/*
  ==============================================================================
    nvs_filters.h
    Created: 6 Dec 2018 12:12:44am
    Author:  Nicholas Solem
  ==============================================================================
*/

/**
 ***TODO:
	-optimize tvap:
		-replace all instances of sin() and cos() with faster versions (pade or lookup table if it remains stable)
**	-remove filter_abstract::z1
	-make lookup tables static
	-make lookup tables constexpr (really in nvs_memoryless)
 */

#pragma once
#include "nvs_memoryless.h"
#include "lookup_tables.h"
#include "nvs_fast_math_impl.h"
#include "nvs_matrix2x2.h"
#include <array>

#include <concepts>
#include <type_traits>

namespace nvs::filters {

template<FloatingPoint float_t>
float_t cutoff_to_g_inaccurate(float_t cutoff, float_t fs_inv){
	return cutoff * fs_inv * 0.5f * (float_t)M_PI;
}
template<FloatingPoint float_t>
float_t cutoff_to_g_slow(float_t cutoff, float_t fs_inv){
	return std::tan(cutoff * fs_inv * (float_t)M_PI);
}

template<FloatingPoint float_t>
float_t cutoff_to_g(float_t cutoff, float_t fs_inv)
{
	auto x = cutoff * fs_inv * (float_t)M_PI;
	return memoryless::math_impl::tan(x);
}

//==================================================================================

enum class mode_e {
	LP = 0,
	HP,
	BP,
	AP,
	NP
};

template<FloatingPoint float_t>
class filter_abstract{
	static_assert(std::is_floating_point<float_t>::value, "filter_abstract type must be floating point");
public:
	virtual ~filter_abstract() = 0;
	//============================================================
	virtual void clear() = 0;
	virtual void setSampleRate(float_t sample_rate){
		_fs_inv = 1.0 / sample_rate;
	}
	float_t getSampleRate() const {
		return 1.0 / _fs_inv;
	}
	virtual void setBlockSize(size_t blockSize){
		_blockSize_inv = static_cast<float_t>(1) / static_cast<float_t>(blockSize);
	}
	//============================================================
	virtual void setCutoffTarget(float_t cutoff_target){
		_cutoffTarget = cutoff_target;
	}
	virtual void setCutoff(float_t cutoff){
		_w_c = cutoff;
		setCutoffTarget(cutoff);
	}
	virtual void updateCutoff(){
		_w_c += (_cutoffTarget - _w_c) * _blockSize_inv;
	}
	virtual float_t getCutoff() const {
		return _w_c;
	}
	virtual void setResonanceTarget(float_t res_target){
		_resonanceTarget = res_target;
	}
	virtual void setResonance(float_t res){
		_q = res;
		setResonanceTarget(res);
	}
	virtual void updateResonance(){
		_q += (_resonanceTarget - _q) * _blockSize_inv;
	}
	virtual float_t getResonance() const {
		return _q;
	}
	virtual void setMode(mode_e mode) = 0;
	//============================================================
	virtual float_t operator()(float_t input) = 0;
	virtual float_t operator()(float_t input, float_t cutoff) = 0;
	virtual float_t operator()(float_t input, float_t cutoff, float_t resonance) = 0;
	
protected:
	float_t _fs_inv;
	float_t _cutoffTarget, _resonanceTarget;
	float_t _w_c, _q;
	float_t _blockSize_inv;
	mode_e _mode {mode_e::LP};
};
template<FloatingPoint float_t>
inline filter_abstract<float_t>::~filter_abstract() { }

/**
 TODO:
**	-get rid of y_n as stored value
	-don't use clamp if unnecessary (and it should be)
 */
template<FloatingPoint float_t>
class onePole   :   public filter_abstract<float_t>
{
public:
	void clear() override {
		v_n = z1 = 0.0;
	}
	//==============================================================================
	/* a one pole filter has no resonance. */
	void setResonanceTarget(float_t) override {}
	void updateResonance() override {}
	//==============================================================================
	
	float_t tpt_lp(float_t x){
		auto const g = cutoff_to_g(this->_w_c, this->_fs_inv);
		auto const G = g / (1.f + g);
		v_n = (x - this->z1) * G;
		auto const y_n = v_n + this->z1;
		this->z1 = y_n + v_n;
		return y_n;
	}
	void setMode(mode_e mode) override {
		if (mode == mode_e::LP)
			this->_mode = mode;
		else	// only 2 modes
			this->_mode = mode_e::HP;
	}
	//==============================================================================
	float_t operator()(float_t input) override {
		auto lp = tpt_lp(input);
		
		if (this->_mode == mode_e::LP)
			return lp;
		else
			return input - lp;
	}
	float_t operator()(float_t input, float_t cutoff) override {
		this->_w_c = cutoff;
		return operator()(input);
	}
	float_t operator()(float_t input, float_t cutoff, float_t) override {
		/*
		 is there a way to warn upon calling a function? but i do not want to deprecate this
		 because it could have some use for generic interfaces
		 */
		return operator()(input, cutoff);
	}
	float_t getState() const {
		return z1;
	}
private:
	float_t v_n {0.0}, z1 {0.0};
};
//==============================================================================

enum class NL_1Pole_Impl {
	Ideal,
	TransistorV2,
	TransistorV3,
	Approx,
	OutputOnly,
	Euler,
	Wrong
};
template<FloatingPoint float_t,
	NL_1Pole_Impl impl=NL_1Pole_Impl::Ideal
>
class onePole_nonlinear_df2   :   public filter_abstract<float_t>
{
public:
	void clear() override {
		z1 = 0.0;
	}
	//=========================================================================
	/* a one pole filter has no resonance. */
	void setResonanceTarget(float_t) override {}
	void updateResonance() override {}
	void setMode(mode_e mode) override {
		if (mode == mode_e::LP)
			this->_mode = mode;
		else if (mode == mode_e::HP)
			this->_mode = mode_e::HP;
		else {
			// assert(false);	// only 2 modes
		}
	}
	//=========================================================================
	float_t tpt_lp_ideal(float_t x) {
		auto const g = cutoff_to_g(this->_w_c, this->_fs_inv);
		auto const G = g / (1.0f + g);  // Still need this for TPT
		using memoryless::math_impl::tanh;

		// We need to solve: u + G*tanh(u) = x + G*tanh(z1)
		// where u is the input to the nonlinearity

		float_t const rhs = x + G * tanh(z1);  // Right-hand side of equation

		// Solve for u using Newton-Raphson (or your preferred method)
		float_t u = u_prev;  // Use previous value as initial guess

		for (int iter = 0; iter < 5; ++iter) {
			float_t tanh_u = tanh(u);
			float_t f = u + G * tanh_u - rhs;
			float_t df = 1.0f + G * (1.0f - tanh_u * tanh_u);

			float_t new_u = u - f / df;
			if (std::abs(new_u - u) < 1e-6f) break;
			u = new_u;
		}
		u_prev = u;  // Store for next iteration

		// Now do TPT integration with the solved input
		float_t tanh_u = tanh(u);
		float_t tanh_z1 = tanh(z1);

		float_t v = G * (tanh_u - tanh_z1);  // TPT difference
		float_t y = v + z1;                   // Output
		z1 = y + v;                           // TPT state update

		return y;
	}
	float_t tpt_approx(float_t x) {
		auto const g = cutoff_to_g(this->_w_c, this->_fs_inv);
		auto const G = g / (1.0f + g);
		using memoryless::math_impl::tanh;

		// Linearize the nonlinearity around the operating point
		float_t tanh_z1 = tanh(z1);
		float_t dtanh_z1 = 1.0f - tanh_z1 * tanh_z1;  // tanh derivative

		// Approximate solution assuming small changes
		float_t u_approx = (x + G * tanh_z1) / (1.0f + G * dtanh_z1);

		// Apply TPT structure
		float_t tanh_u = tanh(u_approx);
		float_t v = G * (tanh_u - tanh_z1);
		float_t y = v + z1;
		z1 = y + v;

		return y;
	}
	float_t tpt_with_output_nonlinearity(float_t x) {
		// Standard linear TPT 1-pole
		auto const g = cutoff_to_g(this->_w_c, this->_fs_inv);
		auto const G = g / (1.0f + g);

		float_t v = (x - z1) * G;
		float_t y_linear = v + z1;
		z1 = y_linear + v;

		// Apply nonlinearity to the output
		using memoryless::math_impl::tanh;
		return tanh(y_linear);  // Nonlinear output
	}
	float_t euler_lp(float_t x) {
		auto const g = cutoff_to_g(this->_w_c, this->_fs_inv);
		using memoryless::math_impl::tanh;

		// Transistor ladder nonlinearity: tanh(x) - tanh(y) drives the integrator
		auto tanh_x = tanh(x);
		auto tanh_y = tanh(z1);  // z1 is the integrator state

		auto integrator_input = g * (tanh_x - tanh_y);
		z1 += integrator_input;  // Simple integration

		return z1;  // Output is the integrator state
	}
	float_t tpt_lp_wrong(float_t x){
		auto const g = cutoff_to_g(this->_w_c, this->_fs_inv);
		auto const x_n = x * g;
		auto const v_n = x_n + z1;
		auto const y_n = memoryless::math_impl::tanh(v_n);
		z1 = x_n + y_n;
		return y_n;
	}
	// Method 2: More accurate transistor differential pair model
	float_t transistor_1pole_v2(float_t x) {
		auto const g = 2 * cutoff_to_g(this->_w_c, this->_fs_inv);
		using memoryless::math_impl::tanh;
		// using juce::dsp::FastMathApproximations::tanh;

		// In a real transistor ladder, we need to solve the implicit equation
		// because the feedback creates a zero-delay loop

		// The equation we're solving:
		// z1_new = z1_old + g * (tanh(x) - tanh(z1_new))
		// Rearranged: z1_new + g*tanh(z1_new) = z1_old + g*tanh(x)

		float_t const rhs = z1 + g * tanh(x);

		// Solve for z1_new using Newton-Raphson
		float_t z1_new = z1;  // Initial guess

		for (int iter = 0; iter < 5; ++iter) {
			float_t tanh_z1_new = tanh(z1_new);
			float_t f = z1_new + g * tanh_z1_new - rhs;
			float_t df = 1.0f + g * (1.0f - tanh_z1_new * tanh_z1_new);

			float_t next_z1 = z1_new - f / df;
			if (std::abs(next_z1 - z1_new) < 1e-6f) break;
			z1_new = next_z1;
		}

		z1 = z1_new;
		return z1;
	}
	float_t transistor_1pole_v3(float_t x) {
		auto const g = 2 * cutoff_to_g(this->_w_c, this->_fs_inv);
		using memoryless::math_impl::tanh;
		// using juce::dsp::FastMathApproximations::tanh;

		// Model the transistor pair more explicitly
		// In analog: I = Is * (tanh(V+/2Vt) - tanh(V-/2Vt))
		// where V+ = input, V- = integrator output

		// Solve for what the integrator output should be
		// given the input and the differential equation

		float_t v_plus = x;      // Input voltage
		float_t v_minus = z1;    // Current integrator state

		// The transconductance creates current: gm * (tanh(v+) - tanh(v-))
		// This current charges the capacitor (integrator)

		// For proper TPT, we need to solve the implicit equation
		float_t target = z1;  // Start with current state

		for (int iter = 0; iter < 3; ++iter) {  // Fewer iterations for real-time
			float_t current = tanh(v_plus) - tanh(target);
			float_t new_target = z1 + g * current;

			if (std::abs(new_target - target) < 1e-5f) break;
			target = new_target;
		}

		z1 = target;
		return z1;
	}
	float_t operator()(float_t x) override {
		auto lp = [this, x]() {
			if constexpr (impl == NL_1Pole_Impl::Ideal) {
				return tpt_lp_ideal(x);
			}
			else if constexpr (impl == NL_1Pole_Impl::TransistorV2) {
				return transistor_1pole_v2(x);
			}
			else if constexpr (impl == NL_1Pole_Impl::TransistorV3) {
				return transistor_1pole_v3(x);
			}
			else if constexpr (impl == NL_1Pole_Impl::Approx) {
				return tpt_approx(x);
			}
			else if constexpr (impl == NL_1Pole_Impl::OutputOnly) {
				return tpt_with_output_nonlinearity(x);
			}
			else if constexpr (impl == NL_1Pole_Impl::Euler) {
				return euler_lp(x);
			}
			else if constexpr (impl == NL_1Pole_Impl::Wrong) {
				return tpt_lp_wrong(x);
			}
		}();
		
		if (this->_mode == mode_e::LP)
			return lp;
		else
			return x - lp;
	}
	float_t operator()(float_t input, float_t cutoff) override {
		this->filter_abstract<float_t>::setCutoff(cutoff);
		return operator()(input);
	}
	float_t operator()(float_t input, float_t cutoff, float_t) override {
		/* is there a way to warn upon calling a function? but i do not want to deprecate this
		 because it could have some use for generic interfaces */
		return operator()(input, cutoff);
	}
	float_t getState() const {
		return z1;
	}
private:
	float_t z1 {0.0};
	float_t u_prev {0.0};
};

template<FloatingPoint float_t>
class fourPole_LP_linear    :   public filter_abstract<float_t>
{
public:
	virtual void setSampleRate(float_t sample_rate) override {
		this->filter_abstract<float_t>::setSampleRate(sample_rate);
		for (auto &pole : _poles){
			pole.setSampleRate(sample_rate);
		}
	}
	virtual void setBlockSize(size_t blockSize) override {
		this->filter_abstract<float_t>::setBlockSize(blockSize);
		for (auto &pole : _poles){
			pole.setBlockSize(blockSize);
		}
	}
	virtual void clear() override {
		for (auto &pole : _poles){
			pole.clear();
		}
	}
	virtual void setCutoffTarget(float_t cutoff_target) override
	{
		this->filter_abstract<float_t>::setCutoffTarget(cutoff_target);
		for (auto &pole : _poles){
			pole.setCutoffTarget(cutoff_target);
		}
	}
	void updateCutoff() override {
		this->filter_abstract<float_t>::updateCutoff();
		for (auto &pole : _poles){
			pole.updateCutoff();
		}
	}
	/*	no need to override updateResonance, setResonanceTarget, because the poles dont care about their resonance value */
	virtual void setMode(mode_e mode) override {
		if (mode == mode_e::LP){
			this->_mode = mode;
			for (auto &p : _poles){
				p.setMode(mode);
			}
		}
		else if (mode == mode_e::HP){
			this->_mode = mode_e::HP;
			for (auto &p : _poles){
				p.setMode(mode);
			}
		}
		else {	// bandpass
			this->_mode = mode_e::BP;
			_poles[0].setMode(mode_e::LP);
			_poles[1].setMode(mode_e::HP);
			_poles[2].setMode(mode_e::LP);
			_poles[3].setMode(mode_e::HP);
		}
	}
	float_t tpt_fourpole(float_t x){
		auto const g = cutoff_to_g(this->_w_c, this->_fs_inv);
		float_t const g2 = g*g;
		float_t const G = g2 * g2;
		float_t const s1 = _poles[0].getState();
		float_t const s2 = _poles[1].getState();
		float_t const s3 = _poles[2].getState();
		float_t const s4 = _poles[3].getState();
		float_t const S = g2*g*s1 + g2*s2 + g*s3 + s4;

		float_t k = this->_q;
		if (this->_mode == mode_e::BP) k *= -1.f;

		// before i was trying to re-emphasize input sig using (x * (1 + this->_q)) instead of x
		float_t u_n = (x - k * S)
					/ (1 + k * G);
		
		float_t y_n = _poles[3]
						(_poles[2]
							(_poles[1]
								(_poles[0]
									(u_n)
								 )
							 )
						 );
		return y_n;
	}
	float_t operator()(float_t input) override {
		return tpt_fourpole(input);
	}
	float_t operator()(float_t input, float_t cutoff) override
	{
		this->_w_c = cutoff;
		for (auto & p : _poles){
			p.setCutoff(cutoff);
		}
		return operator()(input);
	}
	float_t operator()(float_t input, float_t cutoff, float_t reso) override
	{
		this->_q = reso;
		return operator()(input, cutoff);
	}
private:
	std::array<onePole<float_t>, 4> _poles;
};

enum class IterativeMethod {
	FixedPoint = 0,
	NewtonRaphson,
	Bisection,
	None
};
enum class InitialInputGuess {
	Zero = 0,
	Previous
};

template<
	FloatingPoint float_t,
	unsigned int N_iters=16,
	IterativeMethod iterative_method=IterativeMethod::NewtonRaphson,
	bool BisectionConvergenceFallback=true,	// whether to use bisection in case other method doesn't converge
	InitialInputGuess initial_input_guess=InitialInputGuess::Previous,
	class onePole_t=onePole_nonlinear_df2<float_t>
	// class onePole_t=onePole<float_t>
>
class fourPole_LP_nonlinear    :   public filter_abstract<float_t>
{
static_assert(!std::is_same_v<onePole_t, onePole<float_t>>,
		  "OnePole class is not supported - use onePole_nonlinear_df2 instead");
public:
	fourPole_LP_nonlinear() {
		hp.setMode(mode_e::HP);
		lop.setMode(mode_e::LP);
	}
	void setSampleRate(float_t sample_rate) override {
		this->filter_abstract<float_t>::setSampleRate(sample_rate);
		for (auto &pole : _poles){
			pole.setSampleRate(sample_rate);
		}
		hp.setSampleRate(sample_rate);
		// hp.setCutoff(10.f);
		lop.setSampleRate(sample_rate);
		lop.setCutoff(sample_rate / 2.0 - 50.0);
	}
	void setBlockSize(size_t blockSize) override {
		this->filter_abstract<float_t>::setBlockSize(blockSize);
		for (auto &pole : _poles){
			pole.setBlockSize(blockSize);
		}
	}
	void clear() override {
		for (auto &pole : _poles){
			pole.clear();
		}
		hp.clear();
		lop.clear();
	}
	void setCutoff(float_t cutoff) override {
		this->_w_c = cutoff;
		for (auto &pole : _poles){
			pole.setCutoff(cutoff);
		}
		hp.setCutoff(cutoff / 2.0);
		setCutoffTarget(cutoff);
	}
	void setCutoffTarget(float_t cutoff_target) override
	{
		this->filter_abstract<float_t>::setCutoffTarget(cutoff_target);
		for (auto &pole : _poles){
			pole.setCutoffTarget(cutoff_target);
		}
		hp.setCutoffTarget(cutoff_target / 2.0);
	}
	void updateCutoff() override {
		this->filter_abstract<float_t>::updateCutoff();
		for (auto &pole : _poles){
			pole.updateCutoff();
		}
		hp.updateCutoff();
	}
	void setMode(mode_e mode) override {
		// All poles stay in LP mode for the ladder structure
		this->_mode = mode;
		for (auto &p : _poles){
			p.setMode(mode_e::LP);  // Always LP!
		}
	}
	float_t tpt_fourpole(float_t input) {
		constexpr bool use_lop = false;
		constexpr bool use_hp = false;

		auto const g = cutoff_to_g(this->_w_c, this->_fs_inv);
		float_t const g2 = g*g;
		float_t const G = og_way ? g2 * g2
		:							g / (1.0+g);

		float_t const S = [this, G, g, g2]() {
			float_t const s1 = _poles[0].getState();
			float_t const s2 = _poles[1].getState();
			float_t const s3 = _poles[2].getState();
			float_t const s4 = _poles[3].getState();

			return og_way ? g2*g*s1 + g2*s2 + g*s3 + s4		// my way
		:					(G*G*G) / (1.0+g) * s1 +		// way from https://github.com/xodmk/xodVAFilter/blob/master/xodVAFilter.cpp#L36
							(G*G) / (1.0+g) * s2 +
								G / (1.0+g) * s3 +
								1 / (1.0+g) * s4;
		}();
		float_t k = [q = this->_q, w_c = this->_w_c, fs_inv = this->_fs_inv]() {
			[[maybe_unused]] float_t k_og = q;	// probably not accurate

			float_t freq_norm = w_c * fs_inv; // normalized frequency
			[[maybe_unused]] float_t k_freq_dep_lin = q * (1.0f + 2.0f * freq_norm); // increase k at high freq
#pragma message("FURTHER TEST FREQ_DEP K NEAR SELF OSCILLATION POINT")
			[[maybe_unused]] float_t k_freq_dep_pow = q * (1.0f + 1.0 * freq_norm*freq_norm);

			auto _k = k_freq_dep_pow;
			if constexpr (iterative_method == IterativeMethod::None) {
				_k = nvs::memoryless::clamp_high(_k, 4.0);
			}
			return _k;
		}();

		using memoryless::math_impl::tanh;
		bool converged = false;
		if constexpr(initial_input_guess == InitialInputGuess::Zero) {
			u_n = 0.0;	// for newton-raphson, a good guess MIGHT be by using 0 instead of previous. previous is good too. see p.188 of VAFilterDesign_2.1.2.pdf for discussion.
		}	// otherwise it's already set to previous
		if constexpr (iterative_method == IterativeMethod::Bisection) {
			if (use_lop) {
				u_n = lop(u_n);
			}
			return bisection(u_n, input, G, S, k);
		}
		if constexpr (iterative_method==IterativeMethod::None) {
			u_n = og_way ? (input - k*S) / (1.0 + k*G) :
						   (input - k*S) / (1.0 + k*G*G*G*G);
			return process_output(u_n, k);
		}
		for (auto n = 0U; n < N_iters; n++){
			float_t new_u_n {0.0};
			if constexpr (iterative_method == IterativeMethod::FixedPoint) {
				new_u_n = input - k * (G * tanh(u_n) + S);
			}
			else if constexpr (iterative_method == IterativeMethod::NewtonRaphson){	// newton-raphson
				float_t const tanh_u_n = tanh(u_n);
				new_u_n = u_n - (u_n + k * (G * tanh_u_n + S) - input) / (1 + k * G * (1 - tanh_u_n*tanh_u_n));
			}
			// Check convergence
			if (std::abs(new_u_n - u_n) < 1e-6f) {
				converged = true;
				// Debug output occasionally
				// static int counter = 0;
				// if (counter++ % 48000 == 0) {
				// 	std::cout << "Converged in " << n+1 << " iterations" << std::endl;
				// }
				break;
			}
			if (use_lop) {
				u_n = lop(new_u_n);
			} else {
				u_n = new_u_n;
			}
		}
		if (!converged) {
			static int counter = 0;
			if (counter++ % 48000 == 0) {
				std::cout << "Failed to converge!";
			}
			if constexpr (BisectionConvergenceFallback) {
				std::cout << " Using bisection\n";
				return bisection(u_n, input, G, S, k);
			}
			else {std::cout << "\n";}
		}
		return process_output(u_n, k);
	}
	float_t operator()(float_t input) override {
		return tpt_fourpole(input);
	}
	float_t operator()(float_t input, float_t cutoff) override
	{
		this->setCutoff(cutoff);
		return tpt_fourpole(input);
	}
	float_t operator()(float_t input, float_t cutoff, float_t reso) override
	{
		this->_q = reso;
		return operator()(input, cutoff);
	}
private:
	std::array<onePole_t, 4> _poles;
	float_t u_n {0.0}, y_n{0.0};

	static constexpr bool og_way = true;

	onePole<float_t> hp;
	onePole<float_t> lop;
	float_t process_output(float_t u_n, float_t k) {
		using memoryless::math_impl::tanh;

		float_t y0 = /*tanh*/(u_n);
		float_t y1 = _poles[0](u_n);
		float_t y2 = _poles[1](y1);
		float_t y3 = _poles[2](y2);
		float_t y4 = _poles[3](y3);

		if (this->_mode == mode_e::LP) {
			y_n = y4;
		}
		else if (this->_mode == mode_e::BP) {
			y_n = y2 - 2*y3 + y4;
		}
		else {
			y_n = 0.125 * (y0 - 4*y1 + 6*y2 - 4*y3 + y4);
		}

		assert(y_n == y_n);    // not NaN
		assert(!isinf(y_n));

		auto const gain_comp = 1.0f + k * 0.25f;
		return y_n * gain_comp;
	}
	float_t bisection(float_t u_n, float_t input, float_t G, float_t S, float_t k){
	    using memoryless::math_impl::tanh;

	    // Define the function we want to find the root of
	    auto f = [&](float_t u) -> float_t {
	        return u + k * (G * tanh(u) + S) - input;
	    };

	    // Find initial bracketing range [a, b]
	    float_t a, b;

	    // Method from PDF: use bounds based on kG*sgn(u) vs kG*tanh(u)
	    // The line (input - kS) - u intersects kG*sgn(u) at the bracket points
	    float_t x_minus_kS = input - k * S;

	    if (k >= 0) {
	        // For positive k, we can bracket using the fact that:
	        // -kG <= kG*tanh(u) <= kG
	        // So we solve: x_minus_kS - u = ±kG

	        if (G > 0) {
	            a = x_minus_kS + k * G;  // Left bracket
	            b = x_minus_kS - k * G;  // Right bracket

	            // Ensure a < b (swap if needed)
	            if (a > b) std::swap(a, b);

	            // Expand brackets slightly to be safe
	            float_t range = b - a;
	            a -= 0.1f * range;
	            b += 0.1f * range;
	        } else {
	            // Fallback for edge cases
	            a = -10.0f;
	            b = 10.0f;
	        }
	    } else {
	        // For negative k, use a more conservative approach
	        a = -10.0f;
	        b = 10.0f;
	    }

	    // Verify initial bracket is valid
	    float_t fa = f(a);
	    float_t fb = f(b);

	    // If same sign, expand search
	    if (fa * fb > 0) {
	        // Try to find a valid bracket by expanding
	        for (int expand = 0; expand < 5; ++expand) {
	            float_t range = b - a;
	            a -= range;
	            b += range;
	            fa = f(a);
	            fb = f(b);
	            if (fa * fb <= 0) break;
	        }

	        // If still no valid bracket, fall back to previous value or zero
	        if (fa * fb > 0) {
	            u_n = 0.0f; // or keep previous u_n
	            return process_output(u_n, k);
	        }
	    }

	    // Bisection iteration
	    constexpr int MAX_BISECTION_ITERS = 12; // ~60dB precision as mentioned in
	    // constexpr int MAX_BISECTION_ITERS = 24;
	    constexpr float_t TOLERANCE = 2e-6f;

	    for (int n = 0; n < MAX_BISECTION_ITERS; ++n) {
	        float_t c = 0.5f * (a + b);  // Midpoint
	        float_t fc = f(c);

	        // Check convergence
	        if (std::abs(fc) < TOLERANCE || (b - a) * 0.5f < TOLERANCE) {
	            u_n = c;

	            // Debug output occasionally
	            // static int counter = 0;
	            // if (counter++ % 48000 == 0) {
	            //     std::cout << "Bisection converged in " << n+1 << " iterations" << std::endl;
	            // }
	            break;
	        }

	        // Update bracket
	        if (fa * fc < 0) {
	            b = c;
	            fb = fc;
	        } else {
	            a = c;
	            fa = fc;
	        }
	    }
	    return process_output(u_n, k);
	}
};

template<FloatingPoint float_t>
class butterworth2p :   public filter_abstract<float_t>
{
public:
	butterworth2p() {
		clear();
		A.b = 1.0;
		B.b = 1.0;
		D = 0.0;
	}
	void clear() override {
		x.a = x.b = 0.f;
		y = 0.f;
	}
	//============================================================
	virtual void updateCutoff() override {
		this->filter_abstract<float_t>::updateCutoff();
		calcCoefs(this->_w_c);
	}
	// no altering Q for butterworth
	virtual void updateResonance() override {}
	virtual void setResonanceTarget(float_t) override {}
	
	virtual void setMode(mode_e) override {}
	
	void calcCoefs(const float_t cutoff){
		const float_t omega = cutoff_to_g(cutoff, this->_fs_inv);
		const float_t omega2 = omega * omega;
		const float_t twoCosPiOver4xOmega = (float_t)1.414213562373095 * omega;
		const float_t c = 1.f + twoCosPiOver4xOmega + omega2;
		
		const float_t b0 = omega2 / c;
		const float_t b1 = 2.f * b0;
		const float_t a1 = (2.f * (omega2 - 1.f)) / c;
		const float_t a2 = (1.f - twoCosPiOver4xOmega + omega2) / c;
		
		A.c = -a2;
		A.d = -a1;
		C.a = b0 - (a2*b0);
		C.b = b1 - (a1*b0);
		D = b0;
	}
	virtual float_t operator()(float_t x_n) override {
		using namespace nvs_matrix;
		y = vec2::crossProduct(C, x);
		y += D * x_n;
		x = vec2::add(mat2x2::matXvec(A, x), vec2::scale(B, x_n));
		
		return y;
	}
	virtual float_t operator()(float_t x_n, float_t cutoff) override {
		calcCoefs(cutoff);
		return operator()(x_n);
	}
	virtual float_t operator()(float_t x_n, float_t cutoff, float_t) override {
		/* possible to warn for this functions use? */
		return operator()(x_n, cutoff);
	}
	
private:
	nvs_matrix::mat2x2 A;
	nvs_matrix::vec2 B;
	nvs_matrix::vec2 C;
	float_t D;
	nvs_matrix::vec2 x;
	float_t y;
};

//==================================================================================
template<FloatingPoint float_t>
class dcBlock   :   public filter_abstract <float_t>
{
public:
	void clear() override {
		_xz1 = _yz1 = 0.f;
	 }
	////
	virtual void setCutoffTarget(float_t) override {}
	virtual void setCutoff(float_t) override {}
	virtual void updateCutoff() override {}
	virtual void setResonanceTarget(float_t ) override {}
	virtual void setResonance(float_t) override {}
	virtual void updateResonance() override {}
	virtual void setMode(mode_e) override {}

	void setR(float_t R){
		R = nvs::memoryless::clamp<float_t>(R, (float_t)0.f, (float_t)1.f);
		this->_R = R;
	}
	void updateR(float_t R_target, float_t oneOverBlockSize){
		_R += (R_target - _R) * oneOverBlockSize;
	}
	float_t filter(float_t x)
	{// y[i] = x[i] - x[i - 1] + R * y[i - 1];
		_yz1 = x - _xz1 + _R * _yz1;
		_xz1 = x;
		return _yz1;
	 }
	float_t operator()(float_t input) override {
		return filter(input);
	}
	float_t operator()(float_t input, float_t) override {
		return operator()(input);
	}
	float_t operator()(float_t input, float_t, float_t) override {
		return operator()(input);
	}
private:
	float_t _R{0.995f}, _xz1{0.f}, _yz1{0.f};
};


//==================================================================================
namespace {
template<typename>
inline constexpr bool always_false = false;
}

template<FloatingPoint float_t>
class svf_prototype
{
public:
	float_t lp() {return _outputs.lp;}
	float_t bp() {return _outputs.bp;}
	float_t hp() {return _outputs.hp;}
	float_t np() {return _outputs.np;}
protected:
	struct outputs
	{
		float_t lp, bp, hp, np;
	} _outputs = { 0.f, 0.f, 0.f, 0.f };
	struct state
	{
		float_t lp, bp;
	} _state = { 0.f, 0.f };
	
};
// linear state variable filter using 'naive' integrators (i.e., Euler backward difference integration)
template<FloatingPoint float_t>
class svf_lin_naive     :   public filter_abstract<float_t>, svf_prototype<float_t>
{
public:
	//==============================================================================
	void clear() override {
		this->_outputs = {.lp{0.f}, .bp{0.f}, .hp{0.f}, .np{0.f} };
		this->_state = { .lp{0.f}, .bp{0.f} };
	}
	virtual void setResonanceTarget(float_t res_target) override {
		if (res_target > 0.9f){
			res_target = 0.9f;
		}
		this->filter_abstract<float_t>::setResonanceTarget(res_target);
	}
	void filter(float_t input){
		float_t c, d;
		c = 2.f * memoryless::math_impl::sin((float_t)M_PI * this->_w_c * this->_fs_inv);
		d = 2.f * (1.f - pow(this->_q, 0.25f));

		if (c > 0.5f) c = 0.5f;
		if (d > 2.0f) d = 2.f;
		if (d > (2.f/c - (c * 0.5f)))
		 d = 2.f/c - (c * 0.5f);

		this->_outputs.np = input - (d * this->_outputs.bp);
		this->_outputs.lp = this->_outputs.lp + (c * this->_outputs.bp);
		this->_outputs.hp = this->_outputs.np - this->_outputs.lp;
		this->_outputs.bp = this->_outputs.bp + (c * this->_outputs.hp);
	}
	void setMode(mode_e mode) override {
		if (mode == (mode_e::AP)){
			mode = mode_e::NP;
		}
		this->_mode = mode;
	}
	//============================================================
	float_t operator()(float_t input) override {
		filter(input);
		switch (static_cast<int>(this->_mode)) {
			case static_cast<int>(mode_e::LP):
				return this->_outputs.lp;
				break;
			case static_cast<int>(mode_e::HP):
				return this->_outputs.hp;
				break;
			case static_cast<int>(mode_e::BP):
				return this->_outputs.bp;
				break;
			case static_cast<int>(mode_e::NP):
				return this->_outputs.np;
				break;
		}
		return 0.f;
	}
	float_t operator()(float_t input, float_t cutoff) override {
		this->setCutoff(cutoff);
		return operator()(input);
	}
	float_t operator()(float_t inp, float_t cut, float_t res) override {
		this->setResonance(res);
		return operator()(inp, cut);
	}
};

//==================================================================================
/*
 nonlinear state-variable filter using fourth-order runge-kutta
 y[n+1] = y[n] + 1/6 * (k_1 + 2*k_2 + 2*k_3 + k_4)
 h = (1 / fs) / oversample_factor
 
 k_1 = h*f(t_n, y_n)
 k_2 = h*f(t_n + h/2, y_n + k_1/2)
 k_3 = h*f(t_n + h/2, y_n + k_2/2)
 k_4 = h*f(t_n + h, y_n + k_3)
 */
enum class character_e {
	correct = 0,
	wrong = 1
};

template<FloatingPoint float_t>
class svf_nl_rk :   public filter_abstract<float_t>, public svf_prototype<float_t>
{
	using base = filter_abstract<float_t>;
	
public:
	svf_nl_rk() {
		for (int i = 0; i < 2; i++) {
			deriv1[i] = deriv2[i] = deriv3[i] = deriv4[i] = 0.f;
		}
	}
	void setSampleRate(float_t sample_rate) override {
		this->base::setSampleRate(sample_rate);
		_h = 1.f / (_oversample_factor * sample_rate);
	}
	void set_oversample(unsigned int oversample_factor){
		_oversample_factor = oversample_factor;
//		_h = 1.f / (oversample_factor * this->sampleRate);
		_h = this->_fs_inv / oversample_factor;
	}
	void clear() override {
		for (int i = 0; i < 2; i++) {
			deriv1[i] = deriv2[i] = deriv3[i] = deriv4[i] = 0.f;
		}
		this->_outputs = {0.f, 0.f, 0.f, 0.f };
		this->_state = { 0.f, 0.f };
	}
	
	void setResonance(float_t res) override {
		res *= 3.f;
		if (res < 0.5f) res = 0.5f;
		this->base::setResonance(res);
		_resInv = 1.f / res;
	}
	void setResonanceTarget(float_t res_target) override {
		res_target *= 3.f;
		if (res_target < 0.5f)
			res_target = 0.5f;
		this->base::setResonanceTarget(res_target);
	}
	void updateResonance() override {
		this->base::updateResonance();
		_resInv = 1.f / this->_q;
	}
	void setMode(mode_e mode) override {
		if (mode == (mode_e::AP)){
			mode = mode_e::NP;
		}
		this->_mode = mode;
	}
	void setCharacter(character_e character) {
		assert ((character == character_e::wrong) || (character == character_e::correct));
		_character = character;
	}
	void filter(float_t input){
		using namespace nvs::memoryless;
		float_t hp(0), np(0);
		// overwritten states. [0] is bp, [1] is lp.
		float_t tempstate[2];

		auto constexpr twopi = math_impl::two_pi<float_t>();
		auto const res_inv = _resInv * 0.09;
		
		for (unsigned iter = 0; iter < _oversample_factor; iter++) {
			np = input - 2 * res_inv * this->_state.bp;
			hp = np - this->_state.lp;
			deriv1[0] = _h * twopi * this->_w_c * math_impl::tanh(hp);
			deriv1[1] = _h * twopi * this->_w_c * math_impl::tanh(this->_state.bp);
			tempstate[0] = this->_state.bp + deriv1[0] / 2;
			tempstate[1] = this->_state.lp + deriv1[1] / 2;
			
			np = input - 2 * res_inv * tempstate[0];
			hp = np - tempstate[1];
			deriv2[0] = _h * twopi * this->_w_c * math_impl::tanh(hp);
			deriv2[1] = _h * twopi * this->_w_c * math_impl::tanh(tempstate[0]);
			tempstate[0] = this->_state.bp + deriv2[0] / 2;
			tempstate[1] = this->_state.lp + deriv2[1] / 2;
			
			np = input - 2 * res_inv * tempstate[0];
			hp = np - tempstate[1];
			deriv3[0] = _h * twopi * this->_w_c * math_impl::tanh(hp);
			deriv3[1] = _h * twopi * this->_w_c * math_impl::tanh(tempstate[0]);
			tempstate[0] = this->_state.bp + deriv3[0];
			tempstate[1] = this->_state.lp + deriv3[1];
			
			np = input - 2 * res_inv * tempstate[0];
			hp = np - tempstate[1];
			deriv4[0] = _h * twopi * this->_w_c * math_impl::tanh(hp);
			deriv4[1] = _h * twopi * this->_w_c * math_impl::tanh(tempstate[0]);
			this->_state.bp += (1.f/6.f) * (deriv1[0] + 2 * deriv2[0] + 2 * deriv3[0] + deriv4[0]);
			this->_state.lp += (1.f/6.f) * (deriv1[1] + 2 * deriv2[1] + 2 * deriv3[1] + deriv4[1]);
			
			np = input - 2 * res_inv * this->_state.bp;
			hp = np - this->_state.lp;
		}
		this->_outputs.bp = this->_state.bp;
		this->_outputs.lp = this->_state.lp;
		this->_outputs.hp = hp;
		this->_outputs.np = np;
	}
	void filter_wrong(float_t input){
		/*
		 this version made some mistakes (not updating hp and notch at each stage),
		 but led to a very unique chaotic character, so i keep it as an option.
		 */
		using namespace nvs::memoryless;
		float_t hp(0), np(0);
		// overwritten states. [0] is bp, [1] is lp.
		float_t tempstate[2];
		
		auto constexpr twopi = math_impl::two_pi<float_t>();

		for (unsigned iter = 0; iter < _oversample_factor; iter++)
		{
			np = input - 2 * _resInv * this->_state.bp;
			hp = np - this->_state.lp;
			
			auto const hp_for_deriv = math_impl::tanh(hp);
			
			deriv1[0] = _h * twopi * this->_w_c * hp_for_deriv;
			deriv1[1] = _h * twopi * this->_w_c * math_impl::tanh(this->_state.lp);

			tempstate[0] = this->_state.bp + deriv1[0] / 2;
			tempstate[1] = this->_state.lp + deriv1[1] / 2;

			deriv2[0] = _h * twopi * this->_w_c * hp_for_deriv;
			deriv2[1] = _h * twopi * this->_w_c * math_impl::tanh(tempstate[0]);

			tempstate[0] = this->_state.bp + deriv2[0] / 2;
			tempstate[1] = this->_state.lp + deriv2[1] / 2;

			deriv3[0] = _h * twopi * this->_w_c * hp_for_deriv;
			deriv3[1] = _h * twopi * this->_w_c * math_impl::tanh(tempstate[0]);

			tempstate[0] = this->_state.bp + deriv3[0];
			tempstate[1] = this->_state.lp + deriv3[1];

			deriv4[0] = _h * twopi * this->_w_c * hp_for_deriv;
			deriv4[1] = _h * twopi * this->_w_c * math_impl::tanh(tempstate[0]);

			this->_state.bp += (1.f/6.f) * (deriv1[0] + 2 * deriv2[0] + 2 * deriv3[0] + deriv4[0]);
			this->_state.lp += (1.f/6.f) * (deriv1[1] + 2 * deriv2[1] + 2 * deriv3[1] + deriv4[1]);
		}

		this->_outputs.bp = this->_state.bp;
		this->_outputs.lp = this->_state.lp;
		this->_outputs.hp = hp;
		this->_outputs.np = np;
	}
	
	float_t operator()(float_t input) override {
		if (_character == character_e::correct){
			filter(input);
		}
		else if (_character == character_e::wrong){
			filter_wrong(input);
		}
		switch (static_cast<int>(this->_mode)) {
			case static_cast<int>(mode_e::LP):
				return this->_outputs.lp;
				break;
			case static_cast<int>(mode_e::HP):
				return this->_outputs.hp;
				break;
			case static_cast<int>(mode_e::BP):
				return this->_outputs.bp;
				break;
			case static_cast<int>(mode_e::NP):
				return this->_outputs.np;
				break;
		}
		return 0.f;
	}
	float_t operator()(float_t input, float_t cutoff) override {
		this->base::setCutoff(cutoff);
		return operator()(input);
	}
	float_t operator()(float_t inp, float_t cut, float_t res) override {
		setResonance(res);
		return operator()(inp, cut);
	}
private:
	unsigned int _oversample_factor {4};
	float_t _h {1.f / (_oversample_factor*44100.f)}, _resInv {1.f};
	// k_1 through k_4. for each, [0] is bp, [1] is lp.
	float_t deriv1[2], deriv2[2], deriv3[2], deriv4[2];
	
	character_e _character = character_e::correct;
};
//=======================================================================
/*
 time-variant allpass filter
 */
template<FloatingPoint float_t>
class tvap  :   public filter_abstract<float_t>
{
public:
	//===============================================================
	void setSampleRate(float_t sample_rate) override {
		this->filter_abstract<float_t>::setSampleRate(sample_rate);
		dcFilt.setSampleRate(sample_rate);	// doesn't matter for current implementation but may in the future
		lp.setSampleRate(sample_rate);
		lp.setCutoff(sample_rate * 0.125f);
	}
	void setMode(mode_e) override {
		// ignore request, always allpass
		this->_mode = mode_e::AP;
	}
	void clear() override {
		sp->z1 = sp->z2 = sp->fb_proc = 0.f;
	}
	//===============================================================

	void updateCutoff() override {
		update_f_pi();
	}
	void updateResonance() override {
		update_f_b();
	}
	// function aliases just to have more meaningful names
	void update_f_pi (){
		memoryless::clamp_low(this->_cutoffTarget, (float_t)0.f);
		f_pi += (this->_cutoffTarget - f_pi) * this->_blockSize_inv;
		calc_b1();
	}
	void update_f_b(){
		//if (f_b <= 0) f_b = 0.0000000001;
//		constexpr float_t min_f_b = std::numeric_limits<float_t>::
		
		f_b += (this->_resonanceTarget - f_b) * this->_blockSize_inv;
		f_b_to_b0();
	}
	
	void calc_b1(){
		float_t d = -1 * cos((2.f * M_PI * f_pi) * this->_fs_inv);
		float_t c = (tan(M_PI * f_b * this->_fs_inv) - 1.f) / (tan(M_PI * f_b * this->_fs_inv) + 1.f);
		float_t r1 = acos(-1.f * c);
		float_t r2 = acos(-1.f * d);
		b1 = cos(r1) * (1.f + cos(r2));
	}
	void f_b_to_b0(){
		float_t c = (tan(M_PI * f_b * this->_fs_inv) - 1.f) / (tan(M_PI * f_b * this->_fs_inv) + 1.f);
		float_t r1 = acos(-1.f * c);
		b0 = cos(r1);
	}
	
	float_t f_pi2r2(float_t _f_pi){
		float_t d = -1 * cos((2.f * M_PI * _f_pi) * this->_fs_inv);
		float_t r2 = acos(-d);
		return r2;
	}
	float_t f_b2r1(float_t _f_b){
		float_t tmp = tan(M_PI * _f_b * this->_fs_inv);
		float_t c = (tmp - 1.f) / (tmp + 1.f);
		float_t r1 = acos(-c);
		return r1;
	}
	
	float_t filter(float_t x_n) {
		/* float_t _y1 = state.y1;
		float_t _y2 = state.y2;
		float_t _x1 = state.x1;
		float_t _x2 = state.x2;
		float_t y_n = b0 * x_n - b1 * _x1 + _x2 + b1 * _y1 - b0 * _y2;
		state.y2 = _y1;
		state.y1 = y_n;
		state.x2 = _x1;
		state.x1 = x_n;
		return y_n; */
		float_t _r1, _r2, _cr1, _cr2, _sr1, _sr2;
		float_t tmp [3];
		_r1 = f_b2r1(f_b);
		_r2 = f_pi2r2(f_pi);

		_cr1 = cos(_r1);
		_sr1 = memoryless::math_impl::sin(_r1);
		_cr2 = cos(_r2);
		_sr2 = memoryless::math_impl::sin(_r2);
		//tmp[0] = _x_n;
		tmp[1] = _cr2 * state.z1 - _sr2 * state.z2;
		//tmp[2] = _sr2 * _z1 + _cr2 * _z2;

		tmp[0] = _cr1 * x_n - _sr1 * tmp[1];
		tmp[1] = _sr1 * x_n + _cr1 * tmp[1];
		tmp[2] = _sr2 * state.z1 + _cr2 * state.z2;

		state.z1 = tmp[1];
		state.z2 = tmp[2];

		return tmp[0];
	}
	
	// should be in memoryless but got linker error
	float_t unboundSat2(float_t x){
		float_t num = 2.0 * x;
		float_t denom = 1.0 + sqrt(1.0 + abs(4.0 * x));
		return num / denom;
	}
	
	float_t filter_fbmod(float_t x_n, float_t fb_f_pi, float_t fb_f_b){
		float_t _f_pi_n, _f_b_n;
		_f_pi_n = f_pi;
		_f_b_n = f_b;

		fb_f_pi *= f_pi;
		fb_f_pi *= 1000.0;
		fb_f_b *= f_pi;
		fb_f_b *= 1000.0;

		_f_pi_n += unboundSat2(state.fb_proc * fb_f_pi)* 10.0;
		_f_b_n  += unboundSat2(state.fb_proc * fb_f_b) * 10.0;

		// if (_f_b_n < 0.f) _f_b_n = 0.f;
		if (_f_b_n < 0.0)
		_f_b_n = exp(_f_b_n);
		else
		_f_b_n += 1.0;

		float_t highLimit = (1.f / this->_fs_inv) * 0.495f;
		if (_f_pi_n >= highLimit) _f_pi_n = highLimit;
		if (_f_b_n >= highLimit) _f_b_n = highLimit;

		float_t _r1, _r2, _cr1, _cr2, _sr1, _sr2;
		float_t tmp[3];
		_r1 = f_b2r1(_f_b_n);
		_r2 = f_pi2r2(_f_pi_n);

		_cr1 = cos(_r1);
		_sr1 = memoryless::math_impl::sin(_r1);
		_cr2 = cos(_r2);
		_sr2 = memoryless::math_impl::sin(_r2);
		//tmp[0] = _x_n;
		tmp[1] = _cr2 * state.z1 - _sr2 * state.z2;
		//tmp[2] = _sr2 * _z1 + _cr2 * _z2;

		tmp[0] = _cr1 * x_n - _sr1 * tmp[1];
		tmp[1] = _sr1 * x_n + _cr1 * tmp[1];
		tmp[2] = _sr2 * state.z1 + _cr2 * state.z2;

		state.z1 = tmp[1];
		state.z2 = tmp[2];


		float_t fb_filt = dcFilt(tmp[0]);   // used to feed output to modulate control inputs
		fb_filt = lp(fb_filt);
		state.fb_proc = fb_filt;

		return tmp[0];
	}
	float_t operator()(float_t input) override {
		return filter(input);
	}
	float_t operator()(float_t inp, float_t cut) override {
		this->setCutoff(cut);
		return operator()(inp);
	}
	float_t operator()(float_t inp, float_t cut, float_t reso) override {
		this->setResonance(reso);
		return operator()(inp, cut);
	}
protected:
	typedef struct tvapstate {
		//float_t x1, x2, y1, y2;
		float_t z1, z2;
		float_t fb_proc;  // processed fed back output sample
	} _tvapstate;
	_tvapstate state = {.z1 = 0.f, .z2 = 0.f, .fb_proc = 0.f};
	_tvapstate *sp = &state;
	
	dcBlock<float_t> dcFilt;
	onePole<float_t> lp;
	
private:
	float_t f_pi, f_b;
	float_t b0, b1;
};
//===============================================================================
// PIRKLE IMPLEMENTATIONS (not my own work; used only for checking.)
template<FloatingPoint float_t>
class CTPTMoogFilterStage
{
public:
	CTPTMoogFilterStage(){}
	~CTPTMoogFilterStage(){}
protected:
	float_t G;
	float_t scalar;
	float_t sampleRate;
	float_t z1;
public:
	inline void initialize(float_t newSampleRate){
		// save
		sampleRate = newSampleRate;
		z1 = 0;
	}
	void setFc(float_t fc){
		// prewarp the cutoff- these are bilinear-transform filters
		float_t wd = 2 * M_PI * fc;
		float_t fs_inv  = 1 / sampleRate;
		float_t wa = (2 / fs_inv) * tan(wd * fs_inv / 2);
		float_t g  = wa * fs_inv / 2;
		// calculate big G value; see Zavalishin p46 The Art of VA Design
		G = g / (1.0 + g);
	}
	
	float_t doFilterStage(float_t xn){
		float_t const v = (xn - z1) * G;
		float_t const _out = v + z1;
		z1 = _out + v;
		return _out;
	}
	float_t getSampleRate(){
		return sampleRate;
	}
	float_t getStorageRegisterValue(){
		return z1;
	}
};
template <FloatingPoint float_t>
class CTPTMoogLadderFilter
{
//public:
//	CTPTMoogLadderFilter(){}
//	~CTPTMoogLadderFilter(){}
protected:
	CTPTMoogFilterStage<float_t> filter1;
	CTPTMoogFilterStage<float_t> filter2;
	CTPTMoogFilterStage<float_t> filter3;
	CTPTMoogFilterStage<float_t> filter4;
	float_t k; // Q control
	float_t fc; // fc control
public:
	inline void initialize(float_t newSampleRate){
		filter1.initialize(newSampleRate);
		filter2.initialize(newSampleRate);
		filter3.initialize(newSampleRate);
		filter4.initialize(newSampleRate);
	}
	inline void calculateTPTCoeffs(float_t cutoff, float_t Q){
		// 4 sync-tuned filters
		filter1.setFc(cutoff);
		filter2.setFc(cutoff);
		filter3.setFc(cutoff);
		filter4.setFc(cutoff);

		// NOTE: Q is limited to 20 on the UI to prevent blowing up
				// Q=0.707->25==>k=0->4
		k = Q;	//4.0*(Q - 0.707)/(25.0 - 0.707);
		// ours
		fc = cutoff;
	}
	float_t doTPTMoogLPF(float_t xn){	// calculate g
		float_t const wd = 2 * M_PI * fc;
		float_t const fs_inv  = 1 / (float_t)filter1.getSampleRate();
		float_t const wa = (2 / fs_inv) * tan(wd * fs_inv / 2);
		float_t const g  = wa * fs_inv / 2;
		
		float_t const G = g * g * g * g;
		float_t const S = g * g * g * filter1.getStorageRegisterValue() +
							  g * g * filter2.getStorageRegisterValue() +
								  g * filter3.getStorageRegisterValue() +
									  filter4.getStorageRegisterValue();
		
		//uis input to filters, straight from book
		float_t u = (xn - k * S) / (1 + k * G);
		// four cascades using nested functions
		float_t filterOut = filter4.doFilterStage
										(filter3.doFilterStage
											(filter2.doFilterStage
												(filter1.doFilterStage
													(u)
												 )
											 )
										 );
		// output
		return filterOut;
	}
};

//==============================================================================
template<FloatingPoint float_t>
class slewlim
{
public:
	~slewlim()  { }
	//============================================================
	void setSampleRate(float_t sample_rate){
		assert (sample_rate > 0.f);
		fs_inv = 1.f / sample_rate;
	}
	float_t getSampleRate() const {
		return 1.0 / fs_inv;
	}
	void setBlockSize(size_t blockSize){
		_oneOverBlockSize = 1.f / static_cast<float_t>(blockSize);
	}
	//============================================================
	// immediate change
	void setRise(float_t r){
		this->rise = r;
		this->riseInc = (this->fs_inv * 1000.f) / (this->rise);
		setRiseTarget(r);
	}
	void setRiseTarget(float_t rise_target){
		this->_riseTarget = rise_target;
	}
	void updateRise(){ // change over block size
		if (_riseTarget != rise)
		{
			this->rise += (_riseTarget - rise) * _oneOverBlockSize;
			this->riseInc = (this->fs_inv * 1000.f) / (rise);
		}
	}
	// immediate change
	void setFall(float_t f){
		this->fall = f;
		this->fallInc = (this->fs_inv * 1000.f) / (this->fall);
		setFallTarget(f);
	}
	void setFallTarget(float_t fall_target){
		this->_fallTarget = fall_target;
	}
	// change over block size
	void updateFall(){
	if (_fallTarget != fall)
		{
			this->fall += (_fallTarget - fall) * _oneOverBlockSize;
			this->fallInc = (this->fs_inv * 1000.f) / (fall);
		}
	}
	//============================================================
	float_t ASR(float_t gate){
		using namespace nvs::memoryless;
		
		// it was updating the params implicitly before...
//		setRise();
//		setFall();
		
		if (_vOut < gate)
		{
			_vOut += riseInc;
			_vOut = clamp_high<float_t>(_vOut, gate);
		}
		else if (_vOut > gate)
		{
			_vOut -= fallInc;
			_vOut = clamp_low<float_t>(_vOut, gate);
		}
		
		return _vOut;
	}
	
private:
	float_t sampleRate, fs_inv;
	
	float_t _riseTarget, _fallTarget;
	float_t _oneOverBlockSize;
	// 'Inc' variables tell change per sample.
	float_t rise, riseInc, fall, fallInc, _vOut{0.f};
};

}   // namespace nvs::filters
