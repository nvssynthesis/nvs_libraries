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

namespace nvs	{
namespace filters {

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
	auto x = cutoff * fs_inv / 2.0;
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
		y_n = v_n = z1 = 0.0;
	}
	//==============================================================================
	/* a one pole filter has no resonance. */
	virtual void setResonanceTarget(float_t) override {}
	virtual void updateResonance() override {}
	//==============================================================================
	
	float_t tpt_lp(float_t x){
		auto const g = cutoff_to_g(this->_w_c, this->_fs_inv);
		auto const G = g / (1.f + g);
		v_n = (x - this->z1) * G;
		y_n = v_n + this->z1;
		this->z1 = y_n + v_n;
		return nvs::memoryless::clamp<float_t>(y_n, -1000.0, 1000.0);
	 }
	virtual void setMode(mode_e mode) override {
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
	float_t getState(){
		return z1;
	}
private:
	float_t v_n {0.0}, y_n {0.0}, z1 {0.0};
};
//==============================================================================

template<FloatingPoint float_t>
class onePole_nonlinear_df2   :   public filter_abstract<float_t>
{
public:
	void clear() override {
		v_n = z1 = 0.0;
	}
	//=========================================================================
	/* a one pole filter has no resonance. */
	virtual void setResonanceTarget(float_t) override {}
	virtual void updateResonance() override {}
	virtual void setMode(mode_e mode) override {
		if (mode == mode_e::LP)
			this->_mode = mode;
		else	// only 2 modes
			this->_mode = mode_e::HP;
	}
	//=========================================================================
	float_t tpt_lp(float_t x){
		auto const g = cutoff_to_g(this->_w_c, this->_fs_inv);
		auto x_n = x * g;
		v_n = x_n + z1;
		auto y_n = memoryless::math_impl::tanh(v_n);
		z1 = x_n + y_n;
		return y_n;
	}

	float_t operator()(float_t input) override {
		auto lp = tpt_lp(input, {});
		
		if (this->_mode == mode_e::LP)
			return lp;
		else
			return input - lp;
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
	float_t getState(){
		return z1;
	}
private:
	float_t v_n {0.0}, z1 {0.0};
};

/**
 TODO:
	-figure out power loss with non-LP modes
 */
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
		return tpt_fourpole(input, {});
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

// so far quite tame, since the only nonlinearity is in the feedback path. 
/** TODO:
	-fix DC offset that arises after ANY Q has been introduced
	-WHY DID THAT PROBLEM JUST STOP? what did i change?
	-and also, now in highpass mode, it explodes easily, and when you reload the plugin, it remains silent... in reaper at least...
 */
template<
	FloatingPoint float_t,
	unsigned int N_iters=16,
	class onePole_t=onePole_nonlinear_df2<float_t>
>
class fourPole_LP_nonlinear    :   public filter_abstract<float_t>
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
	
	float_t tpt_fourpole(float_t input){
		auto const g = cutoff_to_g(this->_w_c, this->_fs_inv);
		float_t const g2 = g*g;
		float_t const G = g2 * g2;
		float_t const s1 = _poles[0].getState();
		float_t const s2 = _poles[1].getState();
		float_t const s3 = _poles[2].getState();
		float_t const s4 = _poles[3].getState();
		float_t const S = g2*g*s1 + g2*s2 + g*s3 + s4;
		
		float_t k = this->_q;
		if (this->_mode == mode_e::BP) {k *= -1.f;}

		u_n = y_n;  // initial estimation
		for (auto n = 0U; n < N_iters; n++){
			float_t tmp = nvs::memoryless::clamp(u_n, -100.f, 100.f);
			u_n = input - k * (G * memoryless::math_impl::tanh(tmp) + S);
			assert(u_n == u_n);
			assert(!isinf(u_n));
		}
		y_n = _poles[3]
				(_poles[2]
					(_poles[1]
						(_poles[0]
							(u_n)
						 )
					 )
				 );
		assert(y_n == y_n);
		assert(!isinf(y_n));
//		y_n = memoryless::math_impl::tanh(y_n);
		return y_n;
	}
	float_t operator()(float_t input) override {
		return tpt_fourpole(input, {});
	}
	float_t operator()(float_t input, float_t cutoff) override
	{
		this->_w_c = cutoff;
		for (auto & p : _poles){
			p.setCutoff(cutoff);
		}
		return tpt_fourpole(input, {});
	}
	float_t operator()(float_t input, float_t cutoff, float_t reso) override
	{
		this->_q = reso;
		return operator()(input, cutoff);
	}
private:
	std::array<onePole_t, 4> _poles;
	float_t u_n, y_n;
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
		
		for (unsigned iter = 0; iter < _oversample_factor; iter++) {
			np = input - 2 * _resInv * this->_state.bp;
			hp = np - this->_state.lp;
			deriv1[0] = _h * twopi * this->_w_c * math_impl::tanh(hp);
			deriv1[1] = _h * twopi * this->_w_c * math_impl::tanh(this->_state.bp);
			tempstate[0] = this->_state.bp + deriv1[0] / 2;
			tempstate[1] = this->_state.lp + deriv1[1] / 2;
			
			np = input - 2 * _resInv * tempstate[0];
			hp = np - tempstate[1];
			deriv2[0] = _h * twopi * this->_w_c * math_impl::tanh(hp);
			deriv2[1] = _h * twopi * this->_w_c * math_impl::tanh(tempstate[0]);
			tempstate[0] = this->_state.bp + deriv2[0] / 2;
			tempstate[1] = this->_state.lp + deriv2[1] / 2;
			
			np = input - 2 * _resInv * tempstate[0];
			hp = np - tempstate[1];
			deriv3[0] = _h * twopi * this->_w_c * math_impl::tanh(hp);
			deriv3[1] = _h * twopi * this->_w_c * math_impl::tanh(tempstate[0]);
			tempstate[0] = this->_state.bp + deriv3[0];
			tempstate[1] = this->_state.lp + deriv3[1];
			
			np = input - 2 * _resInv * tempstate[0];
			hp = np - tempstate[1];
			deriv4[0] = _h * twopi * this->_w_c * math_impl::tanh(hp);
			deriv4[1] = _h * twopi * this->_w_c * math_impl::tanh(tempstate[0]);
			this->_state.bp += (1.f/6.f) * (deriv1[0] + 2 * deriv2[0] + 2 * deriv3[0] + deriv4[0]);
			this->_state.lp += (1.f/6.f) * (deriv1[1] + 2 * deriv2[1] + 2 * deriv3[1] + deriv4[1]);
			
			np = input - 2 * _resInv * this->_state.bp;
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
	};
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
	void setRise(float_t rise){
		this->rise = rise;
		this->riseInc = (this->fs_inv * 1000.f) / (this->rise);
		setRiseTarget(rise);
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
	void setFall(float_t fall){
		this->fall = fall;
		this->fallInc = (this->fs_inv * 1000.f) / (this->fall);
		setFallTarget(fall);
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

}   // namespace filters
}	// namespace nvs
