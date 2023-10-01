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
		-replace all instances of sin() and cos() with faster versions (lookup table if it remains stable)
	-free cutoff_to_g(), both because it CAN be free and because not all filters need it
	-remove #include <iostream>
**	-remove filter_abstract::z1
	-make lookup tables static
	-make lookup tables constexpr (really in nvs_memoryless)
 */

#pragma once
#include "nvs_memoryless.h"
#include "nvs_matrix2x2.h"
#include <iostream>
#include <array>

#include <concepts>
#include <type_traits>

namespace nvs	{
namespace filters {

template<typename float_t>
float_t cutoff_to_g_inaccurate(float_t cutoff, float_t fs_inv){
	return cutoff * fs_inv * 0.5f * (float_t)PI;
}
template<typename float_t>
float_t cutoff_to_g_slow(float_t cutoff, float_t fs_inv){
	return tan(cutoff * fs_inv * (float_t)PI);
}

template<typename float_t>
struct CutoffToG
{
	static_assert(std::is_floating_point<float_t>::value, "CutoffToG type must be floating point");

	void setSampleRate(float_t sr){
		_fs_inv = static_cast<float_t>(1) / sr;
	}
	float_t operator()(float_t cutoff){
		using namespace nvs::memoryless;
		if ((trig.tan_table) != NULL)
		{
			//float_t wc = TWOPI * cutoff;
			return static_cast<float_t>(trig.tan_LUT(cutoff * _fs_inv / 2.0));
		}
		else
			return 0.f;
	}
private:
	nvs::memoryless::trigTables<float_t> trig;
	float_t _fs_inv;
};

//==================================================================================

enum class mode_e {
	LP = 0,
	HP,
	BP,
	AP,
	NP
};

template<typename float_t>
class filter_abstract{
	static_assert(std::is_floating_point<float_t>::value, "filter_abstract type must be floating point");
public:
	virtual ~filter_abstract() = 0;
	//============================================================
	virtual void clear() = 0;
	virtual void setSampleRate(float_t sample_rate){
		_fs_inv = static_cast<float_t>(1.0) / sample_rate;
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
	nvs::memoryless::trigTables<float_t> trig;
	float_t _fs_inv;
	float_t _cutoffTarget, _resonanceTarget;
	float_t _w_c, _q;
	float_t _blockSize_inv;
	mode_e _mode {mode_e::LP};
};
template<typename float_t>
inline filter_abstract<float_t>::~filter_abstract() { }

/**
 TODO:
**	-get rid of y_n as stored value
	-don't use clamp if unnecessary (and it should be)
 */
template<
	typename float_t,
typename CutoffToGCallable=decltype([](float_t f_c, float_t fs)
				{ return cutoff_to_g_slow(f_c, fs); })
>
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
	
	float_t tpt_lp(float_t x, CutoffToGCallable func){
		auto const g = std::invoke(func, this->_w_c, this->_fs_inv);
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
		auto lp = tpt_lp(input, {});
		
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

template<
	typename float_t,
	typename CutoffToGCallable=decltype([](float_t f_c, float_t fs)
								{ return cutoff_to_g_slow(f_c, fs);})
>
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
	float_t tpt_lp(float_t x, CutoffToGCallable func){
		auto g = std::invoke(func, this->_w_c, this->_fs_inv);
		auto x_n = x * g;
		v_n = x_n + z1;
		auto y_n = tanh(v_n);
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
template<
	typename float_t,
	typename CutoffToGCallable=decltype([](float_t f_c, float_t fs)
							{ return cutoff_to_g_slow(f_c, fs);})
>
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
	float_t tpt_fourpole(float_t x, CutoffToGCallable func){
		float_t const g = std::invoke(func, this->_w_c, this->_fs_inv);
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
	std::array<onePole<float_t, CutoffToGCallable>, 4> _poles;
};

// so far quite tame, since the only nonlinearity is in the feedback path. 
/** TODO:
	-fix DC offset that arises after ANY Q has been introduced
 */
template<
	typename float_t,
	unsigned int N_iters=16,
	typename CutoffToGCallable=decltype([](float_t f_c, float_t fs) { return cutoff_to_g_slow(f_c, fs);}),
	class onePole_t=onePole_nonlinear_df2<float_t, CutoffToGCallable>
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
	
	float_t tpt_fourpole(float_t input, CutoffToGCallable func){
		float_t const g = std::invoke(func, this->_w_c, this->_fs_inv);
		float_t const g2 = g*g;
		float_t const G = g2 * g2;
		float_t const s1 = _poles[0].getState();
		float_t const s2 = _poles[1].getState();
		float_t const s3 = _poles[2].getState();
		float_t const s4 = _poles[3].getState();
		float_t const S = g2*g*s1 + g2*s2 + g*s3 + s4;
		
		float_t k = this->_q;
		if (this->_mode == mode_e::BP) k *= -1.f;

		u_n = y_n;  // initial estimation
		for (auto n = 0U; n < N_iters; n++){
			float_t tmp = nvs::memoryless::clamp(u_n, -100.f, 100.f);
			u_n = input - k * (G * tanh(tmp) + S);
		}
		y_n = _poles[3]
				(_poles[2]
					(_poles[1]
						(_poles[0]
							(u_n)
						 )
					 )
				 );
		
//		y_n = tanh(y_n);
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


/**
 TODO:
*	-replace cutoff_to_g_slow()
 */
template<typename float_t>
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
		const float_t omega = cutoff_to_g_slow(cutoff, this->_fs_inv);
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


template<typename float_t>
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
template<typename float_t>
class svf_lin_naive     :   public filter_abstract<float_t>, svf_prototype<float_t>
{
public:
	//==============================================================================
	svf_lin_naive()
	:   w_c(200.f), R(1.f), resonance(1.f)
	{
		this->sampleRate = 44100.f;
		this->fs_inv = 1.f / this->sampleRate;
		clear();
	}
	void clear(){
		this->_outputs = {0.f, 0.f, 0.f, 0.f };
		this->_state = { 0.f, 0.f };
	}
	void setCutoff(float_t wc){
		this->w_c = wc;
		this->_cutoffTarget = w_c;
	}
	void updateCutoff(float_t cutoff_target, float_t oneOverBlockSize){
		this->w_c += (cutoff_target - this->w_c) * oneOverBlockSize;
	}
	void updateCutoff(){
		this->w_c += (this->_cutoffTarget - this->w_c) * this->_oneOverBlockSize;
	}
	void setResonance(float_t res){
		this->resonance = res;
	 this->R = 1.f / res;
	}
	void updateResonance(float_t res_target, float_t oneOverBlockSize){
		if (res_target > 0.9f)
		res_target = 0.9f;

		this->resonance += (res_target - this->resonance) * oneOverBlockSize;
		this->R = 1.f / res_target;
	}
	void updateResonance(){
		if (this->_resonanceTarget > 0.9f)
		this->_resonanceTarget = 0.9f;

		this->resonance += (this->_resonanceTarget - this->resonance) * this->_oneOverBlockSize;
		this->R = 1.f / this->_resonanceTarget;
	}
	void filter(float_t input){
		float_t c, d;
		c = 2.f * sin(PI * w_c * this->fs_inv);
		d = 2.f * (1.f - pow(resonance, 0.25f));

		if (c > 0.5f) c = 0.5f;
		if (d > 2.0f) d = 2.f;
		if (d > (2.f/c - (c * 0.5f)))
		 d = 2.f/c - (c * 0.5f);

		this->_outputs.np = input - (d * this->_outputs.bp);
		this->_outputs.lp = this->_outputs.lp + (c * this->_outputs.bp);
		this->_outputs.hp = this->_outputs.np - this->_outputs.lp;
		this->_outputs.bp = this->_outputs.bp + (c * this->_outputs.hp);
	}
private:
	float_t w_c, R, resonance;
};
#if UNFINISHED_IMPLIMENTATIONS
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
template<typename float_t>
class svf_nl_rk :   public filter_abstract<float_t>, public svf_prototype<float_t>
{
public:
	svf_nl_rk();
	void setSampleRate(float_t sample_rate);
	void set_oversample(int oversample_factor);
	void clear();
	
	void setCutoff(float_t wc);
	void updateCutoff(float_t cutoff_target, float_t oneOverBlockSize);
	void updateCutoff();
	void setResonance(float_t res);
	void updateResonance(float_t res_target, float_t oneOverBlockSize);
	void updateResonance();
	void filter(float_t input);
private:
	int _oversample_factor;
	float_t h, w_c, R, resonance;
	// k_1 through k_4. for each, [0] is bp, [1] is lp.
	float_t deriv1[2], deriv2[2], deriv3[2], deriv4[2];
};

//==============================================================================
template<typename float_t>
class slewlim
{
public:
	slewlim();
	slewlim(float_t sample_rate);
	~slewlim()  { }
	//============================================================
	void setSampleRate(float_t sample_rate);
	//============================================================
	// immediate change
	void setRise(float_t rise);
	// change over block size
	void setRise();
	void setRise(float_t riseTarget, float_t oneOverBlockSize);
	// immediate change
	void setFall(float_t fall);
	// change over block size
	void setFall();
	void setFall(float_t fallTarget, float_t oneOverBlockSize);
	//============================================================
	float_t ASR(float_t gate);
	
	float_t _riseTarget, _fallTarget;
	float_t _oneOverBlockSize;
private:
	float_t sampleRate, fs_inv;
	
	// 'Inc' variables tell change per sample.
	float_t rise, riseInc, fall, fallInc, _vOut;
};

template<typename float_t>
class dcBlock   :   public filter_abstract <float_t>
{
public:
	dcBlock();
	dcBlock(float_t sampleRate);
	void setSampleRate(float_t sampleRate);
	void clear();
	void updateCutoff(float_t cutoff_target, float_t oneOverBlockSize){ }
	void updateResonance(float_t res_target, float_t oneOverBlockSize){ }
	void updateR(float_t R_target, float_t oneOverBlockSize);
	float_t filter(float_t x);
private:
	float_t R, xz1, yz1;
	float_t sample_rate;
};

//===============================================================================
// PIRKLE IMPLEMENTATIONS (not my own work; used only for checking.)
template<typename float_t>
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
		float_t wd = 2 * PI * fc;
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
template <typename float_t>
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
		float_t const wd = 2 * PI * fc;
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
//==================================================================================
/* 
 time-variant allpass filter
 */
template<typename float_t>
class tvap  :   public filter_abstract<float_t>
{
public:
	tvap() {}
	~tvap() {}
	//==============================================================================
	virtual void setSampleRate(float_t sample_rate);
	
	void clear();
	void updateCutoff(float_t cutoff_target, float_t oneOverBlockSize);
	void updateResonance(float_t res_target, float_t oneOverBlockSize);
	// function aliases just to have more meaningful names
	void update_f_pi (float_t f_pi_target, float_t oneOverBlockSize);
	void update_f_b(float_t f_b_target, float_t oneOverBlockSize);
	
	void calc_b1(void);
	void f_b_to_b0(void);
	
	float_t f_pi2r2(float_t _f_pi);
	float_t f_b2r1(float_t _f_b);
	
	float_t filter(float_t x_n);
	
	// should be in memoryless but got linker error
	inline float_t unboundSat2(float_t x);
	
	float_t filter_fbmod(float_t x_n, float_t fb_f_pi, float_t fb_f_b);
	nvs::memoryless::trigTables<float_t> _LUT;
	
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
#endif	// end of unfinished re-implimentations
}   // namespace filters
}	// namespace nvs
