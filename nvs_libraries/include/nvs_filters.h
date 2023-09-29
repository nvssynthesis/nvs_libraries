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

namespace nvs	{
namespace filters {

template<typename float_t>
float_t cutoff_to_g_inaccurate(float_t cutoff, float_t fs_inv){
	return cutoff * fs_inv * 0.5;
}
template<typename float_t>
float_t cutoff_to_g_slow(float_t cutoff, float_t fs_inv){
	return tan(cutoff * fs_inv * PI);
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
	virtual void updateCutoff(){
		_w_c += (_cutoffTarget - _w_c) * _blockSize_inv;
	}
	virtual void setResonanceTarget(float_t res_target){
		_resonanceTarget = res_target;
	}
	virtual void updateResonance(){
		_q += (_resonanceTarget - _q) * _blockSize_inv;
	}
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
};
template<typename float_t>
inline filter_abstract<float_t>::~filter_abstract() { }

/**
 TODO:
**	-get rid of y_n as stored value
	-don't use clamp if unnecessary (and it should be)
 */
template<typename float_t>
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
	
	float_t tpt_lp(float_t input){
		g = cutoff_to_g_slow(this->_w_c, this->_fs_inv);
		v_n = (input - this->z1) * g / (1.0 + g);
		y_n = v_n + this->z1;
		this->z1 = y_n + v_n;
		return nvs::memoryless::clamp<float_t>(y_n, -1000.0, 1000.0);
	 }
	float_t tpt_lp(float_t input, float_t cutoff){
		g = cutoff_to_g_slow(cutoff, this->_fs_inv);
		v_n = (input - this->z1) * g / (1.0 + g);
		y_n = v_n + this->z1;
		this->z1 = y_n + v_n;
		return nvs::memoryless::clamp<float_t>(y_n, -1000.0, 1000.0);
	}
	
	float_t tpt_hp(float_t input){
		return input - tpt_lp(input);
	}
	float_t tpt_hp(float_t input, float_t cutoff){
		return input - tpt_lp(input, cutoff);
	}
	//==============================================================================
	float_t operator()(float_t input) override {
		return tpt_lp(input);
	}
	float_t operator()(float_t input, float_t cutoff) override {
		return tpt_lp(input, cutoff);
	}
	float_t operator()(float_t input, float_t cutoff, float_t) override {
		/*
		 is there a way to warn upon calling a function? but i do not want to deprecate this
		 because it could have some use for generic interfaces
		 */
		return tpt_lp(input, cutoff);
	}
private:
	float_t v_n {0.0}, y_n {0.0}, z1 {0.0};
	float_t g {0.0};
};
//==============================================================================

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
	
	void calcCoefs(const float_t cutoff){
		const float_t omega = cutoff_to_g_slow(cutoff, this->_fs_inv);
		const float_t omega2 = omega * omega;
		const float_t twoCosPiOver4xOmega = 1.414213562373095 * omega;
		const float_t c = 1.0 + twoCosPiOver4xOmega + omega2;
		
		const float_t b0 = omega2 / c;
		const float_t b1 = 2.0 * b0;
		const float_t a1 = (2.0 * (omega2 - 1.0)) / c;
		const float_t a2 = (1.0 - twoCosPiOver4xOmega + omega2) / c;
		
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

// NOTHING SO FAR.
template<typename float_t>
class onePole_nonlinear   :   public onePole<float_t>
{
public:
private:
};

template<typename float_t>
class fourPole_LP_linear    :   public filter_abstract<float_t>
{
public:
	fourPole_LP_linear();
	fourPole_LP_linear(float_t sample_rate);
	void initialize(float_t sample_rate);
	void updateOneOverBlockSize(float_t oneOverBlockSize);
	void clear();
	void updateCutoff();
	void updateCutoff(float_t cutoff_target, float_t oneOverBlockSize);
	void updateResonance();
	void updateResonance(float_t res_target, float_t oneOverBlockSize);
	
	float_t tpt_fourpole(float_t input);
	float_t tpt_fourpole(float_t input, float_t cutoff);
	
	float_t _resonanceTarget;
private:
	onePole<float_t> H1, H2, H3, H4;
	float_t u_n, y_n, s1, s2, s3, s4, S, y1, y2, y3, y4, g, g_denom, G, k;
	float_t w_c, q;
};

// so far quite tame, since the only nonlinearity is in the feedback path. 
// TODO: convert each onepole into a nonlinear onepole
template<typename float_t>
class fourPole_LP_nonlinear    :   public filter_abstract<float_t>
{
public:
	fourPole_LP_nonlinear();
	fourPole_LP_nonlinear(float_t sample_rate);
	void initialize(float_t sample_rate);
	void updateOneOverBlockSize(float_t oneOverBlockSize);
	void clear();
	void updateCutoff();
	void updateCutoff(float_t cutoff_target, float_t oneOverBlockSize);
	
	void updateResonance();
	void updateResonance(float_t res_target, float_t oneOverBlockSize);
	
	
	float_t tpt_fourpole(float_t input);
	float_t tpt_fourpole(float_t input, float_t cutoff);
	//    float_t _resonanceTarget;
private:
	onePole<float_t> H1, H2, H3, H4;
	nvs::memoryless::trigTables<float_t> tables;
	int iters;
	float_t u_n, y_n, s1, s2, s3, s4, S, y1, y2, y3, y4, g, g_denom, G, k;
	float_t w_c, q;
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
	svf_lin_naive();
	void clear();
	void setCutoff(float_t wc);
	void updateCutoff(float_t cutoff_target, float_t oneOverBlockSize);
	void updateCutoff();
	void setResonance(float_t res);
	void updateResonance(float_t res_target, float_t oneOverBlockSize);
	void updateResonance();
	void filter(float_t input);
private:
	
	float_t w_c, R, resonance;
	
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
	inline void initialize(float_t newSampleRate);
	void setFc(float_t fc);
	
	float_t doFilterStage(float_t xn);
	float_t getSampleRate();
	float_t getStorageRegisterValue();
};
template <typename float_t>
class CTPTMoogLadderFilter
{
public:
	CTPTMoogLadderFilter(){}
	~CTPTMoogLadderFilter(){}
protected:
	CTPTMoogFilterStage<float_t> filter1;
	CTPTMoogFilterStage<float_t> filter2;
	CTPTMoogFilterStage<float_t> filter3;
	CTPTMoogFilterStage<float_t> filter4;
	float_t k; // Q control
	float_t fc; // fc control
public:
	inline void initialize(float_t newSampleRate);
	inline void calculateTPTCoeffs(float_t cutoff, float_t Q);
	float_t doTPTMoogLPF(float_t xn);
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

}   // namespace filters
}	// namespace nvs
