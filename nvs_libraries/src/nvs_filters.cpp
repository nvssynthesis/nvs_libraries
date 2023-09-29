/*
  ==============================================================================
    nvs_filters.cpp
    Created: 6 Dec 2018 12:12:44am
    Author:  Nicholas Solem
  ==============================================================================
*/
//#pragma once
#ifndef NVS_FILTERS_H
#define NVS_FILTERS_H
#include "../include/nvs_filters.h"

/**
 TODO:
	-optimize tvap:
		-replace all instances of sin() and cos() with faster versions (lookup table if it remains stable)
 */


namespace nvs	{
namespace filters {

//==============================================================================
template<typename float_t>
fourPole_LP_linear<float_t>::fourPole_LP_linear()
:   u_n(0.f), s1(0.f), s2(0.f), s3(0.f), s4(0.f), S(0.f), y1(0.f), y2(0.f), y3(0.f), y4(0.f), G(0.f), k(0.f)
{
	float_t def_sampRate = (float_t)44100;
	this->setSampleRate(def_sampRate);
	H1.setSampleRate(def_sampRate);
	H2.setSampleRate(def_sampRate);
	H3.setSampleRate(def_sampRate);
	H4.setSampleRate(def_sampRate);
}
template<typename float_t>
fourPole_LP_linear<float_t>::fourPole_LP_linear(float_t sample_rate)
:   u_n(0.f), s1(0.f), s2(0.f), s3(0.f), s4(0.f), S(0.f), y1(0.f), y2(0.f), y3(0.f), y4(0.f), G(0.f), k(0.f)
{       
	this->setSampleRate(sample_rate);
	H1.setSampleRate(sample_rate);
	H2.setSampleRate(sample_rate);
	H3.setSampleRate(sample_rate);
	H4.setSampleRate(sample_rate);
}
template<typename float_t>
void fourPole_LP_linear<float_t>::initialize(float_t sample_rate)
{
	this->clear();
	this->setSampleRate(sample_rate);
	H1.setSampleRate(sample_rate);
	H2.setSampleRate(sample_rate);
	H3.setSampleRate(sample_rate);
	H4.setSampleRate(sample_rate);
}
template<typename float_t>
void fourPole_LP_linear<float_t>::updateOneOverBlockSize(float_t oneOverBlockSize)
{
	this->_oneOverBlockSize = oneOverBlockSize;
	H1._oneOverBlockSize = oneOverBlockSize;
	H2._oneOverBlockSize = oneOverBlockSize;
	H3._oneOverBlockSize = oneOverBlockSize;
	H4._oneOverBlockSize = oneOverBlockSize;
}
template<typename float_t>
void fourPole_LP_linear<float_t>::clear()
{
	H1.clear();
	H2.clear();
	H3.clear();
	H4.clear();
	s1 = s2 = s3 = s4 = u_n = y_n = S = 0.f;
}
template<typename float_t>
void fourPole_LP_linear<float_t>::updateCutoff()
{
	float_t local_cutoff = this->w_c;
	if (this->_cutoffTarget != local_cutoff)
	{
		this->w_c += (this->_cutoffTarget - local_cutoff) * this->_oneOverBlockSize;
		H1._cutoffTarget = local_cutoff;
		H2._cutoffTarget = local_cutoff;
		H3._cutoffTarget = local_cutoff;
		H4._cutoffTarget = local_cutoff;
		H1.updateCutoff();
		H2.updateCutoff();
		H3.updateCutoff();
		H4.updateCutoff();
	}
}
template<typename float_t>
void fourPole_LP_linear<float_t>::updateCutoff(float_t cutoff_target, float_t oneOverBlockSize)
{
	this->w_c += (cutoff_target - this->w_c) * oneOverBlockSize;
	H1.updateCutoff(cutoff_target, oneOverBlockSize);
	H2.updateCutoff(cutoff_target, oneOverBlockSize);
	H3.updateCutoff(cutoff_target, oneOverBlockSize);
	H4.updateCutoff(cutoff_target, oneOverBlockSize);
}
template<typename float_t>
void fourPole_LP_linear<float_t>::updateResonance()
{
	if (_resonanceTarget != k)
		this->k += (this->_resonanceTarget - this->k) * this->_oneOverBlockSize;
}
template<typename float_t>
void fourPole_LP_linear<float_t>::updateResonance(float_t res_target, float_t oneOverBlockSize)
{
	//this->q += (res_target - this->q) * oneOverBlockSize;
	this->k += (res_target - this->k) * oneOverBlockSize;
}
template<typename float_t>
float_t fourPole_LP_linear<float_t>::tpt_fourpole(float_t input)
{
	g = this->cutoff_to_g(this->w_c);
	float_t g2 = g*g;
	G = g2 * g2;
	s1 = H1.z1;
	s2 = H2.z1;
	s3 = H3.z1;
	s4 = H4.z1;
	S = g2*g*s1 + g2*s2 + g*s3 + s4;
	
	u_n = ((input) * (1 + k) - k * S) / (1 + k * G);
	
	y_n = H4.tpt_lp(H3.tpt_lp(H2.tpt_lp(H1.tpt_lp(u_n))));
	
	return y_n;
}
template<typename float_t>
float_t fourPole_LP_linear<float_t>::tpt_fourpole(float_t input, float_t cutoff)
{
	g = this->cutoff_to_g(cutoff);
	G = g * g * g * g;
	s1 = H1.z1;
	s2 = H2.z1;
	s3 = H3.z1;
	s4 = H4.z1;
	float_t g2 = g*g;
	S = g2*g*s1 + g2*s2 + g*s3 + s4;
	
	u_n = ((input) * (1 + k) - k * S) / (1 + k * G);
	
	y_n = H4.tpt_lp(H3.tpt_lp(H2.tpt_lp(H1.tpt_lp(u_n))));
	
	return y_n;
}

// so far quite tame, since the only nonlinearity is in the feedback path. 
// TODO: convert each onepole into a nonlinear onepole
template<typename float_t>
fourPole_LP_nonlinear<float_t>::fourPole_LP_nonlinear()
:   iters(16),
u_n(0.f), s1(0.f), s2(0.f), s3(0.f), s4(0.f), S(0.f), y1(0.f), y2(0.f), y3(0.f), y4(0.f), G(0.f), k(0.f)
{
	float_t def_sampRate = (float_t)44100;
	this->setSampleRate(def_sampRate);
	H1.setSampleRate(def_sampRate);
	H2.setSampleRate(def_sampRate);
	H3.setSampleRate(def_sampRate);
	H4.setSampleRate(def_sampRate);
}
template<typename float_t>
fourPole_LP_nonlinear<float_t>::fourPole_LP_nonlinear(float_t sample_rate)
:   iters(16),
u_n(0.f), s1(0.f), s2(0.f), s3(0.f), s4(0.f), S(0.f), y1(0.f), y2(0.f), y3(0.f), y4(0.f), G(0.f), k(0.f)
{
	this->setSampleRate(sample_rate);
	H1.setSampleRate(sample_rate);
	H2.setSampleRate(sample_rate);
	H3.setSampleRate(sample_rate);
	H4.setSampleRate(sample_rate);
}
template<typename float_t>
void fourPole_LP_nonlinear<float_t>::initialize(float_t sample_rate)
{
	this->clear();
	this->setSampleRate(sample_rate);
	H1.setSampleRate(sample_rate);
	H2.setSampleRate(sample_rate);
	H3.setSampleRate(sample_rate);
	H4.setSampleRate(sample_rate);
}
template<typename float_t>
void fourPole_LP_nonlinear<float_t>::updateOneOverBlockSize(float_t oneOverBlockSize)
{
	this->_oneOverBlockSize = oneOverBlockSize;
	H1._oneOverBlockSize = oneOverBlockSize;
	H2._oneOverBlockSize = oneOverBlockSize;
	H3._oneOverBlockSize = oneOverBlockSize;
	H4._oneOverBlockSize = oneOverBlockSize;
}
template<typename float_t>
void fourPole_LP_nonlinear<float_t>::clear()
{
	H1.clear();
	H2.clear();
	H3.clear();
	H4.clear();
	s1 = s2 = s3 = s4 = u_n = y_n = S = 0.f;
}
template<typename float_t>
void fourPole_LP_nonlinear<float_t>::updateCutoff()
{
	float_t local_cutoff = this->w_c;
	if (this->_cutoffTarget != local_cutoff)
	{
		this->w_c += (this->_cutoffTarget - local_cutoff) * this->_oneOverBlockSize;
		H1._cutoffTarget = local_cutoff;
		H2._cutoffTarget = local_cutoff;
		H3._cutoffTarget = local_cutoff;
		H4._cutoffTarget = local_cutoff;
		H1.updateCutoff();
		H2.updateCutoff();
		H3.updateCutoff();
		H4.updateCutoff();
	}
}
template<typename float_t>
void fourPole_LP_nonlinear<float_t>::updateCutoff(float_t cutoff_target, float_t oneOverBlockSize)
{
	this->w_c += (cutoff_target - this->w_c) * oneOverBlockSize;
	H1.updateCutoff(cutoff_target, oneOverBlockSize);
	H2.updateCutoff(cutoff_target, oneOverBlockSize);
	H3.updateCutoff(cutoff_target, oneOverBlockSize);
	H4.updateCutoff(cutoff_target, oneOverBlockSize);
}
template<typename float_t>
void fourPole_LP_nonlinear<float_t>::updateResonance()
{
	if (this->_resonanceTarget != k)
		this->k += (this->_resonanceTarget - this->k) * this->_oneOverBlockSize;
}
template<typename float_t>
void fourPole_LP_nonlinear<float_t>::updateResonance(float_t res_target, float_t oneOverBlockSize)
{
	//this->q += (res_target - this->q) * oneOverBlockSize;
	this->k += (res_target - this->k) * oneOverBlockSize;
}
template<typename float_t>
float_t fourPole_LP_nonlinear<float_t>::tpt_fourpole(float_t input)
{
	g = this->cutoff_to_g(this->w_c);
	float_t g2 = g * g;
	G = g2 * g2;
	s1 = H1.z1;
	s2 = H2.z1;
	s3 = H3.z1;
	s4 = H4.z1;
	S = g2*g*s1 + g2*s2 + g*s3 + s4;
	
	u_n = y_n;  // initial estimation
	for (int n = 0; n < iters; n++)
		u_n = input - k * (G * tables.tanh_LUT(u_n) + S);
	
	//u_n = ((input) * (1 + k) - k * S) / (1 + k * G);
	
	y_n = H4.tpt_lp(H3.tpt_lp(H2.tpt_lp(H1.tpt_lp(u_n))));
	
	return y_n;
}
template<typename float_t>
float_t fourPole_LP_nonlinear<float_t>::tpt_fourpole(float_t input, float_t cutoff)
{
	g = this->cutoff_to_g(cutoff);
	float_t g2 = g * g;
	G = g2 * g2;
	s1 = H1.z1;
	s2 = H2.z1;
	s3 = H3.z1;
	s4 = H4.z1;
	S = g2*g*s1 + g2*s2 + g*s3 + s4;
	
	u_n = ((input) * (1 + k) - k * S) / (1 + k * G);
	
	y_n = H4.tpt_lp(H3.tpt_lp(H2.tpt_lp(H1.tpt_lp(u_n))));
	
	return y_n;
}

// linear state variable filter using 'naive' integrators (i.e., Euler backward difference integration)
//==============================================================================
template<typename float_t>
svf_lin_naive<float_t>::svf_lin_naive()
:   w_c(200.f), R(1.f), resonance(1.f)
{
	this->sampleRate = 44100.f;
	this->fs_inv = 1.f / this->sampleRate;
	clear();
}
template<typename float_t>
void svf_lin_naive<float_t>::clear()
{
	this->_outputs = {0.f, 0.f, 0.f, 0.f };
	this->_state = { 0.f, 0.f };
}
template<typename float_t>
void svf_lin_naive<float_t>::setCutoff(float_t wc)
{
	this->w_c = wc;
	this->_cutoffTarget = w_c;
}
template<typename float_t>
void svf_lin_naive<float_t>::updateCutoff(float_t cutoff_target, float_t oneOverBlockSize)
{
	this->w_c += (cutoff_target - this->w_c) * oneOverBlockSize;
}
template<typename float_t>
void svf_lin_naive<float_t>::updateCutoff()
{
	this->w_c += (this->_cutoffTarget - this->w_c) * this->_oneOverBlockSize;
}
template<typename float_t>
void svf_lin_naive<float_t>::setResonance(float_t res)
{
	this->resonance = res;
	this->R = 1.f / res;
}
template<typename float_t>
void svf_lin_naive<float_t>::updateResonance(float_t res_target, float_t oneOverBlockSize)
{
	if (res_target > 0.9f)
		res_target = 0.9f;
	
	this->resonance += (res_target - this->resonance) * oneOverBlockSize;
	this->R = 1.f / res_target;
}
template<typename float_t>
void svf_lin_naive<float_t>::updateResonance()
{
	if (this->_resonanceTarget > 0.9f)
		this->_resonanceTarget = 0.9f;
	
	this->resonance += (this->_resonanceTarget - this->resonance) * this->_oneOverBlockSize;
	this->R = 1.f / this->_resonanceTarget;
}
/*
 void filter(float_t input)
 {
 using namespace nvs_memoryless;
 _outputs.lp = _state.lp + this->w_c * this->fs_inv * _state.bp;
 _outputs.hp = clamp<float_t>((input - (2 * R * _state.bp) - _state.lp), -10.f, 10.f);
 _outputs.bp = _state.bp + this->w_c * this->fs_inv * _outputs.hp;
 
 // update state
 _state.bp = clamp<float_t>(_outputs.bp, -10.f, 10.f);
 _state.lp = clamp<float_t>(_outputs.lp, -10.f, 10.f);
 }
 */
template<typename float_t>
void svf_lin_naive<float_t>::filter(float_t input)
{
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
svf_nl_rk<float_t>::svf_nl_rk()
: _oversample_factor(1), h(0.000022675736961),
w_c(200.0), R(1.0), resonance(1.0)
{
	this->sampleRate = (float_t)44100;
	for (int i = 0; i < 2; i++)
	{
		deriv1[i] = deriv2[i] = deriv3[i] = deriv4[i] = 0.f;
	}
}
template<typename float_t>
void svf_nl_rk<float_t>::setSampleRate(float_t sample_rate)
{
	this->sampleRate = sample_rate;
	this->fs_inv = 1.f / sample_rate;
	this->h = 1.f / (_oversample_factor * sample_rate);
}  
template<typename float_t>
void svf_nl_rk<float_t>::set_oversample(int oversample_factor)
{
	_oversample_factor = oversample_factor;
	h = 1.f / (oversample_factor * this->sampleRate);
}
template<typename float_t>
void svf_nl_rk<float_t>::clear()
{
	for (int i = 0; i < 2; i++)
	{
		deriv1[i] = deriv2[i] = deriv3[i] = deriv4[i] = 0.f;
		this->_outputs = {0.f, 0.f, 0.f, 0.f };
		this->_state = { 0.f, 0.f };
	}
}
template<typename float_t>
void svf_nl_rk<float_t>::setCutoff(float_t wc)
{
	this->w_c = wc;
}
template<typename float_t>
void svf_nl_rk<float_t>::updateCutoff(float_t cutoff_target, float_t oneOverBlockSize)
{
	this->w_c += (cutoff_target - this->w_c) * oneOverBlockSize;
}
template<typename float_t>
void svf_nl_rk<float_t>::updateCutoff()
{
	this->w_c += (this->_cutoffTarget - this->w_c) * this->_oneOverBlockSize;
}
template<typename float_t>
void svf_nl_rk<float_t>::setResonance(float_t res)
{
	if (res < 0.5f) res = 0.5f;
	this->resonance = res;
	this->R = 1.f / res;
}
template<typename float_t>
void svf_nl_rk<float_t>::updateResonance(float_t res_target, float_t oneOverBlockSize)
{
	if (res_target < 0.5f) res_target = 0.5f;
	this->resonance += (res_target - this->resonance) * oneOverBlockSize;
	this->R = 1.f / this->resonance;
}
template<typename float_t>
void svf_nl_rk<float_t>::updateResonance()
{
	if (this->_resonanceTarget < 0.5f) this->_resonanceTarget = 0.5f;
	this->resonance += (this->_resonanceTarget - this->resonance) * this->_oneOverBlockSize;
	this->R = 1.f / this->resonance;
}
template<typename float_t>
void svf_nl_rk<float_t>::filter(float_t input)
{
	using namespace nvs::memoryless;
	float_t hp(0), np(0);
	// overwritten states. [0] is bp, [1] is lp.
	float_t tempstate[2];
	
	for (int iter = 0; iter < _oversample_factor; iter++)
	{
		np = input - 2 * R * this->_state.bp;
		hp = np - this->_state.lp;
		deriv1[0] = h * TWOPI * w_c * this->trig.tanh_LUT(hp);
		deriv1[1] = h * TWOPI * w_c * this->trig.tanh_LUT(this->_state.lp);
		
		// I THINK THE PROBLEM IS WITH SCALING, ORDER OF OPERATIONS REGARDING H
		tempstate[0] = this->_state.bp + deriv1[0] / 2;
		tempstate[1] = this->_state.lp + deriv1[1] / 2;
		
		np = input - 2 * R * this->_state.bp;
		hp = np - this->_state.lp;
		deriv2[0] = h * TWOPI * w_c * this->trig.tanh_LUT(hp); // is this the right move?
		deriv2[1] = h * TWOPI * w_c * this->trig.tanh_LUT(tempstate[0]);
		
		tempstate[0] = this->_state.bp + deriv2[0] / 2;
		tempstate[1] = this->_state.lp + deriv2[1] / 2;
		
		np = input - 2 * R * this->_state.bp;
		hp = np - this->_state.lp;
		deriv3[0] = h * TWOPI * w_c * this->trig.tanh_LUT(hp);
		deriv3[1] = h * TWOPI * w_c * this->trig.tanh_LUT(tempstate[0]);
		
		tempstate[0] = this->_state.bp + deriv3[0];
		tempstate[1] = this->_state.lp + deriv3[1];
		
		np = input - 2 * R * this->_state.bp;
		hp = np - this->_state.lp;
		deriv4[0] = h * TWOPI * w_c * this->trig.tanh_LUT(hp);
		deriv4[1] = h * TWOPI * w_c * this->trig.tanh_LUT(tempstate[0]);
		
		this->_state.bp += (1.f/6.f) * (deriv1[0] + 2 * deriv2[0] + 2 * deriv3[0] + deriv4[0]);
		this->_state.lp += (1.f/6.f) * (deriv1[1] + 2 * deriv2[1] + 2 * deriv3[1] + deriv4[1]);
	}
	
	this->_outputs.bp = this->_state.bp;
	this->_outputs.lp = this->_state.lp;
	this->_outputs.hp = hp;
	this->_outputs.np = np;
}

//==============================================================================
template<typename float_t>
slewlim<float_t>::slewlim()
: _vOut(0.f)
{
	setSampleRate(44100.f);
}
template<typename float_t>
slewlim<float_t>::slewlim(float_t sample_rate)
: _vOut(0.f)
{
	setSampleRate(sample_rate);
}
//============================================================
template<typename float_t>
void slewlim<float_t>::setSampleRate(float_t sample_rate)
{
	this->sampleRate = sample_rate;
	this->fs_inv = 1.f / sample_rate;
}
//============================================================
// immediate change
template<typename float_t>
void slewlim<float_t>::setRise(float_t rise)
{
	this->rise = rise;
	this->riseInc = (this->fs_inv * 1000.f) / (this->rise);
}
// change over block size
template<typename float_t>
void slewlim<float_t>::setRise()
{
	if (_riseTarget != rise)
	{
		this->rise += (_riseTarget - rise) * _oneOverBlockSize;
		this->riseInc = (this->fs_inv * 1000.f) / (rise);
	}
}
template<typename float_t>
void slewlim<float_t>::setRise(float_t riseTarget, float_t oneOverBlockSize)
{
	this->rise += (riseTarget - this->rise) * oneOverBlockSize;
	this->riseInc = (this->fs_inv * 1000.f) / (this->rise);
}
// immediate change
template<typename float_t>
void slewlim<float_t>::setFall(float_t fall)
{
	this->fall = fall;
	this->fallInc = (this->fs_inv * 1000.f) / (this->fall);
}
// change over block size
template<typename float_t>
void slewlim<float_t>::setFall()
{
	if (_fallTarget != fall)
	{
		this->fall += (_fallTarget - fall) * _oneOverBlockSize;
		this->fallInc = (this->fs_inv * 1000.f) / (fall);
	}
}
template<typename float_t>
void slewlim<float_t>::setFall(float_t fallTarget, float_t oneOverBlockSize)
{
	this->fall += (fallTarget - this->fall) * oneOverBlockSize;
	this->fallInc = (this->fs_inv * 1000.f) / (this->fall);
}
//============================================================
template<typename float_t>
float_t slewlim<float_t>::ASR(float_t gate)
{
	using namespace nvs::memoryless;
	
	setRise();
	setFall();
	
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
//=============================================================================
template<typename float_t>
dcBlock<float_t>::dcBlock()   :   dcBlock((float_t)44100.0)
{}
template<typename float_t>
dcBlock<float_t>::dcBlock(float_t sample_rate)
{
	this->setSampleRate((float_t)sample_rate);
	clear();
	R = 0.995f;  //fixed R, should base on sample rate instead
}
template<typename float_t>
void dcBlock<float_t>::setSampleRate(float_t sampleRate){
	this->sample_rate = sampleRate;
	this->fs_inv = 1.f / sampleRate;
}
template<typename float_t>
void dcBlock<float_t>::clear()
{
	xz1 = yz1 = 0.f;
}
template<typename float_t>
void dcBlock<float_t>::updateR(float_t R_target, float_t oneOverBlockSize)
{
	R += (R_target - R) * oneOverBlockSize;
}
template<typename float_t>
float_t dcBlock<float_t>::filter(float_t x)
{// y[i] = x[i] - x[i - 1] + R * y[i - 1];
	yz1 = x - xz1 + R * yz1;
	xz1 = x;
	return yz1;
}
//===============================================================================
// PIRKLE IMPLEMENTATIONS (not my own work; used only for checking.)
template<typename float_t>
inline void CTPTMoogFilterStage<float_t>::initialize(float_t newSampleRate)
{
	// save
	sampleRate = newSampleRate;
	z1 = 0;
}
template<typename float_t>
void CTPTMoogFilterStage<float_t>::setFc(float_t fc)
{
	// prewarp the cutoff- these are bilinear-transform filters
	float_t wd = 2 * PI * fc;
	float_t fs_inv  = 1 / sampleRate;
	float_t wa = (2 / fs_inv) * tan(wd * fs_inv / 2);
	float_t g  = wa * fs_inv / 2;
	// calculate big G value; see Zavalishin p46 The Art of VA Design
	G = g / (1.0 + g);
}
template<typename float_t>
float_t CTPTMoogFilterStage<float_t>::doFilterStage(float_t xn)
{
	float_t v = (xn - z1) * G;
	float_t out = v + z1;
	z1 = out + v;
	return out;
}
template<typename float_t>
float_t CTPTMoogFilterStage<float_t>::getSampleRate()
{
	return sampleRate;
}
template<typename float_t>
float_t CTPTMoogFilterStage<float_t>::getStorageRegisterValue()
{
	return z1;
}
template<typename float_t>
inline void CTPTMoogLadderFilter<float_t>::initialize(float_t newSampleRate)
{
	filter1.initialize(newSampleRate);
	filter2.initialize(newSampleRate);
	filter3.initialize(newSampleRate);
	filter4.initialize(newSampleRate);
}
template<typename float_t>
inline void CTPTMoogLadderFilter<float_t>::calculateTPTCoeffs(float_t cutoff, float_t Q)
{
	// 4 sync-tuned filters
	filter1.setFc(cutoff);
	filter2.setFc(cutoff);
	filter3.setFc(cutoff);
	filter4.setFc(cutoff);
	
	// NOTE: Q is limited to 20 on the UI to prevent blowing up //Q=0.707->25==>k=0->4
	k = 4.0*(Q - 0.707)/(25.0 - 0.707);
	// ours
	fc = cutoff;
}
template<typename float_t>
float_t CTPTMoogLadderFilter<float_t>::doTPTMoogLPF(float_t xn)
{
	// calculate g
	float_t wd = 2 * PI * fc;
	float_t fs_inv  = 1 / (float_t)filter1.getSampleRate();
	float_t wa = (2 / fs_inv) * tan(wd * fs_inv / 2);
	float_t g  = wa * fs_inv / 2;
	float_t G = g * g * g * g;
	float_t S = g * g * g * filter1.getStorageRegisterValue() +
	g * g * filter2.getStorageRegisterValue() +
	g * filter3.getStorageRegisterValue() +
	filter4.getStorageRegisterValue();
	//uis input to filters, straight from book
	float_t u = (xn - k * S) / (1 + k * G);
	// four cascades using nested functions
	float_t filterOut = filter4.doFilterStage(filter3.doFilterStage (filter2.doFilterStage(filter1.doFilterStage(u))));
	// output
	return filterOut;
}
//==================================================================================
/* 
 time-variant allpass filter
 */
template<typename float_t>
void tvap<float_t>::setSampleRate(float_t sample_rate)
{
	this->sampleRate = sample_rate;
	this->fs_inv = 1.f / sample_rate;
	lp.setSampleRate(sample_rate);
	lp.updateCutoff(sample_rate * 0.125f, 1.f);
}    
template<typename float_t>
void tvap<float_t>::clear(){
	sp->z1 = sp->z2 = sp->fb_proc = 0.f;
};
template<typename float_t>
void tvap<float_t>::updateCutoff(float_t cutoff_target, float_t oneOverBlockSize) {
	if (cutoff_target < 0) cutoff_target = 0;
	f_pi += (cutoff_target - f_pi) * oneOverBlockSize;
	calc_b1();
}
template<typename float_t>
void tvap<float_t>::updateResonance(float_t res_target, float_t oneOverBlockSize){
	//if (f_b <= 0) f_b = 0.0000000001;
	f_b += (res_target - f_b) * oneOverBlockSize;
	f_b_to_b0();
}
// function aliases just to have more meaningful names
template<typename float_t>
void tvap<float_t>::update_f_pi (float_t f_pi_target, float_t oneOverBlockSize)
{   updateCutoff(f_pi_target, oneOverBlockSize); }
template<typename float_t>
void tvap<float_t>::update_f_b(float_t f_b_target, float_t oneOverBlockSize)
{   updateResonance(f_b_target, oneOverBlockSize); }

template<typename float_t>
void tvap<float_t>::calc_b1(void) {
	float_t d = -1 * cos((2.f * PI * f_pi) / this->sampleRate);
	float_t c = (tan(PI * f_b / this->sampleRate) - 1.f) / (tan(PI * f_b / this->sampleRate) + 1.f);
	float_t r1 = acos(-1.f * c);
	float_t r2 = acos(-1.f * d);
	b1 = cos(r1) * (1.f + cos(r2));
}
template<typename float_t>
void tvap<float_t>::f_b_to_b0(void) {
	float_t c = (tan(PI * f_b / this->sampleRate) - 1.f) / (tan(PI * f_b / this->sampleRate) + 1.f);
	float_t r1 = acos(-1.f * c);
	b0 = cos(r1);
}

template<typename float_t>
float_t tvap<float_t>::f_pi2r2(float_t _f_pi)
{
	float_t d = -1 * cos((2.f * PI * _f_pi) / this->sampleRate);
	float_t r2 = acos(-d);
	return r2;
}
template<typename float_t>
float_t tvap<float_t>::f_b2r1(float_t _f_b)
{
	float_t tmp = tan(PI * _f_b / this->sampleRate);
	float_t c = (tmp - 1.f) / (tmp + 1.f);
	float_t r1 = acos(-c);
	return r1;
}

template<typename float_t>
float_t tvap<float_t>::filter(float_t x_n) {
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
	_sr1 = sin(_r1);
	_cr2 = cos(_r2);
	_sr2 = sin(_r2);
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
template<typename float_t>
inline float_t tvap<float_t>::unboundSat2(float_t x)
{
	float_t num = 2.0 * x;
	float_t denom = 1.0 + sqrt(1.0 + abs(4.0 * x));
	return num / denom;
}

template<typename float_t>
float_t tvap<float_t>::filter_fbmod(float_t x_n, float_t fb_f_pi, float_t fb_f_b)
{
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
	
	float_t highLimit = (this->sampleRate * 0.495);
	if (_f_pi_n >= highLimit) _f_pi_n = highLimit;
	if (_f_b_n >= highLimit) _f_b_n = highLimit;
	
	float_t _r1, _r2, _cr1, _cr2, _sr1, _sr2;
	float_t tmp[3];
	_r1 = f_b2r1(_f_b_n);
	_r2 = f_pi2r2(_f_pi_n);
	
	_cr1 = cos(_r1);
	_sr1 = sin(_r1);
	_cr2 = cos(_r2);
	_sr2 = sin(_r2);
	//tmp[0] = _x_n;
	tmp[1] = _cr2 * state.z1 - _sr2 * state.z2;
	//tmp[2] = _sr2 * _z1 + _cr2 * _z2;
	
	tmp[0] = _cr1 * x_n - _sr1 * tmp[1];
	tmp[1] = _sr1 * x_n + _cr1 * tmp[1];
	tmp[2] = _sr2 * state.z1 + _cr2 * state.z2;
	
	state.z1 = tmp[1];
	state.z2 = tmp[2];
	
	
	float_t fb_filt = dcFilt.filter(tmp[0]);   // used to feed output to modulate control inputs
	fb_filt = lp.tpt_lp(fb_filt);
	//fb_filt = 0.5f * (fb_filt + _fb_filt_z1);
	state.fb_proc = fb_filt;
	
	return tmp[0];
}

template class filter_abstract<float>;
template class filter_abstract<double>;
template class onePole<float>;
template class onePole<double>;
template class onePole_nonlinear<float>;
template class onePole_nonlinear<double>;
template class fourPole_LP_linear<float>;
template class fourPole_LP_linear<double>;
template class fourPole_LP_nonlinear<float>;
template class fourPole_LP_nonlinear<double>;
template class svf_prototype<float>;
template class svf_prototype<double>;
template class svf_lin_naive<float>;
template class svf_lin_naive<double>;
template class svf_nl_rk<float>;
template class svf_nl_rk<double>;
template class slewlim<float>;
template class slewlim<double>;
template class dcBlock<float>;
template class dcBlock<double>;
template class CTPTMoogFilterStage<float>;
template class CTPTMoogFilterStage<double>;
template class CTPTMoogLadderFilter<float>;
template class CTPTMoogLadderFilter<double>;
template class tvap<float>;
template class tvap<double>;

}	// namespace filters
}	// namespace nvs

#endif // NVS_FILTERS_H
