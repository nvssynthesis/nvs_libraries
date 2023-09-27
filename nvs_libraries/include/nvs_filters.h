/*
  ==============================================================================
    nvs_filters.h
    Created: 6 Dec 2018 12:12:44am
    Author:  Nicholas Solem
  ==============================================================================
*/

/**
 TODO:
	-optimize tvap:
		-replace all instances of sin() and cos() with faster versions (lookup table if it remains stable)
	-remove #include <iostream>
	**-remove filter_abstract::z1
 */ 

#pragma once
#include "nvs_memoryless.h"
#include "nvs_matrix2x2.h"
#include <iostream>

#define USING_EIGEN 0
#if USING_EIGEN	// for state space butterworth
#include <Eigen/Dense>
#endif

namespace nvs	{
namespace filters {

template<typename floatType>
class filter_abstract{
	static_assert(std::is_floating_point<floatType>::value, "filter_abstract type must be floating point");
public:
	virtual ~filter_abstract() = 0;
	//==============================================================================
	virtual void setSampleRate(floatType sample_rate){
		this->sampleRate = sample_rate;
		this->fs_inv = 1.f / sample_rate;
	}
	virtual void clear() = 0;
	virtual void updateCutoff(floatType cutoff_target, floatType oneOverBlockSize) = 0;
	virtual void updateResonance(floatType res_target, floatType oneOverBlockSize) = 0;
	floatType cutoff_to_g(floatType cutoff){ // could be static (or free function)
		//return cutoff * this->fs_inv * 0.5;    // no prewarp
		//return tan(cutoff * fs_inv * PI);      // expensive prewarp
		using namespace nvs::memoryless;
		if ((trig.tan_table) != NULL)
		{
			//floatType wc = TWOPI * cutoff;
			return (floatType)this->trig.tan_LUT(cutoff * fs_inv / 2.f);
		}
		else
			return 0.f;
	}
private:
	nvs::memoryless::trigTables<floatType> trig;
	floatType sampleRate, fs_inv; // why can't these be static?
	floatType _oneOverBlockSize, _cutoffTarget, _resonanceTarget;
	floatType z1; // turn into variable-length array with template<unsigned int>
};
template<typename floatType>
inline filter_abstract<floatType>::~filter_abstract() { }

//==================================================================================

template<typename floatType>
class onePole   :   public filter_abstract<floatType>
{
public:
	onePole();
	onePole(floatType sample_rate);
	void clear();
	//==============================================================================
	void updateCutoff();
	void updateCutoff(floatType cutoff_target, floatType oneOverBlockSize);
	void updateResonance(floatType res_target, floatType oneOverBlockSize);
	//==============================================================================
	floatType tpt_lp(floatType input);
	floatType tpt_lp(floatType input, floatType cutoff);
	floatType tpt_hp(floatType input);
	floatType tpt_hp(floatType input, floatType cutoff);
	//==============================================================================
private:
	floatType g, v_n, y_n;
	floatType w_c;
};
//==============================================================================
#if !USING_EIGEN

template<typename floatType>
class butterworth2p :   public filter_abstract<floatType>
{
public:
	butterworth2p() :   butterworth2p(44100.0){}
	butterworth2p(floatType sample_rate) /*   :   A(2,2), B(2,1), C(1,2), D(1,1), x(2), y(1)*/{
		this->setSampleRate(sample_rate);
		clear();
		A.b = 1.0;
		B.b = 1.0;
		D = 0.0;
		if ((this->trig.tan_table) == NULL){
			std::cout << "TAN TABLE NULL!\n";
		}
	}
	void clear(){
		x.a = x.b = 0.f;
		y = 0.f;
	}
	//==============================================================================
	void updateCutoff() {
		if (this->_cutoffTarget != this->w_c)
			this->w_c += (this->_cutoffTarget - this->w_c) * this->_oneOverBlockSize;
		calcCoefs(this->w_c);
	}
	void updateCutoff(floatType cutoff_target, floatType oneOverBlockSize){
		this->w_c += (cutoff_target - this->w_c) * oneOverBlockSize;
		calcCoefs(this->w_c);
	}
	void updateResonance(floatType res_target, floatType oneOverBlockSize){std::cout << "butterworth has fixed Q\n";}
	
	void calcCoefs(const floatType cutoff){
		//        const floatType fr = this->sampleRate / cutoff;
		//        floatType omega = tan(M_PI / fr);
		const floatType omega = this->trig.tan_LUT(cutoff / this->sampleRate);
		const floatType omega2 = omega * omega;
		//        static const floatType cosPiOver4 = 0.707106781186548;
		const floatType twoCosPiOver4xOmega =  1.414213562373095 * omega;
		const floatType c = 1.0 + twoCosPiOver4xOmega + omega2;
		
		const floatType b0 = omega2 / c;
		const floatType b1 = 2.0 * b0;
		//        const floatType b2 = b0;
		//        const floatType a0 = 1.0;
		const floatType a1 = (2.0 * (omega2 - 1.0)) / c;
		const floatType a2 = (1.0 - twoCosPiOver4xOmega + omega2) / c;
		
		A.c = -a2;
		A.d = -a1;
		C.a = b0 - (a2*b0);
		C.b = b1 - (a1*b0);
		D = b0;
	}
	floatType filter(floatType x_n){
		using namespace nvs_matrix;
		y = vec2::crossProduct(C, x);
		y += D * x_n;
		x = vec2::add(mat2x2::matXvec(A, x), vec2::scale(B, x_n));
		
		return y;
	}
	floatType filter(floatType x_n, floatType cutoff){
		calcCoefs(cutoff);
		return filter(x_n);
	}
	
private:
	nvs_matrix::mat2x2 A;
	nvs_matrix::vec2 B;
	nvs_matrix::vec2 C;
	floatType D;
	nvs_matrix::vec2 x;
	floatType y;
	
	floatType b0, b1, b2, a0, a1, a2;
	floatType w_c;
};
#endif
//==============================================================================
#if USING_EIGEN
template<typename floatType>
class butterworth2p :   public filter_abstract<floatType>
{
public:
	butterworth2p() :   butterworth2p(44100.0){}
	butterworth2p(floatType sample_rate) /*   :   A(2,2), B(2,1), C(1,2), D(1,1), x(2), y(1)*/{
		this->setSampleRate(sample_rate);
		clear();
		A << 0.0, 1.0, -0.0, -0.0;
		B << 0.0, 1.0;
		C << 0.0, 0.0;
		D << 0.0;
		if ((this->trig.tan_table) == NULL){
			std::cout << "TAN TABLE NULL!\n";
		}
	}
	void clear(){
		x[0] = x[1] = 0.f;
		y[0] = 0.f;
	}
	//==============================================================================
	void updateCutoff() {
		if (this->_cutoffTarget != this->w_c)
			this->w_c += (this->_cutoffTarget - this->w_c) * this->_oneOverBlockSize;
		calcCoefs(this->w_c);
	}
	void updateCutoff(floatType cutoff_target, floatType oneOverBlockSize){
		this->w_c += (cutoff_target - this->w_c) * oneOverBlockSize;
		calcCoefs(this->w_c);
	}
	void updateResonance(floatType res_target, floatType oneOverBlockSize){std::cout << "butterworth has fixed Q\n";}    // butterworth has fixed resonance
	
	void calcCoefs(const floatType cutoff){
		//        const floatType fr = this->sampleRate / cutoff;
		//        floatType omega = tan(M_PI / fr);
		const floatType omega = this->trig.tan_LUT(cutoff / this->sampleRate);
		const floatType omega2 = omega * omega;
		//        static const floatType cosPiOver4 = 0.707106781186548;
		const floatType twoCosPiOver4xOmega =  1.414213562373095 * omega;
		const floatType c = 1.0 + twoCosPiOver4xOmega + omega2;
		
		const floatType b0 = omega2 / c;
		const floatType b1 = 2.0 * b0;
		//        const floatType b2 = b0;
		//        const floatType a0 = 1.0;
		const floatType a1 = (2.0 * (omega2 - 1.0)) / c;
		const floatType a2 = (1.0 - twoCosPiOver4xOmega + omega2) / c;
		
		A(1,0) = -a2;
		A(1,1) = -a1;
		
		C(0) = b0 - (a2*b0);
		C(1) = b1 - (a1*b0);
		//        std::cout << C << "\n";
		D(0) = b0;
	}
	floatType filter(floatType x_n){
		y = C*x + D*x_n;
		x = A*x + B*x_n;
		
		return y[0];
	}
	floatType filter(floatType x_n, floatType cutoff){
		calcCoefs(cutoff);
		y = C*x + D*x_n;
		x = A*x + B*x_n;
		
		return y[0];
	}
	
private:
	//    Eigen::MatrixXd A, B, C, D;
	Eigen::Matrix<floatType, 2, 2> A;
	Eigen::Matrix<floatType, 2, 1> B;
	Eigen::Matrix<floatType, 1, 2> C;
	Eigen::Matrix<floatType, 1, 1> D;
	//    Eigen::VectorXd x, y;
	Eigen::Vector<floatType, 2> x;
	Eigen::Vector<floatType, 1> y;
	floatType b0, b1, b2, a0, a1, a2;
	floatType w_c;
};
#endif

// NOTHING SO FAR.
template<typename floatType>
class onePole_nonlinear   :   public onePole<floatType>
{
public:
private:
};

template<typename floatType>
class fourPole_LP_linear    :   public filter_abstract<floatType>
{
public:
	fourPole_LP_linear();
	fourPole_LP_linear(floatType sample_rate);
	void initialize(floatType sample_rate);
	void updateOneOverBlockSize(floatType oneOverBlockSize);
	void clear();
	void updateCutoff();
	void updateCutoff(floatType cutoff_target, floatType oneOverBlockSize);
	void updateResonance();
	void updateResonance(floatType res_target, floatType oneOverBlockSize);
	
	floatType tpt_fourpole(floatType input);
	floatType tpt_fourpole(floatType input, floatType cutoff);
	
	floatType _resonanceTarget;
private:
	onePole<floatType> H1, H2, H3, H4;
	floatType u_n, y_n, s1, s2, s3, s4, S, y1, y2, y3, y4, g, g_denom, G, k;
	floatType w_c, q;
};

// so far quite tame, since the only nonlinearity is in the feedback path. 
// TODO: convert each onepole into a nonlinear onepole
template<typename floatType>
class fourPole_LP_nonlinear    :   public filter_abstract<floatType>
{
public:
	fourPole_LP_nonlinear();
	fourPole_LP_nonlinear(floatType sample_rate);
	void initialize(floatType sample_rate);
	void updateOneOverBlockSize(floatType oneOverBlockSize);
	void clear();
	void updateCutoff();
	void updateCutoff(floatType cutoff_target, floatType oneOverBlockSize);
	
	void updateResonance();
	void updateResonance(floatType res_target, floatType oneOverBlockSize);
	
	
	floatType tpt_fourpole(floatType input);
	floatType tpt_fourpole(floatType input, floatType cutoff);    
	//    floatType _resonanceTarget;
private:
	onePole<floatType> H1, H2, H3, H4;
	nvs::memoryless::trigTables<floatType> tables;
	int iters;
	floatType u_n, y_n, s1, s2, s3, s4, S, y1, y2, y3, y4, g, g_denom, G, k;
	floatType w_c, q;
};

template<typename floatType>
class svf_prototype
{
public:
	floatType lp() {return _outputs.lp;}
	floatType bp() {return _outputs.bp;}
	floatType hp() {return _outputs.hp;}
	floatType np() {return _outputs.np;}
	
protected:
	struct outputs
	{
		floatType lp, bp, hp, np;
	} _outputs = { 0.f, 0.f, 0.f, 0.f };
	struct state
	{
		floatType lp, bp;
	} _state = { 0.f, 0.f };
};

// linear state variable filter using 'naive' integrators (i.e., Euler backward difference integration)
template<typename floatType>
class svf_lin_naive     :   public filter_abstract<floatType>, svf_prototype<floatType>
{
public:
	//==============================================================================
	svf_lin_naive();
	void clear();
	void setCutoff(floatType wc);
	void updateCutoff(floatType cutoff_target, floatType oneOverBlockSize);
	void updateCutoff();
	void setResonance(floatType res);
	void updateResonance(floatType res_target, floatType oneOverBlockSize);
	void updateResonance();
	void filter(floatType input);
private:
	
	floatType w_c, R, resonance;
	
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
template<typename floatType>
class svf_nl_rk :   public filter_abstract<floatType>, public svf_prototype<floatType>
{
public:
	svf_nl_rk();
	void setSampleRate(floatType sample_rate);
	void set_oversample(int oversample_factor);
	void clear();
	
	void setCutoff(floatType wc);
	void updateCutoff(floatType cutoff_target, floatType oneOverBlockSize);
	void updateCutoff();
	void setResonance(floatType res);
	void updateResonance(floatType res_target, floatType oneOverBlockSize);
	void updateResonance();
	void filter(floatType input);
private:
	int _oversample_factor;
	floatType h, w_c, R, resonance;
	// k_1 through k_4. for each, [0] is bp, [1] is lp.
	floatType deriv1[2], deriv2[2], deriv3[2], deriv4[2]; 
};

//==============================================================================
template<typename floatType>
class slewlim
{
public:
	slewlim();
	slewlim(floatType sample_rate);
	~slewlim()  { }
	//============================================================
	void setSampleRate(floatType sample_rate);
	//============================================================
	// immediate change
	void setRise(floatType rise);
	// change over block size
	void setRise();
	void setRise(floatType riseTarget, floatType oneOverBlockSize);
	// immediate change
	void setFall(floatType fall);
	// change over block size
	void setFall();
	void setFall(floatType fallTarget, floatType oneOverBlockSize);
	//============================================================
	floatType ASR(floatType gate);
	
	floatType _riseTarget, _fallTarget;
	floatType _oneOverBlockSize;
private:
	floatType sampleRate, fs_inv;
	
	// 'Inc' variables tell change per sample.
	floatType rise, riseInc, fall, fallInc, _vOut;
};

template<typename floatType>
class dcBlock   :   public filter_abstract <floatType>
{
public:
	dcBlock();
	dcBlock(floatType sampleRate);
	void setSampleRate(floatType sampleRate);
	void clear();
	void updateCutoff(floatType cutoff_target, floatType oneOverBlockSize){ }
	void updateResonance(floatType res_target, floatType oneOverBlockSize){ }
	void updateR(floatType R_target, floatType oneOverBlockSize);
	floatType filter(floatType x);
private:
	floatType R, xz1, yz1;
	floatType sample_rate;
};

//===============================================================================
// PIRKLE IMPLEMENTATIONS (not my own work; used only for checking.)
template<typename floatType>
class CTPTMoogFilterStage
{
public:
	CTPTMoogFilterStage(){}
	~CTPTMoogFilterStage(){}
protected:
	floatType G;
	floatType scalar;
	floatType sampleRate;
	floatType z1;
public:
	inline void initialize(floatType newSampleRate);
	void setFc(floatType fc);
	
	floatType doFilterStage(floatType xn);
	floatType getSampleRate();
	floatType getStorageRegisterValue();
};
template <typename floatType>
class CTPTMoogLadderFilter
{
public:
	CTPTMoogLadderFilter(){}
	~CTPTMoogLadderFilter(){}
protected:
	CTPTMoogFilterStage<floatType> filter1;
	CTPTMoogFilterStage<floatType> filter2;
	CTPTMoogFilterStage<floatType> filter3;
	CTPTMoogFilterStage<floatType> filter4;
	floatType k; // Q control
	floatType fc; // fc control
public:
	inline void initialize(floatType newSampleRate);
	inline void calculateTPTCoeffs(floatType cutoff, floatType Q);
	floatType doTPTMoogLPF(floatType xn);
};
//==================================================================================
/* 
 time-variant allpass filter
 */
template<typename floatType>
class tvap  :   public filter_abstract<floatType>
{
public:
	tvap() {}
	~tvap() {}
	//==============================================================================
	virtual void setSampleRate(floatType sample_rate);
	
	void clear();
	void updateCutoff(floatType cutoff_target, floatType oneOverBlockSize);
	void updateResonance(floatType res_target, floatType oneOverBlockSize);
	// function aliases just to have more meaningful names
	void update_f_pi (floatType f_pi_target, floatType oneOverBlockSize);
	void update_f_b(floatType f_b_target, floatType oneOverBlockSize);
	
	void calc_b1(void);
	void f_b_to_b0(void);
	
	floatType f_pi2r2(floatType _f_pi);
	floatType f_b2r1(floatType _f_b);
	
	floatType filter(floatType x_n);
	
	// should be in memoryless but got linker error
	inline floatType unboundSat2(floatType x);
	
	floatType filter_fbmod(floatType x_n, floatType fb_f_pi, floatType fb_f_b);
	nvs::memoryless::trigTables<floatType> _LUT;
	
protected:
	typedef struct tvapstate {
		//floatType x1, x2, y1, y2;
		floatType z1, z2;
		floatType fb_proc;  // processed fed back output sample
	} _tvapstate;
	_tvapstate state = {.z1 = 0.f, .z2 = 0.f, .fb_proc = 0.f};
	_tvapstate *sp = &state;
	
	dcBlock<floatType> dcFilt;
	onePole<floatType> lp;
private:
	floatType f_pi, f_b;
	floatType b0, b1;
};

}   // namespace filters
}	// namespace nvs
