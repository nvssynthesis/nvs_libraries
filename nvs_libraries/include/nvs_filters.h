/*
  ==============================================================================
    nvs_filters.h
    Created: 6 Dec 2018 12:12:44am
    Author:  Nicholas Solem
  ==============================================================================
*/

#pragma once
#include "nvs_memoryless.h"
#include "nvs_matrix2x2.h"
#include <iostream>

#define USING_EIGEN 0
#if USING_EIGEN // for state space butterworth
#include <Eigen/Dense>
#endif

namespace nvs_filters {
template<typename floatType>
class filter_abstract{
    static_assert(std::is_floating_point<floatType>::value, "filter_abstract type must be floating point");
public:
    virtual ~filter_abstract() = 0;
    //==============================================================================
    virtual void setSampleRate(floatType sample_rate);
    virtual void clear() = 0;
    virtual void updateCutoff(floatType cutoff_target, floatType oneOverBlockSize) = 0;
    virtual void updateResonance(floatType res_target, floatType oneOverBlockSize) = 0;
    floatType cutoff_to_g(floatType cutoff); // could be static
    //==============================================================================
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
/*    onePole_nonlinear()
    {
        this->z1 = 0.0f;
        setSampleRate(44100.f);
    }
    onePole_nonlinear(floatType sample_rate)
    {
        this->z1 = 0.0f;
        setSampleRate(sample_rate);
    }*/
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
    /*
    void filter(floatType input)
    {
        using namespace nvs_memoryless;
        _outputs.lp = _state.lp + this->w_c * this->fs_inv * _state.bp;
        _outputs.hp = clamp<floatType>((input - (2 * R * _state.bp) - _state.lp), -10.f, 10.f);
        _outputs.bp = _state.bp + this->w_c * this->fs_inv * _outputs.hp;
        
        // update state
        _state.bp = clamp<floatType>(_outputs.bp, -10.f, 10.f);
        _state.lp = clamp<floatType>(_outputs.lp, -10.f, 10.f);
    }
    */
    void filter(floatType input);
private:

    floatType w_c, R, resonance;

/*
    nvs_memoryless::trigTables trig;
    floatType sampleRate, fs_inv; // why can't these be static?
    floatType _oneOverBlockSize, _cutoffTarget, _resonanceTarget;
    floatType z1;
*/
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

}   // namespace nvs_filters
//#include "../src/nvs_filters.cpp"

/*
namespace nvs_filters {
template<typename floatType>
void filter_abstract<floatType>::setSampleRate(floatType sample_rate)
{
    this->sampleRate = sample_rate;
    this->fs_inv = 1.f / sample_rate;
}
template<typename floatType>
floatType filter_abstract<floatType>::cutoff_to_g(floatType cutoff)
{
    //return cutoff * this->fs_inv * 0.5;    // no prewarp
    //return tan(cutoff * fs_inv * PI);      // expensive prewarp
    using namespace nvs_memoryless;
    if ((trig.tan_table) != NULL)
    {
        //floatType wc = TWOPI * cutoff;
        return (floatType)this->trig.tan_LUT(cutoff * fs_inv / 2.f);
    }
    else
        return 0.f;
}
//==================================================================================

template<typename floatType>
onePole<floatType>::onePole()
// :   z1(0) // member initializer 'z1' does not name a non-static data member or base class
{
    this->setSampleRate(44100.f);
    this->z1 = 0.0;
}
template<typename floatType>
onePole<floatType>::onePole(floatType sample_rate)
{
    this->setSampleRate(sample_rate);
    this->z1 = 0.0;
}
template<typename floatType>
void onePole<floatType>::clear()
{
    y_n = v_n = this->z1 = 0.0;
}
//==============================================================================
template<typename floatType>
void onePole<floatType>::updateCutoff()
{
    if (this->_cutoffTarget != w_c)
        this->w_c += (this->_cutoffTarget - this->w_c) * this->_oneOverBlockSize;
}
template<typename floatType>
void onePole<floatType>::updateCutoff(floatType cutoff_target, floatType oneOverBlockSize)
{
    this->w_c += (cutoff_target - this->w_c) * oneOverBlockSize;
}
template<typename floatType>
void onePole<floatType>::updateResonance(floatType res_target, floatType oneOverBlockSize)
{}//no resonance for onepole
//==============================================================================
template<typename floatType>
floatType onePole<floatType>::tpt_lp(floatType input)
{
    g = this->cutoff_to_g(this->w_c);
    v_n = (input - this->z1) * g / (1 + g);
    y_n = v_n + this->z1;
    this->z1 = y_n + v_n;
    // cliptest
    return nvs_memoryless::clamp<floatType>(y_n, -1000.f, 1000.f);
}
template<typename floatType>
floatType onePole<floatType>::tpt_lp(floatType input, floatType cutoff)
{
    g = this->cutoff_to_g(cutoff);
    v_n = (input - this->z1) * g / (1 + g);
    y_n = v_n + this->z1;
    this->z1 = y_n + v_n;
        // cliptest
    return nvs_memoryless::clamp<floatType>(y_n, -1000.f, 1000.f);
}
template<typename floatType>
floatType onePole<floatType>::tpt_hp(floatType input)
{
    return input - tpt_lp(input);
}
template<typename floatType>
floatType onePole<floatType>::tpt_hp(floatType input, floatType cutoff)
{
    return input - tpt_lp(input, cutoff);
}

//==============================================================================
template<typename floatType>
fourPole_LP_linear<floatType>::fourPole_LP_linear()
    :   u_n(0.f), s1(0.f), s2(0.f), s3(0.f), s4(0.f), S(0.f), y1(0.f), y2(0.f), y3(0.f), y4(0.f), G(0.f), k(0.f)
{
    floatType def_sampRate = (floatType)44100;
    setSampleRate(def_sampRate);
    H1.setSampleRate(def_sampRate);
    H2.setSampleRate(def_sampRate);
    H3.setSampleRate(def_sampRate);
    H4.setSampleRate(def_sampRate);
}
template<typename floatType>
fourPole_LP_linear<floatType>::fourPole_LP_linear(floatType sample_rate)
    :   u_n(0.f), s1(0.f), s2(0.f), s3(0.f), s4(0.f), S(0.f), y1(0.f), y2(0.f), y3(0.f), y4(0.f), G(0.f), k(0.f)
{
    setSampleRate(sample_rate);
    H1.setSampleRate(sample_rate);
    H2.setSampleRate(sample_rate);
    H3.setSampleRate(sample_rate);
    H4.setSampleRate(sample_rate);
}
template<typename floatType>
void fourPole_LP_linear<floatType>::initialize(floatType sample_rate)
{
    this->clear();
    setSampleRate(sample_rate);
    H1.setSampleRate(sample_rate);
    H2.setSampleRate(sample_rate);
    H3.setSampleRate(sample_rate);
    H4.setSampleRate(sample_rate);
}
template<typename floatType>
void fourPole_LP_linear<floatType>::updateOneOverBlockSize(floatType oneOverBlockSize)
{
    this->_oneOverBlockSize = oneOverBlockSize;
    H1._oneOverBlockSize = oneOverBlockSize;
    H2._oneOverBlockSize = oneOverBlockSize;
    H3._oneOverBlockSize = oneOverBlockSize;
    H4._oneOverBlockSize = oneOverBlockSize;
}
template<typename floatType>
void fourPole_LP_linear<floatType>::clear()
{
    H1.clear();
    H2.clear();
    H3.clear();
    H4.clear();
    s1 = s2 = s3 = s4 = u_n = y_n = S = 0.f;
}
template<typename floatType>
void fourPole_LP_linear<floatType>::updateCutoff()
{
    floatType local_cutoff = this->w_c;
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
template<typename floatType>
void fourPole_LP_linear<floatType>::updateCutoff(floatType cutoff_target, floatType oneOverBlockSize)
{
    this->w_c += (cutoff_target - this->w_c) * oneOverBlockSize;
    H1.updateCutoff(cutoff_target, oneOverBlockSize);
    H2.updateCutoff(cutoff_target, oneOverBlockSize);
    H3.updateCutoff(cutoff_target, oneOverBlockSize);
    H4.updateCutoff(cutoff_target, oneOverBlockSize);
}
template<typename floatType>
void fourPole_LP_linear<floatType>::updateResonance()
{
    if (_resonanceTarget != k)
        this->k += (this->_resonanceTarget - this->k) * this->_oneOverBlockSize;
}
template<typename floatType>
void fourPole_LP_linear<floatType>::updateResonance(floatType res_target, floatType oneOverBlockSize)
{
    //this->q += (res_target - this->q) * oneOverBlockSize;
    this->k += (res_target - this->k) * oneOverBlockSize;
}
template<typename floatType>
floatType fourPole_LP_linear<floatType>::tpt_fourpole(floatType input)
{
    g = cutoff_to_g(this->w_c);
    floatType g2 = g*g;
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
template<typename floatType>
floatType fourPole_LP_linear<floatType>::tpt_fourpole(floatType input, floatType cutoff)
{
    g = cutoff_to_g(cutoff);
    G = g * g * g * g;
    s1 = H1.z1;
    s2 = H2.z1;
    s3 = H3.z1;
    s4 = H4.z1;
    floatType g2 = g*g;
    S = g2*g*s1 + g2*s2 + g*s3 + s4;

    u_n = ((input) * (1 + k) - k * S) / (1 + k * G);

    y_n = H4.tpt_lp(H3.tpt_lp(H2.tpt_lp(H1.tpt_lp(u_n))));

    return y_n;
}

// so far quite tame, since the only nonlinearity is in the feedback path.
// TODO: convert each onepole into a nonlinear onepole
template<typename floatType>
fourPole_LP_nonlinear<floatType>::fourPole_LP_nonlinear()
:   iters(16),
u_n(0.f), s1(0.f), s2(0.f), s3(0.f), s4(0.f), S(0.f), y1(0.f), y2(0.f), y3(0.f), y4(0.f), G(0.f), k(0.f)
{
    floatType def_sampRate = (floatType)44100;
    setSampleRate(def_sampRate);
    H1.setSampleRate(def_sampRate);
    H2.setSampleRate(def_sampRate);
    H3.setSampleRate(def_sampRate);
    H4.setSampleRate(def_sampRate);
}
template<typename floatType>
fourPole_LP_nonlinear<floatType>::fourPole_LP_nonlinear(floatType sample_rate)
:   iters(16),
u_n(0.f), s1(0.f), s2(0.f), s3(0.f), s4(0.f), S(0.f), y1(0.f), y2(0.f), y3(0.f), y4(0.f), G(0.f), k(0.f)
{
    setSampleRate(sample_rate);
    H1.setSampleRate(sample_rate);
    H2.setSampleRate(sample_rate);
    H3.setSampleRate(sample_rate);
    H4.setSampleRate(sample_rate);
}
template<typename floatType>
void fourPole_LP_nonlinear<floatType>::initialize(floatType sample_rate)
{
    this->clear();
    setSampleRate(sample_rate);
    H1.setSampleRate(sample_rate);
    H2.setSampleRate(sample_rate);
    H3.setSampleRate(sample_rate);
    H4.setSampleRate(sample_rate);
}
template<typename floatType>
void fourPole_LP_nonlinear<floatType>::updateOneOverBlockSize(floatType oneOverBlockSize)
{
    this->_oneOverBlockSize = oneOverBlockSize;
    H1._oneOverBlockSize = oneOverBlockSize;
    H2._oneOverBlockSize = oneOverBlockSize;
    H3._oneOverBlockSize = oneOverBlockSize;
    H4._oneOverBlockSize = oneOverBlockSize;
}
template<typename floatType>
void fourPole_LP_nonlinear<floatType>::clear()
{
    H1.clear();
    H2.clear();
    H3.clear();
    H4.clear();
    s1 = s2 = s3 = s4 = u_n = y_n = S = 0.f;
}
template<typename floatType>
void fourPole_LP_nonlinear<floatType>::updateCutoff()
{
    floatType local_cutoff = this->w_c;
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
template<typename floatType>
void fourPole_LP_nonlinear<floatType>::updateCutoff(floatType cutoff_target, floatType oneOverBlockSize)
{
    this->w_c += (cutoff_target - this->w_c) * oneOverBlockSize;
    H1.updateCutoff(cutoff_target, oneOverBlockSize);
    H2.updateCutoff(cutoff_target, oneOverBlockSize);
    H3.updateCutoff(cutoff_target, oneOverBlockSize);
    H4.updateCutoff(cutoff_target, oneOverBlockSize);
}
template<typename floatType>
void fourPole_LP_nonlinear<floatType>::updateResonance()
{
    if (this->_resonanceTarget != k)
        this->k += (this->_resonanceTarget - this->k) * this->_oneOverBlockSize;
}
template<typename floatType>
void fourPole_LP_nonlinear<floatType>::updateResonance(floatType res_target, floatType oneOverBlockSize)
{
    //this->q += (res_target - this->q) * oneOverBlockSize;
    this->k += (res_target - this->k) * oneOverBlockSize;
}
template<typename floatType>
floatType fourPole_LP_nonlinear<floatType>::tpt_fourpole(floatType input)
{
    g = cutoff_to_g(this->w_c);
    floatType g2 = g * g;
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
template<typename floatType>
floatType fourPole_LP_nonlinear<floatType>::tpt_fourpole(floatType input, floatType cutoff)
{
    g = cutoff_to_g(cutoff);
    floatType g2 = g * g;
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
template<typename floatType>
svf_lin_naive<floatType>::svf_lin_naive()
:   w_c(200.f), R(1.f), resonance(1.f)
{
    this->sampleRate = 44100.f;
    this->fs_inv = 1.f / this->sampleRate;
    clear();
}
template<typename floatType>
void svf_lin_naive<floatType>::clear()
{
    this->_outputs = {0.f, 0.f, 0.f, 0.f };
    this->_state = { 0.f, 0.f };
}
template<typename floatType>
void svf_lin_naive<floatType>::setCutoff(floatType wc)
{
    this->w_c = wc;
    this->_cutoffTarget = w_c;
}
template<typename floatType>
void svf_lin_naive<floatType>::updateCutoff(floatType cutoff_target, floatType oneOverBlockSize)
{
    this->w_c += (cutoff_target - this->w_c) * oneOverBlockSize;
}
template<typename floatType>
void svf_lin_naive<floatType>::updateCutoff()
{
    this->w_c += (this->_cutoffTarget - this->w_c) * this->_oneOverBlockSize;
}
template<typename floatType>
void svf_lin_naive<floatType>::setResonance(floatType res)
{
    this->resonance = res;
    this->R = 1.f / res;
}
template<typename floatType>
void svf_lin_naive<floatType>::updateResonance(floatType res_target, floatType oneOverBlockSize)
{
    if (res_target > 0.9f)
        res_target = 0.9f;

    this->resonance += (res_target - this->resonance) * oneOverBlockSize;
    this->R = 1.f / res_target;
}
template<typename floatType>
void svf_lin_naive<floatType>::updateResonance()
{
    if (this->_resonanceTarget > 0.9f)
        this->_resonanceTarget = 0.9f;

    this->resonance += (this->_resonanceTarget - this->resonance) * this->_oneOverBlockSize;
    this->R = 1.f / this->_resonanceTarget;
}

void filter(floatType input)
{
    using namespace nvs_memoryless;
    _outputs.lp = _state.lp + this->w_c * this->fs_inv * _state.bp;
    _outputs.hp = clamp<floatType>((input - (2 * R * _state.bp) - _state.lp), -10.f, 10.f);
    _outputs.bp = _state.bp + this->w_c * this->fs_inv * _outputs.hp;
    
    // update state
    _state.bp = clamp<floatType>(_outputs.bp, -10.f, 10.f);
    _state.lp = clamp<floatType>(_outputs.lp, -10.f, 10.f);
} excluded
*/
/*
template<typename floatType>
void svf_lin_naive<floatType>::filter(floatType input)
{// Erbe version
    floatType c, d;
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

 nonlinear state-variable filter using fourth-order runge-kutta
 y[n+1] = y[n] + 1/6 * (k_1 + 2*k_2 + 2*k_3 + k_4)
 h = (1 / fs) / oversample_factor

 k_1 = h*f(t_n, y_n)
 k_2 = h*f(t_n + h/2, y_n + k_1/2)
 k_3 = h*f(t_n + h/2, y_n + k_2/2)
 k_4 = h*f(t_n + h, y_n + k_3)
 excluded
*/
/*
template<typename floatType>
svf_nl_rk<floatType>::svf_nl_rk()
: _oversample_factor(1), h(0.000022675736961),
 w_c(200.0), R(1.0), resonance(1.0)
{
    this->sampleRate = (floatType)44100;
    for (int i = 0; i < 2; i++)
    {
        deriv1[i] = deriv2[i] = deriv3[i] = deriv4[i] = 0.f;
    }
}
template<typename floatType>
void svf_nl_rk<floatType>::setSampleRate(floatType sample_rate)
{
    this->sampleRate = sample_rate;
    this->fs_inv = 1.f / sample_rate;
    this->h = 1.f / (_oversample_factor * sample_rate);
}
template<typename floatType>
void svf_nl_rk<floatType>::set_oversample(int oversample_factor)
{
    _oversample_factor = oversample_factor;
    h = 1.f / (oversample_factor * this->sampleRate);
}
template<typename floatType>
void svf_nl_rk<floatType>::clear()
{
    for (int i = 0; i < 2; i++)
    {
        deriv1[i] = deriv2[i] = deriv3[i] = deriv4[i] = 0.f;
        this->_outputs = {0.f, 0.f, 0.f, 0.f };
        this->_state = { 0.f, 0.f };
    }
}
template<typename floatType>
void svf_nl_rk<floatType>::setCutoff(floatType wc)
{
    this->w_c = wc;
}
template<typename floatType>
void svf_nl_rk<floatType>::updateCutoff(floatType cutoff_target, floatType oneOverBlockSize)
{
    this->w_c += (cutoff_target - this->w_c) * oneOverBlockSize;
}
template<typename floatType>
void svf_nl_rk<floatType>::updateCutoff()
{
    this->w_c += (this->_cutoffTarget - this->w_c) * this->_oneOverBlockSize;
}
template<typename floatType>
void svf_nl_rk<floatType>::setResonance(floatType res)
{
    if (res < 0.5f) res = 0.5f;
    this->resonance = res;
    this->R = 1.f / res;
}
template<typename floatType>
void svf_nl_rk<floatType>::updateResonance(floatType res_target, floatType oneOverBlockSize)
{
    if (res_target < 0.5f) res_target = 0.5f;
    this->resonance += (res_target - this->resonance) * oneOverBlockSize;
    this->R = 1.f / this->resonance;
}
template<typename floatType>
void svf_nl_rk<floatType>::updateResonance()
{
    if (this->_resonanceTarget < 0.5f) this->_resonanceTarget = 0.5f;
    this->resonance += (this->_resonanceTarget - this->resonance) * this->_oneOverBlockSize;
    this->R = 1.f / this->resonance;
}
template<typename floatType>
void svf_nl_rk<floatType>::filter(floatType input)
{
    using namespace nvs_memoryless;
    floatType hp, np;
    // overwritten states. [0] is bp, [1] is lp.
    floatType tempstate[2];
    
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
template<typename floatType>
slewlim<floatType>::slewlim()
    : _vOut(0.f)
{
    setSampleRate(44100.f);
}
template<typename floatType>
slewlim<floatType>::slewlim(floatType sample_rate)
    : _vOut(0.f)
{
    setSampleRate(sample_rate);
}
//============================================================
template<typename floatType>
void slewlim<floatType>::setSampleRate(floatType sample_rate)
{
    this->sampleRate = sample_rate;
    this->fs_inv = 1.f / sample_rate;
}
//============================================================
// immediate change
template<typename floatType>
void slewlim<floatType>::setRise(floatType rise)
{
    this->rise = rise;
    this->riseInc = (this->fs_inv * 1000.f) / (this->rise);
}
// change over block size
template<typename floatType>
void slewlim<floatType>::setRise()
{
    if (_riseTarget != rise)
    {
        this->rise += (_riseTarget - rise) * _oneOverBlockSize;
        this->riseInc = (this->fs_inv * 1000.f) / (rise);
    }
}
template<typename floatType>
void slewlim<floatType>::setRise(floatType riseTarget, floatType oneOverBlockSize)
{
    this->rise += (riseTarget - this->rise) * oneOverBlockSize;
    this->riseInc = (this->fs_inv * 1000.f) / (this->rise);
}
// immediate change
template<typename floatType>
void slewlim<floatType>::setFall(floatType fall)
{
    this->fall = fall;
    this->fallInc = (this->fs_inv * 1000.f) / (this->fall);
}
// change over block size
template<typename floatType>
void slewlim<floatType>::setFall()
{
    if (_fallTarget != fall)
    {
        this->fall += (_fallTarget - fall) * _oneOverBlockSize;
        this->fallInc = (this->fs_inv * 1000.f) / (fall);
    }
}
template<typename floatType>
void slewlim<floatType>::setFall(floatType fallTarget, floatType oneOverBlockSize)
{
    this->fall += (fallTarget - this->fall) * oneOverBlockSize;
    this->fallInc = (this->fs_inv * 1000.f) / (this->fall);
}
//============================================================
template<typename floatType>
floatType slewlim<floatType>::ASR(floatType gate)
{
    using namespace nvs_memoryless;
    
    setRise();
    setFall();
    
    if (_vOut < gate)
    {
        _vOut += riseInc;
        _vOut = clamp_high<floatType>(_vOut, gate);
    }
    else if (_vOut > gate)
    {
        _vOut -= fallInc;
        _vOut = clamp_low<floatType>(_vOut, gate);
    }
    
    return _vOut;
}
//=============================================================================
template<typename floatType>
dcBlock<floatType>::dcBlock()
{
    dcBlock((floatType)44100.0);
}
template<typename floatType>
dcBlock<floatType>::dcBlock(floatType sample_rate)
{
    this->setSampleRate((floatType) this->sampleRate);
    clear();
    R = 0.995;  //fixed R, should base on sample rate instead
}
template<typename floatType>
void dcBlock<floatType>::setSampleRate(floatType sampleRate){
    this->sample_rate = sampleRate;
    this->fs_inv = 1.f / sampleRate;
}
template<typename floatType>
void dcBlock<floatType>::clear()
{
    xz1 = yz1 = 0.f;
}
template<typename floatType>
void dcBlock<floatType>::updateR(floatType R_target, floatType oneOverBlockSize)
{
    R += (R_target - R) * oneOverBlockSize;
}
template<typename floatType>
floatType dcBlock<floatType>::filter(floatType x)
{
    yz1 = x - xz1 + R * yz1;
    xz1 = x;
    return yz1;
}
//===============================================================================
// PIRKLE IMPLEMENTATIONS (not my own work; used only for checking.)
template<typename floatType>
inline void CTPTMoogFilterStage<floatType>::initialize(floatType newSampleRate)
{
    // save
    sampleRate = newSampleRate;
    z1 = 0;
}
template<typename floatType>
void CTPTMoogFilterStage<floatType>::setFc(floatType fc)
{
    // prewarp the cutoff- these are bilinear-transform filters
    floatType wd = 2 * PI * fc;
    floatType fs_inv  = 1 / sampleRate;
    floatType wa = (2 / fs_inv) * tan(wd * fs_inv / 2);
    floatType g  = wa * fs_inv / 2;
    // calculate big G value; see Zavalishin p46 The Art of VA Design
    G = g / (1.0 + g);
}
template<typename floatType>
floatType CTPTMoogFilterStage<floatType>::doFilterStage(floatType xn)
{
    floatType v = (xn - z1) * G;
    floatType out = v + z1;
    z1 = out + v;
    return out;
}
template<typename floatType>
floatType CTPTMoogFilterStage<floatType>::getSampleRate()
{
    return sampleRate;
}
template<typename floatType>
floatType CTPTMoogFilterStage<floatType>::getStorageRegisterValue()
{
    return z1;
}
template<typename floatType>
inline void CTPTMoogLadderFilter<floatType>::initialize(floatType newSampleRate)
{
    filter1.initialize(newSampleRate);
    filter2.initialize(newSampleRate);
    filter3.initialize(newSampleRate);
    filter4.initialize(newSampleRate);
}
template<typename floatType>
inline void CTPTMoogLadderFilter<floatType>::calculateTPTCoeffs(floatType cutoff, floatType Q)
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
template<typename floatType>
floatType CTPTMoogLadderFilter<floatType>::doTPTMoogLPF(floatType xn)
{
    // calculate g
    floatType wd = 2 * PI * fc;
    floatType fs_inv  = 1 / (floatType)filter1.getSampleRate();
    floatType wa = (2 / fs_inv) * tan(wd * fs_inv / 2);
    floatType g  = wa * fs_inv / 2;
    floatType G = g * g * g * g;
    floatType S = g * g * g * filter1.getStorageRegisterValue() +
    g * g * filter2.getStorageRegisterValue() +
    g * filter3.getStorageRegisterValue() +
    filter4.getStorageRegisterValue();
    //uis input to filters, straight from book
    floatType u = (xn - k * S) / (1 + k * G);
    // four cascades using nested functions
    floatType filterOut = filter4.doFilterStage(filter3.doFilterStage (filter2.doFilterStage(filter1.doFilterStage(u))));
    // output
    return filterOut;
}
//==================================================================================

    time-variant allpass filter


template<typename floatType>
void tvap<floatType>::setSampleRate(floatType sample_rate)
{
    this->sampleRate = sample_rate;
    this->fs_inv = 1.f / sample_rate;
    lp.setSampleRate(sample_rate);
    lp.updateCutoff(sample_rate * 0.125f, 1.f);
}
template<typename floatType>
void tvap<floatType>::clear(){
    sp->z1 = sp->z2 = sp->fb_proc = 0.f;
};
template<typename floatType>
void tvap<floatType>::updateCutoff(floatType cutoff_target, floatType oneOverBlockSize) {
    if (cutoff_target < 0) cutoff_target = 0;
    f_pi += (cutoff_target - f_pi) * oneOverBlockSize;
    calc_b1();
}
template<typename floatType>
void tvap<floatType>::updateResonance(floatType res_target, floatType oneOverBlockSize){
    //if (f_b <= 0) f_b = 0.0000000001;
    f_b += (res_target - f_b) * oneOverBlockSize;
    f_b_to_b0();
}
// function aliases just to have more meaningful names
template<typename floatType>
void tvap<floatType>::update_f_pi (floatType f_pi_target, floatType oneOverBlockSize)
{   updateCutoff(f_pi_target, oneOverBlockSize); }
template<typename floatType>
void tvap<floatType>::update_f_b(floatType f_b_target, floatType oneOverBlockSize)
{   updateResonance(f_b_target, oneOverBlockSize); }

template<typename floatType>
void tvap<floatType>::calc_b1(void) {
    floatType d = -1 * cos((2.f * PI * f_pi) / this->sampleRate);
    floatType c = (tan(PI * f_b / this->sampleRate) - 1.f) / (tan(PI * f_b / this->sampleRate) + 1.f);
    floatType r1 = acos(-1.f * c);
    floatType r2 = acos(-1.f * d);
    b1 = cos(r1) * (1.f + cos(r2));
}
template<typename floatType>
void tvap<floatType>::f_b_to_b0(void) {
    floatType c = (tan(PI * f_b / this->sampleRate) - 1.f) / (tan(PI * f_b / this->sampleRate) + 1.f);
    floatType r1 = acos(-1.f * c);
    b0 = cos(r1);
}

template<typename floatType>
floatType tvap<floatType>::f_pi2r2(floatType _f_pi)
{
    floatType d = -1 * cos((2.f * PI * _f_pi) / this->sampleRate);
    floatType r2 = acos(-d);
    return r2;
}
template<typename floatType>
floatType tvap<floatType>::f_b2r1(floatType _f_b)
{
    floatType tmp = tan(PI * _f_b / this->sampleRate);
    floatType c = (tmp - 1.f) / (tmp + 1.f);
    floatType r1 = acos(-c);
    return r1;
}

template<typename floatType>
floatType tvap<floatType>::filter(floatType x_n) {
     floatType _y1 = state.y1;
    floatType _y2 = state.y2;
    floatType _x1 = state.x1;
    floatType _x2 = state.x2;
    floatType y_n = b0 * x_n - b1 * _x1 + _x2 + b1 * _y1 - b0 * _y2;
    state.y2 = _y1;
    state.y1 = y_n;
    state.x2 = _x1;
    state.x1 = x_n;
    return y_n; *//*
    floatType _r1, _r2, _cr1, _cr2, _sr1, _sr2;
    floatType tmp [3];
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
template<typename floatType>
inline floatType tvap<floatType>::unboundSat2(floatType x)
{
    floatType num = 2.0 * x;
    floatType denom = 1.0 + sqrt(1.0 + abs(4.0 * x));
    return num / denom;
}

template<typename floatType>
floatType tvap<floatType>::filter_fbmod(floatType x_n, floatType fb_f_pi, floatType fb_f_b)
{
    floatType _f_pi_n, _f_b_n;
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

    floatType highLimit = (this->sampleRate * 0.495);
    if (_f_pi_n >= highLimit) _f_pi_n = highLimit;
    if (_f_b_n >= highLimit) _f_b_n = highLimit;

    floatType _r1, _r2, _cr1, _cr2, _sr1, _sr2;
    floatType tmp[3];
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

    
    floatType fb_filt = dcFilt.filter(tmp[0]);   // used to feed output to modulate control inputs
    fb_filt = lp.tpt_lp(fb_filt);
    //fb_filt = 0.5f * (fb_filt + _fb_filt_z1);
    state.fb_proc = fb_filt;

    return tmp[0];
}
} // namespace nvs_filters
                   */
