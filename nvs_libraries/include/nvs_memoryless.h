/*
  ==============================================================================

    nvs_dsp.h
    Created: 5 Dec 2018 9:24:03pm
    Author:  Nicholas Solem

  ==============================================================================
*/

#pragma once
#include <math.h>
#define PI      3.141592653589793
#define TWOPI   6.283185307179586

namespace nvs {
namespace memoryless {

template<typename t>
void metaparamA(t x, t *y)   {
    // y should be a 2-component array
    
    // i fed y with "metap+0"
    if (x < 0.333) {
        y[0] = x * 3;
        y[1] = 0;
    }
    else    {
        y[0] = (x < 0.666) ? 1 : 1 - (x - 0.666) * 3;
        y[1] = (x < 0.666) ? (x - 0.333) * 3 : 1;
    }
}

template <typename t>
void metaparamB(t x, t *y0, t *y1){

}

// linear interpolate. does not check for fade being between 0 and 1.
template<typename t>
t linterp(t x, t y, t fade)
{
    return x - (x * fade) + (y * fade);
}

// unipolar to bipolar (assuming input is [0..1]
template<typename t>
t unibi(t x)
{
    return (x * 2.f) - 1.f;
}
// bipolar to unipolar (assuming input is [-1..1]
template<typename t>
t biuni(t x)
{
    return (x + 1.f) * 0.5f;
}

template<typename t>
t clamp(t x, t min, t max)
{
    if (x > max)
        return max;
    else if (x < min)
        return min;
    else
        return x;
}
template<typename t>
t clamp1(t x)
{
    if (x > 1.f)
        return 1.f;
    else if (x < -1.f)
        return -1.f;
    else
        return x;
}
template<typename t>
t clamp_low(t x, t min)
{
    return (x < min) ? min : x;
}
template<typename t>
t clamp_high(t x, t max)
{
    return (x > max) ? max : x;
}

// based on https://stackoverflow.com/questions/14415753/wrap-value-into-range-min-max-without-division
template<typename t>
t div_wrap(t x, t xmin, t xmax)
{
    return x - (xmax - xmin) * floorf( x / (xmax - xmin));
}

template<typename t>
t wrap(t x, t xMin, t xMax)
{
    return (t)fmodf((fmodf(x - xMin, xMax - xMin) + xMax - xMin) , xMax - xMin) + xMin;
}

template<typename t>
t mod_1(t input)
{
    while (input >= 1.f)
        input -= 1.f;
    while (input < 0.f)
        input += 1.f;
    return input;
}

//    t unboundSat1(t x)
//    {
//        return logf(x + sqrt(x*x + 1.f));
//    }
//    t unboundSat2(t x)
//    {
//        t num = 2.f * x;
//        t denom = 1.f + sqrt(1.f + fabs(4.f * x));
//        return num / denom;
//    }

template<typename t>
t crush(t input, t depth)
{
    //t depth = pow(2, exponent);
    t crushed = int(input * depth) / depth;
    return crushed;
}

template<typename t>
t xOverOnePlusAbsX(t input)
{
    // y = x/(1 + |x|)
    return input / (1 + (input > 0 ? input : -input));
}

//MAKE TABLE INSTEAD! (now have done so!)
// Taylor approximation of sine
// sin(x) =~ x - x^3/3! + x^5/5! - x^7/7!
template<typename t>
t taySin(t input)
{
    t x = input * TWOPI - PI;
    return x - (pow(x, 3) * 0.16666) + (pow(x, 5) * 0.008333) - (pow(x, 7) * 0.0001984);
    //+ (pow(x,9)*0.0000027557); //(use if you want more accuracy)
}

// Taylor approximation of cosine
//cos(x) = 1 − x^2/2 + x^4/24 ...
template<typename t>
t tayCos(t input)
{
    t x = input * TWOPI - PI;
    return 1 - pow(x, 2) * 0.5 + pow(x, 4) * 0.0416666 - pow(x, 6) * 0.0013888 + pow(x,8) * 0.0000248;
}

// linear interpolate from sine to cosine. used in some waveshaping stuff.
template<typename t>
t sin2cos(t x, t phi)
{
    return linterp(taySin(x), tayCos(x), phi);
}

// more of an equal-power interpolation from sine to cosine. used in some waveshaping stuff.
template<typename t>
t sin2cos_ep(t x, t phi)
{
    return (sqrt(1 - phi)) * taySin(x) + (sqrt(phi)) * tayCos(x);
}
// more of an equal-power interpolation from sine to cosine. used in some waveshaping stuff.
/*template<typename t>
t sin2cos_ep(t x, t phi)
{
    return (sqrt(1 - phi)) * taySin(x) + (sqrt(phi)) * tayCos(x);
}*/
template<typename t>
t atan_aprox(t x)
{
    return x / (1 + (0.28125 * x) * (0.28125 * x));
}

template<typename t>
class cosTable {
public:
    cosTable()
    :   m_size(32768) {
        
        for (int i = 0; i < m_size; i++) {
            LUT[i] = cos(((t)(i*4) / (t)m_size) * TWOPI);
        }
        
    }
    t cos(t x) const
    {
        x = x * 0.318309886183791; // divide out PI
        x = wrap<t>(x, 0.0, 2.0); // repeats between 0 and 2PI
        
        t cosx;
        int idx;
        if (x < 0.5){
            idx = int(x * (t)m_size * 4.0);
            cosx = LUT[idx];
        }
        else if (x < 1.) {
            x -= 0.5;
            x = 0.5 - x;
            idx = int(x * (t)m_size * 4.0);
            cosx = -1. * LUT[idx];
        }
        else if (x < 1.5) {
            x -= 1.0;
            idx = int(x * (t)m_size * 4.0);
            cosx = -1.0 * LUT[idx];
        }
        else {
            x -= 1.5;
            x = 0.5 - x;
            idx = int(x * (t)m_size * 4.0);
            cosx = LUT[idx];
        }
        return cosx;
    }
    
private:
    int m_size;
    t LUT[32768];
};
template<typename t>
class trigTables
{
public:
    trigTables()
    : m_size(65536)
    {
        for (int i = 0; i < m_size; i++)
        {
            sin_table[i] = sin(((t)i / (t)m_size) * TWOPI);
            cos_table[i] = cos(((t)i / (t)m_size) * TWOPI);
            tan_table[i] = tan(((t)i / (t)m_size) * PI);
            
            tanh_table[i]= tanh(12 * ((double)(i - (m_size * 0.5)) / (double)m_size));
        }
    }
    t up_sin_LUT(t x) const // unipolar input [0..1] mapped to [-pi..pi]
    {
        x = (clamp<t>(x,0.f,1.f));
        t frac = ceil(x) - x;
        return linterp<t>(sin_table[(int)floor(x * m_size)],
                              sin_table[(int)ceil (x * m_size)],
                              frac);
    }
    t bp_sin_LUT(t x) const // bipolar input [-1..1] mapped to [-pi..pi]
    {
        x = wrap<t>(x, -1.f, 1.f);
//x = biuni<t>(clamp1<t>(x));
        t frac = ceil(x) - x;
        return linterp<t>(sin_table[(int)floor(x * m_size)],
                              sin_table[(int)ceil (x * m_size)],
                              frac);
    }
    t up_cos_LUT(t x) const  // unipolar input [0..1] mapped to [-pi..pi]
    {
        x = (wrap<t>(x,0.f,1.f));
        t frac = ceil(x) - x;
        return linterp<t>(cos_table[(int)floor(x * m_size)],
                              cos_table[(int)ceil (x * m_size)],
                              frac);
    }
    t bp_cos_LUT(t x) const // bipolar input [-1..1] mapped to [-pi..pi]
    {
        x = wrap<t>(x, -1.f, 1.f);
//        x = 0.f;
        t frac = ceil(x) - x;
        int fl = (int)floor(x * m_size);
        int ce = (int)ceil (x * m_size);
        ce = ce + (fl == ce);
        return linterp<t>(cos_table[fl],
                              cos_table[ce],
                              frac);
    }
    
    // send in a value 0 to 1 (or 0 to 0.5 for Nyquist if used for cutoff to freq map), this already has
    // multiplied it by pi.
    t tan_LUT(t x) const //unipolar input for now
    {
        x = clamp<t>(x,0.f,1.f);
        t frac = ceil(x) - x;
        return linterp<t>(tan_table[(int)floor(x * m_size)],
                              tan_table[(int)ceil (x * m_size)],
                              frac);
    }

    // send in a value -6 to 6, because beyond this the output is essentially sign(x)
    double tanh_LUT(double x) const
    {
        x = clamp<double>(x, -12.f, 12.f);
        double x_index = (x + 6.f) * 5461.333333333333333;
        x_index = clamp<double>(x_index, 0.f, 65535.f);
        double frac = ceil(x_index) - x_index;
        return linterp<double>(tanh_table[(int)floor(x_index)],
                              tanh_table[(int)ceil (x_index)],
                              frac);
        //return frac;
    }
    
    const int m_size;
    //bool table_made;
    t tan_table[65536], sin_table[65536], cos_table[65536];
    double tanh_table[65536];
};

}   // namespace memoryless
}   // namespace nvs
