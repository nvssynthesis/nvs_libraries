#pragma once
#include "nvs_filters.h"
#include "nvs_memoryless.h"

namespace nvs_delays {
class delay
{
public:
    delay()
    :   sampleRate(44100.f), T(1.f / 44100.f), delSize(8192)
    {
        isInitialized = false;
        delMask = delSize - 1;
        //delBuff = (float *)malloc(delSize * sizeof(float));
        delBuff = new float[delSize];
        wHead = 0;
        clear();
        isInitialized = true;
    }
    delay(unsigned int _delSize, float sample_rate)
    {
        isInitialized = false;
        // convert delay size to next power of 2
        _delSize = (unsigned int)(_delSize + 0.5);
        unsigned int log2 = 0;
        if (_delSize == 0) log2 =  13; // gives 8192
        if (_delSize != 1)
        {
            while (_delSize > 1)  {
                    _delSize >>= 1;
                    log2++; }
        }
        _delSize = log2;
        delSize =  pow(2, _delSize);
        delBuff = new float[delSize];

        delMask = delSize - 1;
        //delBuff = (float *)malloc(delSize * sizeof(float));
        wHead = 0;
        sampleRate = sample_rate;
        T = 1.f / sample_rate;
        clear(); 
        isInitialized = true;
    }
    ~delay() {delete[] delBuff;}



    void setSampleRate(float fs)
    {
        float ms = (delTimeSamps / sampleRate) * 1000.f;    // maintain former delay time
        sampleRate = (fs > 0.f) ? fs : 44100.f;
        T = 1.f / sampleRate;
        setDelayTimeMS(ms);
    }
    enum interp {
        floor = 0,
        linear,
        cubic
    };
    inline float getEffectiveDelayTimeSamps(){
        return delTimeSamps + 1;// > 1.f ? delTimeSamps : 1.f;
    }
    float (delay::*doTheTick)(float);
    // could also implement by multiplication of 1 or zero with different outputs
    // copy pointer to function and ensure it's really pointing there properly
    // also consider that this MUST be assigned already on initialized 
    float tick(float input) {
        //return (this->*doTheTick)(input);
        return (this->*doTheTick)(input);
    }
    
    float tick_floor(float input) {
        rHead = wHead + (int)getEffectiveDelayTimeSamps();
        rHead &= delMask;
        float y = *(delBuff + rHead);
        
        *(delBuff + wHead) = input;
        wHead--;
        wHead &= delMask;
        return y;
    }
    float tick_linear(float input) {
        rHead = wHead + (int)getEffectiveDelayTimeSamps();
        float frac = getEffectiveDelayTimeSamps() - (float)((int)getEffectiveDelayTimeSamps());
        rHead &= delMask;
        float y = *(delBuff + rHead) * (1.f - frac);
        int rHeadTmp = rHead + 1;
        rHeadTmp &= delMask;
        y += *(delBuff + rHeadTmp) * (frac);
        
        *(delBuff + wHead) = input;
        wHead--;
        wHead &= delMask;
        return y;
    }
    float tick_cubic(float input)   {
        rHead = wHead + (int)getEffectiveDelayTimeSamps();
        float frac = getEffectiveDelayTimeSamps() - (float)((int)getEffectiveDelayTimeSamps());

        float xm1, x0, x1, x2, a, b, c, y;
        xm1 = delBuff[(rHead - 1) & delMask];
        x0  = delBuff[(rHead + 0) & delMask];
        x1  = delBuff[(rHead + 1) & delMask];
        x2  = delBuff[(rHead + 2) & delMask];
        a = (3.f * (x0-x1) - xm1 + x2) / 2.f;
        b = 2.f*x1 + xm1 - (5.f*x0 + x2) / 2.f;
        c = (x1 - xm1) / 2.f;
        y = (((a * frac) + b) * frac + c) * frac + x0;
        
        *(delBuff + wHead) = input;
        wHead--;
        wHead &= delMask;
        return y;
    }
    
    void setInterpolation(int _interpType) {
        this->interpType = _interpType;
        switch (_interpType) {
        case 0: // look at address of doTheTick, make sure that's the address that actually gets returned in the tick function
            doTheTick = &delay::tick_floor;
        case 1:
            doTheTick = &delay::tick_linear;
        case 2:
            doTheTick = &delay::tick_cubic;
        }
    }
    // multiple tick functions for different interpolations,
    // then use function pointer to set tick to that function



    float tap(unsigned int writeOffset)
    {
        unsigned int pos = wHead + writeOffset;
        pos &= pos;
        return *(delBuff + pos);
    }
    void setDelayTimeSamps(float t){
        delTimeSamps = (t <= delSize) ? ((t > 0.f) ? t : 0.f) : delSize;
    }
    void setDelayTimeMS(float t)
    {
        t = t * 0.001f * sampleRate;
        setDelayTimeSamps(t);
    }
//    float getDelayTimeMS(){
//    }

    void updateDelayTimeMS(float target, float oneOverBlockSize)
    {
        target = target * 0.001f * sampleRate; // convert to samples
        float inc = (target - delTimeSamps) * oneOverBlockSize;
        inc = nvs_memoryless::atan_aprox(inc);
        delTimeSamps += inc;
        
    }
    unsigned int getDelaySize() { return delSize; }

    void clear() {
        int n = delSize;
        while (n > 0)
        {
            n -= 1;     // pre-decrement
            *(delBuff + n) = 0.f;
        }
    }
    bool ready(void)
    {
        return isInitialized;
    }
private:
    float delTimeSamps, sampleRate, T; // T is inverse of sample rate
    float *delBuff;
    unsigned int delSize, delMask;
    unsigned int wHead, rHead;
    short int interpType;
    bool isInitialized;
};


class allpass_delay : public delay {
public:
    allpass_delay()
    :   sampleRate(44100.f), T(1.f / 44100.f), delSize(8192), g(0.f), z(0.f)
    {
        isInitialized = false;
        delMask = delSize - 1;
        //delBuff = (float *)malloc(delSize * sizeof(float));
        delBuff = new float[delSize];
        wHead = 0;
        isInitialized = true;
    }
    allpass_delay(unsigned int _delSize, float sample_rate)
    {
        isInitialized = false;
        // convert delay size to next power of 2
        _delSize = (unsigned int)(_delSize + 0.5);
        unsigned int log2 = 0;
        if (_delSize == 0) log2 =  13; // gives 8192
        if (_delSize != 1)
        {
            while (_delSize > 1)  {
                _delSize >>= 1;
                log2++; }
        }
        _delSize = log2;
        delSize = pow(2, _delSize);
        delMask = delSize - 1;
        //delBuff = (float *)malloc(delSize * sizeof(float));
        delBuff = new float[delSize]; // 1 channel, 1 sample of memory
        wHead = 0;
        sampleRate = sample_rate;
        T = 1.f / sample_rate;
        
        isInitialized = true;
    }
    ~allpass_delay() {delete[] delBuff;}
    
    void update_g(float g_target, float oneOverBlockSize)
    {
        g += (g_target - g) * oneOverBlockSize;
    }

    float filter(float x_n)
    {
        float u_n, y_n;
        u_n = x_n + z * g;
        y_n = -1.f * (tick(u_n) + u_n * g);
        z = y_n;
        return y_n;
    }

private:
    float delTimeSamps, sampleRate, T;
    float *delBuff;
    unsigned int delSize, delMask;
    unsigned int wHead, rHead;
    bool isInitialized;
    
    float g, z;
};
}
