/*
  ==============================================================================

    nvs_LFO.h
    Created: 14 Dec 2018 10:48:27am
    Author:  Nicholas Solem

  ==============================================================================
*/

#include "include/nvs_LFO.h"

namespace nvssynthesis_lfo {
  
    template<typename T>
    simple_lfo<T>::simple_lfo() {
        setSampleRate(44100);
        _phase = 0.f;
    }
    template<typename T>
    simple_lfo<T>::simple_lfo(T sample_rate) {
        setSampleRate(sample_rate);
        _phase = 0.f;
    }
    template<typename T>
    void simple_lfo<T>::setSampleRate(T sample_rate) {
        this->sampleRate = sample_rate;
        this->fs_inv = 1.f / sample_rate;
    }
    template<typename T>
    T simple_lfo<T>::phasor() {
        using namespace nvs::memoryless;
        
        _phase += _freq * fs_inv;
        _phase = mod_1<T>(_phase);
        return _phase;
    }
    template<typename T>
    T simple_lfo<T>::saw() {
        return _phase * 2 - 1;
    }
    template<typename T>
    T simple_lfo<T>::square() {
        return _phase > 0.5f ? 1.f : 0.f;
    }
    template<typename T>
    T simple_lfo<T>::tri() {
        return _phase > 0.5f ? (_phase - 0.5) * 2 : (1.f - _phase - 0.5) * 2;
    }
    template<typename T>
    T simple_lfo<T>::sine(){
        using namespace nvs::memoryless;
        T sinewave = shapes.up_cos_LUT(_phase);
        return sinewave;
    }

    // takes in value [0 ... 3] and returns interpolated waveform
    // 0 = tri
    // 1 = saw
    // 2 = square
    // 3 = sine
    template<typename T>
    T simple_lfo<T>::multi(T waveform){
        using namespace nvs::memoryless;
        waveform = clamp<T>(waveform, 0, 3);
        int selector = (int)(floor(waveform));
        T interpolationAmt = mod_1<T>(waveform);

        switch (selector)
        {
            case 0:
            {
                return linterp<T>(tri(), saw(), interpolationAmt);
                break;
            }
            case 1:
            {
                return linterp<T>(saw(), square(), interpolationAmt);
                break;
            }
            case 2:
            {
                return linterp<T>(square(), sine(), interpolationAmt);
                break;
            }
            case 3:
            {
                return sine();
                break;
            }
            default:
            {
                return 0.f;
                break;
            }
        }
    }
}

