/*
  ==============================================================================

    nvs_LFO.h
    Created: 14 Dec 2018 10:48:27am
    Author:  Nicholas Solem

  ==============================================================================
*/

#pragma once
#include <type_traits>
#include "nvs_memoryless.h"

namespace nvs {
namespace lfo {
template<typename T>
class simple_lfo {
	static_assert(std::is_floating_point<T>::value, "simple_lfo type must be floating point");
public:
	simple_lfo();
	simple_lfo(T sample_rate);
	void setSampleRate(T sample_rate);
	
	T phasor(); // use to update state (phase)
	[[nodiscard]] T saw() const;
	[[nodiscard]] T square() const;
	[[nodiscard]] T tri() const;
	[[nodiscard]] T sine() const;
	
	[[nodiscard]] T phasor_offset(T offset) const;
	[[nodiscard]] T saw_offset(T offset) const;
	[[nodiscard]] T square_offset(T offset) const;
	[[nodiscard]] T tri_offset(T offset) const;
	[[nodiscard]] T sine_offset(T offset) const;
	// takes in value [0 ... 3] and returns interpolated waveform
	// 0 = tri
	// 1 = saw
	// 2 = square
	// 3 = sine
	[[nodiscard]] T multi(T waveform) const;
	void reset();
	[[nodiscard]] bool crossedOver() const;
	
	T _freq;
private:
	T sampleRate, fs_inv;
	T _phase, _lastPhase;
	const nvs::memoryless::trigTables<T> shapes;
};
} // namespace lfo
} // namespace nvs

namespace nvs {
namespace lfo {

template<typename T>
simple_lfo<T>::simple_lfo() {
	simple_lfo(44100);
	/*setSampleRate(44100);
	 _phase = 0.0;
	 // better: epsilon for <T>. better still: based on frequency, but that doesn't exist on init
	 _lastPhase = 0.999999;*/

}
template<typename T>
simple_lfo<T>::simple_lfo(T sample_rate) {
	setSampleRate(sample_rate);
	_phase = 0.0;
	// better: epsilon for <T>. better still: based on frequency, but that doesn't exist on init
	_lastPhase = 0.999999;
}
template<typename T>
void simple_lfo<T>::setSampleRate(T sample_rate) {
	this->sampleRate = sample_rate;
	this->fs_inv = 1.0 / sample_rate;
}
template<typename T>
T simple_lfo<T>::phasor() {
	using namespace nvs::memoryless;
	_lastPhase = _phase;
	_phase += _freq * fs_inv;
	_phase = mod_1<T>(_phase);
	return _phase;
}
template<typename T>
T simple_lfo<T>::saw() const {
	return _phase * 2 - 1;
}
template<typename T>
T simple_lfo<T>::square() const {
	return _phase > 0.5f ? 1.f : 0.f;
}
template<typename T>
T simple_lfo<T>::tri() const {
	return _phase > 0.5f ? (_phase - 0.5) * 2 : (1.f - _phase - 0.5) * 2;
}
template<typename T>
T simple_lfo<T>::sine() const {
	using namespace nvs::memoryless;
	T sinewave = shapes.up_cos_LUT(_phase);
	return sinewave;
}
template<typename T>
T inline simple_lfo<T>::phasor_offset(T offset) const {
	return nvs::memoryless::mod_1<T>(_phase + offset);
}
template<typename T>
T simple_lfo<T>::saw_offset(T offset) const {
	return phasor_offset(offset) * 2 - 1;
}
template<typename T>
T simple_lfo<T>::square_offset(T offset) const {
	return phasor_offset(offset) > 0.5f ? 1.f : 0.f;
}
template<typename T>
T simple_lfo<T>::tri_offset(T offset) const {
	T phase = phasor_offset(offset);
	return phase > 0.5f ? (phase - 0.5) * 2 : (1.f - phase - 0.5) * 2;
}
template<typename T>
T simple_lfo<T>::sine_offset(T offset) const {
	using namespace nvs::memoryless;
	T sinewave = this->shapes.up_cos_LUT(phasor_offset(offset));
	return sinewave;
}
// takes in value [0 ... 3] and returns interpolated waveform
// 0 = tri
// 1 = saw
// 2 = square
// 3 = sine
template<typename T>
T simple_lfo<T>::multi(T waveform) const {
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
template<typename T>
void simple_lfo<T>::reset(){
	_phase = 0.0;
}
template<typename T>
bool simple_lfo<T>::crossedOver() const {
	return (_lastPhase > _phase);
}
}	// namespace lfo
}	// namespace nvs
