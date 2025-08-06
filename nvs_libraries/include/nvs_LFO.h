/*
  ==============================================================================

    nvs_LFO.h
    Created: 14 Dec 2018 10:48:27am
    Author:  Nicholas Solem

  ==============================================================================
*/

#pragma once
#include <type_traits>
#include "nvs_fast_math_impl.h"
#include "nvs_gen.h"

namespace nvs::lfo {
template<typename T>
class simple_lfo {
	static_assert(std::is_floating_point<T>::value, "simple_lfo type must be floating point");
public:
	simple_lfo() noexcept;
	simple_lfo(T sample_rate) noexcept;
	void setSampleRate(T sample_rate) noexcept;
	T getSampleRate() const noexcept{ return sampleRate; }
	
	T phasor(); // use to update state (phase)
	[[nodiscard]] T saw() const noexcept;
	[[nodiscard]] T square() const noexcept;
	[[nodiscard]] T tri() const noexcept;
	[[nodiscard]] T sine() const noexcept;
	
	[[nodiscard]] T phasor_offset(T offset) const noexcept;
	[[nodiscard]] T saw_offset(T offset) const noexcept;
	[[nodiscard]] T square_offset(T offset) const noexcept;
	[[nodiscard]] T tri_offset(T offset) const noexcept;
	[[nodiscard]] T sine_offset(T offset) const noexcept;
	// takes in value [0 ... 3] and returns interpolated waveform
	// 0 = tri
	// 1 = saw
	// 2 = square
	// 3 = sine
	[[nodiscard]] T multi(T waveform) const noexcept;
	void reset() noexcept;
	
	
	[[deprecated("Pretty much guaranteed to be buggy in some way; use nvs::gen::ramp2trig instead.")]]
	[[nodiscard]] bool crossedOver() const noexcept;
	
	T _freq {0.0};
private:
	T sampleRate, fs_inv;
	T _phase{};
	
	[[deprecated]] T _lastPhase{};
};
} // namespace nvs::lfo

namespace nvs::lfo {

template<typename T>
simple_lfo<T>::simple_lfo() noexcept {
	simple_lfo(44100);
}
template<typename T>
simple_lfo<T>::simple_lfo(T sample_rate) noexcept {
	setSampleRate(sample_rate);
	_phase = 0.0;
	// better: epsilon for <T>. better still: based on frequency, but that doesn't exist on init
	_lastPhase = 0.999999;
}
template<typename T>
void simple_lfo<T>::setSampleRate(T sample_rate) noexcept {
	assert (0.0 < sample_rate);
	this->sampleRate = sample_rate;
	this->fs_inv = 1.0 / sample_rate;
}
template<typename T>
T simple_lfo<T>::phasor() {
	using namespace nvs::memoryless;
	_lastPhase = _phase;
	_phase += _freq * fs_inv;
	_phase = nvs::memoryless::mspWrap(_phase);
	return _phase;
}
template<typename T>
T simple_lfo<T>::saw() const noexcept {
	return nvs::memoryless::unibi(_phase);
}
template<typename T>
T simple_lfo<T>::square() const noexcept {
	return _phase > 0.5f ? 1.f : -1.f;
}
template<typename T>
T simple_lfo<T>::tri() const noexcept {
	return nvs::memoryless::unibi(_phase > 0.5f ? (_phase - 0.5) * 2 : (1.f - _phase - 0.5) * 2);
}
template<typename T>
T simple_lfo<T>::sine() const noexcept {
	using namespace nvs::memoryless;
	return math_impl::cos(_phase * math_impl::two_pi<T>());
}
template<typename T>
T inline simple_lfo<T>::phasor_offset(T offset) const noexcept {
	return nvs::gen::wrap01(_phase + offset);
}
template<typename T>
T simple_lfo<T>::saw_offset(T offset) const noexcept {
	return nvs::memoryless::unibi(phasor_offset(offset));
}
template<typename T>
T simple_lfo<T>::square_offset(T offset) const noexcept {
	return phasor_offset(offset) > 0.5f ? 1.f : -1.f;
}
template<typename T>
T simple_lfo<T>::tri_offset(T offset) const noexcept {
	T phase = phasor_offset(offset);
	return nvs::memoryless::unibi(phase > 0.5f ? (phase - 0.5) * 2 : (1.f - phase - 0.5) * 2);
}
template<typename T>
T simple_lfo<T>::sine_offset(T offset) const noexcept {
	using namespace nvs::memoryless;
	return math_impl::cos(phasor_offset(offset) * math_impl::two_pi<T>());
}
// takes in value [0 ... 3] and returns interpolated waveform
// 0 = tri
// 1 = saw
// 2 = square
// 3 = sine
template<typename T>
T simple_lfo<T>::multi(T waveform) const noexcept {
	using namespace nvs::memoryless;
	waveform = clamp<T>(waveform, 0, 3);
	int selector = (int)(floor(waveform));
	T interpolationAmt = nvs::gen::wrap01(waveform);
	
	auto const retval = [selector, interpolationAmt, this](){
		switch (selector) {
			case 0: { return linterp<T>(tri(), saw(), interpolationAmt); }
			case 1: { return linterp<T>(saw(), square(), interpolationAmt); }
			case 2: { return linterp<T>(square(), sine(), interpolationAmt); }
			case 3: { return sine(); }
			default: { return 0.f; }
		}
	}();
	assert (retval == retval);
	return retval;
}
template<typename T>
void simple_lfo<T>::reset() noexcept {
	_phase = 0.0;
}
template<typename T>
[[deprecated("Pretty much guaranteed to be buggy in some way; use nvs::gen::ramp2trig instead.")]]
bool simple_lfo<T>::crossedOver() const noexcept {
	return (_lastPhase > _phase);
}

}	// namespace nvs::lfo
