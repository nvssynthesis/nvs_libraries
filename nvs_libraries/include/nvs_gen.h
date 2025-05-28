/*
  ==============================================================================

		nvs_gen.h
    Created: 22 Jun 2023 12:18:02am
    Author:  Nicholas Solem

  ==============================================================================
*/

#pragma once
#include <array>
#include <cmath>
#include <concepts>
#include "nvs_memoryless.h"

namespace nvs {
namespace gen {

template<std::floating_point float_t>
struct phasor {
private:
	float_t phase { 0.f };
	float_t phaseDelta {0.f};	// frequency / samplerate
	float_t _sampleRate {0.f};
public:
	inline void setSampleRate(float_t sr){
		assert (sr > 0.f);
		_sampleRate = sr;
	}
	inline void setPhase(float_t ph){
		phase = ph;
	}
	inline float_t getPhase() const {
		return phase;
	}
	inline void reset(){
		phase = 0.f;
	}
	inline void setPhaseDelta(float_t pd){
		phaseDelta = pd;
	}
	// may be called every sample, so no check for aliasing or divide by 0
	inline void setFrequency(float_t frequency){
		assert(_sampleRate > static_cast<float_t>(0.f));
		phaseDelta = frequency / _sampleRate;
	}
	inline float_t getFrequency() const {
		return phaseDelta * _sampleRate;
	}
	// prefix increment
	phasor& operator++(){
		assert(phaseDelta == phaseDelta);
		phase += phaseDelta;
		while (phase >= 1.f)
			phase -= 1.f;
		while (phase < 0.f)
			phase += 1.f;
		return *this;
	}
	// postfix increment
	phasor operator++(int){
		assert(phaseDelta == phaseDelta);
		phasor old = *this;
		phase += phaseDelta;
		while (phase >= 1.f)
			phase -= 1.f;
		while (phase < 0.f)
			phase += 1.f;
		return old;
	}
};

template <typename T>
struct history {
	T val {0.0};
	T operator()(T x){
		T retVal = val;
		val = x;
		return retVal;
	}
};

template <typename T>
struct accum {
	T val {T(0.0)};
	T operator()(T x){
		val += x;
		return val;
	}
	T operator()(T x, bool shouldReset){
		val += x;
		T retVal = val;
		if (shouldReset){
			reset();
		}
		return retVal;
	}
	void reset(){
		val = static_cast<T>(0);
	}
	void set(float newVal){
		val = newVal;
	}
};

/* the sign of the difference between the current and previous input.
 1 if input is increasing, -1 if decreasing, 0 if unchanging. */
template <typename T>
struct change {
private:
	T histo {static_cast<T>(0)};
public:
	T operator()(T x){
		T val;
		if (x > histo)
			val = 1;//static_cast<T>(1);
		else if (x < histo)
			val = -1;//static_cast<T>(-1);
		else // (x == histo)
			val = 0;//static_cast<T>(0);

		histo = x;
		return val;
	}
};

template <typename T>
struct ramp2trig {
private:
	T histo {static_cast<T>(0)};
	change<T> _change;
public:
	T operator()(T x){
		constexpr T zero = static_cast<T>(0);
		
		T diff = x - histo;
		T sum = x + histo;
		T div = diff / sum;
//		T abs = (div < zero) ? -div : div;
		T abs = std::abs(div);
		T cmp = abs > static_cast<T>(0.5);
		T signOfChange = _change(cmp);
		histo = x;
		return static_cast<T>(signOfChange > zero);
	}
};

/* selects between 2 values based on another (boolean) value.
if true, 1st argument is returned, if false, 2nd argument is returned.
 */
template<typename T>
inline T switcher(bool select, T x1, T x2){
	if (select)
		return x1;
	else
		return x2;
}

/* integral template argument for number of outputs
 lets you choose which output the incoming signal is sent to
 a value of 0 will choose no output
 a value greater than numOuts will go to the last out
 */
template<typename T, size_t numOuts = 1>
std::array<T, numOuts>
gateSelect(size_t outputSelection, T x){
	constexpr T zero = static_cast<T>(0);
	std::array<T, numOuts> outs {zero};
	if (outputSelection == 0)
		return outs;
	outputSelection -= 1;	// now 0 will index 0
	outputSelection = std::min(outputSelection, numOuts - 1);
	
	outs[outputSelection] = x;
	return outs;
}

/* conditionally passes or holds input.
 if control is non-zero, input value is passed.
 otherwise, the previous held input value is passed
 */
template <typename T>
struct latch {
private:
	T histo {static_cast<T>(0)};
public:
	latch(T initial = static_cast<T>(0))
	:	histo(initial){}
	
	T operator()(T x, bool letPass){
		if (letPass){
			histo = x;
			return x;
		} else {
			return histo;
		}
	}
	void set(T x) {	// could use operator(x, false) and ignore return val, but this is more sense-making
		histo = x;
	}
};

template <typename T>
std::array<T, 2> pol2car(T r, T theta){
	std::array<T, 2> cart;
	cart[0] = r * nvs::memoryless::padeCos(theta);
	cart[1] = r * nvs::memoryless::padeSin(theta);
	return cart;
}

template <typename T>
std::array<T, 2> car2pol(T x, T y){
	std::array<T, 2> pol;
	pol[0] = sqrt(x*x + y*y);
	pol[1] = atan2(y, x);
	return pol;
}

template<typename T>
inline T wrap01(T x) {
	T y = x - static_cast<long long int>(x);
	if (x < 0){
		y = static_cast<T>(1) + y;
	}
	return y;
}
template<std::floating_point T>
inline T wrap(T x, T upperLimit) {
	assert(upperLimit > 0);
	x /= upperLimit;
	T y = wrap01(x);
	y *= upperLimit;
	return y;
}
template<std::integral T>
inline T wrap(T x, T upperLimit){
	assert(upperLimit > 0);
	T y = x % upperLimit;
	if (x < 0){
		y += upperLimit;
	}
	return y;
}

//template<typename T>
//T mspWrap(T f){
//	f = (f > std::numeric_limits<int>::max() || f < std::numeric_limits<int>::min()) ? 0. : f;
//	int k = static_cast<int>(f);
//	T val;
//	if (k <= f)
//		val = f-k;
//	else
//		val = f - (k-1);
//	return val;
//}
/*
 potentially assymetrical phasor-to-triangle shaper.
 Unipolar output, [0..1], expecting input also [0..1]; however
 the template 'assumeBoundedInput = false' allows to take care of unbounded inputs.
 If assumeBoundedInput is false, we can either WRAP or CLAMP to bound it, based on
 the wrap boolean.
 
 So keep in mind wrap is only relevant if assumeBoundedInput is false.
 
 Uses the piecewise linear equation:
	  /   x % 1
	 /    -----      if (x % 1) < a
	 |      a
 y = >
	 | 1 - (x % 1)
	 \ -----------   if (x % 1) >= a
	  \   1 - a
 */
#define USING_STD_FMOD 0

template <typename float_t=float, bool assumeBoundedInput=true, bool wrap=false>
float_t triangle(float_t x, float_t skew = 0.5){
	// assume input is bounded 0-1
	if constexpr (!assumeBoundedInput){
		if constexpr (wrap)
		{
#if USING_STD_FMOD
			x = std::fmod(x, 1.f);
			if (x < 0.f)	// built in fmod does not properly take care of negatives
				x *= -1.f;
#else
			x = nvs::memoryless::mspWrap<float_t>(x);
#endif
		} // end wrap
		else {	// clamp
			x = nvs::memoryless::clamp<float_t>(x, 0.f, 1.f);
		}
	}
	skew = nvs::memoryless::clamp<float>(skew, 0.000001, 0.99999);
	if (x < skew)
		return x / skew;
	else
		return (1.f - x) / (1.f - skew);
}

template <typename float_t>
float_t parzen(float_t x){
	x = nvs::memoryless::clamp<float_t>(x, 0.f, 1.f);
	
	float_t x3 = x * x * x;
	float_t oneMinusX = 1.f - x;

	float_t retval = 2.f * x3;
	if (x > 0.5f){
		retval = 1.f - (6.f * oneMinusX * oneMinusX * x);
	}
	return retval;
}

enum class boundsModes_e {
	clamp = 0,
	wrap
};

template <typename T, boundsModes_e b=boundsModes_e::clamp>
T readBuff(T const *data, long index, unsigned long length){
	assert(length > 0);
	if constexpr (b == boundsModes_e::clamp){
		index = nvs::memoryless::clamp<long>(index, 0, length - 1);
	}
	else if constexpr (b == boundsModes_e::wrap){
		index = wrap(index, static_cast<long>(length - 1));
	}
	assert(index >= 0);
	assert(index < length);
	return data[index];
}

enum class interpolationModes_e {
	floor,
	linear,
	hermite
};

// laurent de soras
template<typename float_t>
inline float_t hermite(float_t frac_pos, float_t xm1, float_t x0, float_t x1, float_t x2)
{
   const float_t    c     = (x1 - xm1) * 0.5f;
   const float_t    v     = x0 - x1;
   const float_t    w     = c + v;
   const float_t    a     = w + v + (x2 - x0) * 0.5f;
   const float_t    b_neg = w + a;

   return ((((a * frac_pos) - b_neg) * frac_pos + c) * frac_pos + x0);
}

template <typename T, interpolationModes_e interp=interpolationModes_e::linear, boundsModes_e b=boundsModes_e::clamp>
inline T peek(T const *data, double fracIndex, unsigned long length){
	assert (data != nullptr);
	
	long iidx = static_cast<long>(fracIndex);
	T const x0 = readBuff<T, b>(data, iidx, length);
	
	if constexpr (interp == interpolationModes_e::floor){
		return x0;
	}
	
	double const frac = fracIndex - static_cast<double>(iidx);
	T const x1 = readBuff<T, b>(data, iidx + 1, length);
	
	if constexpr (interp == interpolationModes_e::linear){
		return (1.0 - frac)*x0 + (frac)*x1;
	}
	else if constexpr (interp == interpolationModes_e::hermite){
		T const xm1 = readBuff<T, b>(data, iidx - 1, length);
		T const x2 = readBuff<T, b>(data, iidx + 2, length);
		return hermite(static_cast<T>(frac), xm1, x0, x1, x2);
	}
//	assert(false);
}

#if 0
/*
 Read values from a buffer.
 
Template args:
 -number of output channels, used to determine size of std::array
 -sample data type, normally float
 -container type for each data channel

Construction args:
 -an array of references to buffers (default type vector of floats).
 
Call operator args:
 -sample index to read (no interpolation); indices out of range return zero.
 -channel offset (default 0).
 */
template<size_t numOutputChans, typename float_t = float, typename singleChannelBuffer_t = std::vector<float_t>>
struct peek {
private:
	typedef std::array<singleChannelBuffer_t const *, numOutputChans> multiChannelBuffer_t;
	std::array<singleChannelBuffer_t const *, numOutputChans> buffer;
	
public:
	peek(std::array<singleChannelBuffer_t const *, numOutputChans> ptrsToBuffers)
	:	buffer(ptrsToBuffers){}
	
	void replaceBuffer(std::array<singleChannelBuffer_t const *, numOutputChans> newBuffer) {
		buffer = newBuffer;
	}
	
	template<typename index_t>
	inline index_t clampSample(index_t index, index_t buffLength) const {
		return nvs::memoryless::clamp<index_t>(index, 0, buffLength - 1);
	}
	template<typename index_t>
	inline index_t wrapSample(index_t index, index_t buffLength) const {
		while (index < 0)
			index += buffLength;
		while (index >= buffLength)
			index -= buffLength;
		return index;
	}
	
	template<typename index_t, boundsModes_e boundsMode = boundsModes_e::wrap>
	std::array<float_t, numOutputChans>
	operator()(index_t index, int channel_offset = 0) const {
		std::array<float_t, numOutputChans> channelwiseSamps {};
		
		for (int chan = 0; chan < numOutputChans; ++chan){
			if (buffer[chan] == nullptr) {
				return channelwiseSamps;	// array of zeros
			}
		}
		
		index_t index_0 = static_cast<index_t>(index);
		for (int chan = 0; chan < numOutputChans; ++chan){
			int chan_tmp = nvs::memoryless::clamp<int>(chan + channel_offset, 0, numOutputChans - 1);
			
			singleChannelBuffer_t const& thisChannelBuffer = *(buffer[chan_tmp]);
			
			if (thisChannelBuffer.size() == 0){
				channelwiseSamps[chan] = static_cast<float_t>(0.f);
				continue;
			}
		
			index_t tmp_idx_0;
			if constexpr (boundsMode == boundsModes_e::clamp){
				tmp_idx_0 = clampSample<index_t>(index_0, thisChannelBuffer.size());
			}
			else if constexpr (boundsMode == boundsModes_e::wrap){
				tmp_idx_0 = wrapSample<index_t>(index_0, thisChannelBuffer.size());
			}
			channelwiseSamps[chan] = thisChannelBuffer[tmp_idx_0];
		}
		return channelwiseSamps;
	}
	
	std::array<size_t, numOutputChans>
	getSizes() const {
		std::array<size_t, numOutputChans> sizes;
		for (size_t chan = 0; chan < numOutputChans; ++chan){
			sizes[chan] = (*(buffer[chan])).size();
		}

		return sizes;
	}
};
#endif

}	// namespace gen
}	// namespace nvs

