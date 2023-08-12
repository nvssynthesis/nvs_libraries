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
#include "nvs_libraries/include/nvs_memoryless.h"

namespace nvs {
namespace gen {

struct phasor {
private:
	float phase { 0.f };
	float phaseDelta {0.f};	// frequency / samplerate
	float const &_sampleRate;
public:
	explicit phasor(float const &sampleRate)	:	_sampleRate(sampleRate){}

	inline void setPhase(float ph){
		phase = ph;
	}
	inline float getPhase() const {
		return phase;
	}
	inline void reset(){
		phase = 0.f;
	}
	inline void setPhaseDelta(float pd){
		phaseDelta = pd;
	}
	// may be called every sample, so no check for aliasing or divide by 0
	inline void setFrequency(float frequency){
		phaseDelta = frequency / _sampleRate;
	}
	inline float getFrequency() const {
		return phaseDelta * _sampleRate;
	}
	// prefix increment
	phasor& operator++(){
		phase += phaseDelta;
		while (phase >= 1.f)
			phase -= 1.f;
		while (phase < 0.f)
			phase += 1.f;
		return *this;
	}
	// postfix increment
	phasor operator++(int){
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
	T val;
	T operator()(T x){
		T retVal = val;
		val = x;
		return retVal;
	}
};

template <typename T>
struct accum {
	T val;
	T operator()(T x){
		val += x;
		return val;
	}
	T operator()(T x, bool shouldReset){
		val += x;
		T retVal = val;
		reset(shouldReset);
		return retVal;
	}
	void reset(bool b = true){
		if (b){
			val = static_cast<T>(0);
		}
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
		else if (x == histo)
			val = 0;//static_cast<T>(0);
		else
			assert(false);

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

/* integral argument for number of outputs
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
};


template <typename T>
std::array<T, 2> pol2car(T r, T theta){
	std::array<T, 2> cart;
	cart[0] = r * cos(theta);
	cart[1] = r * sin(theta);
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
inline T wrap(T x, T upperLimit) {
	while (x < static_cast<T>(0))
		x += upperLimit;
	while (x >= upperLimit)
		x -= upperLimit;
	return x;
}

template<typename T>
T mspWrap(T f){
	f = (f > std::numeric_limits<int>::max() || f < std::numeric_limits<int>::min()) ? 0. : f;
	int k = static_cast<int>(f);
	T val;
	if (k <= f)
		val = f-k;
	else
		val = f - (k-1);
	return val;
}
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
		if constexpr (wrap) {
#if USING_STD_FMOD
			x = std::fmod(x, 1.f);
			if (x < 0.f)	// built in fmod does not properly take care of negatives
				x *= -1.f;
#else
			x = mspWrap<float>(x);
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
	
	float_t selector = (x > 0.5f);

	float_t x3 = x * x * x;
	float_t mul5 = 2.f * x3;
	
	float_t rsub = 1.f - x;
	float_t mul_7 = 6.f * rsub * rsub * x;
	float_t rsub_10 = 1.f - mul_7;

	return switcher<float>(static_cast<bool>(selector), rsub_10, mul5);
}



enum class boundsModes_e {
	clamp = 0,
	wrap
};

template <typename T, boundsModes_e b = boundsModes_e::clamp>
T peekBuff(T const *data, long index, unsigned long length){
	if (length < 1){
		return 0.f;
	}
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

}	// namespace gen
}	// namespace nvs

