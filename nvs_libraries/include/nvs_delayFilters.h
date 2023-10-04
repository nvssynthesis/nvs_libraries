#pragma once
#include "nvs_filters.h"
#include "nvs_memoryless.h"

namespace nvs {

constexpr bool is_power_of_2(int v){
	return v && ((v & (v - 1)) == 0);
}

namespace delays {

template<uint32_t _maxDelSize, typename float_t=float>
class Delay {
	static_assert(is_power_of_2(_maxDelSize), "delay size must be power of two");
	static_assert(std::is_floating_point<float_t>::value,
										"float_t type must be floating point");
	using fractionalIdx_t = double;
	using integralIdx_t = decltype(_maxDelSize);
public:
	Delay()	:	delMask(_maxDelSize - 1)
	{
		clear();
		isInitialized = true;
	}
	
	void setSampleRate(double fs)
	{
		fractionalIdx_t ms = (delTimeSamps / sampleRate) * 1000.0;    // maintain former delay time
		sampleRate = (fs > 0.0) ? fs : 44100.0;
		setDelayTimeMS(ms);
	}
	enum class interp {
		floor = 0,
		linear,
		cubic
	};
	inline fractionalIdx_t getEffectiveDelayTimeSamps(){
		return delTimeSamps + 1.0;
	}
	float_t (Delay::*doTheTick)(float_t);
	// could also implement by multiplication of 1 or zero with different outputs
	// copy pointer to function and ensure it's really pointing there properly
	// also consider that this MUST be assigned already on initialized
	float_t tick(float_t input) {
		//return (this->*doTheTick)(input);
		return (this->*doTheTick)(input);
	}
	
	float_t tick_floor(float_t input) {
		rHead = wHead + (integralIdx_t)getEffectiveDelayTimeSamps();
		rHead &= delMask;
		float_t y = _buffer[rHead];
		
		_buffer[wHead] = input;
		wHead--;
		wHead &= delMask;
		return y;
	}
	float_t tick_linear(float_t input) {
		fractionalIdx_t const dTime = getEffectiveDelayTimeSamps();
		integralIdx_t const dTime_floor = static_cast<integralIdx_t>(dTime);
		
		rHead = wHead + dTime_floor;
		rHead &= delMask;
		
		fractionalIdx_t const frac = dTime - static_cast<fractionalIdx_t>(dTime_floor);
		float_t y = _buffer[rHead] * (1.f - frac);
		
		auto rHeadTmp = rHead + 1;
		rHeadTmp &= delMask;
		y += _buffer[rHeadTmp] * (frac);
		
		_buffer[wHead] = input;
		wHead--;
		wHead &= delMask;
		return y;
	}
	float_t tick_cubic(float_t input) {
		rHead = wHead + (int)getEffectiveDelayTimeSamps();
		fractionalIdx_t frac = getEffectiveDelayTimeSamps() -
							(fractionalIdx_t)((int)getEffectiveDelayTimeSamps());
		
		float_t xm1 = _buffer[(rHead - 1) & delMask];
		float_t x0  = _buffer[(rHead + 0) & delMask];
		float_t x1  = _buffer[(rHead + 1) & delMask];
		float_t x2  = _buffer[(rHead + 2) & delMask];
		float_t a = (3.f * (x0-x1) - xm1 + x2) / 2.f;
		float_t b = 2.f*x1 + xm1 - (5.f*x0 + x2) / 2.f;
		float_t c = (x1 - xm1) / 2.f;
		float_t y = (((a * frac) + b) * frac + c) * frac + x0;
		
		_buffer[wHead] = input;
		wHead--;
		wHead &= delMask;
		return y;
	}
	
	void setInterpolation(int _interpType) {
		this->interpType = _interpType;
		switch (_interpType) {
			case 0: // look at address of doTheTick, make sure that's the address that actually gets returned in the tick function
				doTheTick = &Delay::tick_floor;
			case 1:
				doTheTick = &Delay::tick_linear;
			case 2:
				doTheTick = &Delay::tick_cubic;
		}
	}
	// multiple tick functions for different interpolations,
	// then use function pointer to set tick to that function
	
	float tap(unsigned int writeOffset) const {
		integralIdx_t pos = wHead + writeOffset;
		pos &= pos;
		return _buffer[pos];
	}
	void setDelayTimeSamps(double t){
		delTimeSamps = nvs::memoryless::clamp(t, 0.0, (double)_maxDelSize);
	}
	void setDelayTimeMS(double t)
	{
		t = t * 0.001 * sampleRate;
		setDelayTimeSamps(t);
	}
	
	void updateDelayTimeMS(float target, float oneOverBlockSize) {
		target = target * 0.001f * sampleRate; // convert to samples
		float inc = (target - delTimeSamps) * oneOverBlockSize;
		inc = nvs::memoryless::atan_aprox(inc);
		delTimeSamps += inc;
	}
	unsigned int getDelaySize() const { return _maxDelSize; }
	
	void clear() {
		_buffer.fill(0.f);
		/*
		int n = _maxDelSize;
		while (n > 0)
		{
			n -= 1;     // pre-decrement
			_buffer[n] = 0.f;
		}
		 */
	}
	[[deprecated("unnecessary, the object is initialized upon construction")]]
	bool ready() {
		return isInitialized;
	}
private:
	fractionalIdx_t delTimeSamps, sampleRate;
	std::array<float_t, _maxDelSize> _buffer;
	integralIdx_t delMask;
	integralIdx_t wHead {0}, rHead;
	short int interpType;
	[[deprecated("unnecessary, the object is initialized upon construction")]]
	bool isInitialized {false};
};

#if IMPLEMENTED
class AllpassDelay : public Delay {
public:
	AllpassDelay() :   AllpassDelay(8192, 44100.f) {}
	AllpassDelay(unsigned int _delSize, float sample_rate)
	:  Delay(_delSize, sample_rate), g(0.f), z(0.f){}
	
	void update_g(float g_target, float oneOverBlockSize)
	{
		g += (g_target - g) * oneOverBlockSize;
	}
	void set_g(float g_target)
	{
		g = g_target;
	}
	float filter(float x_n)
	{
		float u_n, y_n;
		u_n = x_n + -g * z;
		y_n = z + g * u_n;
		z = tick_cubic(u_n);
		return y_n;
	}
	
private:
	float g, z;
};
#endif

}	// namespace delay
}	// namespace nvs
