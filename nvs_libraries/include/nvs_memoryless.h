/*
  ==============================================================================

    nvs_dsp.h
    Created: 5 Dec 2018 9:24:03pm
    Author:  Nicholas Solem

  ==============================================================================
*/

#pragma once

#include <array>
#include <cmath>
#include <type_traits>
#include <concepts>
//#define USING_SPROUT 0

#ifndef USING_SPROUT
#define USING_SPROUT 1
#endif

#if USING_SPROUT
#include "sprout/math.hpp"
#include "sprout/math/constants.hpp"
#else
#include <numbers>
#endif

namespace nvs::memoryless {

template<typename T>
concept FloatingPoint = std::floating_point<T>;

namespace math_impl {
	template<FloatingPoint T>
	constexpr T exp(T x) noexcept {
#if USING_SPROUT
		return sprout::exp(x);
#else
		return std::exp(x);
#endif
	}

	template<FloatingPoint T>
	constexpr T tan(T x) noexcept {
#if USING_SPROUT
		return sprout::tan(x);
#else
		return std::tan(x);
#endif
	}

	template<FloatingPoint T>
constexpr T atan(T x) noexcept {
#if USING_SPROUT
		return sprout::atan(x);
#else
		return std::atan(x);
#endif
	}

	template<FloatingPoint T>
	constexpr T tanh(T x) noexcept {
#if USING_SPROUT
		return sprout::tanh(x);
#else
		return std::tanh(x);
#endif
	}

	template<FloatingPoint T>
	constexpr T sin(T x) noexcept {
#if USING_SPROUT
		return sprout::sin(x);
#else
		return std::sin(x);
#endif
	}

	template<FloatingPoint T>
	constexpr T cos(T x) noexcept {
#if USING_SPROUT
		return sprout::cos(x);
#else
		return std::cos(x);
#endif
	}

	// Constants with sprout fallback
	template<FloatingPoint T>
	constexpr T pi() noexcept {
#if USING_SPROUT
		return sprout::math::pi<T>();
#else
		return std::numbers::pi_v<T>;
#endif
	}
	template<FloatingPoint T>
	constexpr T halfpi() noexcept {
#if USING_SPROUT
		return sprout::math::pi<T>() * 0.5;
#else
		return std::numbers::pi_v<T> * 0.5;
#endif
	}

	template<FloatingPoint T>
	constexpr T two_pi() noexcept {
#if USING_SPROUT
		return sprout::math::two_pi<T>();
#else
		return 2.0 * std::numbers::pi_v<T>;
#endif
	}
}

// linear interpolate.
template<FloatingPoint T>
constexpr T linterp(T a, T b, T t) noexcept {
	assert ((0.f <= t) && (t <= 1.f)); // t only makes sense in [0..1] range
	return a + t * (b - a);
}

// unipolar to bipolar (assuming input is [0..1]
template<FloatingPoint T>
T unibi(T x) noexcept {
	assert ((0.f <= x) && (x <= 1.f));
    return (x * 2.f) - 1.f;
}
// bipolar to unipolar (assuming input is [-1..1]
template<FloatingPoint T>
T biuni(T x) noexcept {
	assert ((-1.f <= x) && (x <= 1.f));
    return (x + 1.f) * 0.5f;
}

template<FloatingPoint T>
T clamp1(T x) noexcept {
    if (x > 1.f)
        return 1.f;
    else if (x < -1.f)
        return -1.f;
    else
        return x;
}

template<class T> concept Numeric = std::is_arithmetic_v<T>;
template<Numeric T>
T clamp(T x, T min, T max) noexcept {
	return (x < min) ? min : ((x > max) ? max : x);
}
template<Numeric T>
T clamp_low(T x, T min) noexcept {
    return (x < min) ? min : x;
}
template<Numeric T>
T clamp_high(T x, T max) noexcept {
    return (x > max) ? max : x;
}

// based on https://stackoverflow.com/questions/14415753/wrap-value-into-range-min-max-without-division
template<FloatingPoint T>
T div_wrap(T x, T xmin, T xmax)
{
    return x - (xmax - xmin) * floor( x / (xmax - xmin));
}

template<FloatingPoint T>
T wrap(T x, T xMin, T xMax)
{
    return (T)fmod((fmod(x - xMin, xMax - xMin) + xMax - xMin) , xMax - xMin) + xMin;
}

template<FloatingPoint T>
[[deprecated("Non-constant time approach")]]
T mod_1(T input) {
    while (input >= 1.f)
        input -= 1.f;
    while (input < 0.f)
        input += 1.f;
    return input;
}

template<FloatingPoint T>
[[deprecated("Need to use replacement for logf and sqrt")]]
T unboundSat1(T x) {
	return logf(x + sqrt(x*x + 1.f));
}
template<FloatingPoint T>
[[deprecated("Replace sqrt and fabs")]]
T unboundSat2(T x) {
	T num = 2.f * x;
	T denom = 1.f + sqrt(1.f + fabs(4.f * x));
	return num / denom;
}

template<FloatingPoint T>
T crush(T input, T depth) noexcept {
    //T depth = pow(2, exponent);
    T crushed = int(input * depth) / depth;
    return crushed;
}

template<FloatingPoint T>
T xOverOnePlusAbsX(T input) noexcept {
    // y = x/(1 + |x|)
    return input / (1 + (input > 0 ? input : -input));
}

template<FloatingPoint T>
T mspWrap(T f) noexcept
{
	f = (f > std::numeric_limits<int>::max() || f < std::numeric_limits<int>::min()) ? 0. : f;
	int k = static_cast<int>(f);
	T val;
	if (k <= f)
		val = f-k;
	else
		val = f - (k-1);
	return val;
}

template<FloatingPoint T>
T scale(T val, T min, T range) noexcept {
	return (val - min) / range;
}
template<FloatingPoint T>
size_t round(T x) noexcept {
	return x >= 0 ? static_cast<size_t>(x + 0.5) : static_cast<size_t>(x - 0.5);
}

//========================================LUT Stuff=======================================

enum class IndexBoundaryPolicy {
	Clamp,
	Wrap
};

enum class InterpolationType {
	Rounded,
	Linear
};

// base lookup table using CRTP
template<size_t N>
concept PowerOfTwo = N > 0 && (N & (N - 1)) == 0;

template<typename T>
concept HasRange = requires {
	{ T::min } -> std::convertible_to<double>;
	{ T::max } -> std::convertible_to<double>;
} && (T::min < T::max);

template<typename Derived, typename sample_t, typename index_t, HasRange RangeType, size_t Resolution,
	IndexBoundaryPolicy BoundaryPolicy, InterpolationType DefaultInterpType>
requires PowerOfTwo<Resolution> && FloatingPoint<sample_t> && FloatingPoint<index_t>
class LookupTableBase {
public:
	virtual ~LookupTableBase() = default;

	[[nodiscard]]
	virtual constexpr sample_t operator()(index_t x) const noexcept
	{
		if constexpr (DefaultInterpType == InterpolationType::Rounded) { return rounded(x); }
		else if constexpr (DefaultInterpType == InterpolationType::Linear) { return linear(x); }
		else { assert(false); } // No such interpolation type
		return 0.f;
	}
	[[nodiscard]]
	virtual constexpr sample_t rounded(index_t x) const noexcept {
		index_t const normalized = scale(x, static_cast<index_t>(min_x), static_cast<index_t>(x_range));
		index_t const bound_idx = constrainIndexBy01(normalized);
		index_t const continuous_index = [bound_idx]() {
			if constexpr (BoundaryPolicy == IndexBoundaryPolicy::Wrap) {
				return bound_idx * static_cast<index_t>(Resolution);
			}
			else if constexpr (BoundaryPolicy == IndexBoundaryPolicy::Clamp) {
				return bound_idx * static_cast<index_t>(Resolution - 1);
			}
		}();

		// Round to nearest index (handles negative values correctly)
		auto const rounded_index = constrainIndexByReso (nvs::memoryless::round(continuous_index));

		return rounded_lookup(rounded_index);
	}
	[[nodiscard]]
	virtual constexpr sample_t linear(index_t x) const noexcept {
		const index_t normalized = scale(x, static_cast<index_t>(min_x), static_cast<index_t>(x_range));
		return interpolated_lookup(normalized);
	}
	using value_type = sample_t;        // What the table returns
	using index_type = index_t;         // What you pass in for lookup
	using size_type = size_t;           // For consistency with STL

protected:
	static constexpr auto constrainIndexBy01(std::floating_point auto fpIdx) {
		if constexpr (BoundaryPolicy == IndexBoundaryPolicy::Wrap) {
			return mspWrap(fpIdx);
		}
		else if constexpr (BoundaryPolicy == IndexBoundaryPolicy::Clamp) {
			return clamp1(fpIdx);
		}
		else { assert(false); return 0; } //No such boundary policy.
	}
	static constexpr size_t constrainIndexByReso(const size_t idx) {
		if constexpr (BoundaryPolicy == IndexBoundaryPolicy::Wrap) {
			return idx % Resolution;
		}
		else if constexpr (BoundaryPolicy == IndexBoundaryPolicy::Clamp) {
			return clamp_high(idx, Resolution - 1);
		}
		else { assert(false); return 0; } // No such boundary policy.

	}
	constexpr sample_t rounded_lookup(size_t index) const noexcept {
		assert(index < values_.size());
		return values_[index];
	}

	constexpr sample_t interpolated_lookup(index_t normalized_index) const noexcept
	{
		const index_t scaled_index = [normalized_index]() {
			if constexpr (BoundaryPolicy == IndexBoundaryPolicy::Wrap) {
				return normalized_index * static_cast<index_t>(Resolution);
			}
			else if constexpr (BoundaryPolicy == IndexBoundaryPolicy::Clamp) {
				return normalized_index * static_cast<index_t>(Resolution - 1);
			}
		}();

		const auto floored_index = static_cast<size_t>(scaled_index);
		assert (floored_index < values_.size());

		const sample_t frac = scaled_index - static_cast<sample_t>(floored_index);
		const auto next_index = constrainIndexByReso(floored_index + 1);
		return linterp(rounded_lookup(floored_index), rounded_lookup(next_index), frac);
	}
	std::array<sample_t, Resolution> values_{};
	static constexpr double min_x = RangeType::min;
	static constexpr double max_x = RangeType::max;
	static constexpr double x_range = max_x - min_x;
};

template<std::integral auto Min, std::integral auto Max>
requires std::is_same_v<decltype(Min), decltype(Max)>
struct IRange {
	using value_type = decltype(Min);
	static constexpr auto min = Min;
	static constexpr auto max = Max;

	static_assert(Min < Max, "Min must be less than Max");
	static_assert(std::same_as<decltype(Min), decltype(Max)>,
				 "Min and Max must be the same type");
};
template<typename sample_t, typename index_t, HasRange RangeType, size_t Resolution,
	InterpolationType DefaultInterpType>
class ExpTable final : public LookupTableBase<ExpTable<sample_t, index_t, RangeType, Resolution, DefaultInterpType>,
	sample_t, index_t, RangeType, Resolution, IndexBoundaryPolicy::Clamp, DefaultInterpType>
{
	using Base = LookupTableBase<ExpTable, sample_t, index_t, RangeType, Resolution,
		IndexBoundaryPolicy::Clamp, DefaultInterpType>;
	static_assert(std::is_same_v<typename Base::value_type, sample_t>);
	static_assert(std::is_same_v<typename Base::index_type, index_t>);

public:
	constexpr ExpTable() noexcept
	{
		const double increment = Base::x_range / static_cast<double>(Resolution - 1);

		for (size_t i = 0; i < Resolution; ++i) {
			const double x = static_cast<double>(Base::min_x) + (double)i * increment;
			this->values_[i] = static_cast<sample_t>(math_impl::exp(x));
		}
	}
};

template<typename sample_t, typename index_t, HasRange RangeType, size_t Resolution,
InterpolationType DefaultInterpType>
class AtanTable final : public LookupTableBase<AtanTable<sample_t, index_t, RangeType, Resolution, DefaultInterpType>,
	sample_t, index_t, RangeType, Resolution, IndexBoundaryPolicy::Clamp, DefaultInterpType>
{
	using Base = LookupTableBase<AtanTable, sample_t, index_t, RangeType, Resolution,
		IndexBoundaryPolicy::Clamp, DefaultInterpType>;
	static_assert(std::is_same_v<typename Base::value_type, sample_t>);
	static_assert(std::is_same_v<typename Base::index_type, index_t>);

public:
	constexpr AtanTable() noexcept
	{
		const double increment = Base::x_range / static_cast<double>(Resolution - 1);

		for (size_t i = 0; i < Resolution; ++i) {
			const double x = static_cast<double>(Base::min_x) + (double)i * increment;
			this->values_[i] = static_cast<sample_t>(math_impl::atan(x));
		}
	}
};

template<typename sample_t, typename index_t, HasRange RangeType, size_t Resolution,
	InterpolationType DefaultInterpType>
class TanhTable final : public LookupTableBase<TanhTable<sample_t, index_t, RangeType, Resolution, DefaultInterpType>,
		sample_t, index_t, RangeType, Resolution, IndexBoundaryPolicy::Clamp, DefaultInterpType>
{
	using Base = LookupTableBase<TanhTable, sample_t, index_t, RangeType, Resolution,
		IndexBoundaryPolicy::Clamp, DefaultInterpType>;
	static_assert(std::is_same_v<typename Base::value_type, sample_t>);
	static_assert(std::is_same_v<typename Base::index_type, index_t>);

public:
	constexpr TanhTable() noexcept
	{
		const double increment = Base::x_range / static_cast<double>(Resolution - 1);

		for (size_t i = 0; i < Resolution; ++i) {
			const double x = static_cast<double>(Base::min_x) + (double)i * increment;
			this->values_[i] = static_cast<sample_t>(math_impl::tanh(x));
		}
	}
};

struct CosSinRange {
	constexpr static double min {0.0};
	constexpr static double max {math_impl::two_pi<double>()};
};

template<typename sample_t, typename index_t, size_t Resolution, InterpolationType DefaultInterpType>
class CosTable : public LookupTableBase<CosTable<sample_t, index_t, Resolution, DefaultInterpType>,
	sample_t, index_t, CosSinRange, Resolution,
	IndexBoundaryPolicy::Wrap, DefaultInterpType>
{
	using Base = LookupTableBase<CosTable, sample_t, index_t, CosSinRange, Resolution,
		IndexBoundaryPolicy::Wrap, DefaultInterpType>;
	static_assert(std::is_same_v<typename Base::value_type, sample_t>);
	static_assert(std::is_same_v<typename Base::index_type, index_t>);

	static constexpr double two_pi { math_impl::two_pi<double>() };
	static constexpr double inv_two_pi { 1.0 / math_impl::two_pi<double>() };

public:
	static constexpr size_t mask = Resolution - 1; // For power-of-2 wraparound
	constexpr CosTable() noexcept
	{
		const double increment = Base::max_x / static_cast<double>(Resolution);

		for (size_t i = 0; i < Resolution; ++i) {
			const double phase = static_cast<double>(i) * increment;
			this->values_[i] = static_cast<sample_t>(math_impl::cos(phase));
		}
	}
	~CosTable() override = default;

	[[nodiscard]]
	constexpr sample_t operator()(index_t phase) const noexcept override
	{
		// phase = Base::constrainIndexByReso(phase);
		phase = mspWrap(phase * inv_two_pi) * two_pi;
		return Base::operator()(phase);
	}
	[[nodiscard]]
	constexpr sample_t cos(index_t phase) const noexcept
	{
		return operator()(phase);
	}
	[[nodiscard]]
	constexpr sample_t cosRounded(index_t phase) const noexcept {
		phase = constrainIndexByReso(phase);
		return Base::rounded(phase);
	}
	[[nodiscard]]
	constexpr sample_t cosLinear(index_t phase) const noexcept
	{
		phase = mspWrap(phase * inv_two_pi) * two_pi;
		return Base::linear(phase);
	}
	[[nodiscard]]
	inline constexpr sample_t linear(index_t phase) const noexcept override
	{
		return cosLinear(phase);
	}
	[[nodiscard]]
	constexpr sample_t sin(index_t phase) const noexcept
	{
		phase = phase - math_impl::halfpi<index_t>();
		phase = mspWrap(phase * inv_two_pi) * two_pi;
		return Base::operator()(phase);
	}
	[[nodiscard]]
	constexpr sample_t sinRounded(index_t phase) const noexcept
	{
		phase = phase - math_impl::halfpi<index_t>();
		phase = mspWrap(phase * inv_two_pi) * two_pi;
		return Base::rounded(phase);
	}
	[[nodiscard]]
	constexpr sample_t sinLinear(index_t phase) const noexcept
	{
		phase = phase - math_impl::halfpi<index_t>();
		phase = mspWrap(phase * inv_two_pi) * two_pi;
		return Base::linear(phase);
	}
};

struct TanRange {
	constexpr static double min {0.0};
	constexpr static double max {math_impl::pi<double>()};
};
template<typename sample_t, typename index_t, size_t Resolution, InterpolationType DefaultInterpType>
class TanTable final : public LookupTableBase<TanTable<sample_t, index_t, Resolution, DefaultInterpType>,
	sample_t, index_t, TanRange, Resolution,
	IndexBoundaryPolicy::Wrap, DefaultInterpType>
{
	using Base = LookupTableBase<TanTable, sample_t, index_t, TanRange, Resolution,
		IndexBoundaryPolicy::Wrap, DefaultInterpType>;
	static_assert(std::is_same_v<typename Base::value_type, sample_t>);
	static_assert(std::is_same_v<typename Base::index_type, index_t>);

	static constexpr double pi { math_impl::pi<double>() };
	static constexpr double inv_pi { 1.0 / math_impl::pi<double>() };

public:
	static constexpr size_t mask = Resolution - 1; // For power-of-2 wraparound
	constexpr TanTable() noexcept
	{
		const double increment = Base::max_x / static_cast<double>(Resolution);

		for (size_t i = 0; i < Resolution; ++i) {
			const double phase = static_cast<double>(i) * increment;
			this->values_[i] = static_cast<sample_t>(math_impl::tan(phase));
		}
	}
	~TanTable() override = default;

	[[nodiscard]]
	constexpr sample_t operator()(index_t phase) const noexcept override
	{
		// phase = Base::constrainIndexByReso(phase);
		phase = mspWrap(phase * inv_pi) * pi;
		return Base::operator()(phase);
	}
	[[nodiscard]]
	constexpr sample_t rounded(index_t phase) const noexcept override
	{
		phase = mspWrap(phase * inv_pi) * pi;
		return Base::rounded(phase);
	}
	[[nodiscard]]
	inline constexpr sample_t linear(index_t phase) const noexcept override
	{
		phase = mspWrap(phase * inv_pi) * pi;
		return Base::linear(phase);
	}
};

template<FloatingPoint T>
T padeSin(T x)
{
    T xx = x*x;
    T x3 = xx * x;
    T x5 = x3 * xx;
    T num = (551.f * x5)/166320.f - (53.f * x3)/396.f + x;
    T den = (5.f * xx * xx)/11088.f + (13.f * xx)/396.f + 1.f;
    return num / den;
}
template <FloatingPoint T>
T padeCos(T x)
{
    T xx = x * x;
    T x4 = xx * xx;
    T num = (313.f * x4)/15120.f - (115.f * xx)/252.f + 1.f;
    T den = (13.f * x4)/15120.f + (11.f * xx)/252.f + 1.f;
    return num / den;
}

// Taylor approximation of sine
// sin(x) =~ x - x^3/3! + x^5/5! - x^7/7!
template<FloatingPoint T>
T taySin(T input)
{
    T x = input * (math_impl::two_pi) - math_impl::pi;
    return x - (pow(x, 3) * 0.16666) + (pow(x, 5) * 0.008333) - (pow(x, 7) * 0.0001984);
    //+ (pow(x,9)*0.0000027557); //(use if you want more accuracy)
}

// Taylor approximation of cosine
//cos(x) = 1 − x^2/2 + x^4/24 ...
template<FloatingPoint T>
T tayCos(T input)
{
    T x = input * (math_impl::two_pi) - math_impl::pi;
    return 1 - pow(x, 2) * 0.5 + pow(x, 4) * 0.0416666 - pow(x, 6) * 0.0013888 + pow(x,8) * 0.0000248;
}

// linear interpolate from sine to cosine. used in some waveshaping stuff.
template<FloatingPoint T>
T sin2cos(T x, T phi)
{
    return linterp(taySin(x), tayCos(x), phi);
}

// more of an equal-power interpolation from sine to cosine. used in some waveshaping stuff.
template<FloatingPoint T>
T sin2cos_ep(T x, T phi)
{
    return (sqrt(1 - phi)) * taySin(x) + (sqrt(phi)) * tayCos(x);
}

template<FloatingPoint T>
T atan_aprox(T x)
{
    return x / (1 + (0.28125 * x) * (0.28125 * x));
}

}   // namespace nvs::memoryless
