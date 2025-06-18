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

	template<FloatingPoint T>
	constexpr T tan(T x) noexcept {
#if USING_SPROUT
		return sprout::tan(x);
#else
		return std::tan(x);
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

template<typename Derived, typename float_t, HasRange RangeType, size_t Resolution,
	IndexBoundaryPolicy BoundaryPolicy, InterpolationType DefaultInterpType>
requires PowerOfTwo<Resolution> && FloatingPoint<float_t>
class LookupTableBase {
public:
	virtual ~LookupTableBase() = default;

	[[nodiscard]]
	virtual constexpr float_t operator()(float_t x) const noexcept
	{
		if constexpr (DefaultInterpType == InterpolationType::Rounded) { return rounded(x); }
		else if constexpr (DefaultInterpType == InterpolationType::Linear) { return linear(x); }
		else { assert(false); } // No such interpolation type
		return 0.f;
	}
	[[nodiscard]]
	virtual constexpr float operator()(double x) const noexcept
	{
		if constexpr (!std::is_same_v<float_t, float>) {
			static_assert(std::is_same_v<float_t, float>, "This function should only be used when the table is single precision, but we want to index with double-precision.");
			return 0.f;
		}
		if constexpr (DefaultInterpType == InterpolationType::Rounded) { return rounded(x); }
		else if constexpr (DefaultInterpType == InterpolationType::Linear) { return linear(x); }
		else { assert(false); } // No such interpolation type.
		return 0.f;
	}
	[[nodiscard]]
	virtual constexpr float_t rounded(float_t x) const noexcept {
		float_t const normalized = scale(x, static_cast<float_t>(min_x), static_cast<float_t>(x_range));
		float_t const bound_idx = constrainIndexBy01(normalized);
		float_t const continuous_index = [bound_idx]() {
			if constexpr (BoundaryPolicy == IndexBoundaryPolicy::Wrap) {
				return bound_idx * static_cast<float_t>(Resolution);
			}
			else if constexpr (BoundaryPolicy == IndexBoundaryPolicy::Clamp) {
				return bound_idx * static_cast<float_t>(Resolution - 1);
			}
		}();

		// Round to nearest index (handles negative values correctly)
		auto const rounded_index = constrainIndexByReso (nvs::memoryless::round(continuous_index));

		return rounded_lookup(rounded_index);
	}
	[[nodiscard]]
	virtual constexpr float rounded(double x) const noexcept
	{
		if constexpr (!std::is_same_v<float_t, float>) {
			static_assert(std::is_same_v<float_t, float>, "This function should only be used when the table is single precision, but we want to index with double-precision.");
			return 0.f;
		}
		double const normalized = scale(x, min_x, x_range);
		double const bound_idx = constrainIndexBy01(normalized);
		const float continuous_index = [bound_idx]() {
			if constexpr (BoundaryPolicy == IndexBoundaryPolicy::Wrap) {
				return bound_idx * static_cast<double>(Resolution);
			}
			else if constexpr (BoundaryPolicy == IndexBoundaryPolicy::Clamp) {
				return bound_idx * static_cast<double>(Resolution - 1);
			}
		}();
		auto const rounded_index = nvs::memoryless::round(continuous_index);

		return rounded_lookup(rounded_index);
	}
	[[nodiscard]]
	virtual constexpr float_t linear(float_t x) const noexcept {
		const float_t normalized = scale(x, static_cast<float_t>(min_x), static_cast<float_t>(x_range));
		return interpolated_lookup(normalized);
	}
	// this is simply a specialization for when the table type is float but we want double-precision indexing
	[[nodiscard]]
	virtual constexpr float linear(double x) const noexcept
	{
		if constexpr (!std::is_same_v<float_t, float>) {
			static_assert(std::is_same_v<float_t, float>, "This function should only be used when the table is single precision, but we want to index with double-precision.");
			return 0.f;
		}
		const float_t normalized = scale(x, min_x, x_range);
		return interpolated_lookup(normalized);
	}
	using value_type = float_t;
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
	constexpr float_t rounded_lookup(size_t index) const noexcept {
		assert(index < values_.size());
		return values_[index];
	}

	constexpr float_t interpolated_lookup(auto normalized_index) const noexcept
	{
		assert ((0.0 <= normalized_index) && (normalized_index <= 1.0));

		const float scaled_index = [normalized_index]() {
			if constexpr (BoundaryPolicy == IndexBoundaryPolicy::Wrap) {
				return normalized_index * static_cast<float>(Resolution);
			}
			else if constexpr (BoundaryPolicy == IndexBoundaryPolicy::Clamp) {
				return normalized_index * static_cast<float>(Resolution - 1);
			}
		}();

		const auto floored_index = static_cast<size_t>(scaled_index);
		assert (floored_index < values_.size());
		const float frac = scaled_index - static_cast<float>(floored_index);

		const auto next_index = constrainIndexByReso(floored_index + 1);
		return linterp(rounded_lookup(floored_index), rounded_lookup(next_index), frac);
	}
	std::array<float_t, Resolution> values_{};
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
template<typename float_t, HasRange RangeType, size_t Resolution,
	InterpolationType DefaultInterpType>
class ExpTable final : public LookupTableBase<ExpTable<float_t, RangeType, Resolution, DefaultInterpType>,
	float_t, RangeType, Resolution, IndexBoundaryPolicy::Clamp, DefaultInterpType>
{
	using Base = LookupTableBase<ExpTable, float, RangeType, Resolution,
		IndexBoundaryPolicy::Clamp, DefaultInterpType>;

public:
	constexpr ExpTable() noexcept
	{
		const double increment = Base::x_range / static_cast<double>(Resolution - 1);

		for (size_t i = 0; i < Resolution; ++i) {
			const double x = static_cast<double>(Base::min_x) + (double)i * increment;
			this->values_[i] = static_cast<float_t>(math_impl::exp(x));
		}
	}
};
struct CosSinRange {
	constexpr static double min {0.0};
	constexpr static double max {math_impl::two_pi<double>()};
};

template<typename float_t, size_t Resolution, InterpolationType DefaultInterpType>
class CosTable : public LookupTableBase<CosTable<float_t, Resolution, DefaultInterpType>,
	float_t, CosSinRange, Resolution,
	IndexBoundaryPolicy::Wrap, DefaultInterpType>
{
	using Base = LookupTableBase<CosTable, float_t, CosSinRange, Resolution,
		IndexBoundaryPolicy::Wrap, DefaultInterpType>;
	static_assert(std::is_same_v<typename Base::value_type, float_t>);

	static constexpr double two_pi { math_impl::two_pi<double>() };
	static constexpr double inv_two_pi { 1.0 / math_impl::two_pi<double>() };

public:
	static constexpr size_t mask = Resolution - 1; // For power-of-2 wraparound
	constexpr CosTable() noexcept
	{
		const double increment = Base::max_x / static_cast<double>(Resolution);

		for (size_t i = 0; i < Resolution; ++i) {
			const double phase = static_cast<double>(i) * increment;
			this->values_[i] = static_cast<float_t>(math_impl::cos(phase));
		}
	}
	~CosTable() override = default;

	[[nodiscard]]
	constexpr float_t operator()(float_t phase) const noexcept override
	{
		// phase = Base::constrainIndexByReso(phase);
		phase = mspWrap(phase * inv_two_pi) * two_pi;
		return Base::operator()(phase);
	}
	[[nodiscard]]
	constexpr float operator()(double phase) const noexcept override
	{
		phase = mspWrap(phase * inv_two_pi) * two_pi;
		return Base::operator()(phase);
	}
	[[nodiscard]]
	constexpr float_t cos(float_t phase) const noexcept
	{
		return operator()(phase);
	}
	[[nodiscard]]
	constexpr float cos(double phase) const noexcept
	{
		if constexpr (!std::is_same_v<float_t, float>) {
			static_assert(std::is_same_v<float_t, float>, "This function should only be used when the table is single precision, but we want to index with double-precision.");
			return 0.f;
		}
		phase = mspWrap(phase * inv_two_pi) * two_pi;
		return operator()(phase);
	}
	[[nodiscard]]
	constexpr float_t cosRounded(float_t phase) const noexcept {
		phase = constrainIndexByReso(phase);
		return Base::rounded(phase);
	}
	[[nodiscard]]
	constexpr float cosRounded(double phase) const noexcept
	{
		if constexpr (!std::is_same_v<float_t, float>) {
			static_assert(std::is_same_v<float_t, float>, "This function should only be used when the table is single precision, but we want to index with double-precision.");
			return 0.f;
		}
		phase = mspWrap(phase * inv_two_pi) * two_pi;
		return Base::rounded(phase);
	}
	[[nodiscard]]
	constexpr float_t cosLinear(float_t phase) const noexcept
	{
		phase = mspWrap(phase * inv_two_pi) * two_pi;
		return Base::linear(phase);
	}
	[[nodiscard]]
	inline constexpr float_t linear(float_t phase) const noexcept override
	{
		return cosLinear(phase);
	}
	[[nodiscard]]
	constexpr float cosLinear(double phase) const noexcept
	{
		if constexpr (!std::is_same_v<float_t, float>) {
			static_assert(std::is_same_v<float_t, float>, "This function should only be used when the table is single precision, but we want to index with double-precision.");
			return 0.f;
		}
		phase = mspWrap(phase * inv_two_pi) * two_pi;
		return Base::linear(phase);
	}
	[[nodiscard]]
	constexpr float linear(double phase) const noexcept override {
		return cosLinear(phase);
	}
	[[nodiscard]]
	constexpr float_t sin(float_t phase) const noexcept
	{
		phase = phase - math_impl::halfpi<double>();
		phase = mspWrap(phase * inv_two_pi) * two_pi;
		return Base::operator()(phase);
	}
	[[nodiscard]]
	constexpr float sin(double phase) const noexcept
	{
		if constexpr (!std::is_same_v<float_t, float>) {
			static_assert(std::is_same_v<float_t, float>, "This function should only be used when the table is single precision, but we want to index with double-precision.");
			return 0.f;
		}
		phase = phase - math_impl::halfpi<double>();
		phase = mspWrap(phase * inv_two_pi) * two_pi;
		return Base::operator()(phase);
	}
	[[nodiscard]]
	constexpr float_t sinRounded(float_t phase) const noexcept
	{
		phase = phase - math_impl::halfpi<float_t>();
		phase = mspWrap(phase * inv_two_pi) * two_pi;
		return Base::rounded(phase);
	}
	[[nodiscard]]
	constexpr float sinRounded(double phase) const noexcept
	{
		if constexpr (!std::is_same_v<float_t, float>) {
			static_assert(std::is_same_v<float_t, float>, "This function should only be used when the table is single precision, but we want to index with double-precision.");
			return 0.f;
		}
		phase = phase - math_impl::halfpi<float_t>();
		phase = mspWrap(phase * inv_two_pi) * two_pi;
		return Base::rounded(phase);
	}
	[[nodiscard]]
	constexpr float_t sinLinear(float_t phase) const noexcept
	{
		phase = phase - math_impl::halfpi<float_t>();
		phase = mspWrap(phase * inv_two_pi) * two_pi;
		return Base::linear(phase);
	}
	[[nodiscard]]
	constexpr float sinLinear(double phase) const noexcept
	{
		if constexpr (!std::is_same_v<float_t, float>) {
			static_assert(std::is_same_v<float_t, float>, "This function should only be used when the table is single precision, but we want to index with double-precision.");
			return 0.f;
		}
		phase = phase - math_impl::halfpi<float_t>();
		phase = mspWrap(phase * inv_two_pi) * two_pi;
		return Base::linear(phase);
	}
};

#if false // not even slightly tested yet
// unified trig tables with concepts for type safety
template<FloatingPoint T = float, size_t Resolution = 65536>
requires PowerOfTwo<Resolution>
class TrigTables {
private:
	static constexpr T two_pi = math_impl::two_pi<T>();
	static constexpr T pi = math_impl::pi<T>();
	static constexpr size_t mask = Resolution - 1;

	// Structure of Arrays for better cache performance
	alignas(64) std::array<T, Resolution> sin_table_{};
	alignas(64) std::array<T, Resolution> cos_table_{};
	alignas(64) std::array<T, Resolution> tan_table_{};
	alignas(64) std::array<T, Resolution> tanh_table_{};

public:
	constexpr TrigTables() noexcept {
		// Initialize all tables in one pass for better cache usage
		for (size_t i = 0; i < Resolution; ++i) {
			const T phase = (static_cast<T>(i) / static_cast<T>(Resolution)) * two_pi;
			const T half_phase = phase * T{0.5}; // For tan table (0 to π)

			sin_table_[i] = math_impl::sin(phase);
			cos_table_[i] = math_impl::cos(phase);
			tan_table_[i] = math_impl::tan(half_phase);

			// tanh mapping: i ∈ [0, Resolution) → x ∈ [-12, 12)
			const T tanh_input = T{24} * (static_cast<T>(i) / static_cast<T>(Resolution)) - T{12};
			tanh_table_[i] = math_impl::tanh(tanh_input);
		}
	}

	// Unipolar sine lookup [0,1] → [0,2π]
	constexpr T sin_unipolar(T x) const noexcept {
		const T clamped = clamp(x, T{0}, T{1});
		const auto index = static_cast<size_t>(clamped * Resolution) & mask;
		return sin_table_[index];
	}

	// Bipolar sine lookup [-1,1] → [-π,π]
	constexpr T sin_bipolar(T x) const noexcept {
		const T wrapped = wrap(x, T{-1}, T{1});
		const T normalized = (wrapped + 1) * 0.5; // Convert to [0,1]
		return sin_unipolar(normalized);
	}

	// Similar methods for cos
	constexpr T cos_unipolar(T x) const noexcept {
		const T clamped = clamp(x, T{0}, T{1});
		const auto index = static_cast<size_t>(clamped * Resolution) & mask;
		return cos_table_[index];
	}

	constexpr T cos_bipolar(T x) const noexcept {
		const T wrapped = wrap(x, T{-1}, T{1});
		const T normalized = (wrapped + 1) * 0.5;
		return cos_unipolar(normalized);
	}

	// Tan lookup [0,1] → [0,π]
	constexpr T tan_lookup(T x) const noexcept {
		const T clamped = clamp(x, T{0}, T{1});
		const auto index = static_cast<size_t>(clamped * Resolution) & mask;
		return tan_table_[index];
	}

	// Tanh lookup [-12,12]
	constexpr T tanh_lookup(T x) const noexcept {
		const T clamped = clamp(x, T{-12}, T{12});
		const T normalized = (clamped + 12) / 24; // Convert to [0,1]
		const auto index = static_cast<size_t>(normalized * Resolution) & mask;
		return tanh_table_[index];
	}

	// Interpolated versions for higher quality
	constexpr T sin_unipolar_interp(T x) const noexcept {
		const T clamped = clamp(x, T{0}, T{1});
		return interpolated_lookup(sin_table_, clamped);
	}

	constexpr T cos_unipolar_interp(T x) const noexcept {
		const T clamped = clamp(x, T{0}, T{1});
		return interpolated_lookup(cos_table_, clamped);
	}

	constexpr T tanh_lookup_interp(T x) const noexcept {
		const T clamped = clamp(x, T{-12}, T{12});
		const T normalized = (clamped + 12) / 24;
		return interpolated_lookup(tanh_table_, normalized);
	}

	// Access to raw tables if needed for advanced use cases
	constexpr const std::array<T, Resolution>& sin_data() const noexcept { return sin_table_; }
	constexpr const std::array<T, Resolution>& cos_data() const noexcept { return cos_table_; }
	constexpr const std::array<T, Resolution>& tan_data() const noexcept { return tan_table_; }
	constexpr const std::array<T, Resolution>& tanh_data() const noexcept { return tanh_table_; }

private:
	constexpr T interpolated_lookup(const std::array<T, Resolution>& table, T normalized_x) const noexcept {
		const T scaled_index = normalized_x * static_cast<T>(Resolution);
		const auto rounded_index = nvs::memoryless::round(scaled_index) & mask;
		const auto ceil_index = (rounded_index + 1) & mask;
		const T frac = scaled_index - static_cast<T>(rounded_index);

		return linterp(table[rounded_index], table[ceil_index], frac);
	}
};
#endif

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
