//
// Created by Nicholas Solem on 6/9/25.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "catch2/generators/catch_generators.hpp"
#include <catch2/generators/catch_generators_adapters.hpp>

#include <cmath>
#include <type_traits>

#include "../include/nvs_memoryless.h"
#include "catch2/generators/catch_generators_random.hpp"
#include "catch2/generators/catch_generators_range.hpp"


// Include your lookup table headers here
// #include "lookup_tables.hpp"

using Catch::Approx;
using namespace nvs::memoryless;

// Test configuration - you can toggle this
#ifndef USING_SPROUT
#define USING_SPROUT 0  // Set to 1 to test constexpr construction
#endif


// Generic helper function to test accuracy against standard library
template<typename TableType, typename RefFunc, typename TableFunc>
static void test_table_accuracy(const TableType& table, RefFunc ref_func, TableFunc table_func,
                        typename TableType::value_type min_val,
                        typename TableType::value_type max_val,
                        double tolerance, const std::string& method_name,
                        bool use_relative_tolerance = false) {
    constexpr int num_test_points = 100;
    const auto range = max_val - min_val;
    const auto step = range / (num_test_points - 1);

    for (int i = 0; i < num_test_points; ++i) {
        const auto x = min_val + i * step;
        const auto table_result = table_func(table, x);
        const auto expected = static_cast<typename TableType::value_type>(ref_func(x));
        const auto diff = std::abs(expected - table_result);

        INFO("Testing " << method_name << " accuracy at x = " << x << " (point " << i << "/" << num_test_points << ")");
        CAPTURE(x, table_result, expected, tolerance, diff, use_relative_tolerance);

        if (tolerance == 0.0) {
            // For exact comparison, bypass Approx and compare directly
            REQUIRE(table_result == expected);
        }
        else if (use_relative_tolerance) {
            // Use relative tolerance (error proportional to expected value)
            if (std::abs(expected) < std::numeric_limits<typename TableType::value_type>::epsilon()) {
                // For values very close to zero, fall back to absolute tolerance
                REQUIRE(std::abs(table_result - expected) <= tolerance);
            } else {
                // Relative error: |actual - expected| / |expected| <= tolerance
                const auto relative_error = std::abs(table_result - expected) / std::abs(expected);
                CAPTURE(relative_error);
                REQUIRE(relative_error <= tolerance);
            }
        }
        else {
            // Use absolute tolerance (traditional approach)
            // REQUIRE(table_result == Approx(expected).margin(tolerance).epsilon(0.0));
            REQUIRE(table_result == Approx(expected).margin(tolerance));
        }
    }
}

// Enhanced helper for relative tolerance testing
template<typename TableType, typename RefFunc, typename TableFunc>
static void test_table_relative_accuracy(const TableType& table, RefFunc ref_func, TableFunc table_func,
                                 typename TableType::value_type min_val,
                                 typename TableType::value_type max_val,
                                 double relative_tolerance, const std::string& method_name,
                                 double absolute_tolerance_near_zero = 1e-10) {
    constexpr int num_test_points = 100;
    const auto range = max_val - min_val;
    const auto step = range / (num_test_points - 1);

    for (int i = 0; i < num_test_points; ++i) {
        const auto x = min_val + i * step;
        const auto table_result = table_func(table, x);
        const auto expected = static_cast<typename TableType::value_type>(ref_func(x));

        INFO("Testing " << method_name << " relative accuracy at x = " << x << " (point " << i << "/" << num_test_points << ")");

        const auto abs_error = std::abs(table_result - expected);
        const auto abs_expected = std::abs(expected);

        CAPTURE(x, table_result, expected, abs_error, abs_expected, relative_tolerance);

        if (abs_expected < absolute_tolerance_near_zero) {
            // Near zero, use absolute tolerance
            INFO("Using absolute tolerance near zero");
            REQUIRE(abs_error <= absolute_tolerance_near_zero);
        }
        else {
            // Use relative tolerance
            const auto relative_error = abs_error / abs_expected;
            CAPTURE(relative_error);
            INFO("Relative error: " << relative_error << " (limit: " << relative_tolerance << ")");
            REQUIRE(relative_error <= relative_tolerance);
        }
    }
}

// Hybrid tolerance testing - ideal for trigonometric functions
template<typename TableType, typename RefFunc, typename TableFunc>
static void test_table_hybrid_accuracy(const TableType& table, RefFunc ref_func, TableFunc table_func,
                               typename TableType::value_type min_val,
                               typename TableType::value_type max_val,
                               double relative_tolerance, double absolute_tolerance,
                               const std::string& method_name) {
    constexpr int num_test_points = 100;
    const auto range = max_val - min_val;
    const auto step = range / (num_test_points - 1);

    for (int i = 0; i < num_test_points; ++i) {
        const auto x = min_val + i * step;
        const auto table_result = table_func(table, x);
        const auto expected = static_cast<typename TableType::value_type>(ref_func(x));

        INFO("Testing " << method_name << " hybrid accuracy at x = " << x << " (point " << i << "/" << num_test_points << ")");

        const auto abs_error = std::abs(table_result - expected);
        const auto abs_expected = std::abs(expected);

        CAPTURE(x, table_result, expected, abs_error, abs_expected);
        CAPTURE(relative_tolerance, absolute_tolerance);

        // Hybrid approach: pass if EITHER absolute OR relative tolerance is met
        bool passes_absolute = abs_error <= absolute_tolerance;
        bool passes_relative = (abs_expected > absolute_tolerance) ?
                              (abs_error / abs_expected <= relative_tolerance) : true;

        CAPTURE(passes_absolute, passes_relative);

        if (passes_absolute) {
            INFO("Passed absolute tolerance check (error: " << abs_error << " <= " << absolute_tolerance << ")");
        } else if (passes_relative) {
            const auto relative_error = abs_error / abs_expected;
            INFO("Passed relative tolerance check (rel_error: " << relative_error << " <= " << relative_tolerance << ")");
            CAPTURE(relative_error);
        } else {
            const auto relative_error = abs_expected > 0 ? abs_error / abs_expected : std::numeric_limits<double>::infinity();
            INFO("Failed both tolerance checks");
            CAPTURE(relative_error);
            FAIL("Error exceeds both absolute (" << absolute_tolerance << ") and relative (" << relative_tolerance << ") tolerances");
        }

        REQUIRE((passes_absolute || passes_relative));
    }
}

// Convenience wrappers
template<typename TableType, typename RefFunc>
static void test_accuracy(const TableType& table, RefFunc ref_func,
                  typename TableType::value_type min_val,
                  typename TableType::value_type max_val,
                  double tolerance = 0.01) {
    test_table_accuracy(table, ref_func,
                       [](const auto& t, auto x) { return t(x); },
                       min_val, max_val, tolerance, "non-interpolated");
}

template<typename TableType, typename RefFunc>
static void test_interpolated_accuracy(const TableType& table, RefFunc ref_func,
                               typename TableType::value_type min_val,
                               typename TableType::value_type max_val,
                               double tolerance = 0.005) {
    test_table_accuracy(table, ref_func,
                       [](const auto& t, auto x) { return t.linear(x); },
                       min_val, max_val, tolerance, "interpolated");
}

template<typename TableType, typename RefFunc>
static void test_relative_accuracy(const TableType& table, RefFunc ref_func,
                           typename TableType::value_type min_val,
                           typename TableType::value_type max_val,
                           double relative_tolerance = 0.01) {
    test_table_relative_accuracy(table, ref_func,
                                [](const auto& t, auto x) { return t(x); },
                                min_val, max_val, relative_tolerance, "non-interpolated");
}

template<typename TableType, typename RefFunc>
static void test_interpolated_relative_accuracy(const TableType& table, RefFunc ref_func,
                                        typename TableType::value_type min_val,
                                        typename TableType::value_type max_val,
                                        double relative_tolerance = 0.005) {
    test_table_relative_accuracy(table, ref_func,
                                [](const auto& t, auto x) { return t.linear(x); },
                                min_val, max_val, relative_tolerance, "interpolated");
}

template<typename TableType, typename RefFunc>
static void test_hybrid_accuracy(const TableType& table, RefFunc ref_func,
                         typename TableType::value_type min_val,
                         typename TableType::value_type max_val,
                         double relative_tolerance = 0.01,
                         double absolute_tolerance = 1e-6) {
    test_table_hybrid_accuracy(table, ref_func,
                              [](const auto& t, auto x) { return t(x); },
                              min_val, max_val, relative_tolerance, absolute_tolerance, "non-interpolated");
}

template<typename TableType, typename RefFunc>
static void test_interpolated_hybrid_accuracy(const TableType& table, RefFunc ref_func,
                                      typename TableType::value_type min_val,
                                      typename TableType::value_type max_val,
                                      double relative_tolerance = 0.005,
                                      double absolute_tolerance = 1e-6) {
    test_table_hybrid_accuracy(table, ref_func,
                              [](const auto& t, auto x) { return t.interpolated(x); },
                              min_val, max_val, relative_tolerance, absolute_tolerance, "interpolated");
}

TEST_CASE("ExpTable - Basic Properties", "[lookup][exp]") {
    constexpr int min_val = -5;
    constexpr int max_val = 5;
    constexpr size_t resolution = 64;

    REQUIRE(resolution > 0);

    ExpTable<float, IRange<min_val, max_val>, resolution, InterpolationType::Rounded>  const exp_table;

    // SECTION("Table properties") {
    //     REQUIRE(exp_table.min_val == min_val);
    //     REQUIRE(exp_table.max_val == max_val);
    //     REQUIRE(exp_table.size == resolution);
    //     REQUIRE(exp_table.range_val == Approx(max_val - min_val));
    // }

    SECTION("Boundary values") {
        // Test at boundaries
        auto result_min = exp_table(static_cast<float>(min_val));
        auto expected_min = std::exp(static_cast<float>(min_val));
        REQUIRE(result_min == Approx(expected_min).margin(0.1));

        auto result_max = exp_table(static_cast<float>(max_val));
        auto expected_max = std::exp(static_cast<float>(max_val));
        REQUIRE(result_max == Approx(expected_max).margin(1.0)); // exp(5) is large, allow more margin
    }

    SECTION("Out of bounds clamping") {
        // Values below min should clamp to min behavior
        auto result_below = exp_table(static_cast<float>(min_val - 10));
        auto result_min = exp_table(static_cast<float>(min_val));
        REQUIRE(result_below == Approx(result_min).margin(0.01));

        // Values above max should clamp to max behavior
        auto result_above = exp_table(static_cast<float>(max_val + 10));
        auto result_max = exp_table(static_cast<float>(max_val));
        REQUIRE(result_above == Approx(result_max).margin(0.01));
    }
}


TEST_CASE("ExpTable - Accuracy Tests", "[lookup][exp][accuracy]") {
    constexpr int min_val = -3;
    constexpr int max_val = 3;
    constexpr size_t resolution = 1024;

    ExpTable<float, IRange<min_val, max_val>, resolution, InterpolationType::Rounded> const exp_table;

    SECTION("Non-interpolated accuracy (absolute tolerance)") {
        test_accuracy(exp_table, [](float x) { return std::exp(x); }, min_val, max_val, 0.08);
    }
    SECTION("Non-interpolated relative accuracy") {
        test_relative_accuracy(exp_table, [](float x) { return std::exp(x); },
                              static_cast<float>(min_val), static_cast<float>(max_val), 0.006);
    }

    SECTION("Interpolated accuracy (absolute tolerance)") {
        test_interpolated_accuracy(exp_table, [](float x) { return std::exp(x); },
                                  static_cast<float>(min_val), static_cast<float>(max_val), 0.0000001);
    }

    SECTION("Interpolated relative accuracy") {
        test_interpolated_relative_accuracy(exp_table, [](float x) { return std::exp(x); },
                                           static_cast<float>(min_val), static_cast<float>(max_val), 0.000005);
    }

    SECTION("Monotonicity") {
        // Exponential function should be strictly increasing
        constexpr int num_points = 50;
        const float step = (max_val - min_val) / static_cast<float>(num_points - 1);

        float prev_value = exp_table(static_cast<float>(min_val));
        for (int i = 1; i < num_points; ++i) {
            const float x = min_val + i * step;
            const float current_value = exp_table(x);
            REQUIRE(current_value > prev_value);
            prev_value = current_value;
        }
    }
}
TEST_CASE("CosSinTable - Wraparound Behavior", "[lookup][sin][cos][wraparound]") {
    constexpr auto min_val = 0.f;
    constexpr auto twopi = math_impl::two_pi<double>();
    static_assert(std::is_floating_point_v<decltype(min_val)>);
    constexpr auto max_val = twopi;
    constexpr size_t resolution = 1024;

    CosTable<float, resolution, InterpolationType::Rounded> const cos_sin_table;

    SECTION("CosSinTable - Wraparound Behavior - Random Offsets") {
        auto const phi = GENERATE(take(50, Catch::Generators::random(0.0, twopi)));
        auto const multiple = GENERATE(Catch::Generators::range<int>(-100000, 100000, 5321));

        // Capture values for debugging
        CAPTURE(phi, multiple);
        auto const phi_multiple =  (multiple * twopi);
        INFO("Testing phi=" << phi << ", multiple=" << multiple << ", offset=" << phi_multiple);
        {
            float const expected = cos_sin_table.rounded(static_cast<double>(phi));
            float const actual_pos = cos_sin_table.rounded(static_cast<double>(phi + phi_multiple));
            float const actual_neg = cos_sin_table.rounded(static_cast<double>(phi - phi_multiple));

            CAPTURE(expected, actual_pos, actual_neg);

            // Test positive wraparound
            REQUIRE(expected == Approx(actual_pos));

            // Test negative wraparound
            REQUIRE(expected == Approx(actual_neg));
        }
        {
            float const expected = cos_sin_table.linear(static_cast<double>(phi));
            float const actual_pos = cos_sin_table.linear(static_cast<double>(phi + phi_multiple));
            float const actual_neg = cos_sin_table.linear(static_cast<double>(phi - phi_multiple));

            CAPTURE(expected, actual_pos, actual_neg);

            // Test positive wraparound
            REQUIRE(expected == Approx(actual_pos));

            // Test negative wraparound
            REQUIRE(expected == Approx(actual_neg));
        }
    }
}

TEST_CASE("SinCosTable - Accuracy Tests", "[lookup][sin][cos][accuracy]") {
    auto min_val = 0.f;
    static_assert(std::is_floating_point_v<decltype(min_val)>);
    auto max_val = math_impl::two_pi<float>();
    constexpr size_t resolution = 16384;

    CosTable<float, resolution, InterpolationType::Rounded> const cos_sin_table;

    static_assert(std::is_convertible_v<CosTable<float, 16384, InterpolationType::Rounded>::value_type, decltype(min_val)>);

    SECTION("Key points") {
        REQUIRE(cos_sin_table.sin(-M_PI) == Approx(std::sin(-M_PI)).margin(0.00001));
        REQUIRE(cos_sin_table.sin(-M_PI/2) == Approx(std::sin(-M_PI/2)).margin(0.00001));
        auto const a = cos_sin_table.sin(-M_PI/4);
        REQUIRE(a == Approx(std::sin(-M_PI/4)).margin(0.00001));
        REQUIRE(cos_sin_table.sin(0.0f) == Approx(std::sin(0.0f)).margin(0.00001));
        REQUIRE(cos_sin_table.sin(M_PI/4) == Approx(std::sin(M_PI/4)).margin(0.00001));
        REQUIRE(cos_sin_table.sin(M_PI/2) == Approx(std::sin(M_PI/2)).margin(0.00001));
        REQUIRE(cos_sin_table.sin(M_PI) == Approx(std::sin(M_PI)).margin(0.00001));
    }
    SECTION("Misc") {
        float val = 0.1f;

        INFO("error: " << cos_sin_table.sin(val) - std::sin(val));
        REQUIRE(cos_sin_table.sin(val) == Approx(std::sin(val)).margin(0.002));
        val = -0.1f;
        INFO("error: " << cos_sin_table.sin(val) - std::sin(val));
        REQUIRE(cos_sin_table.sin(val) == Approx(std::sin(val)).margin(0.002));
    }

    // WHY COS? Because the DEFAULT function, returned by straight up operator()(float_t), rounded() and linear(), is cos.
    SECTION("Non-interpolated accuracy") {
        test_accuracy(cos_sin_table, [](float x) { return std::cos(x); },0.f, 1.f, 2e-4);
    }
    SECTION("Non-interpolated hybrid accuracy") {
        // Allow 1% relative error OR 1e-6 absolute error (whichever is more lenient)
        test_hybrid_accuracy(cos_sin_table, [](float x) { return std::cos(x); },
                           min_val, max_val, 0.01, 1e-5);
    }
    SECTION("Interpolated hybrid accuracy") {
        test_interpolated_accuracy(cos_sin_table, [](float x) { return std::cos(x); },
                                        min_val, max_val, 2e-7);
    }
}

//
// #if USING_SPROUT
// TEST_CASE("Constexpr Construction", "[lookup][constexpr]") {
//     SECTION("ExpTable constexpr construction") {
//         constexpr ExpTable<-2, 2, float, 256> exp_table{};
//         static_assert(exp_table.size == 256);
//         static_assert(exp_table.min_val == -2);
//         static_assert(exp_table.max_val == 2);
//
//         // Test that we can use it in constexpr context
//         constexpr auto result = exp_table(1.0f);
//         REQUIRE(result == Approx(std::exp(1.0f)).margin(0.1));
//     }
// }
// #endif