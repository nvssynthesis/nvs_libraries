#include "../include/nvs_gen.h"
#include <iostream>
#include <vector>
#include <string>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_vector.hpp"

using namespace nvs::gen;

TEST_CASE("change function tests", "[change]") {
	change<float> c;
	std::vector<float> const input = {0.f, 0.1f, 0.01f, -100.f, -0.5f};
	std::vector<float> const expected_output = {
		0.f, 1.f, -1.f, -1.f, 1.f
	};

	for (auto i = 0; i < input.size(); ++i) {
		auto const out = c(input[i]);
		REQUIRE(out == expected_output[i]);
	}
}

std::vector<float> makePhasorWave(size_t const length, const float freq, const float fs, const float initialPhase=0.f){
	std::vector<float> wave(length);
	phasor<float> ph;
	ph.setSampleRate(fs);
	ph.setFrequency(freq);
	ph.setPhase(initialPhase);
	for (auto &e : wave) {
		e = (ph++).getPhase();
	}
    REQUIRE(wave.size() == length);
	return wave;
}

TEST_CASE("phasor tests", "[phasor]") {
	constexpr size_t N = 10;

	SECTION("nyquist frequency test") {
		std::vector<float> const ph = makePhasorWave(N, 22050.f, 44100.f);
		std::vector<float> const expected_out = {
			0.f, 0.5f, 0.f, 0.5f, 0.f, 0.5f, 0.f, 0.5f, 0.f, 0.5f
		 };

		REQUIRE(ph.size() == expected_out.size());
		for (size_t i = 0; i < ph.size(); ++i) {
			REQUIRE(ph[i] == expected_out[i]);
		}
	}

	SECTION("DC (0Hz) frequency test") {
		std::vector<float> const ph = makePhasorWave(N, 0.f, 44100.f);
		std::vector<float> const expected_out = {
			0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f
		 };

		REQUIRE(ph.size() == expected_out.size());
		REQUIRE(ph.size() == N);

		for (size_t i = 0; i < ph.size(); ++i) {
			INFO("Index: " << i << ", Expected: " << expected_out[i] << ", Actual: " << ph[i]);
			REQUIRE(ph[i] == expected_out[i]);
		}
	}
}
TEST_CASE("ramp2trig tests", "[ramp2trig]") {
	SECTION("basic trigger detection") {
		ramp2trig<float> r2t;

		// Test with manually crafted input that has clear discontinuities
		std::vector<float> const input = {
			0.0f,   // Initial state
			0.0f,   // Initial state
			0.25f,  // Sudden increment
			0.5f,   // Normal increment - no trigger expected
			0.75f,  // Normal increment - no trigger expected
			0.1f,   // RESET/WRAP - trigger expected here!
			0.35f,  // Normal increment after reset
			0.6f,   // Normal increment
			0.85f,  // Normal increment
			0.05f   // RESET/WRAP - trigger expected here!
		};

		std::vector<float> results;
		for (float val : input) {
			results.push_back(r2t(val));
		}

		// First sample behavior might be special, so test from index 1
		REQUIRE(results[0] == 0.f); // No activity
		REQUIRE(results[1] == 0.f); // No activity
		REQUIRE(results[2] == 1.f); // Sudden increment
		REQUIRE(results[3] == 0.f); // Normal increment
		REQUIRE(results[4] == 0.f); // Normal increment
		REQUIRE(results[5] == 1.f); // Reset detected!
		REQUIRE(results[6] == 0.f); // Normal increment after reset
		REQUIRE(results[7] == 0.f); // Normal increment
		REQUIRE(results[8] == 0.f); // Normal increment
		REQUIRE(results[9] == 1.f); // Another reset detected!
	}
}

TEST_CASE("switcher tests", "[switcher]") {
	SECTION("basic switching behavior") {
		std::vector<int> v1 =       {1, 2, 3, 4, 5, 6, 7, 8, 9};
		std::vector<int> v2 =       {9, 8, 7, 6, 5, 4, 3, 2, 1};
		std::vector<int> control =  {0, 10, 0, -1, 0, 5, 0, -100, 0};
		std::vector<int> expected = {9, 2, 7, 4, 5, 6, 3, 8, 1};

		REQUIRE(v1.size() == v2.size());
		REQUIRE(v1.size() == control.size());
		REQUIRE(v1.size() == expected.size());

		for (size_t i = 0; i < v1.size(); ++i) {
			auto out = switcher(control[i], v1[i], v2[i]);
			INFO("Index: " << i << ", Control: " << control[i] << ", V1: " << v1[i] << ", V2: " << v2[i] << ", Expected: " << expected[i] << ", Got: " << out);
			REQUIRE(out == expected[i]);
		}
	}
}

template<typename T>
void print(T x, std::string s = " ") {
	std::cout << x << s;
}
template<typename T, size_t N>
void print(std::array<T, N> xArr, std::string s = " "){
	for (auto e : xArr)
		print(e, s);
}

TEST_CASE("gateSelect tests", "[gateSelect]") {
	SECTION("gate selection with different indices") {
		constexpr size_t N = 3;
		constexpr float val = 10.f;

		for (size_t i = 0; i < N; ++i) {
			std::array<float, N> res = gateSelect<float, N>(i, val);
			std::vector<float> res_vec(res.begin(), res.end());  // Convert to vector

			if (i > 0) {
				std::vector<float> selection(N, 0.f);
				selection[i - 1] = val;
				REQUIRE_THAT(res_vec, Catch::Matchers::Equals(selection));
			} else {
				std::vector<float> const zeros(N, 0.f);
				REQUIRE_THAT(res_vec, Catch::Matchers::Equals(zeros));
			}
		}
	}
}

TEST_CASE("latch tests", "[latch]") {
	SECTION("latch behavior with phasor input") {
		phasor<float> ph;
		ph.setSampleRate(10.f);
		ph.setFrequency(2.f);

		constexpr size_t N = 10;
		std::array<float, N> phRec;
		for (auto &e : phRec) {
			e = (ph++).getPhase();
		}

		std::array<float, N> const expected = {
			0.f, 0.f, 0.4f, 0.6f, 0.8f, 0.f, 0.2f, 0.4f, 0.4f, 0.4f
		};

		latch<float> l;
		bool letPass = false;

		for (size_t i = 0; i < phRec.size(); ++i) {
			// Update letPass based on sample index
			if (i == 2) {        // i+1 == 3 in your original (since you increment before check)
				letPass = true;
			} else if (i == 8) { // i+1 == 9 in your original
				letPass = false;
			}

			float const input = phRec[i];
			float const output = l(input, letPass);
			float const exp = expected[i];

			INFO("Index: " << i << ", Input: " << input << ", LetPass: " << letPass
				 << ", Expected: " << exp << ", Got: " << output);
			REQUIRE(output == exp);
		}
	}
}

TEST_CASE("parzen window tests", "[parzen][windowing]") {
	using Catch::Approx;
	SECTION("boundary conditions") {
		// Half-window should be 0 and 1 at edges (0 and 1)
		REQUIRE(parzen<float>(0.0f) == Approx(0.0f));
		REQUIRE(parzen<float>(1.0f) == Approx(1.0f));
	}

	SECTION("monotonic increasing") {
		// Should be strictly increasing from 0 to 1
		constexpr size_t N = 20;
		std::vector<float> values;

		for (size_t i = 0; i <= N; ++i) {
			float t = static_cast<float>(i) / N;
			values.push_back(parzen<float>(t));
		}

		for (size_t i = 1; i < values.size(); ++i) {
			INFO("Checking monotonicity at index " << i << ": " << values[i-1] << " < " << values[i]);
			REQUIRE(values[i] >= values[i-1]); // Non-decreasing
		}
	}

	SECTION("curvature properties") {
		// Test that it's curved (not linear)
		float quarter = parzen<float>(0.25f);
		float half = parzen<float>(0.5f);
		float three_quarter = parzen<float>(0.75f);

		// For Parzen curvature, early values should be lower than linear
		REQUIRE(quarter < 0.25f);  // Should be below linear interpolation
		REQUIRE(half < 0.5f);      // Should be below linear interpolation

		// But still increasing
		REQUIRE(quarter < half);
		REQUIRE(half < three_quarter);
	}
	SECTION("range constraints") {
		// All values should be in [0,1] range
		for (float t = 0.0f; t <= 1.0f; t += 0.1f) {
			float val = parzen<float>(t);
			INFO("Testing range at t=" << t << ", value=" << val);
			REQUIRE(val >= 0.0f);
			REQUIRE(val <= 1.0f);
		}
	}
	SECTION("creates symmetric parzen from triangle") {
		// Create a triangle wave (0→1→0)
		std::vector<float> const triangle = {0.0f, 0.25f, 0.5f, 0.75f, 1.0f, 0.75f, 0.5f, 0.25f, 0.0f};
		std::vector<float> windowed;

		for (float t : triangle) {
			windowed.push_back(parzen<float>(t));
		}

		// Result should be symmetric
		REQUIRE(windowed[0] == Approx(windowed[8])); // edges
		REQUIRE(windowed[1] == Approx(windowed[7])); // inner points
		REQUIRE(windowed[2] == Approx(windowed[6]));
		REQUIRE(windowed[3] == Approx(windowed[5]));
		// windowed[4] is the peak
	}
}

