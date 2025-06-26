/*
  ==============================================================================

    nvs_math_impl.h
    Created: 25 Jun 2025 6:37:14pm
    Author:  Nicholas Solem

  ==============================================================================
*/

#pragma once
#include "nvs_memoryless.h"
#include "lookup_tables.h"

namespace nvs::memoryless {

namespace math_impl {

template <FloatingPoint T>
static constexpr T tan(T x){
#ifdef TAN_LUT
	if constexpr (std::is_same_v<T, float>){
		return tanTable_f(x);
	}
	else if constexpr (std::is_same_v<T, double>){
		return tanTable_d(x);
	}
#else
	constexpr T inv_half_pi = 0.636619772367581;
	return fastTanPi_2_acc8(inv_half_pi * x);
#endif
}


template <FloatingPoint T>
static constexpr T tanh(T x){
#ifdef TANH_LUT
	if constexpr (std::is_same_v<T, float>){
		return tanhTable_f(x);
	}
	else if constexpr (std::is_same_v<T, double>){
		return tanhTable_d(x);
	}
#else
	return padeTanh(x);
#endif
}


template <FloatingPoint T>
static constexpr T cos(T phase) {
#ifdef COS_LUT
	if constexpr (std::is_same_v<T, float>){
		return cosTable_f.cos(phase);
	}
	else if constexpr (std::is_same_v<T, double>){
		return cosTable_d.cos(phase);
	}
#else
	return padeCos(phase);
#endif
}


template <FloatingPoint T>
static constexpr T sin(T phase) {
#ifdef COS_LUT
	if constexpr (std::is_same_v<T, float>){
		return cosTable_f.sin(phase);
	}
	else if constexpr (std::is_same_v<T, double>){
		return cosTable_d.sin(phase);
	}
#else
	return padeSin(phase);
#endif
}

}	// namespace math_impl

}	// namespace nvs::memoryless
