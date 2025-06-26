//
// Created by Nicholas Solem on 6/19/25.
//
#include "../include/lookup_tables.h"

namespace nvs::memoryless {

#ifdef COS_LUT
alignas(cache_line_size)
const memoryless::CosTable<float, float, COS_SIZE, nvs::memoryless::InterpolationType::Linear> cosTable_f {};
alignas(cache_line_size)
const memoryless::CosTable<double, double, COS_SIZE, nvs::memoryless::InterpolationType::Linear> cosTable_d {};
#endif	// COS_LUT

#ifdef TAN_LUT
alignas(cache_line_size)
const nvs::memoryless::TanTable<float, float, TAN_SIZE, nvs::memoryless::InterpolationType::Linear> tanTable_f {};
alignas(cache_line_size)
const nvs::memoryless::TanTable<double, double, TAN_SIZE, nvs::memoryless::InterpolationType::Linear> tanTable_d {};
#endif	// TAN_LUT

#ifdef TANH_LUT
alignas(cache_line_size)
const nvs::memoryless::TanhTable<float, float, IRange<-20, 20>, TANH_SIZE, nvs::memoryless::InterpolationType::Linear> tanhTable_f {};
alignas(cache_line_size)
const nvs::memoryless::TanhTable<double, double, IRange<-20, 20>, TANH_SIZE, nvs::memoryless::InterpolationType::Linear> tanhTable_d {};
#endif	// TANH_LUT

#ifdef ATAN_LUT
alignas(cache_line_size)
const AtanTable<float, float, IRange<-96, 96>, ATAN_SIZE, InterpolationType::Linear> atanTable_f {};
alignas(cache_line_size)
const AtanTable<double, double, IRange<-96, 96>, ATAN_SIZE, InterpolationType::Linear> atanTable_d {};
#endif	// ATAN_LUT

#ifdef EXP_LUT
alignas(cache_line_size)
const ExpTable<float, float, IRange<-10, 10>, EXP_SIZE, InterpolationType::Rounded> expTable_f {};
alignas(cache_line_size)
const ExpTable<double, double, IRange<-10, 10>, EXP_SIZE, InterpolationType::Rounded> expTable_d {};
#endif // EXP_LUT

}
