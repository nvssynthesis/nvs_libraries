//
// Created by Nicholas Solem on 6/19/25.
//

#ifndef LOOKUP_TABLES_H
#define LOOKUP_TABLES_H
#include "nvs_memoryless.h"

namespace nvs::memoryless {

#ifdef __cpp_lib_hardware_interference_size
    constexpr std::size_t cache_line_size = std::hardware_destructive_interference_size;
#else
    constexpr std::size_t cache_line_size = 64; // Safe assumption for most modern processors
#endif


#ifdef COS_LUT
    // really want 16384
#ifndef COS_SIZE
static constexpr int COS_SIZE {4096};
#endif
    alignas(cache_line_size)
    extern const memoryless::CosTable<float, float, COS_SIZE, nvs::memoryless::InterpolationType::Linear> cosTable_f;
    alignas(cache_line_size)
    extern const memoryless::CosTable<double, double, COS_SIZE, nvs::memoryless::InterpolationType::Linear> cosTable_d;
#endif	// COS_LUT


#ifdef TAN_LUT
    // really want 65536
#ifndef TAN_SIZE
static constexpr int TAN_SIZE {16384};
#endif
    alignas(cache_line_size)
    extern const nvs::memoryless::TanTable<float, float, TAN_SIZE, nvs::memoryless::InterpolationType::Linear> tanTable_f;
    alignas(cache_line_size)
    extern const nvs::memoryless::TanTable<double, double, TAN_SIZE, nvs::memoryless::InterpolationType::Linear> tanTable_d;
#endif	// TAN_LUT


#ifdef TANH_LUT
    // might be fine this small
#ifndef TANH_SIZE
static constexpr int TANH_SIZE {2048};
#endif
    alignas(cache_line_size)
    extern const nvs::memoryless::TanhTable<float, float, IRange<-20, 20>, TANH_SIZE, nvs::memoryless::InterpolationType::Linear> tanhTable_f;
    alignas(cache_line_size)
    extern const nvs::memoryless::TanhTable<double, double, IRange<-20, 20>, TANH_SIZE, nvs::memoryless::InterpolationType::Linear> tanhTable_d;
#endif	// TANH_LUT


#ifdef ATAN_LUT
    // probably want more like 65536. see if we can get away with more restricted range too.
#ifndef ATAN_SIZE
static constexpr int ATAN_SIZE {32768};
#endif
    alignas(cache_line_size)
    extern const AtanTable<float, float, IRange<-96, 96>, ATAN_SIZE, InterpolationType::Linear> atanTable_f;
    alignas(cache_line_size)
    extern const AtanTable<double, double, IRange<-96, 96>, ATAN_SIZE, InterpolationType::Linear> atanTable_d;
#endif // ATAN_LUT


#ifdef EXP_LUT
#ifndef EXP_SIZE
static constexpr int EXP_SIZE{1024};
#endif
    alignas(cache_line_size)
    extern const ExpTable<float, float, IRange<-10, 10>, EXP_SIZE, InterpolationType::Rounded> expTable_f;
    alignas(cache_line_size)
    extern const ExpTable<double, double, IRange<-10, 10>, EXP_SIZE, InterpolationType::Rounded> expTable_d;
#endif	// EXP_LUT

}
#endif //LOOKUP_TABLES_H
