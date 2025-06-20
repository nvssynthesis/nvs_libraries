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
    constexpr size_t cache_line_size = 64; // Safe assumption for most modern processors
#endif

    // really want 16384
    alignas(cache_line_size)
    extern const memoryless::CosTable<float, float, 16384, nvs::memoryless::InterpolationType::Linear> cosTable_f;
    alignas(cache_line_size)
    extern const memoryless::CosTable<double, double, 16384, nvs::memoryless::InterpolationType::Linear> cosTable_d;

    // really want 65536
    alignas(cache_line_size)
    extern const nvs::memoryless::TanTable<float, float, 65536, nvs::memoryless::InterpolationType::Linear> tanTable_f;
    alignas(cache_line_size)
    extern const nvs::memoryless::TanTable<double, double, 65536, nvs::memoryless::InterpolationType::Linear> tanTable_d;

    // might be fine this small
    alignas(cache_line_size)
    extern const nvs::memoryless::TanhTable<float, float, IRange<-20, 20>, 16384, nvs::memoryless::InterpolationType::Linear> tanhTable_f;
    alignas(cache_line_size)
    extern const nvs::memoryless::TanhTable<double, double, IRange<-20, 20>, 16384, nvs::memoryless::InterpolationType::Linear> tanhTable_d;

    // probably want more like 65536. see if we can get away with more restricted range too.
    alignas(cache_line_size)
    extern const AtanTable<float, float, IRange<-96, 96>, 65536, InterpolationType::Linear> atanTable_f;
    alignas(cache_line_size)
    extern const AtanTable<double, double, IRange<-96, 96>, 65536, InterpolationType::Linear> atanTable_d;

    // not sure about range or resolution needed.
    alignas(cache_line_size)
    extern const ExpTable<float, float, IRange<-10, 10>, 1024, InterpolationType::Rounded> expTable_f;
    alignas(cache_line_size)
    extern const ExpTable<double, double, IRange<-10, 10>, 1024, InterpolationType::Rounded> expTable_d;

}

#endif //LOOKUP_TABLES_H
