//
// Created by Nicholas Solem on 6/19/25.
//
#include "../include/lookup_tables.h"

namespace nvs::memoryless {

alignas(cache_line_size)
const memoryless::CosTable<float, float, 16384, nvs::memoryless::InterpolationType::Linear> cosTable_f {};
alignas(cache_line_size)
const memoryless::CosTable<double, double, 16384, nvs::memoryless::InterpolationType::Linear> cosTable_d {};

// really want 65536
alignas(cache_line_size)
const nvs::memoryless::TanTable<float, float, 65536, nvs::memoryless::InterpolationType::Linear> tanTable_f {};
alignas(cache_line_size)
const nvs::memoryless::TanTable<double, double, 65536, nvs::memoryless::InterpolationType::Linear> tanTable_d {};

// might be fine this small
alignas(cache_line_size)
const nvs::memoryless::TanhTable<float, float, IRange<-20, 20>, 16384, nvs::memoryless::InterpolationType::Linear> tanhTable_f {};
alignas(cache_line_size)
const nvs::memoryless::TanhTable<double, double, IRange<-20, 20>, 16384, nvs::memoryless::InterpolationType::Linear> tanhTable_d {};

// probably want more like 65536. see if we can get away with more restricted range too.
alignas(cache_line_size)
const AtanTable<float, float, IRange<-96, 96>, 65536, InterpolationType::Linear> atanTable_f {};
alignas(cache_line_size)
const AtanTable<double, double, IRange<-96, 96>, 65536, InterpolationType::Linear> atanTable_d {};

// not sure about range or resolution needed.
alignas(cache_line_size)
const ExpTable<float, float, IRange<-10, 10>, 1024, InterpolationType::Rounded> expTable_f {};
alignas(cache_line_size)
const ExpTable<double, double, IRange<-10, 10>, 1024, InterpolationType::Rounded> expTable_d {};

}
