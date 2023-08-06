# nvs_libraries

nvs_libraries is a collection of libraries intended for realtime audio dsp. 
this collection grows according to my personal needs, and it is in no way exhaustive.

nvs_filters:
audio filters using a variety of approaches, such as state-space, topology-preserving transform (TPT, trapezoidal integration).
so far, all filters except for perhaps the 'naive' (which uses Euler backward difference integration) and utility (dcBlock) are intended 
as time-variant filters, handling fast modulation well.

nvs_delayFilters:
essentially delay lines. considered as filters because they are equivalent to comb filters, and an allpass delay is also included.

nvs_gen:
some functions/objects adapted from the gen environment in max/msp.

nvs_LFO:
a simple low frequency oscillator with basic waveshapes, simultaneously available but not computed if unwanted.

nvs_logpow:
lookup tables for logarithm and pow.

nvs_matrix2x2:
a simple 2x2 matrix used currently for state space filters. implemented because it was much faster than eigen's (much more flexible) matrices.

nvs_memoryless:
memoryless mapping functions. includes also some (less customizable for now) lookup tables such as sine, tan, and tanh.

interpolation: 
some methods of hermite interpolation, for now. taken from musicdsp.org
