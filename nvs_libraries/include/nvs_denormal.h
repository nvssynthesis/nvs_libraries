/* nvs_denormal.h */
// derived from https://www.musicdsp.org/en/latest/_downloads/d81fc9af8a9fa63332b248772fbb4a54/denormal.pdf
// templated rather than float,
// BUT the precisions here are intended for float (not double or long double)
// TODO: figure out double and long double versions instead

namespace nvs_denormal {
unsigned int rand_state = 1; 

template <typename T>
add_white_noise (T &val)
{
	rand_state = rand_state * 123456UL + 890123UL;
	int mantissa = rand_state & 0x807F0000; // maintain significant bits
	// int mantissa = rand_state & 0x007F0000; // generate positive-only
	int flt_rand = mantissa | 0x1E000000;
	val += *reinterpret_cast <const T *> (&flt_rand);
}

template <typename T>
add_white_noise_buffer(T val_arr [], const T noise_arr[], long n_samps)
{ // requires write of noise to buffer but will then be more efficient
	for (long pos = 0; pos < n_samps; ++pos)
	{
		val_arr[pos] += noise_arr[pos];
	}
}

template <typename T>
void add_dc(T &val) {
	static const T anti_denormal = 1e-20;
	val += anti_denormal;
}

}