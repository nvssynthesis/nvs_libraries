#include "../include/nvs_gen.h"
#include <iostream>
#include <vector>
#include <string>

using namespace nvs::gen;

void change_test(){
	change<float> c;
	std::vector<float> input = {0.f, 0.1f, 0.01f, -100.f, -0.5f};
	std::vector<float> expected_output = {
		0.f, 1.f, -1.f, -1.f, 1.f
	};
	for (auto i = 0; i < input.size(); ++i){
		auto out = c(input[i]);
		assert(out == expected_output[i]);
	}
	std::cout << "change_test success\n";
}

std::vector<float> makePhasorWave(size_t length, float freq, float fs){
	std::vector<float> wave(length);
	phasor ph(fs);
	ph.setFrequency(freq);
	for (auto &e : wave)
		e = (ph++).getPhase();
	assert(wave.size() == length);
	return wave;
}

void phasor_test(){
	constexpr size_t N = 10;
	std::vector<float> ph = makePhasorWave(N, 22050.f, 44100.f);
	std::vector<float> expected_out = {
		0.f, 0.5f, 0.f, 0.5f, 0.f, 0.5f, 0.f, 0.5f, 0.f, 0.5f
	};
	size_t out_idx = 0;
	for (auto &e : ph)
		assert(e == expected_out[out_idx++]);
	std::cout << "phasor test nyquist success\n";
	ph = makePhasorWave(10, 0.f, 44100.f);
	expected_out = {
		0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f
	};
	std::cout << "ABOUT to ACTUALLY phasor test\n";
	assert(ph.size() == expected_out.size());
	assert(ph.size() == N);
	out_idx = 0;
	for (float &val : ph){
		float expected = expected_out[out_idx++];
		if (val != expected){
			std::cout << "expected: " << expected << " val: " << val << '\n';
		}
		assert(val == expected);
	}

	std::cout << "phasor test DC success\n";
}

void ramp2trig_test(){
	size_t len = 10;
	auto wave = makePhasorWave(len, 24000.f, 48000.f);
	std::vector<float> expected = {
		0.f, 1.f, 0.f, 1.f, 0.f, 1.f, 0.f, 1.f, 0.f, 1.f
	};

	ramp2trig<float> r2t;

	size_t idx = 0;
	for (auto &e : wave){
		auto expect = expected[idx++];
		auto out = r2t(e);
		// std::cout << "res: " << out << " expected: " << expect << '\n';
		assert(out == expect);
	}
	std::cout << "ramp2trig test success\n";
}

void switcher_test(){
	std::vector<int> v1 = 		{1,2,3,4,5,6,7,8,9};
	std::vector<int> v2 = 		{9,8,7,6,5,4,3,2,1};
	std::vector<int> control =  {0,10,0,-1,0,5,0,-100,0};	// converts to bool
	std::vector<int> expected = {9,2,7,4,5,6,3,8,1};
	size_t idx = 0;
	for (; idx < v1.size(); ++idx){
		auto out = switcher(control[idx], v1[idx], v2[idx]);
		assert(out == expected[idx]);
	}
	std::cout << "switcher test success\n";
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

void gateSelect_test(){
	constexpr size_t N = 3;

	float val = 10.f;
	for (size_t i = 0; i < N; ++i){
		std::array<float, N> res = gateSelect<float, N> (val, i);
		// print(res, "\n");
		print(i, "   i\n");

		if (i > 0) {
			for (int j = 0; j < N; ++j){
				print(j, "     j\n");
				if (j == (i - 1)){
					assert(res[j] == val);
				}
				else {
					assert(res[j] == 0.f);
				}
			}
		} 
		else {
			for (int j = 0; j < N; ++j){
				assert(res[j] == 0);
			}
		}
	}
	std::cout << "gate select test success\n"; 
}

void latch_test(){
	phasor ph(10.f);	// sample rate of 10
	ph.setFrequency(2.f);
	constexpr size_t N = 10;
	std::array<float, N> phRec;
	for (auto &e : phRec)
		e = (ph++).getPhase();

	std::array<float, N> expected = {
		0.f, 0.f, 0.4f, 0.6f, 0.8f, 0.f, 0.2f, 0.4f, 0.4f, 0.4f
	};
	latch<float> l;
	size_t i = 0;
	bool letPass = false;
	for (float e : phRec){
		auto exp = expected[i++];
		if (i == 3){
			letPass = true;
		} else if (i == 9){
			letPass = false;
		}
		auto v = l(e, letPass);

		if (v != exp)
			std::cout << "v: " << v << " expected: " << exp << '\n';
		assert(v == exp);
	}
	std::cout << "latch test success\n";
}

void parzen_demo(){
	constexpr size_t N = 20;
	std::vector<float> wave = makePhasorWave(N, 90.f, 1000.f);
	for (size_t i = 0; i < N; ++i){
		std::cout << "wave: " << wave[i] << " parzen: " << parzen<float>(wave[i]) << '\n';
	}
}

int main(void){
	switcher_test();
	parzen_demo();

	return 0;
}

