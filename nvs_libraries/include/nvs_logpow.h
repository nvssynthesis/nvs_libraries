//
//  Created by Nicholas Solem on 2/2/23.
//  Copyright © 2023 Nicholas Solem. All rights reserved.
//
#include <cmath>
#include <vector>
#include <climits>
#include <algorithm>
/*
    to use:

    typedef float f_t;  
    f_t xmin = 0.f;     // ~-inf dB, will actually be finitely bound based on eps<T>
    f_t xrange = 4.f;   // ~12 dBFS
    constexpr size_t N ( 0x1UL << 16 ); // 16 bits of floating-point precision
    logTable<f_t> ltab(xmin, xrange, N);
    
    f_t val = 0.5f;     // or any value between [xmin .. xmin+xrange] inclusive
    f_t tabResult = ltab(val);
*/
namespace nvs {
    template <typename T>
    [[deprecated]]
    class abstractBaseTable {
    public:
        abstractBaseTable(T x_min = 0.0, T x_range = 1.0, size_t _N = 8192)
        :   xmin(x_min), xrange(x_range), N(_N){}
        virtual ~abstractBaseTable() = default;
        virtual T operator()(T val) const = 0;
    protected:
        virtual void makeTable() = 0;
        virtual T readTable(T val) const = 0;
        T xmin, xrange;
        size_t N;
        std::vector<T> table;
    };
    template <typename T>
    class logTable  :   public abstractBaseTable<T> {
    public:
        logTable(T x_min = 0.0, T x_range = 1.0, size_t _N = 8192, T _base = 10.0)
        : abstractBaseTable<T>(x_min, x_range, _N), base(_base){
            makeTable();
        }
        T operator()(T val) const override {
            return readTable(val);
        }
    private:
        static constexpr T eps = std::numeric_limits<T>::epsilon();
        const T base;
        void makeTable() override {
            this->table.resize(this->N);
            for (size_t i = 0; i < this->N; i++){
                T x = this->xmin + (this->xrange * ((T)(i + eps) / (T)(this->N-1))); //
                this->table[i] = std::log10(x) / std::log10(base);
            }
        }
        T readTable(T val) const override {
            val = (val - this->xmin);
            val = std::max(val, (T)0);
            val = val * (T)(this->N-1);
            val = val / this->xrange;
            size_t idx0 = static_cast<size_t>(val);
            idx0 = std::max(idx0, (size_t)0);
            T y0 = this->table[idx0]; // return here if no linear interpolation
            size_t idx1 = std::min((idx0 + 1), this->N - 1);
            T y1 = this->table[idx1];
            T frac = val - (T)idx0;
            T y = y0 * (1.0 - frac) + y1 * (frac);
            return y;
        }
    };
    template <typename T>
    class powTable  :   public abstractBaseTable<T> {
    public:
        powTable(T x_min = 0.0, T x_range = 1.0, size_t _N = 8192, T _base = 10.0)
        : abstractBaseTable<T>(x_min, x_range, _N), base(_base) {
            makeTable();
        }
        T operator()(T val) const override {
            return readTable(val);
        }
    private:
        const T base;
        void makeTable() override {
            this->table.resize(this->N);
            for (int i = 0; i < this->N; i++){
                T x =  ((T)i) / ((T)(this->N-1));
                x = x * this->xrange;
                x = x + this->xmin;
                this->table[i] = std::pow(base, (T)x);
            }
        }
        T readTable(T val) const override {
            val = val - this->xmin;
            val = std::max(val, (T)0);
            val = val * (T)(this->N-1);
            val = val / this->xrange;
            size_t idx0 = static_cast<size_t>(val);
//            idx0 = std::max(idx0, (size_t)0);
            T y0 = this->table[idx0];
            size_t idx1 = std::min(idx0 + 1, this->N-1);
            T frac = val - (T)idx0;
            T y1 = this->table[idx1];
            return (1.0 - frac) * y0 + (frac) * y1;
        }
    };
template class logTable<float>;
template class logTable<double>;
template class powTable<float>;
template class powTable<double>;
}   // namespace nvs
