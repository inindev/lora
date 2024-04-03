
#include "types.h"


namespace
{
    long log2(long N)
    {
        long k = N, i = 0;
        while(k) {
            k >>= 1;
            i++;
        }
        return i - 1;
    }

    long rev_bits(long N, long n)
    {
        long p = 0;
        for(long ii=1; ii<=log2(N); ii++) {
            if(n & (1 << (log2(N) - ii))) {
                p |= 1 << (ii - 1);
            }
        }
        return p;
    }

    void rev_arr(cpvarray_t& f1, long N)
    {
        cpvarray_t f2(N);
        for(long ii=0; ii<N; ii++) {
            f2[ii] = f1[rev_bits(N, ii)];
        }
        for(long ii=0; ii<N; ii++) {
            f1[ii] = f2[ii];
        }
    }
}

void fft(cpvarray_t& f, long N)
{
    rev_arr(f, N);

    cpvarray_t W(N/2);
    W[0] = 1;
    W[1] = std::polar(1.0, (-2.0 * M_PI) / N);
    for(long ii= 2; ii<N/2; ii++) {
        W[ii] = pow(W[1], ii);
    }

    long n = 1;
    long a = N / 2;
    for(long jj=0; jj<log2(N); jj++) {
        for(long ii=0; ii<N; ii++) {
            if(!(ii & n)) {
                complex_t temp = f[ii];
                complex_t Temp = W[(ii * a) % (n * a)] * f[ii + n];
                f[ii] = temp + Temp;
                f[ii + n] = temp - Temp;
            }
        }
        n *= 2;
        a = a / 2;
    }
}
