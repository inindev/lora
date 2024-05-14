
#include "types.h"


namespace
{
    const complex_t::value_type tau = 2.0 * std::acos(-1);

    int my_log2(const int v)
    {
        int k = v, i = 0;
        while(k) {
            k >>= 1;
            i++;
        }
        return i - 1;
    }

    int rev_bits(const int mb, const int v)
    {
        const int l2mb = my_log2(mb);
        int p = 0;
        for(int i=1; i<=l2mb; i++) {
            if(v & (1 << (l2mb - i))) {
                p |= 1 << (i - 1);
            }
        }
        return p;
    }

    void dist_arr(cpvarray_t& buf, const int nfft)
    {
        cpvarray_t tmp(nfft);
        for(int i=0; i<nfft; i++) {
            tmp[i] = buf[rev_bits(nfft, i)];
        }
        for(int i=0; i<nfft; i++) {
            buf[i] = tmp[i];
        }
    }
}

void fft(cpvarray_t& buf, const int nfft)
{
    dist_arr(buf, nfft);

    cpvarray_t W(nfft/2);
    W[0] = 1;
    W[1] = std::polar(complex_t::value_type(1.0), -tau / nfft);
    for(int i=2, imax=nfft/2; i<imax; i++) {
        W[i] = pow(W[1], i);
    }

    int n = 1;
    int a = nfft / 2;
    for(int j=0, jmax=my_log2(nfft); j<jmax; j++) {
        for(int i=0; i<nfft; i++) {
            if(!(i & n)) {
                complex_t tmp1 = buf[i];
                complex_t tmp2 = W[(i * a) % (n * a)] * buf[i + n];
                buf[i] = tmp1 + tmp2;
                buf[i + n] = tmp1 - tmp2;
            }
        }
        n *= 2;
        a = a / 2;
    }
}
