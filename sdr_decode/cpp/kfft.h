
#pragma once

#include "types.h"


struct kfft
{
    kfft(void);
    virtual ~kfft(void);

    int init(int nfft, int inverse=0);

    // nfft specified by output buffer size
    void operator()(const cpvarray_t& fin, cpvarray_t& fout);

private:
    int nfft;
    int inverse;
    int factors[2*32];
    complex_t* twiddles;

    void kf_bfly2(complex_t* fout, const int fstride, const int m);
    void kf_bfly3(complex_t* fout, const int fstride, const int m);
    void kf_bfly4(complex_t* fout, const int fstride, const int m);
    void kf_bfly5(complex_t* fout, const int fstride, const int m);
    void kf_bfly_generic(complex_t* fout, const int fstride, const int m, const int p);

    void kf_work(complex_t* fout, const complex_t* fin, const complex_t* fin_end, const int fstride, const int in_stride, const int* factors);
};

