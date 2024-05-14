
// adapted from kissfft: https://github.com/mborgerding/kissfft

#include "kfft.h"


kfft::kfft(void)
  : nfft(0), inverse(0), factors{0}, twiddles(nullptr) { }

kfft::~kfft(void) {
    if(this->twiddles != nullptr) ::free(this->twiddles);
}

int kfft::init(int nfft, int inverse) {
    this->nfft = nfft;
    this->inverse = inverse;

    if(this->twiddles != nullptr) ::free(this->twiddles);
    this->twiddles = reinterpret_cast<complex_t*>(::malloc(sizeof(complex_t) * nfft));
    if(this->twiddles == nullptr) {
        fprintf(stderr, "error: memory allocation failure\n");
        return 1;
    }

    const complex_t::value_type tau = 2.0 * std::acos(-1);
    if(this->inverse) {
        for(int i=0; i<this->nfft; ++i) this->twiddles[i] = std::polar(complex_t::value_type(1.0),  (tau * i) / this->nfft);
    } else {
        for(int i=0; i<this->nfft; ++i) this->twiddles[i] = std::polar(complex_t::value_type(1.0), -(tau * i) / this->nfft);
    }

    // facbuf is populated by p1,m1,p2,m2, ...
    //   where
    //     p[i] * m[i] = m[i-1]
    //     m0 = n
    const complex_t::value_type floor_sqrt = std::floor(std::sqrt(complex_t::value_type(this->nfft)));

    // factor out powers of 4, powers of 2, then any remaining primes
    int* facbuf = this->factors;
    for(int n=this->nfft, p=4; n>1; ) {
        while(n % p) {
            switch(p) {
                case 2: p = 3; break;
                case 4: p = 2; break;
                default: p += 2; break;
            }
            if(p > floor_sqrt) p = n; // no more factors, skip to end
        }

        n /= p;
        *facbuf++ = p;
        *facbuf++ = n;
    }

    return 0;
}

// nfft specified by output buffer size
void kfft::operator()(const cpvarray_t& fin, cpvarray_t& fout) {
    if(this->nfft != fout.size()) init(fout.size());
    kf_work(std::begin(fout), std::begin(fin), std::end(fin), 1, 1, this->factors);
}


//
// private
//

void kfft::kf_bfly2(complex_t* fout, const int fstride, const int m) {
    complex_t* tw1 = this->twiddles;

    complex_t* fout2 = fout + m;

    complex_t t;
    for(int k=m; k>0; --k) {
        t = (*fout2) * (*tw1);

        tw1 += fstride;

        (*fout2) = (*fout) - t;

        (*fout) += t;

        ++fout2;
        ++fout;
    }
}

void kfft::kf_bfly3(complex_t* fout, const int fstride, const int m) {
    const int m2 = 2 * m;
    complex_t epi3 = this->twiddles[fstride * m];
    complex_t* tw1 = this->twiddles;
    complex_t* tw2 = this->twiddles;

    complex_t scratch[5];
    for(int k=m; k>0; --k) {
        scratch[1] = fout[m] * (*tw1);
        scratch[2] = fout[m2] * (*tw2);
        scratch[3] = scratch[1] + scratch[2];
        scratch[0] = scratch[1] - scratch[2];

        tw1 += fstride;
        tw2 += fstride * 2;

        fout[m] = (*fout) - (complex_t::value_type(0.5) * scratch[3]);

        scratch[0] *= epi3.imag();

        (*fout) += scratch[3];

        fout[m2].real(fout[m].real() + scratch[0].imag());
        fout[m2].imag(fout[m].imag() - scratch[0].real());

        fout[m].real(fout[m].real() - scratch[0].imag());
        fout[m].imag(fout[m].imag() + scratch[0].real());

        ++fout;
    }
}

void kfft::kf_bfly4(complex_t* fout, const int fstride, const int m) {
    const int m2 = 2 * m;
    const int m3 = 3 * m;
    complex_t* tw1 = this->twiddles;
    complex_t* tw2 = this->twiddles;
    complex_t* tw3 = this->twiddles;

    complex_t scratch[6];
    for(int k=m; k>0; --k) {
        scratch[0] = fout[m] * (*tw1);
        scratch[1] = fout[m2] * (*tw2);
        scratch[2] = fout[m3] * (*tw3);
        scratch[3] = scratch[0] + scratch[2];
        scratch[4] = scratch[0] - scratch[2];
        scratch[5] = (*fout) - scratch[1];

        (*fout) += scratch[1];

        fout[m2] = (*fout) - scratch[3];

        tw1 += fstride;
        tw2 += fstride * 2;
        tw3 += fstride * 3;

        (*fout) += scratch[3];

        if(this->inverse) {
            fout[m].real(scratch[5].real() - scratch[4].imag());
            fout[m].imag(scratch[5].imag() + scratch[4].real());
            fout[m3].real(scratch[5].real() + scratch[4].imag());
            fout[m3].imag(scratch[5].imag() - scratch[4].real());
        } else {
            fout[m].real(scratch[5].real() + scratch[4].imag());
            fout[m].imag(scratch[5].imag() - scratch[4].real());
            fout[m3].real(scratch[5].real() - scratch[4].imag());
            fout[m3].imag(scratch[5].imag() + scratch[4].real());
        }

        ++fout;
    }
}

void kfft::kf_bfly5(complex_t* fout, const int fstride, const int m) {
    complex_t* twiddles = this->twiddles;
    complex_t ya = twiddles[fstride * m];
    complex_t yb = twiddles[fstride * 2 * m];

    complex_t* fout0 = fout;
    complex_t* fout1 = fout0 + m;
    complex_t* fout2 = fout0 + 2 * m;
    complex_t* fout3 = fout0 + 3 * m;
    complex_t* fout4 = fout0 + 4 * m;

    complex_t* tw = this->twiddles;
    complex_t scratch[13];
    for(int u=0; u<m; ++u) {
        scratch[0] = *fout0;
        scratch[1] = (*fout1) * tw[u*fstride];
        scratch[2] = (*fout2) * tw[2*u*fstride];
        scratch[3] = (*fout3) * tw[3*u*fstride];
        scratch[4] = (*fout4) * tw[4*u*fstride];

        scratch[7] = scratch[1] + scratch[4];
        scratch[10] = scratch[1] - scratch[4];
        scratch[8] = scratch[2] + scratch[3];
        scratch[9] = scratch[2] - scratch[3];

        (*fout0) += scratch[7] + scratch[8];

        scratch[5] = scratch[0] + scratch[7] * ya.real() + scratch[8] * yb.real();
        scratch[6].real( scratch[10].imag() * ya.imag() + scratch[9].imag() * yb.imag());
        scratch[6].imag(-scratch[10].real() * ya.imag() - scratch[9].real() * yb.imag());

        (*fout1) = scratch[5] - scratch[6];
        (*fout4) = scratch[5] + scratch[6];

        scratch[11] = scratch[0] + scratch[7] * yb.real() + scratch[8] * ya.real();
        scratch[12].real(-scratch[10].imag() * yb.imag() + scratch[9].imag() * ya.imag());
        scratch[12].imag( scratch[10].real() * yb.imag() - scratch[9].real() * ya.imag());

        (*fout2) = scratch[11] + scratch[12];
        (*fout3) = scratch[11] - scratch[12];

        ++fout0;
        ++fout1;
        ++fout2;
        ++fout3;
        ++fout4;
    }
}

// perform the butterfly for one stage of a mixed radix fft
void kfft::kf_bfly_generic(complex_t* fout, const int fstride, const int m, const int p) {
    complex_t* twiddles = this->twiddles;
    int Norig = this->nfft;

    complex_t scratch[p];
    complex_t tmp;
    for(int u=0, k=0; u<m; ++u) {
        k = u;
        for(int q1=0; q1<p; ++q1) {
            scratch[q1] = fout[k];
            k += m;
        }

        k = u;
        for(int q1=0; q1<p; ++q1) {
            int twidx = 0;
            fout[k] = scratch[0];

            for(int q=1; q<p; ++q) {
                twidx += fstride * k;
                if(twidx >= Norig) twidx -= Norig;
                tmp = scratch[q] * twiddles[twidx];
                fout[k] += tmp;
            }

            k += m;
        }
    }
}

void kfft::kf_work(complex_t* fout, const complex_t* fin, const complex_t* fin_end, const int fstride, const int in_stride, const int* factors) {
    complex_t* fout_beg = fout;
    const int p = *factors++; // the radix
    const int m = *factors++; // stage's fft length/p
    const complex_t* fout_end = fout + p * m;

    if(m == 1) {
        do{
            *fout = (fin < fin_end) ? *fin : 0;
            fin += fstride * in_stride;
        } while(++fout != fout_end);
    } else {
        do {
            // recursive call:
            // DFT of size m*p performed by doing
            // p instances of smaller DFTs of size m,
            // each one takes a decimated version of the input
            kf_work(fout, fin, fin_end, fstride* p, in_stride, factors);
            fin += fstride * in_stride;
        } while( (fout += m) != fout_end );
    }

    fout = fout_beg;

    // recombine the p smaller DFTs
    switch(p) {
        case 2: kf_bfly2(fout, fstride, m); break;
        case 3: kf_bfly3(fout, fstride, m); break;
        case 4: kf_bfly4(fout, fstride, m); break;
        case 5: kf_bfly5(fout, fstride, m); break;
        default: kf_bfly_generic(fout, fstride, m, p); break;
    }
}

