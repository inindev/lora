
// copyright (C) 2024, John Clark <inindev@gmail.com>

#include "LoraPhy.h"
#include "fft.h"


LoraPhy::LoraPhy(void)
  : sf(0), bw(0), sps(0), fs(0), plen(0), ft_det_bins(0)
{
    printf("default construct\n");
}

LoraPhy::~LoraPhy() {    
    printf("final destruct\n");
}

void LoraPhy::init(const int sf, const int bw, const int plen)
{
    this->sf   = sf;
    this->bw   = bw;
    this->sps  = 2 * pow(2, sf);
    this->fs   = 2 * bw;
    this->plen = plen;      // preamble detect symbols
    this->ft_det_bins = 4;  // fft detect bin variance

    LoraPhy::chirp(this->sps, this->fs, this->bw, false, this->upchirp);
    LoraPhy::chirp(this->sps, this->fs, this->bw, true, this->downchirp);
    //LoraPhy::print_array(this->downchirp, 16);
}

long LoraPhy::detect_preamble(const cpvarray_t& sig, const long pos, const bool invert) {
    if(pos < 0) {
        printf("LoraPhy::detect_preamble - invalid offset (pos = %ld)\n", pos);
        return -1;
    }

    int det_count = 1;  // current sample == 1
    long pos_offs = 0;
    long fbin_last = 0;
    for(long ii=pos; ii < sig.size()-this->sps; ii += this->sps) {
        long fbin = this->dechirp(sig, ii, invert);
        //printf("%ld) detect_preamble - fbin: %ld\n", ii, fbin);

        if(abs(fbin - fbin_last) <= this->ft_det_bins) {
            det_count++;
            if(det_count >= this->plen) {
                // preamble detected, adjust fft bin to zero
                if(invert) {
                    pos_offs = ii + fbin;                // inverted chirp, non-inverted IQ
                } else {
                    pos_offs = ii + (this->sps - fbin);  // non-inverted chirp, inverted IQ
                }
                //printf("start:%ld  fbin:%ld  pos_offs:%ld  *** preamble detected ***\n", pos, fbin, pos_offs);
                return pos_offs;
            }
        } else {
            det_count = 1;
        }

        fbin_last = fbin;
    }

    return -1;  // preamble not detected
}

long LoraPhy::dechirp(const cpvarray_t& sig, const long pos, const bool invert) {
    const cpvarray_t& chp = invert ? this->upchirp : this->downchirp;

    // copy the chirped signal into a 2x buffer for fft proc
    const long fft_len = 2 * this->sps;
    cpvarray_t buf(fft_len);
    for(long ii=0; ii < this->sps; ii++) {
        buf[ii] = chp[ii] * sig[pos + ii];
    }

    fft(buf, fft_len);
    //printf("fft:\n");
    //LoraPhy::print_array(buf, 16);

    // find max fft bin
    long mbin = 0;
    complex_t::value_type mabs = 0.0;
    for(long ii=0; ii < this->sps; ii++) {
        const complex_t::value_type abs = std::abs(buf[ii]) + std::abs(buf[this->sps + ii]);
        if(abs > mabs) {
             mbin = ii;
             mabs = abs;
         }
    }
    //printf("%ld) ** max bin is: %ld  (val: %f)\n", pos, mbin, mabs);

    return mbin;
}

void LoraPhy::chirp(const int sps, const int fs, const int bw, const bool is_downchirp, cpvarray_t& outbuf) {
    const complex_t::value_type f0 = -bw/2.0 * (is_downchirp ? -1 : 1);
    const complex_t::value_type f1 =  bw/2.0 * (is_downchirp ? -1 : 1);
    const complex_t::value_type t1 =  (sps-1) / static_cast<complex_t::value_type>(fs);

    const complex_t::value_type tau = 2.0*std::acos(-1);
    const complex_t::value_type beta = (f1 - f0) / (2.0*t1);

    cpvarray_t buf(sps);
    for(int ii=0; ii<sps; ++ii) {
        const complex_t::value_type t = ii / static_cast<complex_t::value_type>(fs);
        buf[ii] = complex_t(
             cos( tau * (beta * t*t + f0*t) ),
            -cos( tau * (beta * t*t + f0*t + 90.0/360.0) )
        );
    }

    outbuf.swap(buf);
}

int LoraPhy::load(const std::string filename, cpvarray_t& outbuf)
{
    FILE* pf = fopen(filename.c_str(), "rb");
    if(pf == NULL) {
        printf("LoraPhy::load - failed to open binary file for reading\n");
        return -1;
    }

    int rc = fseek(pf, 0L, SEEK_END);
    if(rc != 0) {
        printf("LoraPhy::load - fseek(SEEK_END) failed, rc=%d\n", rc);
        return rc;
    }
    const long filebytes = ftell(pf);

    rc = fseek(pf, 0L, SEEK_SET);
    if(rc != 0) {
        printf("LoraPhy::load - fseek(SEEK_SET) failed, rc=%d\n", rc);
        return rc;
    }

    const long count = (filebytes / sizeof(cpvarray_t::value_type));
    cpvarray_t buf(count);
    const size_t ret = fread(&buf[0], sizeof(cpvarray_t::value_type), count, pf);
    fclose(pf);
    if(ret != count) {
        printf("LoraPhy::load - failed to read binary file (%zu of %ld read)\n", ret, count);
        return -2;
    }

    outbuf.swap(buf);

    return 0;  // success
}

int LoraPhy::save(const std::string filename, const cpvarray_t& buf) {
    FILE* pf = fopen(filename.c_str(), "wb");
    if(pf == NULL) {
        printf("LoraPhy::save - failed to open binary file for writing\n");
        return -1;
    }

    const size_t ret = fwrite(&buf[0], sizeof(cpvarray_t::value_type), buf.size(), pf);
    fclose(pf);
    if(ret != buf.size()) {
        printf("LoraPhy::save - failed to write binary file (%zu of %lu written)\n", ret, buf.size());
        return -1;
    }

    return 0;  // success
}

void LoraPhy::print_array(const cpvarray_t& buf, const long num) {
    const long count = (num==0) ? buf.size() : num;
    for(long ii=0; ii < count; ++ii) {
         printf("  %14.10f %+14.10fi  [%11.8f]\n", buf[ii].real(), buf[ii].imag(), std::abs(buf[ii]));
    }
    printf("  size: %zu\n", buf.size());
}
