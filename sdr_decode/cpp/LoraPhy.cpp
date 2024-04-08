
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

std::tuple<long, int, int> LoraPhy::detect_preamble(const cpvarray_t& sig, const long pos, const bool invert) {
    if(pos < 0) {
        printf("LoraPhy::detect_preamble - invalid offset (pos = %ld)\n", pos);
        return {-1, -1, -1};
    }

    int det_count = 1;  // current sample == 1
    int fbin_last = 0;
    long pos_offs = 0;
    for(long ii=pos; ii < sig.size()-this->sps; ii += this->sps) {
        const int fbin = this->dechirp(sig, ii, invert).first;
        //printf("%4ld) detect_preamble - fbin: %d\n", ii, fbin);

        if(abs(fbin - fbin_last) <= this->ft_det_bins) {
            det_count++;
            if(det_count >= this->plen) {
                // preamble detected, adjust fft bin to zero
                if(invert) {
                    pos_offs = ii + fbin;              // inverted chirp, non-inverted IQ
                } else {
                    pos_offs = ii + this->sps - fbin;  // non-inverted chirp, inverted IQ
                }

                // read network id
                int fbin1 = this->dechirp(sig, pos_offs+this->sps, invert).first;
                int fbin2 = this->dechirp(sig, pos_offs+2*this->sps, invert).first;
                if(invert) {
                    fbin1 = this->sps - fbin1;  // inverted chirp, non-inverted IQ
                    fbin2 = this->sps - fbin2;
                }
                // non-inverted chirp, inverted IQ needs no adjustment

                const int netid1 = fbin1 / 2;
                const int netid2 = fbin2 / 2;

                //printf("  *** preamble detected ***  fbin:%4d  pos_offs:%ld  netid1:%d  netid2:%d\n", fbin, pos_offs, netid1, netid2);
                return {pos_offs, netid1, netid2};
            }
        } else {
            det_count = 1;
        }

        fbin_last = fbin;
    }

    return {-1, -1, -1};  // preamble not detected
}

// look for the start of frame delimiter
sfdinfo_t LoraPhy::detect_sfd(const cpvarray_t& sig, const long pos, const bool invert) {
    if(pos < 0) {
        printf("LoraPhy::detect_sfd - invalid offset (pos = %ld)\n", pos);
        return sfdinfo_t(-1, -1);
    }

    for(long ii=pos; ii < sig.size()-this->sps; ii += this->sps) {
        int mval_up = dechirp(sig, ii, invert).second;
        int mval_dn = dechirp(sig, ii, !invert).second;
        if(mval_dn > mval_up) {
            // downchirp detected
            // look for second downchirp at 1.25 sps
            mval_up = dechirp(sig, ii + 1.25 * this->sps, invert).second;
            mval_dn = dechirp(sig, ii + 1.25 * this->sps, !invert).second;
            if(mval_dn > mval_up) {
                // second downchirp detected
                const long sfd = ii;
                const long hdr = sfd + 2.25 * this->sps;  // the sfd is 2.25 symbols long
                return sfdinfo_t(sfd, hdr);
            }
        }
    }

    return sfdinfo_t(-1, -1);  // sfd not detected
}

bool LoraPhy::decode_header(const cpvarray_t& sig, const long pos, const bool invert) {
    uint32_t symbols[8];
    for(int ii=0; ii<8; ii++) {
        int fbin = this->dechirp(sig, pos+ii*this->sps, invert).first;
        printf("%4d) decode_header - fbin: %d\n", ii, fbin);

        if(invert) {
            fbin = this->sps - fbin;  // inverted chirp, non-inverted IQ
        }
        // non-inverted chirp, inverted IQ needs no adjustment

        symbols[ii] = fbin / 2;
        printf("symbols[%d]: %d\n", ii, symbols[ii]);
    }

    // gray decoding
    uint32_t symbols_g[8];
    for(int ii=0; ii<8; ii++) {
        symbols_g[ii] = (symbols[ii] >> 2) ^ (symbols[ii] >> 3);
        printf("symbols_g[%d]: %d\n", ii, symbols_g[ii]);
    }

    // deinterleave
    const uint32_t bits = this->sf-2;
    uint32_t codewords[bits];
    diag_deinterleave(symbols_g, bits, codewords);
    for(int ii=0; ii<bits; ii++) {
        printf("%4d) cw: %d\n", ii, codewords[ii]);
    }

    return false;
}

// void LoraPhy::gray_decode(const uint32_t symbols[8], uint32_t* symbols_g) {
//     for(int ii=0; ii<8; ii++) {
//         symbols_g[ii] = (symbols[ii] >> 2) ^ (symbols[ii] >> 3);
//         printf("symbols_g[%d]: %d\n", ii, symbols_g[ii]);
//     }
// }

void LoraPhy::diag_deinterleave(const uint32_t symbols_g[8], const uint32_t bits, uint32_t* codewords) {
    std::fill(codewords, codewords+bits, 0);
    for(int ii=0; ii<8; ii++) {                 // 154 0 163 92 0
        // |0  0  1  0  0|: 0  0  1  0  0  <<0  ->  0  0  1  0  0
        //  0 |1  0  1  0 : 0| 1  0  1  0  <<1  ->  1  0  1  0  0
        //  1  0 |0  0  0 : 1  0| 0  0  0  <<2  ->  0  0  0  1  0
        //  0  1  0 |1  0 : 0  1  0| 1  0  <<3  ->  1  0  0  1  0
        //  0  0  1  0 |1 : 0  0  1  0| 1  <<4  ->  1  0  0  1  0
        // |0  0  1  0  0|: 0  0  1  0  0  <<0  ->  0  0  1  0  0
        //  0 |0  0  0  1 : 0| 0  0  0  1  <<1  ->  0  0  0  1  0
        //  0  0 |1  0  1 ; 0  0| 1  0  1  <<2  ->  1  0  1  0  0
        const uint32_t sym = ((symbols_g[ii]<<bits | symbols_g[ii]) << (ii%bits)) >> bits;
        for(int jj=0; jj<bits; jj++) {
            codewords[jj] |= (((sym >> jj) & 1) << ii);
        }
    }
}

chirpval_t LoraPhy::dechirp(const cpvarray_t& sig, const long pos, const bool invert) {
    const cpvarray_t& chp = invert ? this->upchirp : this->downchirp;

    // copy the chirped signal into a 2x buffer for fft proc
    const int fft_len = 2 * this->sps;
    cpvarray_t buf(fft_len);
    for(int ii=0; ii < this->sps; ii++) {
        buf[ii] = chp[ii] * sig[pos + ii];
    }

    fft(buf, fft_len);
    //printf("fft %ld:\n", pos / this->sps);
    //LoraPhy::print_array(buf, 16);

    // find max fft bin
    int mbin = 0;
    complex_t::value_type mabs = 0.0;
    for(int ii=0; ii < this->sps; ii++) {
        const complex_t::value_type abs = std::abs(buf[ii]) + std::abs(buf[this->sps + ii]);
        if(abs > mabs) {
             mbin = ii;
             mabs = abs;
         }
    }
    //printf("%ld) max bin is: %d  (val: %f)\n", pos, mbin, mabs);

    return chirpval_t(mbin, static_cast<int>(mabs));
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

int LoraPhy::load(const std::string filename, cpvarray_t& outbuf, const bool swap_iq)
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

    if(swap_iq) {
        // swap I and Q channels
        outbuf.resize(count);
        for(long ii=0; ii<count; ii++) {
            outbuf[ii] = complex_t(buf[ii].imag(), buf[ii].real());
        }
    } else {
        outbuf.swap(buf);
    }

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
    const long count = (num==0 || num>buf.size()) ? buf.size() : num;
    for(long ii=0; ii < count; ++ii) {
         printf("  %14.10f %+14.10fi  [%11.8f]\n", buf[ii].real(), buf[ii].imag(), std::abs(buf[ii]));
    }
    printf("  size: %zu\n", buf.size());
}
