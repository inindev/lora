
// copyright (C) 2024, John Clark <inindev@gmail.com>

#include "LoraPhy.h"
#include "fft.h"


LoraPhy::LoraPhy(void)
  : sf(0), sr(0), bw(0), sps(0), fs(0), cr(0), use_ldro(0),
    use_crc(0), use_hamming(0), has_header(0), plen(0),
    ft_ratio(0), ft_sym_fct(0), ft_bins(0), get_sample_fcn(nullptr) { }

LoraPhy::~LoraPhy() { }

int LoraPhy::init(const std::string& filename, const int sample_bits, const bool swap_iq,
                  const int sf, const int bw, const int plen) {
    this->sf          = sf;
    this->sr          = 2;
    this->bw          = bw;
    this->sps         = this->sr * pow(2, sf);
    this->fs          = this->sr * bw;
    this->cr          = 1; // 4/5
    this->use_ldro    = 0;
    this->use_crc     = 1;
    this->use_hamming = 1;
    this->has_header  = 1;
    this->plen        = plen; // preamble detect symbols
    this->ft_ratio    = 2;    // fft bins per sps
    this->ft_sym_fct  = this->sr * this->ft_ratio;  // fft symbol factor
    this->ft_bins     = this->ft_ratio * this->sps; // fft bins per symbol

    LoraPhy::chirp(this->sps, this->fs, this->bw, false, this->upchirp);
    LoraPhy::chirp(this->sps, this->fs, this->bw, true, this->downchirp);

    ifs.open(("-"==filename) ? "/dev/stdin" : filename);
    if(!ifs) {
        fprintf(stderr, "error: LoraPhy::init - unable to open input source - filename: %s\n", filename.c_str());
        return 2;
    }

    if(sample_bits == 8) {
        this->get_sample_fcn = swap_iq ? sdr_stream::cu8_qi_sample : sdr_stream::cu8_iq_sample;
    }
    else if(sample_bits == 32) {
        this->get_sample_fcn = swap_iq ? sdr_stream::cf32_qi_sample : sdr_stream::cf32_iq_sample;
    }
    else {
        fprintf(stderr, "error: LoraPhy::init - invalid bits per sample setting: %d\n", sample_bits);
        return 3;
    }

    return 0;
}

std::tuple<int, int, int> LoraPhy::detect_preamble(const bool invert) {
    int det_count = 1;  // current sample == 1
    int fbin_last = 0;
    int pos_adj = 0;
    //int ctr = 0;
    cpvarray_t sig(this->sps);
    for(int rc=0;;) {
        //if(ctr++ % 1000 == 0) printf("kloop# %d\n", ctr/1000);
        if((rc = get_sample(sig, pos_adj))) return {-rc, -1, -1};
        const chirpval_t cv = this->dechirp(sig, 0, invert);
        int fbin = cv.first;
        const int mval = cv.second;
        //printf("detect_preamble(%d) - fbin: %d  mval: %d\n", det_count, fbin, mval);

        int dfbin = std::abs(fbin - fbin_last);
        if((mval > 15) && ((dfbin <= this->ft_ratio) || (dfbin >= (this->ft_bins - this->ft_ratio)))) {
            det_count++;
            pos_adj = 0 - (fbin / static_cast<float>(this->ft_ratio));
            if(pos_adj < -16) {
                fbin = 0;
            } else {
                pos_adj = 0;
            }

            if(det_count >= this->plen) {
                // preamble detected, adjust fft bin to zero
                if(invert) {
                    pos_adj = fbin;                  // inverted chirp, non-inverted IQ
                } else {
                    pos_adj = 0 - (fbin / static_cast<float>(this->ft_ratio));  // non-inverted chirp, inverted IQ
                }

                // consume any extra preamble chirps while looking for network id 1
                int fbin1 = 0;
                while((dfbin <= this->ft_ratio) || (dfbin >= (this->ft_bins - this->ft_ratio))) {
                    // read network id
                    if((rc = get_sample(sig, pos_adj))) return {-rc, -1, -1};
                    fbin1 = this->dechirp(sig, 0, invert).first;
                    dfbin = abs(fbin1 - fbin_last);
                }

                if((rc = get_sample(sig))) return {-rc, -1, -1};
                int fbin2 = this->dechirp(sig, 0, invert).first;

                if(invert) {
                    fbin1 = this->sps - fbin1;  // inverted chirp, non-inverted IQ
                    fbin2 = this->sps - fbin2;
                }
                // non-inverted chirp, inverted IQ needs no adjustment

                const int netid1 = std::round(fbin1 / static_cast<float>(this->ft_sym_fct));
                const int netid2 = std::round(fbin2 / static_cast<float>(this->ft_sym_fct));

                //printf("\n  *** preamble detected ***  fbin:%4d  pos_adj:%d  netid1:%d  netid2:%d\n", fbin, pos_adj, netid1, netid2);
                return {0, netid1, netid2};
            }
        } else {
            det_count = 1;
        }

        fbin_last = fbin;
    }

    return {-1, -1, -1};  // preamble not detected
}

// look for the start of frame delimiter
sfdinfo_t LoraPhy::detect_sfd(const bool invert) {
    // scan up to five samples looking for the sfd
    cpvarray_t sig(this->sps);
    for(int i=0; i<5; i++) {
        if(get_sample(sig)) return sfdinfo_t(-1, -1);
        int mval_up = dechirp(sig, 0, invert).second;
        int mval_dn = dechirp(sig, 0, !invert).second;
        if(mval_dn > mval_up) {
            // downchirp detected
            // look for second downchirp at 1.25 sps
            if(get_sample(sig, 0.25 * this->sps)) return sfdinfo_t(-1, -1);
            mval_up = dechirp(sig, 0, invert).second;
            mval_dn = dechirp(sig, 0, !invert).second;
            if(mval_dn > mval_up) {
                // second downchirp detected
                const int sfd = 0;
                const int hdr = 2.25 * this->sps;  // the sfd is 2.25 symbols long
                return sfdinfo_t(sfd, hdr);
            }
        }
    }

    return sfdinfo_t(-1, -1);  // sfd not detected
}


// <payload_len, cr, crc, is_valid>
std::tuple<uint8_t, uint8_t, uint8_t, bool> LoraPhy::decode_header(const bool invert) {
    cpvarray_t sig(8 * this->sps);
    if(get_sample(sig)) return {0, 0, 0, false};

    u16varray_t symbols(8);
    for(int i=0; i<8; i++) {
        int fbin = this->dechirp(sig, i*this->sps, invert).first;
        // printf("%4d) decode_header - fbin: %d\n", i, fbin);

        if(invert) {
            fbin = this->sps - fbin;  // inverted chirp, non-inverted IQ
        }
        // non-inverted chirp, inverted IQ needs no adjustment

        symbols[i] = std::round(fbin / static_cast<float>(this->ft_sym_fct));
        //printf("symbols[%d]: %d\n", i, symbols[i]);
    }

    // gray decoding
    const u16varray_t& symbols_g = gray_decode(symbols, true);
    // for(int i=0; i<8; i++) {
    //    printf("%4d) symg: %d\n", i, symbols_g[i]);
    // }

    // deinterleave
    const u16varray_t& codewords = diag_deinterleave(symbols_g, 4, this->sf-2); // cr, sf
    // for(int i=0; i<this->sf-2; i++) {
    //    printf("%4d) cw: %d\n", i, codewords[i]);
    // }

    // hamming decode
    const u8varray_t& header = hamming_decode(codewords, 4); // cr

    // parse header
    const uint8_t payload_len = (header[0] << 4) | header[1];
    const uint8_t crc = header[2] & 0x01;
    const uint8_t cr = header[2] >> 1;

    // validate the header checksum
    const uint8_t header_checksum = ((header[3] & 0x01) << 4) | header[4];
    const uint8_t header_checksum_calc = calc_header_csum(header);
    const bool is_valid = (header_checksum == header_checksum_calc);

    return {payload_len, cr, crc, is_valid};
}

u16varray_t LoraPhy::decode_payload(const uint8_t payload_len, const bool invert) {
    const uint8_t symcnt = calc_payload_symbol_count(payload_len);

    cpvarray_t sig(this->sps);
    u16varray_t symbols(symcnt);
    for(uint8_t i=0; i<symcnt; i++) {
        if(get_sample(sig)) {
            fprintf(stderr, "error: decode_payload - unexpected end of data reached - sym_num:%d\n", i);
            return u16varray_t();
        }

        int fbin = this->dechirp(sig, 0, invert).first;
        if(invert) {
            fbin = this->sps - fbin + 1;  // inverted chirp, non-inverted IQ
        }
        // non-inverted chirp, inverted IQ needs no adjustment

        symbols[i] = std::round(fbin / static_cast<float>(this->ft_sym_fct));
        //printf("%d) fbin:%3d  -->  sym:%3d\n", i, fbin, symbols[i]);
    }

    return symbols;
}

uint8_t LoraPhy::calc_payload_symbol_count(const uint8_t payload_len) {
    // see SX1276_Datasheet.pdf p.31
    const int symcnt = std::ceil((2*double(payload_len)-this->sf +7+ 4*this->use_crc-5*(1-this->has_header)) / (this->sf-2*this->use_ldro)) * (this->cr+4);
    if(symcnt < 0) {
        return 0;
    }

    return symcnt;
}

uint8_t LoraPhy::calc_header_csum(const u8varray_t& header) {
    // header checksum
    const uint8_t c4 = (header[0] & 0x08) >> 3 ^ (header[0] & 0x04) >> 2 ^ (header[0] & 0x02) >> 1 ^ (header[0] & 0x01);
    const uint8_t c3 = (header[0] & 0x08) >> 3 ^ (header[1] & 0x08) >> 3 ^ (header[1] & 0x04) >> 2 ^ (header[1] & 0x02) >> 1 ^ (header[2] & 0x01);
    const uint8_t c2 = (header[0] & 0x04) >> 2 ^ (header[1] & 0x08) >> 3 ^ (header[1] & 0x01) ^ (header[2] & 0x08) >> 3 ^ (header[2] & 0x02) >> 1;
    const uint8_t c1 = (header[0] & 0x02) >> 1 ^ (header[1] & 0x04) >> 2 ^ (header[1] & 0x01) ^ (header[2] & 0x04) >> 2 ^ (header[2] & 0x02) >> 1 ^ (header[2] & 0x01);
    const uint8_t c0 = (header[0] & 0x01) ^ (header[1] & 0x02) >> 1 ^ (header[2] & 0x08) >> 3 ^ (header[2] & 0x04) >> 2 ^ (header[2] & 0x02) >> 1 ^ (header[2] & 0x01);
    const uint8_t csum = (c4 << 4) | (c3 << 3) | (c2 << 2) | (c1 << 1) | c0;
    return csum;
}

u16varray_t LoraPhy::gray_decode(const u16varray_t& symbols, const bool use_ldro) {
    const int symcnt = symbols.size();
    u16varray_t symbols_g(symcnt);

    if(use_ldro) {
        for(int i=0; i<symcnt; i++) {
            const uint16_t sym_b2 = (symbols[i] >> 2);
            const uint16_t sym_b3 = (sym_b2 >> 1);
            symbols_g[i] = (sym_b2 ^ sym_b3);
        }
    } else {
        const uint16_t sfb = pow(2, this->sf);
        for(int i=0; i<symcnt; i++) {
            const uint16_t sym_b2 = ((symbols[i] + 0xffff) % sfb);
            const uint16_t sym_b3 = (sym_b2 >> 1);
            symbols_g[i] = (sym_b2 ^ sym_b3);
        }
    }

    return symbols_g;
}

u8varray_t LoraPhy::hamming_decode(const u16varray_t& codewords, const int cr) {
    const int cr_bits = cr + 4;
    u8varray_t result(codewords.size());

    // TODO: support CR 4/7, 4/6, and 4/5 (not needed to decode header)
    if(cr_bits < 8) {
        for(uint8_t i=0; i<codewords.size(); i++) {
            result[i] = codewords[i] & 0x0f;
        }
        return result;
    }

    //           parity    data
    //  cr 4/8  p3p2p1p0 d3d2d1d0
    //  cr 4/7    p2p1p0 d3d2d1d0
    //  cr 4/6      p1p0 d3d2d1d0
    //  cr 4/5        p4 d3d2d1d0
    //
    //    p0 = d0 ^ d1 ^ d2     ^ = xor
    //    p1 = d1 ^ d2 ^ d3
    //    p2 = d0 ^ d1 ^ d3
    //    p3 = d0 ^ d2 ^ d3
    //    p4 = d0 ^ d1 ^ d2 ^ d3   (CR=4/5)
    for(uint8_t i=0; i<codewords.size(); i++) {
        const uint8_t data = codewords[i] & 0x0f;
        const uint8_t parity = (codewords[i] >> 4) & 0x0f;

        // calculate parity
        const uint8_t d0 = data & 0x01;
        const uint8_t d1 = (data >> 1) & 0x01;
        const uint8_t d2 = (data >> 2) & 0x01;
        const uint8_t d3 = (data >> 3) & 0x01;

        const uint8_t p0 = d0 ^ d1 ^ d2;
        const uint8_t p1 = d1 ^ d2 ^ d3;
        const uint8_t p2 = d0 ^ d1 ^ d3;
        const uint8_t p3 = d0 ^ d2 ^ d3;
        //const uint8_t p4 = d0 ^ d1 ^ d2 ^ d3;  // (CR=4/5)

        const uint8_t pcalc = (p3 << 3) | (p2 << 2) | (p1 << 1) | p0;
        if(parity != pcalc) {
            // repair data
            const uint8_t perr = parity ^ pcalc ^ 0x0f;
            // bit 1->8, bit 8->2, bit 2->1
            const uint8_t be = ((perr & 0x01) << 3) | (perr & 0x04) | ((perr & 0x08) >> 2) | ((perr & 0x02) >> 1);
            if((be==1) || (be==2) || (be==4) | (be==8)) {
                result[i] = data ^ be;
                fprintf(stderr, "warn: parity data correction - data elem:%d  bit: %d  repaired data: %d\n", i, be, result[i]);
            } else {
                fprintf(stderr, "error: !!!! unrecoverable parity error !!!! - data elem:%d  bit: %d\n", i, be);
            }
        } else {
            result[i] = data;
        }
    }

    return result;
}

u16varray_t LoraPhy::diag_deinterleave(const u16varray_t& symbols_g, const int cr, const int sf) {
    const int cr_bits = cr + 4;
    const int cw_len = (symbols_g.size() / cr_bits) * sf;
    u16varray_t codewords(cw_len);

    int cwi = 0;
    for(int offs=0; offs<symbols_g.size(); offs+=cr_bits) {
        for(int i=0; i<cr_bits; i++) {              // 154 0 163 92 0
            // |0  0  1  0  0|: 0  0  1  0  0  <<0  ->  0  0  1  0  0
            //  0 |1  0  1  0 : 0| 1  0  1  0  <<1  ->  1  0  1  0  0
            //  1  0 |0  0  0 : 1  0| 0  0  0  <<2  ->  0  0  0  1  0
            //  0  1  0 |1  0 : 0  1  0| 1  0  <<3  ->  1  0  0  1  0
            //  0  0  1  0 |1 : 0  0  1  0| 1  <<4  ->  1  0  0  1  0
            // |0  0  1  0  0|: 0  0  1  0  0  <<0  ->  0  0  1  0  0
            //  0 |0  0  0  1 : 0| 0  0  0  1  <<1  ->  0  0  0  1  0
            //  0  0 |1  0  1 ; 0  0| 1  0  1  <<2  ->  1  0  1  0  0
            const uint8_t sym = ((symbols_g[i+offs] << sf | symbols_g[i+offs]) << (i % sf)) >> sf;
            for(int j=0; j<sf; j++) {
                codewords[j+cwi] |= (((sym >> j) & 1) << i);
            }
        }
        cwi += sf;
    }

    return codewords;
}

u8varray_t LoraPhy::dewhiten(const u8varray_t& bytes, const int len) {
    int blen = len==0 ? bytes.size() : len;
    blen = std::min(len, (int)bytes.size());

    u8varray_t bytes_w(blen);
    for(int i=0; i<blen; i++) {
        bytes_w[i] = bytes[i] ^ WHITENING_SEQ[i];
    }

    return bytes_w;
}

chirpval_t LoraPhy::dechirp(const cpvarray_t& sig, const int pos, const bool invert) {
    const cpvarray_t& chp = invert ? this->upchirp : this->downchirp;

    // copy the chirped signal into a 2x buffer for fft processing
    cpvarray_t buf(2*this->ft_bins);
    for(int i=0; i < this->sps; i++) {
        buf[i] = chp[i] * sig[pos + i];
    }

    fft(buf, 2*this->ft_bins);

    // find max fft bin
    int mbin = 0;
    complex_t::value_type mval = 0.0;
    for(int i=0; i < this->ft_bins; i++) {
        const complex_t::value_type val = std::abs(buf[i]) + std::abs(buf[this->ft_bins + i]);
        if(val > mval) {
             mbin = i;
             mval = val;
         }
    }
    //printf("%d) max bin is: %d  (val: %f)\n", pos, mbin, mval);

    return chirpval_t(mbin, static_cast<int>(mval));
}

void LoraPhy::chirp(const int sps, const int fs, const int bw, const bool is_downchirp, cpvarray_t& outbuf) {
    const complex_t::value_type f0 = -bw/2.0 * (is_downchirp ? -1 : 1);
    const complex_t::value_type f1 =  bw/2.0 * (is_downchirp ? -1 : 1);
    const complex_t::value_type t1 =  (sps-1) / static_cast<complex_t::value_type>(fs);

    const complex_t::value_type tau = 2.0*std::acos(-1);
    const complex_t::value_type beta = (f1 - f0) / (2.0*t1);

    cpvarray_t buf(sps);
    for(int i=0; i<sps; ++i) {
        const complex_t::value_type t = i / static_cast<complex_t::value_type>(fs);
        buf[i] = complex_t(
             cos( tau * (beta * t*t + f0*t) ),
            -cos( tau * (beta * t*t + f0*t + 90.0/360.0) )
        );
    }

    outbuf.swap(buf);
}

int LoraPhy::get_sample(cpvarray_t& buf, const int skip) {
    int rc = ifs.rdstate();
    if(rc != 0) {
        if(!ifs.eof()) {
            fprintf(stderr, "error: LoraPhy::get_sample - cannot read from input stream - rc: %d\n", rc);
        }
        return rc;
    }

    this->get_sample_fcn(ifs, buf, skip);

    rc = ifs.rdstate();
    if(rc != 0) {
        if(!ifs.eof()) {
            fprintf(stderr, "error: LoraPhy::get_sample - read from input stream failed - rc: %d\n", rc);
        }
        return rc;
    }

    return 0;
}

int LoraPhy::load(const std::string filename, cpvarray_t& outbuf, const bool swap_iq) {
    FILE* pf = fopen(filename.c_str(), "rb");
    if(pf == NULL) {
        fprintf(stderr, "error: load - failed to open binary file for reading\n");
        return -1;
    }

    int rc = fseek(pf, 0L, SEEK_END);
    if(rc != 0) {
        fprintf(stderr, "error: load - fseek(SEEK_END) failed, rc=%d\n", rc);
        return rc;
    }
    const long filebytes = ftell(pf);

    rc = fseek(pf, 0L, SEEK_SET);
    if(rc != 0) {
        fprintf(stderr, "error: load - fseek(SEEK_SET) failed, rc=%d\n", rc);
        return rc;
    }

    const long count = (filebytes / sizeof(cpvarray_t::value_type));
    cpvarray_t buf(count);
    const size_t ret = fread(&buf[0], sizeof(cpvarray_t::value_type), count, pf);
    fclose(pf);
    if(ret != count) {
        fprintf(stderr, "error: load - failed to read binary file (%zu of %ld read)\n", ret, count);
        return -2;
    }

    if(swap_iq) {
        // swap I and Q channels
        for(cpvarray_t::value_type& cpx : outbuf) {
            cpx = complex_t(cpx.imag(), cpx.real());
        }
    } else {
        outbuf.swap(buf);
    }

    return 0;  // success
}

int LoraPhy::save(const std::string filename, const cpvarray_t& buf) {
    FILE* pf = fopen(filename.c_str(), "wb");
    if(pf == NULL) {
        fprintf(stderr, "error: save - failed to open binary file for writing\n");
        return -1;
    }

    const size_t ret = fwrite(&buf[0], sizeof(cpvarray_t::value_type), buf.size(), pf);
    fclose(pf);
    if(ret != buf.size()) {
        fprintf(stderr, "error: save - failed to write binary file (%zu of %lu written)\n", ret, buf.size());
        return -1;
    }

    return 0;  // success
}

void LoraPhy::print_array(const cpvarray_t& buf, const int num) {
    const int count = (num==0 || num>buf.size()) ? buf.size() : num;
    for(int i=0; i < count; ++i) {
         printf("  %14.10f %+14.10fi  [%11.8f]\n", buf[i].real(), buf[i].imag(), std::abs(buf[i]));
    }
    printf("  size: %zu\n", buf.size());
}

