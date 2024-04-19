
// copyright (C) 2024, John Clark <inindev@gmail.com>

#include "LoraPhy.h"
#include "fft.h"


LoraPhy::LoraPhy(void)
  : sf(0), bw(0), sps(0), fs(0), cr(0),
    use_ldro(0), use_crc(0), use_hamming(0), has_header(0),
    plen(0), ft_det_bins(0) { }

LoraPhy::~LoraPhy() { }

void LoraPhy::init(const int sf, const int bw, const int plen) {
    this->sf          = sf;
    this->bw          = bw;
    this->sps         = 2 * pow(2, sf);
    this->fs          = 2 * bw;
    this->cr          = 1; // 4/5
    this->use_ldro    = 0;
    this->use_crc     = 1;
    this->use_hamming = 1;
    this->has_header  = 1;
    this->plen        = plen; // preamble detect symbols
    this->ft_det_bins = 4;    // fft detect bin variance

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

// <payload_len, cr, crc, is_valid>
std::tuple<uint8_t, uint8_t, uint8_t, bool> LoraPhy::decode_header(const cpvarray_t& sig, const long pos, const bool invert) {
    u16varray_t symbols(8);
    for(int ii=0; ii<8; ii++) {
        int fbin = this->dechirp(sig, pos+ii*this->sps, invert).first;
        //printf("%4d) decode_header - fbin: %d\n", ii, fbin);

        if(invert) {
            fbin = this->sps - fbin;  // inverted chirp, non-inverted IQ
        }
        // non-inverted chirp, inverted IQ needs no adjustment

        symbols[ii] = fbin / 2;
        //printf("symbols[%d]: %d\n", ii, symbols[ii]);
    }

    // gray decoding
    const u16varray_t& symbols_g = gray_decode(symbols, true);

    // deinterleave
    const u16varray_t& codewords = diag_deinterleave(symbols_g, 4, this->sf-2); // cr, sf
    // for(int ii=0; ii<this->sf-2; ii++) {
    //    printf("%4d) cw: %d\n", ii, codewords[ii]);
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

u16varray_t LoraPhy::decode_payload(const cpvarray_t& sig, const long pos_hdr, const uint8_t payload_len, const bool invert) {
    const uint8_t symcnt = calc_payload_symbol_count(payload_len);

    u16varray_t symbols(symcnt);
    long pos = pos_hdr + 8*this->sps;
    for(uint8_t ii=0; ii<symcnt; ii++) {
        if(pos > (sig.size() - this->sps)) {
            printf("error: unexpected end of data reached - pos:%ld  sym_num:%d\n", pos, ii);
            return u16varray_t();
        }

        int fbin = this->dechirp(sig, pos, invert).first;
        if(invert) {
            fbin = this->sps - fbin + 1;  // inverted chirp, non-inverted IQ
        }
        // non-inverted chirp, inverted IQ needs no adjustment

        symbols[ii] = fbin / 2;
        //printf("%d) pos:%4ld  fbin:%3d  -->  sym:%3d\n", ii, pos, fbin, symbols[ii]);

        pos = pos + this->sps;
    }

    return symbols;
}

uint8_t LoraPhy::calc_payload_symbol_count(const uint8_t payload_len) {
    // see SX1276_Datasheet.pdf p.31
    const int symcnt = std::ceil((2*double(payload_len)-this->sf +7+ 4*this->use_crc-5*this->has_header) / (this->sf-2*this->use_ldro)) * (this->cr+4);
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
                printf("parity data correction - data elem:%d  bit: %d  repaired data: %d\n", i, be, result[i]);
            } else {
                printf(" !!!! unrecoverable parity error !!!! - data elem:%d  bit: %d\n", i, be);
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
        for(int ii=0; ii<cr_bits; ii++) {                 // 154 0 163 92 0
            // |0  0  1  0  0|: 0  0  1  0  0  <<0  ->  0  0  1  0  0
            //  0 |1  0  1  0 : 0| 1  0  1  0  <<1  ->  1  0  1  0  0
            //  1  0 |0  0  0 : 1  0| 0  0  0  <<2  ->  0  0  0  1  0
            //  0  1  0 |1  0 : 0  1  0| 1  0  <<3  ->  1  0  0  1  0
            //  0  0  1  0 |1 : 0  0  1  0| 1  <<4  ->  1  0  0  1  0
            // |0  0  1  0  0|: 0  0  1  0  0  <<0  ->  0  0  1  0  0
            //  0 |0  0  0  1 : 0| 0  0  0  1  <<1  ->  0  0  0  1  0
            //  0  0 |1  0  1 ; 0  0| 1  0  1  <<2  ->  1  0  1  0  0
            const uint8_t sym = ((symbols_g[ii+offs] << sf | symbols_g[ii+offs]) << (ii % sf)) >> sf;
            for(int jj=0; jj<sf; jj++) {
                codewords[jj+cwi] |= (((sym >> jj) & 1) << ii);
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

int LoraPhy::load(const std::string filename, cpvarray_t& outbuf, const bool swap_iq) {
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
