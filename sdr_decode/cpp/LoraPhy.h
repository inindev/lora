
// copyright (C) 2024, John Clark <inindev@gmail.com>

#include <string>
#include <stdint.h>
#include <fstream>

#include "types.h"
#include "sdr_stream.h"

typedef std::pair<int, int> chirpval_t;
typedef std::pair<int, int> sfdinfo_t;


#define WHITENING_SEQ (uint8_t[]){ \
        0xff, 0xfe, 0xfc, 0xf8, 0xf0, 0xe1, 0xc2, 0x85, 0x0b, 0x17, 0x2f, 0x5e, 0xbc, 0x78, 0xf1, 0xe3, \
        0xc6, 0x8d, 0x1a, 0x34, 0x68, 0xd0, 0xa0, 0x40, 0x80, 0x01, 0x02, 0x04, 0x08, 0x11, 0x23, 0x47, \
        0x8e, 0x1c, 0x38, 0x71, 0xe2, 0xc4, 0x89, 0x12, 0x25, 0x4b, 0x97, 0x2e, 0x5c, 0xb8, 0x70, 0xe0, \
        0xc0, 0x81, 0x03, 0x06, 0x0c, 0x19, 0x32, 0x64, 0xc9, 0x92, 0x24, 0x49, 0x93, 0x26, 0x4d, 0x9b, \
        0x37, 0x6e, 0xdc, 0xb9, 0x72, 0xe4, 0xc8, 0x90, 0x20, 0x41, 0x82, 0x05, 0x0a, 0x15, 0x2b, 0x56, \
        0xad, 0x5b, 0xb6, 0x6d, 0xda, 0xb5, 0x6b, 0xd6, 0xac, 0x59, 0xb2, 0x65, 0xcb, 0x96, 0x2c, 0x58, \
        0xb0, 0x61, 0xc3, 0x87, 0x0f, 0x1f, 0x3e, 0x7d, 0xfb, 0xf6, 0xed, 0xdb, 0xb7, 0x6f, 0xde, 0xbd, \
        0x7a, 0xf5, 0xeb, 0xd7, 0xae, 0x5d, 0xba, 0x74, 0xe8, 0xd1, 0xa2, 0x44, 0x88, 0x10, 0x21, 0x43, \
        0x86, 0x0d, 0x1b, 0x36, 0x6c, 0xd8, 0xb1, 0x63, 0xc7, 0x8f, 0x1e, 0x3c, 0x79, 0xf3, 0xe7, 0xce, \
        0x9c, 0x39, 0x73, 0xe6, 0xcc, 0x98, 0x31, 0x62, 0xc5, 0x8b, 0x16, 0x2d, 0x5a, 0xb4, 0x69, 0xd2, \
        0xa4, 0x48, 0x91, 0x22, 0x45, 0x8a, 0x14, 0x29, 0x52, 0xa5, 0x4a, 0x95, 0x2a, 0x54, 0xa9, 0x53, \
        0xa7, 0x4e, 0x9d, 0x3b, 0x77, 0xee, 0xdd, 0xbb, 0x76, 0xec, 0xd9, 0xb3, 0x67, 0xcf, 0x9e, 0x3d, \
        0x7b, 0xf7, 0xef, 0xdf, 0xbf, 0x7e, 0xfd, 0xfa, 0xf4, 0xe9, 0xd3, 0xa6, 0x4c, 0x99, 0x33, 0x66, \
        0xcd, 0x9a, 0x35, 0x6a, 0xd4, 0xa8, 0x51, 0xa3, 0x46, 0x8c, 0x18, 0x30, 0x60, 0xc1, 0x83, 0x07, \
        0x0e, 0x1d, 0x3a, 0x75, 0xea, 0xd5, 0xaa, 0x55, 0xab, 0x57, 0xaf, 0x5f, 0xbe, 0x7c, 0xf9, 0xf2, \
        0xe5, 0xca, 0x94, 0x28, 0x50, 0xa1, 0x42, 0x84, 0x09, 0x13, 0x27, 0x4f, 0x9f, 0x3f, 0x7f, }

struct LoraPhy
{
    LoraPhy(void);
    virtual ~LoraPhy(void);

    int init(const std::string& filename, const int sample_bits, const bool swap_iq, const int sf, const int bw, const int plen);

    std::tuple<int, int, int> detect_preamble(const bool invert=false);
    sfdinfo_t detect_sfd(const bool invert=false);

    std::tuple<uint8_t, uint8_t, uint8_t, bool> decode_header(const bool invert=false);
    u16varray_t decode_payload(const uint8_t payload_len, const bool invert=false);

    uint8_t calc_payload_symbol_count(const uint8_t payload_len);
    uint8_t calc_header_csum(const u8varray_t& header);
    u16varray_t gray_decode(const u16varray_t& symbols, const bool use_ldro=false);

    u8varray_t hamming_decode(const u16varray_t& codewords, const int cr);
    u16varray_t diag_deinterleave(const u16varray_t& symbols_g, const int cr, const int sf);
    u8varray_t dewhiten(const u8varray_t& bytes, const int len=0);

    chirpval_t dechirp(const cpvarray_t& sig, const int pos=0, const bool invert=false);
    static void chirp(const int sps, const int fs, const int bw, const bool is_downchirp, cpvarray_t& outbuf);

    int get_sample(cpvarray_t& buf, const int skip=0);

    static int load(const std::string filename, cpvarray_t& buf, const bool swap_iq=false);
    static int save(const std::string filename, const cpvarray_t& buf);
    static void print_array(const cpvarray_t& buf, const int num=0);


private:
    int sf;           // spreading factor (7,8,9,10,11,12)
    int sr;           // sample rate: 2x
    int bw;           // bandwidth (125e3, 250e3, 500e3)
    int sps;          // samples per symbol
    int fs;           // sampling frequency (2x bandwidth)
    int cr;           // coding rate: (1:4/5 2:4/6 3:4/7 4:4/8)
    int use_ldro;     // low data rate optimization
    int use_crc;      // calculate cyclic redundancy check
    int use_hamming;  // hamming error detection and correction
    int has_header;   // data has header
    int plen;         // chrips in preamble
    int ft_ratio;     // fft bins per sps
    int ft_sym_fct;   // fft symbol factor
    int ft_bins;      // fft bins per symbol

    cpvarray_t upchirp;
    cpvarray_t downchirp;

    std::ifstream ifs;
    std::function<void(std::ifstream& ifs, cpvarray_t&, int)> get_sample_fcn;
};

