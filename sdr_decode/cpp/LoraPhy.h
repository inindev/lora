
// copyright (C) 2024, John Clark <inindev@gmail.com>

#include <string>
#include <stdint.h>

#include "types.h"

typedef std::pair<int, int> chirpval_t;
typedef std::pair<long, long> sfdinfo_t;


struct LoraPhy
{
    LoraPhy(void);
    virtual ~LoraPhy(void);

    void init(const int sf, const int bw, const int plen);

    std::tuple<long, int, int> detect_preamble(const cpvarray_t& sig, const long pos=0, const bool invert=false);
    sfdinfo_t detect_sfd(const cpvarray_t& sig, const long pos=0, const bool invert=false);
    bool decode_header(const cpvarray_t& sig, const long pos, const bool invert=false);
    void diag_deinterleave(const uint32_t symbols_g[8], const uint32_t bits, uint32_t* codewords);

    chirpval_t dechirp(const cpvarray_t& sig, const long pos=0, const bool invert=false);
    static void chirp(const int sps, const int fs, const int bw, const bool is_downchirp, cpvarray_t& outbuf);

    static int load(const std::string filename, cpvarray_t& buf, const bool swap_iq=false);
    static int save(const std::string filename, const cpvarray_t& buf);
    static void print_array(const cpvarray_t& buf, const long num=0);


private:
    int sf;    // spreading factor (7,8,9,10,11,12)
    int bw;    // bandwidth (125e3, 250e3, 500e3)
    int sps;   // samples per symbol
    int fs;    // sampling frequency (2x bandwidth)
    int plen;  // chrips in preamble
    int ft_det_bins;  // fft detect bin variance
    cpvarray_t upchirp;
    cpvarray_t downchirp;
};
