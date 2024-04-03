
// copyright (C) 2024, John Clark <inindev@gmail.com>

#include <string>

#include "types.h"


struct LoraPhy
{
    LoraPhy(void);
    virtual ~LoraPhy(void);

    void init(const int sf, const int bw, const int plen);

    long detect_preamble(const cpvarray_t& sig, const long pos=0, const bool invert=false);
    long dechirp(const cpvarray_t& sig, const long pos=0, const bool invert=false);

    static void chirp(const int sps, const int fs, const int bw, const bool is_downchirp, cpvarray_t& outbuf);

    static int load(const std::string filename, cpvarray_t& buf);
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
