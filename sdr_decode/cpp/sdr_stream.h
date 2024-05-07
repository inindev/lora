
#include <iostream>
#include <fstream>
#include <memory>
#include <functional>

#include "types.h"


namespace sdr_stream
{
    static constexpr complex_t::value_type cu82f[256] = {
        (-127/128.0f), (-126/128.0f), (-125/128.0f), (-124/128.0f), (-123/128.0f), (-122/128.0f), (-121/128.0f), (-120/128.0f), (-119/128.0f), (-118/128.0f), (-117/128.0f), (-116/128.0f), (-115/128.0f), (-114/128.0f), (-113/128.0f), (-112/128.0f),
        (-111/128.0f), (-110/128.0f), (-109/128.0f), (-108/128.0f), (-107/128.0f), (-106/128.0f), (-105/128.0f), (-104/128.0f), (-103/128.0f), (-102/128.0f), (-101/128.0f), (-100/128.0f), ( -99/128.0f), ( -98/128.0f), ( -97/128.0f), ( -96/128.0f),
        ( -95/128.0f), ( -94/128.0f), ( -93/128.0f), ( -92/128.0f), ( -91/128.0f), ( -90/128.0f), ( -89/128.0f), ( -88/128.0f), ( -87/128.0f), ( -86/128.0f), ( -85/128.0f), ( -84/128.0f), ( -83/128.0f), ( -82/128.0f), ( -81/128.0f), ( -80/128.0f),
        ( -79/128.0f), ( -78/128.0f), ( -77/128.0f), ( -76/128.0f), ( -75/128.0f), ( -74/128.0f), ( -73/128.0f), ( -72/128.0f), ( -71/128.0f), ( -70/128.0f), ( -69/128.0f), ( -68/128.0f), ( -67/128.0f), ( -66/128.0f), ( -65/128.0f), ( -64/128.0f),
        ( -63/128.0f), ( -62/128.0f), ( -61/128.0f), ( -60/128.0f), ( -59/128.0f), ( -58/128.0f), ( -57/128.0f), ( -56/128.0f), ( -55/128.0f), ( -54/128.0f), ( -53/128.0f), ( -52/128.0f), ( -51/128.0f), ( -50/128.0f), ( -49/128.0f), ( -48/128.0f),
        ( -47/128.0f), ( -46/128.0f), ( -45/128.0f), ( -44/128.0f), ( -43/128.0f), ( -42/128.0f), ( -41/128.0f), ( -40/128.0f), ( -39/128.0f), ( -38/128.0f), ( -37/128.0f), ( -36/128.0f), ( -35/128.0f), ( -34/128.0f), ( -33/128.0f), ( -32/128.0f),
        ( -31/128.0f), ( -30/128.0f), ( -29/128.0f), ( -28/128.0f), ( -27/128.0f), ( -26/128.0f), ( -25/128.0f), ( -24/128.0f), ( -23/128.0f), ( -22/128.0f), ( -21/128.0f), ( -20/128.0f), ( -19/128.0f), ( -18/128.0f), ( -17/128.0f), ( -16/128.0f),
        ( -15/128.0f), ( -14/128.0f), ( -13/128.0f), ( -12/128.0f), ( -11/128.0f), ( -10/128.0f), (  -9/128.0f), (  -8/128.0f), (  -7/128.0f), (  -6/128.0f), (  -5/128.0f), (  -4/128.0f), (  -3/128.0f), (  -2/128.0f), (  -1/128.0f), (   0/128.0f),
        (   1/128.0f), (   2/128.0f), (   3/128.0f), (   4/128.0f), (   5/128.0f), (   6/128.0f), (   7/128.0f), (   8/128.0f), (   9/128.0f), (  10/128.0f), (  11/128.0f), (  12/128.0f), (  13/128.0f), (  14/128.0f), (  15/128.0f), (  16/128.0f),
        (  17/128.0f), (  18/128.0f), (  19/128.0f), (  20/128.0f), (  21/128.0f), (  22/128.0f), (  23/128.0f), (  24/128.0f), (  25/128.0f), (  26/128.0f), (  27/128.0f), (  28/128.0f), (  29/128.0f), (  30/128.0f), (  31/128.0f), (  32/128.0f),
        (  33/128.0f), (  34/128.0f), (  35/128.0f), (  36/128.0f), (  37/128.0f), (  38/128.0f), (  39/128.0f), (  40/128.0f), (  41/128.0f), (  42/128.0f), (  43/128.0f), (  44/128.0f), (  45/128.0f), (  46/128.0f), (  47/128.0f), (  48/128.0f),
        (  49/128.0f), (  50/128.0f), (  51/128.0f), (  52/128.0f), (  53/128.0f), (  54/128.0f), (  55/128.0f), (  56/128.0f), (  57/128.0f), (  58/128.0f), (  59/128.0f), (  60/128.0f), (  61/128.0f), (  62/128.0f), (  63/128.0f), (  64/128.0f),
        (  65/128.0f), (  66/128.0f), (  67/128.0f), (  68/128.0f), (  69/128.0f), (  70/128.0f), (  71/128.0f), (  72/128.0f), (  73/128.0f), (  74/128.0f), (  75/128.0f), (  76/128.0f), (  77/128.0f), (  78/128.0f), (  79/128.0f), (  80/128.0f),
        (  81/128.0f), (  82/128.0f), (  83/128.0f), (  84/128.0f), (  85/128.0f), (  86/128.0f), (  87/128.0f), (  88/128.0f), (  89/128.0f), (  90/128.0f), (  91/128.0f), (  92/128.0f), (  93/128.0f), (  94/128.0f), (  95/128.0f), (  96/128.0f),
        (  97/128.0f), (  98/128.0f), (  99/128.0f), ( 100/128.0f), ( 101/128.0f), ( 102/128.0f), ( 103/128.0f), ( 104/128.0f), ( 105/128.0f), ( 106/128.0f), ( 107/128.0f), ( 108/128.0f), ( 109/128.0f), ( 110/128.0f), ( 111/128.0f), ( 112/128.0f),
        ( 113/128.0f), ( 114/128.0f), ( 115/128.0f), ( 116/128.0f), ( 117/128.0f), ( 118/128.0f), ( 119/128.0f), ( 120/128.0f), ( 121/128.0f), ( 122/128.0f), ( 123/128.0f), ( 124/128.0f), ( 125/128.0f), ( 126/128.0f), ( 127/128.0f), ( 128/128.0f),
    };

    inline void cu8_iq_sample(std::ifstream& ifs, cpvarray_t& buf, const int skip) {
        int offs = 0;
        if(skip < 0) {  // partial buffer read
            offs = ((0 - skip) % buf.size());
            std::memcpy(&buf[0], &buf[buf.size() - offs], (offs * sizeof(cpvarray_t::value_type)));
        }

        else if(skip > 0) {  // move the stream forward by the specified amount
            ifs.ignore(skip * sizeof(uint16_t));
        }

        uint16_t u16a[buf.size() - offs];
        ifs.read((char*)u16a, sizeof(u16a));
        for(int i=0; i<(sizeof(u16a)/sizeof(uint16_t)); i++) {
            buf[i+offs] = complex_t(
                cu82f[u16a[i] & 0xff], // real
                cu82f[u16a[i] >> 8]    // imag
            );
        }
    }

    // i & q swapped
    inline void cu8_qi_sample(std::ifstream& ifs, cpvarray_t& buf, const int skip) {
        int offs = 0;
        if(skip < 0) {  // partial buffer read
            offs = ((0 - skip) % buf.size());
            std::memcpy(&buf[0], &buf[buf.size() - offs], (offs * sizeof(cpvarray_t::value_type)));
        }

        else if(skip > 0) {  // move the stream forward by the specified amount
            ifs.ignore(skip * sizeof(uint16_t));
        }

        uint16_t u16a[buf.size() - offs];
        ifs.read((char*)u16a, sizeof(u16a));
        for(int i=0; i<(sizeof(u16a)/sizeof(uint16_t)); i++) {
            buf[i+offs] = complex_t(
                cu82f[u16a[i] >> 8],  // imag
                cu82f[u16a[i] & 0xff] // real
            );
        }
    }

    inline void cf32_iq_sample(std::ifstream& ifs, cpvarray_t& buf, const int skip) {
        int offs = 0;
        if(skip < 0) {  // partial buffer read
            offs = ((0 - skip) % buf.size());
            std::memcpy(&buf[0], &buf[buf.size() - offs], (offs * sizeof(cpvarray_t::value_type)));
        }

        else if(skip > 0) {  // move the stream forward by the specified amount
            ifs.ignore(skip * sizeof(cpvarray_t::value_type));
        }

        ifs.read(reinterpret_cast<char*>(&buf[offs]), (buf.size() - offs) * sizeof(cpvarray_t::value_type));
    }

    // i & q swapped
    inline void cf32_qi_sample(std::ifstream& ifs, cpvarray_t& buf, const int skip) {
        int offs = 0;
        if(skip < 0) {  // partial buffer read
            offs = ((0 - skip) % buf.size());
            std::memcpy(&buf[0], &buf[buf.size() - offs], (offs * sizeof(cpvarray_t::value_type)));
        }

        else if(skip > 0) {  // move the stream forward by the specified amount
            ifs.ignore(skip * sizeof(cpvarray_t::value_type));
        }

        ifs.read(reinterpret_cast<char*>(&buf[offs]), (buf.size() - offs) * sizeof(cpvarray_t::value_type));

        // swap I and Q channels
        for(int i=offs; i<buf.size(); i++) {
            buf[i] = complex_t(buf[i].imag(), buf[i].real());
        }
    }
}

