
// Copyright (c) 2024, John Clark <inindev@gmail.com>

#pragma once

#include <cstring>
#include <iostream>
#include <fstream>
#include <memory>
#include <functional>

#include "types.h"


namespace sdr_stream
{
    static constexpr complex_t::value_type cu82f[256] = {
        (-127.4/128), (-126.4/128), (-125.4/128), (-124.4/128), (-123.4/128), (-122.4/128), (-121.4/128), (-120.4/128), (-119.4/128), (-118.4/128), (-117.4/128), (-116.4/128), (-115.4/128), (-114.4/128), (-113.4/128), (-112.4/128),
        (-111.4/128), (-110.4/128), (-109.4/128), (-108.4/128), (-107.4/128), (-106.4/128), (-105.4/128), (-104.4/128), (-103.4/128), (-102.4/128), (-101.4/128), (-100.4/128), ( -99.4/128), ( -98.4/128), ( -97.4/128), ( -96.4/128),
        ( -95.4/128), ( -94.4/128), ( -93.4/128), ( -92.4/128), ( -91.4/128), ( -90.4/128), ( -89.4/128), ( -88.4/128), ( -87.4/128), ( -86.4/128), ( -85.4/128), ( -84.4/128), ( -83.4/128), ( -82.4/128), ( -81.4/128), ( -80.4/128),
        ( -79.4/128), ( -78.4/128), ( -77.4/128), ( -76.4/128), ( -75.4/128), ( -74.4/128), ( -73.4/128), ( -72.4/128), ( -71.4/128), ( -70.4/128), ( -69.4/128), ( -68.4/128), ( -67.4/128), ( -66.4/128), ( -65.4/128), ( -64.4/128),
        ( -63.4/128), ( -62.4/128), ( -61.4/128), ( -60.4/128), ( -59.4/128), ( -58.4/128), ( -57.4/128), ( -56.4/128), ( -55.4/128), ( -54.4/128), ( -53.4/128), ( -52.4/128), ( -51.4/128), ( -50.4/128), ( -49.4/128), ( -48.4/128),
        ( -47.4/128), ( -46.4/128), ( -45.4/128), ( -44.4/128), ( -43.4/128), ( -42.4/128), ( -41.4/128), ( -40.4/128), ( -39.4/128), ( -38.4/128), ( -37.4/128), ( -36.4/128), ( -35.4/128), ( -34.4/128), ( -33.4/128), ( -32.4/128),
        ( -31.4/128), ( -30.4/128), ( -29.4/128), ( -28.4/128), ( -27.4/128), ( -26.4/128), ( -25.4/128), ( -24.4/128), ( -23.4/128), ( -22.4/128), ( -21.4/128), ( -20.4/128), ( -19.4/128), ( -18.4/128), ( -17.4/128), ( -16.4/128),
        ( -15.4/128), ( -14.4/128), ( -13.4/128), ( -12.4/128), ( -11.4/128), ( -10.4/128), (  -9.4/128), (  -8.4/128), (  -7.4/128), (  -6.4/128), (  -5.4/128), (  -4.4/128), (  -3.4/128), (  -2.4/128), (  -1.4/128), (  -0.4/128),
        (   0.6/128), (   1.6/128), (   2.6/128), (   3.6/128), (   4.6/128), (   5.6/128), (   6.6/128), (   7.6/128), (   8.6/128), (   9.6/128), (  10.6/128), (  11.6/128), (  12.6/128), (  13.6/128), (  14.6/128), (  15.6/128),
        (  16.6/128), (  17.6/128), (  18.6/128), (  19.6/128), (  20.6/128), (  21.6/128), (  22.6/128), (  23.6/128), (  24.6/128), (  25.6/128), (  26.6/128), (  27.6/128), (  28.6/128), (  29.6/128), (  30.6/128), (  31.6/128),
        (  32.6/128), (  33.6/128), (  34.6/128), (  35.6/128), (  36.6/128), (  37.6/128), (  38.6/128), (  39.6/128), (  40.6/128), (  41.6/128), (  42.6/128), (  43.6/128), (  44.6/128), (  45.6/128), (  46.6/128), (  47.6/128),
        (  48.6/128), (  49.6/128), (  50.6/128), (  51.6/128), (  52.6/128), (  53.6/128), (  54.6/128), (  55.6/128), (  56.6/128), (  57.6/128), (  58.6/128), (  59.6/128), (  60.6/128), (  61.6/128), (  62.6/128), (  63.6/128),
        (  64.6/128), (  65.6/128), (  66.6/128), (  67.6/128), (  68.6/128), (  69.6/128), (  70.6/128), (  71.6/128), (  72.6/128), (  73.6/128), (  74.6/128), (  75.6/128), (  76.6/128), (  77.6/128), (  78.6/128), (  79.6/128),
        (  80.6/128), (  81.6/128), (  82.6/128), (  83.6/128), (  84.6/128), (  85.6/128), (  86.6/128), (  87.6/128), (  88.6/128), (  89.6/128), (  90.6/128), (  91.6/128), (  92.6/128), (  93.6/128), (  94.6/128), (  95.6/128),
        (  96.6/128), (  97.6/128), (  98.6/128), (  99.6/128), ( 100.6/128), ( 101.6/128), ( 102.6/128), ( 103.6/128), ( 104.6/128), ( 105.6/128), ( 106.6/128), ( 107.6/128), ( 108.6/128), ( 109.6/128), ( 110.6/128), ( 111.6/128),
        ( 112.6/128), ( 113.6/128), ( 114.6/128), ( 115.6/128), ( 116.6/128), ( 117.6/128), ( 118.6/128), ( 119.6/128), ( 120.6/128), ( 121.6/128), ( 122.6/128), ( 123.6/128), ( 124.6/128), ( 125.6/128), ( 126.6/128), ( 127.6/128),
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

    inline void cs16_iq_sample(std::ifstream& ifs, cpvarray_t& buf, const int skip) {
        int offs = 0;
        if(skip < 0) {  // partial buffer read
            offs = ((0 - skip) % buf.size());
            std::memcpy(&buf[0], &buf[buf.size() - offs], (offs * sizeof(cpvarray_t::value_type)));
        }

        else if(skip > 0) {  // move the stream forward by the specified amount
            ifs.ignore(skip * sizeof(uint32_t));
        }

        uint32_t u32a[buf.size() - offs];
        ifs.read((char*)u32a, sizeof(u32a));
        for(int i=0; i<(sizeof(u32a)/sizeof(uint32_t)); i++) {
            buf[i+offs] = complex_t(
                (int16_t(u32a[i] & 0xffff) / 32768.0), // real
                (int16_t(u32a[i] >> 16) / 32768.0)     // imag
            );
        }
    }

    inline void cs16_qi_sample(std::ifstream& ifs, cpvarray_t& buf, const int skip) {
        int offs = 0;
        if(skip < 0) {  // partial buffer read
            offs = ((0 - skip) % buf.size());
            std::memcpy(&buf[0], &buf[buf.size() - offs], (offs * sizeof(cpvarray_t::value_type)));
        }

        else if(skip > 0) {  // move the stream forward by the specified amount
            ifs.ignore(skip * sizeof(uint32_t));
        }

        uint32_t u32a[buf.size() - offs];
        ifs.read((char*)u32a, sizeof(u32a));
        for(int i=0; i<(sizeof(u32a)/sizeof(uint32_t)); i++) {
            buf[i+offs] = complex_t(
                (int16_t(u32a[i] >> 16) / 32768.0),   // imag
                (int16_t(u32a[i] & 0xffff) / 32768.0) // real
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

