
#include "LoraPhy.h"


int do_work(LoraPhy& phy, std::ifstream& ifs, const bool swap_iq, const bool invert) {
    const auto [ppos, netid1, netid2] = phy.detect_preamble(ifs, invert);
    if(ppos < 0) {
        fprintf(stderr, "error: preamble not detected\n");
        return 3;
    }
    printf("netid1: %d  netid2: %d\n", netid1, netid2);

    const auto [pos_sfd, pos_hdr] = phy.detect_sfd(ifs, invert);
    if(pos_sfd < 0) {
        fprintf(stderr, "error: sfd not detected\n");
        return 4;
    }
    // printf("sfd pos: %ld  hdr pos: %ld\n", pos_sfd, pos_hdr);

    // cpvarray_t sig(256 * 8);
    // phy.get_sample(ifs, sig);
    const auto [payload_len, cr, crc, is_valid] = phy.decode_header(ifs, invert);
    if(!is_valid) {
        fprintf(stderr, "error: header is invalid\n");
        return 5;
    }
    printf("header is valid - payload_len: %d  cr: %d  crc: %d\n", payload_len, cr, crc);

    const u16varray_t& payload_symbols = phy.decode_payload(ifs, payload_len, invert);
    // for(int i=0; i<payload_symbols.size(); i++) {
    //     printf("%d) payload_symbols:%3d\n", i, payload_symbols[i]);
    // }

    const u16varray_t& payload_symbols_g = phy.gray_decode(payload_symbols, false);
    // for(int i=0; i<payload_symbols_g.size(); i++) {
    //     printf("%d) payload_symbols_g:%3d\n", i, payload_symbols_g[i]);
    // }

    const u16varray_t& codewords = phy.diag_deinterleave(payload_symbols_g, 1, 7);  // cr, sf
    // for(int i=0; i<codewords.size(); i++) {
    //     printf("%d) codewords: %d\n", i, codewords[i]);
    // }

    // hamming decode
    const u8varray_t& nibbles = phy.hamming_decode(codewords, 1); // cr
    // for(int i=0; i<nibbles.size(); i++) {
    //     printf("%d) nibbles: %d\n", i, nibbles[i]);
    // }

    const int blen = nibbles.size() / 2;
    u8varray_t bytes(blen);
    for(int i=0; i<blen; i++) {
        const int n = i * 2;
        bytes[i] = (nibbles[n+1] << 4) | nibbles[n];
    }

    const u8varray_t& payload = phy.dewhiten(bytes, payload_len);
    printf("payload_len: %zu\n", payload.size());
    printf("payload: ");
    for(int i=0; i<payload.size(); i++) {
        printf("%02x", payload[i]);
    }
    printf("\n");
    printf("csum1: %d\n", bytes[payload_len]);
    printf("csum2: %d\n", bytes[payload_len+1]);

    return 0;
}

// https://github.com/rxseger/rx_tools
// rx_sdr -g12 -f 910300000 -s 250000 -F CF32 /tmp/lora.raw
// rx_sdr -g12 -f 910300000 -s 250000 -F CF32 - | ./loraphy -
int main(int argc, char** argv) {
    if(argc < 2) {
        fprintf(stderr, "error: no input source specified\n");
        return 1;
    }

    // test iq / chirp inversion
    const bool test_inv = false;
    const bool swap_iq = !test_inv ? false: true;
    const bool invert  = !test_inv ? false: !swap_iq;

    const std::string filename(argv[1]);
    std::ifstream ifs(("-"==filename) ? "/dev/stdin" : filename);
    if(!ifs) {
        fprintf(stderr, "error: unable to open input source - filename: %s\n", filename.c_str());
        return 2;
    }

    LoraPhy phy;
    phy.init(7, 125e3, 8);

    for(;;) {
        int rc = do_work(phy, ifs, swap_iq, invert);
        if(rc != 0) {
            return rc;
        }
    }

    return 0;
}

