
#include "LoraPhy.h"


// rx_sdr -g12 -f 910300000 -s 250000 -F CF32 /tmp/lora.raw
int main(void) {
    // test iq / chirp inversion
    const bool test_inv = false;
    const bool swap_iq = !test_inv ? false: true;
    const bool invert  = !test_inv ? false: !swap_iq;

    LoraPhy phy;
    phy.init(7, 125e3, 8);

    cpvarray_t sig;
    int rc = LoraPhy::load("/tmp/lora.raw", sig, swap_iq);
    if(rc != 0)
    {
        return rc;
    }
    // LoraPhy::print_array(sig, 16);

    const auto [ppos, netid1, netid2] = phy.detect_preamble(sig, 0, invert);
    if(ppos < 0) {
        printf("preamble not detected\n");
        return -1;
    }
    printf("netid1: %d  netid2: %d\n", netid1, netid2);

    const auto [pos_sfd, pos_hdr] = phy.detect_sfd(sig, ppos, invert);
    if(pos_sfd < 0) {
        printf("sfd not detected\n");
        return -2;
    }
    //printf("sfd pos: %ld  hdr pos: %ld\n", pos_sfd, pos_hdr);

    const auto [payload_len, cr, crc, is_valid] = phy.decode_header(sig, pos_hdr, invert);
    if(!is_valid) {
        printf("header is invalid\n");
        return -3;
    }
    printf("header is valid - payload_len: %d  cr: %d  crc: %d\n", payload_len, cr, crc);


    const u16varray_t& payload_symbols = phy.decode_payload(sig, pos_hdr, payload_len, invert);
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
        printf("%02x ", payload[i]);
    }
    printf("\n");
    printf("csum1: %d\n", bytes[payload_len]);
    printf("csum2: %d\n", bytes[payload_len+1]);

    return 0;
}
