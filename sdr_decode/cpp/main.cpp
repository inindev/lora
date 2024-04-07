
#include "LoraPhy.h"
#include "temp.h"


int main(void) {

    // test iq / chirp inversion
    const bool test_inv = true;
    const bool swap_iq = !test_inv ? false: true;
    const bool invert  = !test_inv ? false: !swap_iq;

    LoraPhy phy;
    phy.init(7, 125e3, 8);

    cpvarray_t sig;
    int rc = LoraPhy::load("./lora.raw", sig, swap_iq);
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
    printf("detect_preamble - pos: %ld  netid1: %d  netid2: %d\n", ppos, netid1, netid2);

    sfdinfo_t sfd_info = phy.detect_sfd(sig, ppos, invert);
    if(sfd_info.first < 0) {
        printf("sfd not detected\n");
        return -2;
    }
    printf("sfd pos: %ld  hdr pos: %ld\n", sfd_info.first, sfd_info.second);

    return 0;
}

