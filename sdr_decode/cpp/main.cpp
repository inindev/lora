
#include "LoraPhy.h"
#include "temp.h"


int main(void) {

    // test iq / chirp inversion
    const bool test_inv = false;

    LoraPhy phy;
    phy.init(7, 125e3, 8);

    cpvarray_t sig;
    int rc = load_float("./lora.raw", sig, test_inv);
    if(rc != 0)
    {
        return rc;
    }
    // LoraPhy::print_array(sig, 16);

    long pos = phy.detect_preamble(sig, 0, !test_inv);
    printf("detect_preamble - pos: %ld\n", pos);

    sfdinfo_t sfd_info = phy.detect_sfd(sig, pos, !test_inv);
    printf("sfd pos: %ld  hdr pos: %ld\n", sfd_info.first, sfd_info.second);

    return 0;
}

