
#include "LoraPhy.h"


int main(void) {
    LoraPhy phy;
    phy.init(7, 125e3, 8);

    cpvarray_t sig;
    int rc = LoraPhy::load("./lora.raw", sig);
    if(rc != 0)
    {
        return rc;
    }
    // LoraPhy::print_array(sig, 16);

    long pos = phy.detect_preamble(sig, 0, true);
    printf("detect_preamble - pos: %ld\n", pos);

    return 0;
}

