

def genovec_3freq():
    # src/plink_common.cpp:8024
    pass

def get_final_mask(sample_ct, BITCT2):
    # lib/plink_common.hpp:2605

    uii = sample_ct % BITCT2
    if uii != 0:
        return (1 << (2 * uii)) - 1
    else:
        return -1