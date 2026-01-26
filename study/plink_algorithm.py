

def genovec_3freq(sample_ctl2):
    # src/plink_common.cpp:8024
    acc_even = 0
    acc_odd = 0
    acc_and = 0

    cur_decr = 120
    sample_ctl2 -= sample_ctl2 % 12

    return sample_ctl2
    pass

def get_final_mask(sample_ct, BITCT2):
    # lib/plink_common.hpp:2605

    uii = sample_ct % BITCT2
    if uii != 0:
        return (1 << (2 * uii)) - 1
    else:
        return -1