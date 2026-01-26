from pandas_plink import read_plink
import pandas as pd
from data_structs import *
import plink_algorithm
from functools import cached_property

BITCT = 64
BITCT2 = BITCT / 2
VEC_BYTES = 16
VEC_BITS = VEC_BYTES * 8
VEC_WORDS = VEC_BITS / BITCT

BITCT_TO_VECCT = lambda val: (((val) + (VEC_BITS - 1)) / VEC_BITS)
BITCT_TO_ALIGNED_WORDCT = lambda val: VEC_WORDS * BITCT_TO_VECCT(val)
BITCT_TO_WORDCT = lambda val: (((val) + (BITCT - 1)) / BITCT)
QUATERCT_TO_VECCT = lambda val: (((val) + ((VEC_BITS / 2) - 1)) / (VEC_BITS / 2))

QUATERCT_TO_WORDCT = lambda val: (((val) + (BITCT2 - 1)) / BITCT2)
QUATERCT_TO_ALIGNED_WORDCT = lambda val: (VEC_WORDS * QUATERCT_TO_VECCT(val))


class Genotype:

    def __init__(self, plink_files, base_file, cov_file):
        self.bim, self.fam, self.G = read_plink(plink_files, verbose=False)
        self.base_data = pd.read_csv(base_file, sep="\t", compression="gzip")
        self.covariate_data = pd.read_csv(cov_file, sep=" ")
        self.genotype_filename = f"{plink_files}.bed"

        self.num_samples = len(self.fam) 

        # inc/storage.hpp:193
        self.clumping_info = Clumping()

    @property
    def m_sort_by_p_index(self) -> list[int]:
        """
        Sort by: chr ↑ → p_value ↑ → loc ↑ → rs ↑
        """
        # inc/genotype.hpp:139
        idx = self.bim.merge(
                    self.base_data[["SNP", "P", "CHR", "BP"]], 
                    left_on="snp", 
                    right_on="SNP", 
                    how="left")\
                .sort_values(by=["CHR", "P", "BP", "SNP"])\
                .index.to_list()                # src/snp.cpp:27
        
        return idx

    def get_chrom_boundary(self) -> list[SnpRange]:
        # src/genotype.cpp:1171
        dta = self.bim.loc[:,"chrom"].value_counts()\
                        .reset_index()\
                        .astype({"chrom": int})\
                        .sort_values(by="chrom")
        dta["second"] = dta["count"].cumsum()
        dta["first"] = dta["second"].shift(1, fill_value=0)
        dta_list = [SnpRange(start, end) for start, end in dta.loc[:,["first", "second"]].values.tolist()]

        return dta_list

    @cached_property
    def m_max_window_size(self):
        # src/genotype.cpp:99 & src/genotype.cpp:111
        # src/genotype.cpp:83
        m_max_window_size = 0
        prev_chr = 0
        low_bound = 0
        prev_loc = 0
        cur_dist = 0

        for i_dx, (chr, pos) in enumerate(self.bim[["chrom", "pos"]].values):
            ## Get max gap between chromosomes (number of snp)
            if (prev_chr != chr):
                prev_chr = chr
                prev_loc = pos
                low_bound = i_dx
            snp_distance = i_dx - low_bound
            if (m_max_window_size < snp_distance): 
                m_max_window_size = snp_distance
            
            ## Lift up low_bound by limiting clumping distance (position distance)
            # src/genotype.cpp:101
            cur_dist = pos - prev_loc
            while(cur_dist > self.clumping_info["distance"] and low_bound < i_dx):
                snp_distance = i_dx - low_bound
                self.m_existed_snps[low_bound].m_clump_info.up_bound = i_dx

                low_bound+=1
                prev_loc = self.bim["pos"][low_bound]
                cur_dist = pos - prev_loc
            
            # src/genotype.cpp:111
            if (m_max_window_size < snp_distance): 
                m_max_window_size = snp_distance

            self.m_existed_snps[i_dx].m_clump_info.low_bound = low_bound
            self.m_existed_snps[i_dx].m_clump_info.up_bound = len(self.m_existed_snps)

        
        return m_max_window_size
    
    @cached_property
    def m_existed_snps(self):
        df = self.bim.merge(
                    self.base_data[["SNP", "P", "CHR", "BP"]], 
                    left_on="snp", 
                    right_on="SNP", 
                    how="left")\
                .reset_index()\
                .loc[:,["index", "chrom", "pos", "P"]]
        return [SNP(row.index, row.chrom, row.pos, row.P) for row in df.itertuples()]
    
    def update_index_tot(self):
        # inc/genotype.hpp:1118
        pass

    def get_r2(self, window_data_ptr):
        # inc/genotype.hpp:1149

        is_x = False
        counts = [0]*18
        freq11 = 0
        freq11_expected = 0
        freq1x = 0
        freq2x = 0
        freqx1 = 0
        freqx2 = 0
        dxx = 0

        # inc/genotype.hpp:1165
        plink_algorithm.genovec_3freq(counts, window_data_ptr)


        dxx = freq11 - freq11_expected
        r2 = dxx * dxx / (freq11_expected * freq2x * freqx2)


        pass

    def genotype_file_read(self, byte_pos, read_size):
        # inc/memoryread.hpp:14
        with open(self.genotype_filename, "rb") as f:
            f.seek(byte_pos)
            result = f.read(read_size)
        return result

    def read_genotype(self, snp, m_unfiltered_sample_ct):
        # inc/binaryplink.hpp:129
        final_mask = plink_algorithm.get_final_mask(self.num_samples, BITCT2)
        unfiltered_sample_ct4 = (m_unfiltered_sample_ct + 3) // 4           # 1 byte -> 8 bits -> 4 samples
        byte_pos = 3 + unfiltered_sample_ct4 * snp.index

        # inc/binaryplink.hpp:141
        snp.m_genotype_storage = self.genotype_file_read(byte_pos, unfiltered_sample_ct4)
        pass

    def clumping(self):
        # src/genotype.cpp:1194
        # src/genotype.cpp:1274
        min_r2 = 0.1
        m_founder_ct = len(self.fam)
        m_unfiltered_sample_ct = len(self.fam)
        founder_ctv3 = BITCT_TO_ALIGNED_WORDCT(m_founder_ct)
        founder_ctl2 = QUATERCT_TO_WORDCT(m_founder_ct)

        founder_ctsplit = 3 * founder_ctv3
        founder_ctv2 = QUATERCT_TO_ALIGNED_WORDCT(m_founder_ct)
        unfiltered_sample_ctl = BITCT_TO_WORDCT(m_unfiltered_sample_ct)

        unfiltered_sample_ctv2 = 2 * unfiltered_sample_ctl
        index_data = [0] * int(3 * founder_ctsplit + founder_ctv3)
        index_tots = [0] * 6

        founder_include2 = [0] * int(founder_ctv2)

        snp_ranges = self.get_chrom_boundary()            # src/genotype.cpp:1200

        num_snp_in_chr = sum([i.second - i.first for i in snp_ranges])
        max_snp_in_chr = 0

        # src/genotype.cpp:1300
        max_size = max(self.m_max_window_size, num_snp_in_chr * 0.01)
        r2 = -1
        num_processed = 0
        prev_processed = 0
        local_progress = 0
        prev_progress = 0
        local_num_core = 0
        for snp_range in snp_ranges:
            for i_snp in range(snp_range.first, snp_range.second):
                # src/genotype.cpp:1326
                core_snp_idx = self.m_sort_by_p_index[i_snp]
                core_snp = self.m_existed_snps[core_snp_idx]

                if (core_snp.clumped) or (core_snp.p_value > self.clumping_info["pvalue"]): 
                    continue

                # src/genotype.cpp:1330
                # inc/snp.hpp:351
                clump_start_idx = core_snp.m_clump_info.low_bound
                clump_end_idx = core_snp.m_clump_info.up_bound
                
                # src/genotype.cpp:1334
                for clump_idx in range(clump_start_idx, core_snp_idx):
                    clump_snp = self.m_existed_snps[clump_idx]

                    if (clump_snp.clumped) or (clump_snp.p_value > self.clumping_info["pvalue"]):
                        continue

                    # src/genotype.cpp:1341
                    # inc/snp.hpp:455
                    if (clump_snp.m_genotype_storage == 0):
                        # src/genotype.cpp:1344
                        self.read_genotype(clump_snp, m_founder_ct, m_unfiltered_sample_ct)
                        pass
                    pass
                
                # src/genotype.cpp:1360
                self.update_index_tot()     # Just relating to memory manage 

                for clump_idx in range(clump_start_idx, core_snp_idx):
                    clump_snp = self.m_existed_snps[clump_idx]

                    if(clump_snp.clumped or (clump_snp.p_value > self.clumping_info["pvalue"])):
                        continue
                    
                    # src/genotype.cpp:1374
                    r2 = self.get_r2(window_data_ptr)

                    pass



            pass

        pass

gt = Genotype(
        "tests/data/EUR.QC",
        "tests/data/Height.QC.gz",
        "tests/data/EUR.covariate"
    )

# print(gt.m_max_window_size)
# gt.clumping()
print(gt.G.compute()[1:3,:5])