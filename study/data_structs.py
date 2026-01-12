from dataclasses import dataclass, field


@dataclass(frozen=True)
class SnpRange:
    first: int 
    second: int

@dataclass
class Clumping:     # inc/storage.hpp:193
    r2 = 0.1
    proxy = 0.0
    pvalue = 1
    distance = 250000
    no_clump = False
    use_proxy = False
    provided_distance = False

@dataclass
class SNPClump:     # inc/storage.hpp:126
    low_bound: int = 0
    up_bound: int = 0
    max_flag_idx: int = 0
    clumped: bool = False

@dataclass
class SNP:
    chr: int 
    pos: int
    p_value: float 
    m_clump_info: SNPClump = field(default_factory=SNPClump)
    clumped: bool = False
    m_genotype_storage: int = 0