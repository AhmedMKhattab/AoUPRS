# utils.py

import hail as hl

def calculate_effect_allele_count(mt):
    """
    Calculate the effect allele count from the given MatrixTable.
    
    :param mt: Hail MatrixTable.
    :return: Expression to compute effect allele count.
    """
    effect_allele = mt.prs_info['effect_allele']
    non_effect_allele = mt.prs_info['noneffect_allele']
        
    ref_allele = mt.alleles[0]

    # Create a set of alternate alleles using hl.set
    alt_alleles_set = hl.set(mt.alleles[1:].map(lambda allele: allele))

    is_effect_allele_ref = ref_allele == effect_allele
    is_effect_allele_alt = alt_alleles_set.contains(effect_allele)
    is_non_effect_allele_ref = ref_allele == non_effect_allele
    is_non_effect_allele_alt = alt_alleles_set.contains(non_effect_allele)

    return hl.case() \
        .when(mt.GT.is_hom_ref() & is_effect_allele_ref, 2) \
        .when(mt.GT.is_hom_var() & is_effect_allele_alt, 2) \
        .when(mt.GT.is_het() & is_effect_allele_ref, 1) \
        .when(mt.GT.is_het() & is_effect_allele_alt, 1) \
        .default(0)

def calculate_effect_allele_count_na_hom_ref(vds):
    """
    Calculate the effect allele count from the given VariantDataset (VDS), handling NA and homozygous reference cases.
    
    :param vds: Hail VariantDataset.
    :return: Expression to compute effect allele count.
    """
    effect_allele = vds.prs_info['effect_allele']
    non_effect_allele = vds.prs_info['noneffect_allele']
        
    ref_allele = vds.alleles[0]

    # Create a set of alternate alleles using hl.set
    alt_alleles_set = hl.set(vds.alleles[1:].map(lambda allele: allele))

    is_effect_allele_ref = ref_allele == effect_allele
    is_effect_allele_alt = alt_alleles_set.contains(effect_allele)
    is_non_effect_allele_ref = ref_allele == non_effect_allele
    is_non_effect_allele_alt = alt_alleles_set.contains(non_effect_allele)

    return hl.case() \
        .when(hl.is_missing(vds.GT) & is_effect_allele_ref, 2) \
        .when(hl.is_missing(vds.GT) & is_effect_allele_alt, 0) \
        .when(vds.GT.is_hom_ref() & is_effect_allele_ref, 2) \
        .when(vds.GT.is_hom_var() & is_effect_allele_alt, 2) \
        .when(vds.GT.is_het() & is_effect_allele_ref, 1) \
        .when(vds.GT.is_het() & is_effect_allele_alt, 1) \
        .default(0)
    
# Add more utility functions as needed