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

# def calculate_effect_allele_count_na_hom_ref(vds):
#     """
#     Calculate the effect allele count from the given VariantDataset (VDS), handling NA and homozygous reference cases.
    
#     :param vds: Hail VariantDataset.
#     :return: Expression to compute effect allele count.
#     """
#     effect_allele = vds.prs_info['effect_allele']
#     non_effect_allele = vds.prs_info['noneffect_allele']
        
#     ref_allele = vds.alleles[0]

#     # Create a set of alternate alleles using hl.set
#     alt_alleles_set = hl.set(vds.alleles[1:].map(lambda allele: allele))

#     is_effect_allele_ref = ref_allele == effect_allele
#     is_effect_allele_alt = alt_alleles_set.contains(effect_allele)
#     is_non_effect_allele_ref = ref_allele == non_effect_allele
#     is_non_effect_allele_alt = alt_alleles_set.contains(non_effect_allele)

#     return hl.case() \
#         .when(hl.is_missing(vds.GT) & is_effect_allele_ref, 2) \
#         .when(hl.is_missing(vds.GT) & is_effect_allele_alt, 0) \
#         .when(vds.GT.is_hom_ref() & is_effect_allele_ref, 2) \
#         .when(vds.GT.is_hom_var() & is_effect_allele_alt, 2) \
#         .when(vds.GT.is_het() & is_effect_allele_ref, 1) \
#         .when(vds.GT.is_het() & is_effect_allele_alt, 1) \
#         .default(0)

def calculate_effect_allele_count_na_hom_ref(vds):
    """
    Computes the effect allele count per sample for a given variant in a Hail VariantDataset (VDS).
    Supports both VDS v7 (GT-based) and VDS v8 (LGT + LA-based) formats.

    This function is essential for PRS calculation. It resolves each sample's genotype into the
    number of effect alleles (0, 1, or 2), accounting for multiallelic sites and missing data.
    
    If the genotype is missing, the function assumes a homozygous reference genotype. This is consistent
    with how the All of Us VDS sparsely encodes genotype calls—only non-reference genotypes are stored,
    and missing values are interpreted as homozygous reference. For further discussion of this assumption
    and its implications for polygenic risk score calculation, see:

        Khattab et al., “AoUPRS: A Cost-Effective and Versatile PRS Calculator for the All of Us Program”,
        BMC Genomics 2025. https://link.springer.com/article/10.1186/s12864-025-11693-9
    """
    effect_allele = vds.prs_info['effect_allele']
    non_effect_allele = vds.prs_info['noneffect_allele']

    # Always row-scoped (safe for v7 and v8)
    ref_allele = vds.row.alleles[0]
    alt_alleles_set = hl.set(vds.row.alleles[1:])

    # Boolean flags to identify where the effect allele lies
    is_effect_allele_ref = ref_allele == effect_allele
    is_effect_allele_alt = alt_alleles_set.contains(effect_allele)

    if "GT" in vds.entry:  # v7 — dense GT genotypes
        return (
            hl.case()
              .when(hl.is_missing(vds.GT) & is_effect_allele_ref, 2)
              .when(hl.is_missing(vds.GT) & is_effect_allele_alt, 0)
              .when(vds.GT.is_hom_ref() & is_effect_allele_ref, 2)
              .when(vds.GT.is_hom_var() & is_effect_allele_alt, 2)
              .when(vds.GT.is_het() & (is_effect_allele_ref | is_effect_allele_alt), 1)
              .default(0)
        )
        
    else:  # v8 — reconstruct GT from LGT + LA
        lgt = vds.entry.LGT
        la = hl.or_else(vds.entry.LA, hl.empty_array(hl.tint32))

        alleles = hl.or_missing(
            hl.is_defined(lgt) & hl.is_defined(la),
            hl.array([
                vds.row.alleles[hl.or_else(la[lgt[0]], 0)],
                vds.row.alleles[hl.or_else(la[lgt[1]], 0)]
            ])
        )

        missing_expr = hl.case() \
            .when(hl.is_missing(alleles) & is_effect_allele_ref, 2) \
            .when(hl.is_missing(alleles) & is_effect_allele_alt, 0)

        effect_allele_count = hl.len(
            hl.or_else(alleles, hl.empty_array(hl.tstr)).filter(lambda a: a == effect_allele)
        )

        return missing_expr.default(effect_allele_count)
    
    
# Add more utility functions as needed
