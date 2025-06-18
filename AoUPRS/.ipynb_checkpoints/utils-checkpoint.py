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

    Parameters
    ----------
    vds : hail.vds.VariantDataset
        A Hail VariantDataset containing PRS metadata in `vds.prs_info` and either `GT` (v7) or `LGT + LA` (v8).

    Returns
    -------
    hail.expr.Int32Expression
        An integer Hail expression representing the effect allele dosage (0–2) for each sample.
    """
    # Retrieve the effect and non-effect alleles from PRS weights
    effect_allele = vds.prs_info['effect_allele']
    non_effect_allele = vds.prs_info['noneffect_allele']

    # Reference allele is always the first allele in the global allele list
    ref_allele = vds.alleles[0]
    alt_alleles_set = hl.set(vds.alleles[1:])

    # Boolean flags to identify where the effect allele lies
    is_effect_allele_ref = ref_allele == effect_allele
    is_effect_allele_alt = alt_alleles_set.contains(effect_allele)

    if 'GT' in vds.entry:
        # VDS v7 — use GT (global genotypes) directly
        gt_alleles = hl.or_else(vds.GT, hl.null(hl.tcall))
        alleles = hl.or_missing(
            hl.is_defined(gt_alleles),
            gt_alleles.alleles().map(lambda i: vds.alleles[i])  # map index → actual allele
        )
    else:
        # VDS v8 — use LGT + LA to reconstruct genotypes from local alleles
        lgt = vds.entry.LGT
        la = hl.or_else(vds.entry.LA, hl.empty_array(hl.tint32))

        # Manually map local genotype indices → global alleles via LA → vds.alleles
        alleles = hl.or_missing(
            hl.is_defined(lgt) & hl.is_defined(la),
            hl.array([
                vds.row.alleles[hl.or_else(la[lgt[0]], 0)],
                vds.row.alleles[hl.or_else(la[lgt[1]], 0)]
            ])
        )

    # Handle cases where genotype is missing:
    # If missing and effect allele is reference → assume 2 copies
    # If missing and effect allele is alternate → assume 0 copies
    missing_expr = hl.case() \
        .when(hl.is_missing(alleles) & is_effect_allele_ref, 2) \
        .when(hl.is_missing(alleles) & is_effect_allele_alt, 0)

    # Final effect allele dosage: count how many genotype alleles match the effect allele
    effect_allele_count = hl.len(
        hl.or_else(alleles, hl.empty_array(hl.tstr)).filter(lambda a: a == effect_allele)
    )

    return missing_expr.default(effect_allele_count)
    
    
# Add more utility functions as needed