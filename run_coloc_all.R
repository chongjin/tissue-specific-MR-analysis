args = commandArgs(trailingOnly=TRUE)

# 5 combinations of IV and coloc datasets;
# 2 tissues with expression for coloc;
# 22 chromosomes.
if (length(args) == 3) {
  coloc_indices = as.numeric(args[1])
  tissue_indices = as.numeric(args[2])
  chrs = as.numeric(args[3])
}

# get exposure GWAS
library(tidyverse)
library(TwoSampleMR)
library(ieugwasr)
library(coloc)
# ao = available_outcomes()

IV_dataset = c("Savage", "Hill", "Okbay", "OkbaySavage", "OkbaySavage")[coloc_indices]
coloc_dataset =  c("Savage", "Savage", "Okbay", "Savage", "Okbay")[coloc_indices]
tissues = c("BrainMeta_cis_eqtl_summary", "cage_eqtl_data")[tissue_indices]

eQTL_folder = "data"
analysis_folder = sprintf("analysis_IV_%s_coloc_%s_%s", IV_dataset, coloc_dataset, tissues)
if(!file.exists(analysis_folder)) {
  dir.create(analysis_folder)
}

minp_QTL_cutoff = 1e-4  # The existence of a small QTL p-value as a prerequisite for coloc analysis
cis_cutoff = 2e6    # unused; 2e6 as the cutoff between cis- and trans- eQTL

# instrumental variables from GWAS
exposure_ld_clumped_file = file.path(analysis_folder, "IV_dataset.csv")
if (!file.exists(exposure_ld_clumped_file)) {
  # exposure_dat = extract_instruments("ebi-a-GCST006250", r2=0.01, kb=10000, p1=5e-8, p2=5e-8)
  # # > a = extract_instruments("ebi-a-GCST006250", r2=0.01, kb=10000, p1=1e-4)
  # # > dim(a)
  # # [1] 1489   15
  # # > b = extract_instruments("ebi-a-GCST006250", r2=0.01, kb=10000, p1=5e-3)
  # # > dim(b)
  # # [1] 6180   15
  # for (chr in 1:22) {
  #   exposure_dat_per_chr = exposure_dat[exposure_dat$chr.exposure == chr, ]
  #   cat(sprintf("chr%s: %s\n", chr, min(diff(sort(exposure_dat_per_chr$pos.exposure)))))
  # }
  # write_csv(exposure_dat, exposure_ld_clumped_file)
  # #
  # # Note: the LD clumping is equivalent to:
  # # # devtools::install_github("explodecomputer/plinkbinr")
  # Using the Intelligence GWAS form Hill et al. based on the multi-trait method
  # # # local LD clumping. See https://mrcieu.github.io/ieugwasr/articles/local_ld.html
  if (IV_dataset == "Savage") {
    exposure_dat = extract_instruments("ebi-a-GCST006250", r2=0.01, kb=10000, p1=1e-5)
    # > a = extract_instruments("ebi-a-GCST006250", r2=0.01, kb=10000, p1=1e-5)
    # > dim(a)
    # [1] 709    15
    # > a = extract_instruments("ebi-a-GCST006250", r2=0.01, kb=10000, p1=1e-4)
    # > dim(a)
    # [1] 1489   15
    # > b = extract_instruments("ebi-a-GCST006250", r2=0.01, kb=10000, p1=5e-3)
    # > dim(b)
    # [1] 6180   15
  } else if (IV_dataset == "Hill") {
    raw_exp_data = readxl::read_xlsx(file.path(eQTL_folder, "41380_2017_1_MOESM2_ESM.xlsx"), skip=1)
    exposure_dat_renamed = dplyr::rename(raw_exp_data,
                                         SNP='Tagged SNPs',
                                         IndSNP='Independent significant SNPs',
                                         beta='Beta',
                                         se='SE of tagged SNP',
                                         eaf='MAF of tagged SNP',
                                         effect_allele='effect allele',
                                         other_allele='non-effect allele',
                                         pval='P-value of tagged SNP',
                                         chr='Chromosome',
                                         pos='bp position of the tagged SNPs')
    clumped_exp_data = filter(exposure_dat_renamed,SNP==IndSNP)
    clumped_exp_data1 = mutate(clumped_exp_data,samplesize=248482,id='anderson.2020')
    exposure_dat = format_data(clumped_exp_data1)
    
  } else if (IV_dataset == "Okbay") {
    raw_exp_data = read_tsv(file.path(eQTL_folder, "Okbay_2022", "EA4_additive_p1e-5_clumped.txt"))
    exposure_dat_renamed = dplyr::rename(raw_exp_data,
                                         SNP=rsID,
                                         beta=Beta,
                                         se=SE,
                                         eaf=EAF_HRC,
                                         effect_allele=Effect_allele,
                                         other_allele=Other_allele,
                                         pval=P,
                                         chr=Chr,
                                         pos=BP)
    clumped_exp_data = mutate(exposure_dat_renamed, rsid=SNP, samplesize=3037499, id='EA4_additive')
    clumped_exp_data = ld_clump(clumped_exp_data, clump_r2=0.001) # If r2=0.01, there are 2950 vs. 746 SNPs after clumping
    exposure_dat = format_data(clumped_exp_data) 
  } else if (IV_dataset == "OkbaySavage") {
    raw_exp_data = read_tsv(file.path(eQTL_folder, "Okbay_2022", "EA4_additive_p1e-5_clumped.txt"))
    exposure_dat_Okbay = dplyr::rename(raw_exp_data,
                                         SNP=rsID,
                                         beta=Beta,
                                         se=SE,
                                         eaf=EAF_HRC,
                                         effect_allele=Effect_allele,
                                         other_allele=Other_allele,
                                         pval=P,
                                         chr=Chr,
                                         pos=BP)
    exposure_dat_Savage = extract_instruments("ebi-a-GCST006250", r2=0.01, kb=10000, p1=1e-5)
    exposure_dat_Savage = dplyr::rename(exposure_dat_Savage,
                                        SNP=SNP,
                                        beta=beta.exposure,
                                        se=se.exposure,
                                        eaf=eaf.exposure,
                                        effect_allele=effect_allele.exposure,
                                        other_allele=other_allele.exposure,
                                        pval=pval.exposure,
                                        chr=chr.exposure,
                                        pos=pos.exposure)

    common_columns = c("SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "pval", "chr", "pos")
    clumped_exp_data1 = rbind(exposure_dat_Okbay[,common_columns], exposure_dat_Savage[,common_columns])
    clumped_exp_data1 = mutate(clumped_exp_data1, rsid=SNP, samplesize=3037499, id='EA4_additive')
    clumped_exp_data2 = ld_clump(clumped_exp_data1, clump_r2=0.001) # If r2=0.01, there are 2950 vs. 746 SNPs after clumping
    exposure_dat = format_data(clumped_exp_data2)
    
    # set the IVs by taking the SNPs separately from Savage and Okbay out of "exposure_dat"
    # exposure_dat_Okbay_filtered = exposure_dat_Okbay[exposure_dat_Okbay$SNP %in% exposure_dat$SNP, ]
    
    if (coloc_dataset == "Okbay") {
      exposure_full = read_tsv("data/Okbay_2022/EA4_additive_excl_23andMe.txt.gz")
      exposure_dat_Okbay_filtered = exposure_full[match(exposure_dat$SNP, exposure_full$rsID) , ]
      exposure_dat_Okbay_filtered = dplyr::rename(exposure_dat_Okbay_filtered,
                                         SNP=rsID,
                                         beta=Beta,
                                         se=SE,
                                         eaf=EAF_HRC,
                                         effect_allele=Effect_allele,
                                         other_allele=Other_allele,
                                         pval=P,
                                         chr=Chr,
                                         pos=BP)
      exposure_dat_Okbay_filtered$Phenotype = "Okbay_2022"
      exposure_dat = format_data(exposure_dat_Okbay_filtered)
      
    } else if (coloc_dataset == "Savage") {
      exposure_full = read_tsv("data/29942086-GCST006250-EFO_0004337-Build37.f.tsv.gz")
      exposure_dat_Savage_filtered = exposure_full[match(exposure_dat$SNP, exposure_full$variant_id) , ]
      exposure_dat_Savage_filtered = dplyr::rename(exposure_dat_Savage_filtered,
                                    SNP=variant_id,
                                    beta=beta,
                                    se=standard_error,
                                    eaf=eaf_ref,
                                    effect_allele=effect_allele,
                                    other_allele=other_allele,
                                    pval=p_value,
                                    chr=chromosome,
                                    pos=base_pair_location)
      exposure_dat_Savage_filtered$Phenotype = "Savage_2018"
      exposure_dat = format_data(exposure_dat_Savage_filtered)
    }
  }
  write_csv(exposure_dat, exposure_ld_clumped_file)
} else {
  exposure_dat = read_csv(exposure_ld_clumped_file)
}

# Education GWAS
# https://thessgac.com/papers/14
# coloc: EA4_additive_excl_23andMe.txt.gz
# IV: EA4_additive_p1e-5_clumped.txt

# Read Education GWAS summary statistics. We will use it in coloc.abf().
max_bp = 2e5  # The value used in the 2014 coloc paper and the Leyden et al. 2022 paper
if (coloc_dataset == "Okbay") {
  exposure_full = read_tsv("data/Okbay_2022/EA4_additive_excl_23andMe.txt.gz")
  N_exposures = 765283
} else if (coloc_dataset == "Savage") {
  exposure_full = read_tsv("data/29942086-GCST006250-EFO_0004337-Build37.f.tsv.gz")
  N_exposures = max(exposure_full$n_analyzed, na.rm=TRUE)
  exposure_full = dplyr::rename(exposure_full, 
                                Chr=chromosome,
                                BP=base_pair_location,
                                SE=standard_error,
                                Beta=beta,
                                EAF_HRC=eaf_ref)
}

# exposure_dat = ld_clump(
#   exposure_dat_before_clump,
#   clump_kb = 10000,
#   clump_r2 = 0.01,
#   clump_p  = 5e-08,
#   plink_bin = plinkbinr::get_plink_exe(),
#   bfile = "1000G_reference/EUR"
# )

# load gene coordinates
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene

# Run coloc package and select the SNPs colocalized with Intelligence GWAS using
# multi-tissue eQTL data.
# The eQTL data SNPs are in GRCh37
# For each tissue, find the key. 
# For example, "ForBOX_CSF_afterQC_featureFile.csv".
# With the EntrezGeneSymbol, we can get the gene locations.
# Each expression corresponds to a separate file of summary statistics.
number_of_leading_SNPs = nrow(exposure_dat)

for (tissue in tissues) {
  
  # get geneid_keys, a table of probes
  if (tissue == "cage_eqtl_data") {
    expression_table = read_tsv("data/cage_eqtl_data/CAGE.sparse.summary", skip=9)
  } else if (tissue == "BrainMeta_cis_eqtl_summary") {
    expression_tables = list()
    for (chr in 1:22) {
      expression_tables[[chr]] = read_tsv(sprintf("data/BrainMeta_cis_eqtl_summary/BrainMeta_cis_eQTL_chr%d.summary", chr), skip=9)
    }
    expression_table = do.call(rbind, expression_tables)
  }
  
  expression_table$ProbeID = sapply(str_split(expression_table$`{ProbeID, ProbeChr, ProbeBP}`,
                                              pattern = regex("\\,|\\{|\\}")),
                                    function(x) x[2])
  expression_table$Chr = sapply(str_split(expression_table$`[Chr,cis_startBP,cis_endBP,NumSNPs]`,
                                          pattern = regex("\\,|\\[|\\]")),
                                function(x) as.numeric(x[2]))
  expression_table$cis_startBP = sapply(str_split(expression_table$`[Chr,cis_startBP,cis_endBP,NumSNPs]`,
                                                  pattern = regex("\\,|\\[|\\]")),
                                        function(x) as.numeric(x[3]))
  expression_table$cis_endBP = sapply(str_split(expression_table$`[Chr,cis_startBP,cis_endBP,NumSNPs]`,
                                                pattern = regex("\\,|\\[|\\]")),
                                      function(x) as.numeric(x[4]))
  minp_QTL_matrix = coloc_matrix = matrix(NA, nrow=nrow(expression_table), ncol=number_of_leading_SNPs)
  
  rownames(minp_QTL_matrix) = rownames(coloc_matrix) = expression_table$ProbeID
  colnames(minp_QTL_matrix) = colnames(coloc_matrix) = exposure_dat$SNP
  expression_table = expression_table[expression_table$Chr %in% 1:22, ]
  
  expression_table = expression_table[expression_table$Chr %in% chrs, ]
  
  for (j in seq_len(nrow(expression_table))) {
    
    cat(date(), sprintf("Processing coloc of expression %d of chr %d in tissue %s\n", j, chrs, tissue))
    
    # Get a list of nearby lead SNPs and check whether the lead and the expression is in proximity.
    candidate_lead_SNPs = which((exposure_dat$chr.exposure == expression_table$Chr[j]) &
                                (exposure_dat$pos.exposure >= expression_table$cis_startBP[j]) &
                                (exposure_dat$pos.exposure <= expression_table$cis_endBP[j]))
    
    if (length(candidate_lead_SNPs) > 0) {
    
      f = paste0("temp_", chrs, tissue, coloc_indices)
      if (tissue == "cage_eqtl_data") {
        system(
          sprintf("/home/chongjin/mtmr/smr_v1.3.1_linux_x86_64_static/smr --beqtl-summary data/cage_eqtl_data/CAGE.sparse --query 1 --snp-chr %d --from-snp-kb %d --to-snp-kb %d --probe %s --out %s",
                  expression_table$Chr[j],
                  max(0, floor((expression_table$cis_startBP[j] - max_bp) / 1000)),
                  ceiling((expression_table$cis_endBP[j] + max_bp) / 1000),
                  expression_table$ProbeID[j],
                  f))
      
      } else if (tissue == "BrainMeta_cis_eqtl_summary") {
        system(
          sprintf("/home/chongjin/mtmr/smr_v1.3.1_linux_x86_64_static/smr --beqtl-summary data/BrainMeta_cis_eqtl_summary/BrainMeta_cis_eQTL_chr%d --query 1 --snp-chr %d --from-snp-kb %d --to-snp-kb %d --probe %s --out %s",
                  expression_table$Chr[j],
                  expression_table$Chr[j],
                  max(0, floor((expression_table$cis_startBP[j] - max_bp) / 1000)),
                  ceiling((expression_table$cis_endBP[j] + max_bp) / 1000),
                  expression_table$ProbeID[j],
                  f))
      }
      if (!file.exists(paste0(f, ".txt"))) {
         next
      }
      QTLs = read_tsv(paste0(f, ".txt"))
      unlink(paste0(f, ".txt"))
      for (i in candidate_lead_SNPs) {
        
        # SNPs from the pre-LD clumping Intelligence GWAS at most 1m base pairs apart:
        exposure_local_SNPs = which((exposure_full$Chr == exposure_dat$chr.exposure[i]) &
                                    (abs(exposure_full$BP - exposure_dat$pos.exposure[i]) < max_bp) &
                                    !is.na(exposure_full$SE))
  
        QTL_local_SNPs = which((QTLs$Chr == exposure_dat$chr.exposure[i]) &
                               (abs(QTLs$BP - exposure_dat$pos.exposure[i]) < max_bp) &
                               !is.na(QTLs$SE))
        
        # length mismatch
        
        # For each expression, judging if there are any significant eQTLs within the region:
        # We proceed to select the intersection of
        # all the SNPs in the 1mb window from both the (not already LD clumped) GWAS and
        # the eQTL summary statistics of the expression we've found.
        common_SNP_positions = intersect(exposure_full$BP[exposure_local_SNPs],
                                         QTLs$BP[QTL_local_SNPs])
        
        exposure_local_common_SNPs = exposure_local_SNPs[match(common_SNP_positions,
                                                               exposure_full$BP[exposure_local_SNPs])]
        QTL_local_common_SNPs = QTL_local_SNPs[match(common_SNP_positions,
                                                     QTLs$BP[QTL_local_SNPs])]
        minp_QTL_matrix[j, i] = min(QTLs$p[QTL_local_common_SNPs], na.rm=TRUE)
        
        
        if (minp_QTL_matrix[j, i] < minp_QTL_cutoff) {
          # For a lead SNP, for all the associated expressions, run coloc:
          # colocalization using coloc::coloc.abf()
          # Take the cutoff of PPA4 > 0.8 as strong evidence of colocalization
          # between any of the Intelligence GWAS SNPs and the expression,
          # as in Leyden et. al. "proximal gene within a 200kb window".
          
          # supply a list of beta, varbeta, type="quant" to coloc.abf()
          # sample size; either do "--show-n" or check https://yanglab.westlake.edu.cn/software/smr/#DataResource
          if (tissue == "cage_eqtl_data") {
            N_QTLs = 2765
          } else if (tissue == "BrainMeta_cis_eqtl_summary") {
            N_QTLs = 2865
          }
          
          QTLs_coloc = list(beta = QTLs$b[QTL_local_common_SNPs],
                            varbeta = QTLs$SE[QTL_local_common_SNPs]^2,
                            N = N_QTLs,
                            type = "quant")
          
          exposures_coloc = list(beta = exposure_full$Beta[exposure_local_common_SNPs], 
                                 varbeta = exposure_full$SE[exposure_local_common_SNPs]^2,
                                 N = N_exposures,
                                 type = "quant")
          
          # Note that we do not have MAF for the QTLs; we are using the same MAF for both of QTL and the exposure:
          
          coloc_result = coloc.abf(QTLs_coloc, exposures_coloc, MAF = exposure_full$EAF_HRC[exposure_local_common_SNPs])
          
          coloc_matrix[j, i] = coloc_result$summary["PP.H4.abf"]
        }
        
        # Add to the SNP list for each tissue if we have evidence of colocalization.
      }
    }
  }
  # write_csv(as.data.frame(minp_QTL_matrix), file.path(analysis_folder, sprintf("minp_QTL_matrix_expression_%d_%s.csv", chrs, tissue)))
  # write_csv(as.data.frame(coloc_matrix), file.path(analysis_folder, sprintf("coloc_matrix_expression_%d_%s.csv", chrs, tissue)))
  
  write_csv(data.frame(SNP = exposure_dat$SNP, is_coloc = apply(coloc_matrix, 2, function(x) {any(na.omit(x) > 0.8)})),
            file.path(analysis_folder, sprintf("coloc_vector_expression_%s.csv", chrs)))
}
