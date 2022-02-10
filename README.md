# SIgHeterosis
The directory contains all essential scripts to perform the analyses in the manuscript. 
  1_RNAseqPipeline: Raw RNA-sequencing data in fastq format are mapped to reference.
  2_Postalignment_processing: Reads mapped to targets that may introduce artefacts in abundance estimate (e.g. chloroplast and rDNA clusters) were removed and the Transcript Per Million (TPM) re-estimated. Outlier samples are also dropped.
  3_find_Non-additive_genes.r: Genes whose expression deviate systematically from mid-parent value are identified (SHB1 data set).
  4_SHB1clusterHvI_BayesianLMMspline.r: clusters of non-additive genes were modeled for the relationships between gene expression level and rosette size of plants (SHB1 data set).
  5_SHB2lm_additive_MPHpredict.r: Genes whose expression in F1 are systematically similar to that of Mid-Parent Value are identified (SHB2 data set).
  6_SHB2Call_BTHresponsive_genes.r: Genes whose expression changes significantly in response to BTH treatment are identified (SHB2 data set).
  7_lmmTreatment_responsive_genes_bsplinecluster.r: BTH-responsive genes identified in the last step are clustered based on how their expression deviation from parental mean in F1s is associated with rosette size heterosis.
  8_MPHpredictCluster_Wilcoxontest.r: individual clusters acquired from the previous steps are tested for whether the two most extreme decile of F1 expression deviation is associated with significantly different rosette size heterosis.
  9_Subgenelist_Geneticdistance.r: whether inter-parental genetic (Hamming) distance is associated with hybrid performance. Various inter-parental genetic distances are calculated by using common SNPs in the whole genome, and subset of SNPs in the gene body and flanking regions of a) common additive genes, b) genes whose non-additive expression in F1 are negatively associated with heterosis, and c) genes whose non-additive expression in F1 are positively associated with heterosis.
  10_Rareallele_xpn.r, 10A_Rareallele_xpn_PosNeggenelist.r, 10B_Rareallele_xpn_additivegenelist.r: whether rare allele count in regulatory region of genes affects the gene expression, and whether gene expression in turn associates with plant growth. It is also tested whether the association differs between hybrids and inbreds using various (additive and non-additive) gene lists.
  
  For detailed explanation, see #README sections of each script.
  