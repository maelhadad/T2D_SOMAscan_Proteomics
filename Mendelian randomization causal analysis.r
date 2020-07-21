# Two-samples bi-directional Mendelian randomization: 
## Using summary statistics of published GWAS studies of type 2 diabetes (T2D)
## and GWAS studies of SOMAscan measured proteins, we applied the causal analysis
## two-samples Mendelian randomization:

### Install the TwoSampleMR R package: full documentation could be accessed here:
### https://mrcieu.github.io/TwoSampleMR/articles/introduction.html
devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)

################################################################################
## Download the summary statistics data:
### Download the T2D GWAS by Xue et al. 2018, Nature Communications
### Title: "Genome-wide association analyses identify 143 risk variants 
### and putative regulatory mechanisms for type 2 diabetes"
### at that link: "http://cnsgenomics.com/data/t2d/Xue_et_al_T2D_META_Nat_Commun_2018.gz"

### Download the SOMAscan measured protein GWAS study by Sun et al., Nature, 2018
### Title: "Genomic atlas of the human plasma proteome" 
### at that link: "http://www.phpc.cam.ac.uk/ceu/proteins/"

### Download the SOMAscan measured protein GWAS study by Suhre et al., Nature Communications, 2017
### Title: "Connecting genetic risk to disease end points through the human blood plasma proteome" 
### at that link: "http://metabolomics.helmholtz-muenchen.de/pgwas/"

################################################################################
### First direction: T2D --> Proteins: here T2D would be the exposure in the MR analysis.

#### Load the Xue GWAS summary stats from the download directory: 
t2d_gw_all <- readr::read_delim(file = gzfile("Xue_et_al_T2D_META_Nat_Commun_2018.gz"),
                                delim = " ",col_names = T)

#### Step 1: Subset to the GWAS significant results:
t2d_gw_sig <- t2d_gw_all[t2d_gw_all$P<5e-08,]
t2d_gw_sig$pheno <- "T2D"

t2d_exp_data<-format_data(dat = t2d_gw_sig,type = "exposure",snp_col = "SNP",
                          beta_col = "b",se_col = "se",eaf_col = "frq_A1",
                          effect_allele_col = "A1",other_allele_col = "A2",
                          pval_col = "P",samplesize_col = "N",chr_col = "CHR",
                          pos_col = "BP",phenotype_col = "pheno")

###  Step 2: Clump the results to get a list of indepenedent SNPs:
t2d_exp_data<-clump_data(t2d_exp_data)

###  Step 3: Checking if there are any ambiguous palindromic SNPs:
### Replace ambiguous and palindromic SNPs with non ambiguous ones:
apsnps<-t2d_exp_data[TwoSampleMR:::check_palindromic(A1 = t2d_exp_data[,"effect_allele.exposure"],
                                                     A2 = t2d_exp_data[,"other_allele.exposure"]),"SNP"]

### Step 4: Read protein summary statistics from Sun et al.: 
### (Use Suhre et al. for proteins not found in Sub et al.)
protein_sumstat <- protein_sumstat[protein_sumstat$SNP %in% t2d_exp_data$SNP,]

### Step 5: for SNPs missing in protein summary statistics extract LD proxies
### use data extracted from the online tool: "https://ldlink.nci.nih.gov/?tab=ldproxy"
snps_for_lookup <- t2d_exp_data$SNP[t2d_exp_data$SNP %in% protein_sumstat$SNP]

### add LD proxies if available.

### Step 6: format the protein summary statistics:
#### if protein_sumstat from Sun et al.
protein_sumstat$N<-3301
outcome_dat <- format_data(dat = protein_sumstat,type = "outcome",id_col = "ID",
                         snps = exposure_dat$SNP,snp_col = "SNPID",
                         beta_col = "Effect",se_col = "StdErr",eaf_col = "eaf",
                         effect_allele_col = "Allele1",other_allele_col = "Allele2",
                         pval_col = "Pval",chr_col = "chromosome",
                         pos_col = "position",samplesize_col = "N")
#### if protein_sumstat from Suhre et al.
outcome_dat <- format_data(dat = protein_sumstat,type = "outcome",id_col = "ID",
                         snps = exposure_dat$SNP,snp_col = "SNP",se_col = "SE",
                         beta_col = "BETA",z_col = "STAT",eaf_col = "eaf",
                         effect_allele_col = "A1",other_allele_col = "A2",
                         pval_col = "P",chr_col = "CHR",pos_col = "BP",samplesize_col = "NMISS")

### Step 7: Harmonize data.
dat<-harmonise_data(exposure_dat = t2d_exp_data,outcome_dat = outcome_dat)

### Step 8: Perform MR.
mr_r<-mr(dat)

### Step 9: Sensitivity analyses:
heterogeneity<-mr_heterogeneity(dat)
pleiotropy<-mr_pleiotropy_test(dat)
res_single<-mr_singlesnp(dat)
res_loo<-mr_leaveoneout(dat)

################################################################################
### Second direction: Protein --> T2D: here the protein would be the exposure in the MR analysis.
### Apply the same steps starting with the extraction of the IVs for the protein of interest.

#### for proteins with no IVs in Sun et al., IVs were extracted from Emilsson et al.,Science, 2018
#### results of which could be downloaded at: "www.sciencemag.org/content/361/6404/769/suppl/DC1"

#### Step 1: Subset to the GWAS significant results:
protein_gw_sig <- protein_gw_sig[protein_gw_sig$P<5e-08,]
protein_gw_sig$pheno <- "protein"
protein_gw_sig$N <- 3301
protein_exp_data<-format_data(dat = protein_sumstat,type = "outcome",id_col = "ID",
                              snps = exposure_dat$SNP,snp_col = "SNPID",
                              beta_col = "Effect",se_col = "StdErr",eaf_col = "eaf",
                              effect_allele_col = "Allele1",other_allele_col = "Allele2",
                              pval_col = "Pval",chr_col = "chromosome",
                              pos_col = "position",samplesize_col = "N")

#### Step 2: subset to only SNPs in cis i.e. within 1 MB:
library("biomaRt")
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",GRCh=37)

##### Extract info on the protein's Gene:
info<-sapply(unique(protein_exp_data[,"KORA.EntrezGeneSymbol"]),function(i){
  getBM(attributes=c('ensembl_gene_id','gene_biotype','hgnc_symbol','chromosome_name',
                     'start_position','end_position'),
        filters = 'hgnc_symbol', values =i, mart = ensembl)
},USE.NAMES = T,simplify = F)

info<-plyr::rbind.fill(info)

info$min<-pmax(info$start_position-1000000,1)
info$max<-info$end_position+1000000

for(i in 1:nrow(protein_exp_data)){
  s<-protein_exp_data[i,"SNP"]
  g<-protein_exp_data[i,"KORA.EntrezGeneSymbol"]
  info2<-info[info$hgnc_symbol==g,]
  samechr<-ifelse(protein_exp_data[i,"Chr"] == info2[,"chromosome_name"],"Yes","No")
  
  if(samechr=="Yes"){
    protein_exp_data[i,"cis_trans"]<-
      ifelse(protein_exp_data[i,"Position"] > info2$min & 
               protein_exp_data[i,"Position"] < info2$max,"cis","trans")
  }else{
    protein_exp_data[i,"cis_trans"]<-"trans"
  }
}

protein_exp_data <- protein_exp_data[protein_exp_data$cis_trans %in% "cis",]

#### Step 3: Clump the results to get a list of indepenedent SNPs:
protein_exp_data<-clump_data(protein_exp_data)

##### Step 4: Checking if there are any ambiguous palindromic SNPs:
### Replace ambiguous and palindromic SNPs with non ambiguous ones:
apsnps<-protein_exp_data[TwoSampleMR:::check_palindromic(A1 = protein_exp_data[,"effect_allele.exposure"],
                                                     A2 = protein_exp_data[,"other_allele.exposure"]),"SNP"]

### Step 5: Read t2d summary statistics from Xue et al.: 
t2d_sumstat <- t2d_sumstat[t2d_sumstat$SNP %in% protein_exp_data$SNP,]

### Step 6: for SNPs missing in protein summary statistics extract LD proxies
### use data extracted from the online tool: "https://ldlink.nci.nih.gov/?tab=ldproxy"
snps_for_lookup <- t2d_sumstat$SNP[t2d_sumstat$SNP %in% protein_exp_data$SNP]

### add LD proxies if available.

### Step 7: format the t2d summary statistics:
outcome_dat <- format_data(dat = outcome_dat,type = "exposure",snp_col = "SNP",
                             beta_col = "b",se_col = "se",eaf_col = "frq_A1",
                             effect_allele_col = "A1",other_allele_col = "A2",
                             pval_col = "P",samplesize_col = "N",chr_col = "CHR",
                             pos_col = "BP",phenotype_col = "pheno")

### Step 8: Harmonize data.
dat<-harmonise_data(exposure_dat = protein_exp_data,outcome_dat = outcome_dat)

### Step 9: Perform MR.
mr_r<-mr(dat)

### Step 10: Sensitivity analyses:
heterogeneity<-mr_heterogeneity(dat)
pleiotropy<-mr_pleiotropy_test(dat)
res_single<-mr_singlesnp(dat)
res_loo<-mr_leaveoneout(dat)
