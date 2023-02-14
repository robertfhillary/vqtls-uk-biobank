##########################################
######### Variant Annotation #############
##########################################

# Load requisite libraries
library(tidyverse)
library(data.table)

# Set working directory 
setwd("/UKB_PPP_pQTLs/")

# Read in variant file 
b2=readRDS("/rhillary/dis_independent_annotated.rds") 
# Prepare list for annotation 
names(b2)[1]="ID"
query_vep_res = list()
# Loop stage 
for(i in seq(22)){
  query_vep_res[[i]] = fread(cmd = paste0("grep -v '^##' analysis/interim/VEP_annot/annotations/ukb_imp_chr", i, "_mac5_info03_b0_7_b38VEP.txt")) %>% 
    semi_join(b2, by = c("#Uploaded_variation" = "ID"))
# Print to denote completion 
  print(i)
}

# Tidy and format final file 
query_vep_res2 = bind_rows(query_vep_res)
names(query_vep_res2)[1]="SNP"
query_vep_res2$SNP=as.character(query_vep_res2$SNP)
anno=query_vep_res2[which(!duplicated(query_vep_res2$SNP)),]
# Merge back in with original file to ensure they have annotations 
b3=merge(b2,anno,by.x="ID",by.y="SNP")

# Obtain rsid information
chr_rsid_map = fread("/UKBB_imputed_hg38/RSID_MAP.tsv.gz")
chr_rsid_map$ID=paste(chr_rsid_map$ID,"v1",sep=":")
chr2=chr_rsid_map[which(chr_rsid_map$ID %in% b3$SNP),]
# Merge in with annotated information 
b4=merge(b3, chr2, by.x="SNP", by.y="ID")
