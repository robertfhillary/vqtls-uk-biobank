#################################################################
######## Reverse Look-Up Analysis: meQTL to vQTL ################
#################################################################

# Load requisite libraries 
library(data.table)

# Read in pQTLs or main effect QTLs 
qtl=read.csv("/rhillary/inputs/main_effect_qtls.csv")

# Subset to MAF > 0.05 
qtls1=qtl[-which(qtl$A1FREQ..discovery. < 0.05),]
qtls1=qtls1[-which(qtls1$A1FREQ..discovery. > 0.95),] #8878 qtls in 1355 proteins 

# Extract vQTL files and set working directory 
setwd("/rhillary/residualised_vqtls/")
loop=list.files(".", ".vqtl")

# Tidy main effect/pQTL file 
names(qtls1)[2]="Chr"
# Identify and fix problematic file names 
qtls1[grep("IL12A_IL12B_P29459_P29460_OID21327_v1",qtls1$UKBPPP.ProteinID),"UKBPPP.ProteinID"]="IL12B_P29459_OID21327_v1"
qtls1[grep("CKMT1A_CKMT1B_P12532_OID20721_v1",qtls1$UKBPPP.ProteinID),"UKBPPP.ProteinID"]="CKMT1B_P12532_OID20721_v1"
qtls1[grep("DEFA1_DEFA1B_P59665_OID20344_v1",qtls1$UKBPPP.ProteinID),"UKBPPP.ProteinID"]="DEFA1B_P59665_OID20344_v1"
qtls1[grep("DEFB4A_DEFB4B_O15263_OID21373_v1",qtls1$UKBPPP.ProteinID),"UKBPPP.ProteinID"]="DEFB4B_O15263_OID21373_v1"
qtls1[grep("EBI3_IL27_Q14213_Q8NEV9_OID21389_v1",qtls1$UKBPPP.ProteinID),"UKBPPP.ProteinID"]="IL27_Q14213_OID21389_v1"
qtls1[grep("FUT3_FUT5_P21217_Q11128_OID21013_v1",qtls1$UKBPPP.ProteinID),"UKBPPP.ProteinID"]="FUT5_P21217_OID21013_v1"
qtls1[grep("LGALS7_LGALS7B_P47929_OID21406_v1",qtls1$UKBPPP.ProteinID),"UKBPPP.ProteinID"]="LGALS7B_P47929_OID21406_v1"
qtls1[grep("MICB_MICA_Q29980_Q29983_OID20593_v1",qtls1$UKBPPP.ProteinID),"UKBPPP.ProteinID"]="MICA_Q29980_OID20593_v1"

# Set up naming system common to both files
qtls1$query=paste(paste0("chr",qtls1$Chr),qtls1$UKBPPP.ProteinID,sep="_")

# Create list to store output
list1=list()
# Loop through each unique protein/chromosome combination and extract vQTL summary statistics 
for(i in unique(qtls1$query)){
# Subset file to just that protein and chromosome 
tmp=qtls1[qtls1$query %in% i,]	
# Read in vqtl file 
file=loop[grep(i,loop)]
# Check if file is missing - will suggest mismatched names 
if(is.na(file)){print(paste(i, "missing"))}else{
# Read in vQTL file and subset to SNPs in main effect GWAS
file1=as.data.frame(fread(file))
file2=file1[which(file1$SNP %in% tmp$Variant.ID..CHROM.GENPOS..hg37..A0.A1.imp.v1.),]
# Subset main effect QTL file to vQTLs to ensure consistency 
tmp2=tmp[which(tmp$Variant.ID..CHROM.GENPOS..hg37..A0.A1.imp.v1. %in% file2$SNP),]
# Add in vQTL summary stats 
if(nrow(tmp2)==0){ 
NULL 
} else { 
# Match order of SNPs 
ids=file2$SNP 
tmp2=tmp2[match(ids,tmp2$Variant.ID..CHROM.GENPOS..hg37..A0.A1.imp.v1.),]
tmp2$vQTL_beta=file2$beta
tmp2$vQTL_se=file2$se
tmp2$vQTL_p = file2$P
# Store output 
list1[[i]]=tmp2
}
}
print(i)
}

# Combine outputs and format 
list2=as.data.frame(do.call("rbind",list1))

# Tabulate information
length(which(list2$vQTL_p < 0.05)) 
length(which(list2$vQTL_p < 5e-8)) 
length(which(list2$vQTL_p < 3.4e-11)) 

# Identify final list 
overlap=list2[which(list2$vQTL_p < 3.4e-11),]


