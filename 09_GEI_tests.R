###############################################
######## GEI association tests ################
###############################################

# Load requisite libraries 
library(data.table)

# Create functions for downstream analyses 
# Create function for outlier removal - 5SD
outlierID <- function(x, cut=5) {
  xx <- scale(x)
  retval <- ifelse (abs(xx) >= cut, TRUE, FALSE)
  retval
}
outlierTrim <- function(x, cut=5) {
  id <- outlierID(x, cut=cut)
  retval <- ifelse(id, NA, x)
}

# Create function for rank-based inverse normal transformation 
INT <- function(x) qnorm((rank(x, na.last="keep") - 0.5) / sum(!is.na(x)))


###########################################
##### Phenotye preparation stage ##########
###########################################

# Read in phenotypes 
phen=as.data.frame(fread("/rhillary/all_quant_20200904.pheno.v2.tsv.gz"))
# Subset to Olink individuals (discovery set first)
ids=fread("/rhillary/discovery_phenotypes_patched/ACE2_Q9BYF1_OID20105_v1_Cardiometabolic.phen",header=T)
phen=phen[phen$IND %in% ids$FID,]

# Add in education phenotype - coded as to UKB standard 
ed=as.data.frame(fread("/home/rhillary/education_separate.csv"))
ed$V1=NULL
sub1 <- ed[,c(1,2)]
sub1 <- sub1[which(sub1[,2] %in% 1),]
names(sub1)[2] <- 'X'
sub2 <- ed[,c(1,8)]
sub2 <- sub2[which(sub2[,2] %in% 1),]
names(sub2)[2] <- 'X'
sub3 <- ed[,c(1,14)]
sub3 <- sub3[which(sub3[,2] %in% 1),]
names(sub3)[2] <- 'X'
sub4 <- ed[,c(1,20)]
sub4 <- sub4[which(sub4[,2] %in% 1),]
names(sub4)[2] <- 'X'
bind <- rbind(sub1, sub2)
bind <- rbind(bind, sub3)
bind <- rbind(bind, sub4)
names(bind)[1] <- 'SampleID'
bind <- bind[which(bind$SampleID %in% phen[,1]),]
length(unique(bind[,1])) 
names(ed) <- c('SampleID', 'Educ') # 16,216
names(ed)[1] <- 'SampleID'
ed <- ed[which(ed$SampleID %in% phen[,1]),]
ed$Edu <- ifelse(ed$SampleID %in% bind$SampleID, 1, 0)
ed <- ed[c(1,26)]
names(ed)="IND"


# Add in additional lifestyle covariates 
covs=as.data.frame(fread("/rhillary/covs_separate.csv"))
covs$f.21001.0.0=NULL
covs[which(covs$f.20116.0.0%in%-3),"f.20116.0.0"]=NA
covs[which(covs$f.1558.0.0%in%-3),"f.1558.0.0"]=NA
covs=covs[which(covs[,1] %in% phen$IND),]
covs$f.189.0.0=INT(covs$f.189.0.0)
names(covs)[1]="IND"


# Add in age and sex at baseline 
f=as.data.frame(fread("/rhillary/combined_covars_v1.tsv"))
f=f[,c(1,3,4)]
f=f[which(f[,1] %in% phen$IND),]
f$age=INT(f$age)
names(f)[1]="IND"
covs1=merge(f,covs,by="IND")
covs1=merge(ed,covs1,by="IND")

# Format files and add in season of blood draw 
names(covs1)=c("IND","Education (College Degree)", "Age", "Sex", "Deprivation (Townsend)", "Alcohol Intake Frequency", "Smoking Status")
ses=as.data.frame(fread("ukb_ppp_blood_collection_season_from_janssen.csv"))
a=as.data.frame(fread("/UKB_PPP_pQTLs/analysis/input/maps/olink_sample_map_projectinfo_v2.tsv") )
ses1=merge(ses, a[,c("UKBPPP_SampleID","App_26041")], by="UKBPPP_SampleID")
ses1$UKBPPP_SampleID=NULL
names(ses1)[3]="FID"
covs1=merge(covs1, ses1[,c("FID", "blood_season")], by.x="IND",by.y="FID")
names(covs1)[ncol(covs1)]="Season of Blood Draw"

# Remove phenotypes with too much missingess (>20%)
cols=which(apply(phen,2,function(a) (length(which(is.na(a)))/length(a)>0.20))==TRUE)
phen1=phen[,-which(names(phen) %in% names(cols))]

# Remove phenotypes with insufficient heterogeneity 
het=as.data.frame(apply(phen1, 2, function(a) length(unique(a[!is.na(a)]))))
names(het)[1]="count"
het$col=row.names(het)
het=het[order(het$count),]
phen2=phen1[,which(!names(phen1) %in% het[het$count <=1, "col"])]

# Remove extreme outliers >5 SDs
phen2[,2:ncol(phen2)]=apply(phen2[,2:ncol(phen2)], 2, function(a) outlierTrim(a))

# INVT step of continuous variables 
phen2[,2:ncol(phen2)]=apply(phen2[,2:ncol(phen2)], 2, function(a) INT(a))

# Remove additional columns post QC 
phen2$f_21_0_0=NULL
phen2$f_23_0_0=NULL

# Save out file 
fwrite(phen2, "/rhillary/interaction_phenos.txt",row.names=F)
phen3=merge(phen2,covs1,by="IND")
phen3$Sex=NULL
phen3$Age=NULL 
fwrite(phen3, "/rhillary/interaction_phenos_new.txt",row.names=F)

###############################
## QTLs - vQTL and meQTL ######
###############################

############################################################################
#### PREPARE FILES FOR GEI TESTS - GET QTLs and THEIR PROTEIN(S) ###########
############################################################################

### Subset recoded dosgae plink files to just SNPs in vQTLs or main effect QTLs ###

# Get vQTLs 
vqtls=readRDS("/rhillary/outputs/dis_independent_annotated.rds")
# Get pQTLs or main effect QTLs 
pqtls=read.csv("/rhillary/inputs/main_effect_qtls.csv")
# Subset pQTLs to MAF > 5% 
pqtls=pqtls[which(pqtls$A1FREQ..discovery. > 0.05),]
# Extract SNPs and Chromosome  
vqtls=vqtls[,c(1,2)]
names(vqtls)=c("SNP","Chr")
vqtls$Type="vqtl"
pqtls=pqtls[,c(1,2)]
names(pqtls)=c("SNP","Chr")
pqtls$Type="pqtl"
# Combine for query
tmp=rbind(vqtls,pqtls)
# Loop through each chromsome to extract SNPs of interest
for(i in 1:22){
snps=tmp[which(tmp$Chr %in% i),]
# Read in each vQTL and its associated phenotype 
x=as.data.frame(fread(paste0("/rhillary/recode/chr",i,".raw")))
names(x)=gsub("_.*", "", names(x))
x2=x[,c(1,2,which(names(x) %in% snps$SNP))]
x2$IID=NULL
# Find missing SNPs 
missed=snps[which(!snps$SNP %in% names(x2)),]
fwrite(missed, paste0("/rhillary/interaction_chrs/missing_chr",i,".txt"),row.names=F)
# Save out new file 
fwrite(x2, paste0("/rhillary/interaction_chrs/chr",i,".txt"),row.names=F)
rm(x)
rm(x2)
gc()
# Print to denote completion
print(i)
} 

# Combine files 
setwd("/rhillary/interaction_chrs/")
loop=list.files("/rhillary/interaction_chrs/", ".")
a=lapply(loop,fread)
a1=as.data.frame(do.call("cbind", a))
ids=a1[,1]
a2=a1[,which(!names(a1) %in% names(a1)[grep("FID", names(a1))])]
a2$FID=ids 
fwrite(a2, "/rhillary/interaction_chrs/complete_snps.txt", row.names=F)



# Combine proteins with respective associated SNP as phenotype file 
# Prepare vqtl and pqtl data 
vqtls=readRDS("/rhillary/outputs/dis_independent_annotated.rds")
# pQTLs or main effect QTLs 
qtls1=read.csv("/rhillary/inputs/main_effect_qtls.csv")
qtls1=qtls1[which(qtls1$A1FREQ..discovery. > 0.05),]
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
# Extract SNPs and Chromosome  and filename
vqtls=vqtls[,c(1,2,17)]
names(vqtls)=c("SNP","Chr","file")
vqtls$Type="vqtl"
qtls1=qtls1[,c(1,2,ncol(qtls1))]
names(qtls1)=c("SNP","Chr","file")
qtls1$Type="pqtl"
## Additional QC 
vqtls$compare=sub(".*?_", "", vqtls$file)
qtls1$compare=sub(".*?_", "", qtls1$file)
maps=as.data.frame(fread("/rhillary/olink_protein_map_1.5k_cleaned.tsv"))
maps=maps[-which(duplicated(maps$UKBPPP_ProteinID)),]
maps$UKBPPP.ProteinID=gsub(":","_",maps$UKBPPP_ProteinID)
maps[grep("IL12A_IL12B_P29459_P29460_OID21327_v1",maps$UKBPPP.ProteinID),"UKBPPP.ProteinID"]="IL12B_P29459_OID21327_v1"
maps[grep("CKMT1A_CKMT1B_P12532_OID20721_v1",maps$UKBPPP.ProteinID),"UKBPPP.ProteinID"]="CKMT1B_P12532_OID20721_v1"
maps[grep("DEFA1_DEFA1B_P59665_OID20344_v1",maps$UKBPPP.ProteinID),"UKBPPP.ProteinID"]="DEFA1B_P59665_OID20344_v1"
maps[grep("DEFB4A_DEFB4B_O15263_OID21373_v1",maps$UKBPPP.ProteinID),"UKBPPP.ProteinID"]="DEFB4B_O15263_OID21373_v1"
maps[grep("EBI3_IL27_Q14213_Q8NEV9_OID21389_v1",maps$UKBPPP.ProteinID),"UKBPPP.ProteinID"]="IL27_Q14213_OID21389_v1"
maps[grep("FUT3_FUT5_P21217_Q11128_OID21013_v1",maps$UKBPPP.ProteinID),"UKBPPP.ProteinID"]="FUT5_P21217_OID21013_v1"
maps[grep("LGALS7_LGALS7B_P47929_OID21406_v1",maps$UKBPPP.ProteinID),"UKBPPP.ProteinID"]="LGALS7B_P47929_OID21406_v1"
maps[grep("MICB_MICA_Q29980_Q29983_OID20593_v1",maps$UKBPPP.ProteinID),"UKBPPP.ProteinID"]="MICA_Q29980_OID20593_v1"
maps$compare=paste(maps$UKBPPP.ProteinID, maps$Panel,sep="_")
maps1=maps[maps$UKBPPP.ProteinID %in% unique(qtls1$compare),]
qtls2=merge(qtls1, maps1[,c("UKBPPP.ProteinID", "compare")], by.x="compare", by.y="UKBPPP.ProteinID")
qtls2$file=paste0(qtls2$compare.y,".phen")
qtls=qtls2[,c(2,3,4,5)]
# Find which ones do not make final file
missed=as.data.frame(fread("/rhillary/interaction_chrs/complete_snps.txt"))
qtls=qtls[which(qtls$SNP %in% names(missed)),] # 8857 



# Loop stage to combine vQTLs and their proteins 
for(i in 1:nrow(vqtls)){
tmp=vqtls[i,]
file=tmp[,"compare"]
snp=tmp[,"SNP"]
x=as.data.frame(fread("/rhillary/interaction_chrs/complete_snps.txt"))
y=read.table(paste0("/rhillary/discovery_phenotypes_patched/",file),header=T)
x2=x[,c(ncol(x),which(names(x) %in% snp))]
x3=merge(x2,y[,c(1,3)],by="FID")
# Create name 
nm1=gsub(".phen","",tmp$compare)
names(x3)[3]=nm1
# Write out file 
nm=paste(i, snp, file,sep="_")
fwrite(x3, paste0("/rhillary/interaction_vqtls/", nm),row.names=F)
# Print to denote completion 
print(i)
} 


# Loop stage to combine main effect or pQTLs and their proteins 
for(i in 1:nrow(qtls)){
tmp=qtls[i,]
file=tmp[,"file"]
snp=tmp[,"SNP"]
x=as.data.frame(fread("/scratch/rhillary/interaction_chrs/complete_snps.txt"))
y=read.table(paste0("/rhillary/discovery_phenotypes_patched/",file),header=T)
x2=x[,c(ncol(x),which(names(x) %in% snp))]
x3=merge(x2,y[,c(1,3)],by="FID")
# Create name
nm1=gsub(".phen","",tmp$file)
names(x3)[3]=nm1
# Write out file 
nm=paste(i, snp, file,sep="_")
fwrite(x3, paste0("/rhillary/interaction_pqtls/", nm),row.names=F)
# Print to denote completion 
print(i)
} 
