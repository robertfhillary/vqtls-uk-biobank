########################################################################################################
#### Prepare Proteins adjusting for top-associated phenotypes in GEI tests (sensitivity analysis) ######
########################################################################################################

# Load requisite libraries 
library(data.table)

# Create functions 
# Create function for outlier removal - 4SD
outlierID <- function(x, cut=5) {
  xx <- scale(x)
  retval <- ifelse (abs(xx) >= cut, TRUE, FALSE)
  retval
}
outlierTrim <- function(x, cut=5) {
  id <- outlierID(x, cut=cut)
  retval <- ifelse(id, NA, x)
}

# Rank-inverse normal based transformation 
INT <- function(x) qnorm((rank(x, na.last="keep") - 0.5) / sum(!is.na(x)))

## DISCOVERY 
# Read in the raw phenotype file 
prot=as.data.frame(fread("/UKB_PPP_pQTLs/analysis/input/dat_wide_internal_use.csv"))
# Read in information about the proteins i.e. which panel are they on
info=as.data.frame(fread("/UKB_PPP_pQTLs/analysis/input/maps/olink_protein_map_1.5k_v1.tsv"))
# Loop through different panels with different covariate files
vars=c("Cardiometabolic","Inflammation","Neurology","Oncology")
# Subset to proteins with conditional GEI associations 
cond=fread("/rhillary/outputs/cond_gei_final2023.csv")
# Align protein formats 
cond$Protein=gsub("_",":", cond$Protein)
# Remove everything after "v1"
cond$Protein= sub(":[^:]+$", "", cond$Protein)
cond[which(cond$Protein %in% "CKMT1B:P12532:OID20721:v1"),"Protein"]="CKMT1A_CKMT1B:P12532:OID20721:v1"
info=info[which(!duplicated(info$UKBPPP_ProteinID)),]
info=info[which(info$UKBPPP_ProteinID %in% cond$Protein),]
# Read in GEI phenotypes 
phenos=as.data.frame(fread("/rhillary/interaction_phenos_new.txt"))
phenos=phenos[which(!duplicated(phenos$IND)),]
# Subset to those with conditionally significant associations with vQTL-protein pairs 
phenos1=phenos[,c(1,which(names(phenos) %in% unique(cond$Phenotype)))]
names(phenos1)[1]="FID"
# Loop stage 
for(i in 1:nrow(info)){
# Extract protein of interest 
tmp.df=info[i,] 
# Read in discovery covariate file 
cov=as.data.frame(fread(paste0("/UKB_PPP_pQTLs/analysis/input/UKBPPP_GWAS_RUN1/discovery_covars_", tmp.df$Panel, "_v1.tsv")))
# Subset protein dataframe to only those proteins on panel being queried in this stage of the loop 
queries=info[which(info$Panel %in% i),]
prot1=prot[,c(1,which(names(prot) %in% tmp.df$UKBPPP_ProteinID))]
# Create temporary dataframe with individuals common to raw phenotype and covariate files
tmp=merge(prot1,cov,by.x="App_26041",by.y="FID")
# Subset conditionally significant phenotypes to those involving the protein 
cond.df=cond[which(cond$Protein %in% tmp.df$UKBPPP_ProteinID),]
phenos.tmp=phenos1[,c(1,which(names(phenos1) %in% cond.df$Phenotype))]
# Extract names of phenotypes for regressions later 
phens=names(phenos.tmp)[2:ncol(phenos.tmp)]
# Merge phenotype(s) in with protein and covariate file 
tmp=merge(tmp,phenos.tmp,by.x="App_26041",by.y="FID")
# Prepare variables for regression 
names(tmp)=gsub(" ", "_", names(tmp))
names(tmp)=gsub("\\(", "", names(tmp))
names(tmp)=gsub("\\)", "", names(tmp))
tmp$sex=as.factor(tmp$sex)
tmp$Batch=as.factor(tmp$Batch)
tmp$ukb_centre_fct=as.factor(tmp$ukb_centre_fct)
tmp$array_fct=as.factor(tmp$array_fct)
tmp[,2]=INT(tmp[,2])
names(tmp)[2]="prot"
tmp$IID=NULL 
# Extract variable names 
xs=paste(names(tmp)[3:ncol(tmp)],collapse="+")
form1=as.formula(sprintf(paste("prot","~",xs)))
model <- lm(form1,data=tmp,na.action=na.exclude)
# Regression step 
resids1=as.data.frame(resid(model))
# Tidy up dataframe 
names(resids1)[1]="resids"
resids1$FID=tmp[,1]
resids1$IID=tmp[,1]
resids1=resids1[,c(2,3,1)]
# Scale to mean zero and unit variance 
resids1$resids=scale(resids1$resids)
# Prepare name for writing out to match other files used in analyses 
nm=names(prot1)[2]
nm=gsub(":", "_", nm)
nm1=paste0(nm,"_",tmp.df$Panel)
# Check length of file name - needs to be 5 
if(length(unlist(strsplit(nm1, "_"))) > 5){ 
nm2=gsub("^.*?_","",nm1)
names(resids1)[3]=nm2
write.table(resids1, paste0("/home/rhillary/discovery_gei_phenotypes_sensitivity/", nm2, ".phen"), quote = F, row.names = F, sep = ' ')
} else 
{ 
## Write out file in format required for OSCA 
write.table(resids1, paste0("/home/rhillary/discovery_gei_phenotypes_sensitivity/", nm1, ".phen"), quote = F, row.names = F, sep = ' ')
names(resids1)[3]=nm1}
# Print to denote completion
print(i)
print(phens)
} 

# Prepare files to help with bash scripts 
# We need CHR and Protein for each vQTL involved 
cond$test=paste(cond$SNP, cond$Protein, sep="_")
vqtl=readRDS("/rhillary/outputs/dis_independent_annotated.rds")
vqtl$compare1=gsub("_",":", vqtl$compare)
vqtl$compare1= sub(":[^:]+$", "", vqtl$compare1)
vqtl[which(vqtl$compare1 %in% "CKMT1B:P12532:OID20721:v1"),"compare1"]="CKMT1A_CKMT1B:P12532:OID20721:v1"
vqtl$test=paste(vqtl$SNP, vqtl$compare1, sep="_")
vqtls1=vqtl[vqtl$test %in% cond$test,]
# Get CHRs 
chrs=as.data.frame(vqtls1[,c("Chr")])
names(chrs)[1]="chr"
phenos=as.data.frame(vqtls1[,c("compare")])
names(phenos)[1]="phen"
phenos$phen=paste0("/rhillary/discovery_gei_phenotypes_sensitivity/", phenos$phen, ".phen")
# Write out tables 
write.table(chrs,"/rhillary/gei_sensitivity_chrs.txt", row.names=F, col.names=F, quote = F)
write.table(phenos,"/rhillary/gei_sensitivity_phenos.txt", row.names=F, col.names=F, quote = F)



