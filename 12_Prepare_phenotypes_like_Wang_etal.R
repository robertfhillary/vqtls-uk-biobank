###################################################################################
#### Prepare Phenotypes as recommended by Wang et al. (sensitivity analysis) ######
###################################################################################

# Load requisite libraries 
library(data.table)

# Define outlier removal function - 5 SD
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

## DISCOVERY 
# Read in the raw phenotype file 
prot=as.data.frame(fread("/UKB_PPP_pQTLs/analysis/input/dat_wide_internal_use.csv"))
# Read in information about the proteins i.e. which panel are they on
info=as.data.frame(fread("/UKB_PPP_pQTLs/analysis/input/maps/olink_protein_map_1.5k_v1.tsv"))
# Loop through different panels with different covariate files
vars=c("Cardiometabolic","Inflammation","Neurology","Oncology")
# Loop stage 
for(i in vars){ 
# Read in discovery covariate file 
cov=as.data.frame(fread(paste0("/UKB_PPP_pQTLs/analysis/input/UKBPPP_GWAS_RUN1/discovery_covars_", i, "_v1.tsv")))
# Subset protein dataframe to only those proteins on panel being queried in this stage of the loop 
queries=info[which(info$Panel %in% i),]
prot1=prot[,c(1,which(names(prot) %in% queries$UKBPPP_ProteinID))]
# Create temporary dataframe with individuals common to raw phenotype and covariate files
tmp=merge(prot1,cov,by.x="App_26041",by.y="FID")
# Extract column ids of the proteins so that we can use these in a regression loop 
ids=grep("v1", names(tmp))
# Also find index of panel-specific covariate - helps in regression step 
covind=which(names(tmp) %in% paste0(i,"_tbms"))
# Loop through the proteins in question 
for(j in ids){ 
# Regression step - want to regress each protein on all covariates stated 
tmp1=tmp[tmp$sex%in%1,]
tmp2=tmp[tmp$sex%in%0,]
resids1=as.data.frame(resid(lm(tmp1[,j] ~ age+tmp1[,covind]+UKBPC_1+UKBPC_2+UKBPC_3+UKBPC_4+UKBPC_5+UKBPC_6+UKBPC_7+UKBPC_8+UKBPC_9+UKBPC_10+UKBPC_11+UKBPC_12+UKBPC_13+UKBPC_14+UKBPC_15+UKBPC_16+UKBPC_17+UKBPC_18+UKBPC_19+UKBPC_20+factor(Batch)+factor(ukb_centre_fct)+factor(array_fct),data=tmp1,na.action=na.exclude)))
resids2=as.data.frame(resid(lm(tmp2[,j] ~ age+tmp2[,covind]+UKBPC_1+UKBPC_2+UKBPC_3+UKBPC_4+UKBPC_5+UKBPC_6+UKBPC_7+UKBPC_8+UKBPC_9+UKBPC_10+UKBPC_11+UKBPC_12+UKBPC_13+UKBPC_14+UKBPC_15+UKBPC_16+UKBPC_17+UKBPC_18+UKBPC_19+UKBPC_20+factor(Batch)+factor(ukb_centre_fct)+factor(array_fct),data=tmp2,na.action=na.exclude)))
# Tidy up dataframe 
names(resids1)[1]="resids"
resids1$FID=tmp1[,1]
resids1$IID=tmp1[,1]
resids1=resids1[,c(2,3,1)]
names(resids2)[1]="resids"
resids2$FID=tmp2[,1]
resids2$IID=tmp2[,1]
resids2=resids2[,c(2,3,1)]
# Remove outliers that are >5 SD from mean 
resids1$resids=outlierTrim(resids1$resids)
resids2$resids=outlierTrim(resids2$resids)
# Scale to mean zero and unit variance 
resids1$resids=scale(resids1$resids)
resids2$resids=scale(resids2$resids)
resids=rbind(resids1,resids2)
# Prepare name for writing out to match other files used in analyses 
nm=names(tmp1)[j]
nm=gsub(":", "_", nm)
nm1=paste0(nm,"_",i)
# Check length of file name - needs to be 5 
if(length(unlist(strsplit(nm1, "_"))) > 5){ 
nm2=gsub("^.*?_","",nm1)
write.table(resids, paste0("/rhillary/discovery_wangtransformed/", nm2, ".phen"), quote = F, row.names = F, sep = ' ')
print(nm2)
} else 
{ 
## Write out file in format required for OSCA 
write.table(resids, paste0("/rhillary/discovery_wangtransformed/", nm1, ".phen"), quote = F, row.names = F, sep = ' ')
print(nm1)
}
}
# Print to denote completion
print(i)
} 


## REPLICATION 
# Read in the raw phenotype file 
prot=as.data.frame(fread("/UKB_PPP_pQTLs/analysis/input/dat_wide_internal_use.csv"))
# Read in information about the proteins i.e. which panel are they on
info=as.data.frame(fread("/UKB_PPP_pQTLs/analysis/input/maps/olink_protein_map_1.5k_v1.tsv"))
# Loop through different panels with different covariate files
vars=c("Cardiometabolic","Inflammation","Neurology","Oncology")
# Loop stage 
for(i in vars){ 
# Read in replication covariate file 
cov=as.data.frame(fread(paste0("/UKB_PPP_pQTLs/analysis/input/UKBPPP_GWAS_RUN1/replication_covars_", i, "_v1.tsv")))
# Subset protein dataframe to only those proteins on panel being queried in this stage of the loop 
queries=info[which(info$Panel %in% i),]
prot1=prot[,c(1,which(names(prot) %in% queries$UKBPPP_ProteinID))]
# Create temporary dataframe with individuals common to raw phenotype and covariate files
tmp=merge(prot1,cov,by.x="App_26041",by.y="FID")
# Extract column ids of the proteins so that we can use these in a regression loop 
ids=grep("v1", names(tmp))
# Also find index of panel-specific covariate - helps in regression step 
covind=which(names(tmp) %in% paste0(i,"_tbms"))
# Loop through the proteins in question 
for(j in ids){ 
tmp1=tmp[tmp$sex%in%1,]
tmp2=tmp[tmp$sex%in%0,]
# Regression step - want to regress each protein on all covariates stated 
resids1=as.data.frame(resid(lm(tmp1[,j] ~ age+tmp1[,covind]+UKBPC_1+UKBPC_2+UKBPC_3+UKBPC_4+UKBPC_5+UKBPC_6+UKBPC_7+UKBPC_8+UKBPC_9+UKBPC_10+UKBPC_11+UKBPC_12+UKBPC_13+UKBPC_14+UKBPC_15+UKBPC_16+UKBPC_17+UKBPC_18+UKBPC_19+UKBPC_20+factor(Batch)+factor(ukb_centre_fct)+factor(array_fct),data=tmp1,na.action=na.exclude)))
resids2=as.data.frame(resid(lm(tmp2[,j] ~ age+tmp2[,covind]+UKBPC_1+UKBPC_2+UKBPC_3+UKBPC_4+UKBPC_5+UKBPC_6+UKBPC_7+UKBPC_8+UKBPC_9+UKBPC_10+UKBPC_11+UKBPC_12+UKBPC_13+UKBPC_14+UKBPC_15+UKBPC_16+UKBPC_17+UKBPC_18+UKBPC_19+UKBPC_20+factor(Batch)+factor(ukb_centre_fct)+factor(array_fct),data=tmp2,na.action=na.exclude)))
# Tidy up dataframe 
names(resids1)[1]="resids"
resids1$FID=tmp1[,1]
resids1$IID=tmp1[,1]
resids1=resids1[,c(2,3,1)]
names(resids2)[1]="resids"
resids2$FID=tmp2[,1]
resids2$IID=tmp2[,1]
resids2=resids2[,c(2,3,1)]
# Remove outliers that are >5 SD from mean 
resids1$resids=outlierTrim(resids1$resids)
resids2$resids=outlierTrim(resids2$resids)
# Scale to mean zero and unit variance 
resids1$resids=scale(resids1$resids)
resids2$resids=scale(resids2$resids)
resids=rbind(resids1,resids2)
# Prepare name for writing out to match other files used in analyses 
nm=names(tmp)[j]
nm=gsub(":", "_", nm)
nm1=paste0(nm,"_",i)
# Check length of file name - needs to be 5 
if(length(unlist(strsplit(nm1, "_"))) > 5){ 
nm2=gsub("^.*?_","",nm1)
write.table(resids, paste0("/rhillary/replication_wangtransformed/", nm2, ".phen"), quote = F, row.names = F, sep = ' ')
print(nm2)
} else 
{ 
## Write out file in format required for OSCA 
write.table(resids, paste0("/rhillary/replication_wangtransformed/", nm1, ".phen"), quote = F, row.names = F, sep = ' ')
print(nm1)
}
}
# Print to denote completion
print(i)
} 


