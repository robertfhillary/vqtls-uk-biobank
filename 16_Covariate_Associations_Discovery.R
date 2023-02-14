#####################################
###### COVARIATE ASSOCIATIONS #######
#####################################

# Load requisite libraries 
library(data.table)

# Read in protein data 
prot=as.data.frame(fread("/rhillary/dat_wide_internal_use.csv"))
names(prot)[1]="FID"

# Read in discovery and replication covariates 
cov=readRDS("/UKB_PPP_pQTLs/analysis/input/covars_nonscaled.rds")
cov[,c(5,9:28)] = apply(cov[,c(5,9:28)],2,as.numeric)
names(cov)[1]="FID"
names(cov)[9:28]=paste0("PC",1:20)

# Merge in correct ID format 
prot1=merge(prot,cov,by="FID")
# Add in remaining covariates 
extra=as.data.frame(fread("/rhillary/combined_covars_v1.tsv"))
extra1=extra[,c("FID","sample_selection","Batch")]
prot1=merge(prot1,extra1,by="FID",all.x=T)

# Subset to discovery individuals only 
dis=as.data.frame(fread("/rhillary/discovery_phenotypes_patched/AARSD1_Q9BTE6_OID21311_v1_Oncology.phen"))
prot1=prot1[which(prot1$FID %in% dis$FID),]
# Read in protein and panel information
info=as.data.frame(fread("/UKB_PPP_pQTLs/analysis/input/maps/olink_protein_map_1.5k_v1.tsv"))
info=info[which(!duplicated(info$UKBPPP_ProteinID)),]
# Prepare dictionaries of variable names 
# quantitative for linear regression 
vars=c("age","age2","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")
# binary variables for linear regression
bins=c("sex")
# interaction variables 
inter=c("age2_sex","age_sex")
# ANOVE variables 
anovs=c("Batch","ukb_centre_char","array_char")

#########################################################
########### 1. LINEAR REGRESSION ########################
#########################################################
# Create matrix to store output 
mat=as.data.frame(matrix(nrow=1472,ncol=(((length(vars)+1)*3)+5)))
names(mat)[1:5]=c("Protein","Panel","UniProt","OlinkID","Assay")
# Loop through proteins for regression 
for(i in 3:1474){
# Get protein name 
nm=names(prot)[i] 
# Get panel of the protein 
inf=info[which(info$UKBPPP_ProteinID %in% nm),]
panel=inf$Panel
pan=paste0(panel,"_tbms")
# Store protein and panel name 
mat[(i-2),1]=as.character(nm)
mat[(i-2),2]=as.character(panel) 
# Store UniProt, OlinkID and Assay name 
mat[(i-2),3]=as.character(inf$UniProt)
mat[(i-2),4]=as.character(inf$OlinkID) 
mat[(i-2),5]=as.character(inf$Assay)
# Add in panel-specific covariate to the covariates 
vars1=c(vars,pan)
# Loop through variables for linear regression
for(j in vars1){  
# Get index of variables 
ind=which(vars1%in%j)
# Store coefficients 
mat[(i-2),((3*(ind-0))+3)]=signif(summary(lm(scale(prot1[,i]) ~ scale(prot1[,which(names(prot1) %in% j)])))$coefficients[2,1],2)
mat[(i-2),((3*(ind-0))+4)]=signif(summary(lm(scale(prot1[,i]) ~ scale(prot1[,which(names(prot1) %in% j)])))$coefficients[2,2],2) 
mat[(i-2),((3*(ind-0))+5)]=signif(summary(lm(scale(prot1[,i]) ~ scale(prot1[,which(names(prot1) %in% j)])))$coefficients[2,4],2)
# Format column names 
names(mat)[((3*(ind-0))+3)]=paste0(j,".Beta")
names(mat)[((3*(ind-0))+4)]=paste0(j,".SE")
names(mat)[((3*(ind-0))+5)]=paste0(j,".P")
}
# Print to denote completion
print(i)
} 
# Store output for combining later 
linear1=mat

#########################################################
############## 2 . BINARY VARIABLES #####################
#########################################################
# Create new matrix for storing outputs with binary variables 
mat=as.data.frame(matrix(nrow=1472,ncol=(((length(bins)+0)*3)+5)))
names(mat)[1:5]=c("Protein","Panel","UniProt","OlinkID","Assay")
# Loop through proteins for regression 
for(i in 3:1474){
# Get protein name 
nm=names(prot)[i] 
# Get panel of the protein 
inf=info[which(info$UKBPPP_ProteinID %in% nm),]
panel=inf$Panel
# Store protein and panel name 
mat[(i-2),1]=as.character(nm)
mat[(i-2),2]=as.character(panel) 
# Store UniProt, OlinkID and Assay name 
mat[(i-2),3]=as.character(inf$UniProt)
mat[(i-2),4]=as.character(inf$OlinkID) 
mat[(i-2),5]=as.character(inf$Assay)
# Loop through variables
for(j in bins){  
# Get index of variables 
ind=which(bins%in%j)
# Store coefficients 
mat[(i-2),((3*(ind-0))+3)]=signif(summary(lm(scale(prot1[,i]) ~ prot1[,which(names(prot1) %in% j)]))$coefficients[2,1],2)
mat[(i-2),((3*(ind-0))+4)]=signif(summary(lm(scale(prot1[,i]) ~ prot1[,which(names(prot1) %in% j)]))$coefficients[2,2],2) 
mat[(i-2),((3*(ind-0))+5)]=signif(summary(lm(scale(prot1[,i]) ~ prot1[,which(names(prot1) %in% j)]))$coefficients[2,4],2)
# Format column names 
names(mat)[((3*(ind-0))+3)]=paste0(j,".Beta")
names(mat)[((3*(ind-0))+4)]=paste0(j,".SE")
names(mat)[((3*(ind-0))+5)]=paste0(j,".P")
}
# Print to denote completion
print(i)
} 
# Store output for combining later 
bins1=mat


#########################################################
############## 3 . ANOVA VARIABLES ######################
#########################################################
# Create new matrix for storing outputs with binary variables 
mat=as.data.frame(matrix(nrow=1472,ncol=(((length(anovs)+0)*2)+5)))
names(mat)[1:5]=c("Protein","Panel","UniProt","OlinkID","Assay")
# Loop through proteins for regression 
for(i in 3:1474){
# Get protein name 
nm=names(prot)[i] 
# Get panel of the protein 
inf=info[which(info$UKBPPP_ProteinID %in% nm),]
panel=inf$Panel
# Store protein and panel name 
mat[(i-2),1]=as.character(nm)
mat[(i-2),2]=as.character(panel) 
# Store UniProt, OlinkID and Assay name 
mat[(i-2),3]=as.character(inf$UniProt)
mat[(i-2),4]=as.character(inf$OlinkID) 
mat[(i-2),5]=as.character(inf$Assay)
# Loop through variables
for(j in anovs){  
# Get index of variables 
ind=which(anovs%in%j)
# Store coefficients 
mat[(i-2),((2*(ind-0))+4)]=signif(summary(aov(scale(prot1[,i]) ~ factor(prot1[,which(names(prot1) %in% j)])))[[1]][1,4],2)
mat[(i-2),((2*(ind-0))+5)]=signif(summary(aov(scale(prot1[,i]) ~ factor(prot1[,which(names(prot1) %in% j)])))[[1]][1,5],2)
# Format column names 
names(mat)[((2*(ind-0))+4)]=paste0(j,".F")
names(mat)[((2*(ind-0))+5)]=paste0(j,".P")
}
# Print to denote completion
print(i)
} 
# Store output for combining later 
anovs1=mat

#########################################################
############ 4 . INTERACTION VARIABLES ##################
#########################################################
# Create new matrix for storing outputs with binary variables 
mat=as.data.frame(matrix(nrow=1472,ncol=(((length(inter)+0)*3)+5)))
names(mat)[1:5]=c("Protein","Panel","UniProt","OlinkID","Assay")
# Loop through proteins for regression 
for(i in 3:1474){
# Get protein name 
nm=names(prot)[i] 
# Get panel of the protein 
inf=info[which(info$UKBPPP_ProteinID %in% nm),]
panel=inf$Panel
# Store protein and panel name 
mat[(i-2),1]=as.character(nm)
mat[(i-2),2]=as.character(panel) 
# Store UniProt, OlinkID and Assay name 
mat[(i-2),3]=as.character(inf$UniProt)
mat[(i-2),4]=as.character(inf$OlinkID) 
mat[(i-2),5]=as.character(inf$Assay)
# Age*Sex interaction  
mat[(i-2),6]=signif(summary(lm(scale(prot1[,i]) ~ scale(prot1$age)*prot1$sex))$coefficients[4,1],2)
mat[(i-2),7]=signif(summary(lm(scale(prot1[,i]) ~ scale(prot1$age)*prot1$sex))$coefficients[4,2],2)
mat[(i-2),8]=signif(summary(lm(scale(prot1[,i]) ~ scale(prot1$age)*prot1$sex))$coefficients[4,4],2)
# Format column names 
names(mat)[6]=paste0("Age_Sex",".Beta")
names(mat)[7]=paste0("Age_Sex",".SE")
names(mat)[8]=paste0("Age_Sex",".P")
# Age2*Sex interaction  
mat[(i-2),9]=signif(summary(lm(scale(prot1[,i]) ~ scale(prot1$age2)*prot1$sex))$coefficients[4,1],2)
mat[(i-2),10]=signif(summary(lm(scale(prot1[,i]) ~ scale(prot1$age2)*prot1$sex))$coefficients[4,2],2)
mat[(i-2),11]=signif(summary(lm(scale(prot1[,i]) ~ scale(prot1$age2)*prot1$sex))$coefficients[4,4],2)
# Format column names 
names(mat)[9]=paste0("Age2_Sex",".Beta")
names(mat)[10]=paste0("Age2_Sex",".SE")
names(mat)[11]=paste0("Age2_Sex",".P")
# Print to denote completion
print(i)
} 
# Store output for combining later 
inter1=mat

# Combine all 
total=cbind(linear1,inter1[,c(6:11)],bins1[,c(6,7,8)],anovs1[,c(6:13)])
# Save out file 
fwrite(total,"/rhillary/covariate_assocations_discovery_patched.csv",row.names=F)

