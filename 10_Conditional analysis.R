##########################################################
####### STEPWISE CONDITIONAL GEI ANALYSES ################
##########################################################

#### Just first stage shown for brevity as an example but same principle applies at each stage #####

# Load requisite libraries
library(data.table)

# Set working directory 
setwd("/rhillary/interaction_vqtls/")

# Read in GEI results 
gei=read.csv("/rhillary/outputs/gei_updated.csv")
gei[gei$Protein %in% "CKMT1A_CKMT1B_P12532_OID20721_v1", "Protein"] = "CKMT1B.P12532.OID20721.v1"
# Read in necessary input files 
phenos=as.data.frame(fread("/rhillary/interaction_phenos_new.txt"))
# Additional QC 
names(phenos)[1]="FID"
phenos$f_3089_0_0=NULL
phenos$f_3090_0_0=NULL
phenos$f_30230_0_0=NULL


# Extract files in directory 
loop=list.files(".",".")
# Set up list to store outputs 
list1=list()
gei$test=paste(gei$Protein, gei$SNP, sep = "_")
# Loop through proteins 
for(i in unique(gei$test)){ 
# Subset to protein of interest 
tmp=gei[gei$test %in% i,]
# Remove duplicates 
tmp=tmp[which(!duplicated(tmp$Phenotype)),]
# Ensure most significant interaction is ranked first 
tmp=tmp[order(tmp$Interaction_P) ,]
# If only one interaction present, skip 
if(nrow(tmp)==1){ NULL }else { 
top.pheno=tmp[1,"Phenotype"]
# Loop through phenotypes
output=as.data.frame(matrix(nrow=(nrow(tmp)-1),ncol=7))
names(output)=c("Protein","SNP","Phenotype","Conditioned_On", "Interaction_Beta","Interaction_SE","Interaction_P") 
for(j in 2:nrow(tmp)){ 
tmp1=tmp[j,]
# Extract file name 
test=paste(tmp1$SNP, tmp1$Protein,sep="_")
# Read in SNP and protein file 
file=as.data.frame(fread(loop[grep(test,loop)]))
# Merge with phenotypes 
phenos1=merge(phenos, file, by="FID")
# Extract phenotype of interest 
phen=tmp1[1,"Phenotype"]
# Extract SNP 
snp=tmp1[1,"SNP"]
# Extract protein name 
prot=tmp1[1,"Protein"]
# Run interaction now including the most significant phenotype as a covariate 
mod=lm(phenos1[,prot] ~ phenos1[,snp]*phenos1[,phen] + phenos1[,snp]*phenos1[,top.pheno])
summary(mod)$coefficients[5,c(1,2,4)]
output[(j-1),1]=prot
output[(j-1),2]=snp
output[(j-1),3]=phen
output[(j-1),4]=top.pheno
# Store coefficients
output[(j-1),5:7]=summary(mod)$coefficients[5,c(1,2,4)]
} 
# Tidy up and store output 
output=as.data.frame(output)
list1[[i]]=output
}
# Print to denote completion
print(paste0("Completed ", ind=which(unique(gei$test) %in% i), "/", length(unique(gei$test)), " Proteins"))
} 

# Combine final list 
list2=as.data.frame(do.call("rbind",list1))
list3=list2[which(list2$Interaction_P < 5.4e-6),]
gei1=gei[which(!duplicated(gei$test)),]



############# REPEAT UNTIL NO PHENOTYPES REMAIN #################
