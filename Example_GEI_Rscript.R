## Extract taskid from SLURM job 
taskid <- commandArgs(trailingOnly=TRUE)
print(as.numeric(taskid))
taskid=as.numeric(taskid)

# Load requisite libraries 
library(data.table)

# Read in phenotypes 
phenos=as.data.frame(fread("/rhillary/interaction_phenos_new.txt"))
# QC 
names(phenos)[1]="FID"
phenos$f_3089_0_0=NULL
phenos$f_3090_0_0=NULL
phenos$f_30230_0_0=NULL
# Extract protein files along with their vQTL info 
loop=list.files("/rhillary/interaction_vqtls/", ".")
# Loop through each vqtl file 
tmp=as.data.frame(fread(paste0("/rhillary/interaction_vqtls/",loop[[taskid]])))
# Combine datasets 
tmp1=merge(tmp,phenos,by="FID")
# Create output matrix 
output=as.data.frame(matrix(nrow=length(4:ncol(tmp1)),ncol=6))
names(output)=c("Protein","SNP","Phenotype","Interaction_Beta","Interaction_SE","Interaction_P")
for(j in 4:ncol(tmp1)){ 
summary(lm(tmp1[,3] ~ tmp1[,2]*tmp1[,j]))
# Store variable names
output[(j-3),1]=names(tmp1)[3]
output[(j-3),2]=names(tmp1)[2]
output[(j-3),3]=names(tmp1)[j]
# Store coefficients
output[(j-3),4:6]=summary(lm(tmp1[,3] ~ tmp1[,2]*tmp1[,j]))$coefficients[4,c(1,2,4)]
}
# Combine outputs and save out 
output1=as.data.frame(output)
write.table(output1, paste0("/rhillary/interaction_vqtls_output/", loop[[taskid]]), row.names=F)


