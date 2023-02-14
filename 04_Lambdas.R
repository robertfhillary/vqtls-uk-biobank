############################################################
########### Lambdas - vQTL association studies #############
############################################################

# Load requisite libraries 
library(data.table)

# Set working directory 
setwd("/rhillary/residualised_vqtls/")
# Extract results files 
loop=list.files(".",".vqtl")
loop1=sub(".*?_", "", loop)
# Ensure 22 chromosomes are present 
loop2=as.data.frame(table(loop1))
loop3=loop2[which(loop2$Freq == 22),]
loop1=loop1[which(loop1 %in% loop3$loop1),]
loop1=unique(loop1)
# Set up dataframe to store results 
mat=as.data.frame(matrix(nrow=1472,ncol=7))
names(mat)=c("Protein","Lambda","Genome.Wide.Significant","Bonferroni.Significant","UniProt","Array","Filename")
# Loop step 
for(i in loop1){ 
# Get index of file in loop
ind=which(loop1%in%i)
# Extract file name 
file=gsub(".vqtl","",i)
sp=strsplit(file,split="_")[[1]]
# Extract all the files for this protein 
tmp.files=loop[grep(i, loop)]
# Set up list to store output 
list1=list()
for(j in tmp.files){ 
# Read in file 
tmp=as.data.frame(fread(j))
list1[[i]]=tmp
# Store results 
} 
# Prepare file for calculation 
list2=as.data.frame(do.call("rbind",list1))
# Genomic Inflation Factor 
chisq <- qchisq(1-list2$P,1)
lambda=median(chisq)/qchisq(0.5,1)
# Get number of genome-wide significant hits 
hits1=length(which(list2$P < 5e-8))
# Get number of Bonferroni-corrected significant hits 
hits2=length(which(list2$P < 3.417635e-11))
# Remove file for cleaning purposes
rm(list1)
rm(list2)
gc()
# Store outputs 
# Protein 
mat[ind,1]=as.character(sp[[1]])  
# Lambda 
mat[ind,2]=lambda 
# No. of significant hits at P<5e-8
mat[ind,3]=hits1
# No. of significant hits at P<3.4e-11
mat[ind,4]=hits2
# Uniprot 
mat[ind,5]=as.character(sp[[2]])
# Array
mat[ind,6]=as.character(sp[[5]])
# Filename 
mat[ind,7]=file
# Print to denote completion 
print(ind)
}
# Save out final file 
fwrite(mat,"/rhillary/outputs/lambdas.csv",row.names=F)