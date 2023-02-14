###################################################
########### Collate discovery results #############
###################################################

# Load requisite libraries 
library(data.table)

# Set working directory 
setwd("/residualised_vqtls/")
# Extract results files 
loop=list.files(".",".vqtl")
# Set up list to store significant results 
res=list()
# Loop step 
for(i in loop){ 
# Extract file name 
file=gsub(".vqtl","",i)
sp=strsplit(file,split="_")[[1]]
# Read in file 
tmp=as.data.frame(fread(i))
# Store filename in results - ease for later
tmp$Protein=as.character(sp[[2]])
tmp$UniProt=as.character(sp[[3]])
tmp$Array=as.character(sp[[6]])
tmp$file=paste0(file,".phen")
# Extract significant hits 
tmp1=tmp[tmp$P<5e-8,]
# If there are no significant hits, then skip file
if(nrow(tmp1)==0){ 
NULL} else { 
# If there are significant results, then store them 
res[[i]] <- tmp1
# Find index of the file in the overall loop and print to denote completion 
ind=which(loop %in% i)
print(ind)
}
}

# Set up files for replication testing - unique CHRs and Proteins 
res1=as.data.frame(do.call("rbind",res))
saveRDS(res1, "/rhillary/discovery_p5e8.rds")
# Subset to Bonferroni correction 
res1=res1[res1$P < 3.417635e-11,]
saveRDS(res1, "/rhillary/discovery_p3.4e11.rds")
# Summarise outputs  
length(unique(res1$UniProt))
# Tidy up file
res1$file=sub(".*?_", "", res1$file)
# Get unique value per model 
res1$model=paste(res1$Chr,res1$Protein,sep="_") 
res2=res1[-which(duplicated(res1$model)),]
# Extract chromosomes and phenotype files to run replication analyses on 
res2=res2[,c("Chr","file")]
res2$Chr=paste0("chr",res2$Chr)
chrs=as.data.frame(res2[,c("Chr")])
names(chrs)[1]="chr"
phenos=as.data.frame(res2[,c("file")])
names(phenos)[1]="phen"
phenos$phen=paste0("/rhillary/replication_phenotypes_patched/", phenos$phen)
# Write out tables 
write.table(chrs,"/rhillary/replication_chrs.txt", row.names=F, col.names=F, quote = F)
write.table(phenos,"/rhillary/wang_phenos.txt", row.names=F, col.names=F, quote = F)
