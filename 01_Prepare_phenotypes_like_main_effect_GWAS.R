########################################################################################
#### Prepare Phenotypes as in  Main Effect GWAS Comparison (main strategy) #############
########################################################################################

# Load requisite libraries
library(data.table)

# Set working directory 
setwd("/UKB_PPP_pQTLs/analysis/")

# Read in discovery phenotype file 
phen1=as.data.frame(fread("discovery_pheno_v1.tsv"))

# Create phenotype files for OSCA 
# Extract protein names 
names.phen <- names(phen1)
# Remove FID and IID from loop - only want to loop through protein names 
names.phen=names.phen[-c(1,2)]
# Start loop 
for(i in names.phen){
# Create temp file for each protein
tmp <- phen1[,c("FID","IID",as.character(i))]
# Convert infinite values to NA
tmp[is.infinite(tmp[,3]),3]<-NA
# Tidy up name
i2=gsub(":", "_", i)
# Check length of file name - needs to be 5 for convenience of OSCA models 
if(length(unlist(strsplit(i2, "_"))) > 5){ 
i3=gsub("^.*?_","",i2)
write.table(tmp, paste0("/rhillary/phenotypes/", i3, ".phen"), quote = F, row.names = F, sep = ' ')
print(i3)
} else 
{ 
# Write out file in format required for OSCA 
write.table(tmp, paste0("/rhillary/phenotypes/", i2, ".phen"), quote = F, row.names = F, sep = ' ')
print(i2)
} 
}


# Replication phenotypes 
# Read in replication phenotype file 
phen1=as.data.frame(fread("replication_pheno_v1.tsv"))

# Create phenotype files for OSCA 
# Extract protein names 
names.phen <- names(phen1)
# Remove FID and IID from loop - only want to loop through protein names 
names.phen=names.phen[-c(1,2)]
# Start loop 
for(i in names.phen){
# Create temp file for each protein
tmp <- phen1[,c("FID","IID",as.character(i))]
# Convert infinite values to NA
tmp[is.infinite(tmp[,3]),3]<-NA
# Tidy up name
i2=gsub(":", "_", i)
# Check length of file name - needs to be 5 
if(length(unlist(strsplit(i2, "_"))) > 5){ 
i3=gsub("^.*?_","",i2)
write.table(tmp, paste0("/rhillary/replication_phenotypes/", i3, ".phen"), quote = F, row.names = F, sep = ' ')
print(i3)
} else 
{ 
## Write out file in format required for OSCA 
write.table(tmp, paste0("/rhillary/replication_phenotypes/", i2, ".phen"), quote = F, row.names = F, sep = ' ')
print(i2)
} 
}
