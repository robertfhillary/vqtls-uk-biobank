##########################################################################################
######## Summary of output from Wang et al. sensitivity analyses  ########################
##########################################################################################

# Read in vQTL results from main analytical strategy 
res=readRDS("/rhillary/outputs/dis_independent_annotated.rds")
# Set up results columns 
res$wang_P=NA
res$wang_beta=NA
res$wang_se=NA
res$wang_A1=NA 

# Loop through each GWAS file 
for(i in 1:nrow(res)){
# Extract file of interest 
gw=res[i,"file"]
# Extract SNP 
snp=res[i,"SNP"]
# Check file 
gw=gsub(".phen",".vqtl",gw)
# Subset to GWAS file of interest 
tmp=as.data.frame(fread(paste0("/rhillary/residualised_vqtls/",gw)))
# Read in replication version 
tmp1=as.data.frame(fread(paste0("/rhillary/wang_vqtls/",gw)))
# Subset files to SNP being queried 
tmp=tmp[tmp$SNP %in% snp,]
tmp1=tmp1[tmp1$SNP %in% snp,]
if(nrow(tmp1)==nrow(tmp)){ 
# Add in replication stats 
res[i,"wang_P"]=tmp1$P
res[i,"wang_beta"]=tmp1$beta
res[i,"wang_se"]=tmp1$se
res[i,"wang_A1"]=tmp1$A1
# Save out main effect QTLs + print to denote completion
print(i)
} else { 
print("FAILED")
}
}

# Post-formatting 
# Replicated? 
res$Remains=0
res[which(res$wang_P < 0.05),"Remains"]=1
# Sign concordance 
res$dis_sign=ifelse(res$beta<0,"+","-")
res$rep_sign=ifelse(res$wang_beta<0,"+","-")
res1=res[which(res$A1==res$wang_A1),]
res2=res[which(!res$A1==res$wang_A1),]
res3=res[which(is.na(res$wang_A1)),]
#res2$replication_beta=res2$replication_beta*-1
fin=rbind(res1,res2,res3)
fin$dis_sign=ifelse(fin$beta<0,"+","-")
fin$rep_sign=ifelse(fin$wang_beta<0,"+","-")
fin2=fin[,c("SNP","Chr","bp","rsid","Assay","compare","beta","se","P","wang_beta","wang_se","wang_P","Remains")]

