res=readRDS("/home/rhillary/outputs/dis_independent_annotated.rds")
# Subset to just those in conditional GEIs
cond=as.data.frame(fread("/home/rhillary/outputs/cond_gei_final2023.csv"))
cond$test=paste(cond$SNP, cond$Protein,sep="_")
res$test=paste(res$SNP, res$compare,sep="_")
res=res[res$test %in% cond$test,]
# Set up results columns 
res$gei_P=NA
res$gei_beta=NA
res$gei_se=NA
res$gei_A1=NA 

# Loop through each GWAS file 
for(i in 1:nrow(res)){
# Extract file of interest 
gw=res[i,"file"]
# Extract SNP 
snp=res[i,"SNP"]
# Check file 
gw=gsub(".phen",".vqtl",gw)
gw1=gsub("chr","",gw)
# Subset to GWAS file of interest 
tmp=as.data.frame(fread(paste0("/scratch/rhillary/residualised_vqtls/",gw)))
# Read in replication version 
tmp1=as.data.frame(fread(paste0("/scratch/rhillary/gei_sensitivity/",gw1)))
# Subset files to SNP being queried 
tmp=tmp[tmp$SNP %in% snp,]
tmp1=tmp1[tmp1$SNP %in% snp,]
if(nrow(tmp1)==nrow(tmp)){ 
# Add in replication stats 
res[i,"gei_P"]=tmp1$P
res[i,"gei_beta"]=tmp1$beta
res[i,"gei_se"]=tmp1$se
res[i,"gei_A1"]=tmp1$A1
# Save out main effect QTLs + print to denote completion
print(i)
} else { 
print("FAILED")
}
}

# Post-formatting 
# Replicated? 
res$Remains=0
res[which(res$gei_P < 0.05),"Remains"]=1
# Sign concordance 
res$dis_sign=ifelse(res$beta<0,"+","-")
res$rep_sign=ifelse(res$gei_beta<0,"+","-")
res1=res[which(res$A1==res$gei_A1),]
res2=res[which(!res$A1==res$gei_A1),]
res3=res[which(is.na(res$gei_A1)),]
#res2$replication_beta=res2$replication_beta*-1
fin=rbind(res1,res2,res3)
fin$dis_sign=ifelse(fin$beta<0,"+","-")
fin$rep_sign=ifelse(fin$gei_beta<0,"+","-")
fin2=fin[,c("SNP","Chr","bp","rsid","Assay","compare","beta","se","P","gei_beta","gei_se","gei_P","Remains")]
write.csv(fin2, "/home/rhillary/outputs/wang_vqtls.csv",row.names=F)
