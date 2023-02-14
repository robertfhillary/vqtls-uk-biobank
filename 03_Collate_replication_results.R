#####################################################
########### Collate replication results #############
#####################################################

# Load requisite libraries 
library(data.table)
library(ggplot2)

res=readRDS("/outputs/dis_independent_annotated.rds")
# Set up results columns 
res$replication_P=NA
res$replication_beta=NA
res$replication_se=NA
res$replication_A1=NA 

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
tmp1=as.data.frame(fread(paste0("/rhillary/replication_vqtls/",gw)))
# Subset files to SNP being queried 
tmp=tmp[tmp$SNP %in% snp,]
tmp1=tmp1[tmp1$SNP %in% snp,]
if(nrow(tmp1)==nrow(tmp)){ 
# Add in replication stats 
res[i,"replication_P"]=tmp1$P
res[i,"replication_beta"]=tmp1$beta
res[i,"replication_se"]=tmp1$se
res[i,"replication_A1"]=tmp1$A1
# Save out main effect QTLs + print to denote completion
print(i)
} else { 
print("FAILED")
}
}

# Post-formatting 
# Replicated? 
res$Replicated=0
res[which(res$replication_P < 0.05),"Replicated"]=1
# Sign concordance 
res$dis_sign=ifelse(res$beta<0,"+","-")
res$rep_sign=ifelse(res$replication_beta<0,"+","-")
res1=res[which(res$A1==res$replication_A1),]
res2=res[which(!res$A1==res$replication_A1),]
res3=res[which(is.na(res$replication_A1)),]
res2$replication_beta=res2$replication_beta*-1
fin=rbind(res1,res2,res3)
fin$dis_sign=ifelse(fin$beta<0,"+","-")
fin$rep_sign=ifelse(fin$replication_beta<0,"+","-")
res4=res[which(is.na(res$replication_A1)),]

data=fin3
data$xmin=data$beta-(data$se*1.96)
data$xmax=data$beta+(data$se*1.96)
data$ymin=data$replication_beta-(data$replication_se*1.96)
data$ymax=data$replication_beta+(data$replication_se*1.96)

pdf("/rhillary/dissemination2/replication.pdf",width=8,height=7)
print(ggplot(data =data,aes(x = beta,y =replication_beta)) + 
    geom_point(color='darkturquoise') + 
    geom_errorbar(aes(ymin = ymin,ymax = ymax),color='darkturquoise') + 
    geom_errorbarh(aes(xmin = xmin,xmax = xmax),color='darkturquoise') + xlab("Discovery - Betas [95% CI]") + ylab("Replication - Betas [95% CI]") +
     geom_vline(xintercept = 0,linetype="dotted") + geom_hline(yintercept=0,linetype="dotted") + scale_x_continuous(limits=c(-0.7,0.7),breaks=c(-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)) + scale_y_continuous(limits=c(-0.7,0.7),breaks=c(-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)) + theme(axis.text = element_text(size=13.5)) + theme_scientific() +  geom_smooth(method='lm', formula= y~x,se=F, size=0.5,color="#E69F00") + annotate(x = 0.42, y = 0.18, label = expression(paste(italic("r"), " = 0.96")), size = 4.7, geom = "text"))
 dev.off()
