#############################################
######## GEI DISCORDANCE ANALYSIS ###########
#############################################

# Load requisite libraries
library(data.table)

# Create functions
# Create function for outlier removal - 5SD
outlierID <- function(x, cut=5) {
  xx <- scale(x)
  retval <- ifelse (abs(xx) >= cut, TRUE, FALSE)
  retval
}
outlierTrim <- function(x, cut=5) {
  id <- outlierID(x, cut=cut)
  retval <- ifelse(id, NA, x)
}

# Rank-based inverse normal transformation 
INT <- function(x) qnorm((rank(x, na.last="keep") - 0.5) / sum(!is.na(x)))

########################################
######## NON-MAIN EFFECTS ##############
########################################

# Read in vQTLs and define those that lack main effects 
vqtl=as.data.frame(fread("/rhillary/outputs/dis_independent_annotated_main_effect.csv"))
vqtl=vqtl[which(!vqtl$main_effect_P <3.4e-11),]

# Read in conditional GEIs 
cond=as.data.frame(fread("/rhillary/outputs/cond_gei_final2023.csv"))

# Align files 
vqtl$test=paste(vqtl$SNP,vqtl$compare,sep="_")
cond$test=paste(cond$SNP,cond$Protein,sep="_")
vqtl1=vqtl[which(vqtl$test %in% cond$test),]
cond1=cond[which(cond$test %in% vqtl$test),]

# Set up list of files to extract 
setwd("/rhillary/interaction_vqtls/")
loop=list.files(".", ".")

# Read in phenotypes 
new.phenos=as.data.frame(fread("/rhillary/interaction_phenos_new.txt"))
phenos=new.phenos[,c(1,which(names(new.phenos) %in% cond1$Phenotype))]

# Set up output matrix 
output=as.data.frame(matrix(nrow=(nrow(cond1)*3),ncol=11))
names(output)=c("SNP","rsid","Protein","Analyte","Phenotype","Phenotype_Name","Group","Tertile","Beta","SE","P")

# Loop through each unique vqtl association 
for(i in 1:nrow(cond1)){ 
# Select first set of files 
cond.tmp=cond1[i,]
# Read in snp 
tmp=vqtl1[which(vqtl1$test%in%cond.tmp$test),]
# Extract and read in associated protein 
file=as.data.frame(fread(loop[grep(tmp$test, loop)]))
# Loop through phenotypes 
phen.tmp=phenos[,c(1,which(names(phenos)%in%cond.tmp$Phenotype))]
file1=merge(file,phen.tmp,by.x="FID",by.y="IND")
# Remove outliers 
file1[,4]=outlierTrim(file1[,4])
# Obtain tertiles of phenotypes 
vTert = quantile(file1[,4], c(0:3/3),na.rm=T)
# Classify values
file1$tert = with(file1, 
               cut(file1[,4], 
                   vTert, 
                   include.lowest = T, 
                   labels = c("Low", "Medium", "High")))
# Split into each tertile 
tmp1=file1[file1$tert%in%"Low",]
tmp2=file1[file1$tert%in%"Medium",]
tmp3=file1[file1$tert%in%"High",]
# Extract coefficients 
output[(((i-1)*3)+1),c(9,10,11)]=summary(lm(tmp1[,3] ~ tmp1[,2]))$coefficients[2,c(1,2,4)]
output[(((i-1)*3)+2),c(9,10,11)]=summary(lm(tmp2[,3] ~ tmp2[,2]))$coefficients[2,c(1,2,4)]
output[(((i-1)*3)+3),c(9,10,11)]=summary(lm(tmp3[,3] ~ tmp3[,2]))$coefficients[2,c(1,2,4)]
# Extract remaining column names
# Tertile  
output[(((i-1)*3)+1),8]=unique(tmp1$tert)
output[(((i-1)*3)+2),8]=unique(tmp2$tert)
output[(((i-1)*3)+3),8]=unique(tmp3$tert)
# Phenotype ID, Name and Group 
output[(((i-1)*3)+1),c(5,6,7)]=c(cond.tmp$Phenotype,cond.tmp$Phenotype_Name,cond.tmp$Group)
output[(((i-1)*3)+2),c(5,6,7)]=c(cond.tmp$Phenotype,cond.tmp$Phenotype_Name,cond.tmp$Group)
output[(((i-1)*3)+3),c(5,6,7)]=c(cond.tmp$Phenotype,cond.tmp$Phenotype_Name,cond.tmp$Group)
# SNP, rsid, Protein and Analyte 
output[(((i-1)*3)+1),c(1,2,3,4)]=c(tmp$SNP,tmp$rsid,tmp$compare,tmp$Assay)
output[(((i-1)*3)+2),c(1,2,3,4)]=c(tmp$SNP,tmp$rsid,tmp$compare,tmp$Assay)
output[(((i-1)*3)+3),c(1,2,3,4)]=c(tmp$SNP,tmp$rsid,tmp$compare,tmp$Assay)
# Print to denote completion
print(i)
} 

# Code whether effects are concordant or discordant within each sign 
output$Sign=ifelse(output$Beta >0, "+", "-")
output$Test=paste(output$Protein, output$Phenotype, output$SNP,sep= "_")
# Define results 
res=as.data.frame(table(output$Sign, output$Test) )
res=res[which(!res$Freq ==0),]
res=res[which(!res$Freq ==3),]
mat1=output[which(output$Test %in% unique(res$Var2)),]
mat2=output[which(!output$Test %in% mat1$Test),]
mat1$Sign="Discordant"
mat2$Sign="Concordant"
total=rbind(mat1,mat2)
# Write out file 
fwrite(total,"/rhillary/outputs/nonmaineffect_discordance.csv",row.names=F)


########################################
######## MAIN EFFECTS ##################
########################################

# Read in vQTLs and define those that lack main effects 
vqtl=as.data.frame(fread("/rhillary/outputs/dis_independent_annotated_main_effect.csv"))
vqtl=vqtl[which(vqtl$main_effect_P <3.4e-11),]

# Read in conditional GEIs 
cond=as.data.frame(fread("/rhillary/outputs/cond_gei_final2023.csv"))

# Align files 
vqtl$test=paste(vqtl$SNP,vqtl$compare,sep="_")
cond$test=paste(cond$SNP,cond$Protein,sep="_")
vqtl1=vqtl[which(vqtl$test %in% cond$test),]
cond1=cond[which(cond$test %in% vqtl$test),]

# Set up list of files to extract 
setwd("/rhillary/interaction_vqtls/")
loop=list.files(".", ".")


# Read in phenotypes 
new.phenos=as.data.frame(fread("/rhillary/interaction_phenos_new.txt"))
phenos=new.phenos[,c(1,which(names(new.phenos) %in% cond1$Phenotype))]
phenos$Phenotype

# Remove non-continuous phenotypes 
cond2=cond1[which(!cond1$Phenotype%in%c("Alcohol Intake Frequency","Smoking Status","Season of Blood Draw")),]
# Set up output matrix 
output=as.data.frame(matrix(nrow=nrow(cond2),ncol=11))
names(output)=c("SNP","rsid","Protein","Analyte","Phenotype","Phenotype_Name","Group","Tertile","Beta","SE","P")

# Loop through each unique vqtl association 
for(i in 1:nrow(cond2)){ 
# Select first set of files 
cond.tmp=cond2[i,]
# Read in snp 
tmp=vqtl1[which(vqtl1$test%in%cond.tmp$test),]
# Extract and read in associated protein 
file=as.data.frame(fread(loop[grep(tmp$test, loop)]))
# Loop through phenotypes 
phen.tmp=phenos[,c(1,which(names(phenos)%in%cond.tmp$Phenotype))]
file1=merge(file,phen.tmp,by.x="FID",by.y="IND")
# Remove outliers 
file1[,4]=outlierTrim(file1[,4])
# Obtain tertiles of phenotypes 
vTert = quantile(file1[,4], c(0:3/3),na.rm=T)
# Classify values
file1$tert = with(file1, 
               cut(file1[,4], 
                   vTert, 
                   include.lowest = T, 
                   labels = c("Low", "Medium", "High")))
# Split into each tertile 
tmp1=file1[file1$tert%in%"Low",]
tmp2=file1[file1$tert%in%"Medium",]
tmp3=file1[file1$tert%in%"High",]
# Extract coefficients 
output[(((i-1)*3)+1),c(9,10,11)]=summary(lm(tmp1[,3] ~ tmp1[,2]))$coefficients[2,c(1,2,4)]
output[(((i-1)*3)+2),c(9,10,11)]=summary(lm(tmp2[,3] ~ tmp2[,2]))$coefficients[2,c(1,2,4)]
output[(((i-1)*3)+3),c(9,10,11)]=summary(lm(tmp3[,3] ~ tmp3[,2]))$coefficients[2,c(1,2,4)]
# Extract remaining column names
# Tertile  
output[(((i-1)*3)+1),8]=unique(tmp1$tert)
output[(((i-1)*3)+2),8]=unique(tmp2$tert)
output[(((i-1)*3)+3),8]=unique(tmp3$tert)
# Phenotype ID, Name and Group 
output[(((i-1)*3)+1),c(5,6,7)]=c(cond.tmp$Phenotype,cond.tmp$Phenotype_Name,cond.tmp$Group)
output[(((i-1)*3)+2),c(5,6,7)]=c(cond.tmp$Phenotype,cond.tmp$Phenotype_Name,cond.tmp$Group)
output[(((i-1)*3)+3),c(5,6,7)]=c(cond.tmp$Phenotype,cond.tmp$Phenotype_Name,cond.tmp$Group)
# SNP, rsid, Protein and Analyte 
output[(((i-1)*3)+1),c(1,2,3,4)]=c(tmp$SNP,tmp$rsid,tmp$compare,tmp$Assay)
output[(((i-1)*3)+2),c(1,2,3,4)]=c(tmp$SNP,tmp$rsid,tmp$compare,tmp$Assay)
output[(((i-1)*3)+3),c(1,2,3,4)]=c(tmp$SNP,tmp$rsid,tmp$compare,tmp$Assay)
# Print to denote completion
print(i)
} 

# Code whether effects are concordant or discordant within each sign 
output$Sign=ifelse(output$Beta >0, "+", "-")
output$Test=paste(output$Protein, output$Phenotype, output$SNP,sep= "_")
# Define results 
res=as.data.frame(table(output$Sign, output$Test) )
res=res[which(!res$Freq ==0),]
res=res[which(!res$Freq ==3),]
mat1=output[which(output$Test %in% unique(res$Var2)),]
mat2=output[which(!output$Test %in% mat1$Test),]
mat1$Sign="Discordant"
mat2$Sign="Concordant"
total1=rbind(mat1,mat2)

# Smoking status 
cond2=cond1[which(cond1$Phenotype%in%c("Smoking Status")),]
# Set up output matrix 
output=as.data.frame(matrix(nrow=(nrow(cond2)*3),ncol=11))
names(output)=c("SNP","rsid","Protein","Analyte","Phenotype","Phenotype_Name","Group","Tertile","Beta","SE","P")

# Loop through each unique vqtl association 
for(i in 1:nrow(cond2)){ 
# Select first set of files 
cond.tmp=cond2[i,]
# Read in snp 
tmp=vqtl1[which(vqtl1$test%in%cond.tmp$test),]
# Extract and read in associated protein 
file=as.data.frame(fread(loop[grep(tmp$test, loop)]))
# Loop through phenotypes 
phen.tmp=phenos[,c(1,which(names(phenos)%in%cond.tmp$Phenotype))]
file1=merge(file,phen.tmp,by.x="FID",by.y="IND")
# Remove outliers 
# Split into each tertile 
tmp1=file1[file1[,4]%in%0,]
tmp2=file1[file1[,4]%in%1,]
tmp3=file1[file1[,4]%in%2,]
# Extract coefficients 
output[(((i-1)*3)+1),c(9,10,11)]=summary(lm(tmp1[,3] ~ tmp1[,2]))$coefficients[2,c(1,2,4)]
output[(((i-1)*3)+2),c(9,10,11)]=summary(lm(tmp2[,3] ~ tmp2[,2]))$coefficients[2,c(1,2,4)]
output[(((i-1)*3)+3),c(9,10,11)]=summary(lm(tmp3[,3] ~ tmp3[,2]))$coefficients[2,c(1,2,4)]
# Extract remaining column names
# Tertile  
output[(((i-1)*3)+1),8]="Never"
output[(((i-1)*3)+2),8]="Former"
output[(((i-1)*3)+3),8]="Current"
# Phenotype ID, Name and Group 
output[(((i-1)*3)+1),c(5,6,7)]=c(cond.tmp$Phenotype,cond.tmp$Phenotype_Name,cond.tmp$Group)
output[(((i-1)*3)+2),c(5,6,7)]=c(cond.tmp$Phenotype,cond.tmp$Phenotype_Name,cond.tmp$Group)
output[(((i-1)*3)+3),c(5,6,7)]=c(cond.tmp$Phenotype,cond.tmp$Phenotype_Name,cond.tmp$Group)
# SNP, rsid, Protein and Analyte 
output[(((i-1)*3)+1),c(1,2,3,4)]=c(tmp$SNP,tmp$rsid,tmp$compare,tmp$Assay)
output[(((i-1)*3)+2),c(1,2,3,4)]=c(tmp$SNP,tmp$rsid,tmp$compare,tmp$Assay)
output[(((i-1)*3)+3),c(1,2,3,4)]=c(tmp$SNP,tmp$rsid,tmp$compare,tmp$Assay)
} 
# Code whether effects are concordant or discordant within each sign 
output$Sign=ifelse(output$Beta >0, "+", "-")
output$Test=paste(output$Protein, output$Phenotype, output$SNP,sep= "_")
# Define results 
res=as.data.frame(table(output$Sign, output$Test) )
res=res[which(!res$Freq ==0),]
res=res[which(!res$Freq ==3),]
mat1=output[which(output$Test %in% unique(res$Var2)),]
mat2=output[which(!output$Test %in% mat1$Test),]
mat1$Sign="Discordant"
mat2$Sign="Concordant"
total2=rbind(mat1,mat2)

# Alcohol Intake Frequency 
cond2=cond1[which(cond1$Phenotype%in%c("Alcohol Intake Frequency")),]
# Set up output matrix 
output=as.data.frame(matrix(nrow=(nrow(cond2)*6),ncol=11))
names(output)=c("SNP","rsid","Protein","Analyte","Phenotype","Phenotype_Name","Group","Tertile","Beta","SE","P")

# Loop through each unique vqtl association 
for(i in 1:nrow(cond2)){ 
# Select first set of files 
cond.tmp=cond2[i,]
# Read in snp 
tmp=vqtl1[which(vqtl1$test%in%cond.tmp$test),]
# Extract and read in associated protein 
file=as.data.frame(fread(loop[grep(tmp$test, loop)]))
# Loop through phenotypes 
phen.tmp=phenos[,c(1,which(names(phenos)%in%cond.tmp$Phenotype))]
file1=merge(file,phen.tmp,by.x="FID",by.y="IND")
# Remove outliers 
# Split into each group 
tmp1=file1[file1[,4]%in%1,]
tmp2=file1[file1[,4]%in%2,]
tmp3=file1[file1[,4]%in%3,]
tmp4=file1[file1[,4]%in%4,]
tmp5=file1[file1[,4]%in%5,]
tmp6=file1[file1[,4]%in%6,]
# Extract coefficients 
output[(((i-1)*3)+1),c(9,10,11)]=summary(lm(tmp1[,3] ~ tmp1[,2]))$coefficients[2,c(1,2,4)]
output[(((i-1)*3)+2),c(9,10,11)]=summary(lm(tmp2[,3] ~ tmp2[,2]))$coefficients[2,c(1,2,4)]
output[(((i-1)*3)+3),c(9,10,11)]=summary(lm(tmp3[,3] ~ tmp3[,2]))$coefficients[2,c(1,2,4)]
output[(((i-1)*3)+4),c(9,10,11)]=summary(lm(tmp4[,3] ~ tmp4[,2]))$coefficients[2,c(1,2,4)]
output[(((i-1)*3)+5),c(9,10,11)]=summary(lm(tmp5[,3] ~ tmp5[,2]))$coefficients[2,c(1,2,4)]
output[(((i-1)*3)+6),c(9,10,11)]=summary(lm(tmp6[,3] ~ tmp6[,2]))$coefficients[2,c(1,2,4)]
# Extract remaining column names
# Tertile  
output[(((i-1)*3)+1),8]="Daily"
output[(((i-1)*3)+2),8]="Three-Four-Times-Week"
output[(((i-1)*3)+3),8]="Once-Twice-Week"
output[(((i-1)*3)+4),8]="One-Three-Times-Month"
output[(((i-1)*3)+5),8]="Special-Occasions-Only"
output[(((i-1)*3)+6),8]="Never"
# Phenotype ID, Name and Group 
output[(((i-1)*3)+1),c(5,6,7)]=c(cond.tmp$Phenotype,cond.tmp$Phenotype_Name,cond.tmp$Group)
output[(((i-1)*3)+2),c(5,6,7)]=c(cond.tmp$Phenotype,cond.tmp$Phenotype_Name,cond.tmp$Group)
output[(((i-1)*3)+3),c(5,6,7)]=c(cond.tmp$Phenotype,cond.tmp$Phenotype_Name,cond.tmp$Group)
output[(((i-1)*3)+4),c(5,6,7)]=c(cond.tmp$Phenotype,cond.tmp$Phenotype_Name,cond.tmp$Group)
output[(((i-1)*3)+5),c(5,6,7)]=c(cond.tmp$Phenotype,cond.tmp$Phenotype_Name,cond.tmp$Group)
output[(((i-1)*3)+6),c(5,6,7)]=c(cond.tmp$Phenotype,cond.tmp$Phenotype_Name,cond.tmp$Group)
# SNP, rsid, Protein and Analyte 
output[(((i-1)*3)+1),c(1,2,3,4)]=c(tmp$SNP,tmp$rsid,tmp$compare,tmp$Assay)
output[(((i-1)*3)+2),c(1,2,3,4)]=c(tmp$SNP,tmp$rsid,tmp$compare,tmp$Assay)
output[(((i-1)*3)+3),c(1,2,3,4)]=c(tmp$SNP,tmp$rsid,tmp$compare,tmp$Assay)
output[(((i-1)*3)+4),c(1,2,3,4)]=c(tmp$SNP,tmp$rsid,tmp$compare,tmp$Assay)
output[(((i-1)*3)+5),c(1,2,3,4)]=c(tmp$SNP,tmp$rsid,tmp$compare,tmp$Assay)
output[(((i-1)*3)+6),c(1,2,3,4)]=c(tmp$SNP,tmp$rsid,tmp$compare,tmp$Assay)
} 
# Code whether effects are concordant or discordant within each sign 
output$Sign=ifelse(output$Beta >0, "+", "-")
output$Test=paste(output$Protein, output$Phenotype, output$SNP,sep= "_")
# Define results 
res=as.data.frame(table(output$Sign, output$Test) )
res=res[which(!res$Freq ==0),]
res=res[which(!res$Freq ==6),]
mat1=output[which(output$Test %in% unique(res$Var2)),]
mat2=output[which(!output$Test %in% mat1$Test),]
mat1$Sign="Discordant"
mat2$Sign="Concordant"
total3=rbind(mat1,mat2)

# Season of Blood Draw
cond2=cond1[which(cond1$Phenotype%in%c("Season of Blood Draw")),]
# Set up output matrix 
output=as.data.frame(matrix(nrow=(nrow(cond2)*2),ncol=11))
names(output)=c("SNP","rsid","Protein","Analyte","Phenotype","Phenotype_Name","Group","Tertile","Beta","SE","P")

# Loop through each unique vqtl association 
for(i in 1:nrow(cond2)){ 
# Select first set of files 
cond.tmp=cond2[i,]
# Read in snp 
tmp=vqtl1[which(vqtl1$test%in%cond.tmp$test),]
# Extract and read in associated protein 
file=as.data.frame(fread(loop[grep(tmp$test, loop)]))
# Loop through phenotypes 
phen.tmp=phenos[,c(1,which(names(phenos)%in%cond.tmp$Phenotype))]
file1=merge(file,phen.tmp,by.x="FID",by.y="IND")
# Remove outliers 
# Split into each tertile 
tmp1=file1[file1[,4]%in%"Winter/Spring",]
tmp2=file1[file1[,4]%in%"Summer/Autumn",]
# Extract coefficients 
output[(((i-1)*3)+1),c(9,10,11)]=summary(lm(tmp1[,3] ~ tmp1[,2]))$coefficients[2,c(1,2,4)]
output[(((i-1)*3)+2),c(9,10,11)]=summary(lm(tmp2[,3] ~ tmp2[,2]))$coefficients[2,c(1,2,4)]
# Extract remaining column names
# Tertile  
output[(((i-1)*3)+1),8]="Winter/Spring"
output[(((i-1)*3)+2),8]="Summer/Autumn"
# Phenotype ID, Name and Group 
output[(((i-1)*3)+1),c(5,6,7)]=c(cond.tmp$Phenotype,cond.tmp$Phenotype_Name,cond.tmp$Group)
output[(((i-1)*3)+2),c(5,6,7)]=c(cond.tmp$Phenotype,cond.tmp$Phenotype_Name,cond.tmp$Group)
# SNP, rsid, Protein and Analyte 
output[(((i-1)*3)+1),c(1,2,3,4)]=c(tmp$SNP,tmp$rsid,tmp$compare,tmp$Assay)
output[(((i-1)*3)+2),c(1,2,3,4)]=c(tmp$SNP,tmp$rsid,tmp$compare,tmp$Assay)
} 
# Code whether effects are concordant or discordant within each sign 
output$Sign=ifelse(output$Beta >0, "+", "-")
output$Test=paste(output$Protein, output$Phenotype, output$SNP,sep= "_")
# Define results 
res=as.data.frame(table(output$Sign, output$Test) )
res=res[which(!res$Freq ==0),]
res=res[which(!res$Freq ==2),]
mat1=output[which(output$Test %in% unique(res$Var2)),]
mat2=output[which(!output$Test %in% mat1$Test),]
mat1$Sign="Discordant"
mat2$Sign="Concordant"
total4=rbind(mat1,mat2)

# Combine outputs 
total5=rbind(total1,total2,total3,total4)
fwrite(total5,"/rhillary/outputs/maineffect_discordance.csv",row.names=F)

# Summary 
final=rbind(total,total5)
fwrite(final,"/rhillary/outputs/discordance_all_qtls.csv",row.names=F)
final1=final[which(!duplicated(final$Test)),]
# Ask if main effect qtl 
final1$Main_Effect=0
final1$test=paste(final1$SNP,final1$Protein,sep="_")
final1[final1$test%in%vqtl$test,"Main_Effect"]=1
final1=final1[,c("SNP","rsid","Protein","Analyte","Sign","Main_Effect","test")]
fwrite(final1,"/rhillary/outputs/discordance_summary.csv",row.names=F)
# Info for paper 
non=final1[final1$Main_Effect==0,]
on=final1[final1$Main_Effect==1,]


