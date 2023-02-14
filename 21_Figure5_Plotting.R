##########################################
######## FIGURE 5 ########################
##########################################

# Load requisite libraries 
library(data.table)
library(ggplot2)
library(artyfarty)
library(forcats)
library(RColorBrewer)
library(grafify)
library(cowplot)

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

#################################
##### FIGURE 6A #################
################################# 

# Extract GEI for PAEP and combine with phenotype data 
setwd("/rhillary/interaction_vqtls/")
loop=list.files(".", ".")
new.phenos=as.data.frame(fread("/rhillary/all_quant_20200904.pheno.v2.tsv.gz"))
names(new.phenos)[1]="FID"
tmp=as.data.frame(fread(paste0("/rhillary/interaction_vqtls/",loop[633])))
vqtl=readRDS("/home/rhillary/outputs/dis_independent_annotated.rds")
tmp1=merge(tmp,new.phenos,by="FID")
# Define tertiles 
vTert = quantile(tmp1$f_21002_0_0, c(0:3/3),na.rm=T)
# Classify values
tmp1$tert = with(tmp1, 
               cut(f_21002_0_0, 
                   vTert, 
                   include.lowest = T, 
                   labels = c("Low", "Medium", "High")))
# Prepare extra info for plots 
dat=as.data.frame(fread("/rhillary/dat_wide_internal_use.csv"))
dat1=dat[,c(1,grep("PAEP",names(dat)))]
names(dat1)=c("FID","PAEP.2")
tmp1=merge(tmp1,dat1,by="FID")
names(tmp1)[3]="PAEP"
tmp1=tmp1[which(!is.na(tmp1$tert)),]
tmp1=tmp1[which(!is.na(tmp1[,2])),]
tmp1=tmp1[which(!is.na(tmp1[,3])),]

# Get summary data 
tmp1$supergroup=paste(tmp1[,2],tmp1$tert,sep="_")
data=as.data.frame(tapply(tmp1[,3], tmp1$supergroup, mean))
data1=as.data.frame(tapply(tmp1[,3], tmp1$supergroup, function(x)sd(x)/sqrt(length(x))))
data2=cbind(data,data1)
names(data2)=c("mean","se")
data2$HCI=data2$mean+(data2$se*1.96)
data2$LCI=data2$mean-(data2$se*1.96)
data2$group=row.names(data2)
gsub("_.*", "", data2$group)
data2$SNP=gsub("_.*", "", data2$group)
data2$tert=gsub("*._", "", data2$group)

# Get info for plot 
a=table(tmp1[,2])[1]
b=table(tmp1[,2])[2]
c=table(tmp1[,2])[3]
a1=vqtl[which(vqtl$SNP %in% names(tmp1)[2]),"A1"]
a2=vqtl[which(vqtl$SNP %in% names(tmp1)[2]),"A2"]
rsid=vqtl[which(vqtl$SNP %in% names(tmp1)[2]), "rsid"]
prot=vqtl[which(vqtl$SNP %in% names(tmp1)[2]),"Assay"]
g1=range(tmp1[which(tmp1$tert%in%"Low"),"f_21002_0_0"],na.rm=T)
g2=range(tmp1[which(tmp1$tert%in%"Medium"),"f_21002_0_0"],na.rm=T)
g3=range(tmp1[which(tmp1$tert%in%"High"),"f_21002_0_0"],na.rm=T)

# Get ranges of tertiles 
g1=paste0("<",g1[2])
g3=paste0(">", g3[1])
g2=paste0(g2[1], "-", g2[2])

# Prepare plot 
data2$tert=factor(data2$tert, levels=c("Low", "Medium", "High"))
# Plot stage 
pdf("/rhillary/dissemination2/interaction_paep_transform2.pdf",height=7,width=10)
ggplot(data2, aes(x=factor(SNP),y=mean)) +
  geom_point(
    aes(color = tert),
    position = position_dodge(0.3), size =5) +   geom_errorbar(
    aes(ymin = LCI, ymax = HCI, color = tert),
    position = position_dodge(0.3), width = 0.4,
    ) +
    scale_x_discrete(name=paste0("Genotype - ",rsid), labels = c(paste0(a2,"/",a2, "\n", "(n=",a,")"), paste0(a2,"/",a1, "\n", "(n=",b,")") , paste0(a1,"/",a1, "\n", "(n=",c,")"))) + 
     scale_color_brewer(palette = "Set2", labels = c(g1,g2,g3), name="Body Weight \n(kg) Tertile") + ylab(paste0("Mean Transformed ", prot, "\n levels [95 %CI]")) + 
        theme_scientific() + scale_y_continuous(limits=c(-1.5,1.5),breaks=c(-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5)) + theme(legend.position="right",axis.title =element_text(size=16), axis.text=element_text(size=15), legend.text=element_text(size=15),legend.title=element_text(size=15)) 
 dev.off()


#################################
##### FIGURE 6B #################
################################# 

# Male and female specific plots 
### FEMALES #####

# Read in files 
demo=as.data.frame(fread("UKB_PPP_pQTLs/analysis/input/baseline_covars.tsv"))
demo=demo[,c("f.eid", "f.31.0.0", "f.21003.0.0")]
names(demo)=c("FID","sex","age")
tmp1=merge(tmp,demo,by="FID")
# Subset to females 
tmp1=tmp1[which(tmp1$sex%in%0),]

# Read in vQTL data 
vqtl=readRDS("/rhillary/outputs/dis_independent_annotated.rds")
# Merge datasets 
tmp1=merge(tmp1,new.phenos,by="FID")

# Define tertiles 
vTert = quantile(tmp1$f_21002_0_0, c(0:3/3),na.rm=T)
# Classify values
tmp1$tert = with(tmp1, 
               cut(f_21002_0_0, 
                   vTert, 
                   include.lowest = T, 
                   labels = c("Low", "Medium", "High")))

# Obtain remaining data needed for plot 
dat=as.data.frame(fread("/rhillary/dat_wide_internal_use.csv"))
dat1=dat[,c(1,grep("PAEP",names(dat)))]
names(dat1)=c("FID","PAEP.2")
tmp1=merge(tmp1,dat1,by="FID")
names(tmp1)[3]="PAEP"
tmp1=tmp1[which(!is.na(tmp1$tert)),]
tmp1=tmp1[which(!is.na(tmp1[,2])),]
tmp1=tmp1[which(!is.na(tmp1[,3])),]

# Get summary data 
tmp1$supergroup=paste(tmp1[,2],tmp1$tert,sep="_")
data=as.data.frame(tapply(tmp1[,3], tmp1$supergroup, mean))
data1=as.data.frame(tapply(tmp1[,3], tmp1$supergroup, function(x)sd(x)/sqrt(length(x))))
data2=cbind(data,data1)
names(data2)=c("mean","se")
data2$HCI=data2$mean+(data2$se*1.96)
data2$LCI=data2$mean-(data2$se*1.96)
data2$group=row.names(data2)
gsub("_.*", "", data2$group)
data2$SNP=gsub("_.*", "", data2$group)
data2$tert=gsub("*._", "", data2$group)

# Get info needed for plot 
a=table(tmp1[,2])[1]
b=table(tmp1[,2])[2]
c=table(tmp1[,2])[3]
a1=vqtl[which(vqtl$SNP %in% names(tmp1)[2]),"A1"]
a2=vqtl[which(vqtl$SNP %in% names(tmp1)[2]),"A2"]
rsid=vqtl[which(vqtl$SNP %in% names(tmp1)[2]), "rsid"]
prot=vqtl[which(vqtl$SNP %in% names(tmp1)[2]),"Assay"]
g1=range(tmp1[which(tmp1$tert%in%"Low"),"f_21002_0_0"],na.rm=T)
g2=range(tmp1[which(tmp1$tert%in%"Medium"),"f_21002_0_0"],na.rm=T)
g3=range(tmp1[which(tmp1$tert%in%"High"),"f_21002_0_0"],na.rm=T)

# Get ranges of tertiles 
g1=paste0("<",g1[2])
g3=paste0(">", g3[1])
g2=paste0(g2[1], "-", g2[2])

# Prepare for plot 
data2$tert=factor(data2$tert, levels=c("Low", "Medium", "High"))
data2$label="Females"
  my_col = brewer.pal(n = 5, "Blues")[3:5]

# Female-specific plot 
pdf("/rhillary/dissemination2/interaction_paep_transform_females.pdf",height=7,width=8)
plotx1=ggplot(data2, aes(x=factor(SNP),y=mean)) +
  geom_point(
    aes(color = tert),
    position = position_dodge(0.3), size =5) +   geom_errorbar(
    aes(ymin = LCI, ymax = HCI, color = tert),
    position = position_dodge(0.3), width = 0.4,
    ) +
    scale_x_discrete(name=paste0("Genotype - ",rsid), labels = c(paste0(a2,"/",a2, "\n", "(n=",a,")"), paste0(a2,"/",a1, "\n", "(n=",b,")") , paste0(a1,"/",a1, "\n", "(n=",c,")"))) + 
     scale_color_manual(values = my_col, labels = c(g1,g2,g3), name="Body Weight \n(kg) Tertile") + ylab(paste0("Mean Transformed ", prot, "\n levels [95% CI]")) +
        theme_scientific() + scale_y_continuous(limits=c(-1.5,1.5),breaks=c(-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5)) + theme(legend.position="right",axis.title =element_text(size=13), axis.text=element_text(size=13), legend.text=element_text(size=13),legend.title=element_text(size=13))+facet_wrap(factor(data2$label))+ theme(strip.background =element_rect(fill="grey94"),strip.text.x = element_text(size = 14.5)) 
 print(plotx1)
 dev.off()


### MALES #####

# Repeat as above 
demo=as.data.frame(fread("/UKB_PPP_pQTLs/analysis/input/baseline_covars.tsv"))
demo=demo[,c("f.eid", "f.31.0.0", "f.21003.0.0")]
names(demo)=c("FID","sex","age")
tmp2=merge(tmp,demo,by="FID")
# Subset to males 
tmp2=tmp2[which(tmp2$sex%in%1),]

# Read in and merge datasets 
vqtl=readRDS("/rhillary/outputs/dis_independent_annotated.rds")
tmp2=merge(tmp2,new.phenos,by="FID")

# Define tertiles 
vTert = quantile(tmp2$f_21002_0_0, c(0:3/3),na.rm=T)
# Classify values
tmp2$tert = with(tmp2, 
               cut(f_21002_0_0, 
                   vTert, 
                   include.lowest = T, 
                   labels = c("Low", "Medium", "High")))

# Prepare for plot 
daty=as.data.frame(fread("/rhillary/dat_wide_internal_use.csv"))
dat1y=daty[,c(1,grep("PAEP",names(daty)))]
names(dat1y)=c("FID","PAEP.2")
tmp2=merge(tmp2,dat1y,by="FID")
names(tmp2)[3]="PAEP"
tmp2=tmp2[which(!is.na(tmp2$tert)),]
tmp2=tmp2[which(!is.na(tmp2[,2])),]
tmp2=tmp2[which(!is.na(tmp2[,3])),]

# Get summary data
tmp2$supergroup=paste(tmp2[,2],tmp2$tert,sep="_")
datay=as.data.frame(tapply(tmp2[,3], tmp2$supergroup, mean))
data1y=as.data.frame(tapply(tmp2[,3], tmp2$supergroup, function(x)sd(x)/sqrt(length(x))))
data2y=cbind(datay,data1y)
names(data2y)=c("mean","se")
data2y$HCI=data2y$mean+(data2y$se*1.96)
data2y$LCI=data2y$mean-(data2y$se*1.96)
data2y$group=row.names(data2y)
gsub("_.*", "", data2y$group)
data2y$SNP=gsub("_.*", "", data2y$group)
data2y$tert=gsub("*._", "", data2y$group)

# Get info needed for plot 
ay=table(tmp2[,2])[1]
by=table(tmp2[,2])[2]
cy=table(tmp2[,2])[3]
a1y=vqtl[which(vqtl$SNP %in% names(tmp2)[2]),"A1"]
a2y=vqtl[which(vqtl$SNP %in% names(tmp2)[2]),"A2"]
rsidy=vqtl[which(vqtl$SNP %in% names(tmp2)[2]), "rsid"]
proty=vqtl[which(vqtl$SNP %in% names(tmp2)[2]),"Assay"]
g1y=range(tmp2[which(tmp2$tert%in%"Low"),"f_21002_0_0"],na.rm=T)
g2y=range(tmp2[which(tmp2$tert%in%"Medium"),"f_21002_0_0"],na.rm=T)
g3y=range(tmp2[which(tmp2$tert%in%"High"),"f_21002_0_0"],na.rm=T)

# Get ranges of tertiles 
g1y=paste0("<",g1y[2])
g3y=paste0(">", g3y[1])
g2y=paste0(g2y[1], "-", g2y[2])

# Prepare for plot 
data2y$tert=factor(data2y$tert, levels=c("Low", "Medium", "High"))
data2y$label="Males"
  my_coly = brewer.pal(n = 5, "Greens")[3:5]

# Male-specific plots 
pdf("/rhillary/dissemination2/interaction_paep_transform_males.pdf",height=7,width=9)
plotx2=ggplot(data2y, aes(x=factor(SNP),y=mean)) +
  geom_point(
    aes(color = tert),
    position = position_dodge(0.3), size =5) +   geom_errorbar(
    aes(ymin = LCI, ymax = HCI, color = tert),
    position = position_dodge(0.3), width = 0.4,
    ) +
    scale_x_discrete(name=paste0("Genotype - ",rsidy), labels = c(paste0(a2y,"/",a2y, "\n", "(n=",ay,")"), paste0(a2y,"/",a1y, "\n", "(n=",by,")") , paste0(a1y,"/",a1y, "\n", "(n=",cy,")"))) + 
     scale_color_manual(values = my_coly, labels = c(g1y,g2y,g3y), name="Body Weight \n(kg) Tertile") + ylab(paste0("Mean Transformed ", proty, "\n levels [95% CI]")) +
        theme_scientific() + scale_y_continuous(limits=c(-1.5,1.5),breaks=c(-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5)) + theme(legend.position="right",axis.title =element_text(size=13), axis.text=element_text(size=13), legend.text=element_text(size=13),legend.title=element_text(size=13))+facet_wrap(factor(data2y$label))+ theme(strip.background =element_rect(fill="grey94"),strip.text.x = element_text(size = 14.5)) 
 print(plotx2)
 dev.off()

# Combine female- and male-specific plots 
pdf("/home/rhillary/dissemination2/interaction_paep_transform_sexes.pdf",height=7,width=13)
plot1=plot_grid(plotx1,plotx2,nrow=1,ncol=2)
print(plot1)
dev.off()

#################################
##### FIGURE 6C #################
################################# 

# Read in phenotype data 
phenos=as.data.frame(fread("/rhillary/all_quant_20200904.pheno.v2.tsv.gz"))
# QC 
names(phenos)[1]="FID"
phenos$f_3089_0_0=NULL
phenos$f_3090_0_0=NULL
phenos$f_30230_0_0=NULL
# Extract results 
loop=list.files("/rhillary/interaction_vqtls_threeterm/",".")
# Extract PAEP result specifically from age-by-sex-genotype analyses 
taskid=grep("PAEP",loop)
tmp=as.data.frame(fread(paste0("/scratch/rhillary/interaction_vqtls_threeterm/",loop[[taskid]])))
# Prepare for further analyses 
test=paste(names(tmp)[3],names(tmp)[2],sep="_")
# Subset phenotype file to only handle phenotypes with conditionally significant interaction
cond=as.data.frame(fread("/rhillary/outputs/cond_gei_final2023.csv"))
# Format further 
cond$Protein <- sub("_[^_]+$", "", cond$Protein)
cond$Protein=gsub("_","\\.",cond$Protein)
cond$test=paste(cond$Protein, cond$SNP, sep="_")
cond1=cond[which(cond$test %in% test),]
phenos1=phenos[,c(1,which(names(phenos) %in% c(cond1$Phenotype,"f_21002_0_0")))]
# Now loop through each phenotype to run interaction test
tmp1=merge(tmp,phenos1,by="FID")

# Now, based on demographics we need to create 12 groups at first
a=as.data.frame(fread("/camhpc/home/bsun/UKB_PPP_pQTLs/analysis/input/baseline_covars.tsv"))
demo=a[,c("f.eid", "f.31.0.0", "f.21003.0.0")]
names(demo)=c("FID","sex","age")
demo1=demo[which(demo$FID %in% phenos1$FID),]

# Groups 1-2
grp=demo1[which(demo1$age <= 50),]
grp1=grp[grp$sex %in% 1,]
grp2=grp[grp$sex %in% 0,]
# Groups 3-4
grp=demo1[which(demo1$age > 50 & demo1$age <= 60),]
grp3=grp[grp$sex %in% 1,]
grp4=grp[grp$sex %in% 0,]
# Groups 5-6
grp=demo1[which(demo1$age > 60 & demo1$age <= 70),]
grp5=grp[grp$sex %in% 1,]
grp6=grp[grp$sex %in% 0,]


# Create new column based on groups 
tmp1$group=NA 
tmp1[which(tmp1$FID %in% grp1$FID),"group"]="Males 40-50yr"
tmp1[which(tmp1$FID %in% grp2$FID),"group"]="Females 40-50yr"
tmp1[which(tmp1$FID %in% grp3$FID),"group"]="Males 50-60yr"
tmp1[which(tmp1$FID %in% grp4$FID),"group"]="Females 50-60yr"
tmp1[which(tmp1$FID %in% grp5$FID),"group"]="Males 60-70yr"
tmp1[which(tmp1$FID %in% grp6$FID),"group"]="Females 60-70yr"


# Set up for plots 
vqtl=readRDS("/rhillary/outputs/dis_independent_annotated.rds")
myplots=list()
for(i in unique(tmp1$group)){ 
tmp2=tmp1[which(as.character(tmp1$group)%in%i),]
# Split into tertiles 
vTert = quantile(tmp2$f_21002_0_0, c(0:3/3),na.rm=T)
# classify values
tmp2$tert = with(tmp2, 
               cut(f_21002_0_0, 
                   vTert, 
                   include.lowest = T, 
                   labels = c("Low", "Medium", "High")))

names(tmp2)[3]="PAEP"
tmp2=tmp2[which(!is.na(tmp2$tert)),]
tmp2=tmp2[which(!is.na(tmp2[,2])),]
tmp2=tmp2[which(!is.na(tmp2[,3])),]

# Get Summary data 
tmp2$supergroup=paste(tmp2[,2,],tmp2$tert,sep="_")
data=as.data.frame(tapply(tmp2[,3], tmp2$supergroup, mean))
data1=as.data.frame(tapply(tmp2[,3], tmp2$supergroup, function(x)sd(x)/sqrt(length(x))))
data2=cbind(data,data1)
names(data2)=c("mean","se")
data2$HCI=data2$mean+(data2$se*1.96)
data2$LCI=data2$mean-(data2$se*1.96)
data2$group=row.names(data2)
gsub("_.*", "", data2$group)
data2$SNP=gsub("_.*", "", data2$group)
data2$tert=gsub("*._", "", data2$group)

# Get info for plots 
a=table(tmp2[,2])[1]
b=table(tmp2[,2])[2]
c=table(tmp2[,2])[3]
a1=vqtl[which(vqtl$SNP %in% names(tmp2)[2]),"A1"]
a2=vqtl[which(vqtl$SNP %in% names(tmp2)[2]),"A2"]
rsid=vqtl[which(vqtl$SNP %in% names(tmp2)[2]),"rsid"]
prot=vqtl[which(vqtl$SNP %in% names(tmp2)[2]),"Assay"]
g1=range(tmp2[which(tmp2$tert%in%"Low"),"f_21002_0_0"],na.rm=T)
g2=range(tmp2[which(tmp2$tert%in%"Medium"),"f_21002_0_0"],na.rm=T)
g3=range(tmp2[which(tmp2$tert%in%"High"),"f_21002_0_0"],na.rm=T)

# Get ranges of tertiles 
g1=paste0("<",g1[2])
g3=paste0(">", g3[1])
g2=paste0(g2[1], "-", g2[2])

# Set up for facet portion 
data2$overall = i
# Set up for plotting 
data2$tert=factor(data2$tert, levels=c("Low", "Medium", "High"))
if(i %in% c("Males 40-50yr","Males 50-60yr","Males 60-70yr")) { 
	my_col = brewer.pal(n = 5, "Greens")[3:5]
p1=ggplot(data2, aes(x=factor(SNP),y=mean)) +
  geom_point(
    aes(color = tert),
    position = position_dodge(0.3), size =3.3) +   geom_errorbar(
    aes(ymin = LCI, ymax = HCI, color = tert),
    position = position_dodge(0.3), width = 0.4,
    ) +
    scale_x_discrete(name=paste0(rsid), labels = c(paste0(a2,"/",a2, "\n", "(n=",a,")"), paste0(a2,"/",a1, "\n", "(n=",b,")") , paste0(a1,"/",a1, "\n", "(n=",c,")"))) + 
     scale_colour_manual(values = my_col, labels = c(g1,g2,g3), name="Body Weight \n (kg) Tertile") + ylab(paste0("Mean Transformed ",prot, "\n levels [95% CI]")) + ylim(limits=c(-2.5,2.5))+ 
   theme_scientific() +  theme(legend.position="right",axis.title =element_text(size=13), axis.text=element_text(size=13), legend.text=element_text(size=11),legend.title=element_text(size=13))+facet_wrap(~overall)+ theme(strip.background =element_rect(fill="grey94"),strip.text.x = element_text(size = 14))

} else { 
	my_col = brewer.pal(n = 5, "Blues")[3:5]
p1=ggplot(data2, aes(x=factor(SNP),y=mean)) +
  geom_point(
    aes(color = tert),
    position = position_dodge(0.3), size =3.3) +   geom_errorbar(
    aes(ymin = LCI, ymax = HCI, color = tert),
    position = position_dodge(0.3), width = 0.4,
    ) +
    scale_x_discrete(name="", labels = c(paste0(a2,"/",a2, "\n", "(n=",a,")"), paste0(a2,"/",a1, "\n", "(n=",b,")") , paste0(a1,"/",a1, "\n", "(n=",c,")"))) + 
     scale_colour_manual(values=my_col, labels = c(g1,g2,g3), name="Body Weight \n (kg) Tertile") + ylab(paste0("Mean Transformed ",prot, "\n levels [95% CI]"))  +
   theme_scientific() +  theme(legend.position="right",axis.title =element_text(size=13), axis.text=element_text(size=13), legend.text=element_text(size=11),legend.title=element_text(size=13)) +facet_wrap(~overall)+ theme(strip.background =element_rect(fill="grey94"),strip.text.x = element_text(size = 14))

}
myplots[[i]]=p1
} 

# Combine plots - just females for Figure 5B
pdf("/home/rhillary/dissemination2/PAEP_agesex_noresid.pdf",width=15,height=6)
p1.1=myplots[[1]]
p2=myplots[[2]]
p3=myplots[[3]]
p4=myplots[[4]]
p5=myplots[[5]]
p6=myplots[[6]]

plot1=plot_grid(p6,p1.1,p2,nrow=1,ncol=3)
print(plot1)
dev.off()


#################################
##### FIGURE S1 #################
################################# 

# Collate age*sex specific conditional GEIs 
setwd("/rhillary/interaction_vqtls_threeway_output2/")
loop=list.files(".", ".")
final=lapply(loop,function(x) fread(x,header=T))
## Combine final files 
final1=as.data.frame(do.call("rbind",final))

# Test PAEP and body weight 
tmp=final1[grep("PAEP", final1$Protein),]
tmp=tmp[tmp$Phenotype %in% "f_21002_0_0",]
tmp[grep("F",tmp$group),"Sex"]="Females" 
tmp[grep("M",tmp$group),"Sex"]="Males"
# Obtain Z scores 
tmp$Interaction_Z=tmp$Interaction_Beta/tmp$Interaction_SE
# Set up plot 
tmp$group=factor(tmp$group,levels=rev(c("F-50","M-50","F-50.60","M-50.60","F-60","M-60")))
# Plot stage 
 inter=ggplot(tmp, aes(x=factor(group),y=Interaction_Beta)) +  geom_point(
    aes(color = Sex),alpha =0.7,
    position = position_dodge(0.3), size =4) +   geom_errorbar(
    aes(ymin = Interaction_Beta-Interaction_SE, ymax = Interaction_Beta+Interaction_SE, color = Sex),
    position = position_dodge(0.3), width = 0.4,
    ) + scale_y_continuous(limits=c(-0.10,0.25)) + scale_color_manual(values = c("Females" = "dodgerblue", "Males" = "forestgreen"))+
 theme_scientific()+scale_x_discrete(name="Demographic Group", labels=rev(c("M 40-50yr","F 40-50yr","M 50-60yr","F 50-60yr","M 60-70yr","F 60-70yr"))) + ylab("Interaction Beta [SE]") + theme(legend.title=element_text(hjust=0.5), axis.text=element_text(size=12), axis.title=element_text(size=13)) + coord_flip()
# Finish plot 
inter1=inter + geom_text(aes(label=ifelse(tmp$Interaction_P < 0.001, "***", ifelse(tmp$Interaction_P<0.05,"*",''))),y=0.22,vjust=0.83,size=7,col="gray15")
# Save out plot 
pdf("/rhillary/dissemination2/PAEP_interaction_coefficients.pdf",width=7,height=6)
print(inter1)
dev.off()


#################################
##### FIGURE S2 #################
################################# 

# Read in datasets
phenos=as.data.frame(fread("/rhillary/all_quant_20200904.pheno.v2.tsv.gz"))
# QC
names(phenos)[1]="FID"
phenos$f_3089_0_0=NULL
phenos$f_3090_0_0=NULL
phenos$f_30230_0_0=NULL
# Extract PAEP
loop=list.files("/rhillary/interaction_vqtls_threeterm/",".")
taskid=grep("PAEP",loop)
tmp=as.data.frame(fread(paste0("/rhillary/interaction_vqtls_threeterm/",loop[[taskid]])))
test=paste(names(tmp)[3],names(tmp)[2],sep="_")
# Tidy up file 
cond=as.data.frame(fread("/rhillary/outputs/cond_gei_final2023.csv"))
cond$Protein <- sub("_[^_]+$", "", cond$Protein)
cond$Protein=gsub("_","\\.",cond$Protein)
cond$test=paste(cond$Protein, cond$SNP, sep="_")
cond1=cond[which(cond$test %in% test),]
phenos1=phenos[,c(1,which(names(phenos) %in% c(cond1$Phenotype,"f_21002_0_0")))]
tmp1=merge(tmp,phenos1,by="FID")

# Now, based on demographics we need to create 12 distinct groups according to sex and 5-year age bins 
a=as.data.frame(fread("/UKB_PPP_pQTLs/analysis/input/baseline_covars.tsv"))
demo=a[,c("f.eid", "f.31.0.0", "f.21003.0.0")]
names(demo)=c("FID","sex","age")
demo1=demo[which(demo$FID %in% phenos1$FID),]

# Groups 1-2
grp=demo1[which(demo1$age <= 45),]
grp1=grp[grp$sex %in% 1,]
grp2=grp[grp$sex %in% 0,]
# Groups 3-4
grp=demo1[which(demo1$age > 45 & demo1$age <= 50),]
grp3=grp[grp$sex %in% 1,]
grp4=grp[grp$sex %in% 0,]
# Groups 5-6
grp=demo1[which(demo1$age > 50 & demo1$age <= 55),]
grp5=grp[grp$sex %in% 1,]
grp6=grp[grp$sex %in% 0,]
# Groups 7-8
grp=demo1[which(demo1$age > 55 & demo1$age <= 60),]
grp7=grp[grp$sex %in% 1,]
grp8=grp[grp$sex %in% 0,]
# Groups 9-10
grp=demo1[which(demo1$age > 60 & demo1$age <= 65),]
grp9=grp[grp$sex %in% 1,]
grp10=grp[grp$sex %in% 0,]
# Groups 11-12
grp=demo1[which(demo1$age > 65 & demo1$age <= 70),]
grp11=grp[grp$sex %in% 1,]
grp12=grp[grp$sex %in% 0,]

# Create new column based on groups 
tmp1$group=NA 
tmp1[which(tmp1$FID %in% grp1$FID),"group"]="Males 40-45yr"
tmp1[which(tmp1$FID %in% grp2$FID),"group"]="Females 40-45yr"
tmp1[which(tmp1$FID %in% grp3$FID),"group"]="Males 45-50yr"
tmp1[which(tmp1$FID %in% grp4$FID),"group"]="Females 45-50yr"
tmp1[which(tmp1$FID %in% grp5$FID),"group"]="Males 50-55yr"
tmp1[which(tmp1$FID %in% grp6$FID),"group"]="Females 50-55yr"
tmp1[which(tmp1$FID %in% grp7$FID),"group"]="Males 55-60yr"
tmp1[which(tmp1$FID %in% grp8$FID),"group"]="Females 55-60yr"
tmp1[which(tmp1$FID %in% grp9$FID),"group"]="Males 60-65yr"
tmp1[which(tmp1$FID %in% grp10$FID),"group"]="Females 60-65yr"
tmp1[which(tmp1$FID %in% grp11$FID),"group"]="Males 65-70yr"
tmp1[which(tmp1$FID %in% grp12$FID),"group"]="Females 65-70yr"


# Read in vQTL data 
vqtl=readRDS("/home/rhillary/outputs/dis_independent_annotated.rds")
# Prepare list to store plots 
myplots=list()
# Loop through each group 
for(i in unique(tmp1$group)){ 
# Subset to this group 
tmp2=tmp1[which(as.character(tmp1$group)%in%i),]
# Define tertiles for this group 
vTert = quantile(tmp2$f_21002_0_0, c(0:3/3),na.rm=T)
# Classify values
tmp2$tert = with(tmp2, 
               cut(f_21002_0_0, 
                   vTert, 
                   include.lowest = T, 
                   labels = c("Low", "Medium", "High")))

# Prepare data for plot 
names(tmp2)[3]="PAEP"
tmp2=tmp2[which(!is.na(tmp2$tert)),]
tmp2=tmp2[which(!is.na(tmp2[,2])),]
tmp2=tmp2[which(!is.na(tmp2[,3])),]

# Get summary data 
tmp2$supergroup=paste(tmp2[,2,],tmp2$tert,sep="_")
data=as.data.frame(tapply(tmp2[,3], tmp2$supergroup, mean))
data1=as.data.frame(tapply(tmp2[,3], tmp2$supergroup, function(x)sd(x)/sqrt(length(x))))
data2=cbind(data,data1)
names(data2)=c("mean","se")
data2$HCI=data2$mean+(data2$se*1.96)
data2$LCI=data2$mean-(data2$se*1.96)
data2$group=row.names(data2)
gsub("_.*", "", data2$group)
data2$SNP=gsub("_.*", "", data2$group)
data2$tert=gsub("*._", "", data2$group)

# Get info needed for plot 
a=table(tmp2[,2])[1]
b=table(tmp2[,2])[2]
c=table(tmp2[,2])[3]
a1=vqtl[which(vqtl$SNP %in% names(tmp2)[2]),"A1"]
a2=vqtl[which(vqtl$SNP %in% names(tmp2)[2]),"A2"]
rsid=vqtl[which(vqtl$SNP %in% names(tmp2)[2]),"rsid"]
prot=vqtl[which(vqtl$SNP %in% names(tmp2)[2]),"Assay"]
g1=range(tmp2[which(tmp2$tert%in%"Low"),"f_21002_0_0"],na.rm=T)
g2=range(tmp2[which(tmp2$tert%in%"Medium"),"f_21002_0_0"],na.rm=T)
g3=range(tmp2[which(tmp2$tert%in%"High"),"f_21002_0_0"],na.rm=T)

# Get ranges of tertiles 
g1=paste0("<",g1[2])
g3=paste0(">", g3[1])
g2=paste0(g2[1], "-", g2[2])

# Get info for facet portion 
data2$overall = i

# Prepare data for plot 
data2$tert=factor(data2$tert, levels=c("Low", "Medium", "High"))
# Define male groups 
if(i %in% c("Males 40-45yr","Males 45-50yr","Males 50-55yr","Males 55-60yr","Males 60-65yr","Males 65-70yr")) { 
    my_col = brewer.pal(n = 5, "Greens")[3:5]
p1=ggplot(data2, aes(x=factor(SNP),y=mean)) +
  geom_point(
    aes(color = tert),
    position = position_dodge(0.3), size =3.3) +   geom_errorbar(
    aes(ymin = LCI, ymax = HCI, color = tert),
    position = position_dodge(0.3), width = 0.4,
    ) +
    scale_x_discrete(name=paste0(rsid), labels = c(paste0(a2,"/",a2, "\n", "(n=",a,")"), paste0(a2,"/",a1, "\n", "(n=",b,")") , paste0(a1,"/",a1, "\n", "(n=",c,")"))) + 
     scale_colour_manual(values = my_col, labels = c(g1,g2,g3), name="Body Weight \n(kg) Tertile") + ylab(paste0("Mean Transformed ",prot, "\n levels [95% CI]")) + ylim(limits=c(-2.5,2.5))+ 
   theme_scientific() +  theme(legend.position="right",axis.title =element_text(size=13), axis.text=element_text(size=13), legend.text=element_text(size=11),legend.title=element_text(size=13))+facet_wrap(~overall)+ theme(strip.background =element_rect(fill="grey94"),strip.text.x = element_text(size = 14))

} else { # Define female groups 
    my_col = brewer.pal(n = 5, "Blues")[3:5]
p1=ggplot(data2, aes(x=factor(SNP),y=mean)) +
  geom_point(
    aes(color = tert),
    position = position_dodge(0.3), size =3.3) +   geom_errorbar(
    aes(ymin = LCI, ymax = HCI, color = tert),
    position = position_dodge(0.3), width = 0.4,
    ) +
    scale_x_discrete(name="", labels = c(paste0(a2,"/",a2, "\n", "(n=",a,")"), paste0(a2,"/",a1, "\n", "(n=",b,")") , paste0(a1,"/",a1, "\n", "(n=",c,")"))) + 
     scale_colour_manual(values=my_col, labels = c(g1,g2,g3), name="Body Weight \n(kg) Tertile") + ylab(paste0("Mean Transformed ",prot, "\n levels [95% CI]"))  +
   theme_scientific() +  theme(legend.position="right",axis.title =element_text(size=13), axis.text=element_text(size=13), legend.text=element_text(size=11),legend.title=element_text(size=13)) +facet_wrap(~overall)+ theme(strip.background =element_rect(fill="grey94"),strip.text.x = element_text(size = 14))

}
# Store output 
myplots[[i]]=p1
} 

# Combine plots 
pdf("/rhillary/dissemination2/PAEP_agesex_noresid_larger.pdf",width=15,height=15)
p1.1=myplots[[1]]
p2=myplots[[2]]
p3=myplots[[3]]
p4=myplots[[4]]
p5=myplots[[5]]
p6=myplots[[6]]
p7=myplots[[7]]
p8=myplots[[8]]
p9=myplots[[9]]
p10=myplots[[10]]
p11=myplots[[11]]
p12=myplots[[12]]

plot1=plot_grid(p10,p12,p1.1,p3,p9,p2,nrow=2,ncol=3)
plot2=plot_grid(p5,p11,p7,p4,p8,p6,nrow=2,ncol=3)
parta=ggarrange(plot1,plot2, nrow=2)
print(parta)
dev.off()
