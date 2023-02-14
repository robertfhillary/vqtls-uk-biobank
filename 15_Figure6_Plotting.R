##########################################
######## FIGURE 6 ########################
##########################################

# Load requisite libraries 
library(data.table)
library(ggplot2)
library(artyfarty)
library(RColorBrewer)

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

# Rank-inverse based normal transformation
INT <- function(x) qnorm((rank(x, na.last="keep") - 0.5) / sum(!is.na(x)))


######################################################################
######## FIGURE 6B - 6A is from Biorender.com ########################
######################################################################

# Set up files 
setwd("/rhillary/interaction_vqtls/")
loop=list.files(".", ".")
new.phenos=as.data.frame(fread("/rhillary/all_quant_20200904.pheno.v2.tsv.gz"))
names(new.phenos)[1]="FID"
tmp=as.data.frame(fread(paste0("/rhillary/interaction_vqtls/",loop[66])))
vqtl=readRDS("/rhillary/outputs/dis_independent_annotated.rds")
tmp1=merge(tmp,new.phenos,by="FID")
tmp1=tmp1[which(!is.na(tmp1[,2])),]
a1=vqtl[which(vqtl$SNP %in% names(tmp1)[2]),"A1"]
a2=vqtl[which(vqtl$SNP %in% names(tmp1)[2]),"A2"]
rsid=vqtl[which(vqtl$SNP %in% names(tmp1)[2]),"rsid"]
prot=vqtl[which(vqtl$SNP %in% names(tmp1)[2]),"Assay"]

# Define tertiles of monocyte count 
vTert = quantile(tmp1$f_30130_0_0, c(0:3/3),na.rm=T)

# Classify values
tmp1$tert = with(tmp1, 
               cut(f_30130_0_0, 
                   vTert, 
                   include.lowest = T, 
                   labels = c("Low", "Medium", "High")))

# Prepare files for plotting 
dat=as.data.frame(fread("/rhillary/dat_wide_internal_use.csv"))
dat1=dat[,c(1,grep("FLT3LG",names(dat)))]
names(dat1)=c("FID","FLT3LG")
tmp1=merge(tmp1,dat1,by="FID")
names(tmp1)[3]="FLT3LG"
tmp1=tmp1[which(!is.na(tmp1$tert)),]
tmp1=tmp1[which(!is.na(tmp1[,2])),]
tmp1=tmp1[which(!is.na(tmp1[,3])),]

# Get summary data across genotype groups 
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
g1=range(tmp1[which(tmp1$tert%in%"Low"),"f_30130_0_0"],na.rm=T)
g2=range(tmp1[which(tmp1$tert%in%"Medium"),"f_30130_0_0"],na.rm=T)
g3=range(tmp1[which(tmp1$tert%in%"High"),"f_30130_0_0"],na.rm=T)

# Get ranges of tertiles 
g1=paste0("<",g1[2])
g3=paste0(">", g3[1])
g2=paste0(g2[1], "-", g2[2])

# Prep for plot 
data2$tert=factor(data2$tert, levels=c("Low", "Medium", "High"))
# Plot stage 
pdf("/rhillary/dissemination2/interaction_ftl3lg2.pdf",height=7,width=11)
ggplot(data2, aes(x=factor(SNP),y=mean)) +
  geom_point(
    aes(color = tert),
    position = position_dodge(0.3), size =5) +   geom_errorbar(
    aes(ymin = LCI, ymax = HCI, color = tert),
    position = position_dodge(0.3), width = 0.4,
    ) +
    scale_x_discrete(name=paste0("Genotype - ",rsid), labels = c(paste0(a2,"/",a2, "\n", "(n=",a,")"), paste0(a2,"/",a1, "\n", "(n=",b,")") , paste0(a1,"/",a1, "\n", "(n=",c,")"))) + 
     scale_color_brewer(palette = "Set2", labels = c(g1,g2,g3), name=expression(atop(paste("Monocyte Count"), paste("\n (", 10^9, "cells/Litre) Tertile")))) + ylab(paste0("Mean Transformed ", prot, "\n levels [95% CI]")) + 
        theme_scientific() + scale_y_continuous(limits=c(-0.5,0.5),breaks=c(-0.5,-0.25,0,0.25,0.5)) + theme(legend.position="right",axis.title =element_text(size=16), axis.text=element_text(size=15), legend.text=element_text(size=15),legend.title=element_text(size=15)) + theme(legend.title.align=0)
 dev.off()


############################
###### FIGURE 6C ###########
############################

# Correlations between monocyte count and FLT3LG by genotype group
tmper=tmp1
tmper$f_30130_0_0=outlierTrim(tmper$f_30130_0_0)
tmper=tmper[which(!is.na(tmper$f_30130_0_0)),]
tmper=tmper[which(!is.na(tmper$FLT3LG)),]
names(tmper)[603]="FLT3LG.2"
tmper$label=NA
a=table(tmper[,2])[1]
b=table(tmper[,2])[2]
c=table(tmper[,2])[3]
tmper[tmper[,2]%in%0,"label"]=paste0(a2,"/",a2, "\n", "(n=",a,")")
tmper[tmper[,2]%in%1,"label"]=paste0(a2,"/",a1, "\n", "(n=",b,")")
tmper[tmper[,2]%in%2,"label"]=paste0(a1,"/",a1, "\n", "(n=",c,")")

# Plot Stage 
pdf("/rhillary/dissemination2/correlation_monocyte_counts.pdf",height=6.5,width=8)
c1=ggplot(tmper) + geom_point(size=1,aes(y=FLT3LG, x=f_30130_0_0, color=factor(tmper[,2]))) + scale_colour_manual(values = brewer.pal(4, "Purples")[2:4]) + 
 geom_smooth(aes(y=FLT3LG, x=f_30130_0_0, color=factor(tmper[,2])), method="lm",se=F,color="gray17") + scale_y_continuous(limits=c(-5,5),,breaks=c(-5,-2.5,0,2.5,5)) + ylab("Transformed FLT3LG levels") + labs(x = expression(paste("Monocyte Count (", 10^9, "cells/Litre)"))) +theme_scientific() + facet_wrap(~label) + theme(legend.position="none",axis.title =element_text(size=15), axis.text=element_text(size=13)) +   theme(strip.background =element_rect(fill="white"))
print(c1)
dev.off() 


############################
###### FIGURE 6D ###########
############################

# Correlations between FLT3 and FLT3LG by genotype group
# Merge in FLT3 
flt3=as.data.frame(fread("/rhillary/discovery_phenotypes_patched/FLT3_P36888_OID21272_v1_Oncology.phen"))
flt3$IID=NULL 
names(flt3)[2]="FLT3"
tmp1=merge(tmp1,flt3,by="FID")
tmper=tmp1
tmper=tmper[which(!is.na(tmper$FLT3)),]
tmper=tmper[which(!is.na(tmper$FLT3LG)),]
a=table(tmper[,2])[1]
b=table(tmper[,2])[2]
c=table(tmper[,2])[3]
tmper$label=NA
tmper[tmper[,2]%in%0,"label"]=paste0(a2,"/",a2, "\n", "(n=",a,")")
tmper[tmper[,2]%in%1,"label"]=paste0(a2,"/",a1, "\n", "(n=",b,")")
tmper[tmper[,2]%in%2,"label"]=paste0(a1,"/",a1, "\n", "(n=",c,")")

# Plot Stage 
pdf("/rhillary/dissemination2/correlation_proteins.pdf",height=6.5,width=8)
c1=ggplot(tmper, aes(y=FLT3LG, x=FLT3, group=factor(tmper[,2]))) + geom_point(size=1,aes(color=factor(tmper[,2]))) +  scale_colour_manual(values = brewer.pal(4, "Purples")[2:4]) + 
 geom_smooth(method="lm",se=F,color="gray17") + ylab("Transformed FLT3LG levels") +  xlab("Transformed FLT3 levels") + scale_y_continuous(limits=c(-5,5),breaks=c(-5,-2.5,0,2.5,5)) +scale_x_continuous(limits=c(-5,5),,breaks=c(-5,-2.5,0,2.5,5)) + theme_scientific() + facet_wrap(~label) + theme(legend.position="none",axis.title =element_text(size=15), axis.text=element_text(size=13)) +   theme(strip.background =element_rect(fill="white"))
print(c1)
dev.off() 


############################
###### FIGURE 6E ###########
############################
# Is monocyte count different by each genotype group? 
tmper=tmp1
tmper$f_30130_0_0=outlierTrim(tmper$f_30130_0_0)
tmper=tmper[which(!is.na(tmper$f_30130_0_0)),]
a=table(tmper[,2])[1]
b=table(tmper[,2])[2]
c=table(tmper[,2])[3]
names(tmper)[603]="FLT3LG.2"

# Plot Stage 
pdf("/rhillary/dissemination2/monocytes.pdf",height=6,width=7)
ggplot(tmper, aes(x=factor(tmper[,2]), y=f_30130_0_0, group=factor(tmper[,2]))) + geom_violin(aes(color=factor(tmper[,2])),width=0.4,position = position_dodge(0.9)) + geom_boxplot(aes(color=factor(tmper[,2])),width=0.16,outlier.shape = 1, outlier.size = 1.5,position = position_dodge(0.9)) +  scale_colour_manual(values = brewer.pal(5, "Purples")[3:5]) + 
 theme_scientific() + theme(legend.position="none") + labs(y = expression(paste("Monocyte Count (", 10^9, "cells/Litre)"))) +  theme(strip.background =element_rect(fill="white")) + scale_x_discrete(name=paste0("Genotype - ",rsid), labels = c(paste0(a2,"/",a2, "\n", "(n=",a,")"), paste0(a2,"/",a1, "\n", "(n=",b,")") , paste0(a1,"/",a1, "\n", "(n=",c,")"))) + 
theme(axis.title =element_text(size=15), axis.text=element_text(size=14))
dev.off() 


####################################################################
###############  FIGURE S3 - DIFFERENT BLOOD CELL TYPES ############
####################################################################

# Prepare files - need a different one for each plot 
names(tmp1)[3]="FLT3LG"
names(tmp1)[603]="FLT3LG.2"
tmper=tmp1
tmper=tmper[which(!is.na(tmper$FLT3LG)),]

###############
# Neutrophils #
###############
tmper1=tmper
tmper1$f_30140_0_0=outlierTrim(tmper1$f_30140_0_0)
tmper1=tmper1[which(!is.na(tmper1$f_30140_0_0)),]
a=table(tmper1[,2])[1]
b=table(tmper1[,2])[2]
c=table(tmper1[,2])[3]
tmper1$label=NA
tmper1[tmper1[,2]%in%0,"label"]=paste0(a2,"/",a2, "\n", "(n=",a,")")
tmper1[tmper1[,2]%in%1,"label"]=paste0(a2,"/",a1, "\n", "(n=",b,")")
tmper1[tmper1[,2]%in%2,"label"]=paste0(a1,"/",a1, "\n", "(n=",c,")")
# Plot 
p1=ggplot(tmper1, aes(y=FLT3LG, x=f_30140_0_0, group=factor(tmper1[,2]))) + geom_point(size=1,aes(color=factor(tmper1[,2]))) +  scale_colour_manual(values = brewer.pal(4, "Blues")[2:4]) + 
 geom_smooth(method="lm",se=F,color="grey55") + scale_y_continuous(limits=c(-5,5)) + ylab("Transformed FLT3LG levels") + labs(x = expression(paste("Neutrophil Count (", 10^9, "cells/Litre)"))) +theme_scientific() + facet_wrap(~label) + theme(legend.position="none",axis.title =element_text(size=13), axis.text=element_text(size=12)) +   theme(strip.background =element_rect(fill="white"))

#############
# Basophils #
#############
tmper2=tmper
tmper2$f_30160_0_0=outlierTrim(tmper2$f_30160_0_0)
tmper2=tmper2[which(!is.na(tmper2$f_30160_0_0)),]
a=table(tmper2[,2])[1]
b=table(tmper2[,2])[2]
c=table(tmper2[,2])[3]
tmper2$label=NA
tmper2[tmper2[,2]%in%0,"label"]=paste0(a2,"/",a2, "\n", "(n=",a,")")
tmper2[tmper2[,2]%in%1,"label"]=paste0(a2,"/",a1, "\n", "(n=",b,")")
tmper2[tmper2[,2]%in%2,"label"]=paste0(a1,"/",a1, "\n", "(n=",c,")")
# Plot 
p2=ggplot(tmper2, aes(y=FLT3LG, x=f_30160_0_0, group=factor(tmper2[,2]))) + geom_point(size=1,aes(color=factor(tmper2[,2]))) +  scale_colour_manual(values = brewer.pal(4, "Blues")[2:4]) + 
 geom_smooth(method="lm",se=F,color="grey55") + scale_y_continuous(limits=c(-5,5)) + ylab("Transformed FLT3LG levels") + labs(x = expression(paste("Basophil Count (", 10^9, "cells/Litre)"))) +theme_scientific() + facet_wrap(~label) + theme(legend.position="none",axis.title =element_text(size=13), axis.text=element_text(size=12)) +   theme(strip.background =element_rect(fill="white"))

###############
# Eosinophils #
###############
tmper3=tmper
tmper3$f_30150_0_0=outlierTrim(tmper3$f_30150_0_0)
tmper3=tmper3[which(!is.na(tmper3$f_30150_0_0)),]
a=table(tmper3[,2])[1]
b=table(tmper3[,2])[2]
c=table(tmper3[,2])[3]
tmper3$label=NA
tmper3[tmper3[,2]%in%0,"label"]=paste0(a2,"/",a2, "\n", "(n=",a,")")
tmper3[tmper3[,2]%in%1,"label"]=paste0(a2,"/",a1, "\n", "(n=",b,")")
tmper3[tmper3[,2]%in%2,"label"]=paste0(a1,"/",a1, "\n", "(n=",c,")")
# Plot 
p3=ggplot(tmper3, aes(y=FLT3LG, x=f_30150_0_0, group=factor(tmper3[,2]))) + geom_point(size=1,aes(color=factor(tmper3[,2]))) +  scale_colour_manual(values = brewer.pal(4, "Blues")[2:4]) + 
 geom_smooth(method="lm",se=F,color="grey55") + scale_y_continuous(limits=c(-5,5)) + ylab("Transformed FLT3LG levels") + labs(x = expression(paste("Eosinophil Count (", 10^9, "cells/Litre)"))) +theme_scientific() + facet_wrap(~label) + theme(legend.position="none",axis.title =element_text(size=13), axis.text=element_text(size=12)) +   theme(strip.background =element_rect(fill="white"))

#################
# Reticulocytes #
#################
tmper4=tmper
tmper4$f_30250_0_0=outlierTrim(tmper4$f_30250_0_0)
tmper4=tmper4[which(!is.na(tmper4$f_30250_0_0)),]
a=table(tmper4[,2])[1]
b=table(tmper4[,2])[2]
c=table(tmper4[,2])[3]
tmper4$label=NA
tmper4[tmper4[,2]%in%0,"label"]=paste0(a2,"/",a2, "\n", "(n=",a,")")
tmper4[tmper4[,2]%in%1,"label"]=paste0(a2,"/",a1, "\n", "(n=",b,")")
tmper4[tmper4[,2]%in%2,"label"]=paste0(a1,"/",a1, "\n", "(n=",c,")")
# Plot 
p4=ggplot(tmper4, aes(y=FLT3LG, x=f_30250_0_0, group=factor(tmper4[,2]))) + geom_point(size=1,aes(color=factor(tmper4[,2]))) +  scale_colour_manual(values = brewer.pal(4, "Blues")[2:4]) + 
 geom_smooth(method="lm",se=F,color="grey55") + scale_y_continuous(limits=c(-5,5)) + ylab("Transformed FLT3LG levels") + labs(x = expression(paste("Reticulocyte Count (", 10^12, "cells/Litre)"))) +theme_scientific() + facet_wrap(~label) + theme(legend.position="none",axis.title =element_text(size=13), axis.text=element_text(size=12)) +   theme(strip.background =element_rect(fill="white"))

###############
# Lymphocytes #
###############  
tmper5=tmper
tmper5$f_30120_0_0=outlierTrim(tmper5$f_30120_0_0)
tmper5=tmper5[which(!is.na(tmper5$f_30120_0_0)),]
a=table(tmper5[,2])[1]
b=table(tmper5[,2])[2]
c=table(tmper5[,2])[3]
tmper5$label=NA
tmper5[tmper5[,2]%in%0,"label"]=paste0(a2,"/",a2, "\n", "(n=",a,")")
tmper5[tmper5[,2]%in%1,"label"]=paste0(a2,"/",a1, "\n", "(n=",b,")")
tmper5[tmper5[,2]%in%2,"label"]=paste0(a1,"/",a1, "\n", "(n=",c,")")
# Plot 
p5=ggplot(tmper5, aes(y=FLT3LG, x=f_30120_0_0, group=factor(tmper5[,2]))) + geom_point(size=1,aes(color=factor(tmper5[,2]))) +  scale_colour_manual(values = brewer.pal(4, "Blues")[2:4]) + 
 geom_smooth(method="lm",se=F,color="grey55") + scale_y_continuous(limits=c(-5,5)) + ylab("Transformed FLT3LG levels") + labs(x = expression(paste("Lymphocyte Count (", 10^9, "cells/Litre)"))) +theme_scientific() + facet_wrap(~label) + theme(legend.position="none",axis.title =element_text(size=13), axis.text=element_text(size=12)) +   theme(strip.background =element_rect(fill="white"))

##############
# Leukocytes #
##############  
tmper6=tmper
tmper6$f_30000_0_0=outlierTrim(tmper6$f_30000_0_0)
tmper6=tmper6[which(!is.na(tmper6$f_30000_0_0)),]
a=table(tmper6[,2])[1]
b=table(tmper6[,2])[2]
c=table(tmper6[,2])[3]
tmper6$label=NA
tmper6[tmper6[,2]%in%0,"label"]=paste0(a2,"/",a2, "\n", "(n=",a,")")
tmper6[tmper6[,2]%in%1,"label"]=paste0(a2,"/",a1, "\n", "(n=",b,")")
tmper6[tmper6[,2]%in%2,"label"]=paste0(a1,"/",a1, "\n", "(n=",c,")")
# Plot 
p6=ggplot(tmper6, aes(y=FLT3LG, x=f_30000_0_0, group=factor(tmper6[,2]))) + geom_point(size=1,aes(color=factor(tmper6[,2]))) +  scale_colour_manual(values = brewer.pal(4, "Blues")[2:4]) + 
 geom_smooth(method="lm",se=F,color="grey55") + scale_y_continuous(limits=c(-5,5)) + ylab("Transformed FLT3LG levels") + labs(x = expression(paste("Leukocyte Count (", 10^9, "cells/Litre)"))) +theme_scientific() + facet_wrap(~label) + theme(legend.position="none",axis.title =element_text(size=13), axis.text=element_text(size=12)) +   theme(strip.background =element_rect(fill="white"))

#########################
# Red blood cell counts #
#########################
tmper7=tmper
tmper7$f_30010_0_0=outlierTrim(tmper7$f_30010_0_0)
tmper7=tmper7[which(!is.na(tmper7$f_30010_0_0)),]
a=table(tmper7[,2])[1]
b=table(tmper7[,2])[2]
c=table(tmper7[,2])[3]
tmper7$label=NA
tmper7[tmper7[,2]%in%0,"label"]=paste0(a2,"/",a2, "\n", "(n=",a,")")
tmper7[tmper7[,2]%in%1,"label"]=paste0(a2,"/",a1, "\n", "(n=",b,")")
tmper7[tmper7[,2]%in%2,"label"]=paste0(a1,"/",a1, "\n", "(n=",c,")")
# Plot 
p7=ggplot(tmper7, aes(y=FLT3LG, x=f_30010_0_0, group=factor(tmper7[,2]))) + geom_point(size=1,aes(color=factor(tmper7[,2]))) +  scale_colour_manual(values = brewer.pal(4, "Blues")[2:4]) + 
 geom_smooth(method="lm",se=F,color="grey55") + scale_y_continuous(limits=c(-5,5)) + ylab("Transformed FLT3LG levels") + labs(x = expression(paste("Erythrocyte Count (", 10^12, "cells/Litre)"))) +theme_scientific() + facet_wrap(~label) + theme(legend.position="none",axis.title =element_text(size=13), axis.text=element_text(size=12)) +   theme(strip.background =element_rect(fill="white"))

###################
# Platelet counts #
###################
tmper8=tmper
tmper8$f_30080_0_0=outlierTrim(tmper8$f_30080_0_0)
tmper8=tmper8[which(!is.na(tmper8$f_30080_0_0)),]
a=table(tmper8[,2])[1]
b=table(tmper8[,2])[2]
c=table(tmper8[,2])[3]
tmper8$label=NA
tmper8[tmper8[,2]%in%0,"label"]=paste0(a2,"/",a2, "\n", "(n=",a,")")
tmper8[tmper8[,2]%in%1,"label"]=paste0(a2,"/",a1, "\n", "(n=",b,")")
tmper8[tmper8[,2]%in%2,"label"]=paste0(a1,"/",a1, "\n", "(n=",c,")")
# Plot 
p8=ggplot(tmper8, aes(y=FLT3LG, x=f_30080_0_0, group=factor(tmper8[,2]))) + geom_point(size=1,aes(color=factor(tmper8[,2]))) +  scale_colour_manual(values = brewer.pal(4, "Blues")[2:4]) + 
 geom_smooth(method="lm",se=F,color="grey55") + scale_y_continuous(limits=c(-5,5)) + ylab("Transformed FLT3LG levels") + labs(x = expression(paste("Platelet Count (", 10^9, "cells/Litre)"))) +theme_scientific() + facet_wrap(~label) + theme(legend.position="none",axis.title =element_text(size=13), axis.text=element_text(size=12)) +   theme(strip.background =element_rect(fill="white"))

# Combine all 
pdf("/rhillary/dissemination2/correlation_allbloodcells_counts.pdf",height=14,width=10)
g1=plot_grid(p4,p5,p2,p3,p1,p6,p7,p8, nrow=4,ncol=2, labels =c("A","B","C","D","E","F","G","H"))
print(g1)
dev.off()





