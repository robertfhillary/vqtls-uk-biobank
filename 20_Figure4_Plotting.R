########################
###### Figure 4 ########
########################

# Load requisite libraries 
library(data.table)
library(ggplot2)
library(cowplot)
library(forcats)
library(RColorBrewer)
library(grafify)
library(artyfarty)
library(circlize)


# Create functions necessary for Fig 4C 
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


####################################
########## PLOT FIGURE 4B ##########
####################################

## 4B is considered first for flow of code

# Read in gei results 
phen1=as.data.frame(fread("/rhillary/outputs/cond_gei_final2023.csv"))

# Merge in their panel and annotation information for plotting 
dis=readRDS("/rhillary/outputs/dis_independent_annotated.rds")
dis=dis[which(!duplicated(dis$compare)),]
dis=dis[,c("compare","Array")]
names(dis)[1]="Protein"
phen1=merge(phen1,dis,by="Protein")

# Tabulate by array 
table1=table(phen1$Array, phen1$Group)
table2=t(table1)
res=apply(table2, 2, function(x) x/sum(x))
res1=as.data.frame(res)
res1$Group=row.names(res1)

# Heatmap of Z scores - Prepare and calculate 
phen1$Z=phen1$Interaction_Beta/phen1$Interaction_SE
# Split by array 
cardio=phen1[phen1$Array %in% "Cardiometabolic",]
inflam=phen1[phen1$Array %in% "Inflammation",]
onc=phen1[phen1$Array %in% "Oncology",]
neuro=phen1[phen1$Array %in% "Neurology",]
# Tidy - this protein has two assigned to it  
phen1[grep("CKMT1A",phen1$Assay),"Assay"]="CKMT1B"

# Prepare data for heatmap 
p1<-ggplot(aes(y=Phenotype.y, x=Assay, fill=Z), data=rbind(cardio,inflam,onc,neuro)) + geom_tile(col="black") +  theme_minimal() +scale_fill_gradient2(high="#D7191C", mid="white", low="#2C7BB6")+ theme(axis.text.x=element_text(angle = -70, hjust = 0)) + xlab("") + ylab("Phenotype")  
p2= p1 +    facet_grid(.~factor(Array), scales = "free", switch = "x", space = "free_x") +  theme(strip.placement = "outside")
phen2=as.data.frame(table(phen1$Assay))
phen3=phen2[which(phen2$Freq >=3),]
phen.plot=phen1[which(phen1$Assay %in% phen3$Var1),]

# Plot stage - IL6 is removed as it has four distinct analytes that each have a GEI (so looks like 4 for IL6) but only one protein associated so doesn't meet criteria 
phen.plot=phen.plot[which(!phen.plot$Assay %in% "IL6"),]
pdf("/rhillary/dissemination2/heatmap_point.pdf",height=8.2,width=12.5)
p1<-ggplot(aes(y=Phenotype_Name, x=Assay, col=Z), data=phen.plot) +   geom_point(aes(col=Z), size =3.4) +  theme_bw() +labs(colour="G*E \nZ-score")+scale_color_gradient(low="#a123cd",high="#13ed6b")+ theme(axis.text.x=element_text(angle = -49, hjust = 0),axis.title=element_text(size=12.5)) + xlab("Protein") + ylab("Phenotype") + ggtitle("Proteins with 3 or more conditionally significant interactions") +   theme(plot.title = element_text(hjust = 0.5, size=14.5),axis.title=element_text(size=15),axis.text.x=element_text(size=11),axis.text.y=element_text(size=10)) 
p2= p1 +     scale_y_discrete(limits=rev) +   facet_grid(.~factor(Array), scales = "free", switch = "y", space = "free_x") +  theme(strip.placement = "outside") +   theme(strip.background =element_rect(fill="whitesmoke"), strip.text.x=element_text(size=10)) 
print(p2)
dev.off()


####################################
########## PLOT FIGURE 4A ##########
####################################

mat=phen1[,c("Assay","Group","Z")]
names(mat)=c("from","to","value")
circos.clear()
pdf("/rhillary/dissemination2/chord_diagram_p16.pdf", height = 13, width =15.5)
par(cex = 0.9, mar = c(0, 0, 0, 0))
chordDiagram(mat, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), col = "black")
}, bg.border = NA)
dev.off()


####################################
########## PLOT FIGURE 4C ##########
####################################

# Reset working directory and obtain LDLR interaction from results 
setwd("/rhillary/interaction_vqtls/")
loop=list.files(".", ".")
new.phenos=as.data.frame(fread("/rhillary/all_quant_20200904.pheno.v2.tsv.gz"))
names(new.phenos)[1]="FID"
tmp=as.data.frame(fread(paste0("/rhillary/interaction_vqtls/",loop[207])))
vqtl=readRDS("/rhillary/outputs/dis_independent_annotated.rds")
vqtl=vqtl[which(vqtl$Protein%in%"LDLR"),]
# Merge phenotype + interaction data together 
tmp1=merge(tmp,new.phenos,by="FID")

# Obtain tertiles of phenotypes 
vTert = quantile(tmp1$f_30780_0_0, c(0:3/3),na.rm=T)

# Classify values according to tertile 
tmp1$tert = with(tmp1, 
               cut(f_30780_0_0, 
                   vTert, 
                   include.lowest = T, 
                   labels = c("Low", "Medium", "High")))

# Tidy up dataframe and format for plot 
dat=as.data.frame(fread("/home/rhillary/dat_wide_internal_use.csv"))
dat1=dat[,c(1,grep("LDLR",names(dat)))]
names(dat1)=c("FID","LDLR.2")
tmp1=merge(tmp1,dat1,by="FID")
names(tmp1)[3]="LDLR"
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

# Obtain info needed for plot 
a=table(tmp1[,2])[1]
b=table(tmp1[,2])[2]
c=table(tmp1[,2])[3]
a1=vqtl[which(vqtl$SNP %in% names(tmp1)[2]),"A1"]
a2=vqtl[which(vqtl$SNP %in% names(tmp1)[2]),"A2"]
rsid=vqtl[which(vqtl$SNP %in% names(tmp1)[2]), "rsid"]
prot=vqtl[which(vqtl$SNP %in% names(tmp1)[2]),"Assay"]

# Get ranges within each tertile 
g1=range(tmp1[which(tmp1$tert%in%"Low"),"f_30780_0_0"],na.rm=T)
g2=range(tmp1[which(tmp1$tert%in%"Medium"),"f_30780_0_0"],na.rm=T)
g3=range(tmp1[which(tmp1$tert%in%"High"),"f_30780_0_0"],na.rm=T)

# Tidy values 
g1=paste0("<",signif(g1[2],2))
g3=paste0(">", signif(g3[1],2))
g2=paste0(signif(g2[1],2), "-", signif(g2[2],2))
data2$tert=factor(data2$tert, levels=c("Low", "Medium", "High"))

# Plot stage 
pdf("/rhillary/dissemination2/interaction_LDLR_transform.pdf",height=7,width=11)
ggplot(data2, aes(x=factor(SNP),y=mean)) +
  geom_point(
    aes(color = tert),
    position = position_dodge(0.3), size =5) +   geom_errorbar(
    aes(ymin = LCI, ymax = HCI, color = tert),
    position = position_dodge(0.3), width = 0.4,
    ) +
    scale_x_discrete(name=paste0("Genotype - ",rsid), labels = c(paste0(a2,"/",a2, "\n", "(n=",a,")"), paste0(a2,"/",a1, "\n", "(n=",b,")") , paste0(a1,"/",a1, "\n", "(n=",c,")"))) + 
     scale_color_brewer(palette = "Set2", labels = c(g1,g2,g3), name="LDL Cholesterol (mmol/L) \nTertile") + ylab(paste0("Mean Transformed ", prot, "\n levels [95% CI]")) + 
        theme_scientific() + scale_y_continuous(limits=c(-1.0,1),breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1)) + theme(legend.position="right",axis.title =element_text(size=16), axis.text=element_text(size=15), legend.text=element_text(size=15),legend.title=element_text(size=15)) 
 dev.off()


