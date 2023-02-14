########################
###### Figure 3 ########
########################

# Load requisite libraries 
library(data.table)
library(ggplot2)
library(cowplot)
library(forcats)
library(RColorBrewer)
library(grafify)
library(artyfarty)

#################################
####### FIGURE 3A ###############
#################################

# Read in vQTL data 
vqtl=as.data.frame(fread("/rhillary/outputs/dis_independent_annotated_main_effect.csv"))
big.me=vqtl[grep("FLT3LG", vqtl$Assay),]
## both 
chr=big.me[,"Chr"]
file=big.me[,"compare"]
snp=big.me[,"SNP"]
a1=big.me[,"A1"]
a2=big.me[,"A2"]
rsid=big.me[,"rsid"]
x=as.data.frame(fread(paste0("/rhillary/recode/chr",chr,".raw")))
names(x)=gsub("_.*", "", names(x))
y=read.table(paste0("/rhillary/discovery_phenotypes_patched/",file, ".phen"),header=T)
x2=x[,c(1,2,which(names(x) %in% snp))]
x3=merge(x2,y[,c(1,3)],by="FID")
x3=x3[which(!is.na(x3[,4])),]
x3=x3[which(!is.na(x3[,3])),]
a=length(which(x3[,3]==0))
b=length(which(x3[,3]==1))
c=length(which(x3[,3]==2))

# Plot stage 
pdf("/rhillary/dissemination2/Fig3A_FLT3LG.pdf",width=6.2,height=5.8)
plot2= ggplot(x3,aes(factor(x3[,3]),x3[,4],color=factor(x3[,3])))+  geom_jitter(width = 0.03)+geom_boxplot()+
scale_x_discrete(name=paste0("Genotype - ",rsid), labels = c(paste0(a2,"/",a2, "\n", "(n=",a,")"), paste0(a2,"/",a1, "\n", "(n=",b,")") , paste0(a1,"/",a1, "\n", "(n=",c,")"))) + 
ylab(paste0("Transformed ",big.me$Assay, " levels")) + ggtitle("vQTL Effect Only") + guides(fill="none") + theme_scientific() + theme(plot.title = element_text(hjust = 0.5),axis.text=element_text(size=12),axis.title=element_text(size=13))
plot2=plot2+ scale_color_brewer(palette = "Set2") + guides(color=FALSE)+scale_y_continuous(limits=c(-5,5),breaks=c(-5,-4,-3,-2,-1,0,1,2,3,4,5))
plot4=plot2 +annotate(geom="text", x=3, y=4.7, label = "atop(vQTL~italic(p) == 2.1~x10^-16,Main~effect~italic(p) == 0.95)", parse = TRUE,size=4)
print(plot4)
dev.off()


#################################
####### FIGURE 3B ###############
#################################

# vQTL and main effect QTL together 
vqtl=as.data.frame(fread("/rhillary/outputs/dis_independent_annotated_main_effect.csv"))
big.me1=vqtl[grep("ACP6", vqtl$Assay),]
## both 
chr1=big.me1[,"Chr"]
file1=big.me1[,"compare"]
snp1=big.me1[,"SNP"]
a11=big.me1[,"A1"]
a21=big.me1[,"A2"]
rsid1=big.me1[,"rsid"]
x1=as.data.frame(fread(paste0("/scratch/rhillary/recode/chr",chr1,".raw")))
names(x1)=gsub("_.*", "", names(x1))
y1=read.table(paste0("/home/rhillary/main_effect_phenotypes/",file1, ".phen"),header=T)
x21=x1[,c(1,2,which(names(x1) %in% snp1))]
x31=merge(x21,y1[,c(1,3)],by="FID")
x31=x31[which(!is.na(x31[,4])),]
x31=x31[which(!is.na(x31[,3])),]
a1=length(which(x31[,3]==0))
b1=length(which(x31[,3]==1))
c1=length(which(x31[,3]==2))

# Plot stage 
pdf("/rhillary/dissemination2/Fig3B_ACP6.pdf",width=6.2,height=5.8)
plot2.1= ggplot(x31,aes(factor(x31[,3]),x31[,4],color=factor(x31[,3])))+  geom_jitter(width = 0.03)+geom_boxplot()+
scale_x_discrete(name=paste0("Genotype - ",rsid1), labels = c(paste0(a21,"/",a21, "\n", "(n=",a1,")"), paste0(a21,"/",a11, "\n", "(n=",b1,")") , paste0(a11,"/",a11, "\n", "(n=",c1,")"))) + 
ylab(paste0("Transformed ",big.me1$Assay, " levels")) + ggtitle("vQTL and Main Effect QTL") + guides(fill="none") + theme_scientific() + theme(plot.title = element_text(hjust = 0.5),axis.text=element_text(size=12),axis.title=element_text(size=13))
plot2.1=plot2.1+ scale_color_brewer(palette = "Set2") + guides(color=FALSE)+scale_y_continuous(limits=c(-5,5),breaks=c(-5,-4,-3,-2,-1,0,1,2,3,4,5))
plot3=plot2.1 +annotate(geom="text", x=3, y=4.7, label = "atop(vQTL~italic(p) == 1.0~x10^-300, Main~effect~italic(p) == 1.0~x10^-300)", parse = TRUE,size=4)
print(plot3)
dev.off()


###################################
####### FIGURE 3C+D ###############
###################################

# Read in vqtl data 
vqtl=as.data.frame(fread("/rhillary/outputs/dis_independent_annotated_main_effect.csv"))
# Read in pqtl data 
pqtl=as.data.frame(fread("/rhillary/outputs/meQTL_2_vQTL.csv"))
pqtl=pqtl[-which(pqtl$Annotated.gene.consequence %in% ""),]

# Reformat vQTL annotations 
vqtl[vqtl$Consequence %in% c("3_prime_UTR_variant"), "Consequence"]="3' UTR"
vqtl[vqtl$Consequence %in% c("5_prime_UTR_variant","splice_region_variant,5_prime_UTR_variant"), "Consequence"]="5' UTR"
vqtl[vqtl$Consequence %in% c("non_coding_transcript_exon_variant"," missense_variant,splice_region_variant","missense_variant", "missense_variant,NMD_transcript_variant", "missense_variant","missense_variant,stop_retained_variant", "stop_gained", "synonymous_variant"), "Consequence"]="Exon"
vqtl[vqtl$Consequence %in% c("intergenic_variant"), "Consequence"]="Intergenic"
vqtl[vqtl$Consequence %in% c("intron_variant", "intron_variant,NMD_transcript_variant", "intron_variant,non_coding_transcript_variant", "splice_region_variant,intron_variant"), "Consequence"]="Intron"
vqtl[vqtl$Consequence %in% c("frameshift_variant,splice_region_variant,intron_variant", "splice_region_variant,non_coding_transcript_exon_variant","splice_region_variant", "missense_variant,splice_region_variant"), "Consequence"]="Splice Region"
       
# Reformat pQTL annotations 
pqtl[pqtl$Annotated.gene.consequence %in% c("3_prime_UTR_variant","3_prime_UTR_variant&NMD_transcript_variant"), "Annotated.gene.consequence"]="3' UTR"
pqtl[pqtl$Annotated.gene.consequence %in% c("5_prime_UTR_variant","5_prime_UTR_variant&NMD_transcript_variant"), "Annotated.gene.consequence"]="5' UTR"
#pqtl[pqtl$Annotated.gene.consequence %in% c("downstream_gene_variant"), "Annotated.gene.consequence"]="Downstream"
pqtl[pqtl$Annotated.gene.consequence %in% c("frameshift_variant","frameshift_variant&splice_region_variant&intron_variant","inframe_deletion","inframe_insertion","missense_variant","missense_variant&NMD_transcript_variant","missense_variant&splice_region_variant","non_coding_transcript_exon_variant","stop_gained","stop_lost","synonymous_variant"), "Annotated.gene.consequence"]="Exon"
pqtl[pqtl$Annotated.gene.consequence %in% c("intergenic_variant","upstream_gene_variant", "downstream_gene_variant", "regulatory_region_variant"), "Annotated.gene.consequence"]="Intergenic"
pqtl[pqtl$Annotated.gene.consequence %in% c("intron_variant","intron_variant&NMD_transcript_variant","intron_variant&non_coding_transcript_variant"), "Annotated.gene.consequence"]="Intron"
#pqtl[pqtl$Annotated.gene.consequence %in% c("upstream_gene_variant"), "Annotated.gene.consequence"]="Upstream"
pqtl[pqtl$Annotated.gene.consequence %in% c("splice_polypyrimidine_tract_variant&intron_variant","splice_acceptor_variant","splice_acceptor_variant&splice_polypyrimidine_tract_variant&coding_sequence_variant&intron_variant&NMD_transcript_variant","splice_donor_5th_base_variant&intron_variant","splice_donor_5th_base_variant&intron_variant&non_coding_transcript_variant","splice_donor_region_variant&intron_variant","splice_donor_region_variant&intron_variant&NMD_transcript_variant","splice_donor_region_variant&intron_variant&non_coding_transcript_variant", "splice_donor_variant","splice_donor_variant&non_coding_transcript_variant"," splice_polypyrimidine_tract_variant&intron_variant","splice_polypyrimidine_tract_variant&intron_variant&NMD_transcript_variant","splice_polypyrimidine_tract_variant&intron_variant&non_coding_transcript_variant","splice_polypyrimidine_tract_variant&splice_region_variant&intron_variant","splice_polypyrimidine_tract_variant&splice_region_variant&intron_variant&non_coding_transcript_variant","splice_region_variant&intron_variant","splice_region_variant&intron_variant&non_coding_transcript_variant","splice_region_variant&non_coding_transcript_exon_variant","splice_region_variant&synonymous_variant"), "Annotated.gene.consequence"]="Splice Region"
pqtl[pqtl$Annotated.gene.consequence %in% c("TF_binding_site_variant"), "Annotated.gene.consequence"]="TF Binding Site"

# Tabulate information
t1=as.data.frame(table(vqtl$Consequence))
t2=as.data.frame(table(pqtl$Annotated.gene.consequence))
t1$Prop=t1$Freq/sum(t1$Freq)
t2$Prop=t2$Freq/sum(t2$Freq)
missing=as.data.frame(as.character(t2[which(!t2$Var1 %in% t1$Var1),"Var1"]))
names(missing)[1]="Var1"
missing$Freq=NA
missing$Prop=0 
t1=rbind(t1,missing)

t1$Var1=as.character(t1$Var1)
t2$Var1=as.character(t2$Var1)
t1=t1[order(t1$Var1),]
t2=t2[order(t2$Var1),]
t1$Prop=t1$Prop*100
t2$Prop=t2$Prop*100

# Plots 
pdf("/rhillary/dissemination2/Fig3C.pdf",width=6.2,height=5.8)
p1=ggplot(data=t1, aes(x=Var1, y=Prop)) +
  geom_bar(stat="identity",aes(fill=Var1)) + xlab("Annotated Variant Consequence") + ylab("Proportion (%)") + theme_scientific()+theme(axis.text=element_text(size=10.5), axis.title=element_text(size=13)) + scale_fill_grafify(palette="okabe_ito")+theme(legend.position="none") + theme(axis.text.x = element_text(angle = 45,hjust=1))+ geom_text(aes(label=paste0(round(Prop,1), " %")), position=position_dodge(width=0.9), vjust=-0.25)
print(p1)
dev.off()
pdf("/rhillary/dissemination2/Fig3D.pdf",width=6.2,height=5.8)
p2=ggplot(data=t2, aes(x=Var1, y=Prop)) +
  geom_bar(stat="identity",aes(fill=Var1)) + xlab("Annotated Variant Consequence") + ylab("Proportion (%)") + theme_scientific()+theme(axis.text=element_text(size=10.5), axis.title=element_text(size=13)) + scale_fill_grafify(palette="okabe_ito")+theme(legend.position="none") + theme(axis.text.x = element_text(angle = 45, hjust=1))+ geom_text(aes(label=paste0(round(Prop,1), " %")), position=position_dodge(width=0.9), vjust=-0.25)
print(p2)
dev.off()
              

# Plot all four panels together 
tiff("/rhillary/dissemination/fig3.tiff", units="in", width=10, height=10, res=300)
plot_grid(plot4,plot3,p1,p2,nrow=2,ncol=2,labels=c("A","B","C","D"),label_size = 20)                                                                                               
dev.off()                  


                                         
                                                                                    
                                                       
                                                                                               
                                                                                                   
                        