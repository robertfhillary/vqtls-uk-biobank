##########################################
######## FIGURE 2A-C #####################
##########################################

# Load requisite libraries 
library(pacman)
p_load(tidyverse,
       data.table, magrittr, tools, ggthemes)


# Read in olink information file 
olink_prot_map4 = as.data.frame(fread("/rhillary/olink_protein_map_1.5k_cleaned.tsv"))
olink_prot_map4$UKBPPP.ProteinID=gsub("\\:","_",olink_prot_map4$UKBPPP_ProteinID)
# Fix proteins with names that alternate across files 
olink_prot_map4[grep("LGALS7_LGALS7B_P47929_OID21406_v1",olink_prot_map4$UKBPPP.ProteinID),"compare"]="LGALS7B_P47929_OID21406_v1_Oncology"
olink_prot_map4[grep("MICB_MICA_Q29980_Q29983_OID20593_v1",olink_prot_map4$UKBPPP.ProteinID),"compare"]="MICA_Q29980_OID20593_v1_Inflammation"
# Read in loci file 
pQTLS_ALL_3_4e11_top_loci2 = as.data.frame(read_rds("/rhillary/outputs/discovery_p3.4e11_clumped.rds"))
# Clean file 
maps1=olink_prot_map4[,c("compare", "chr", "gene_start","Assay")]
maps1=maps1[which(!duplicated(maps1$compare)),]
olink_prot_map4=olink_prot_map4[which(!duplicated(olink_prot_map4$compare)),]
# Combine files to help with plotting step 
pQTLS_ALL_3_4e11_top_loci2=merge(pQTLS_ALL_3_4e11_top_loci2,maps1,by="compare")

# Tidy further and add in Cis/Trans information 
pQTLS_ALL_3_4e11_top_loci2$Type=ifelse(pQTLS_ALL_3_4e11_top_loci2$Chr==pQTLS_ALL_3_4e11_top_loci2$chr,"Cis", "Trans")
pQTLS_ALL_3_4e11_top_loci2$diff=NA
cis=pQTLS_ALL_3_4e11_top_loci2[which(pQTLS_ALL_3_4e11_top_loci2$Type%in%"Cis"),]
trans=pQTLS_ALL_3_4e11_top_loci2[which(pQTLS_ALL_3_4e11_top_loci2$Type%in%"Trans"),]
cis$diff=abs(cis$bp-as.numeric(cis$gene_start))
cis$Type = ifelse(cis$diff <=1e6, "Cis","Trans")
pQTLS_ALL_3_4e11_top_loci2=rbind(cis,trans)


####################################
########## PLOT FIGURE 2A ##########
####################################

pQTLS_ALL_3_4e11_top_loci_matrixplot = pQTLS_ALL_3_4e11_top_loci2 %>% 
  left_join(olink_prot_map4 %>% distinct(compare, chr, gene_start), by = c("compare"))

# Obtain chromosome lengths for y-axis
chr_lengths = fread("UKB_PPP_pQTLs/analysis/input/chr_lengths_hg38.txt") %>% select(1,2) %>% setNames(c("CHR", "CHR_length")) %>% 
  as_tibble() %>% 
  mutate(CHR = gsub("X", "23", CHR))
pQTLS_ALL_3_4e11_top_loci_matrixplot$log10P=-log10(pQTLS_ALL_3_4e11_top_loci_matrixplot$P)
lead_3_4e11_matrixplot = pQTLS_ALL_3_4e11_top_loci_matrixplot %>% data.frame() %>% 
  filter(!(Chr != chr.x & Type == "Cis")) %>% 
  group_by(SNP, compare) %>% 
  filter(Type == min(Type)) %>% 
  ungroup %>% 
  dplyr::rename(CHR = Chr, POS = bp, logP = log10P, CIS_TRANS = Type, UNIPROT = compare,
                CHR_gene = chr.x, start_position = gene_start.x) %>% 
  mutate(CHR_gene = gsub("X", "23", CHR_gene)) %>% 
  left_join(chr_lengths, by = "CHR") %>%
  # mutate(logP2 = ifelse(logP > 300, 300,logP),
  #        trunc = ifelse(logP > 300, "diamond-open","circle-open")) %>%
  mutate(pos_plot = as.character(POS/CHR_length)) %>%
  mutate(pos_plot = gsub("0[.]",".",pos_plot)) %>%
  mutate(pos_plot2 = as.numeric(paste0(CHR, pos_plot))) %>%
  mutate(prot_index = 1:n()) %>%
  # mutate(logP2 = ifelse(logP > 300, 300,logP),
  #        trunc = ifelse(logP > 300, "diamond-open","circle-open")) %>%
  # mutate(P_META_chr = conv_log_bs(-logP)) %>%
  left_join(chr_lengths %>%
              dplyr::select(CHR_gene = CHR, CHR_length_gene = CHR_length),
            by = "CHR_gene") %>%
  mutate(start_pos_plot = as.character(as.numeric(start_position)/CHR_length_gene)) %>%
  mutate(start_pos_plot = gsub("0[.]",".",start_pos_plot)) %>%
  mutate(start_pos_plot2 = as.numeric(paste0(CHR_gene, start_pos_plot))) %>%
  ungroup


  lead_3_4e11_matrixplotgg = lead_3_4e11_matrixplot %>% 
  arrange(desc(CIS_TRANS))
  lead_3_4e11_matrixplotgg2 = lead_3_4e11_matrixplotgg %>% 
  filter(!is.na(start_position))
  # mutate(CHR = ifelse(CHR %in% c(23,24) & CIS_TRANS == "zlimits",1,CHR))
lead_3_4e11_matrixplotgg2$CHR_gene = factor(lead_3_4e11_matrixplotgg2$CHR_gene, levels = c(23:1), labels = c("X",22:1))
lead_3_4e11_matrixplotgg2$CHR = as.integer(lead_3_4e11_matrixplotgg2$CHR)

######## FINAL PLOT ###########
lead_3_4e11_matrixplotgg2$P[lead_3_4e11_matrixplotgg2$P==0]=1e-300
pdf(paste0("/rhillary/dissemination2/olink_disc_matrixplot_2D_v2.pdf"),
    width = 8, height = 6)
cols <- c("red", "dark blue", NA)
p_gg = ggplot(data = lead_3_4e11_matrixplotgg2 %>% 
                mutate(abs_zscore = abs(beta/se)) %>% 
                mutate(logP2 = ifelse(logP > 20, 20, logP)),
              aes(x = pos_plot2, y = start_pos_plot2, 
                  colour = CIS_TRANS)) + 
  geom_point(aes(size = logP2), shape = 19, alpha = 0.8, size = 1) + 
  scale_size("p-value", range=c(0.1, 1), breaks=c(11, 13, 15, 20), labels = c("1e-11", "1e-13","1e-15", "<1e-20")) + 
  facet_grid(CHR_gene~CHR, scales = "free", space="free", switch=c("both")) + 
  theme_bw() + 
  labs(y= paste("Protein Position")) +
  labs(x= paste("vQTL Position"))  + 
  scale_color_npg(name=NULL, 
                     labels=c(bquote(italic(Cis)), bquote(italic(Trans)), "")) + 
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) + 
  theme(panel.spacing = unit(0, "lines")) + 
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.2),
        strip.text.y = element_text(angle=180)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(legend.position="right") + 
  theme(text = element_text(size = 12)) + theme(legend.text=element_text(size=11))
print(p_gg)
dev.off()





####################################
########## PLOT FIGURE 2B ##########
####################################

# Read in file and tabulate frequency of vQTL per Protein
a=as.data.frame(table(lead_3_4e11_matrixplotgg2$Protein))
a2=as.data.frame(table(a$Freq))
# Plot stage 
cbPalette <- c("#599ad3", "#f1595f", "#79c36a", "#9e66ab","#f9a65a","#D3D3D3")
pdf("/rhillary/dissemination2/prop_variants_per_protein.pdf",width=7.8,height=7.5)
bar = ggplot(data = a2, aes(a2$Var1, a2$Freq)) + geom_bar(stat="identity", fill = cbPalette) + xlab("Number of Independent Loci") + ylab("Number of Proteins")  + theme_scientific()
bar1=bar + geom_text(aes(label=a2$Freq), size=5, vjust=-0.3) + theme(axis.title.x = element_text(size=13.5), axis.title.y = element_text(size=13.5), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))
print(bar1)
dev.off()



####################################
########## PLOT FIGURE 2C ##########
####################################

# Format file for plot
plots=lead_3_4e11_matrixplotgg2
plots$freq1=log10(plots$freq)
plots$beta1=abs(plots$beta)
plots$beta1=log(plots$beta1)
cols <- c("red", "dark blue", NA)
plots$logP=-log10(plots$P)
plots$logP2 = ifelse(plots$logP > 20, 20, plots$logP) # in keeping with previous study 

# Plot stage 
pdf("/rhillary/dissemination2/mafs.pdf",width=8.8,height=7)
mafs = ggplot(data = plots, aes(plots$freq1, plots$beta1,color=factor(plots$CIS_TRANS))) + geom_point(aes(size = logP2),shape = 19, alpha = 0.8, size = 1) +
  scale_color_npg(name=NULL, 
                     labels=c(bquote(italic(Cis)), bquote(italic(Trans)), "")) + labs(y=expression(paste("log","[Beta]",sep=""))) + labs(x=expression(paste(log[10], "[MAF]",sep="")))
mafs1=mafs + theme_scientific() +
  theme(legend.position="right") + scale_x_continuous(limits=c(-1.30,-0.25),breaks=c(-1.25,-1,-0.75,-0.5,-0.25))  + scale_y_continuous(limits=c(-3,-0.5),breaks=c(-3,-2.5,-2,-1.5,-1,-0.5)) + 
  theme(axis.text = element_text(size = 11),axis.title=element_text(size=13)) + theme(legend.text=element_text(size=12.5))
mafs1+geom_smooth(method="lm",se=F)
dev.off()



##############################################################################################
######## FIGURE 2D - for convenience, plotted as part of Script_03 - vQTL replication ########
##############################################################################################


