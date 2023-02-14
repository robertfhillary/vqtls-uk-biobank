######################################################
######## Main Effect Look-Up Analysis ################
######################################################

# Obtain meQTL results 
list_of_traits = readRDS("/UKB_PPP_pQTLs/analysis/interim/UKBPPP_GWAS_v1/discovery_rs2_paths_details_patched.rds")
# Obtain vQTLs for lookup 
# Read in discovery vQTLs 
res=readRDS("/rhillary/outputs/dis_independent_annotated.rds")
# Prepare phenotype names 
res$compare1=gsub("_", ":", res$compare)
res$compare1=sub(":[^:]+$", "",res$compare1)
# Fix protein names that may alternate between files 
res[grep("MICA_Q29980_OID20593_v1_Inflammation",res$compare),"compare1"]="MICB_MICA:Q29980_Q29983:OID20593:v1"
res[grep("LGALS7B_P47929_OID21406_v1_Oncology",res$compare),"compare1"]="LGALS7_LGALS7B:P47929:OID21406:v1"
res[grep("IL12B_P29459_OID21327_v1",res$compare),"compare1"]="IL12A_IL12B:P29459_P29460:OID21327:v1"
res[grep("CKMT1B_P12532_OID20721_v1",res$compare),"compare1"]="CKMT1A_CKMT1B:P12532:OID20721:v1"
res[grep("DEFA1B_P59665_OID20344_v1",res$compare),"compare1"]="DEFA1_DEFA1B:P59665:OID20344:v1"
res[grep("IL27_Q14213_OID21389_v1",res$compare),"compare1"]="EBI3_IL27:Q14213_Q8NEV9:OID21389:v1"
res[grep("FUT5_P21217_OID21013_v1",res$compare),"compare1"]="FUT3_FUT5:P21217_Q11128:OID21013:v1"
# Ensure entire set is present 
length(which(res$compare1 %in% list_of_traits$UKBPPP_ProteinID))
names(res)[1]="SNP"


# Loop through each GWAS file 
list3=list()
for(gw in unique(res$compare1)){
# Subset to GWAS file of interest 
tmp=list_of_traits[grep(gw,list_of_traits$UKBPPP_ProteinID),]
# Check 22 chromosomes are present 
print(nrow(tmp))
# Create list to store output from GWAS files (chromosome by chromosome)
list1=list()
# Loop through each chromosome 
for(i in seq(nrow(tmp))){
# Read in each chromosome and save out 
cwas_res_fil = fread(paste0(tmp$result_paths[i]))
cwas=cwas_res_fil[,c(1,2,3,6,7,8,10,11,13)]
list1[[i]]=cwas
}
# Combine chromosomes into one GWAS summary file for main effect QTLs 
list2=as.data.frame(do.call("rbind",list1))
# Subset our vQTL file to just the GWAS phenotype being queried 
tmp1=res[res$compare1 %in% gw,]
# Subset pQTL or main effect QTLs to same SNPs 
list2.1=list2[list2$ID %in% tmp1$SNP,]
# Tidy up pQTL file 
list2.1$P = 10^-list2.1$LOG10P
# Match order of SNPs and combine summ stats 
ids=tmp1$SNP
list2.1=list2.1[match(ids,list2.1$ID),]
tmp1$main_effect_P=list2.1$P
tmp1$main_effect_beta=list2.1$BETA
tmp1$main_effect_se=list2.1$SE
# Save out main effect QTLs + print to denote completion
list3[[gw]]=tmp1
print(gw)
} 
# Collate final set 
vqtl=as.data.frame(do.call("rbind",list3))