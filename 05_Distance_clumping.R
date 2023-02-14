###########################################
########### LD-based Clumping #############
###########################################

# Load requisite libraries 
library(pacman)
p_load(tidyverse,
       data.table, magrittr, tools, ggthemes,bitops)


# Read in Olink information
olink_protein_map = fread("/rhillary/olink_protein_map_1.5k_cleaned.tsv")
olink_protein_map$UKBPPP.ProteinID=gsub("\\:","_",olink_protein_map$UKBPPP_ProteinID)
# Fix proteins with names that alternate across files 
olink_protein_map[grep("LGALS7_LGALS7B_P47929_OID21406_v1",olink_protein_map$UKBPPP.ProteinID),"compare"]="LGALS7B_P47929_OID21406_v1_Oncology"
olink_protein_map[grep("MICB_MICA_Q29980_Q29983_OID20593_v1",olink_protein_map$UKBPPP.ProteinID),"compare"]="MICA_Q29980_OID20593_v1_Inflammation"


# Preprocess data
tmp=read_rds("/rhillary/outputs/discovery_p3.4e11.rds")
tmp$log10P=-log10(tmp$P)
tmp$compare=sub(".*?_", "", tmp$file)
tmp$compare=gsub(".phen","",tmp$compare)
#saveRDS(tmp,"/rhillary/outputs/discovery_p3.4e11_anno.rds")


# Identify all vQTLs at  Bonferroni threshold 
pQTLS_ALL_3_4e11_dis = read_rds("/rhillary/outputs/discovery_p3.4e11_anno.rds") %>% 
 filter(log10P > -log10(3.4e-11))
pQTLS_ALL_3_4e11_dis %<>% 
  mutate(HLA = case_when((Chr == 6 & bp > 25.5*1000000 & bp < 34.0*1000000) ~ 1, TRUE ~ 0)) %>% 
  mutate(position_start = ifelse(Chr == 6 & between(bp, 25500000,34000000), 25500000,bp - 1000000), 
         position_end = ifelse(Chr == 6 & between(bp, 25500000,34000000), 34000000,bp + 1000000)) %>% 
  mutate(position_start = ifelse(position_start < 0, 0, position_start)) %>% 
  arrange(desc(log10P))
pQTLS_ALL_3_4e11_dis_chromreg = pQTLS_ALL_3_4e11_dis %>% distinct(compare, Chr)


# Clumping stage 
clumped_regions = list()
 for(i in seq(nrow(pQTLS_ALL_3_4e11_dis_chromreg))){
   print(i)
   
   clumped_regions[[i]] = pQTLS_ALL_3_4e11_dis %>% 
     filter(compare %in% pQTLS_ALL_3_4e11_dis_chromreg$compare[i] & Chr %in% pQTLS_ALL_3_4e11_dis_chromreg$Chr[i]) %>% 
     GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T) %>% GenomicRanges::reduce() %>% 
     as_tibble() %>% select(Chr = seqnames,region_start = start, region_end = end) %>% 
     mutate(region_size = region_end-region_start) %>% 
     mutate(region_index = paste0(i,"_", 1:n())) %>% 
          mutate(compare = pQTLS_ALL_3_4e11_dis_chromreg$compare[i])
 }
 clumped_regions2 = bind_rows(clumped_regions)


# Select only lead signal(s) within predefined clumping window (1Mb here)
pQTLS_ALL_3_4e11_dis2 = pQTLS_ALL_3_4e11_dis %>% mutate(Chr = as.character(Chr)) %>% 
  left_join(clumped_regions2 %>% mutate(Chr = as.character(Chr)), by = c("compare", "Chr")) %>% 
  mutate(row_select = ifelse(bp > region_start & bp < region_end, 1, 0)) %>% 
  filter(row_select %in% 1) %>% 
  mutate(distance_to_midreg = abs(bp - (region_end + region_start)/2))

pQTLS_ALL_3_4e11_top_loci = pQTLS_ALL_3_4e11_dis2 %>% 
  group_by(compare, Chr, region_index) %>% 
  filter(log10P == max(log10P)) %>% 
  filter(NMISS == max(NMISS)) %>% 
  filter(distance_to_midreg == min(distance_to_midreg)) %>% 
  filter(bp == min(bp)) %>% ungroup()

#saveRDS(pQTLS_ALL_3_4e11_top_loci, "/rhillary/outputs/discovery_p3.4e11_clumped.rds")

  pQTLS_ALL_3_4e11_top_loci2 = pQTLS_ALL_3_4e11_top_loci %>% 
  left_join(olink_protein_map, by = "compare") %>% 
  mutate(cis_trans = ifelse(Chr == chr & (region_start < gene_end & region_end > gene_start), "cis", "trans")) %>% 
  select(Chr:row_select, Assay, Panel, UniProt2, HGNC.symbol, ensembl_id, cis_trans) %>% distinct() %>% 
  group_by(SNP, compare) %>% 
  filter(cis_trans == min(cis_trans)) %>% 
  mutate(UniProt2 = paste(unique(UniProt2), collapse = ","),
        HGNC.symbol =  paste(unique(HGNC.symbol), collapse = ","),
        ensembl_id =  paste(unique(ensembl_id), collapse = ",")) %>% 
  ungroup() %>% 
  distinct()

saveRDS(pQTLS_ALL_3_4e11_top_loci2, "/rhillary/outputs/discovery_p3.4e11_patched.rds")
