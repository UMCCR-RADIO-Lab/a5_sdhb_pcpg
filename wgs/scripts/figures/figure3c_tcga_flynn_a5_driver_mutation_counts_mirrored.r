library(ggplot2)
library(patchwork)
################
# Data Loaders #
################

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)

source("/g/data/pq08/projects/ppgl/a5/wgs/scripts/data_loaders/wgs_dataloaders.r")
data_loader_somatic_variants(quickload = T)
data_loader_sv_gridsslinx(quickload = T)

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

###############
# Public data #
###############

tcga_mutation <- readr::read_tsv("/g/data/pq08/projects/ppgl/public_data/wgs_wes/tcga/data_mutations.txt")
tcga_cgi <- readr::read_tsv("/g/data/pq08/projects/ppgl/public_data/wgs_wes/tcga/tcga_pcpg_cgi_ouput.tsv")

tcga_mutation <- tcga_mutation %>% 
  left_join(tcga_cgi %>%  
              mutate(Chromosome=paste0("chr", CHROMOSOME)) %>% 
              dplyr::select(Chromosome, Start_Position=pos, `CGI-Oncogenic Prediction`), relationship = "one-to-one" )

flynn_mutation <- readr::read_tsv("/g/data/pq08/projects/ppgl/public_data/wgs_wes/flynn/data_mutations.txt")
flynn_cgi <- readr::read_tsv("/g/data/pq08/projects/ppgl/public_data/wgs_wes/tcga/tcga_pcpg_cgi_ouput.tsv")

flynn_mutation <- flynn_mutation %>% 
  left_join(flynn_cgi %>%  
              mutate(Chrom=paste0("chr", CHROMOSOME)) %>% 
              dplyr::select(Chrom, Position=pos, `CGI-Oncogenic Prediction`), relationship = "one-to-one" )

tcga_flynn <- bind_rows(tcga_mutation %>%  mutate(Sample=stringr::str_sub(Tumor_Sample_Barcode, start = 1, end = 12)) %>%  dplyr::select(Sample, Hugo_Symbol, Consequence), 
                        flynn_mutation %>%  dplyr::select(Sample=SampleName, Hugo_Symbol=GeneHGNC, Consequence))

################
# Driver Genes #
################

cgi_drivers <- bind_rows(a5_somatic_variants_keep %>% 
                           filter(A5_ID %in% a5_anno$A5_ID, A5_ID != "E167-1", Tumour_AF > 0.1)  %>% 
                           ungroup() %>%  
                           dplyr::select(Hugo_Symbol=PCGR_SYMBOL, `CGI-Oncogenic Prediction`=CGI.Oncogenic.Prediction),
                         tcga_cgi %>%  dplyr::select(Hugo_Symbol=`CGI-Gene`, `CGI-Oncogenic Prediction`),
                         flynn_cgi %>%  dplyr::select(Hugo_Symbol=`CGI-Gene`, `CGI-Oncogenic Prediction`)
                         ) %>% 
  filter(`CGI-Oncogenic Prediction` %in% c("driver (boostDM: non-tissue-specific model)", "driver (oncodriveMUT)", "driver (oncodriveMUT)/passenger (oncodriveMUT)"))

cgi_drivers <- c(unique(cgi_drivers$Hugo_Symbol), "TERT")

cgi_drivers <- intersect(cgi_drivers, 
                         (a5_somatic_variants_keep %>% 
                           filter(A5_ID %in% a5_anno$A5_ID, A5_ID != "E167-1", Tumour_AF > 0.1)  %>% 
                            pull(PCGR_SYMBOL)))

##########
# Filter #
##########

a5_somatic_variants_keep_drivers <- a5_somatic_variants_keep %>% 
  filter(A5_ID %in% a5_anno$A5_ID, A5_ID != "E167-1", Tumour_AF > 0.1) %>%
  mutate(SampleName=gsub("-.$", "", A5_ID)) %>% 
  ungroup() %>% 
  filter(PCGR_SYMBOL %in% cgi_drivers) %>% 
  dplyr::select(SampleName, Hugo_Symbol=PCGR_SYMBOL, Consequence) %>%  distinct()

A5_gridss_keep_drivers <- bind_rows(A5_gridss_keep %>% dplyr::select(A5_ID, Hugo_Symbol=GeneStartName, GeneDisrupted=GeneStartDisrupted),
          A5_gridss_keep %>% dplyr::select(A5_ID, Hugo_Symbol=GeneStartName, GeneDisrupted=GeneStartDisrupted)) %>% 
  filter((Hugo_Symbol %in% cgi_drivers & GeneDisrupted) | (Hugo_Symbol %in% c("ATRX","TERT"))) %>% distinct() %>% mutate(Consequence = "Structural variant")  %>%
  mutate(SampleName=gsub("-T0.$", "", A5_ID))

tcga_flynn_drivers <- tcga_flynn %>%  filter(Hugo_Symbol %in% cgi_drivers)


#########
# Count #
#########

plot_data <- bind_rows(a5_somatic_variants_keep_drivers %>% mutate(source="A5", total_samples=length(unique(a5_anno$`Patient ID`))),
                       A5_gridss_keep_drivers %>% mutate(source="A5", total_samples=length(unique(a5_anno$`Patient ID`))),
                       tcga_flynn_drivers %>% mutate(source="TCGA/Flynn", total_samples = length(unique(tcga_flynn$Sample)))) %>% 
  group_by(source, total_samples, Hugo_Symbol, Consequence) %>% 
  dplyr::count() %>% 
  mutate(pcnt=(n/total_samples)*100) %>% 
    pivot_wider(id_cols = c(Hugo_Symbol,Consequence), names_from = source, values_from = c(pcnt,n)) %>%  
  mutate(across(.cols=c(pcnt_A5, `pcnt_TCGA/Flynn`,  n_A5, `n_TCGA/Flynn`), .fns=~replace_na(.x,0))) %>% 
  group_by(Hugo_Symbol) %>% 
  mutate(total_pcnt_a5=sum(pcnt_A5)) %>% 
  arrange(desc(total_pcnt_a5)) %>% 
  filter(pcnt_A5 > 0)  %>% 
  mutate(Hugo_Symbol=factor(Hugo_Symbol, levels=unique(.$Hugo_Symbol)))
  
plot_data <- plot_data %>% mutate(Consequence=case_when( 
  grepl("missense_variant", Consequence) ~ "Missense", 
  Consequence == "frameshift_variant" ~ "Frameshift", 
  Consequence == "stop_gained" ~ "Stop gained",
  grepl("splice_acceptor_variant", Consequence) ~ "Splice acceptor",
  grepl("splice_donor_variant", Consequence) ~ "Splice donor",
  grepl("splice_region_variant", Consequence) ~ "Splice region",
  Consequence == "inframe_deletion" ~ "In-frame deletion",
  Consequence == "upstream_gene_variant" ~ "Promoter mutation",
  .default = Consequence
))
  
############
# Plotting #
############

ggplot(plot_data, aes(x=Hugo_Symbol, y=n_A5, fill=Consequence)) + geom_col() + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(6,6,0,6)) +
  scale_fill_manual(values=genomic_alteration_cols) +
ggplot(plot_data, aes(x=Hugo_Symbol, y=`n_TCGA/Flynn`, fill=Consequence)) + 
  geom_col() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=45, vjust=0.5,hjust=1),
        plot.margin = margin(0,6,6,6)) + 
  scale_fill_manual(values=genomic_alteration_cols) +
  scale_y_reverse() + 
  plot_layout(nrow=2, guides="collect")
