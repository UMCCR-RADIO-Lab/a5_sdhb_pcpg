## TCGA/Comete/A5 SuperSet

#Run 
#"/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Code/Combine TCGA Comete A5 small RNA-Seq.R"
# Lines 1 - 146

annotation.smallrna <- annotation
batch_removed.smallrna <- batch_removed


#Run 
#"/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Code/Combine PPGL methylation values.R"
# Lines 1 - 159

annotation.methylation <- annotation
batch_removed.methylation <- batch_removed

#Run 
#"/data/cephfs/punim0010/projects/flynna/A5//WTS/tcga_a5_integration/tcga_a5_integration.r"
# Lines 1 - 113

annotation.rnaseq <- annotation
batch_removed.rnaseq <- batch_removed

## Homogenise TCGA Labelling
colnames(batch_removed.rnaseq) <- gsub("(TCGA-..-....-...-..).+","\\1", colnames(batch_removed.rnaseq))
colnames(batch_removed.smallrna) <- gsub("(TCGA-..-....-...-..).+","\\1", colnames(batch_removed.smallrna))

complete_samples <- intersect(intersect(colnames(batch_removed.smallrna), colnames(batch_removed.rnaseq)), colnames(batch_removed.methylation))

batch_removed.rnaseq <- batch_removed.rnaseq[,complete_samples]
batch_removed.smallrna <- batch_removed.smallrna[,complete_samples]
batch_removed.methylation <- batch_removed.methylation[,complete_samples]

annotation.superset <- bind_rows(annotation.methylation, annotation.rnaseq, annotation.smallrna) %>% 
  dplyr::select(Barcode, Cluster, Dataset, new_naming) %>% distinct() %>% 
  filter(Barcode %in% complete_samples) %>% arrange(Barcode)

#Combine datasets
superset <- bind_rows(data.frame(batch_removed.rnaseq, check.names = F), 
                      data.frame(batch_removed.smallrna, check.names = F), 
                      data.frame(batch_removed.methylation, check.names = F))

#Remove features with missing values
superset <- superset[!apply(superset, 1, anyNA),]

set.seed(42)

umap_config <- umap.defaults
umap_config$n_neighbors=10
umap_config$spread=3

# Plot a UMAP as well
umap <- umap(t(superset), config = umap_config)

to_plot_umap <- data.frame(umap$layout)%>%
  rownames_to_column("Barcode")%>%
  left_join(annotation.superset)

write_delim(to_plot_umap, 
            file =  "/data/cephfs/punim0010/projects/flynna/A5/tcga_a5_superset/plots/umap_coord_superset_tcga+a5_allsamples_nn10_seed42_spread3.txt", 
            delim = "\t")

ggplot(data = to_plot_umap %>% mutate(Label=ifelse(Dataset %in% c("A5"),Barcode,NA)), aes(x = X1, y = X2, colour = new_naming,shape= Dataset, label = Label)) + #
  geom_point()+
  geom_text_repel(color ="black",nudge_y = 0.5)+
  blank_theme+
  guides(label= F)+
  labs(shape ="", x= "UMAP 1", y = "UMAP 2", colour = "Subtype")+
  scale_color_manual(values = subtpye_cols)+
  theme(aspect.ratio=1)
