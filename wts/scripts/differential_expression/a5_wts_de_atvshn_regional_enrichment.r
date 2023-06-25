##########################################
# This script generates plots            #
# focused on co-located clusters of      #
# genes found to be DE between           #
# abdo-thoracic and head and neck PPGL   #
#                                        #
# Author: Aidan Flynn                    #
# Date: 07/06/2023                       #
#                                        #
##########################################

library(ggplot2)
library(ggrepel)
library(patchwork)

setwd("/g/data/pq08/projects/ppgl")


#################
# Run DE script #
#################

#Performs DE
# creates globals:
# - EnsIds
# - data_loader_ensgid_to_chr()
# - wts_top_tables[]
source("./a5/wts/scripts/differential_expression/a5_wts_differential_expression.r")

######################
# Load cytoband data #
######################

#Cytobands - HG38
if(!file.exists("/g/data/pq08/reference/GRCh38/cytobands_hg38.txt.gz")) {
  download.file(url = "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz",
                destfile ="/g/data/pq08/reference/GRCh38/cytobands_hg38.txt.gz")
}

cytoband <- read.delim("/g/data/pq08/reference/GRCh38/cytobands_hg38.txt.gz", header=F)
colnames(cytoband) <- c("chr","start","end","cytoband","stain")
cytoband <- cytoband %>% filter(chr %in% paste0("chr", c(1:22,"X","Y")))

cytoband_gr <- GenomicRanges::makeGRangesFromDataFrame(cytoband, 
                                                       seqnames.field ="chr", 
                                                       ignore.strand = T, 
                                                       keep.extra.columns = T)

###########################################
# Annotate DE TopTable with cytoband data #
###########################################

top_table <- wts_top_tables[["Parasympathetic_vs_Sympathetic"]] 

top_table <- top_table %>% 
  separate(Gene, into=c("ensembl_gene_id", "gene_symbol"), 
           sep="_", 
           extra = "merge", remove = F) %>% 
  mutate(ensembl_gene_id=gsub("[.][0-9]+$", "", ensembl_gene_id))

data_loader_ensgid_to_chr(EnsIds = EnsIds, 
                          use_cache = T, 
                          update_cache = F, 
                          ensembl_mirror = "https://www.ensembl.org")

ensgid_to_chr <- ensgid_to_chr %>% 
  filter(!is.na(chromosome_name)) %>% 
  mutate(chromosome_name=paste0("chr",chromosome_name))

ensgid_to_chr_gr <-   GenomicRanges::makeGRangesFromDataFrame(ensgid_to_chr, 
                                                          seqnames.field ="chromosome_name", 
                                                          start.field = "start_position",
                                                          end.field = "start_position",
                                                          ignore.strand = T,keep.extra.columns = T)

hits <- findOverlaps(subject = ensgid_to_chr_gr, query = cytoband_gr)
mcols(hits)[["cytoband"]] <- mcols(cytoband_gr)[["cytoband"]][queryHits(hits)]
mcols(hits)[["cytoband_start_pos"]] <- start(cytoband_gr)[queryHits(hits)]
mcols(hits)[["cytoband_end_pos"]] <- end(cytoband_gr)[queryHits(hits)]
mcols(hits)[["arm"]] <- stringr::str_extract(string = mcols(hits)[["cytoband"]], pattern = "p|q")

ensgid_to_chr[subjectHits(hits),"cytoband"] <- mcols(hits)[["cytoband"]]

top_table <- top_table %>% inner_join(ensgid_to_chr) 

top_table <- top_table %>% 
  mutate(chromosome_name = factor(chromosome_name, 
                                  levels=paste0("chr", c(1:22,"X","Y")))) %>% 
  filter(!is.na(chromosome_name))

###################################
# Identify significant gene peaks #
###################################

plot.data <- top_table %>% group_by(chromosome_name, cytoband) %>% 
  dplyr::summarise(start=min(start_position), nGenes_cytoband=n()) %>% 
  left_join(top_table %>% 
              filter(abs(logFC) > 1, adj.P.Val < 0.05) %>% 
              group_by(chromosome_name, cytoband) %>% 
              dplyr::count(name="nGenes_DE"), by=c("chromosome_name", "cytoband")) %>% 
  mutate(nGenes_DE=replace_na(nGenes_DE,0)) %>% 
  arrange(chromosome_name, start) %>% 
  mutate(region=factor(paste0(chromosome_name, cytoband))) %>% 
  mutate(chromosome_name=factor(chromosome_name, levels=paste0("chr", c(1:22,"X","Y")))) %>% 
  mutate(ratio=nGenes_DE/nGenes_cytoband) %>% 
  filter(!is.na(chromosome_name), chromosome_name!="chrY")


ggplot(plot.data, aes(x=region, y=ratio)) + 
  geom_col() + 
  facet_wrap("chromosome_name", scales = "free_x") + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.5,hjust=1)) 


ggplot(top_table, aes(x=start_position, y=-log10(adj.P.Val), color=logFC)) + 
  geom_point(size=1) + 
  geom_hline(yintercept = -log10(0.05), linetype=2, color="red") +
  facet_wrap("chromosome_name", scales="free_x") +
  scale_color_gradientn(colours = c("blue", "green", "grey", "orange", "red"), 
                        values = c(0, 0.25, 0.5, 0.8, 1)) 


##################################
# Extract significant gene peaks #
##################################

peak_genes <- top_table %>%  filter(abs(logFC) > 1, adj.P.Val < 0.05) %>% 
  filter(
    (chromosome_name=="chr7"  & start_position > 2.65*10^7  & start_position < 2.72*10^7) |
      (chromosome_name=="chr12"  & start_position > 5.38*10^7 & start_position < 5.41*10^7) | 
      (chromosome_name=="chr17" & start_position > 3.355*10^7 & start_position < 3.5*10^7) | 
      (chromosome_name=="chr17" & start_position > 4.82*10^7 & start_position < 4.9*10^7) |
      (chromosome_name=="chr20" & start_position > 6.1*10^7 & start_position < 6.2*10^7))  %>% 
  dplyr::select(Gene, gene_symbol, gene_biotype, chromosome_name, start_position, end_position, cytoband)


#####################
# Prepare plot data #
#####################

########
# Join per sample expression 
########

plot_data <- peak_genes %>% inner_join(a5_wts_lcpm_list[["SDHB"]] %>% 
  as_tibble(rownames = "ensgid_symbol") %>% 
  pivot_longer(cols=-ensgid_symbol, 
               names_to = "A5_ID", 
               values_to = "log2cpm"), 
  by=c("Gene"="ensgid_symbol"))

########
# Join clinical annotation
########

plot_data <- plot_data %>% inner_join(a5_anno %>% dplyr::select(A5_ID, differential_group_anatomy))


########
# Z-scale expression
########

plot_data <- plot_data %>%  group_by(gene_symbol) %>% mutate(log2cpm_z= (log2cpm-mean(log2cpm))/sd(log2cpm))

########
# Replace contig names with Ensembl gene IDs for genes for un-annotated genes
########

plot_data <- plot_data %>% 
  mutate(gene_symbol=ifelse(grepl("^AC[0-9]", gene_symbol), 
                            stringr::str_extract(string = Gene, pattern = "ENSG[0-9]+"), 
                            gene_symbol))

########
# Factorise gene symbols in chromosome positional order
########

plot_data <- plot_data %>% 
  arrange(chromosome_name, start_position) %>% 
  mutate(gene_symbol=factor(gene_symbol, levels=unique(.$gene_symbol)))

########
# Create cytoband label
########

plot_data <- plot_data %>% mutate(cytoband=paste0(chromosome_name,cytoband))

########
# Make Plot friendly names for clinical values
########

plot_data <- plot_data %>% 
  mutate(`Tumour Location`=dplyr::recode(differential_group_anatomy, 
                                         "Abdominal_Thoracic"="Abdominal/Thoracic",
                                         "Head_Neck"="Head and neck"))

##################
# Generate Plots #
##################

########
# Gene expression by anatomical type for each region
########

gg_expr <- list()
for (current_cytoband in c("chr7p15.2", "chr12q13.13", "chr17q12", "chr17q21.32", "chr20q13.33"))
{
  current_plot_data <- plot_data  %>% 
    filter(cytoband == current_cytoband, 
           `Tumour Location` %in% c("Abdominal/Thoracic","Head and neck")) %>% 
    arrange(start_position) %>% 
    mutate(gene_symbol=factor(gene_symbol, levels=unique(.$gene_symbol)))
  
  gg_expr[[current_cytoband]] <- 
    ggplot(data=current_plot_data, 
           mapping=aes(x=gene_symbol, y=log2cpm_z, 
                       # fill=`Tumour Location`,
                       color=`Tumour Location`)) + 
    geom_boxplot(outlier.size = 0.5) +
    geom_vline(xintercept = seq(1.5,length(levels(current_plot_data$gene_symbol)),1),
               linetype=2, color="grey") +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5,hjust=1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) + 
    # scale_fill_manual(values=c(`Abdominal/Thoracic`=location_cols[["Head_neck"]],
    #                            `Head and neck`=location_cols[["Extraadrenal"]])) +
    scale_color_manual(values=c(`Abdominal/Thoracic`=location_cols[["Head_neck"]],
                               `Head and neck`=location_cols[["Extraadrenal"]])) +
    ylab(bquote('CPM ('~log[2]~' z-scaled)')) +
    facet_wrap("cytoband", scales="free_x", ncol=1) + 
    xlab("") 

}


########
# Manhattan style plot (chr. position by -log10(p)) from limma top-table
########

gg_manhattan <- ggplot(top_table %>% filter(chromosome_name %in% c("chr7","chr12","chr17")), aes(x=start_position, y=-log10(adj.P.Val), color=logFC)) + 
  geom_point(size=1) + 
  geom_hline(yintercept = -log10(0.05), linetype=2, color="red") +
  facet_wrap("chromosome_name", scales="free_x") +
  scale_x_continuous(labels = \(x) paste(x/10^6, "Mb")) +
  scale_color_gradientn(colours = c("blue", "green", "grey", "orange", "red"), 
                        values = c(0, 0.25, 0.5, 0.8, 1),
                        name= bquote(Log[2]~" FC")) +
  ylab(bquote('-'~log[10]~'(adjusted p-value)')) + 
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust=1,
                                   hjust=1),
        axis.title.x = element_blank()) 
  
################
# Format plots #
################

layout_design <- "AABBB\nCCCCC\nDDDEE"

gg_expr[["chr7p15.2"]] +  #AA
  gg_expr[["chr12q13.13"]] + ylab("") + #BBB
  gg_manhattan + #CCCC
  gg_expr[["chr17q12]"]] + #DDD
  gg_expr[["chr17q21.32"]]  +  ylab("") + #EE
  plot_layout(guides="collect", 
              design = layout_design, 
              heights = c(1,1.5,1)) 



