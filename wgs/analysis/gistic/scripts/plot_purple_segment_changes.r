# Make a heatmap of segments from PURPLE
library(tidyverse)
library(RCircos)

# This script plots copy number changes as determined by PURPLE across the A5 cohort

blank_theme <- theme_bw()+
  theme(panel.grid=element_blank(),
        strip.background = element_blank())

# Samples that were either wrong or had the wrong genotype
to_remove <- c("E124-1", "E145-1","E191-1", "E181-1", "E167-1", "E154-1")

# Read in purple outputs
purple_outputs <- read_csv("~/Documents/Projects/Year_2018/A5study/WGS/Bioinformatics/Processed_data/Purple/All_purple_cnv_combined.csv") %>%
  # Remove redundannt samples 
  filter(! sample %in% to_remove)

# Back to processing for the coverage plot 
plot.data_2 <- purple_outputs %>%
  mutate(sample_number_start = as.numeric(as.factor(purple_outputs$sample)))%>%
  mutate(sample_number_end = as.numeric(as.factor(purple_outputs$sample)) +1 )%>%
  mutate(copyNumber = as.numeric(copyNumber))%>%
  # Set negative CNs to 0.1
  mutate(copyNumber = replace(copyNumber,copyNumber <= 0.1,0.1))%>%
  # Set a copy number of 2 to 0
  mutate(copyNumber = log2(copyNumber)-1)%>%
  dplyr::select(sample,sample_number_start,sample_number_end, chrom = X.chromosome,region_start = start, region_end = end, copyNumber)%>%
  filter(chrom != "Y")%>%
  mutate(chrom = paste0("chr",chrom))%>%
  group_by(chrom)%>%
  # Get the max chromosome sizes
  mutate(max_end = max(region_end))%>%
  ungroup()

# Correct ordering of chromosomes 
chr_order <- paste0("chr",c(as.character(1:22), "X"))

# Get the max value from the data to dtermine chromosome ends
overall_chr_sizes <- plot.data_2 %>%
  group_by(chrom)%>%
  dplyr::summarise(length = max(region_end)- min(region_end))%>%
  mutate(chrom = factor(chrom, levels = chr_order))%>%
  arrange(chrom)%>%
  mutate(length= as.numeric(length))%>%
  # The max for the x axis seems to be way off so reducing to 
  mutate(length = replace(length, chrom == "chrX", 51304565))%>%
  mutate(length = replace(length, chrom == "chr20", 59128982))%>%
  # Cumulative sum of chromosome length
  mutate(cumu_sum_chr = cumsum(length))


# Subtract chr1 from the cumulative sums to make the plot start at 0
chr_join <- dplyr::select(overall_chr_sizes,chrom,cumu_sum_chr)%>%
  mutate(cumu_sum_chr = cumu_sum_chr - cumu_sum_chr[[1]])%>%          
  mutate(chrom = as.character(chrom))

# Join the plot data and the info about where to move the chrs to plot them
plot.data_3 <- left_join(plot.data_2, chr_join)%>%
  mutate(region_start_plot = region_start + cumu_sum_chr)%>%
  mutate(region_end_plot = region_end + cumu_sum_chr)%>%
  # Make sure chrs have correct order
  mutate(chrom = factor(chrom, levels = chr_order))
# Read in the clinical data output to plot over the heatmap ----


samps <- unique(plot.data_3$sample)  
  
y_axis <- data.frame(sample = c(samps,  "Gistic broad regions", "Gistic focal regions"), 
                     pos = 1:(length(samps)+2)+0.5)

# Move up the gistic labels
y_axis$pos[length(y_axis$pos)] <- y_axis$pos[length(y_axis$pos)] +0.5

y_axis$pos[length(y_axis$pos)-1] <- y_axis$pos[length(y_axis$pos)-1] +3.5

# Add on an endpoint value that is not present in the cumsums
break_mid <- c(unique(plot.data_3$cumu_sum_chr), max(plot.data_3$region_end_plot))

# Get the midpoints between breaks for the x axis
midpoints <- numeric()
for (i in 1:length(break_mid)){
  if(i == length(break_mid)){
    next
  }
  midpont <- median(c(break_mid[[i]], break_mid[[i+1]]))
  midpoints <- c(midpoints, midpont)
}

x_axis <- data.frame(chrs = chr_join$chrom, midpoints, stringsAsFactors = F)

# Read in the cyto bands for hg-19
cyto_bands_hg19 <- read_tsv("/Users/adpattison/Documents/Projects/Year_2018/A5study/WGS/Bioinformatics/Raw_data/Genome_annotations/cytoBand_hg_19.txt",col_names = F)%>%
  dplyr::rename(chrom = "X1")

# Joined to the cumsum chr lengths to get regions for plotting
joined_cyto <- inner_join(cyto_bands_hg19,chr_join)%>%
  mutate(cyto_start = X2 + cumu_sum_chr, cyto_end = X3 + cumu_sum_chr)

cyto_summary <- joined_cyto%>%
  group_by(chrom)%>%
  summarise(min = min(cyto_start),max = max(cyto_end))

cyto_lines <- data.frame(lines = c(cyto_summary$min, max(cyto_summary$max)))

# Make data frames to plot p and q as separate layers  
p_frame <- filter(joined_cyto, grepl("p", X4))%>%
  # Add a 'join column to join arms with the gisitc output. 
  mutate(Arm = paste0(chrom, substr(X4, 1, 1)))

q_frame <- filter(joined_cyto, grepl("q", X4))%>%
  # Add a 'join column to join arms with the gisitc output. 
  mutate(Arm = paste0(chrom, substr(X4, 1, 1)))

# Read in gistic broad chnages to get significance for cytoband loss/gain

broad_signif_amp <- read_tsv("/Volumes/mdhs-clinical/5100/UMCCR/Lab-Tothill/Restricted/Projects/PPGL/A5/Bioinformatics/WGS/Bioinformatics_tools_outputs/GISTIC/A5_non_redundant_set/broad_significance_results.txt")%>%
  dplyr::select(Arm,`Amp frequency`, `Amp q-value`)%>%
  filter(`Amp q-value` <= 0.05)%>%
  mutate(Type = "Amp")

colnames(broad_signif_amp) <- c('Arm', 'Freq' , 'FDR', 'Type')

broad_signif_del <- read_tsv("/Volumes/mdhs-clinical/5100/UMCCR/Lab-Tothill/Restricted/Projects/PPGL/A5/Bioinformatics/WGS/Bioinformatics_tools_outputs/GISTIC/A5_non_redundant_set/broad_significance_results.txt")%>%
  dplyr::select(Arm, `Del frequency`, `Del q-value`)%>%
  filter(`Del q-value` <= 0.05)%>%
  mutate(Type = "Del")

colnames(broad_signif_del) <- c('Arm', 'Freq' , 'FDR', 'Type')

combined_arm <- bind_rows(list(broad_signif_amp,broad_signif_del))%>%
  mutate(Type = as.numeric(gsub(" .*", "", `Type`) == "Amp"))%>%
  mutate(Type = replace(Type,Type ==0 , -1))%>%
  mutate(Arm = paste0("chr", Arm))

p_gisitc <- inner_join(p_frame, combined_arm)

# 1q had both signifcant AMP and DEl so I went with the AMP as it was more common 
q_gisitc <- inner_join(q_frame, combined_arm)

filter <- q_gisitc$Arm == "chr1q" & q_gisitc$Type == -1

q_gisitc <- q_gisitc[!filter,]

# Read in the gistic output to plot over the heatmap ----
processed_gistic_info <- read_csv("~/Documents/Projects/Year_2018/A5study/WGS/Bioinformatics/Processed_data/GISTIC_outputs/Processed_gisitc_oututs/gisitic_positional_info.csv")

important_gisitc_cols <- processed_gistic_info %>%
  select(`Unique Name`, Descriptor, `Region Limits`, `q values`, chrom =chr, start, end)%>%
  mutate(chrom = paste0("chr",chrom))

# Plot gistic values at the right point
joined_gistic <- inner_join (important_gisitc_cols, chr_join)%>%
  mutate(start_plot = start + cumu_sum_chr, end_plot = end + cumu_sum_chr)%>%
  mutate(amp_del = as.numeric(gsub(" .*", "", `Unique Name`) == "Amplification"))%>%
  mutate(amp_del = replace(amp_del,amp_del ==0 , -1))

gisitic_y_axis <- length(unique(plot.data_3$sample))
  
aberrations_plot <- ggplot(plot.data_3)

aberrations_plot + 
  geom_rect(aes(xmin=region_start_plot, xmax=region_end_plot, ymin=sample_number_start, ymax=sample_number_end, fill=copyNumber))+ 
  theme_bw() +
  scale_fill_gradient2(low = "darkblue", high = "darkred",midpoint = 0)+
  geom_vline(data = cyto_lines, aes(xintercept = lines, alpha= 0.1),linetype="dotted")+
  labs(x = 'Chromosome',y = 'Sample', fill = "Total copy\nnumber")+
  scale_y_continuous(breaks= y_axis$pos,labels = y_axis$sample)+
 geom_rect(data = joined_gistic, aes(xmin=start_plot, xmax=end_plot, ymin=gisitic_y_axis+2, ymax=gisitic_y_axis+4, fill = amp_del))+
  scale_x_continuous(breaks= x_axis$midpoints, labels= x_axis$chrs)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size=20))+
  geom_rect(data = q_frame, aes(xmin=cyto_start, xmax=cyto_end, ymin=0, ymax=1), fill = "purple")+
  geom_rect(data = p_frame, aes(xmin=cyto_start, xmax=cyto_end, ymin=0, ymax=1), fill = "grey")+
  geom_rect(data = q_gisitc, aes(xmin=cyto_start, xmax=cyto_end, ymin=gisitic_y_axis+4, ymax=gisitic_y_axis+6,fill = Type , alpha = Freq))+
  geom_rect(data = p_gisitc, aes(xmin=cyto_start, xmax=cyto_end, ymin=gisitic_y_axis+4, ymax=gisitic_y_axis+6,fill = Type ,alpha = Freq))+
  guides(alpha = F)+
  ggsave("/Users/adpattison/Documents/Projects/Year_2018/A5study/Results/Plots/CNV/A5_all_CNV.pdf",
         width = 30, height =24)

# Make a plot to go next to the CNV plot that shows TERT/ATRX status

running_summary <- read_csv("~/Documents/Projects/Year_2018/A5study/Clinical_data/A5_key_stats_summary.csv")%>%
  select(A5_ID, `Metastasis driver`)


ordering <- data.frame (A5_ID = unique(plot.data_3$sample))%>%
  left_join(running_summary)%>%
  mutate(`Metastasis driver` = replace(`Metastasis driver`, is.na(`Metastasis driver`), 'Unknown'))%>%
  mutate(A5_ID = factor(A5_ID, levels = unique(A5_ID)))

ggplot(data = ordering , aes(x = A5_ID, y = 1, fill = `Metastasis driver`)) + 
  geom_tile(width=0.7, height=0.7,size = 0.7)+
  coord_flip()+  
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.margin = unit(c(1,1,1,1), "cm"))+
  labs(x = "Sample (Non-redundant)", y = "")+
  scale_fill_manual(values = c ("orange","red", "white"))+
  blank_theme+
  ggsave("/Users/adpattison/Documents/Projects/Year_2018/A5study/Results/Plots/CNV/A5_all_CNV_annotation.pdf")
  









