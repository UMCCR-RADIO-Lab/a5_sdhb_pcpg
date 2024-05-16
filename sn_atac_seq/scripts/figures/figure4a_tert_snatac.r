library(Seurat)
library(Signac)

###############
# Data loader #
###############

source("a5/sn_atac_seq/scripts/data_loaders/a5_snatac_dataloader.r")
data_loader_a5_snatac(quickload = T)

#########
# Paths #
#########

MACS_out_dir="/g/data/pq08/projects/ppgl/a5/sn_atac_seq/analysis/macs2_via_signac"
gencode_gtf="/g/data/pq08/reference/GRCh38/gencode/gencode.v36.primary_assembly.annotation.gtf"
spark_script="/g/data/pq08/projects/ppgl/a5/software/spark/spark.py"
spark_out="/g/data/pq08/projects/ppgl/a5/sn_atac_seq/results/figures/tert_region_mac2_fragments_spark2"

########
# MACS #
########

peaks <- CallPeaks(
  object = a5_snatac,
  assay = "cellranger",
  group.by = "cluster",
  cleanup = FALSE,
  macs2.path = "/g/data/pq08/projects/ppgl/a5/software/macs2/bin/macs2",
  additional.args = "--bdg",
  outdir = MACS_out_dir
)

#########
# Spark #
#########

samples <- c("E156-1","E188-1","E201-1","E200-1","E123-1","E197-1","E166-1","Chromaffin_cells","Adrenocortical_cells","Endothelial_cells","Fibroblasts")

run_script = file.path(dirname(spark_script),"run_spark.sh")

script_lines <- c(
  "module load htslib/1.16",
  paste("MACS_out_dir",MACS_out_dir,sep="="),
  paste("gencode_gtf", gencode_gtf,sep="="),
  paste("spark_script",spark_script,sep="="),
  paste("tabix -p bed", paste0("$MACS_out_dir/", samples, "_treat_pileup.bdg.gz")),
  "python3 ${spark_script} \\",
  "-pr chr5:1242692-1305544 \\",
  paste("-cf", paste(paste0("$MACS_out_dir/", samples, "_treat_pileup.bdg.gz"), collapse=" "), "\\"),
  "-gtf ${gencode_gtf} \\",
  paste("-gl", paste(samples, collapse=" "), "\\"),
  "-dg TERT \\",
  "-gs yes \\",
  paste("--output_name", spark_out)
)

readr::write_lines(script_lines, run_script)
system2(command = "bash", args = run_script)
file.remove(run_script)
