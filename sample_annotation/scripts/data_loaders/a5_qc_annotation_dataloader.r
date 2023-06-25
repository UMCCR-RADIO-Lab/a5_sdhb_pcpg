library(googlesheets4)
library(dplyr)

#Function to read the clinical data from the google sheet or offline-cache
data_loader_a5_qc_anno <- function(google_account=NULL, use_cache=FALSE, offline_cache="/g/data/pq08/projects/ppgl/a5/offline_cache/a5_qc_annotation.tsv")
{
  if(exists("a5_qc")) {
    message("WARNING: Object a5_qc already exists in the global environment and will be overwritten")
  }
  
  if(use_cache)
  {
    message("Loading clinical data from offline cache...")
    a5_qc <- read.delim(offline_cache, sep="\t", header=T, check.names = F)
  } else {
    if(is.null(google_account)) { stop("You must provide a google credential with access to the clinical annotation sheet")}
    gs4_auth(google_account)
    message("Fetching annotation from google sheet ...")
    a5_qc <- read_sheet("1hnXdXI29KvvuLxsaTBID1-EbE7mTSk6bG05HcSFfhgo", col_types = "c", sheet = "QC_metrics")
    message("Updating offline cache ...")
    write.table(a5_qc, offline_cache, sep="\t", row.names = F)
  }
  

  assign(x = "a5_qc", value = a5_qc, envir = globalenv())
  message("Loaded QC annotation into the global environment as a5_qc")
}
message("Created data loader function data_loader_a5_qc_anno()")
