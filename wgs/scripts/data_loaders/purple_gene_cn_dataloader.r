library(dplyr)
library(tidyr)

fetchPurpleGeneCNFileNames <- function(purple_cn_file_dir)
{
  A5_CN.files <- list.files(purple_cn_file_dir, pattern = ".purple.cnv.gene.tsv", full.names = T, recursive = T)
  names(A5_CN.files) <- gsub("__.+$","",basename(A5_CN.files))
  return(A5_CN.files)
}

readPurpleGeneCN <- function(FileNames, Combine=T)
{
  if (Combine)
  {
    return(bind_rows(lapply(FileNames, read.delim), .id="A5_ID"))
  } else {
    return(lapply(FileNames, read.delim))
  }
}


