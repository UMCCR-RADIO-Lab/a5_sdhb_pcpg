#Helper function to give classifications to CNA events based on various metrics
#Classes of CNA events

cn_event_types <- c("Diploid/Haploid-X", "Loss",  "Subclonal Loss", "Minor Subclonal Loss",  
                    "Loss + Subclonal CNLOH", "CNLOH", "Gain", "Subclonal Gain", "Gain+LOH", "WGD", "WGD+Gain",
                    "Hom. Del.", "Other", "Chromothripsis")

classify_cna_event <- function(minorAlleleCopyNumber,majorAlleleCopyNumber, Gender, chromosome, copyNumber, mean_copyNumber)
{
  dplyr::case_when(
    (majorAlleleCopyNumber >= 1.85 &
       majorAlleleCopyNumber <= 2.15) &
      (minorAlleleCopyNumber >= 1.85 &
         minorAlleleCopyNumber <= 2.15) &
      (mean_copyNumber >= 3) ~ "WGD",
    (majorAlleleCopyNumber >= 2.85) &
      (minorAlleleCopyNumber >= 1.85 &
         minorAlleleCopyNumber <= 2.15) &
      (mean_copyNumber >= 3) ~ "WGD+Gain",
    (majorAlleleCopyNumber >= 0.5 &
       majorAlleleCopyNumber <= 1.15) &
      (minorAlleleCopyNumber >= 0.90 &
         minorAlleleCopyNumber <= 1.15) ~ "Diploid/Haploid-X",
    (majorAlleleCopyNumber >= 0.5 &
       majorAlleleCopyNumber <= 1.15) &
      (minorAlleleCopyNumber <= 0.25) &
      (chromosome != "chrX" |
         (chromosome == "chrX" & Gender == "female")) ~ "Loss",
    (majorAlleleCopyNumber >= 0.97 &
       majorAlleleCopyNumber <= 1.85) &
      (minorAlleleCopyNumber <= 0.25) &
      (chromosome != "chrX" |
         (chromosome == "chrX" &
            Gender == "female")) ~ "Loss + Subclonal CNLOH",
    (majorAlleleCopyNumber >= 0.5 &
       majorAlleleCopyNumber <= 1.15) &
      (minorAlleleCopyNumber <= 0.25) &
      (chromosome == "chrX" &
         Gender == "male") ~ "Diploid/Haploid-X",
    (majorAlleleCopyNumber >= 0.5 &
       majorAlleleCopyNumber <= 1.15) &
      (minorAlleleCopyNumber <= 0.7) &
      (chromosome != "chrX" |
         (chromosome == "chrX" &
            Gender == "female")) ~ "Subclonal Loss",
    (majorAlleleCopyNumber >= 0.5 &
       majorAlleleCopyNumber <= 1.15) &
      (minorAlleleCopyNumber <= 0.93) &
      (!(chromosome %in% c("chrX", "chr19")) |
         (chromosome == "chrX" &
            Gender == "female")) ~ "Minor Subclonal Loss",
    (majorAlleleCopyNumber >= 1.85 &
       majorAlleleCopyNumber <= 2.15) &
      (minorAlleleCopyNumber <= 0.25) ~ "CNLOH",
    (majorAlleleCopyNumber >= 1.85) &
      (minorAlleleCopyNumber >= 0.85) ~ "Gain",
    #& minorAlleleCopyNumber <= 1.15
    (majorAlleleCopyNumber >= 1.15) &
      (minorAlleleCopyNumber >= 0.85) ~ "Subclonal Gain",
    #& minorAlleleCopyNumber <= 1.15
    (majorAlleleCopyNumber >= 1.85) &
      (minorAlleleCopyNumber <= 0.85) ~ "Gain+LOH",
    copyNumber <= 0.5 ~ "Hom. Del.",
    TRUE ~ "Other"
  )
}
