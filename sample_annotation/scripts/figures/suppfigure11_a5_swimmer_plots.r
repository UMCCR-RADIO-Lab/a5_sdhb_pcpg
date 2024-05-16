library(tidyverse)
library(swimplot)

setwd("/g/data/pq08/projects/ppgl")

################
# Data loaders #
################

source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)


########################
# Function Definitions #
########################

#########
# Date Parsing
#########

date_accuracy <- function(string_date)
{
  case_when(
    is.na(string_date) ~ "Uncertain",
    grepl("^[0-9]{4}$", string_date) ~ "Year",
    grepl("^[0-9]{1,2}/[0-9]{4}$", string_date) ~ "Month",
    grepl("^[0-9]{1,2}/[0-9]{1,2}/[0-9]{4}$", string_date) ~ "Day" 
  ) %>% return()
}

pad_and_format_date <- function(string_date, date_accuracy)
{
  if (!all(grepl("([0-9]{1,2}/)?([0-9]{1,2}/)?[0-9]{4}", string_date))){
    stop("Unrecognised date format:", toString(string_date), ". Needed DD/MM/YYYY or MM/YYYY or YYYY.")
  }
  
  case_when(date_accuracy == "Year" ~ dmy(paste0("01/01/", string_date), quiet = T),
            date_accuracy == "Month" ~ dmy(paste0("01/", string_date),quiet = T),
            date_accuracy == "Day" ~ dmy(string_date, quiet = T)
  ) %>% return()
}

#Reformat dates and fill in missing dates for treatments
conform_treatment_dates <- function(start_date, end_date)
{ 
  
  start_date_accuracy = date_accuracy(start_date)
  
  end_date_accuracy = date_accuracy(end_date)
  
  start_date_formatted = pad_and_format_date(start_date, start_date_accuracy)
  
  if(is.na(start_date_formatted)) { stop("Failed to parse start date:", start_date) }
  
  end_date_formatted = case_when(
    end_date_accuracy == "Uncertain" & start_date_accuracy == "Year" ~ start_date_formatted + years(1),
    end_date_accuracy == "Uncertain" & start_date_accuracy == "Month" ~ start_date_formatted + months(1),
    end_date_accuracy == "Uncertain" & start_date_accuracy == "Day" ~ start_date_formatted + months(1),
    end_date_accuracy == "Year"  ~ dmy(paste0("01/01/", end_date), quiet = T),
    end_date_accuracy == "Month" ~ dmy(paste0("01/", end_date), quiet = T),
    end_date_accuracy == "Day" ~ dmy(end_date, quiet = T),
  ) 
  
  if(is.na(end_date_formatted)) { stop("Failed to parse or guess end date:", end_date) }
  
  return(
    list(start_date = start_date_formatted, 
         start_date_accuracy = start_date_accuracy, 
         end_date_accuracy = end_date_accuracy,
         end_date = end_date_formatted)
  )
}

translate_dates_to_days <- function(swimmer_events) 
{
  
  zero_reference <- swimmer_events[["presentation_timepoint"]] 
  
  swimmer_events[["presentation_timepoint"]] <- 0
  
  swimmer_events[["followup_end_date"]] <- as.numeric(difftime(time1 = swimmer_events[["followup_end_date"]], 
                                                               time2 = zero_reference, 
                                                               units = "days"))
  
  if("death_timepoint" %in% names(swimmer_events)) {
    swimmer_events[["death_timepoint"]] <- 
      as.numeric(difftime(time1 = swimmer_events[["death_timepoint"]], 
                          time2 = zero_reference, 
                          units = "days"))
  }
  
  if("metastasis_free_duration" %in% names(swimmer_events)) {
    if (!is.na(swimmer_events[["metastasis_free_duration"]] ))
    {
      swimmer_events[["metastasis_free_duration"]] <- swimmer_events[["metastasis_free_duration"]] * 30.4 #Average month
      #Handle margin of error in conversion of follow up time and metastasis free survival
      if(swimmer_events[["metastasis_free_duration"]] > 0  & swimmer_events[["followup_end_date"]] > 0){
        ratio <- swimmer_events[["metastasis_free_duration"]] /  swimmer_events[["followup_end_date"]]
        if (ratio > 0.98 & ratio < 1.03) { 
          swimmer_events[["metastasis_free_duration"]] <-  swimmer_events[["followup_end_date"]] }
      }
    }
  }
  
  swimmer_events[["sample_timepoint"]][["sx_date"]] <- 
    as.numeric(difftime(time1 = swimmer_events[["sample_timepoint"]][["sx_date"]], 
                        time2 = zero_reference, 
                        units = "days"))
  
  if ("additional_surgical_events" %in% names(swimmer_events)) {
    swimmer_events[["additional_surgical_events"]][["start_date"]] <- 
      as.numeric(difftime(time1 = swimmer_events[["additional_surgical_events"]][["start_date"]], 
                          time2 = zero_reference, 
                          units = "days"))
    swimmer_events[["additional_surgical_events"]][["end_date"]] <- 
      as.numeric(difftime(time1 = swimmer_events[["additional_surgical_events"]][["end_date"]], 
                          time2 = zero_reference, 
                          units = "days"))
  }
  
  if ("treament_timepoint" %in% names(swimmer_events)) {
    swimmer_events[["treament_timepoint"]][["start_date"]] <- 
      as.numeric(difftime(time1 = swimmer_events[["treament_timepoint"]][["start_date"]], 
                          time2 = zero_reference, 
                          units = "days"))
    swimmer_events[["treament_timepoint"]][["end_date"]] <- 
      as.numeric(difftime(time1 = swimmer_events[["treament_timepoint"]][["end_date"]], 
                          time2 = zero_reference, 
                          units = "days"))
  }
  
  return(swimmer_events)
}

#########
# Annotation Parsing
#########

#Parse treatment events
parse_treatments <- function(sample_treament_annotation)
{
  treament_timepoint <- sample_treament_annotation %>% 
    separate_rows(non_surgical_therapy, sep = ";") %>% 
    separate(col = non_surgical_therapy, 
             into=c("date", "treatment_type"), sep=":") %>% 
    mutate(date = stringr::str_trim(date))
  
  treament_timepoint <- treament_timepoint %>% 
    separate(col = date, 
             into=c("start_date", "end_date"), sep="-")
  
  treatments <- tibble(start_date=ymd(), 
                       start_date_accuracy=vector(mode="character"), 
                       end_date=ymd(),
                       end_date_accuracy=vector(mode="character"))
  
  
  treatments_formatted <- map2(.x = treament_timepoint[["start_date"]], 
                               .y= treament_timepoint[["end_date"]], 
                               .f = conform_treatment_dates) %>% 
    map(.f= \(x) add_case(.data = treatments, !!!x)) %>% 
    bind_rows() %>% 
    mutate(event="Tx") 
  
  treatments_formatted[["treatment_type"]] <- treament_timepoint[["treatment_type"]]
  
  #remove duplicated cause by repeating data in the clinical table
  treatments_formatted <- treatments_formatted %>% distinct()
  
  return(treatments_formatted)
  
}

#Parse treatment events
parse_additional_surgeries <- function(additional_surgery_annotation)
{
  surgical_timepoint <- additional_surgery_annotation %>% 
    separate_rows(surgical_events_excluding_profiled, sep = ";") %>% 
    separate(col = surgical_events_excluding_profiled, 
             into=c("date", "surgery"), sep=":") %>% 
    mutate(date = stringr::str_trim(date))
  
  surgeries_formatted <- surgical_timepoint %>% 
    mutate(
      start_date_accuracy = date_accuracy(date),
      start_date = pad_and_format_date(date, start_date_accuracy), 
      end_date = start_date,
      end_date_accuracy = start_date_accuracy) %>% 
    dplyr::select(-date, -PublicationID)
  
  #remove duplicated cause by repeating data in the clinical table
  surgeries_formatted <- surgeries_formatted %>% distinct()
  
  return(surgeries_formatted)
  
}

anno_to_swimmer <- function(patient_anno)
{
  
  swimmer_events <- list()
  
  swimmer_events[["presentation_timepoint"]] <- 
    patient_anno %>% 
   dplyr::select(`Patient ID`, dx_date) %>% 
    distinct() %>% 
    pull(dx_date) 
  
  swimmer_events[["presentation_timepoint"]] <- 
    pad_and_format_date(swimmer_events[["presentation_timepoint"]],
                        date_accuracy(swimmer_events[["presentation_timepoint"]]))
  
  followup_duration <- patient_anno %>% pull(Post_diagnosis_follow_up_months) %>% unique()
  if(is.na(followup_duration) | followup_duration == "No data")
  {
    followup_duration <- 0
  } else {
    followup_duration <- as.numeric(followup_duration)
  }
  
  vital_duration <- patient_anno %>% pull(Time_to_death_months) %>% unique() 
  if(is.na(vital_duration) | vital_duration == "Not applicable")
  {
    vital_duration <- followup_duration 
  } else 
  {
    vital_duration <- as.numeric(vital_duration)
  }
  
  
  swimmer_events[["followup_end_date"]] <- swimmer_events[["presentation_timepoint"]] + months(floor(as.numeric(followup_duration)))
  
  swimmer_events[["metastasis_free_duration"]] <- patient_anno %>% pull(Metastasis_free_survival_months) %>% unique()
  if(is.na(swimmer_events[["metastasis_free_duration"]]) | swimmer_events[["metastasis_free_duration"]]  == "No data")
  {
    swimmer_events[["metastasis_free_duration"]]  <- followup_duration 
  } else 
  {
    swimmer_events[["metastasis_free_duration"]]  <- as.numeric(swimmer_events[["metastasis_free_duration"]])
  }
  
  if (!is.na(vital_duration) & is.numeric(vital_duration) & unique(patient_anno$`Patient is alive`) == "No") {
    swimmer_events[["death_timepoint"]] <- swimmer_events[["presentation_timepoint"]] + months(as.numeric(vital_duration))
  }
  
  swimmer_events[["sample_timepoint"]] <- patient_anno %>% 
   dplyr::select(PublicationID, sx_date, is_primary_or_met) %>% rowwise() %>% 
    mutate(sx_date = pad_and_format_date(sx_date,
                                         date_accuracy(sx_date)))
  
  additional_surgery_annotation <- patient_anno %>% 
   dplyr::select(PublicationID, surgical_events_excluding_profiled) %>% 
    filter(!is.na(surgical_events_excluding_profiled)) 
  
  if(nrow(additional_surgery_annotation) > 0) {
    swimmer_events[["additional_surgical_events"]] <- parse_additional_surgeries(additional_surgery_annotation)
  }
  
  sample_treament_annotation <- patient_anno %>% 
   dplyr::select(PublicationID, non_surgical_therapy) %>% 
    filter(!is.na(non_surgical_therapy)) 
  
  if(nrow(sample_treament_annotation) > 0) {
    swimmer_events[["treament_timepoint"]] <- parse_treatments(sample_treament_annotation)
    
    #clamp imprecise end dates at follow-up end date
    swimmer_events[["treament_timepoint"]] <- swimmer_events[["treament_timepoint"]] %>% 
      mutate(end_date=dplyr::if_else(difftime(time1 = end_date, 
                                              time2 = swimmer_events[["followup_end_date"]]) < 0, 
                                     end_date, 
                                     swimmer_events[["followup_end_date"]]))
    
  }
  
  return(swimmer_events)
  
  
}

########
# Plot preparation 
########


event_to_plottable <- function(patient_swimmer_events)
{
  
  plot_table <- tibble(event_start=vector(mode = "numeric"), 
                       event_end=vector(mode = "numeric"), 
                       event_class=vector(mode = "character"),
                       event_name=vector(mode = "character"), 
                       event_annotation=vector(mode = "character"))
  
  plot_table <- add_row(plot_table, event_start=patient_swimmer_events[["presentation_timepoint"]],
                        event_end=patient_swimmer_events[["presentation_timepoint"]],
                        event_class="Observation",
                        event_name="Presentation",
                        event_annotation="Presentation")
  
  if ("death_timepoint" %in% names(patient_swimmer_events))
  {
    plot_table <- add_row(plot_table, 
                          event_start=patient_swimmer_events[["death_timepoint"]],
                          event_end=patient_swimmer_events[["death_timepoint"]],
                          event_class="Observation",
                          event_name="Death",
                          event_annotation="Death")
  }
  
  plot_table <- add_row(plot_table, 
                        event_start=patient_swimmer_events[["followup_end_date"]],
                        event_end=patient_swimmer_events[["followup_end_date"]],
                        event_class="Observation",
                        event_name="End Follow-up",
                        event_annotation="End Follow-up")
  
  plot_table <- add_row(plot_table, 
                        event_start=patient_swimmer_events[["sample_timepoint"]][["sx_date"]],
                        event_end=patient_swimmer_events[["sample_timepoint"]][["sx_date"]],
                        event_class="Surgery",
                        event_name=patient_swimmer_events[["sample_timepoint"]][["PublicationID"]],
                        event_annotation=patient_swimmer_events[["sample_timepoint"]][["is_primary_or_met"]])
  
  if ("treament_timepoint" %in% names(patient_swimmer_events))
  {
    plot_table <- add_row(plot_table, 
                          event_start=patient_swimmer_events[["treament_timepoint"]][["start_date"]],
                          event_end=patient_swimmer_events[["treament_timepoint"]][["end_date"]],
                          event_class="Treatment",
                          event_name=patient_swimmer_events[["treament_timepoint"]][["treatment_type"]],
                          event_annotation=patient_swimmer_events[["treament_timepoint"]][["start_date_accuracy"]])
  }
  
  if ("additional_surgical_events" %in% names(patient_swimmer_events))
  {
    plot_table <- add_row(plot_table, 
                          event_start=patient_swimmer_events[["additional_surgical_events"]][["start_date"]],
                          event_end=patient_swimmer_events[["additional_surgical_events"]][["end_date"]],
                          event_class="Surgery",
                          event_name=patient_swimmer_events[["additional_surgical_events"]][["surgery"]],
                          event_annotation="Surgery not profiled")
  }
  
  return(plot_table)
  
}

######################
# Process Annotation #
######################

a5_anno_use <- a5_anno %>% 
  group_by(`Patient ID`) %>% 
  filter(Exclude == "N") %>% 
  #mutate(multisample=n()>1) %>% 
  #filter(multisample | !is.na(non_surgical_therapy)) %>% 
  #filter(any(tumour_metastasised =="Yes")) %>% 
  #filter(!(non_surgical_therapy %in% c("Interferon","SSA"))) %>% 
  #filter(A5_ID != "E168-1") %>% #Unreliable information 
  mutate(tumour_metastasised = ifelse(A5_ID == "E159-2", "Yes", tumour_metastasised)) %>% 
  mutate(Metastasis_free_survival_months=case_when(
    tumour_metastasised=="No" ~ Overall_Survival, 
    tumour_metastasised == "Short follow up" ~ Post_diagnosis_follow_up_months, 
    Metastasis_free_survival_months == "No Data" ~ NA,
    TRUE ~ Metastasis_free_survival_months)) %>% 
  mutate(tumour_metastasised = ifelse(A5_ID == "E159-2", "No", tumour_metastasised)) %>% 
  mutate(`Date of resection (DD/MM/YYYY)` = ifelse(`Date of resection (DD/MM/YYYY)` == "NA" | is.na(`Date of resection (DD/MM/YYYY)`), 
                                                   `Date of first pathologic diagnosis of PPGL (MM/YYYY)`, 
                                                   `Date of resection (DD/MM/YYYY)`)) %>% 
  mutate(Post_diagnosis_follow_up_months = gsub("<","",Post_diagnosis_follow_up_months)) %>% 
  mutate(Overall_Survival = gsub("<","",Overall_Survival)) %>% 
  mutate(`Date of first pathologic diagnosis of PPGL (MM/YYYY)` = ifelse(A5_ID=="E204-1", `Date of resection (DD/MM/YYYY)`, `Date of first pathologic diagnosis of PPGL (MM/YYYY)`))



a5_anno_use <- a5_anno_use %>% 
  filter(Exclude != "Y") %>% 
  dplyr::rename("dx_date"=`Date of first pathologic diagnosis of PPGL (MM/YYYY)`,
                "sx_date"= `Date of resection (DD/MM/YYYY)`) %>% 
  group_by(`Patient ID`) %>% 
  group_split() %>% 
  purrr::set_names(purrr::map_chr(., ~.x[["Patient ID"]][1]))


swimmer_events <- map(.x = a5_anno_use, .f = anno_to_swimmer)

swimmer_events <- map(.x = swimmer_events, .f=translate_dates_to_days)


swimmer_base <- purrr::map2(
  .x = swimmer_events, 
  .y = names(swimmer_events), 
  .f = function(patient_anno, patient_id) {
    
    last_sx <- patient_anno$sample_timepoint %>% ungroup() %>% 
      slice_max(order_by = sx_date) %>% pull(sx_date) %>%  unique()
    
    return_table = NULL
    if(patient_id == "E165") { 
      #making an assumption that follow up duration was relative to surgery for this case
      return_table <- data.frame(patient_id=patient_id, 
                                 clinical_behaviour=c("Limited data","Metastatic Disease"), 
                                 end=c(3650, (last_sx + patient_anno$followup_end_date)))
    } else if(patient_anno$followup_end_date == 0) { 
      last_sx <- patient_anno$sample_timepoint %>%
        slice_max(order_by = sx_date) %>% pull(sx_date)
      return_table <- data.frame(patient_id=patient_id, clinical_behaviour="Short follow up", end=last_sx)
    } else if(patient_anno$metastasis_free_duration == 0) { 
      return_table <- data.frame(patient_id=patient_id, clinical_behaviour="Metastatic Disease", end=patient_anno$followup_end_date)
    } else {
      return_table <- data.frame(patient_id=patient_id, 
                                 clinical_behaviour=c("No Metastatic Disease", "Metastatic Disease"), 
                                 end=c(patient_anno$metastasis_free_duration, patient_anno$followup_end_date))
    }
    
    if(!is.null(patient_anno$death_timepoint) && patient_anno$death_timepoint > patient_anno$followup_end_date) { 
      return_table <- bind_rows(return_table,
                                data.frame(patient_id=patient_id, 
                                           clinical_behaviour="Limited data", 
                                           end=patient_anno$death_timepoint))
    }
    
    return(return_table)
    
  }) %>%  bind_rows()

swimmer_sx <- purrr::map2(
  .x = swimmer_events, 
  .y = names(swimmer_events), 
  .f = function(patient_anno, patient_id) {
    
    return_table <- tibble(patient_id = vector(mode = "character"), 
                           event = vector(mode = "character"),
                           time = vector(mode = "numeric"))
    
    if ("death_timepoint" %in% names(patient_anno))
    {
      return_table <- return_table %>% add_row(tibble(patient_id = patient_id, event = "Death", time = patient_anno$death_timepoint))
    }
    
    if ("sample_timepoint" %in% names(patient_anno)) {
      if (nrow(patient_anno[["sample_timepoint"]]) > 0)
      {
        return_table <- return_table %>% 
          add_row(tibble(patient_id = patient_id, 
                         event = paste("Surgery -", patient_anno[["sample_timepoint"]][["is_primary_or_met"]]), 
                         time = patient_anno[["sample_timepoint"]][["sx_date"]]))
      }}
    
    if ("additional_surgical_events" %in% names(patient_anno)) {
      if (nrow(patient_anno[["additional_surgical_events"]]) > 0)
      {
        return_table <- return_table %>% 
          add_row(tibble(patient_id = patient_id, 
                         event = paste("Surgery - Not profiled"), 
                         time = patient_anno[["additional_surgical_events"]][["start_date"]]))
      }}
    
    return(return_table)
    
  }) %>%  bind_rows() %>% as.data.frame()


swimmer_tx <- purrr::map2(
  .x = swimmer_events, 
  .y = names(swimmer_events), 
  .f = function(patient_anno, patient_id) {
    
    return_table <- tibble(patient_id = vector(mode = "character"), 
                           treatment = vector(mode = "character"),
                           quantity = vector(mode = "character"),
                           time_start = vector(mode = "numeric"),
                           time_end = vector(mode = "numeric"),
                           start_date_accuracy = vector(mode = "character"),
                           end_date_accuracy = vector(mode = "character"))
    
    if ("treament_timepoint" %in% names(patient_anno)) {
      if (nrow(patient_anno[["treament_timepoint"]]) > 0){
        TOI <- patient_anno[["treament_timepoint"]] %>% filter(grepl("MIBG|CVD|TMZ|Luta",treatment_type)) %>% 
          separate(col = treatment_type, into = c("treatment","quantity"), sep=" - ")
        if (nrow(TOI) > 0) {
          return_table <- return_table %>% 
            add_row(tibble(patient_id = patient_id, 
                           treatment = TOI$treatment, 
                           quantity = TOI[["quantity"]],
                           time_start = TOI[["start_date"]],
                           time_end = TOI[["end_date"]],
                           start_date_accuracy = TOI[["start_date_accuracy"]],
                           end_date_accuracy = TOI[["end_date_accuracy"]]))
        }
      }
    }
    
    
  }) %>%  bind_rows() %>% as.data.frame()

##################
# Generate Plots #
##################

swimmer_plot <- swimmer_plot(df=swimmer_base,id='patient_id',end='end',name_fill='clinical_behaviour', width=.8)


swimmer_plot <- swimmer_plot +
  swimmer_lines(df_lines=swimmer_tx,id='patient_id',start =
                  'time_start',end='time_end',name_col='treatment', size=1.5 ) #name_size = "start_date_accuracy"

swimmer_plot <- swimmer_plot + swimmer_points(df_points= swimmer_sx,
                                              id='patient_id',
                                              time='time',
                                              name_shape ='event',
                                              size=1.5,
                                              fill='white',
                                              col='black')

swimmer_plot <- swimmer_plot + 
  scale_fill_manual(name="Clinical behaviour", values=c(`No Metastatic Disease`="#98ca8cff", `Metastatic Disease`="#de5353ff",  "Limited data"="#e2e0d9ff", "Short follow up"="#f4e764ff")) +
  scale_color_manual(name="Treatment", values=c(CVD="#f4e764ff", `I131-MIBG`="#cadae8ff", `CAPE/TMZ`="#d6a097ff", Lutate="#efecdcff")) +
  scale_shape_manual(name="Event", values=c("Death"=4,"Surgery - Metastasis"=8,"Surgery - Primary"=16, "Surgery - Recurrent"=3, "Surgery - Not profiled"=1)) +
  scale_size_manual(values=c("Day"=3, "Month"=1.5, "Year"=0.6))

swimmer_plot
