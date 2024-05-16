library(tidyverse)
library(lubridate)
library(ggrepel)
library(patchwork)

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
    select(`Patient ID`, dx_date) %>% 
    distinct() %>% 
    pull(dx_date) 
  
  swimmer_events[["presentation_timepoint"]] <- 
    pad_and_format_date(swimmer_events[["presentation_timepoint"]],
                        date_accuracy(swimmer_events[["presentation_timepoint"]]))
  
  followup_duration <- patient_anno %>% pull(Post_diagnosis_follow_up_months) %>% unique()
  vital_duration <- patient_anno %>% pull(Time_to_death_months) %>% unique() %>% as.numeric()
  
  swimmer_events[["followup_end_date"]] <- swimmer_events[["presentation_timepoint"]] + months(as.numeric(followup_duration))
  
  if (!is.na(vital_duration) & is.numeric(vital_duration)) {
    swimmer_events[["death_timepoint"]] <- swimmer_events[["presentation_timepoint"]] + months(as.numeric(vital_duration))
  }
  
  swimmer_events[["sample_timepoint"]] <- patient_anno %>% 
    select(PublicationID, sx_date, is_primary_or_met) %>% rowwise() %>% 
    mutate(sx_date = pad_and_format_date(sx_date,
                                         date_accuracy(sx_date)))
  
  additional_surgery_annotation <- patient_anno %>% 
    select(PublicationID, surgical_events_excluding_profiled) %>% 
    filter(!is.na(surgical_events_excluding_profiled)) 
  
  if(nrow(additional_surgery_annotation) > 0) {
    swimmer_events[["additional_surgical_events"]] <- parse_additional_surgeries(additional_surgery_annotation)
  }
  
  sample_treament_annotation <- patient_anno %>% 
    select(PublicationID, non_surgical_therapy) %>% 
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
  mutate(multisample=n()>1) %>% 
  filter(multisample | !is.na(non_surgical_therapy)) %>% 
  filter(!(non_surgical_therapy %in% c("Interferon","SSA"))) %>% 
  filter(A5_ID != "E168-1") #Unreliable information

a5_anno_use <- a5_anno_use %>% 
  filter(Exclude != "Y") %>% 
  dplyr::rename("dx_date"=`Date of first pathologic diagnosis of PPGL (MM/YYYY)`,
                "sx_date"= `Date of resection (DD/MM/YYYY)`) %>% 
  group_by(`Patient ID`) %>% 
  group_split() %>% 
  purrr::set_names(purrr::map_chr(., ~.x[["Patient ID"]][1]))


swimmer_events <- map(.x = a5_anno_use, .f = anno_to_swimmer)

swimmer_events <- map(.x = swimmer_events, .f=translate_dates_to_days)


##################
# Generate Plots #
##################

plot_tables <- map(.x = swimmer_events, .f = event_to_plottable)


swimmer_plots <- list()
for (patient in names(plot_tables))
{
  
  plot_data <- plot_tables[[patient]]
  
  plot_data <- plot_data %>% 
    mutate(event_class = factor(event_class, 
                                levels = intersect(c("Treatment","Surgery","Observation"), plot_data$event_class)))
  
  #remove end_observation timepoint if the same as death timepoint
  if("Death" %in% plot_data$event_name)
  {
    death_time = plot_data %>% filter(event_name=="Death") %>% pull(event_end)
    obs_end_time = plot_data %>% filter(event_name=="End Follow-up") %>% pull(event_end)
    if(obs_end_time==death_time) {
      plot_data <- plot_data %>% filter(event_name != "End Follow-up")
    }
  }
  
  #Add offsets for overlapping events
  # plot_data <- plot_data %>% 
  #   arrange(event_class, event_start, event_end) %>% 
  #   mutate(offset=ifelse(lag(event_class)==event_class & event_start-lag(event_end, default = 0) <= (max(event_end)*0.1), 0.3,0)) %>% 
  #   mutate(offset=ifelse(lead(offset, default = 0) != 0, -0.3, offset))
  # plot_data$offset <- plot_data$offset * rep(c(-1,1), ceiling(length(plot_data$offset)/2))[1:length(plot_data$offset)]
  # 
  plot_data <- plot_data %>% 
    arrange(event_class, event_start, event_end) %>% 
    mutate(overlap=ifelse(lag(event_class)==event_class & event_start-lag(event_end, default = 0) <= (max(event_end)*0.1) |
                           lead(event_class)==event_class & event_end-lead(event_start, default = 0) <= (max(event_end)*0.1), "overlap","no_overlap")) %>%
    mutate(overlap=replace_na(overlap, "no_overlap")) %>% 
   group_by(event_class) %>% 
    mutate(offset=rep(c(-0.3,0.3), ceiling(n()/2))[1:n()],
           offset=ifelse(overlap=="overlap", offset,0))

  swim <- ggplot()  
  
  
  #Observations
  swim <- swim +  geom_point(data=plot_data %>% filter(event_class=="Observation"),
                             mapping=aes(x=event_start, 
                                         y=event_class, 
                                         shape=event_name), size=2) 
  
  swim <- swim +  geom_segment(data=plot_data %>% filter(event_class=="Observation") %>% 
                                 mutate(event_start=min(event_start), event_end=max(event_end)) %>% 
                                 slice_head(n=1),
                               mapping=aes(x=event_start, 
                                           xend=event_end, 
                                           y=event_class, 
                                           yend=event_class)) 
  
  #Surgery
  swim <- swim + geom_point(data=plot_data %>% filter(event_class=="Surgery"),
                            mapping=aes(x=event_start, 
                                        y=as.numeric(event_class)+offset)) 
  
  swim <- swim + ggrepel::geom_text_repel(data=plot_data %>% filter(event_class=="Surgery") %>% 
                             mutate(event_name=gsub("E...-","",event_name),
                                    event_name=gsub(" - ","\n",event_name)),
                           mapping=aes(x=event_start, 
                                       y=as.numeric(event_class)+offset, 
                                       label=event_name),
                           nudge_y=-0.2) 
    
  if("Treatment" %in% plot_data$event_class)
  {
  #Treatment(Date imprecise)
  swim <- swim + geom_point(data=plot_data %>% 
                              filter(event_class=="Treatment", 
                                     event_annotation=="Year"),
                            mapping=aes(x=(event_end-((event_end-event_start)/2)), 
                                        y=as.numeric(event_class) + offset)) 
  
  swim <- swim + geom_segment(data=plot_data %>% 
                                filter(event_class=="Treatment", 
                                       event_annotation=="Year"),
                              mapping=aes(x=event_start, 
                                          xend=event_end, 
                                          y=as.numeric(event_class) + offset,
                                          yend=as.numeric(event_class) + offset), 
                              linetype="21") 

  #Treatment(Date precise)
  swim <- swim + geom_point(data=plot_data %>% 
                              filter(event_class=="Treatment", 
                                     event_annotation != "Year") %>% 
                              pivot_longer(cols = c(event_start, event_end), 
                                           names_to = "timepoint", 
                                           values_to = "event_start"),
                            mapping=aes(x=event_start , y=as.numeric(event_class)+offset)) 
  
  swim <- swim + 
    geom_segment(data=plot_data %>% filter(event_class=="Treatment", 
                                           event_annotation != "Year"),
                       mapping=aes(x=event_start, 
                                   xend=event_end, 
                                   y=as.numeric(event_class)+offset,
                                   yend=as.numeric(event_class)+offset))
  
  
  #Treatment label 
  swim <- swim + 
    ggrepel::geom_text_repel(data=plot_data %>% filter(event_class=="Treatment") %>% 
                     mutate(event_name=gsub("[(].+[)]","", event_name)),
                   mapping=aes(x=(event_end-((event_end-event_start)/2)), 
                               y=as.numeric(event_class)+offset,
                               label=event_name),
                   nudge_y=-0.2)   
  }
  
  swim <- swim + 
    geom_rect(data=plot_data %>% 
                group_by(event_class) %>% 
                     mutate(max_offset=max(offset)) %>%  
                     ungroup() %>% 
                     transmute(event_class=event_class,
                               max_offset=max_offset,
                               event_start=min(event_start)-(0.05 * max(event_end)), 
                               event_end=max(event_end)+(0.05 * max(event_end))) %>% 
                     distinct() %>% 
                     mutate(ymin=seq(0.5, nlevels(event_class) - 0.5),
                            ymin=dplyr::if_else("Treatment" %in% plot_data$event_class & event_class=="Treatment", ymin-max_offset, ymin),
                            ymax=seq(1.5, nlevels(event_class) + 0.5)),
                   aes(xmin=event_start, xmax=event_end, ymin=ymin, ymax=ymax), 
              alpha=c(0.3,0.1,0.3)[1:nlevels(plot_data$event_class)])
  
  swimmer_plots[[patient]] <- 
    swim + theme_bw() + 
    theme(rect=element_blank(), 
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    theme(panel.grid = element_blank(), 
          axis.line.x.bottom = element_line()) + 
    scale_y_discrete(drop=F) +
    xlab("Days") + ggtitle(patient)
}

# all_plots <- Reduce(x=swimmer_plots, f="/")
# ggsave(plot = all_plots, filename = "./a5/sample_annotation/figures/swimmer_plots.pdf", 
#        width = 12, 
#        height = 
#          60,
#        units = "in",
#        limitsize = F)
