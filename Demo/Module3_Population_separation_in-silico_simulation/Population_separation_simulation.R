library(tidyverse)
library(fs)
library(furrr)
library(ggplot2)
library(patchwork)
library(yaml)
library(optparse)
library(dplyr)

### Set to TRUE if you want to run import and processing even if file already exists
force_redo=TRUE
tracks_provided=NULL

### Checks if being run in GUI (e.g. Rstudio) or command line
if (interactive()) {
  BGT_dir = paste0(dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path))),"/")
  pars = yaml.load_file("config_template.yml")
} else {
  BGT_dir <- "~/BGT"
  option_list = list(
    make_option(c("-c", "--config"), type="character", default=NULL, help="Path to the BGT config file", metavar="character"),
    make_option(c("-f", "--force_redo"), action="store_true", default=FALSE, help="Force the pipeline to re-import data even if files exists"),
    make_option(c("-t", "--tracks_rds"), type="character", default=NULL, help="(Optional) Path to RDS file containing processed T cell track data", metavar="character")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  if (is.null(opt$config)){
    print_help(opt_parser)
    stop("Config file -c|--config, must be supplied", call.=FALSE)
  }
  pars = yaml.load_file(opt$config)
  force_redo=opt$force_redo
  tracks_provided=opt$tracks_rds
}
# Set path to you config manually
# pars = yaml.load_file("config_template.yml")
output_dir= paste0(pars$BGT_dir, pars$output_dir)

dir.create(output_dir, recursive=TRUE)

n_timepoints = pars$n_timepoints
first_timepoint_interval_starts = pars$first_timepoint_interval_starts
first_timepoint_interval_ends = pars$first_timepoint_interval_ends
reoccuring_time = pars$reoccuring_time

master_clust_Live <- readRDS(pars$classified_tcell_track_data_filepath_rds)
time_intervals <- list()

for (i in 1:n_timepoints) {
  interval_start <- first_timepoint_interval_starts + (i - 1) * reoccuring_time
  interval_end <- first_timepoint_interval_ends + (i - 1) * reoccuring_time
  time_intervals[[i]] <- c(interval_start, interval_end)
}

t1 = c(first_timepoint_interval_starts:first_timepoint_interval_ends)
t2 = c(interval_start:interval_end)

color_palette <- c("darkolivegreen3",
                   "seagreen3",
                   "forestgreen",
                   "dodgerblue",
                   "cyan1",
                   "indianred",
                   "firebrick",
                   "brown1")

Behavioral_signatures <- c(
  "Static",
  "Lazy",
  "Slow-scanner",
  "Medium-scanner",
  "Super-scanner",
  "Tickler",
  "Engager",
  "Super-engager")  



# Generate timepoint ranges
time_intervals <- list()
for (i in 1:n_timepoints) {
  interval_start <- first_timepoint_interval_starts + (i - 1) * reoccuring_time
  interval_end <- first_timepoint_interval_ends + (i - 1) * reoccuring_time
  time_intervals[[i]] <- c(interval_start, interval_end)
}

# Set up parallel processing
plan(multisession, workers = 4)
# Process data and plot for each time interval
plots <- future_map(time_intervals, function(interval) {
  subset_data <- master_clust_Live %>%
    filter(Time >= interval[1], Time <= interval[2])
  
  cd8_data <- subset_data %>% filter(str_detect(tcell_line, "CD8"), cluster != "1")
  cd4_data <- subset_data %>% filter(str_detect(tcell_line, "CD4"), cluster != "1")
  
  plot_cd8 <- ggplot(cd8_data, aes(x = factor(contact), fill = factor(cluster))) +
    geom_bar(position = "fill") +
    scale_fill_manual(name = "Behavioral\nsignature", labels = Behavioral_signatures, values=color_palette) +
    scale_x_discrete(labels=c("No Contact", "Contact"))+
    labs(y = "Percentage", x = "Contact", title = sprintf("CD8 Engagement between timepoint %s", paste(interval, collapse = "-"))) +
    theme_minimal() +
    coord_flip() +
    facet_wrap(~tcell_line)
  print(plot_cd8)
  
  plot_cd4 <- ggplot(cd4_data, aes(x = factor(contact), fill = factor(cluster))) +
    geom_bar(position = "fill") +
    scale_fill_manual(name = "Behavioral\nsignature", labels = Behavioral_signatures, values=color_palette) +
    scale_x_discrete(labels=c("No Contact", "Contact"))+
    labs(y = "Percentage", x = "Contact", title = sprintf("CD4 Engagement between timepoint %s", paste(interval, collapse = "-"))) +
    theme_minimal() +
    coord_flip() +
    facet_wrap(~tcell_line)
  print(plot_cd4)
  
  list(CD8 = plot_cd8, CD4 = plot_cd4)
})

if (n_timepoints == 1) {
  # Save plots
  walk2(plots, seq_along(plots), ~{
    ggsave(filename = file.path(output_dir, sprintf("plot_cd8_timepoint_%s.png", .y)), plot = .x$CD8, width = 8, height = 6)
    ggsave(filename = file.path(output_dir, sprintf("plot_cd4_timepoint_%s.png", .y)), plot = .x$CD4, width = 8, height = 6)
  })
  
  # Engagement frequency analysis and CSV output
  engagement_analysis_and_output <- function(cell_type) {
    engagement <- master_clust_Live %>% filter(str_detect(tcell_line, cell_type))
    behav <- engagement %>%
      group_by(cluster, ifelse(contact == 1, "engagerd", "Non_engaged")) %>%
      summarize(frequency = n(), .groups = 'drop') %>%
      mutate(percentage = frequency / sum(frequency) * 100)
    
    write.csv(behav, file.path(output_dir, paste0(cell_type, "_engagement_behavior_freq.csv")), row.names = FALSE)
  }
  
  engagement_analysis_and_output("CD8")
  engagement_analysis_and_output("CD4")
  
  print("Engagement frequency analysis completed and results saved to CSV.")
  print("Simulation completed.")
  dev.off()
  
} else if (n_timepoints == 2) {
  
  # Analysis 1
  ### calculate the frequency of behaviors per conditions for CD8 only:
  
  ### Simulate the selection at t2 hours of T cells that were or were not in contact with organoids at that timepoint
  master_clust_Live2_CD4_CD8<-filter(master_clust_Live, grepl("CD4|CD8", tcell_line))
  
  ### Select cells at timepoint t2  and remove cluster 1 since in the real experiment dead cells are excluded by FACS
  master_clust_Live2_CD8<-subset(master_clust_Live2_CD4_CD8, Time%in%t2 & grepl("CD8", tcell_line) & !cluster=="1") 
  ### For each trackID select the last timepoint in the t2 hours range
  master_clust_Live2_CD8<-master_clust_Live2_CD8%>%group_by(TrackID)%>%filter(Time==max(Time)) 
  
  ## Plot the engager and the non engager populations and their cluster behavior distribution
  Plot_Eng_CD8<-ggplot(master_clust_Live2_CD8, aes(x=factor(contact), fill=as.factor(cluster))) + 
    geom_bar(position="fill")+
    # scale_x_discrete(breaks=c(0,1), minor_breaks = NULL, n.breaks=2, labels = c("True", "False"))+
    labs(y = "Percentage", x = "Contact", title = paste0("Engager & non-engager CD8 \nbetween timepoint ", interval_start, " and ", interval_end))+
    scale_x_discrete(labels=c("No Contact", "Contact"))+
    scale_fill_manual(name = "Behavioral\nsignature", labels = Behavioral_signatures, values=color_palette)+
    theme_minimal()+
    theme(aspect.ratio = 0.2, axis.text.y = element_text())+  
    coord_flip() +
    # facet_grid(interaction(exp_nr, well, organoid_line)  ~ tcell_line) #For plotting all wells separately
    facet_grid(organoid_line  ~ tcell_line)+
    facet_wrap(~tcell_line)
  Plot_Eng_CD8
  ggsave(filename = paste0(output_dir, "Engager & non-engager CD8 between timepoint.png"), plot = Plot_Eng_CD8, width = 8, height = 6)
  
  ### Now plot the super engagers and the never engagers. These cells under go two washing and separation steps:
  ### 1) at t1 hours they are separated into  contacting organoids and non-contacting organoids and then put back in culture. For the non-contacting organoids condition new organoids are added.
  ### 2) at t2 hours they are separated again into contact organoids and non=contacting organoids T cells.
  
  ### Orverall we recover two populations: Super-engagers> T cells that were in contact with organoids at t1 and t2 hours
  ### Never engagers> T cells that were not in contact with organoids at t1 hours and also not at t2 hours.
  
  ### Simulate the selection at t1 hours of T cells that were in contact with an orgnaoid
  master_clust_Live1_EN_CD8<-subset(master_clust_Live2_CD4_CD8, Time%in%t1& grepl("CD8", tcell_line)& contact==1)
  master_clust_Live1_EN_CD8<-master_clust_Live1_EN_CD8%>%group_by(TrackID)%>%filter(Time==min(Time))
  ### Simulate the selection at t1 hours of T cells that were not in contact with an orgnaoid
  master_clust_Live1_NEN_CD8<-subset(master_clust_Live2_CD4_CD8, Time%in%t1& grepl("CD8", tcell_line) &contact==0)
  master_clust_Live1_NEN_CD8<-master_clust_Live1_NEN_CD8%>%group_by(TrackID)%>%filter(Time==min(Time))
  
  ###Simulate the separation at t2 hours of T cells in contact with organoid or not, based on the previously selected cells
  master_clust_Live2_SEN_CD8<-subset(master_clust_Live2_CD8, contact==1 &TrackID%in%master_clust_Live1_EN_CD8$TrackID)
  master_clust_Live2_NEN_CD8<-subset(master_clust_Live2_CD8, contact==0 &TrackID%in%master_clust_Live1_NEN_CD8$TrackID)
  master_clust_Live2_NEN_SEN_CD8<-rbind(master_clust_Live2_NEN_CD8,master_clust_Live2_SEN_CD8)
  master_clust_Live2_NEN_SEN_CD8<-subset(master_clust_Live2_NEN_SEN_CD8, !cluster=="1") ## remove dying cell cluster, since will be deleted by FACS
  
  Plot_SEng_CD8<-ggplot(master_clust_Live2_NEN_SEN_CD8, aes(fill=as.factor(cluster), x=as.factor(contact))) + 
    geom_bar( position="fill")+ 
    labs(y = "Percentage", x = "Contact", title = "Super-engager & Never-engager CD8")+
    scale_x_discrete(labels=c("No Contact", "Contact"))+
    scale_fill_manual(name = "Behavioral\nsignature", labels = Behavioral_signatures, values=color_palette)+
    theme_minimal()+
    theme(aspect.ratio = 0.2, axis.text.y = element_text())+  
    coord_flip() +
    # facet_grid(interaction(exp_nr, well, organoid_line)  ~ tcell_line) #For plotting all wells separately
    facet_grid(organoid_line  ~ tcell_line)
  Plot_SEng_CD8
  ggsave(filename = paste0(output_dir, "Super-engager & Never-engager CD8.png"), plot = Plot_SEng_CD8, width = 8, height = 6)
  
  
  #### Repeat the same for the CD4 cells
  ### calculate the frequency of behaviors per conditions for CD4 only:
  ### Simulate the selection at t2 of T cells that were or were not in contact with organoids at that timepoint
  master_clust_Live2_CD4<-subset(master_clust_Live2_CD4_CD8, Time%in%t2 & grepl("CD4", tcell_line) & !cluster=="1")
  master_clust_Live2_CD4<-master_clust_Live2_CD4%>%group_by(TrackID)%>%filter(Time==max(Time))
  
  Plot_Eng_CD4<-ggplot(master_clust_Live2_CD4, aes(fill=as.factor(cluster), x=as.factor(contact))) + 
    geom_bar( position="fill")+ 
    # scale_x_discrete(breaks=c(0,1), minor_breaks = NULL, n.breaks=2, labels = c("True", "False"))+
    labs(y = "Percentage", x = "Contact", title = paste0("Engager & non-engager CD4 \nbetween timepoint ", interval_start, " and ", interval_end))+
    scale_x_discrete(labels=c("No Contact", "Contact"))+
    scale_fill_manual(name = "Behavioral\nsignature", labels = Behavioral_signatures, values=color_palette)+
    theme_minimal()+
    theme(aspect.ratio = 0.2, axis.text.y = element_text())+  
    coord_flip() +
    # facet_grid(interaction(exp_nr, well, organoid_line)  ~ tcell_line) #For plotting all wells separately
    facet_grid(organoid_line  ~ tcell_line)
  Plot_Eng_CD4
  ggsave(filename = paste0(output_dir, "Engager & non-engager CD4 between timepoint.png"), plot = Plot_Eng_CD4, width = 8, height = 6)
  
  #Analysis 2
  ### Never engagers> T cells that were not in contact with organoids at t1 and also not at t2.
  ### Simulate the selection at t1 of T cells that were in contact with an orgnaoid
  master_clust_Live1_EN_CD4<-subset(master_clust_Live2_CD4_CD8, Time%in%t1& grepl("CD4", tcell_line)& contact==1)
  master_clust_Live1_EN_CD4<-master_clust_Live1_EN_CD4%>%group_by(TrackID)%>%filter(Time==min(Time))
  ### Simulate the selection at t1 of T cells that were not in contact with an orgnaoid
  master_clust_Live1_NEN_CD4<-subset(master_clust_Live2_CD4_CD8, Time%in%t1 & grepl("CD4", tcell_line) &contact==0)
  master_clust_Live1_NEN_CD4<-master_clust_Live1_NEN_CD4%>%group_by(TrackID)%>%filter(Time==min(Time))
  
  ###Simulate the separation at t2 of T cells in contact with organoid or not, based on the previously selected cells
  master_clust_Live2_SEN_CD4<-subset(master_clust_Live2_CD4, contact==1 &TrackID%in%master_clust_Live1_EN_CD4$TrackID)
  master_clust_Live2_NEN_CD4<-subset(master_clust_Live2_CD4, contact==0 &TrackID%in%master_clust_Live1_NEN_CD4$TrackID)
  master_clust_Live2_NEN_SEN_CD4<-rbind(master_clust_Live2_NEN_CD4,master_clust_Live2_SEN_CD4)
  master_clust_Live2_NEN_SEN_CD4<-subset(master_clust_Live2_NEN_SEN_CD4, !cluster=="1") ## remove dying cell cluster, since will be deleted by FACS
  
  Plot_SEng_CD4<-ggplot(master_clust_Live2_NEN_SEN_CD4, aes(fill=as.factor(cluster), x=as.factor(contact))) + 
    geom_bar( position="fill")+ 
    labs(y = "Percentage", x = "Contact", title = "Super-engager & Never-engager CD4")+
    scale_x_discrete(labels=c("No Contact", "Contact"))+
    scale_fill_manual(name = "Behavioral\nsignature", labels = Behavioral_signatures, values=color_palette)+
    theme_minimal()+
    theme(aspect.ratio = 0.2, axis.text.y = element_text())+  
    coord_flip() +
    # facet_grid(interaction(exp_nr, well, organoid_line)  ~ tcell_line) #For plotting all wells separately
    facet_grid(organoid_line  ~ tcell_line)
  Plot_SEng_CD4
  ggsave(filename = paste0(output_dir, "Engager & non-engager CD4 between timepoint.png"), plot = Plot_Eng_CD4, width = 8, height = 6)
  
  # Analysis 3
  ########################### CD8 ########################################################
  ### calculare the proportions for combination with scRNA seq:
  ### create a condition for each experimental condition
  master_clust_Live2_CD8$engagement<-ifelse(master_clust_Live2_CD8$contact==1, "engager","nonengager")
  master_clust_Live2_NEN_SEN_CD8$engagement<-ifelse(master_clust_Live2_NEN_SEN_CD8$contact==1,"super-engaged","never-engaged")
  ## join both datasets:
  CD8_engagement<-rbind(master_clust_Live2_CD8,master_clust_Live2_NEN_SEN_CD8)
  library(dplyr)
  CD8_behav<-CD8_engagement %>% group_by(cluster,engagement) %>%summarize(frec = n(), .groups = "drop")
  engagement_n<-CD8_engagement %>% group_by(engagement) %>%summarize(total_n = n(), .groups = "drop")
  CD8_behav<-left_join(CD8_behav,engagement_n)
  CD8_behav$cluster_prop<-CD8_behav$frec/CD8_behav$total_n
  colnames(CD8_behav)=c("behavioral_cluster", "exp_condition", "count", "total_n", "cluster_proportion")
  
  CD8_engagement$engagement <- factor(CD8_engagement$engagement, levels = c("never-engaged", "nonengager", "engager", "super-engaged"))
  Plot_x<-ggplot(CD8_engagement, aes(fill=as.factor(cluster), x=engagement)) + 
    geom_bar( position="fill")+ 
    labs(y = "Percentage", x = "Contact", title = "Proportions of the four different\n types of engagers CD8")+
    # scale_x_discrete(breaks=labels=c("No Contact", "Contact"))+
    scale_fill_manual(name = "Behavioral\nsignature", labels = Behavioral_signatures, values=color_palette)+
    theme_minimal()+
    theme(aspect.ratio = 0.2, axis.text.y = element_text())+  
    coord_flip() +
    # facet_grid(interaction(exp_nr, well, organoid_line)  ~ tcell_line) #For plotting all wells separately
    facet_grid(organoid_line  ~ tcell_line)
  Plot_x
  ggsave(filename = paste0(output_dir, "Proportions of the four different types of engagers_CD8.png"), plot = Plot_x, width = 8, height = 6)
  
  ### save data to csv same folder where the pseudotime clustering from Farid
  write.csv(CD8_behav,paste0(output_dir, "CD8_engagement_behavior_freq.csv"), row.names = FALSE)
  
  #### Calculate the mean contact time per condition for information (used in the paper)
  ## SE mean contact time:
  master_clust_Live2_CD8_se_pop<-subset(CD8_engagement, engagement=="super-engaged")
  ###calculate the mean engagement time for each signature:
  mean_cont_SE<- master_clust_Live2_CD4_CD8%>%filter(TrackID%in%master_clust_Live2_CD8_se_pop$TrackID)%>%group_by(cluster)%>%
    summarise(contact_per_h= mean(contact)*60)
  mean_cont_SE$exp<-"super-engaged"
  
  
  ## E mean contact time:
  master_clust_Live2_CD8_e_pop<-subset(CD8_engagement, engagement=="engaged")
  ###calculate the mean engagement time for each signature:
  mean_cont_E<- master_clust_Live2_CD4_CD8%>%filter(TrackID%in%master_clust_Live2_CD8_e_pop$TrackID)%>%group_by(cluster)%>%
    summarise(contact_per_h= mean(contact)*60)
  mean_cont_E$exp<-"engaged"
  
  
  ## NonE mean contact time:
  master_clust_Live2_CD8_noe_pop<-subset(CD8_engagement, engagement=="non-engaged")
  ###calculate the mean engagement time for each signature:
  mean_cont_NoE<- master_clust_Live2_CD4_CD8%>%filter(TrackID%in%master_clust_Live2_CD8_noe_pop$TrackID)%>%group_by(cluster)%>%
    summarise(contact_per_h= mean(contact)*60)
  mean_cont_NoE$exp<-"non-engaged"
  
  
  ## No mean contact time:
  master_clust_Live2_CD8_NE_pop<-subset(CD8_engagement, engagement=="never-engaged")
  ###calculate the mean engagement time for each signature:
  mean_cont_NE<- master_clust_Live2_CD4_CD8%>%filter(TrackID%in%master_clust_Live2_CD8_NE_pop$TrackID)%>%group_by(cluster)%>%
    summarise(contact_per_h= mean(contact)*60)
  
  mean_cont_NE$exp<-"never-engaged"
  
  mean_cont<-rbind(mean_cont_NE,mean_cont_NoE,mean_cont_E,mean_cont_SE)
  
  #saveRDS(mean_cont, paste0(output_dir,"min_contact_per_hour_13T_per_exp.rds"))
  
  
  ########################### CD4 ########################################################
  ### calculate the proportions for combination with scRNA seq:
  ### create a condition for each experimental condition
  master_clust_Live2_CD4$engagement<-ifelse(master_clust_Live2_CD4$contact==1, "engager","nonengager")
  master_clust_Live2_NEN_SEN_CD4$engagement<-ifelse(master_clust_Live2_NEN_SEN_CD4$contact==1,"super-engaged","never-engaged")
  
  ## join both datasets:
  CD4_engagement<-rbind(master_clust_Live2_CD4,master_clust_Live2_NEN_SEN_CD4)
  CD4_behav<-CD4_engagement %>% group_by(cluster,engagement) %>%summarize(frec = n(), .groups = "drop")
  engagement_n<-CD4_engagement %>% group_by(engagement) %>%summarize(total_n = n(), .groups = "drop")
  CD4_behav<-left_join(CD4_behav,engagement_n)
  CD4_behav$cluster_prop<-CD4_behav$frec/CD4_behav$total_n
  colnames(CD4_behav)=c("behavioral_signature", "exp_condition", "count", "total_n", "cluster_proportion")
  
  CD4_engagement$engagement <- factor(CD4_engagement$engagement, levels = c("never-engaged", "nonengager", "engager", "super-engaged"))
  
  Plot_y<-ggplot(CD4_engagement, aes(fill=as.factor(cluster), x=engagement)) + 
    geom_bar(position="fill")+ 
    labs(y = "Percentage", x = "Contact", title = "Proportions of the four different\n types of engagers CD4")+
    # scale_x_discrete(breaks=labels=c("No Contact", "Contact"))+
    scale_fill_manual(name = "Behavioral\nsignature", labels = Behavioral_signatures, values=color_palette)+
    theme_minimal()+
    theme(aspect.ratio = 0.2, axis.text.y = element_text())+  
    coord_flip() +
    facet_wrap(~tcell_line)+
    # facet_grid(interaction(exp_nr, well, organoid_line)  ~ tcell_line) #For plotting all wells separately
    facet_grid(organoid_line  ~ tcell_line)
  Plot_y
  ggsave(filename = paste0(output_dir, "Proportions of the four different types of engagers_CD4.png"), plot = Plot_y, width = 8, height = 6)
  ### save data to csv same folder where the pseudotime clustering from Farid
  write.csv(CD4_behav,paste0(output_dir,"CD4_engagement_behavior_freq.csv"), row.names = FALSE)
  
  #pdf(file = paste0(output_dir,"Engager_Super_ENG_proportions_CD4_CD8.pdf"))
  #Plot_x
  #Plot_y
  
  print("Engagement frequency analysis completed and results saved to CSV.")
  print("Simulation completed.")
}

