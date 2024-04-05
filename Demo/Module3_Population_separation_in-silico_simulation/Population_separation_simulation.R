library(tidyverse)
library(fs)
library(furrr)
library(ggplot2)
library(patchwork)
library(yaml)


### Set to TRUE if you want to run import and processing even if file already exists
force_redo=TRUE
tracks_provided=NULL

### Checks if being run in GUI (e.g. Rstudio) or command line
if (interactive()) {
  ### !!!!!! Change the path to the BGT_config file here if running the code in RStudio !!!!!!
  ### Demo path
  BGT_dir = paste0(dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path))),"/")
  pars = yaml.load_file("~/BGT/Demo/Module3_Population_separation_in-silico_simulation/config_template.yml")

  ### For your own file, uncomment following line and add own path to the BEHAV3D_config.yml
  # pars = yaml.load_file("")
  
} else {
  option_list = list(
    make_option(c("-c", "--config"), type="character", default=NULL, 
                help="Path to the BGT config file", metavar="character"),
    make_option(c("-f", "--force_redo"), action="store_true", default=FALSE, 
                help="Force the pipeline to re-import data even if files exists"),
    make_option(c("-t", "--tracks_rds"), type="character", default=NULL, 
                help="(Optional) Path to RDS file containing processed T cell track data", metavar="character")
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

### Setting data directory (if specified) and creating output directories
output_dir=paste0(BGT_dir, "/Results/Module3_Population_seperation_in-silico_simulation/")
dir.create(output_dir, recursive=TRUE)


n_timepoints = pars$n_timepoints
first_timepoint_interval_starts = pars$first_timepoint_interval_starts
first_timepoint_interval_ends = pars$first_timepoint_interval_ends
reoccuring_time = pars$reoccuring_time


 
# Set up parallel processing
plan(multisession, workers = 4)

# Load data
master_clust_Live <- readRDS(pars$classified_tcell_track_data_filepath_rds)

# Generate timepoint ranges
time_intervals <- list()
for (i in 1:n_timepoints) {
  interval_start <- first_timepoint_interval_starts + (i - 1) * reoccuring_time
  interval_end <- first_timepoint_interval_ends + (i - 1) * reoccuring_time
  time_intervals[[i]] <- c(interval_start, interval_end)
}
  
  # Define color palette
color_palette <- c("gold", "darkolivegreen", "seagreen", "forestgreen", "dodgerblue", "cyan", "indianred", "firebrick", "brown")
  
  # Process data and plot for each time interval
plots <- future_map(time_intervals, function(interval) {
  subset_data <- master_clust_Live %>%
    filter(Time >= interval[1], Time <= interval[2])
    
    cd8_data <- subset_data %>% filter(str_detect(tcell_line, "CD8"), cluster != "1")
    cd4_data <- subset_data %>% filter(str_detect(tcell_line, "CD4"), cluster != "1")
    
    plot_cd8 <- ggplot(cd8_data, aes(x = factor(contact), fill = factor(cluster))) +
      geom_bar(position = "fill") +
      scale_fill_manual(values = color_palette) +
      labs(y = "Percentage", x = "Contact", title = sprintf("CD8 Engagement between timepoint %s", paste(interval, collapse = "-"))) +
      theme_minimal() +
      coord_flip() +
      facet_wrap(~tcell_line)
    print(plot_cd8)
    
    plot_cd4 <- ggplot(cd4_data, aes(x = factor(contact), fill = factor(cluster))) +
      geom_bar(position = "fill") +
      scale_fill_manual(values = color_palette) +
      labs(y = "Percentage", x = "Contact", title = sprintf("CD4 Engagement between timepoint %s", paste(interval, collapse = "-"))) +
      theme_minimal() +
      coord_flip() +
      facet_wrap(~tcell_line)
    print(plot_cd4)
    
    list(CD8 = plot_cd8, CD4 = plot_cd4)
  })
  
# Save plots
walk2(plots, seq_along(plots), ~{
  ggsave(filename = file.path(output_dir, sprintf("plot_cd8_timepoint_%s.png", .y)), plot = .x$CD8, width = 8, height = 6)
  ggsave(filename = file.path(output_dir, sprintf("plot_cd4_timepoint_%s.png", .y)), plot = .x$CD4, width = 8, height = 6)
})
  
# Engagement frequency analysis and CSV output
engagement_analysis_and_output <- function(cell_type) {
  engagement <- master_clust_Live %>% filter(str_detect(tcell_line, cell_type))
  behav <- engagement %>%
    group_by(cluster, ifelse(contact == 1, "engager", "nonengager")) %>%
    summarize(frequency = n(), .groups = 'drop') %>%
    mutate(percentage = frequency / sum(frequency) * 100)
    
  write.csv(behav, file.path(output_dir, paste0(cell_type, "_engagement_behavior_freq.csv")), row.names = FALSE)
  }
  
engagement_analysis_and_output("CD8")
engagement_analysis_and_output("CD4")
  
print("Engagement frequency analysis completed and results saved to CSV.")
print("Simulation completed.")

