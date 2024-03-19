#In this in-silico simulation for T cell population separation into two predominant experimental engagement states enagegers and non engagers,
# which is collectively called as t cell contacts with the Tumor organoid. 
# The purpose of this simulation is to recreate a T cell and tumor organoid co culture experiment and,
# for the different user defined time intervals find out the behavioral signatures frequencies of the t cells in engaged and non engaged states.

population_separation_simulation <- function(working_directory_path, classified_tcell_track_data_filepath_rds, n_timepoints, first_timepoint_intervals, reoccuring_time) {
  library(tidyverse) 
  library(fs)
  library(furrr)
  library(ggplot2)
  library(patchwork)
  
  # Set up parallel processing
  plan(multisession, workers = 4)
  
  # Create necessary directories
  dirs <- c("Behaviour_signature_simulation_only", "Results", "Probability_Map", "Metadata")
  dirs_paths <- file.path(working_directory_path, dirs)
  walk(dirs_paths, ~if (!dir_exists(.x)) dir_create(.x))
  
  # Load data: Classified T cell Track data From Behav3D
  master_clust_Live <- readRDS(classified_tcell_track_data_filepath_rds)
  
  # From the %T cell engagement vs Time graph we identified the first targeted time point to find out
  # the Behavioral signatures frequency for each Experimental Engagement population. 
  # After determining the first time point you can setup the reoccurring time, which will automatically calculated the the next time point interval for which the following estimation is to be done, as it was done for the first time point.  
  # Generate time point ranges
  time_intervals <- list()
  for (i in 1:n_timepoints) {
    interval_start <- first_timepoint_intervals[1] + (i - 1) * reoccuring_time
    interval_end <- first_timepoint_intervals[2] + (i - 1) * reoccuring_time
    time_intervals[[i]] <- c(interval_start, interval_end)
  }
  
  # Define color palette
  color_palette <- c("gold", "darkolivegreen", "seagreen", "forestgreen", "dodgerblue", "cyan", "indianred", "firebrick", "brown")
  Behavioral_signatures_scale_values <- c("Dead", "Static", "Lazy", "Slow scanner", "Medium scanner", "Super scanner", "Tickler", "Engager", "Super-engager")
  
  # Process data and plot for each time interval
  plots <- future_map(time_intervals, function(interval) {
    subset_data <- master_clust_Live %>%
      filter(Time >= interval[1], Time <= interval[2])
    
    cd8_data <- subset_data %>% filter(str_detect(tcell_line, "CD8"), cluster != "1")
    cd4_data <- subset_data %>% filter(str_detect(tcell_line, "CD4"), cluster != "1")
   
    cd8_data$cluster <- factor(cd8_data$cluster, levels = 1:9, labels = Behavioral_signatures_scale_values)
    cd4_data$cluster <- factor(cd4_data$cluster, levels = 1:9, labels = Behavioral_signatures_scale_values)
     
    plot_cd8 <- ggplot(cd8_data, aes(x = factor(contact), fill = factor(cluster))) +
      geom_bar(position = "fill") +
      scale_fill_manual(values = color_palette, name="Behavioral Signatures") +
      labs(y = "Percentage", x = "CD8 T-cell Contact", title = sprintf("CD8 Engagement between timepoint %s minutes", paste(interval, collapse = "-"))) +
      theme_minimal() +
      coord_flip() +
      facet_wrap(~tcell_line)
    print(plot_cd8)
    
    plot_cd4 <- ggplot(cd4_data, aes(x = factor(contact), fill = factor(cluster))) +
      geom_bar(position = "fill") +
      scale_fill_manual(values = color_palette, name="Behavioral Signatures") +
      labs(y = "Percentage", x = "CD4 T-cell Contact", title = sprintf("CD4 Engagement between timepoint %s minutes", paste(interval, collapse = "-"))) +
      theme_minimal() +
      coord_flip() +
      facet_wrap(~tcell_line)
    print(plot_cd4)
    
    list(CD8 = plot_cd8, CD4 = plot_cd4)
  })
  
  # Save plots
  walk2(plots, seq_along(plots), ~{
    ggsave(filename = file.path(working_directory_path, "Results", sprintf("plot_cd8_timepoint_%s.png", .y)), plot = .x$CD8, width = 8, height = 6)
    ggsave(filename = file.path(working_directory_path, "Results", sprintf("plot_cd4_timepoint_%s.png", .y)), plot = .x$CD4, width = 8, height = 6)
  })
  
  # Engagement frequency analysis and CSV output
  engagement_analysis_and_output <- function(cell_type) {
    engagement <- master_clust_Live %>% filter(str_detect(tcell_line, cell_type))
    behav <- engagement %>%
      group_by(cluster, ifelse(contact == 1, "engager", "nonengager")) %>%
      summarize(frequency = n(), .groups = 'drop') %>%
      mutate(percentage = frequency / sum(frequency) * 100)
    
    write.csv(behav, file.path(working_directory_path, paste0("/Results/", cell_type, "_engagement_behavior_freq.csv")), row.names = FALSE)
  }
  
  engagement_analysis_and_output("CD8")
  engagement_analysis_and_output("CD4")
  
  print("Engagement frequency analysis completed and results saved to CSV.")
  print("Simulation completed.")
}

# Example usage of the function
population_separation_simulation("~/Desktop/working folder/NaPo_bgt/13T_processing/Test",
                                 "~/Desktop/working folder/NaPo_bgt/13T_processing/classified_tcell_track_data.rds",
                                 4,
                                 c(86, 95),
                                 30)
