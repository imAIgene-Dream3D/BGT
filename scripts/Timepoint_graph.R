library(dplyr)
library(ggplot2)
library(yaml)
library(optparse)


### Set to TRUE if you want to run import and processing even if file already exists
force_redo=TRUE
tracks_provided=NULL

BGT_dir <- NULL

### Checks if being run in GUI (e.g. Rstudio) or command line
if (interactive()) {
  ### !!!!!! Change the path to the BGT_config file here if running the code in RStudio !!!!!!
  ### Demo path
  BGT_dir = paste0(dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path))),"/")
  pars = yaml.load_file("~/BGT/Demo/Module1_Evalution_of_Super_engager_population_dynamics_in_co_culture/config_template.yml")
} else {
  ### Define a default BGT_dir for non-interactive sessions if not set
  ### You may adjust this path as necessary for your non-interactive environment
  if(is.null(BGT_dir)) {
    BGT_dir <- "~/BGT/"  # Adjust this path as necessary
  }
  
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
output_dir=paste0(BGT_dir, "/Results/Module1_Evaluation_of_Super_engager_population_dynamics_in_co_culture/")
dir.create(output_dir, recursive=TRUE)

imaging_time = pars$imaging_time
# Load data
master_clust_Live <- readRDS(pars$classified_tcell_track_data_filepath_rds)

# Data manipulations
#Group data by T cell line , cluster and time, then calculate  the count of records in each group
# Subsequently, regroups by t cell line and time to calculate the percentage of each group over the total count per timepoint.
Raw_per_overall<-master_clust_Live %>%
  group_by(tcell_line,cluster, Time) %>%
  summarise(n = n()) %>% 
  filter(n() > 5) %>%
  group_by(tcell_line, Time) %>%
  mutate(perc = n*100 / sum(n))

# Time conversion
Raw_per_overall$Time <- Raw_per_overall$Time/30

# filter out groups with fewer than 5 cells
Raw_per_overall <- Raw_per_overall %>% filter(n>5) ##at least 5 cells


# filter for imaging time, remove all the points before imaging time
Raw_per_overall <- Raw_per_overall %>%
  filter(Time >= imaging_time) 


### subset cluster 9:
Raw_per_overall <- subset(Raw_per_overall,cluster=="9")

# Further subset the data for specific t cell lines
Raw_per_overall <- subset(Raw_per_overall,tcell_line %in% c("CD4_TEG","CD8_TEG"))

imaging_time <- min(Raw_per_overall$Time) 

# Assuming Raw_per_overall is already loaded and prepared
max_points_time <- Raw_per_overall %>%
  group_by(tcell_line) %>%
  summarise(MaxTime = Time[which.max(perc)], MaxPerc = max(perc)) %>%
  ungroup()

mean_max_time <- mean(max_points_time$MaxTime)

# Calculate the window
window_start <- mean_max_time - (15 / 30)
window_end <- mean_max_time + (15 / 30)

# Filter data within this window
windowed_data <- Raw_per_overall %>%
  filter(Time >= window_start & Time <= window_end)


Per2 <- ggplot(Raw_per_overall, aes(Time, perc, group = tcell_line, color = as.factor(tcell_line))) +
  geom_smooth(size = 1, span = 0.6) +
  geom_rect(aes(xmin = window_start, xmax = window_end, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.005, inherit.aes = FALSE) +
  #geom_line() +
  geom_vline(xintercept = c(window_start, window_end), color = "black") +
  theme_bw() + 
  ylab("% T cells in cluster 9 (super-engagers)") +
  xlab("Time in Co-culture(Hours)") +
  scale_x_continuous(breaks = c(imaging_time, seq(from = ceiling(imaging_time), to = max(Raw_per_overall$Time), by = 1))) + # Ensure this line is properly continued with a '+'
  scale_fill_manual(values = c(
    "CD4_TEG" = "darkolivegreen3",
    "CD8_TEG" = "dodgerblue" 
  ))+
  scale_color_manual(values = c(
    CD4_TEG = "darkolivegreen3", 
    CD8_TEG = "dodgerblue" 
   
  )) +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 30, l = 10, unit = "pt"),  
    axis.text.x = element_text(vjust = 0.5, size=10), 
    axis.text.y = element_text(size = 10), 
    axis.title.y = element_text(size = 12), 
    axis.title.x = element_text(size = 12), 
    legend.text = element_text(size = 10),
    aspect.ratio = 1
  ) +
  labs(color = "T Cell Type", 
       fill = "T Cell Type"  
  )+
  ggtitle("T-cell engagement at different time points ") 
ggsave(filename = paste0(output_dir, "T-cell_engagementVsTime.png"), width = 8, height = 6)
write.csv(Raw_per_overall, file.path(paste0(output_dir, "T-cell_engagementVsTime.csv")), row.names = FALSE)
Per2
