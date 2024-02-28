#### Inputs
imaging_time = 1
tenT_data = '~/Desktop/working folder/NaPo_bgt/tcell_behavior/results/classified_tcell_track_data.rds'
thirteenT_data = "~/Desktop/working folder/NaPo_bgt/13T_processing/classified_tcell_track_data.rds"


library(dplyr)
library(ggplot2)
#read data
master_clust_Live6 = readRDS(file=thirteenT_data)

# Data manipulations
#Group data by T cell line , cluster and time, then calculate  the count of records in each group
# Subsequently, regroups by t cell line and time to calculate the percentage of each group over the total count per timepoint.
Raw_per_overall<-master_clust_Live6 %>%
  group_by(tcell_line,cluster, Time) %>%
  summarise(n = n()) %>% 
  group_by(tcell_line, Time) %>%
  mutate(perc = n*100 / sum(n))

# Time conversion
Raw_per_overall$Time <- Raw_per_overall$Time/30

# filter out groups with fewer than 5 cells
Raw_per_overall <- Raw_per_overall %>%
  filter(n>5) ##at least 5 cells

# filter for imaging time, remove all the points before imaging time
Raw_per_overall <- Raw_per_overall %>%
  filter(Time >= imaging_time) 


### subset cluster 9:
Raw_per_overall <- subset(Raw_per_overall,cluster=="9")

# Further subset the data for specific t cell lines
Raw_per_overall <- subset(Raw_per_overall,tcell_line %in% c("CD4_TEG","CD8_TEG"))

# Assuming imaging_time is defined. If not, replace it with an actual start time from your dataset.
imaging_time <- min(Raw_per_overall$Time) # Example definition; adjust as necessary

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
  geom_smooth(size = 1, span = 0.5) +
  geom_rect(aes(xmin = window_start, xmax = window_end, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.005, inherit.aes = FALSE) +
  #geom_line() +
  geom_vline(xintercept = c(window_start, window_end), color = "black") +
  theme_bw() + # Duplicate theme_bw() removed
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
    # Add more mappings for each tcell_line as needed
  )) +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 30, l = 10, unit = "pt"),  # Adjust the right margin to ensure there's enough space
    axis.text.x = element_text(vjust = 0.5, size=10), 
    axis.text.y = element_text(size = 10), 
    axis.title.y = element_text(size = 12), 
    axis.title.x = element_text(size = 12), 
    legend.text = element_text(size = 10),
    aspect.ratio = 1
  ) +
  labs(color = "T Cell Type",  # Set the legend title for color
       fill = "T Cell Type"  # If you also have a fill legend, set its title (optional)
  )+
  ggtitle("Frequency of T-cells in the co-culture experiment") # Corrected typo in the comment

Per2
