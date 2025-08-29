# - - - Masterscript for FlowCode analysis  
# with log2 transformation, and accounting for clone skewing, and minimum detection skewing, 
# written for 3 sets of CD8 data, where the particular tissue of interest was the pancreas
# can easily be modified to input other data
# requires percell, pct_matrix, and thresholds .csv files from the debarcoding app output

# I've written a guide to this with explanations and examples

# The script measures phenotypic change by barcode, to process and make graphs of the debarcoded data
# I've added in my code to visually check if there is clone expansion which may pose an issue/skew data.
# if that shows there is, then
# code is altered to remove the overexpanded cd8 clone from affected mouse
# by an automated method for identifying affected clones with editable criteria
# there are also options that change the design and features of the graphs, like how many rows of heatmaps and if you want a log2fc threshold on volcanos etc
 
#option for non parametric testing has been added, for volcano plots and heatmaps when comparing tissues

# last edited 28th August (Natalie) - tidying up the 20250806 version


#new since last version  (for the 20250806)
# it saved the set_guides wrong, I've fixed it
# Also the ability to group tissues by lymphoid, non lymphoid, and gut, and then compare between them as volcanos, is now there
# fixed the shapiro section from crashing
# added the t-test depletion/enrichment heatmap in a tissue of choice, comparing the activated and naive groups, to those groups in the spleen
# and barcharts showing the depletion of guides in a tissue of interest as compared to the spleen 
# alphabetical option to organise target guides in heatmap
# min_det of all tissues to be twice spleen is now an option
# added the neat layout for many sets in the clone heatmaps

# Remember to manually change the titles of graphs to describe what your graph shows
# For different FlowCode experimental designs and numbers of libraries, this example usage of the script needs to be adjusted, as per the comments





required.packages <- c("Rtsne", "ggplot2", "RColorBrewer", "dplyr", "emmeans",  
                       "EmbedSOM", "tidyr", "coda", "scattermore",
                       "data.table", "ggrepel", "pheatmap", "parallelly",
                       "purrr", "lsa", "stringr", "ineq" )

for (req.package in required.packages){
  if(!requireNamespace(req.package, quietly=TRUE)){
    install.packages(req.package, repos='http://cran.us.r-project.org')
  }
}

bioconductor.packages <- c("ConsensusClusterPlus")

for (req.package in bioconductor.packages){
  if(!requireNamespace(req.package, quietly=TRUE)){
    BiocManager::install(req.package)
  }
}

invisible( lapply( c(required.packages, bioconductor.packages), library, character.only = TRUE ) )

# tip: set the seed to be today's date in 8-digit format.
# this will name the output files as well as setting clustering reproducibility
flowcode.seed <- 20250806
flowcode.threads <- parallelly::availableCores() - 1
data.table::setDTthreads(flowcode.threads)

# import all perCell data
# you can do it set by set (if you have multiple libraries) and then combine it at this stage 
perCellData<-read.csv("CD8_combined_perCellData.csv")

#heres the set by set

#importing each set
#set1<-read.csv("CD8_Set1_FlowcodeDecoder_20250303_perCell_data.csv")
#set2<-read.csv("CD8_Set2_FlowcodeDecoder_20250303_perCell_data.csv")
#set3<-read.csv("CD8_Set3_FlowcodeDecoder_20250303_perCell_data.csv")

#combining all sets to a perCell data
#perCellData <- rbind(set1,set2,set3)

#cleaning up the perCell Data
non.procode.cells <- c( "no procode signal", "1 or 2 unexpected signals","3 or more unexpected signals")

perCellData$Procode_combination[perCellData$Procode_combination==""] <- "Untransduced"
perCellData$Id[perCellData$Id %in% non.procode.cells] <- "Other"

# add source and individual columns
# this will need to be modified based on your naming convention
perCellData <- perCellData %>% 
  separate(sample, remove = TRUE, sep = "[ _]", into = c("X1", "tissue", "X3", "X4", "X5", "X6",
                                                         "mouse", "X8", "X9" ) ) %>%
  select(-X1, -X3, -X4, -X5, -X6, -X8, -X9)

#  remove non-working guides here from the perCell ###-------------
# you must look at the pct_matrix and see which guides have 0s for all samples
# then remove them 

# Filter gene names that appear less than 20 times
# this step does not identify guides that completely failed, that have 0 occurrences
# you may want to alter this threshold, for example if a guide has 30 but none in the spleen, you are limited in future comparisons with it
gene_counts <- table(perCellData$Id)

rare_genes <- names(gene_counts[gene_counts < 20])

cat("The following genes have below 20 cells across all tissues and mice:\n")
print(rare_genes)


# input guides you consider non-working, the <20 ones and manually identified zero ones
perCellData <- perCellData %>%
  filter(!str_detect(Id, "APC|Celsr2|Ctnna2|Cd69|Itga3|Siglece"))

any(str_detect(perCellData$Id, "APC|Celsr2|Ctnna2|Cd69|Itga3|Siglece"))




# import in each set's pct matrix, will need to be processed seperately then combined after
if ( PREPARE_PCT_MATRIX <- TRUE ) {


#uses the pct_matrix
# ----------------  set1
set1<- read.csv("CD8_Set1_FlowcodeDecoder_20250303_pct_matrix.csv")

set1 <- set1 %>% 
  separate(sample, remove = TRUE, sep = "[ _]", into = c("X1", "tissue", "X3", "X4", "X5",
                                                         "X6", "mouse", "X8", "X9") ) %>%
  select(-X1, -X3, -X4, -X5, -X6, -X8, -X9)

set1 <- set1 %>% select(-no.procode.signal, -X1.or.2.unexpected.signals, -X3.or.more.unexpected.signals ) 

# remove non-working guides here 
set1 <- set1 %>% select( -APC, -Celsr2, -Ctnna2, -Cd69)

# correct frequencies to be among transduced cells only
set1$total <- rowSums(set1[,-c(1,2)])

## Important: remove any rows with total = 0 ###
set1 <- set1 %>% filter(total != 0)

set1 <- set1 %>% mutate(across(where(is.numeric), ~ . * 100 / total))
anyNA(set1)

set1 <- set1 %>% 
  pivot_longer(cols = -c(tissue, mouse, total), names_to = "target_gene", values_to = "freq" )

#get the info you need for later
# Save unique mouse entries into set1_mice
set1_mice <- unique(set1$mouse)

# Save all column names except the first as its ('mouse') into set1_guides
# version 20250714 fixed this line
set1_guides <- unique(set1$target_gene)




#  --- repeat for set 2
set2<- read.csv("CD8_Set2_FlowcodeDecoder_20250303_pct_matrix.csv")

set2 <- set2 %>% 
  separate(sample, remove = TRUE, sep = "[ _]", into = c("X1", "tissue", "X3", "X4", "X5",
                                                         "X6", "mouse", "X8", "X9") ) %>%
  select(-X1, -X3, -X4, -X5, -X6, -X8, -X9)

set2 <- set2 %>% select(-no.procode.signal, -X1.or.2.unexpected.signals, -X3.or.more.unexpected.signals ) 

# remove non-working guides here 
set2 <- set2 %>% select( -Itga3)

# correct frequencies to be among transduced cells only
set2$total <- rowSums(set2[,-c(1,2)])

## Important: remove any rows with total = 0 ###
set2 <- set2 %>% filter(total != 0)

set2 <- set2 %>% mutate(across(where(is.numeric), ~ . * 100 / total))
anyNA(set2)

set2 <- set2 %>% 
  pivot_longer(cols = -c(tissue, mouse, total), names_to = "target_gene", values_to = "freq" )


#get the info you need
# Save unique mouse entries into set2_mice
set2_mice <- unique(set2$mouse)

# Save all column names except the first as its 'mouse' into set2_guides
set2_guides <- unique(set2$target_gene)


# ---  repeat for set 3
set3<- read.csv("CD8_Set3_FlowcodeDecoder_20250303_pct_matrix.csv")

set3 <- set3 %>% 
  separate(sample, remove = TRUE, sep = "[ _]", into = c("X1", "tissue", "X3", "X4", "X5",
                                                         "X6", "mouse", "X8", "X9") ) %>%
  select(-X1, -X3, -X4, -X5, -X6, -X8, -X9)

set3 <- set3 %>% select(-no.procode.signal, -X1.or.2.unexpected.signals, -X3.or.more.unexpected.signals ) 

# remove non-working guides here 
set3 <- set3 %>% select( -Siglece)

# correct frequencies to be among transduced cells only
set3$total <- rowSums(set3[,-c(1,2)])

## Important: remove any rows with total = 0 ###
set3 <- set3 %>% filter(total != 0)

set3 <- set3 %>% mutate(across(where(is.numeric), ~ . * 100 / total))
anyNA(set3)

set3 <- set3 %>% 
  pivot_longer(cols = -c(tissue, mouse, total), names_to = "target_gene", values_to = "freq" )


#get the info you need
# Save unique mouse entries into set3_mice
set3_mice <- unique(set3$mouse)

# Save all column names except the first as its 'mouse' into set3_guides
set3_guides <- unique(set3$target_gene)


#combine them
pct_matrix <- rbind(set1,set2,set3)
}
#### end of pct_matrix making section #####


#### clone heatmaps ####

#Change this between TRUE or FALSE depending on if you want the section to run.
if ( PLOT_HEATMAPS_TO_DETECT_TCR_EXPANSION <- TRUE ) {
  
#### Optional - Makes heatmap-style graphs to visually identify if there is clonal expansion in samples ####
#  --- now plots the graphs, this is pancreas and spleen, you can change however you'd like
# whichever tissues you want to identify potentially affected clones. 
# if you have multiple sets there is a neater option that gets rid of the blank space

### option 1 - for if you have one set, or if you have multiple sets and dont mind blank space on the heatmap ###
#pancreas

pancreas_pcts <- pct_matrix %>%
  filter(tissue=="Pancreas")

# Pivot to wide format
pan_matrix <- pancreas_pcts %>%
  select(-tissue) %>%
  select(-total) %>%
  pivot_wider(names_from = target_gene, values_from = freq, values_fill = 0)

# Convert to long format for plotting
pan_long <- pan_matrix %>%
  pivot_longer(cols = -mouse, names_to = "target_gene", values_to = "cell_count")

# Plot the heatmap
pan <- ggplot(pan_long, aes(x = target_gene, y = mouse, fill = cell_count)) +
  geom_tile(color = "grey85") +
  scale_fill_gradient(low = "white", high = "purple4") +
  labs(title = "Heatmap of CD8 pct_matrix (Pancreas) per Guide and Mouse",
       x = "Target Gene",
       y = "Mouse",
       fill = "pct frequency") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 12)
  )

# Display the plot
print(pan)

# Save the plot
ggsave("CD8_Pancreas_pct_heatmap.png",
       plot = pan,
       width = 10, height = 6, dpi = 300, bg = "white")


# spleen 

spleen_pcts <- pct_matrix %>%
  filter(tissue=="Spleen")

# Pivot to wide format
spl_matrix <- spleen_pcts %>%
  select(-tissue) %>%
  select(-total) %>%
  pivot_wider(names_from = target_gene, values_from = freq, values_fill = 0)

# Convert to long format for plotting
spl_long <- spl_matrix %>%
  pivot_longer(cols = -mouse, names_to = "target_gene", values_to = "cell_count")

# Plot the heatmap
spl <- ggplot(spl_long, aes(x = target_gene, y = mouse, fill = cell_count)) +
  geom_tile(color = "grey85") +
  scale_fill_gradient(low = "white", high = "maroon4") +
  labs(title = "Heatmap of CD8 pct_matrix (Spleen) per Guide and Mouse",
       x = "Target Gene",
       y = "Mouse",
       fill = "pct frequency") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 12)
  )

# Display the plot
print(spl)

# Save the plot
ggsave("CD8_Spleen_pct_heatmap.png",
       plot = spl,
       width = 10, height = 6, dpi = 300, bg = "white")

#### option 2 - if you have multiple sets and want them presented neatly ####
pancreas_pcts <- pct_matrix %>%
  filter(tissue=="Pancreas")

#split into each set
set1_pcts <- pancreas_pcts[pancreas_pcts$mouse %in% set1_mice, ]
set2_pcts <- pancreas_pcts[pancreas_pcts$mouse %in% set2_mice, ]
set3_pcts <- pancreas_pcts[pancreas_pcts$mouse %in% set3_mice, ]

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# --- Step 1: Clean and combine datasets ---
set1_clean <- set1_pcts 
set2_clean <- set2_pcts
set3_clean <- set3_pcts 

# Combine all for global color scale
all_pcts <- bind_rows(set1_clean, set2_clean, set3_clean)

# --- Step 2: Get global min and max for consistent color scale ---
global_min <- min(all_pcts$freq, na.rm = TRUE)
#global_max <- max(all_pcts$freq, na.rm = TRUE)

#you can set the max of the colour scale manually to match for example a different t cell type's heatmap you made
# or you can just make it the max pct from this data
global_max <- 85

# --- Step 3: Function to create heatmap ---

create_heatmap <- function(data, title = NULL, show_x_label = TRUE) {
  pan_matrix <- data %>%
    select(-tissue) %>%
    pivot_wider(names_from = target_gene, values_from = freq, values_fill = 0)
  
  pan_long <- pan_matrix %>%
    pivot_longer(cols = -mouse, names_to = "target_gene", values_to = "cell_count")
  
  p <- ggplot(pan_long, aes(x = target_gene, y = mouse, fill = cell_count)) +
    geom_tile(color = "grey85") +
    scale_fill_gradient(low = "white", high = "#048077", limits = c(global_min, global_max)) +
    labs(
      title = title,
      x = if (show_x_label) "Target Gene" else NULL,
      y = "Mouse",
      fill = "% frequency"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      text = element_text(size = 12)
    )
  
  return(p)
}


# --- Step 4: Create individual heatmaps ---
pan1 <- create_heatmap(set1_clean, title = "Sets: CD8 pct_matrix (Pancreas)", show_x_label = FALSE)
pan2 <- create_heatmap(set2_clean, title = NULL, show_x_label = FALSE)
pan3 <- create_heatmap(set3_clean, title = NULL, show_x_label = TRUE)

# --- Step 5: Stack and save ---
combined_plot <- pan1 / pan2 / pan3

combined_plot

ggsave("sets_heatmap.png", plot = combined_plot, width = 25, height = 20, units = "cm", dpi = 300)

#### end of option 2 ####

### for any option continue from here ##
# ---  bar chart of total cell counts per guide, adds up from all the mice.
#may be informative  in terms of cell counts rather than percentages

pancreas_perCells <-perCellData %>% 
  filter (tissue == "Pancreas") %>%
  filter (Id != "Other")

spleen_perCells <-perCellData %>% 
  filter (tissue == "Spleen") %>%
  filter (Id != "Other")


# # -------  plot pancreas

# Count number of cells per Id
cell_counts <- pancreas_perCells %>%
  count(Id)

# Create the plot
cell_plot <- ggplot(cell_counts, aes(x = Id, y = n)) +
  geom_bar(stat = "identity", fill = "orchid") + # Customize bar color
  labs(title = "Number of CD8s per FlowCode in the Pancreas",
       x = "Guide",
       y = "Number of CD8s") +
  theme_minimal() +
  theme(
    text = element_text(size = 14), # Customize text size
    axis.text.x = element_text(angle = 45, hjust = 1, size=7),
    axis.line = element_line(color = "black") # Add axis lines
    
  )

# Display the plot
print(cell_plot)

# Save the plot as a PNG
ggsave("CD8s_per_id_pancreas.png", plot = cell_plot, width = 10, height = 4, dpi = 300, bg = "white")


# --------  plot spleen
cell_counts <- spleen_perCells %>%
  count(Id)

# Create the plot
cell_plot <- ggplot(cell_counts, aes(x = Id, y = n)) +
  geom_bar(stat = "identity", fill = "turquoise3") + # Customize bar color
  labs(title = "Number of CD8s per FlowCode in the Spleen",
       x = "Guide",
       y = "Number of CD8s") +
  theme_minimal() +
  theme(
    text = element_text(size = 14), # Customize text size
    axis.text.x = element_text(angle = 45, hjust = 1, size=7),
    axis.line = element_line(color = "black") # Add axis lines
    
  )

# Display the plot
print(cell_plot)

# Save the plot as a PNG
ggsave("CD8s_per_id_spleen.png", plot = cell_plot, width = 10, height = 4, dpi = 300, bg = "white")

}
 

#### end of clone graphs ####




#### autothreshold remove clones ####

# ----- here is an added part as compared to the original code -----
#automated detection and removal of clones expanded in the pancreas, can be altered to whichever tissue you want.
# can be used if the heatmap graph above suggests that there was clonal expansion that takes over a sample, that you want to remove
# this block of code does require the above heatmap code to be run, at least the part that loads the pct_matrices

#if you want to skip this, start using the code again at "rejoin code here if not removing clones"

if ( REMOVE_BIG_TCR_CLONES_AUTOMATICALLY <- FALSE ) {

# however in instances where distribution is non normal, or clonal expansion is posing an issue that this is insufficient to solve,
# it is likely better to use non-paramtetric testing - instead of comparing pcts of each guide per mouse, you rank them and compare the rankings
# there is code a little later in the script that checks which method is more appropriate for your data
# and there is the code to run non-parametric statistics


perCellData <- perCellData %>%
  filter (Id != "Other")

# calculate median values for each guide in the pancreas


pancreas_perCells <-perCellData %>% 
  filter (tissue == "Pancreas")

spleen_perCells <-perCellData %>% 
  filter (tissue == "Spleen")

pancreas_summary <- pancreas_perCells %>%
  group_by(mouse, Id) %>%
  summarise(pan_count = n(), .groups = "drop")

pancreas_summary <- pancreas_summary %>% 
  group_by(mouse) %>% 
  mutate( total = sum( pan_count ))

pancreas_summary <- pancreas_summary %>% 
  group_by(mouse) %>% 
  mutate( freq = pan_count / total* 100 )

pancreas_summary <- pancreas_summary %>% 
  group_by(Id) %>% 
  mutate( median_freq = median(freq) )

pancreas_summary <- pancreas_summary %>% 
  group_by(Id) %>% 
  mutate( sd_freq = sd(freq) )

#the selection criteria for guides to remove, can be adjusted as needed
to_remove <- pancreas_summary %>% 
  group_by(Id) %>% 
  filter( freq > median_freq + sd_freq) %>% # could be 3*median
  filter( pan_count > 200 ) %>% # could leave this line out
  filter(freq>10) %>%
  print( n=50 )

# this checks if a guide is affected in more than 2 mice
remove <- to_remove %>%
  select(mouse,Id)

id_counts <- table(remove$Id)

# Check if any Id appears more than 2 times
more_than_two <- id_counts[id_counts > 2]
more_than_two

# Print the Ids that appear more than twice
if (length(more_than_two) > 0) {
  cat("The following Ids appear more than 2 times:\n")
  print(more_than_two)
} else {
  cat("No Id appears more than 2 times.\n")}
  
#remove them from the list of guides to remove
ids_to_remove <- names(more_than_two)

remove <- remove %>%
  filter(!(Id %in% ids_to_remove))


#this is now the list of things to remove.
#remove them from the perCell Data
# this removes affected guides in every tissue in the affected mouse
perCellData <- perCellData %>%
  anti_join(remove, by = c("mouse", "Id"))


#--- make this into a pct_matrix

# Step 1: Count the number of cells per tissue, mouse, and Id
cell_counts <- perCellData %>%
  group_by(tissue, mouse, Id) %>%
  summarise(count = n(), .groups = "drop")

# Step 2: Pivot to wide format to create a matrix with Ids as columns
count_matrix <- cell_counts %>%
  pivot_wider(names_from = Id, values_from = count, values_fill = 0)

# Step 3: Calculate row-wise totals
count_matrix <- count_matrix %>%
  rowwise() %>%
  mutate(total = sum(c_across(-c(tissue, mouse)))) %>%
  ungroup()

count_matrix <- count_matrix %>% filter(total != 0)

# Step 4: Convert counts to percentages
percentage_matrix <- count_matrix %>%
  mutate(across(-c(tissue, mouse, total), ~ (.x / total) * 100))

# Step 5: remove samples below 20 cells total, then remove total column
# Subset the rows with fewer than 20 cells
removed_rows <- count_matrix[count_matrix$total < 20, ]

# Print each removed sample with mouse and tissue info and cell count
cat("Removed samples:\n")
apply(removed_rows, 1, function(row) {
  cat(row["mouse"], row["tissue"], "—", row["total"], "cells\n")
})

# Then apply the filter and drop the 'total' column
percentage_matrix <- percentage_matrix[count_matrix$total >= 20, ] %>%
  select(-total)

#check each row sums to 100
row_sums <- rowSums(percentage_matrix[ , -(1:2)], na.rm = TRUE)

# View the result
print(row_sums)

#change format to what the code expects
percentage_matrix <- percentage_matrix %>% 
  pivot_longer(cols = -c(tissue, mouse), names_to = "target_gene", values_to = "freq" )


#this code has an entry for every gene in every mouse, disregarding the sets, that many of the
# zeros just mean that gene was not in that mouse as it wasn't in that set.

#so
#use the mouse names and guide names per set as previously saved to remove these
# this is unecessary if you have just one set


# Filter out rows for set1
filtered_matrix <- percentage_matrix %>%
  filter(!(freq == 0 & mouse %in% set1_mice & !(target_gene %in% set1_guides)))

# Filter out rows for set2
filtered_matrix <- filtered_matrix %>%
  filter(!(freq == 0 & mouse %in% set2_mice & !(target_gene %in% set2_guides)))

# Filter out rows for set3
filtered_matrix <- filtered_matrix %>%
  filter(!(freq == 0 & mouse %in% set3_mice & !(target_gene %in% set3_guides)))

#for subsequent code to work
pct_matrix <- filtered_matrix

#check the rows are approx equal to (set1's mice x tissues x guides) + (set2's mice x tissues x guides) etc
dim(pct_matrix)

#checking it worked, can choose a guide, the output should be the correct number of mice in that set

checks <- pct_matrix %>%
   filter(target_gene == "Ccr4") %>%
   distinct(mouse) %>%
   nrow()

checks
}


#### end ####



# rejoin code here if not removing clones
# use this
# calculate minimum detection frequency per sample from perCellData

if ( ASSIGNING_MIN_VALUES <- TRUE ) {
  tissue_counts <- perCellData %>%
    group_by( tissue, mouse, Id ) %>%
    filter( Id != "Other") %>%
    summarize( count = n()) %>%
    mutate(pct = round(count/sum(count)*100, 2)) 
  
  
  tissue_mins <- tissue_counts %>%
    group_by( tissue, mouse ) %>%
    summarize( total = sum(count) ) %>%
    mutate( min_det = 100/total )

# calculate half-min frequency to replace zeros in percent matrix

tissue_mins$half_min <- tissue_mins$min_det/2

#optional - for checking differences in minimum detection frequencies between a reference tissue and tissue of interest
#this is pancreas and spleen
#sort of for checking what may cause skew and how best to compensate or correct in a way that is scientifically sound


tissue_mins_sep <- tissue_mins

# first check if there is a difference between the mins of different tissues which you think is sufficent to cause issues
# Filter for Pancreas and rename column
pan_half_mins <- tissue_mins_sep %>%
  filter(tissue == "Pancreas") %>%
  select(mouse, half_min) %>%
  rename(pan_half_min = half_min)

# Filter for Spleen and rename column
spl_half_mins <- tissue_mins_sep %>%
  filter(tissue == "Spleen") %>%
  select(mouse, half_min) %>%
  rename(spl_half_min = half_min)

# Print average values
mean_pan <- mean(pan_half_mins$pan_half_min, na.rm = TRUE)
mean_spl <- mean(spl_half_mins$spl_half_min, na.rm = TRUE)

print(paste("Average pan_half_min:", mean_pan))
print(paste("Average spl_half_min:", mean_spl))


# Combine dataframes by mouse and compare values
combined_mins <- pan_half_mins %>%
  inner_join(spl_half_mins, by = "mouse") %>%
  mutate(pan_greater_than_spl = pan_half_min > spl_half_min)
}


#### option 1 - make the half mins for non-spleen tissues = twice the mouses spleen half-mins ####
#  has been recommended to just use as standard for all non-spleen tissues

all_mins <- tissue_mins_sep

spl_half_mins <- spl_half_mins %>%
  ungroup() %>%
  select(-tissue)


all_mins <- all_mins %>%
  left_join(spl_half_mins %>% select(mouse, spl_half_min), by = "mouse")


all_mins <- all_mins %>%
  mutate(half_min = if_else(tissue != "Spleen", 2 * spl_half_min, half_min))

tissue_mins <- all_mins
  

#### option 2 - make the pancreas half_mins  = twice that mouses spleen half_mins ####
# taking the tissue_mins_sep made in the previous section and replacing the pan values with x2 spl
combined_mins$pan_half_min <- 2*(combined_mins$spl_half_min)

# First, create a lookup table from combined_mins
pan_lookup <- combined_mins %>%
  select(mouse, pan_half_min)

# Update half_min in tissue_mins_sep for Pancreas rows
tissue_mins_sep <- tissue_mins_sep %>%
  left_join(pan_lookup, by = "mouse") %>%
  mutate(half_min = if_else(tissue == "Pancreas", pan_half_min, half_min)) %>%
  select(-pan_half_min)  # Remove the temporary column

#make the change fit back into original code
tissue_mins <- tissue_mins_sep



#### or option 3 -leaving the minimum detections as they are ####

pan_lookup <- combined_mins %>%
  select(mouse, pan_half_min)

# Update half_min in tissue_mins_sep for Pancreas rows
tissue_mins_sep <- tissue_mins_sep %>%
  left_join(pan_lookup, by = "mouse") %>%
  mutate(half_min = if_else(tissue == "Pancreas", pan_half_min, half_min)) %>%
  select(-pan_half_min)  # Remove the temporary column

#i think this line could have been missing from the previous version
tissue_mins <- tissue_mins_sep

#back to normal code (all options rejoin here)

tissue_mins <- tissue_mins %>% unite(col = sample, tissue, mouse)

pct_matrix_corrected <- pct_matrix
pct_matrix_corrected <- pct_matrix_corrected %>% unite(col = sample, tissue, mouse)

for (s in unique(pct_matrix_corrected$sample)) {
  min_temp <- filter(tissue_mins, sample == s)
  pct_temp <- filter(pct_matrix_corrected, sample == s)
  
  if (nrow(min_temp) == 1) {
    pct_temp$freq[pct_temp$freq == 0] <- min_temp$half_min
    pct_matrix_corrected$freq[pct_matrix_corrected$sample %in% unique(pct_temp$sample)] <- pct_temp$freq
  }
}



# Cell number cutoff has been done previously if you made the pct_matrix from the perCell Data
pct_matrix_corrected <- pct_matrix_corrected %>% select(-total)

pct_matrix_corr <- left_join(pct_matrix_corrected, tissue_mins, by = "sample")

#if it hasn't, filter out the <20 samples
#sometimes, if it cant find the total column, the relevant total column is called total.y

filtered_out <- pct_matrix_corr %>% filter(total <= 20)

cat("Samples filtered out due to total <= 20:\n")
print(unique(filtered_out$sample)) 

pct_matrix_corr <- pct_matrix_corr %>% filter(total > 20)

pct_matrix_corr <- pct_matrix_corr %>% 
  separate(sample, remove = TRUE, sep = "_", into = c("tissue", "mouse")) 
  #select(-total, -min_det, -half_min)

anyNA(pct_matrix_corr)




#back to normal code

# before you run any t test you should need to check your input data is normally distributed.
# it probably isn't, so you then check if the log2 transformation (used here) makes it sufficiently normal
# Vaclav has written code that both of these things under "Shapiro normality histograms", and I've modified it a bit 
# and put it here but it takes at least a few minutes to run so you don't need to rerun it after you've initially checked it

####  -------   Shapiro Normality Hitograms  -----------   ####

#needs the pct matrix corr to by in tissue abundance or activation subset format depending on what you want
# modify as needed
if (RUN_SHAPIRO_TESTING <-  TRUE) {

# requires packages patchwork and readr
library(patchwork)
library(readr)

comparison.results <- data.frame()
l_hist <- list()
l_hist_log <- list()

# code for tissue abundance. its the same but the naming changes
graphs.dir <- paste0("./", flowcode.seed, "_graphs/")

if( !file.exists(graphs.dir) ){
  dir.create(graphs.dir)
}

comparison.dir <- "tissue_vs_Spleen/"
if( !file.exists( file.path(graphs.dir, comparison.dir) ) ){
  dir.create(file.path(graphs.dir, comparison.dir) )
}
reference_tissue <- "Spleen"



dir.create(file.path(graphs.dir, comparison.dir, "per_gene_histograms"), recursive = TRUE, showWarnings = FALSE)


#sometimes this runs into issues if due to clone removal or other reasons, for a gene and tissue, there are less than 3 datapoint remaining
#you can figure out which genes and tissues are the issue and remove them from the dataset that plots this graph, if its a few
#make sure you do not remove them from the dataset used for subsequent analysis.

#### does not need to be in masterscript
# but would make the code work better
#added 20250711

# option1 - for when no guide has <3 mice remaining after filtering
pct_matrix_corr_smaller <- pct_matrix_corr

### option2 for when one or more guides have<3 mice remaining after filtering
#as that causes the Shapiro code to fail

#I made it just an example tissue and Spleen as a control, it can be anything and spleen.
# or you could do all your tissues
pct_matrix_corr_smaller <- pct_matrix_corr %>%
  filter(tissue %in% c("Pancreas", "Spleen")) 


# Split by tissue
pancreas_df <- pct_matrix_corr_smaller %>% filter(tissue == "Pancreas")
spleen_df   <- pct_matrix_corr_smaller %>% filter(tissue == "Spleen")

# Count target_gene occurrences
pancreas_counts <- pancreas_df %>%
  count(target_gene, name = "pancreas_count")

spleen_counts <- spleen_df %>%
  count(target_gene, name = "spleen_count")

# Join counts and replace missing with 0
gene_counts <- full_join(pancreas_counts, spleen_counts, by = "target_gene") %>%
  mutate(across(everything(), ~replace_na(.x, 0)))

# Find genes with < 3 in either tissue
low_count_genes <- gene_counts %>%
  filter(pancreas_count < 3 | spleen_count < 3)

# Print genes and counts
print(low_count_genes)

# Remove those genes from the dataframe to check stats
pct_matrix_corr_smaller <- pct_matrix_corr_smaller %>%
  filter(!target_gene %in% low_count_genes$target_gene)


#### end

#for all options this now runs the shapiro testing
for (t in unique(pct_matrix_corr_smaller$tissue)) {
  
  if(t != "Spleen"){
    tissue.dir <- paste0(graphs.dir, comparison.dir, t, "/")
    
    if( !file.exists(tissue.dir)){
      dir.create(tissue.dir)
    }
    
    toi = pct_matrix_corr_smaller %>% filter(tissue == t)
    
    spleen_toi = pct_matrix_corr_smaller %>% filter(tissue == reference_tissue)
    
    resdf = data.frame(target_gene = sort(unique(toi$target_gene)), log2FC = 0, pvalue = 1, diff_shap_pval = 1, logdiff_shap_pval = 1 )
    
    for (i in 1:nrow(resdf)) {
      g = resdf$target_gene[i]
      
      xvdf = (toi %>% filter(target_gene == g) %>% arrange(mouse))
      yvdf = spleen_toi %>% filter(target_gene == g)%>% filter(mouse %in% xvdf$mouse) %>% arrange(mouse)
      xv = (xvdf %>% filter(mouse %in% yvdf$mouse))$freq
      yv = yvdf$freq
      
      if (length(c(xv,yv)) >2){
        # VACLAV ADDED
        if ( VACLAV_ADDED <- TRUE ) {
          
          dd <- data.frame( 
            diff = (xv + 0.001) - (yv + 0.001),
            logdiff = log2((xv + 0.001)) - log2((yv + 0.001))) 
          
          resdf$diff_shap_pval[i] <- shapiro.test( dd$diff )$p.value
          resdf$logdiff_shap_pval[i] <- shapiro.test( dd$logdiff )$p.value
          
          l_hist[[ i ]] <- dd %>% 
            ggplot(aes(diff)) + 
            geom_histogram(aes(y = ..density..), bins = 10, fill = "#D6B996") +
            geom_density(color = "red") + 
            labs( x = "Cell Prop in gr1 - Cell Prop in gr2", 
                  y = "Density / Normalized count",
                  title = paste0(g, ", Normality Shapiro's P=", 
                                 round(resdf$diff_shap_pval[i], 4))) +
            theme_bw(); l_hist[[ i ]]
          
          l_hist_log[[ i ]] <- dd %>% 
            ggplot(aes(logdiff)) + 
            geom_histogram(aes(y = ..density..), bins = 10, fill = "#DAB7E2") +
            geom_density(color = "red") + 
            labs( x = "log2(Cell Prop in gr1) - log2(Cell Prop in gr2)", 
                  y = "Density / Normalized count",
                  title = paste0(g, ", Normality Shapiro's P=", 
                                 round(resdf$logdiff_shap_pval[i], 4))) +
            theme_bw(); l_hist_log[[ i ]]
          
          gg_hist <- l_hist[[ i ]] + l_hist_log[[ i ]]
          
          
          
          ggsave(paste0(graphs.dir, comparison.dir, "per_gene_histograms/normality_", t, "_vs_", reference_tissue, "_", g, ".pdf"), 
                 plot = gg_hist, height = 10, width = 20, units = "cm")
          
        }
        
        resdf$log2FC[i] = mean(log2((xv + 0.001) / (yv + 0.001)))
        
        res <- t.test(x = log2(xv + 0.001), y = log2(yv + 0.001), 
                      paired = TRUE, alternative = "two.sided")
        resdf$pvalue[i] = res$p.value
        #removed excess code making here
      }
    }
    
    
    
    resdf$tissue <- rep(t)
    
    # Apply multiple hypothesis testing correction
    resdf$pvalue <- p.adjust(resdf$pvalue, method = "BH")
    
    maxmlog10pvalue <- max(-log10(resdf$pvalue))
    maxlog2FC <- max(abs(resdf$log2FC))
    
    resdf$dist0 <- sqrt(resdf$log2FC^2 + (-log10(resdf$pvalue) * maxlog2FC / maxmlog10pvalue)^2)
    resdf <- resdf %>% arrange(-dist0)
    
    comparison.results <- rbind(comparison.results, resdf)
    write.csv(resdf, paste0(tissue.dir, t, "_vs_", reference_tissue, "_pvals.csv"))
    
    # VACLAV ADDED:
    write.csv(resdf, paste0(graphs.dir, comparison.dir, t, "_vs_", reference_tissue, "_pvals_shap.csv"))
    write.csv(resdf %>% select( -diff_shap_pval, -logdiff_shap_pval), 
              paste0(graphs.dir, comparison.dir, t, "_vs_", reference_tissue, "_pvals.csv"))
    
  }
  
}

#set ti- to match the tissue in the data you import for this, that you are interested in
ti<-"Pancreas"

#the second line is for the title, so copy paste the file origin, or at least what is true
file_path <- paste0(graphs.dir, comparison.dir, ti, "_vs_Spleen_pvals_shap.csv")
dd_shap <- read.csv(file_path)
file_input <- paste0(ti, " vs Spleen")

# if the above line can't find the file, just import it directly by name
#dd_shap <- read.csv("")


gg_norm_pvals_orig <- dd_shap %>% 
  ggplot( aes( x = diff_shap_pval ) ) +
  geom_histogram(fill = "#966427") + 
  labs(x = "P values of Shapiro test for all genes", 
       y = "Count",
       title = paste0("Differences in Proportions")) +
  annotate(   "text", 
              x = Inf, y = Inf, 
              label = file_input, 
              hjust = 1.1, vjust = 2, 
              size = 0.5, fontface = "italic",
              colour = "lightgrey")+
  theme_bw(); gg_norm_pvals_orig

gg_norm_pvals_log <- dd_shap %>% 
  ggplot( aes( x = logdiff_shap_pval ) ) +
  geom_histogram(fill = "#9D4FAB") + 
  labs(x = "P values of Shapiro test for all genes", 
       y = "Count",
       title = paste0("Differences in log2(Proportions)")) +
  theme_bw(); gg_norm_pvals_log

gg_norm_pvals <- gg_norm_pvals_orig + gg_norm_pvals_log; gg_norm_pvals
ggsave(paste0(graphs.dir, comparison.dir, "norm_pval_distribution_", ti, "_vs_", reference_tissue, ".pdf"), 
       plot = gg_norm_pvals, height = 10, width = 20, units = "cm")

}

#### end ####
 

#### START OF FROM_VACLAV ####
# - RANKS + NONPARAMETRIC STATISTICS

# Note, that for the non-parametric approach some of the previous steps
# are not needed (e.g., dealing with minimal values).
# But they are kept for keeping minimum changes.


# Check if there are no missing data
if ( CHECK_DATA <- TRUE ) {  
  pct_matrix_corr %>% 
    group_by( tissue, mouse ) %>% 
    summarise( n_genes_per_mouse_tissue = n() ) %>%
    pivot_wider(
      names_from = tissue,
      values_from = n_genes_per_mouse_tissue
    ) %>% as.data.frame()  
  # %>% View

}

# Do we deal with TCR expansion and incomparable distributions?
# ie tests if the "remove some of the most expanded clones" method works well enough
# or non-parametric testing would be better
if ( ASSESS_UNEQUAL_DISTRIBUTIONS__TCR_CLONES <- TRUE ) {
  # install.packages("ineq")
  library(ineq)
  
  # Let's assess the if the distribution of proportions of cells 
  # with given target_gene and compare that across organs.
  # Measure: Gini coefficient (0...equal, 1...absolutelly unequal)
  # Particularly high (e.g., relatively among organs)  Gini coeff.
  # values will strongly suggest presence of TCR expansion.
  # Significant difference between Gini coeff will requiere non-parametric
  # testing of ranks
  pct_matrix_corr %>% 
    group_by( mouse, tissue ) %>%
    summarise( 
      gini_score = ineq( freq, type = "Gini" ) ) %>% 
    ggplot( aes( x = tissue, y = gini_score ) ) +
    geom_boxplot( aes( colour = tissue ) ) +
    labs( title = "Nothing removed" )
  ggsave( file = "gini_removed_0_per_mouse.png", width = 5, height = 5 )
  
  # Q: Will removing of the most expanded TCR per mouse help?
  # A: No.  
  pct_matrix_corr %>% 
    group_by( mouse, tissue ) %>%  
    filter( min_rank( desc( freq ) ) != 1 ) %>%
    summarise( 
      gini_score = ineq( freq, type = "Gini" ) ) %>% 
    ggplot( aes( x = tissue, y = gini_score ) ) +
    geom_boxplot( aes( colour = tissue ) ) +
    labs( title = "1 most frequent freqency PER MOUSE removed" )
  ggsave( file = "gini_removed_1_per_mouse.png", width = 5, height = 5 )
  
  # Q: Will removing of two most expanded TCR per mouse help?
  pct_matrix_corr %>% 
    group_by( mouse, tissue ) %>% 
    filter( min_rank( desc( freq ) ) > 2 ) %>% 
    summarise( 
      gini_score = ineq( freq, type = "Gini" ) ) %>% 
    ggplot( aes( x = tissue, y = gini_score ) ) +
    geom_boxplot( aes( colour = tissue ) ) +
    labs( title = "2 most frequent freqency PER MOUSE removed" )
  ggsave( file = "gini_removed_2_per_mouse.png", width = 5, height = 5 )
  
  pct_matrix_corr %>% 
    group_by( mouse, tissue ) %>% 
    filter( min_rank( desc( freq ) ) > 3 ) %>% 
    summarise( 
      gini_score = ineq( freq, type = "Gini" ) ) %>% 
    ggplot( aes( x = tissue, y = gini_score ) ) +
    geom_boxplot( aes( colour = tissue ) ) +
    labs( title = "3 most frequent freqency PER MOUSE removed" )
  ggsave( file = "gini_removed_3_per_mouse.png", width = 5, height = 5 )
  
  pct_matrix_corr %>% 
    group_by( mouse, tissue ) %>% 
    filter( min_rank( desc( freq ) ) > 4 ) %>% 
    summarise( 
      gini_score = ineq( freq, type = "Gini" ) ) %>% 
    ggplot( aes( x = tissue, y = gini_score ) ) +
    geom_boxplot( aes( colour = tissue ) ) +
    labs( title = "4 most frequent freqency PER MOUSE removed" )
  ggsave( file = "gini_removed_4_per_mouse.png", width = 5, height = 5 )
  # The TCR expansion is present even in much smaller scales. No top clone  
  # removal can correct it.
  
  # The consequence for proportions of genes.
  pct_matrix_corr %>% 
    group_by( target_gene, tissue ) %>% 
    filter( min_rank( desc( freq ) ) > 1 ) %>% 
    summarise( 
      gini_score = ineq( freq, type = "Gini" ) ) %>% 
    ggplot( aes( x = tissue, y = gini_score ) ) +
    geom_boxplot( aes( colour = tissue ) ) +
    labs( title = "1 most frequent freqency PER GENE removed" )
  ggsave( file = "gini_removed_1_per_gene.png", width = 5, height = 5 )
  
  pct_matrix_corr %>% 
    group_by( target_gene, tissue ) %>% 
    filter( min_rank( desc( freq ) ) > 2 ) %>% 
    summarise( 
      gini_score = ineq( freq, type = "Gini" ) ) %>% 
    ggplot( aes( x = tissue, y = gini_score ) ) +
    geom_boxplot( aes( colour = tissue ) ) + 
    labs( title = "2 most frequent freqency PER GENE removed" )
  ggsave( file = "gini_removed_2_per_gene.png", width = 5, height = 5 )
}

# Based on the analysis above, determine if you want to ranks and non-parametric test for Pancreas
# set this to TRUE or FALSE
# if you set it to FALSE then statistical tests will run on comparing percentage frequencies of the guides,
# if you set it to TRUE then the statistical tests will compare the guides non-parametrically, by ranks
USE_RANKS <- FALSE

# Translate ranks back to spleen values, i.e. assuming there is 
# the same distribution/TCR clonality in tissue and spleen - assuming this 
# (and thus removing the signal of excessive TCR expansion), 
# we want to measure the effect size in terms of logFC as in other comparisons
# it is an assumption, but should be reasonably accurate
# if you set it to FALSE, the code will output graphs with the x axis of change in ranks, not logfc
RANK_BASED_LOGFC <- FALSE


if ( USE_RANKS ) {
  pct_matrix_rank <- pct_matrix_corr %>% 
    group_by( tissue, mouse ) %>% 
    mutate( rank_raw = rank( freq, ties.method = "average" ) ) %>% 
    mutate( mean_rank_raw = mean( rank_raw ) ) %>% 
    ungroup() %>%  
    
    # Normalize to deal with unequal number of mice per gene
    mutate( max_mean_rank_raw = max( mean_rank_raw ) ) %>% 
    mutate( rank = rank_raw * ( max_mean_rank_raw / mean_rank_raw ) )
}

#this is only an option if use_ranks is set to TRUE
if ( CHECK_DATA_REVERT_TO_SPLEEN_DISTRIB <- TRUE ) {
  pct_matrix_rank %>% 
    group_by( tissue, target_gene ) %>% 
    mutate( med_freq = median( freq ),
            med_rank = median( rank ) ) %>%
    ungroup() %>% 
    ggplot( aes( x = med_rank, y = med_freq ) ) +
    geom_point() +
    geom_smooth( method = loess ) +
    facet_wrap( vars( tissue ) )
  
  
  d_spleen_ranks <- pct_matrix_rank %>% 
    filter( tissue == "Spleen" ) %>%
    group_by( tissue, target_gene ) %>% 
    mutate( med_freq = median( freq ),
            med_rank = median( rank ) ) %>%
    ungroup() 
  
  d_pancreas_ranks <- pct_matrix_rank %>% 
    filter( tissue == "Pancreas" ) %>%
    group_by( tissue, target_gene ) %>% 
    mutate( med_freq = median( freq ),
            med_rank = median( rank ) ) %>%
    ungroup() 
  
  d_spleen_ranks %>% 
    ggplot( aes( x = med_rank, y = med_freq ) ) +
    geom_point() +
    geom_smooth( method = loess, span = 0.4 ) +
    facet_wrap( vars( tissue ) )
  
  
  fit_spleen <- smooth.spline( 
    x = d_spleen_ranks$med_rank, y = d_spleen_ranks$med_freq, spar = 0.8 )
  rank_to_value <- function( xval ) predict( fit_spleen, x = xval )$y
  
  plot( d_pancreas_ranks$rank, rank_to_value( d_pancreas_ranks$rank ) )
}



## Simple tissue vs Spleen freq comparison----------------
graphs.dir <- paste0("./", flowcode.seed, "_graphs/")

if( !file.exists(graphs.dir) ){
  dir.create(graphs.dir)
}

comparison.dir <- "tissue_vs_Spleen/"
if( !file.exists( file.path(graphs.dir, comparison.dir) ) ){
  dir.create(file.path(graphs.dir, comparison.dir) )
}

reference_tissue <- "Spleen"

# this performs the chosen statistical test and outputs boxplots and volcanos,
# it depends on if you have set USE_RANKS to TRUE or FALSE
# within the volcano code you can adjust p value significance level and log2fc thresholds

# VACLAV: The code for non-parametric approach is copied from below and modified to
# allow for USE_RANKS by using "flag"
# STILL FROM_VACLAV
if ( PLOT_BOXPLOTS_AND_VOLCANO_FOR_RANKS <- TRUE ) {
  
  
  comparison.results <- data.frame()
  
  if ( USE_RANKS ) { pct_matrix_corr <- pct_matrix_rank } 
  for (t in unique(pct_matrix_corr$tissue)) {
    
    if(t != "Spleen"){
      tissue.dir <- paste0(graphs.dir, comparison.dir, t, "/")
      
      if( !file.exists(tissue.dir)){
        dir.create(tissue.dir)
      }
      
      toi = pct_matrix_corr %>% filter(tissue == t)
      
      spleen_toi = pct_matrix_corr %>% filter(tissue == reference_tissue)
      
      resdf = data.frame(target_gene = sort(unique(toi$target_gene)), 
                         log2FC = 0, pvalue_t = 1, pvalue = 1,
                         median_rank_diff = 0, pvalue_w = 1,
                         rank_based_logfc = 0 )
      
      for (i in 1:nrow(resdf)) {
        g = resdf$target_gene[i]
        xvdf = (toi %>% filter(target_gene == g) %>% arrange(mouse))
        yvdf = spleen_toi %>% filter(target_gene == g)%>% filter(mouse %in% xvdf$mouse) %>% arrange(mouse)
        
        xv = (xvdf %>% filter(mouse %in% yvdf$mouse))$freq
        yv = yvdf$freq
        
        if (length(c(xv,yv)) >2){
          
          # Parametric testing
          resdf$log2FC[i] = mean(log2((xv + 0.001) / (yv + 0.001)))
          res_t <- t.test(x = log2(xv + 0.001), y = log2(yv + 0.001),
                          paired = TRUE, alternative = "two.sided")
          resdf$pvalue_t[i] = res_t$p.value
          
          ggplot(data.frame(toi = xv, spleen = yv, id = 1:length(xv)) %>% 
                   pivot_longer(cols = -id, names_to = "group", values_to = "freq")) +
            scale_x_discrete(name = "", labels = c(reference_tissue, t)) +
            scale_y_log10() +
            scale_fill_manual(values = c("darkgrey", "firebrick1")) +
            ggtitle(paste0(g, ": ", t, " vs ", reference_tissue, 
                           "\n log2FC = ", round(resdf$log2FC[i], 3), 
                           "; t-test p = ", round(res_t$p.value, 3))) +
            geom_boxplot(aes(x = group, y = freq, fill = group), alpha = 0.2) +
            geom_point(aes(x = group, y = freq)) +
            geom_line(aes(x = group, y = freq, group = id)) +
            theme_bw() +
            theme(legend.position = "none")
          
          ggsave(paste0(tissue.dir, t, "_vs_", reference_tissue, "_", g, "_prop.pdf"), 
                 height = 10, width = 10, units = "cm")
          
          # Nonparametric testing
          if ( USE_RANKS ) {
            xv_r = (xvdf %>% filter(mouse %in% yvdf$mouse))$rank
            yv_r = yvdf$rank
            resdf$median_rank_diff[i] = median(xv_r - yv_r)
            res_w <- wilcox.test(x = xv_r + 0.001, y = yv_r + 0.001, 
                                 paired = TRUE, alternative = "two.sided")
            resdf$pvalue_w[i] = res_w$p.value
            
            ggplot(data.frame(toi = xv_r, spleen = yv_r, id = 1:length(xv_r)) %>% 
                     pivot_longer(cols = -id, names_to = "group", values_to = "rank")) +
              scale_x_discrete(name = "", labels = c(reference_tissue, t)) +
              # scale_y_log10() +
              scale_fill_manual(values = c("darkgrey", "firebrick1")) +
              ggtitle(paste0(g, ": ", t, " vs ", reference_tissue, 
                             "\n Rank diff = ", round(resdf$median_rank_diff[i], 3), 
                             "; Wilcox p = ", round(res_w$p.value, 3))) +
              geom_boxplot(aes(x = group, y = rank, fill = group), alpha = 0.2) +
              geom_point(aes(x = group, y = rank)) +
              geom_line(aes(x = group, y = rank, group = id)) +
              theme_bw() +
              theme(legend.position = "none")
            
            ggsave(paste0(tissue.dir, t, "_vs_", reference_tissue, "_", g, "_rank.pdf"), 
                   height = 10, width = 10, units = "cm")
            
          }
          if ( USE_RANKS & RANK_BASED_LOGFC ) {
            resdf$rank_based_logfc[i] = log2( 
              rank_to_value( median( xv_r ) ) / rank_to_value( median( yv_r ) ) )
            
          }
        }
      }
      
      
      
      resdf$tissue <- rep(t)
      
      # Apply multiple hypothesis testing correction
      if ( USE_RANKS ) {
        resdf$pvalue <- p.adjust(resdf$pvalue_w, method = "BH")  # FDR, Benjamini-Hochberg
      } else {
        resdf$pvalue <- p.adjust(resdf$pvalue_t, method = "BH")  # FDR, Benjamini-Hochberg
      }
      maxmlog10pvalue <- max(-log10(resdf$pvalue))
      maxlog2FC <- max(abs(resdf$log2FC))
      
      resdf$dist0 <- sqrt(resdf$log2FC^2 + (-log10(resdf$pvalue) * maxlog2FC / maxmlog10pvalue)^2)
      resdf <- resdf %>% arrange(-dist0)
      
      comparison.results <- rbind(comparison.results, resdf)
      write.csv(resdf, paste0(tissue.dir, t, "_vs_", reference_tissue, 
                              ifelse( USE_RANKS, "_ranks", "" ), ".csv"))
      
    }
    
  }
  
  
  
  if ( VG_VOLCANO_ON_RANKS_NORM <- TRUE ) {
    
    tissues <- c( "Pancreas", "LN", "Blood" )
    reference_subset <- "Spleen"
    
    for (s in tissues) {
      subset.results.vs.res <- read.csv(
        file.path(
          paste0( flowcode.seed, "_graphs" ), "tissue_vs_Spleen", s,
          paste0(s, "_vs_", reference_tissue, 
                 ifelse( USE_RANKS, "_ranks", "" ), ".csv") ) )
      
      
      # option 1 - normal plotting of everything with the log2fc threshold
      subset.results <- subset.results.vs.res
      if ( USE_RANKS ) {
        subset.results$x_value <- subset.results$median_rank_diff
        x_label <- "Rank difference (across 55 genes)"    
        # if less than 55 genes the ranks are already normalized. VG prefers this 
        # instead of using some relative values as this directly looks like ranks.
      } else {
        subset.results$x_value <- subset.results$log2FC
        x_label <- "log2FC"
      }
      
      y_max <- max(-log10(subset.results$pvalue))
      x_max <- max( abs( subset.results$x_value ) )
      
      options(ggrepel.max.overlaps = Inf)
      
      plot.figure <- paste0(s, "_vs_", reference_subset, "_CD8_Volcano p=0.01_para.png" )
      
      resdf <- subset.results %>% filter( tissue == s )
      resdf <- resdf %>% arrange(-dist0)
      
      # Store the ggplot object in a variable
      p <- ggplot(resdf) +
        ggtitle(paste(s, "vs", reference_subset, "CD8_Volcano p=0.01")) +
        labs( x = x_label ) +
        scale_x_continuous(limits = c(-x_max - 0.5, x_max + 0.5)) +
        scale_y_continuous(limits = c(0, y_max + 1.5)) +
        geom_point(data = resdf, aes(x = x_value, y = -log10(pvalue))) +
        geom_hline(yintercept = -log10(0.01), color = "firebrick1") +
        geom_vline(xintercept = c(-1, 1), color = "red", linetype = "dashed") +
        #normally I set this to the pval so significant genes are labelled, choose your thresholds
        geom_text_repel(data = resdf %>% filter(pvalue < 0.01, x_value < -1 | x_value > 1), 
                        aes(x = x_value, y = -log10(pvalue), label = target_gene),
                        nudge_x = -0.5) +
        theme_bw(); p
      
      # Pass the ggplot object instead of filename
      ggsave(filename = plot.figure, plot = p, 
             height = 10, width = 10, units = "cm", dpi = 300) 
      
      # translate ranks back to spleen values, i.e. assuming there is 
      # the same distribution/TCR clonality in tissue and spleen
      if ( RANK_BASED_LOGFC ) {
        x_max <- max( abs( subset.results$rank_based_logfc ) )
        plot.figure_2 <- paste0(s, "_vs_", reference_subset, 
                                "_CD8_Volcano p=0.05_rank_logfc.png" )
        p_rank_logfc <- ggplot(resdf) +
          ggtitle(paste(s, "vs", reference_subset, "CD8_Volcano p=0.05")) +
          labs( x = "Rank-based logFC" ) +
          scale_x_continuous(limits = c(-x_max - 0.5, x_max + 0.5)) +
          scale_y_continuous(limits = c(0, y_max + 1.5)) +
          geom_point(data = resdf, aes(x = rank_based_logfc, y = -log10(pvalue))) +
          geom_hline(yintercept = -log10(0.05), color = "firebrick1") +
          geom_vline(xintercept = c(-1, 1), color = "red", linetype = "dashed") +
          #normally I set this to the pval so significant genes are labelled, choose your thresholds
          #match them to the red lines
          geom_text_repel(data = resdf %>% filter(pvalue < 0.05, rank_based_logfc < -1 | rank_based_logfc > 1),
                          aes(x = rank_based_logfc, y = -log10(pvalue), 
                              label = target_gene),
                          nudge_x = -0.5) +
          theme_bw(); p_rank_logfc
        
        # Pass the ggplot object instead of filename
        ggsave(filename = plot.figure_2, plot = p_rank_logfc, 
               height = 10, width = 10, units = "cm", dpi = 300) 
      }
    }
    
  }

}

# END OF FROM_VACLAV










#plots heatmaps corresponding to the volcanos just made previously
# they correspond if you input in the same p value and log2fc thresholds
# should work for if you are using or not using ranks

if ( HEATMAPS_CLONE_AUTOREMOVED <- F ) {
  #  Remove the column 'distr0' from comparison.results
  comparison.results <- comparison.results %>% select(-dist0)
  
  #export out the comparison results as a csv for use in other scripts
  #change naming to be descriptive
  write.csv(comparison.results, "CD8_comparison_results.csv", row.names = FALSE)
  
  #taken from my script "log_abundance_volcano_and_heatmap"
  #replacing the standard heatmap and volcano
  
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(ggplot2)
  
  
  heatmap.input <- read.csv("CD8_comparison_results.csv")
  
  #to format the result of the statistical test for a heatmap
  # you can alter the p value threshold as you like
  
  #option 1 - when you have used parametric testing 
  if (USE_RANKS == FALSE) {
    # heatmap specific code, to only show colours for significant p<0.01 values
    heatmap.input[heatmap.input$pvalue>0.01,]$log2FC <- 0
    #and the log2fc threshold - to match any log2fc threshold your volcano has, optional
    heatmap.input[heatmap.input$log2FC > -1 & heatmap.input$log2FC < 1, ]$log2FC <- 0
    
    heatmap.input <- heatmap.input %>% select(-pvalue)
  }
  
  #or option 2 - when you have used non parametric testing ( ranks)
  if (USE_RANKS == TRUE){
    heatmap.input <-  heatmap.input %>% select(target_gene,pvalue, rank_based_logfc,tissue)
    # heatmap specific code, to only show colours for significant p<0.05 or whatever you want, values
    heatmap.input[heatmap.input$pvalue>0.05,]$rank_based_logfc <- 0
    #and the log2fc threshold - to match any log2fc threshold your volcano has, optional
    heatmap.input[heatmap.input$rank_based_logfc > -1 & heatmap.input$rank_based_logfc < 1, ]$rank_based_logfc <- 0
    heatmap.input <- heatmap.input %>% select(-pvalue)
    #for the subsequent code to work, this is renamed, though its not an accurate name
    heatmap.input <- heatmap.input %>%rename(log2FC = rank_based_logfc)
  }
  
  #both options continue from here
  #  Aggregating duplicate values using mean
  
  heatmap.input <- heatmap.input %>%
    group_by(target_gene, tissue) %>%
    summarise(log2FC = mean(log2FC)) %>%
    ungroup()
  
  # Assign group labels to tissues
  
  group.tissues <- c( "Pancreas")
  
  group.blood <- c("Blood")
  
  group.lymphoid <- c("Spleen", "LN")
  
  # Define order of tissues based on groups
  tissue_order <- c(group.blood, group.lymphoid, group.tissues)
  
  
  # Calculate the average log2FC for each target_gene across all tissues
  average_log2FC_between_tissues <- aggregate(log2FC ~ target_gene, data = heatmap.input, FUN = mean)
  
  # Rank the target genes based on the average log2FC
  # This will naturally order from most negative to most positive
  rank_log2FC <- rank(average_log2FC_between_tissues$log2FC)
  
  # Create a data frame with target_gene and rank_log2FC
  ranked_genes <- data.frame(target_gene = average_log2FC_between_tissues$target_gene, 
                             rank_log2FC = rank_log2FC,
                             avg_log2FC = average_log2FC_between_tissues$log2FC)
  
  ### option 1 - genes shown in order of rank_log2FC ###
  # This will order from most negative to most positive log2FC
  ranked_genes <- ranked_genes[order(ranked_genes$avg_log2FC), ]
  
  
  ### option 2 - genes shown in alphabetical order ###
  # better for when you have lots of heatmaps and want to compare results for a guide between them
  ranked_genes <- ranked_genes %>%
    arrange(target_gene) %>%
    mutate(rank_log2FC = row_number())
  
  # all options rejoin here - reshape data into wide format: one row per gene, one column per tissue
  heatmap.input <- heatmap.input %>% pivot_wider(names_from = tissue, values_from = log2FC)
  heatmap.rownames <- heatmap.input$target_gene
  heatmap.input <- heatmap.input %>% select(-target_gene)
  heatmap.input <- as.matrix(heatmap.input)
  rownames(heatmap.input) <- heatmap.rownames
  
  
  
  palette.length <- 50
  heatmap.colors <- colorRampPalette(c("blue", "white", "red"))(palette.length)
  heatmap.breaks <- c(seq(min(heatmap.input), 0, length.out=ceiling(palette.length/2) + 1), 
                      seq(max(heatmap.input)/palette.length, max(heatmap.input), length.out=floor(palette.length/2)))
  
  #------ Option 1 - long heatmap
  # Create one long heatmap with clustered columns based on designated groups
  # all these options are just design choices
  
  gene.impact.heatmap <- pheatmap(
    t(heatmap.input[ranked_genes$target_gene, order(match(colnames(heatmap.input), tissue_order))]),  # Transpose the input matrix, order rows by the calculated ranks, and order columns based on tissue_order
    color = heatmap.colors, 
    breaks = heatmap.breaks,
    cluster_rows = FALSE,  # Do not perform hierarchical clustering on genes
    cluster_cols = FALSE,  # Do not perform hierarchical clustering on tissues
    border_color = c("grey", "grey")  # Add grey borders between rows and columns
  )
  
  ggsave(gene.impact.heatmap, filename = "CD8_P0.01_heatmap.png", height = 10, 
         width = 60, units = "cm")
  
  ggsave(gene.impact.heatmap, filename = "CD8_P0.01_heatmap.pdf", height = 10, 
         width = 60, units = "cm")
  
  #-------- Option2 - 2 rows of heatmaps -----------
  # OR split the heatmap into two rows for legibility if you have a lot of guides
  
  library(gridExtra)
  
  # Split the genes into two groups
  half_length <- ceiling(length(ranked_genes$target_gene) / 2)
  genes_group1 <- ranked_genes$target_gene[1:half_length]
  genes_group2 <- ranked_genes$target_gene[(half_length + 1):length(ranked_genes$target_gene)]
  
  #might need this line if you get breaks arent unique error
  heatmap.breaks <- unique(heatmap.breaks)
  
  # Create the first heatmap
  gene.impact.heatmap1 <- pheatmap(
    t(heatmap.input[genes_group1, order(match(colnames(heatmap.input), tissue_order))]),  # Transpose the input matrix, order rows by the calculated ranks, and order columns based on tissue_order
    color = heatmap.colors, 
    breaks = heatmap.breaks,
    cluster_rows = FALSE,  # Do not perform hierarchical clustering on genes
    cluster_cols = FALSE,  # Do not perform hierarchical clustering on tissues
    border_color = c("#898989", "#898989")  # Add grey borders between rows and columns
  )
  
  # Create the second heatmap
  gene.impact.heatmap2 <- pheatmap(
    t(heatmap.input[genes_group2, order(match(colnames(heatmap.input), tissue_order))]),  # Transpose the input matrix, order rows by the calculated ranks, and order columns based on tissue_order
    color = heatmap.colors, 
    breaks = heatmap.breaks,
    cluster_rows = FALSE,  # Do not perform hierarchical clustering on genes
    cluster_cols = FALSE,  # Do not perform hierarchical clustering on tissues
    border_color = c("#898989", "#898989"),  # Add black borders between rows and columns
    legend = FALSE  # Remove the legend
  )
  
  # Arrange the two heatmaps vertically
  combined_heatmap <- grid.arrange(gene.impact.heatmap1$gtable, gene.impact.heatmap2$gtable, ncol = 1)
  
  # Save the combined heatmap
  # REMEMBER to change the name as to not overwrite any you previously made
  ggsave(combined_heatmap, filename = "CD8_P0.01_2rows_heatmap.png", height = 10, 
         width = 30, units = "cm")
  
  #-------- Option3 - 3 rows of heatmaps with log2fc threshold -----------
  # the threshold is optional, my volcano plots for this data had the log2fc threshold so I added it here to match.
  
  library(gridExtra)
  
  # Calculate the length of each group
  third_length <- ceiling(length(ranked_genes$target_gene) / 3)
  genes_group1 <- ranked_genes$target_gene[1:third_length]
  genes_group2 <- ranked_genes$target_gene[(third_length + 1):(2 * third_length)]
  genes_group3 <- ranked_genes$target_gene[(2 * third_length + 1):length(ranked_genes$target_gene)]
  
  # Might need this line if you get breaks aren't unique error
  heatmap.breaks <- unique(heatmap.breaks)
  
  # Create the first heatmap
  gene.impact.heatmap1 <- pheatmap(
    t(heatmap.input[genes_group1, order(match(colnames(heatmap.input), tissue_order))]),  # Transpose the input matrix, order rows by the calculated ranks, and order columns based on tissue_order
    color = heatmap.colors, 
    breaks = heatmap.breaks,
    cluster_rows = FALSE,  # Do not perform hierarchical clustering on genes
    cluster_cols = FALSE,  # Do not perform hierarchical clustering on tissues
    border_color = c("#898989", "#898989"),  # Add grey borders between rows and columns
    fontsize = 16,  
    fontsize_row = 14,  
    fontsize_col = 14,  
    #advised to change this to an accurate descriptor
    main = "Heatmap rank-based logFC. Threshold: p < 0.05 and logFC > 1 or < -1" 
  )
  
  # Create the second heatmap
  gene.impact.heatmap2 <- pheatmap(
    t(heatmap.input[genes_group2, order(match(colnames(heatmap.input), tissue_order))]),  # Transpose the input matrix, order rows by the calculated ranks, and order columns based on tissue_order
    color = heatmap.colors, 
    breaks = heatmap.breaks,
    cluster_rows = FALSE,  # Do not perform hierarchical clustering on genes
    cluster_cols = FALSE,  # Do not perform hierarchical clustering on tissues
    border_color = c("#898989", "#898989"),  # Add grey borders between rows and columns
    legend = FALSE,  # Remove the legend
    fontsize = 16,  
    fontsize_row = 14,  
    fontsize_col = 14  
  )
  
  # Create the third heatmap
  gene.impact.heatmap3 <- pheatmap(
    t(heatmap.input[genes_group3, order(match(colnames(heatmap.input), tissue_order))]),  # Transpose the input matrix, order rows by the calculated ranks, and order columns based on tissue_order
    color = heatmap.colors, 
    breaks = heatmap.breaks,
    cluster_rows = FALSE,  # Do not perform hierarchical clustering on genes
    cluster_cols = FALSE,  # Do not perform hierarchical clustering on tissues
    border_color = c("#898989", "#898989"),  # Add grey borders between rows and columns
    legend = FALSE,  
    fontsize = 16,  
    fontsize_row = 14,  
    fontsize_col = 14  
  )
  
  # Arrange the three heatmaps vertically
  combined_heatmap <- grid.arrange(gene.impact.heatmap1$gtable, gene.impact.heatmap2$gtable, gene.impact.heatmap3$gtable, ncol = 1)
  
  # Save the combined heatmap
  # REMEMBER to change the name as to not overwrite any you previously made
  ggsave(combined_heatmap, filename = "heatmap_rank-based_logFC.png", height = 15, 
         width = 30, units = "cm")
  
}

















#### optional - looking into min_det, if that's the issue ####
# what are the sums of pct_matrix_corr, as due to min_dets some will be above 100

#add it to a column to pct_matrix_corr

pct_matrix_corr <- pct_matrix_corr %>%
  group_by(tissue, mouse) %>%
  mutate(sum_pct = sum(freq, na.rm = TRUE)) %>%
  ungroup()


pan_sums <- pct_matrix_corr %>%
  filter(tissue == "Pancreas") %>%
  select(mouse, sum_pct) %>%
  distinct() 

# Arrange in descending order and print each entry
pan_sums %>%
  arrange(desc(sum_pct)) %>%
  rowwise() %>%
  do({
    cat("Mouse:", .$mouse, "- sum_pct:", .$sum_pct, "\n")
    data.frame()  # return an empty data frame to satisfy do()
  })

#### end ####





#back to original code
## Gene tissue profile analysis ####
if ( PLOTS_ORIGINAL_PROFILE_ANALYSIS <- TRUE ) {
  
  tissue.profile.dir <- "./gene_tissue_profile/"
  
  if( !file.exists( file.path(graphs.dir, tissue.profile.dir) ) ){
    dir.create(file.path(graphs.dir, tissue.profile.dir) )
  }
  
  #if you have used ranks and non parametric testing, you can replace freq with rank to see how rank changes across tissues
  #if you want to look at that
  
  for(g in unique(pct_matrix_corr$target_gene)){
    plot.figure <- paste0(g, "_gene_tissue_profile.pdf")
    
    ggplot(pct_matrix_corr %>% filter(target_gene == g))+
      ggtitle(g)+
      scale_y_log10()+
      geom_point(aes(x = tissue, y = freq))+
      geom_line(aes(x = tissue, y = freq, group = mouse, color = mouse))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    ggsave(file.path(graphs.dir, tissue.profile.dir, plot.figure), height = 10, width = 20, units = "cm")
    
  }
}



# subset analysis within tissues-----------------------------------------
# this doesn't have the option to use ranks, to use non parametric testing yet

# assign thresholds--------- the below line only works when you only have one threshold file in your file
#thresholds <- read.csv(file = list.files(pattern = "thresholds.csv"))

#if the 3 thresholds files have the same values for cd62l cd44 and cd69 pick any to use
thresholds<- read.csv("CD8_Set1_FlowcodeDecoder_20250303_thresholds.csv")

# assign each cell as "naive", "resident" or "activated" based on:
# CD62L+CD44-, CD69+, CD62L-CD44+ or remainder

t.CD62L <- thresholds[grep("CD62L", thresholds[,1]),2]
t.CD44 <- thresholds[grep("CD44", thresholds[,1]),2]
t.CD69 <- thresholds[grep("CD69", thresholds[,1]),2]

#if you want to do this analysis for just one tissue, such as spleen, use this line here
#perCellData <- perCellData %>% filter(tissue == "Spleen")

perCellData <- perCellData %>% mutate(subset = case_when(CD62L >= t.CD62L & CD44 <= t.CD44 ~ 'Naive',
                                                         CD69 >= t.CD69 ~ 'Resident',
                                                         CD69 <= t.CD69 & CD44 >= t.CD44 ~ 'Activated'))
perCellData$subset[is.na(perCellData$subset)] <- "Activated"
anyNA(perCellData$subset)

# set as factor
perCellData$subset <- factor(perCellData$subset)

# assign tissue groups 
# if you have lots of tissues and you want to group them
# or if you want to compare all tissues separately, skip this chunk, and change subsequent instances of tissue.group to tissue

unique(perCellData$tissue)

group.tissues <- c( "Brain", "Skin", "Testes", "FRT", "Kidney", "WAT" ,
                    "Lung", "Pancreas", "Liver", "Muscle" )

group.galt <- c("PP", "IEL", "LPL")

group.blood <- c("Blood")

group.bm <- c("BM")

group.lymphoid <- c( "LN", "MLN")

perCellData$tissue.group[perCellData$tissue %in% group.tissues] <- "Tissue"
perCellData$tissue.group[perCellData$tissue %in% group.blood] <- "Blood"
perCellData$tissue.group[perCellData$tissue %in% group.bm] <- "BM"
perCellData$tissue.group[perCellData$tissue %in% group.galt] <- "GALT"
perCellData$tissue.group[perCellData$tissue %in% group.lymphoid] <- "Lymphoid"
anyNA(perCellData$tissue.group)

# make a pct matrix by subset info
# both options need this 

subset_pct_matrix <- perCellData %>%
  group_by( tissue.group, mouse, subset, Id ) %>%
  filter( Id != "Other") %>%
  summarize( count = n()) %>%
  mutate(pct = round(count/sum(count)*100, 2)) %>%
  complete(Id, fill = list(count = 0, pct = 0))

#### Heatmap of guide frequency in grouped tissues - Emily ####
#optional 

tissue_order<-c(group.blood, group.bm, group.lymphoid, group.galt, group.tissues)

# for any p>0.05, set log2FC to be 0
heatmap.input <- comparison.results
heatmap.input[heatmap.input$pvalue>0.05,]$log2FC <- 0
heatmap.input <- heatmap.input %>% select(-pvalue)

# reshape data into wide format: one row per gene, one column per tissue
heatmap.input <- heatmap.input %>% pivot_wider(names_from = tissue, values_from = log2FC)
heatmap.rownames <- heatmap.input$target_gene
heatmap.input <- heatmap.input %>% select(-target_gene)
heatmap.input <- as.matrix(heatmap.input)
rownames(heatmap.input) <- heatmap.rownames

#reorder heatmap input based on tissue order

heatmap.input.ordered <- heatmap.input[, match(tissue_order, colnames(heatmap.input))] #REORDER THE HEATMAP INPUT DATAFRAME TO MATCH THE COLUMNS DEFINED BY THE VECTOR TISSUE_ORDER

palette.length <- 50
heatmap.colors <- colorRampPalette(c("blue", "white", "red"))(palette.length)
heatmap.breaks <- c(seq(min(heatmap.input), 0, length.out=ceiling(palette.length/2) + 1), 
                    seq(max(heatmap.input)/palette.length, max(heatmap.input), length.out=floor(palette.length/2)))


#gene.impact.heatmap <- pheatmap(
# t(heatmap.input[ranked_genes$target_gene, order(match(colnames(heatmap.input), tissue_order))]), #COPIED FRPOM MAGDA
# color = heatmap.colors, 
#breaks = heatmap.breaks,
#cluster_rows = FALSE, cluster_cols = FALSE)

gene.impact.heatmap <- pheatmap(t(heatmap.input.ordered), color = heatmap.colors, #MADE FROM NEW ORDERED HEATMAP DATA FRAME
                                breaks = heatmap.breaks,
                                cluster_rows = FALSE, cluster_cols = FALSE)

ggsave(gene.impact.heatmap, filename = paste0(graphs.dir, "Frequency_heatmap_TISSUEORDER.png"), height = 10, 
       width = 20, units = "cm")

#### end ####


#### Optional, if you want the volcanos and subsequent graphs to group tissues, instead of every tissue vs spleen ####
# assign tissue groups
unique(perCellData$tissue)

group.tissues <- c( "Brain", "Skin", "Testes", "FRT", "Kidney", "WAT" ,
                    "Lung", "Pancreas", "Liver", "Muscle" )

group.galt <- c("IEL", "PP", "LPL")

group.blood <- c("Blood")

group.bm <- c("BM")

group.lymphoid <- c("Spleen", "aLN", "MLN") # this has the spleen as compared to before


perCellData$tissue.group[perCellData$tissue %in% group.tissues] <- "Tissue"
perCellData$tissue.group[perCellData$tissue %in% group.blood] <- "Blood"
perCellData$tissue.group[perCellData$tissue %in% group.bm] <- "BM"
perCellData$tissue.group[perCellData$tissue %in% group.galt] <- "GALT"
perCellData$tissue.group[perCellData$tissue %in% group.lymphoid] <- "Lymphoid"
anyNA(perCellData$tissue.group)

subset_pct_matrix <- perCellData %>%
  group_by( tissue.group, mouse, subset, Id ) %>%
  filter( Id != "Other") %>%
  summarize( count = n()) %>%
  mutate(pct = round(count/sum(count)*100, 2)) %>%
  complete(Id, fill = list(count = 0, pct = 0))
#### end ####

## within Tissue group subset comparison----------------

subset.dir <- "./tissue_subsets/"
if( !file.exists( file.path(graphs.dir, subset.dir) ) ){
  dir.create(file.path(graphs.dir, subset.dir) )
}

#for comparisons between all the groups, the following bit of code needs to be run twice,
# up until writing the comparison csv
# once with the below line and "if(s != "Resident"){" line set to Resident, and once with Naive
# remember to change the title of the csv file you save to Resident and Naive
reference_subset <- "Resident"

# filter to Tissue group
tissue_subset_matrix <- subset_pct_matrix %>% filter(tissue.group == "Tissue")

# for subset != Resident, run comparison as below (only 2 comparisons)

# set up collection of p.vals and log2FC
### --- Checked to be edited the same way the first paired t test above was edited
### --- Also needed for exporting subset.results for other graphs you make
subset.results <- data.frame()

reference.subset.matrix <- tissue_subset_matrix %>% filter(subset == reference_subset)

for (s in unique(tissue_subset_matrix$subset)) {
  
  if(s != "Resident"){
    each.subset.dir <- paste0(graphs.dir, "tissue_subsets/", s, "/")
    
    if( !file.exists(each.subset.dir)){
      dir.create(each.subset.dir)
    }
    
    soi = tissue_subset_matrix %>% filter(subset == s)
    
    resdf = data.frame(target_gene = sort(unique(soi$Id)), log2FC = 0, pvalue = 1 )
    
    for (i in 1:nrow(resdf)) {
      g = resdf$target_gene[i]
      xvdf = (soi %>% filter(Id == g) %>% arrange(mouse))
      yvdf = reference.subset.matrix %>% filter(Id == g)%>% filter(mouse %in% xvdf$mouse) %>% arrange(mouse)
      xv = (xvdf %>% filter(mouse %in% yvdf$mouse))$pct
      yv = yvdf$pct
      
      if (length(c(xv,yv)) >2){
        resdf$log2FC[i] = mean(log2((xv + 0.001) / (yv + 0.001)))
        
        res <- t.test(x = log2(xv + 0.001), y = log2(yv + 0.001), 
                      paired = TRUE, alternative = "two.sided")
        resdf$pvalue[i] = res$p.value
        
        ggplot(data.frame(soi = xv, resident = yv, id = 1:length(xv)) %>% 
                 pivot_longer(cols = -id, names_to = "group", values_to = "freq")) +
          scale_x_discrete(name = "", labels = c(reference_subset, s)) +
          scale_y_log10() +
          scale_fill_manual(values = c("darkgrey", "firebrick1")) +
          ggtitle(paste(g, ":", s, "vs", reference_subset, "\n log2FC =", round(resdf$log2FC[i], 3), "; pvalue =", round(res$p.value, 3))) +
          geom_boxplot(aes(x = group, y = freq, fill = group), alpha = 0.2) +
          geom_point(aes(x = group, y = freq)) +
          geom_line(aes(x = group, y = freq, group = id)) +
          theme_bw() +
          theme(legend.position = "none")
        
        ggsave(paste0(each.subset.dir, s, "_vs_", reference_subset, "_", g, ".pdf"), height = 10, width = 10, units = "cm")
      }
    }
    
    resdf$subset <- rep(s)
    
    # Apply multiple hypothesis testing correction
    resdf$pvalue <- p.adjust(resdf$pvalue, method = "BH")
    
    maxmlog10pvalue <- max(-log10(resdf$pvalue))
    maxlog2FC <- max(abs(resdf$log2FC))
    
    resdf$dist0 <- sqrt(resdf$log2FC^2 + (-log10(resdf$pvalue) * maxlog2FC / maxmlog10pvalue)^2)
    resdf <- resdf %>% arrange(-dist0)
    
    subset.results <- rbind(subset.results, resdf)
    write.csv(resdf, paste0(each.subset.dir, s, "_vs_", reference_subset, "_pvals.csv"))
    
  }
  
}
write.csv(subset.results, file = "cd8.subset.results.vs.res.csv")

#plotting the graphs by activation subset
subset.results.vs.res <- read.csv("cd8.subset.results.vs.res.csv")
subset.results.vs.nai <- read.csv("cd8.subset.results.vs.nai.csv")


subset.results.vs.nai <- subset.results.vs.nai %>% select (-X)
subset.results.vs.res <- subset.results.vs.res %>% select (-X)





# heatmap of gene impacts by subset---------
#### HEATMAP using the vs naive ####
#change p value if you want


heatmap.input <- subset.results.vs.res

#this sets all p values over 0.01 to white, to be represnted as not significant
heatmap.input[heatmap.input$pvalue>0.01,]$log2FC <- 0
heatmap.input <- heatmap.input %>% select(-pvalue)


# Example: Aggregating duplicate values using mean
heatmap.input <- heatmap.input %>%
  group_by(target_gene, subset) %>%
  summarise(log2FC = mean(log2FC)) %>%
  ungroup()



# Calculate the average log2FC for each target_gene across all tissues
average_log2FC_between_subsets <- aggregate(log2FC ~ target_gene, data = heatmap.input, FUN = mean)

# Rank the target genes based on the average log2FC
# This will naturally order from most negative to most positive
rank_log2FC <- rank(average_log2FC_between_subsets$log2FC)

# Create a data frame with target_gene and rank_log2FC
ranked_genes <- data.frame(target_gene = average_log2FC_between_subsets$target_gene, 
                           rank_log2FC = rank_log2FC,
                           avg_log2FC = average_log2FC_between_subsets$log2FC)

# Sort the ranked genes data frame based on rank_log2FC
# This will order from most negative to most positive log2FC
ranked_genes <- ranked_genes[order(ranked_genes$avg_log2FC), ]


# reshape data into wide format: one row per gene, one column per subsets
heatmap.input <- heatmap.input %>% pivot_wider(names_from = subset, values_from = log2FC)
heatmap.rownames <- heatmap.input$target_gene
heatmap.input <- heatmap.input %>% select(-target_gene)
heatmap.input <- as.matrix(heatmap.input)
rownames(heatmap.input) <- heatmap.rownames

heatmap.input[is.na(heatmap.input)] <- 0

palette.length <- 50
heatmap.colors <- colorRampPalette(c("blue", "white", "red"))(palette.length)
heatmap.breaks <- c(seq(min(heatmap.input, na.rm = TRUE), 0, length.out=ceiling(palette.length/2) + 1), 
                    seq(max(heatmap.input, na.rm = TRUE)/palette.length, max(heatmap.input, na.rm = TRUE), length.out=floor(palette.length/2)))


# Define order of tissues based on groups
subset_order <- c("Activated","Resident")

#code i added cos "breaks" are not unique
heatmap.breaks <- unique(c(seq(min(heatmap.input), 0, length.out = ceiling(palette.length / 2) + 1), 
                           seq(max(heatmap.input) / palette.length, max(heatmap.input), length.out = floor(palette.length / 2))))

#----------- Option1 - Create heatmap with clustered columns based on designated groups
gene.impact.heatmap <- pheatmap(
  t(heatmap.input[ranked_genes$target_gene, order(match(colnames(heatmap.input), subset_order))]),  # Transpose the input matrix, order rows by the calculated ranks, and order columns based on tissue_order
  color = heatmap.colors, 
  breaks = heatmap.breaks,
  cluster_rows = FALSE,  # Do not perform hierarchical clustering on genes
  cluster_cols = FALSE,  # Do not perform hierarchical clustering on tissues
  border_color = c("black", "black")  # Add black borders between rows and columns
)

ggsave(gene.impact.heatmap, filename = "Activation_CD8_P0.01_heatmap.png", height = 10, 
       width = 60, units = "cm")

ggsave(gene.impact.heatmap, filename = "Homing_CD8_P0.01_heatmap.pdf", height = 10, 
       width = 60, units = "cm")

#------------ Option2 - 2 rows of heatmaps

library(gridExtra)

# Split the genes into two groups
half_length <- ceiling(length(ranked_genes$target_gene) / 2)
genes_group1 <- ranked_genes$target_gene[1:half_length]
genes_group2 <- ranked_genes$target_gene[(half_length + 1):length(ranked_genes$target_gene)]

# Create the first heatmap
gene.impact.heatmap1 <- pheatmap(
  t(heatmap.input[genes_group1, order(match(colnames(heatmap.input), subset_order))]),  # Transpose the input matrix, order rows by the calculated ranks, and order columns based on tissue_order
  color = heatmap.colors, 
  breaks = heatmap.breaks,
  cluster_rows = FALSE,  
  cluster_cols = FALSE,  
  border_color = c("#898989", "#898989") 
)

# Create the second heatmap
gene.impact.heatmap2 <- pheatmap(
  t(heatmap.input[genes_group2, order(match(colnames(heatmap.input), subset_order))]),  # Transpose the input matrix, order rows by the calculated ranks, and order columns based on tissue_order
  color = heatmap.colors, 
  breaks = heatmap.breaks,
  cluster_rows = FALSE,  
  cluster_cols = FALSE,  
  border_color = c("#898989", "#898989"),  
  legend = FALSE  # Remove the legend
)

# Arrange the two heatmaps vertically
combined_heatmap <- grid.arrange(gene.impact.heatmap1$gtable, gene.impact.heatmap2$gtable, ncol = 1)

# Save the combined heatmap
# REMEMBER to change the name as to not overwrite any you previously made
ggsave(combined_heatmap, filename = "CD8_P0.01_heatmap_VS_RES.png", height = 8, 
       width = 35, units = "cm")


#### end of heatmap ####



#### volcanos, activation subsets ####
#volcanos---------------------------------------Volcanos
#p=0.01 but can change here, requires changing 4 instances of the number, in the function and in the titles for clarity



#use the correct vs data, and reference subset
#---so change as needed
subset.results<-subset.results.vs.nai
reference_subset<- "Naive"

# Calculate maximum values for scaling the plot
maxmlog10pvalue <- max(-log10(subset.results$pvalue), na.rm = TRUE)
maxlog2FC <- max(abs(subset.results$log2FC), na.rm = TRUE)

# Compute distance from (0,0) for label prioritization
subset.results$dist0 <- sqrt(subset.results$log2FC^2 + 
                               (-log10(subset.results$pvalue) * maxlog2FC / maxmlog10pvalue)^2)

options(ggrepel.max.overlaps = Inf)

#option 1 - no log2fc threshold

for (s in unique(subset.results$subset)) {
  
  plot.figure <- paste0(s, "_vs_", reference_subset, "_volcano p=0.01.pdf" )
  
  resdf <- subset.results %>% filter(subset == s) %>% arrange(-dist0)
  
  plot<-ggplot(resdf) +
    ggtitle(paste(s, "vs", reference_subset, "Volcano Plot p=0.01")) +
    scale_x_continuous(limits = c(-maxlog2FC - 0.5, maxlog2FC + 0.5)) +
    scale_y_continuous(limits = c(0, maxmlog10pvalue + 1)) +
    geom_point(aes(x = log2FC, y = -log10(pvalue)), color = "black") +
    geom_hline(yintercept = -log10(0.01), color = "firebrick1") +
    geom_text_repel(data = resdf %>% filter(pvalue < 0.01), 
                    aes(x = log2FC, y = -log10(pvalue), label = target_gene),
                    nudge_x = -0.5,
                    segment.color = "#AAAAAC",
                    segment.size = 0.5) +
    theme_bw()
  
  
  ggsave(filename = plot.figure, plot = plot, height = 10, width = 10, units = "cm")
}


#option 2 - yes log2fc threshold

subset.results<-subset.results.vs.nai
reference_subset<- "Naive"

# Calculate maximum values for scaling the plot
maxmlog10pvalue <- max(-log10(subset.results$pvalue), na.rm = TRUE)
maxlog2FC <- max(abs(subset.results$log2FC), na.rm = TRUE)

# Compute distance from (0,0) for label prioritization
subset.results$dist0 <- sqrt(subset.results$log2FC^2 + 
                               (-log10(subset.results$pvalue) * maxlog2FC / maxmlog10pvalue)^2)

options(ggrepel.max.overlaps = Inf)


for (s in unique(subset.results$subset)) {
  
  plot.figure <- paste0(s, "_vs_", reference_subset, "_volcano_log2fc.pdf")
  
  resdf <- subset.results %>% filter(subset == s) %>% arrange(-dist0)
  
  p <- ggplot(resdf) +
    ggtitle(paste(s, "vs", reference_subset, "Volcano log2fc=1 p=0.01")) +
    scale_x_continuous(limits = c(-maxlog2FC - 0.5, maxlog2FC + 0.5)) +
    scale_y_continuous(limits = c(0, maxmlog10pvalue + 1)) +
    geom_point(aes(x = log2FC, y = -log10(pvalue)), color = "black") +
    geom_hline(yintercept = -log10(0.01), color = "firebrick1") +
    geom_vline(xintercept = c(-1, 1), color = "red", linetype = "dashed") +
    geom_text_repel(data = resdf %>% filter(pvalue < 0.01 & (log2FC <= -1 | log2FC >= 1)), 
                    aes(x = log2FC, y = -log10(pvalue), label = target_gene),
                    nudge_x = -0.5,
                    segment.color = "#AAAAAC",
                    segment.size = 0.5) +
    theme_bw()
  
  ggsave(filename = plot.figure, plot = p, height = 10, width = 10, units = "cm")
}

#option 3 - only plots the genes of interest, regardless of log2fc and p values

subset.results<-subset.results.vs.nai
reference_subset<- "Naive"

# Compute distance from (0,0) for label prioritization
subset.results$dist0 <- sqrt(subset.results$log2FC^2 + 
                               (-log10(subset.results$pvalue) * maxlog2FC / maxmlog10pvalue)^2)

options(ggrepel.max.overlaps = Inf)

# List of genes to label, edit however you want
genes_to_label <- c("Itga4", "Itga9", "Cxcr3", "Itgb1", "Vcam1")

for (s in unique(subset.results$subset)) {
  
  plot.figure <- paste0(s, "_vs_", reference_subset, "_volcano_hits.pdf")
  
  resdf <- subset.results %>% filter(subset == s) %>% arrange(-dist0)
  
  p <- ggplot(resdf) +
    ggtitle(paste(s, "vs", reference_subset, "Volcano hits")) +
    scale_x_continuous(limits = c(-maxlog2FC - 0.5, maxlog2FC + 0.5)) +
    scale_y_continuous(limits = c(0, maxmlog10pvalue + 1)) +
    geom_point(data = resdf %>% filter(target_gene %in% genes_to_label), 
               aes(x = log2FC, y = -log10(pvalue)), color = "black") +
    geom_hline(yintercept = -log10(0.01), color = "firebrick1") +
    geom_vline(xintercept = c(-1, 1), color = "red", linetype = "dashed") +
    geom_text_repel(data = resdf %>% filter(target_gene %in% genes_to_label), 
                    aes(x = log2FC, y = -log10(pvalue), label = target_gene),
                    nudge_x = -0.5,
                    segment.color = "#AAAAAC",
                    segment.size = 0.5) +
    theme_bw()
  
  ggsave(filename = plot.figure, plot = p, height = 10, width = 10, units = "cm")
}


#option 4 - plots all genes, only labels genes of interest. should change code to change datapoint colour
subset.results<-subset.results.vs.nai
reference_subset<- "Naive"

# Compute distance from (0,0) for label prioritization
subset.results$dist0 <- sqrt(subset.results$log2FC^2 + 
                               (-log10(subset.results$pvalue) * maxlog2FC / maxmlog10pvalue)^2)

options(ggrepel.max.overlaps = Inf)

# List of genes to label
genes_to_label <- c("Itga4", "Itga9", "Cxcr3", "Itgb1", "Vcam1")

for (s in unique(subset.results$subset)) {
  
  plot.figure <- paste0(s, "_vs_", reference_subset, "_volcano_opt4.pdf")
  
  resdf <- subset.results %>% filter(subset == s) %>% arrange(-dist0)
  
  p <- ggplot(resdf) +
    ggtitle(paste(s, "vs", reference_subset, "Volcano opt4")) +
    scale_x_continuous(limits = c(-maxlog2FC - 0.5, maxlog2FC + 0.5)) +
    scale_y_continuous(limits = c(0, maxmlog10pvalue + 1)) +
    geom_point(aes(x = log2FC, y = -log10(pvalue)), color = "black") +
    geom_hline(yintercept = -log10(0.01), color = "firebrick1") +
    geom_vline(xintercept = c(-1, 1), color = "red", linetype = "dashed") +
    geom_text_repel(data = resdf %>% filter(target_gene %in% genes_to_label), 
                    aes(x = log2FC, y = -log10(pvalue), label = target_gene),
                    nudge_x = -0.5,
                    segment.color = "#AAAAAC",
                    segment.size = 0.5) +
    theme_bw()
  
  ggsave(filename = plot.figure, plot = p, height = 10, width = 10, units = "cm")
}


####end of volcano code####





#### heatmap for pancreas vs spleen, within activated and resident subsets ####
# t-test to detect guides that are significantly enriched or depleted in the pancreas in comparison to the reference tissue the spleen
# this is done within the populations of activated and resident to indicate which guides are important for which subset

# its parametric, does not use ranks

#the percelldata needs to have been assigned naive/activated/resident previously
perCellData <- perCellData %>% filter (Id!="Other")


# the tissue t test needs data in the format pct matrix corr so subset it

naive_perCells <- perCellData %>% filter (subset=="Naive")
activated_perCells <- perCellData %>% filter (subset=="Activated")
resident_perCells <- perCellData %>% filter (subset=="Resident")

#--- activated code ----
# make the pct matrix from the percell data
# Step 1: Count the number of cells per tissue, mouse, and Id
cell_counts <- activated_perCells %>%
  group_by(tissue, mouse, Id) %>%
  summarise(count = n(), .groups = "drop")

count_matrix <- cell_counts %>%
  pivot_wider(names_from = Id, values_from = count, values_fill = 0)

count_matrix <- count_matrix %>%
  rowwise() %>%
  mutate(total = sum(c_across(-c(tissue, mouse)))) %>%
  ungroup()

count_matrix <- count_matrix %>% filter(total != 0)

activated_pct_matrix <- count_matrix %>%
  mutate(across(-c(tissue, mouse, total), ~ (.x / total) * 100))

#remove samples below 20 cells total, then remove total column
# Subset the rows with fewer than 20 cells
removed_rows <- count_matrix[count_matrix$total < 20, ]

# Print each removed sample with mouse and tissue info and cell count
cat("Removed samples:\n")
apply(removed_rows, 1, function(row) {
  cat(row["mouse"], row["tissue"], "—", row["total"], "cells\n")
})

activated_pct_matrix <- activated_pct_matrix[count_matrix$total >= 20, ] %>%
  select(-total)

#check each row sums to 100
row_sums <- rowSums(activated_pct_matrix[ , -(1:2)], na.rm = TRUE)
# View the result
print(row_sums)

#change format to what the code expects
activated_pct_matrix <- activated_pct_matrix %>% 
  pivot_longer(cols = -c(tissue, mouse), names_to = "target_gene", values_to = "freq" )

#if you have multiple sets in the experiment (>56 guides), you need to run this code
# Filter out rows for set1
filtered_matrix <- activated_pct_matrix %>%
  filter(!(freq == 0 & mouse %in% set1_mice & !(target_gene %in% set1_guides)))

# Filter out rows for set2
filtered_matrix <- filtered_matrix %>%
  filter(!(freq == 0 & mouse %in% set2_mice & !(target_gene %in% set2_guides)))

# Filter out rows for set3
filtered_matrix <- filtered_matrix %>%
  filter(!(freq == 0 & mouse %in% set3_mice & !(target_gene %in% set3_guides)))

#for subsequent code to work
activated_pct_matrix <- filtered_matrix

dim(activated_pct_matrix)

#tissue mins

tissue_counts <- activated_perCells %>%
  group_by( tissue, mouse, Id ) %>%
  filter( Id != "Other") %>%
  summarize( count = n()) %>%
  mutate(pct = round(count/sum(count)*100, 2)) 

tissue_mins <- tissue_counts %>%
  group_by( tissue, mouse ) %>%
  summarize( total = sum(count) ) %>%
  mutate( min_det = 100/total )

# calculate half-min frequency to replace zeros in percent matrix

tissue_mins$half_min <- tissue_mins$min_det/2


#back to normal code
tissue_mins <- tissue_mins %>% unite(col = sample, tissue, mouse)


pct_matrix_corrected <- activated_pct_matrix
pct_matrix_corrected <- pct_matrix_corrected %>% unite(col = sample, tissue, mouse)


for (s in unique(pct_matrix_corrected$sample)) {
  min_temp <- filter(tissue_mins, sample == s)
  pct_temp <- filter(pct_matrix_corrected, sample == s)
  
  if (nrow(min_temp) == 1) {
    pct_temp$freq[pct_temp$freq == 0] <- min_temp$half_min
    pct_matrix_corrected$freq[pct_matrix_corrected$sample %in% unique(pct_temp$sample)] <- pct_temp$freq
  }
}

# Cell number cutoff has been done previously if you made the pct_matrix from the perCell Data

pct_matrix_corr <- left_join(pct_matrix_corrected, tissue_mins, by = "sample")

#if it hasn't, filter out the <20 samples

filtered_out <- pct_matrix_corr %>% filter(total <= 20)

cat("Samples filtered out due to total <= 20:\n")
print(unique(filtered_out$sample)) 

pct_matrix_corr <- pct_matrix_corr %>% filter(total > 20)

pct_matrix_corr <- pct_matrix_corr %>% separate(sample, remove = TRUE, sep = "_", into = c("tissue", "mouse"))  %>%
  select(-total, -min_det, -half_min)

anyNA(pct_matrix_corr)

activated_pct_matrix_corr <- pct_matrix_corr %>%
  filter(tissue %in% c("Spleen", "Pancreas"))


#now run it through a paired t test, pancreas compared to spleen
# Simple tissue vs Spleen freq comparison----------------
graphs.dir <- paste0("./", flowcode.seed, "_graphs/")

if( !file.exists(graphs.dir) ){
  dir.create(graphs.dir)
}

comparison.dir <- "pancreas_vs_Spleen/"
if( !file.exists( file.path(graphs.dir, comparison.dir) ) ){
  dir.create(file.path(graphs.dir, comparison.dir) )
}

reference_tissue <- "Spleen"


#  normal version using mean log2fc, set up collection of p.vals and log2FC
comparison.results <- data.frame()

for (t in unique(activated_pct_matrix_corr$tissue)) {
  
  if(t != "Spleen"){
    tissue.dir <- paste0(graphs.dir, comparison.dir, t, "/")
    
    if( !file.exists(tissue.dir)){
      dir.create(tissue.dir)
    }
    
    toi = activated_pct_matrix_corr %>% filter(tissue == t)
    
    spleen_toi = activated_pct_matrix_corr %>% filter(tissue == reference_tissue)
    
    resdf = data.frame(target_gene = sort(unique(toi$target_gene)), log2FC = 0, pvalue = 1 )
    
    for (i in 1:nrow(resdf)) {
      g = resdf$target_gene[i]
      xvdf = (toi %>% filter(target_gene == g) %>% arrange(mouse))
      yvdf = spleen_toi %>% filter(target_gene == g)%>% filter(mouse %in% xvdf$mouse) %>% arrange(mouse)
      xv = (xvdf %>% filter(mouse %in% yvdf$mouse))$freq
      yv = yvdf$freq
      
      if (length(c(xv,yv)) >2){
        resdf$log2FC[i] = mean(log2((xv + 0.001) / (yv + 0.001)))
        
        res <- t.test(x = log2(xv + 0.001), y = log2(yv + 0.001), 
                      paired = TRUE, alternative = "two.sided")
        resdf$pvalue[i] = res$p.value
        
        ggplot(data.frame(toi = xv, spleen = yv, id = 1:length(xv)) %>% 
                 pivot_longer(cols = -id, names_to = "group", values_to = "freq")) +
          scale_x_discrete(name = "", labels = c(reference_tissue, t)) +
          scale_y_log10() +
          scale_fill_manual(values = c("darkgrey", "firebrick1")) +
          ggtitle(paste(g, ":", t, "vs", reference_tissue, "\n log2FC =", round(resdf$log2FC[i], 3), "; pvalue =", round(res$p.value, 3))) +
          geom_boxplot(aes(x = group, y = freq, fill = group), alpha = 0.2) +
          geom_point(aes(x = group, y = freq)) +
          geom_line(aes(x = group, y = freq, group = id)) +
          theme_bw() +
          theme(legend.position = "none")
        
        ggsave(paste0(tissue.dir, t, "_vs_", reference_tissue, "_", g, ".pdf"), height = 10, width = 10, units = "cm")
      }
    }
    
    
    
    resdf$tissue <- rep(t)
    
    # Apply multiple hypothesis testing correction
    resdf$pvalue <- p.adjust(resdf$pvalue, method = "BH")
    
    maxmlog10pvalue <- max(-log10(resdf$pvalue))
    maxlog2FC <- max(abs(resdf$log2FC))
    
    resdf$dist0 <- sqrt(resdf$log2FC^2 + (-log10(resdf$pvalue) * maxlog2FC / maxmlog10pvalue)^2)
    resdf <- resdf %>% arrange(-dist0)
    
    comparison.results <- rbind(comparison.results, resdf)
    write.csv(resdf, paste0(tissue.dir, t, "_vs_", reference_tissue, "_pvals.csv"))
    
  }
  
}




#  Remove the column 'distr0' from comparison.results
comparison.results <- comparison.results %>% select(-dist0)

#export out the comparison results as a csv for use in other scripts
#change naming to be descriptive
write.csv(comparison.results, "CD8_activated_comparison_results.csv", row.names = FALSE)

#save it
activated.comparison.results <- comparison.results







#### ------ segment for resident ----- ####

# make the pct matrix from the percell data
# Step 1: Count the number of cells per tissue, mouse, and Id
cell_counts <- resident_perCells %>%
  group_by(tissue, mouse, Id) %>%
  summarise(count = n(), .groups = "drop")

count_matrix <- cell_counts %>%
  pivot_wider(names_from = Id, values_from = count, values_fill = 0)

count_matrix <- count_matrix %>%
  rowwise() %>%
  mutate(total = sum(c_across(-c(tissue, mouse)))) %>%
  ungroup()

count_matrix <- count_matrix %>% filter(total != 0)

resident_pct_matrix <- count_matrix %>%
  mutate(across(-c(tissue, mouse, total), ~ (.x / total) * 100))

#remove samples below 20 cells total, then remove total column
# Subset the rows with fewer than 20 cells
removed_rows <- count_matrix[count_matrix$total < 20, ]

# Print each removed sample with mouse and tissue info and cell count
cat("Removed samples:\n")
apply(removed_rows, 1, function(row) {
  cat(row["mouse"], row["tissue"], "—", row["total"], "cells\n")
})

resident_pct_matrix <- resident_pct_matrix[count_matrix$total >= 20, ] %>%
  select(-total)

#change format to what the code expects
resident_pct_matrix <- resident_pct_matrix %>% 
  pivot_longer(cols = -c(tissue, mouse), names_to = "target_gene", values_to = "freq" )


#minimum detection, will be from the subset of interest

tissue_counts <- resident_perCells %>%
  group_by( tissue, mouse, Id ) %>%
  filter( Id != "Other") %>%
  summarize( count = n()) %>%
  mutate(pct = round(count/sum(count)*100, 2)) 

tissue_mins <- tissue_counts %>%
  group_by( tissue, mouse ) %>%
  summarize( total = sum(count) ) %>%
  mutate( min_det = 100/total )

# calculate half-min frequency to replace zeros in percent matrix

tissue_mins$half_min <- tissue_mins$min_det/2


#back to normal code
tissue_mins <- tissue_mins %>% unite(col = sample, tissue, mouse)

pct_matrix_corrected <- resident_pct_matrix
pct_matrix_corrected <- pct_matrix_corrected %>% unite(col = sample, tissue, mouse)
#added
#pct_matrix_corrected <- pct_matrix_corrected %>% rename (freq=pct)

for (s in unique(pct_matrix_corrected$sample)) {
  min_temp <- filter(tissue_mins, sample == s)
  pct_temp <- filter(pct_matrix_corrected, sample == s)
  
  if (nrow(min_temp) == 1) {
    pct_temp$freq[pct_temp$freq == 0] <- min_temp$half_min
    pct_matrix_corrected$freq[pct_matrix_corrected$sample %in% unique(pct_temp$sample)] <- pct_temp$freq
  }
}

# Cell number cutoff has been done previously if you made the pct_matrix from the perCell Data

pct_matrix_corr <- left_join(pct_matrix_corrected, tissue_mins, by = "sample")

#if it hasn't, filter out the <20 samples

filtered_out <- pct_matrix_corr %>% filter(total <= 20)

cat("Samples filtered out due to total <= 20:\n")
print(unique(filtered_out$sample)) 

pct_matrix_corr <- pct_matrix_corr %>% filter(total > 20)

pct_matrix_corr <- pct_matrix_corr %>% separate(sample, remove = TRUE, sep = "_", into = c("tissue", "mouse"))  %>%
  select(-total, -min_det, -half_min)

anyNA(pct_matrix_corr)

resident_pct_matrix_corr <- pct_matrix_corr %>%
  filter(tissue %in% c("Spleen", "Pancreas"))


#now run it through a paired t test, pancreas compared to spleen
# Simple tissue vs Spleen freq comparison----------------
graphs.dir <- paste0("./", flowcode.seed, "_graphs/")

if( !file.exists(graphs.dir) ){
  dir.create(graphs.dir)
}

comparison.dir <- "pancreas_vs_Spleen/"
if( !file.exists( file.path(graphs.dir, comparison.dir) ) ){
  dir.create(file.path(graphs.dir, comparison.dir) )
}

reference_tissue <- "Spleen"


#  normal version using mean log2fc, set up collection of p.vals and log2FC
comparison.results <- data.frame()

for (t in unique(resident_pct_matrix_corr$tissue)) {
  
  if(t != "Spleen"){
    tissue.dir <- paste0(graphs.dir, comparison.dir, t, "/")
    
    if( !file.exists(tissue.dir)){
      dir.create(tissue.dir)
    }
    
    toi = resident_pct_matrix_corr %>% filter(tissue == t)
    
    spleen_toi = resident_pct_matrix_corr %>% filter(tissue == reference_tissue)
    
    resdf = data.frame(target_gene = sort(unique(toi$target_gene)), log2FC = 0, pvalue = 1 )
    
    for (i in 1:nrow(resdf)) {
      g = resdf$target_gene[i]
      xvdf = (toi %>% filter(target_gene == g) %>% arrange(mouse))
      yvdf = spleen_toi %>% filter(target_gene == g)%>% filter(mouse %in% xvdf$mouse) %>% arrange(mouse)
      xv = (xvdf %>% filter(mouse %in% yvdf$mouse))$freq
      yv = yvdf$freq
      
      if (length(c(xv,yv)) >2){
        resdf$log2FC[i] = mean(log2((xv + 0.001) / (yv + 0.001)))
        
        res <- t.test(x = log2(xv + 0.001), y = log2(yv + 0.001), 
                      paired = TRUE, alternative = "two.sided")
        resdf$pvalue[i] = res$p.value
        
        ggplot(data.frame(toi = xv, spleen = yv, id = 1:length(xv)) %>% 
                 pivot_longer(cols = -id, names_to = "group", values_to = "freq")) +
          scale_x_discrete(name = "", labels = c(reference_tissue, t)) +
          scale_y_log10() +
          scale_fill_manual(values = c("darkgrey", "firebrick1")) +
          ggtitle(paste(g, ":", t, "vs", reference_tissue, "\n log2FC =", round(resdf$log2FC[i], 3), "; pvalue =", round(res$p.value, 3))) +
          geom_boxplot(aes(x = group, y = freq, fill = group), alpha = 0.2) +
          geom_point(aes(x = group, y = freq)) +
          geom_line(aes(x = group, y = freq, group = id)) +
          theme_bw() +
          theme(legend.position = "none")
        
        ggsave(paste0(tissue.dir, t, "_vs_", reference_tissue, "_", g, ".pdf"), height = 10, width = 10, units = "cm")
      }
    }
    
    
    
    resdf$tissue <- rep(t)
    
    # Apply multiple hypothesis testing correction
    resdf$pvalue <- p.adjust(resdf$pvalue, method = "BH")
    
    maxmlog10pvalue <- max(-log10(resdf$pvalue))
    maxlog2FC <- max(abs(resdf$log2FC))
    
    resdf$dist0 <- sqrt(resdf$log2FC^2 + (-log10(resdf$pvalue) * maxlog2FC / maxmlog10pvalue)^2)
    resdf <- resdf %>% arrange(-dist0)
    
    comparison.results <- rbind(comparison.results, resdf)
    write.csv(resdf, paste0(tissue.dir, t, "_vs_", reference_tissue, "_pvals.csv"))
    
  }
  
}




#  Remove the column 'distr0' from comparison.results
comparison.results <- comparison.results %>% select(-dist0)

#export out the comparison results as a csv for use in other scripts
#change naming to be descriptive
write.csv(comparison.results, "CD8_resident_comparison_results.csv", row.names = FALSE)

#save it
resident.comparison.results <- comparison.results



#### ----- naive segment ----- #####

#i did run this but you may get no result, or nothing useful, as the cell numbers can be small
# if you want to run it you can copy the previous section for activated or resident and change those occurances to naive
# would require some adjustments to the code to plot the heatmap


#### ---- plot the heatmap---- ####

#if you need to load it
resident.comparison.results <- read.csv("CD8_resident_comparison_results.csv")
activated.comparison.results <- read.csv("CD8_activated_comparison_results.csv")

#trying to be neater
resident.results <- resident.comparison.results %>% mutate(subset="Resident")
activated.results <- activated.comparison.results %>% mutate(subset="Activated")

#reformatting the data frames for heatmap
# heatmap specific code

# - activated data
# uses a p value of 0.01, feel free to adjust
activated.results[activated.results$pvalue>0.01,]$log2FC <- 0
activated.results <- activated.results %>% select(-pvalue)

activated.results <- activated.results %>% 
  rename (Activated=log2FC) %>%
  select(-tissue,-subset)

#optional line - to only keep the genes of interest
activated.results <- activated.results %>% filter (target_gene %in% c("Itga4", "Itga9", "Itgb1", "Cxcr3", "Vcam1"))


# - resident data
# uses a p value of 0.01, feel free to adjust
resident.results[resident.results$pvalue>0.01,]$log2FC <- 0
resident.results <- resident.results %>% select(-pvalue)

resident.results <- resident.results %>% 
  rename (Resident=log2FC) %>%
  select(-tissue,-subset)

#optional line - to only keep the genes of interest
resident.results <- resident.results %>% filter (target_gene %in% c("Itga4", "Itga9", "Itgb1", "Cxcr3", "Vcam1"))

#join them
heatmap.input <- full_join(resident.results, activated.results, by = "target_gene")
# Set target_gene as row names
rownames(heatmap.input) <- heatmap.input$target_gene

# Remove the target_gene column now that it's used as row names
heatmap.input <- heatmap.input[, !names(heatmap.input) %in% "target_gene"]


# Step 1: Calculate average log2FC across Resident and Activated
average_log2FC <- rowMeans(heatmap.input, na.rm = TRUE)

# Step 2: Create ranked gene list
ranked_genes <- data.frame(
  target_gene = rownames(heatmap.input),
  avg_log2FC = average_log2FC
) %>%
  arrange(avg_log2FC)

# Step 3: Set up heatmap color palette and breaks
palette.length <- 50
heatmap.colors <- colorRampPalette(c("blue", "white", "red"))(palette.length)

heatmap.breaks <- c(
  seq(min(heatmap.input, na.rm = TRUE), 0, length.out = ceiling(palette.length / 2) + 1),
  seq(max(heatmap.input, na.rm = TRUE) / palette.length, max(heatmap.input, na.rm = TRUE), length.out = floor(palette.length / 2))
)

#if heatmap.breaks causes an error for not being unique
heatmap.breaks <- unique(heatmap.breaks)


#------ Option 1 - long heatmap
# Create one long heatmap with clustered columns based on designated groups

library(pheatmap)
library(grid)
library(gridExtra)

# Open a PNG device
png("CD8_subset_depletion_heatmap_p0.01.png", height = 5, width = 10, units = "cm", res = 300)

# Create the heatmap without a title
heatmap <- pheatmap(
  t(heatmap.input[ranked_genes$target_gene, ]),
  color = heatmap.colors,
  breaks = heatmap.breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey",
  main = "",               # Suppress default title
  fontsize = 10,
  angle_col = 45 
)

# Add relevant title 
grid.text(
  "p0.01 CD8 significant depletion of guides in panc compared to spl",
  x = 0.5, y = 0.98, gp = gpar(fontsize = 9)
)

# Close the device to save the file
dev.off()


#### end####

#### bar chart of depletion in a tissue comapred to the spleen ####


#option1 to focus on a few genes
small_pct_matrix <- pct_matrix_corr %>% filter (target_gene %in% c("Itga4", "Itga9", "Itgb1", "Cxcr3", "Vcam1"))
pan_pcts<-small_pct_matrix %>% filter (tissue =="Pancreas")
spl_pcts<-small_pct_matrix %>% filter (tissue =="Spleen")

#option2 for all genes

pan_pcts<-pct_matrix_corr %>% filter (tissue =="Pancreas")
spl_pcts<-pct_matrix_corr %>% filter (tissue =="Spleen")


#all options rejoin here
#average pancreas data

pan_pcts <- pan_pcts %>%
  group_by(target_gene) %>%
  summarise(pan_freq = mean(freq))

#average spleen data

spl_pcts <- spl_pcts %>%
  group_by(target_gene) %>%
  summarise(spl_freq = mean(freq))

#bind data

combined_pcts <- pan_pcts %>%
  left_join(spl_pcts, by = "target_gene")


combined_pcts <- combined_pcts %>%
  mutate(norm_freq = (pan_freq / spl_freq)*100)

combined_pcts <- combined_pcts %>%
  mutate(depletion_pct=100-norm_freq)





# plot barchart


# Load necessary library
library(ggplot2)

# Order the dataframe by norm_freq in ascending order

combined_pcts <- combined_pcts %>%
  arrange(desc(depletion_pct))

#remove unwanted guides
genes_to_remove <- c("Itgb2l", "Cd22", "Ctnna2", "Itga10")

# Filter the dataframe
combined_pcts <- combined_pcts[!combined_pcts$target_gene %in% genes_to_remove, ]



# Create the bar chart
p <- ggplot(combined_pcts, aes(x = reorder(target_gene, -depletion_pct), y = depletion_pct)) +
  geom_bar(stat = "identity", fill = "blue3") +
  theme_minimal() +
  theme(
    text = element_text(size = 16, color = "black"), # Customize font size and color
    axis.title = element_text(size = 18, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.line = element_line(size = 0.5, colour = "black") # Add black axis lines
  ) +
  labs(
    title = "Abundance in pancreas normalized to spleen",
    x = "Target Gene",
    y = "Relative Depletion (%)"
    
  ) +
  scale_y_continuous(limits = c(0, 100)) # Set y-axis limits

print(p)



# Save the plot as a PNG file
ggsave(filename = "normalized_frequency_plot.png", plot = p, height = 17, width = 15, units = "cm", dpi = 300, bg="white")
















## tissue group analysis---------------------




### phenotypic divergence analysis by Oliver-----------------------------

# an alternative to this part is flowcytoscript to make UMAP and tSNEs instead
# which is found on github

# select channels to use for dimensionality reduction
# tip: check pdf markers plot to look for FRET artefacts
# exclude those channels, plus any irrelevant ones
## channels to exclude:AF, BB660, ivCD45, CCR2, CD25
## run with channels: 2:3,7,10:15,18:22,26:27,29,31:33
## first, try to fix errors crudely by removing hypernegative events
## bring up values below 250 to 250

# remove untransduced cells here
perCellData <- perCellData %>% dplyr::filter(Id != "Other")

unwanted.channels <- c("SSC", "FSC", "Id", "sample", "abovethreshold", 
                       "combination", "tissue", "mouse", "subset")
channel.choices <- colnames(perCellData)[!grepl( paste0(unwanted.channels, collapse = "|"),
                                                 colnames(perCellData))]
channels.for.dr <- multi.menu(channel.choices, "Select the channels to use for visualization and clustering: ")










channels.for.dr <- channel.choices[channels.for.dr]
flowCode.dr.data <- perCellData[,channels.for.dr]
flowCode.dr.data.corr <- flowCode.dr.data
flowCode.dr.data.corr[flowCode.dr.data.corr<250] <- 250

# SOM-based clustering
cluster.n <- 10
set.seed( flowcode.seed ) 
flow.som <- EmbedSOM::SOM(flowCode.dr.data.corr, xdim = 24, 
                          ydim = 24, batch = TRUE,
                          parallel = TRUE, threads = flowcode.threads )
flow.som.mapping <- flow.som$mapping[ , 1 ]
flow.som.codes <- flow.som$codes
consensus.cluster <- ConsensusClusterPlus( t( flow.som.codes ),
                                           maxK = cluster.n, reps = 100, pItem = 0.9, 
                                           pFeature = 1,
                                           clusterAlg = "hc", innerLinkage = "average", 
                                           finalLinkage = "average",
                                           distance = "euclidean", 
                                           seed = flowcode.seed )

flow.som.event.cluster <- consensus.cluster[[ cluster.n ]]$consensusClass[ flow.som.mapping ]
flow.som.cluster.rank <- 1 + cluster.n - rank( table( flow.som.event.cluster ), ties.method = "last" )
flow.som.event.cluster <- flow.som.cluster.rank[ flow.som.event.cluster ]
names( flow.som.event.cluster ) <- NULL

# add clusters as column to data
perCellData$Cluster <- factor( flow.som.event.cluster )

# cluster distribution by flowcode per tissue
cluster.fig.name <- "flowcode_cluster_"
cluster.fig.dir <- "./cluster_histograms/"
dir.create(cluster.fig.dir)

cluster.distr <- perCellData %>% group_by( mouse, tissue, Id ) %>% 
  count(Cluster) %>%
  mutate(freq = round(100*n/sum(n),2))

cluster.distr.summary <- cluster.distr %>% group_by( tissue, Id, Cluster ) %>%
  summarise( mean = mean(freq), sd = sd(freq))

cluster.distr.summary$sd[is.na(cluster.distr.summary$sd)] <- 0

for(t in unique(cluster.distr.summary$tissue)){
  tissue.dist <- cluster.distr.summary %>% dplyr::filter(tissue == t)
  
  tissue.clust.plot <- ggplot(tissue.dist, aes(x = Cluster, y = mean, fill = Id ))+
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes( ymin = mean, ymax = mean + sd ),
                  width = 0.2 ) +
    facet_wrap(~Id) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_discrete(name = "FlowCode") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 5))
  
  ggsave(paste0(cluster.fig.dir, cluster.fig.name, t, ".jpg"), tissue.clust.plot,
         width = 10, height = 6)
}



# cosine similarity assessment----------------

dir.create("./volcano_plots")

complete.cluster.distr <- cluster.distr %>% 
  group_by(mouse, tissue) %>% 
  complete(Id, Cluster)

complete.cluster.distr[is.na(complete.cluster.distr)] <- 0

# summarize cluster distribution by tissue for reference
cluster.distr.tissue <- complete.cluster.distr %>% group_by(mouse, tissue, Cluster) %>%
  summarise( mean = mean(freq))

# calculate cosine similarity to means for each tissue
similarity.results.df <- complete.cluster.distr %>%
  inner_join(cluster.distr.tissue, by = c("mouse", "tissue", "Cluster")) %>%
  group_by(mouse, tissue, Id) %>%
  summarize(cosine.similarity = cosine(freq, mean), .groups = 'drop') %>%
  ungroup()

colnames(similarity.results.df)[4] <- "cosine"

# convert to divergence (1/cosine.similarity)
similarity.results.df$cosine <- 1/similarity.results.df$cosine

# set up collection of p.vals and mean divergence

cosine.results <- data.frame()

# compare to median divergence
#note - this has a t test on percell derived data but already has a log10, I've left it as it was, also I don't really use this code section
for (t in unique(similarity.results.df$tissue)){
  
  toi = similarity.results.df %>% filter(tissue == t)
  
  resdf = data.frame(target_gene = sort(unique(toi$Id)), divergence = 0, pval = 1, 
                     tissue = rep(t, length(unique(toi$Id))) )
  
  for(i in 1:nrow(resdf)){
    g = resdf$target_gene[i]
    xvdf = (toi %>% filter(Id == g) %>% arrange(mouse))
    xv = xvdf$cosine
    yvdf = (toi %>% filter(mouse %in% xvdf$mouse))
    yv = (yvdf %>% group_by(mouse) %>% summarise(median = median(cosine)))$median
    
    if(length(xv) >1){
      resdf$divergence[i] = mean(xv)
      
      res = t.test(x = log10(xv+0.1), 
                   y = log10(yv+0.1), 
                   paired = FALSE, 
                   alternative = "two.sided")
      resdf$pval[i] = res$p.value
    }
    
  }
  cosine.results = rbind(cosine.results, resdf)
}

cosine.results$divergence[cosine.results$divergence==0] <- 1

min.divergence = min(cosine.results$divergence)
max.divergence = max(cosine.results$divergence)
maxlog10.pval = max(-log10(cosine.results$pval))

for (t in unique(cosine.results$tissue)){
  plotdf <- cosine.results %>% dplyr::filter(tissue == t)
  
  plotdf$dist0 = sqrt(plotdf$divergence^2+(-log10(plotdf$pval)*max.divergence/maxlog10.pval)^2 )
  
  plotdf = plotdf %>% arrange(-dist0)
  
  ggplot()+
    ggtitle(paste(t, "Phenotypic divergence volcano"))+
    scale_x_continuous(limits = c(min.divergence, max.divergence))+
    scale_y_continuous(limits = c(0, maxlog10.pval))+
    geom_point(data = plotdf, aes(x = divergence, y=-log10(pval)))+
    geom_hline(yintercept = -log10(0.05), color = "firebrick1")+
    geom_text_repel(data = plotdf %>% filter(pval < 0.05), 
                    aes(x = divergence, y = -log10(pval), label = target_gene),
                    nudge_x = -0.5) +
    theme_bw()
  ggsave(paste0("./volcano_plots/", t, "_phenotype_volcano.pdf"),
         height = 10, width = 10, units = "cm")
  
}



# heatmap of phenotype impacts by tissue---------

# for any p>0.05, set log2FC to be 0
pheno.heatmap.input <- cosine.results

#to manually remove a FRET-errored guide
pheno.heatmap.input <- pheno.heatmap.input %>% filter(target_gene != "Ctnna2")

pheno.heatmap.input[pheno.heatmap.input$pval>0.05,]$divergence <- 0
pheno.heatmap.input <- pheno.heatmap.input %>% select(-pval)

# reshape data into wide format: one row per gene, one column per tissue
pheno.heatmap <- pheno.heatmap.input  %>% 
  pivot_wider(names_from = tissue, values_from = divergence)
pheno.rownames <- pheno.heatmap$target_gene
pheno.heatmap <- pheno.heatmap %>% select(-target_gene)
pheno.heatmap <- as.matrix(pheno.heatmap)

rownames(pheno.heatmap) <- pheno.rownames
pheno.heatmap[is.na(pheno.heatmap)] <- 0

pheno.heatmap.colors <- colorRampPalette(c("white", "red"))(palette.length)
pheno.heatmap.breaks <- c( seq(max(pheno.heatmap)/palette.length, max(pheno.heatmap), length.out=floor(palette.length)))


pheno.heatmap.plot <- pheatmap(t(pheno.heatmap), color = pheno.heatmap.colors, 
                               breaks = pheno.heatmap.breaks,
                               fontsize_col = 9,
                               cluster_rows = FALSE, cluster_cols = FALSE )

ggsave(pheno.heatmap.plot, filename = paste0(graphs.dir, "Pheno_heatmap.png"), height = 10, 
       width = 50, units = "cm")



