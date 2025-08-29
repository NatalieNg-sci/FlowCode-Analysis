# Comparing 2 groups of mice - volcano plot - FlowCode Analysis - Natalie
# parametric testing, with the option of non parametric testing if you want

# to compare more groups to each other, change these group names and input files as desired


#---  adapting for copd vs healthy ---
# meaning for two separate groups of mice, a non paired t test. and for lung normalized to spleen

#last edited 29th August, Natalie Ng - to make presentable


# --- Load packages

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
flowcode.seed <- 20250819
flowcode.threads <- parallelly::availableCores() - 1
data.table::setDTthreads(flowcode.threads)





#### Process Data ####
# --- Assign which groups you are comparing
# probably best to have group1 as the group of interest and group2 as the reference one
# regardless, group2 will be used as reference.

#also i am going to edit this as suitable for the COPD data
# for T cell type vs t cell type functionality use that version of the code

Group1 <- "Elastase"
Group2 <-  "PBSCTRL"

# --- Load the input data
#this code allows for multiple sets, you can change it if you just have one library

#load the percell data
perCellData <- read.csv("CD8_FlowcodeDecoder_202508181516_perCell_data.csv")

non.procode.cells <- c( "no procode signal", "1 or 2 unexpected signals","3 or more unexpected signals")

perCellData$Procode_combination[perCellData$Procode_combination==""] <- "Untransduced"
perCellData$Id[perCellData$Id %in% non.procode.cells] <- "Other"

# add source and individual columns
# will need to be changed depending on how you named your samples
perCellData <- perCellData %>% 
  separate(sample, remove = TRUE, sep = "[ _]", into = c("X1", "tissue", "status", "X4", "mouse", "X6",
                                                         "X7", "X8" ) ) %>%
  select(-X1, -X4, -X6, -X7, -X8)

#removing untransduced cells becuase they aren't needed
perCellData <- perCellData %>% filter (Id != "Other")

#annotate the data
perCellData <- perCellData %>% 
  rename(group = status)

# remove guides that didnt work in the experiment, such as 0 transduction or too few cells 
perCellData <- perCellData %>%
  filter(!str_detect(Id, "Celsr2"))

any(str_detect(perCellData$Id, "Celsr2"))

#segment
Group1_perCells <- perCellData %>%
  filter(group==Group1)

Group2_perCells <- perCellData %>%
  filter(group==Group2)






# --- load the pct data 

# -- load Group1
# read in CD8 percent matrix files-- --- --- ---

pct_matrix <- read.csv ("CD8_FlowcodeDecoder_202508181516_pct_matrix.csv")


pct_matrix <- pct_matrix %>% 
  separate(sample, remove = TRUE, sep = "[ _]", into = c("X1", "tissue", "status", "X4", "mouse", "X6",
                                                         "X7", "X8" ) ) %>%
  select(-X1, -X4, -X6, -X7, -X8)

pct_matrix <- pct_matrix %>% select(-no.procode.signal, -X1.or.2.unexpected.signals, -X3.or.more.unexpected.signals ) 

# remove non-working guides here
pct_matrix <- pct_matrix %>% select(-Celsr2)

# calculate totals per row
text_cols <- c("tissue", "mouse","status")
pct_matrix$total <- rowSums(pct_matrix[ , setdiff(names(pct_matrix), text_cols)], na.rm = TRUE)

## remove any rows with total = 0 ###
pct_matrix <- pct_matrix %>% filter(total != 0)

#correct frequencies to be among transduced cells only
pct_matrix <- pct_matrix %>% mutate(across(where(is.numeric)) *100/ total)
anyNA(pct_matrix)

#annotate the data with the group
pct_matrix <- pct_matrix %>% 
  rename(group = status)

Group1_pct_matrix <- pct_matrix %>%
  filter(group == Group1) %>%
  select(tissue, mouse, group, everything())  # Reorder columns

#into the required format
Group1_pct_matrix <- Group1_pct_matrix %>% 
  pivot_longer(cols = -c(tissue, mouse, group, total), names_to = "target_gene", values_to = "freq" )

#group2
Group2_pct_matrix <- pct_matrix %>%
  filter(group == Group2) %>%
  select(tissue, mouse, group, everything())  # Reorder columns

#into the required format
Group2_pct_matrix <- Group2_pct_matrix %>% 
  pivot_longer(cols = -c(tissue, mouse, group, total), names_to = "target_gene", values_to = "freq" )

#bind the pcts for checking rank-based log2fc later
tissue_pct_matrix <- bind_rows(Group1_pct_matrix,Group2_pct_matrix)

####  end - Process Data section ####


#### Next section ####
# Note I've left it as pancreas, for actual COPD use just change it to Lung, its fine

#### all options should run this - normalization to spleen
# need to remove any samples with <20 cells

tissue_counts <- perCellData %>%
  group_by( tissue, mouse, Id, group ) %>%
  filter( Id != "Other") %>%
  summarize( count = n()) %>%
  mutate(pct = round(count/sum(count)*100, 2),  .groups = "drop") 

tissue_mins <- tissue_counts %>%
  group_by(tissue, mouse, group) %>%
  summarize(total = sum(count), .groups = "drop") %>%
  mutate(min_det = 100 / total)

# calculate half-min frequency to replace zeros in percent matrix

tissue_mins$half_min <- tissue_mins$min_det/2

#this prints info on what, if anything, it will remove later from pct_matrix
low_cell_samples <- tissue_mins %>% filter(total <= 20)

if (nrow(low_cell_samples) > 0) {
  cat("Samples removed due to low cell number:\n")
  apply(low_cell_samples, 1, function(row) {
    cat(row["mouse"], row["tissue"], "—", row["total"], "cells\n")
  })
}



# --- Now to normalize the freq of the tissue of interest to the spleen, both groups
#makes a function to do so
norm_to_spleen <- function(group_matrix, group_name) {
  spleen_pcts <- group_matrix %>%
    filter(tissue == "Spleen") %>%
    na.omit()
  
  Lung_pcts <- group_matrix %>%
    filter(tissue == "Lung") %>%
    na.omit()
  
  norm_pcts <- Lung_pcts %>%
    select(mouse, tissue, target_gene, freq_Lung = freq) %>%
    inner_join(
      spleen_pcts %>%
        select(mouse, target_gene, freq_spleen = freq),
      by = c("mouse", "target_gene")
    ) %>%
    mutate(group = group_name) 
  
  #but if the spleen freq=0, you get infinity as norm_freq. so replace 0 with half_min
  # Join half_min from tissue_mins to norm_pcts by mouse, group, and tissue
  # i think the min dets here are just from Lung. need spleen too
  norm_pcts <- norm_pcts %>%
    left_join(
      tissue_mins %>% select(mouse, group, tissue, half_min),
      by = c("mouse", "group", "tissue")
    ) %>%
    mutate(freq_Lung= ifelse(freq_Lung == 0, half_min, freq_Lung)
    )
  
  # Join spleen-based half_min values
  norm_pcts <- norm_pcts %>%
    left_join(
      tissue_mins %>%
        filter(tissue == "Spleen") %>%
        select(mouse, group, spl_half_min = half_min),
      by = c("mouse", "group")
    )
  
  norm_pcts <- norm_pcts %>%
    mutate(
      freq_spleen   = ifelse(freq_spleen == 0, spl_half_min, freq_spleen))
  
  
  norm_pcts <- norm_pcts %>%
    mutate(norm_freq = freq_Lung / freq_spleen)
  
  
  
  return(norm_pcts)
}

# Applies the function to Group1 and Group2
norm_Group1_pcts <- norm_to_spleen(Group1_pct_matrix, Group1)
norm_Group2_pcts <- norm_to_spleen(Group2_pct_matrix, Group2)

#and bind the two groups together
pct_matrix <- bind_rows(norm_Group1_pcts,norm_Group2_pcts)

#as the norm_freq is of the Lung, annotate it
pct_matrix<-pct_matrix %>%
  mutate(tissue = "Lung")

#remove NAs
pct_matrix <- pct_matrix %>%
  na.omit(data)







#add on the total cell count per sample column
pct_matrix_corr <- left_join(pct_matrix, tissue_mins, by = c("tissue", "mouse", "group"))
#and filter by it
pct_matrix_corr <- pct_matrix_corr %>% filter(total>20)
#tidy the dataframe
pct_matrix_corr <- pct_matrix_corr %>% select(-freq_Lung, -freq_spleen, -total)

pct_matrix_corr <- pct_matrix_corr %>% select(-half_min.x, -half_min.y, -spl_half_min, -min_det)
#### end of normalization to spleen section ####


#### parametric - between two separate groups of mice ####
# Run the comparison code, look at the treg vs tconv as exemplar
# Remember UNPAIRED

graphs.dir <- paste0("./", flowcode.seed, "_graphs/")

if( !file.exists(graphs.dir) ){
  dir.create(graphs.dir)
}

comparison.dir <- paste0(Group1, "_vs_" ,Group2 ,"/")
if( !file.exists( file.path(graphs.dir, comparison.dir) ) ){
  dir.create(file.path(graphs.dir, comparison.dir) )
}

# set up collection of p.vals and log2FC
comparison.results <- data.frame()

reference_subset <- Group2  # Specify your reference subset

### altered to not loop through tissues
#this needs modification , remove t 
# as compared to the parametric compare t cells loop I stopped it cycling through tissues
# and made it unpaired, so it doesnt need to organise mice to match
# nd made it work when there are different numbers of mic in each group
# by changing slightly the way it calculates log2fc
  

subset_data <- pct_matrix_corr %>% filter(group == Group1)
ref_data <- pct_matrix_corr %>% filter(group == Group2)

resdf = data.frame(target_gene = sort(unique(subset_data$target_gene)), log2FC = 0, pvalue = 1 )
# Initialize an empty results data frame
comparison.results <- data.frame()

# Loop over each gene in resdf
for (i in 1:nrow(resdf)) {
  g <- resdf$target_gene[i]
  
  # Filter data for the current gene
  subset_xvdf <- subset_data %>% filter(target_gene == g)
  ref_yvdf <- ref_data %>% filter(target_gene == g)
  
  xv <- subset_xvdf$norm_freq
  yv <- ref_yvdf$norm_freq
  
  if (length(c(xv, yv)) > 2) {
    # Compute log2 fold change
    log2fc <- log2(mean(xv + 0.001) / mean(yv + 0.001))
    
    # Perform t-test
    res <- t.test(log2(xv + 0.001), log2(yv + 0.001), paired = FALSE)
    pval <- res$p.value
    
    # Store results in a temporary data frame
    temp_result <- data.frame(
      target_gene = g,
      log2FC = log2fc,
      pvalue = pval
    )
    
    # Append to master results
    comparison.results <- rbind(comparison.results, temp_result)
    
    # Plotting
    plot_df <- bind_rows(
      data.frame(norm_freq = xv, group = Group1),
      data.frame(norm_freq = yv, group = reference_subset)
    )
    
    p <- ggplot(plot_df) +
      scale_x_discrete(name = "", labels = c(Group1, reference_subset)) +
      scale_y_log10() +
      scale_fill_manual(values = c("darkgrey", "firebrick1")) +
      ggtitle(paste(g, ":", Group1, "vs", reference_subset, 
                    "\n log2FC =", round(log2fc, 3), 
                    "; pvalue =", round(pval, 3))) +
      geom_boxplot(aes(x = group, y = norm_freq, fill = group), alpha = 0.2) +
      geom_point(aes(x = group, y = norm_freq), position = position_jitter(width = 0.2)) +
      theme_bw() +
      theme(legend.position = "none")
    
    ggsave(paste0(graphs.dir, comparison.dir, Group1, "_vs_", reference_subset, "_", g, ".pdf"), plot = p, height = 10, width = 10, units = "cm")
  }
}

# Apply multiple hypothesis correction
comparison.results$pvalue <- p.adjust(comparison.results$pvalue, method = "BH")

# Calculate distance metric
maxmlog10pvalue <- max(-log10(comparison.results$pvalue))
maxlog2FC <- max(abs(comparison.results$log2FC))

comparison.results$dist0 <- sqrt(comparison.results$log2FC^2 + 
                                   (-log10(comparison.results$pvalue) * maxlog2FC / maxmlog10pvalue)^2)

# Sort by distance
comparison.results <- comparison.results %>% arrange(-dist0)

# Save results
write.csv(comparison.results, paste0(graphs.dir, comparison.dir, Group1, "_vs_", reference_subset, "_pvals.csv"))



# #### cache in case the above section deosnt work####
# subset_data <- pct_matrix_corr %>% filter(subset == s)
# ref_data <- pct_matrix_corr %>% filter(subset == reference_subset)
# 
# resdf = data.frame(target_gene = sort(unique(subset_data$target_gene)), log2FC = 0, pvalue = 1 )
# 
# for (i in 1:nrow(resdf)) {
#   g = resdf$target_gene[i]
#   
#   subset_xvdf = subset_data %>% filter(target_gene == g) 
#   ref_yvdf = ref_data %>% filter(target_gene == g) 
#   
#   xv = subset_xvdf$norm_freq
#   yv = ref_yvdf$norm_freq
#   
#   if (length(c(xv,yv)) > 2){
#     resdf$log2FC[i] = mean(log2((xv + 0.001) / (yv + 0.001)))
#     
#     res = t.test(x = log2(xv+0.001),
#                  y = log2(yv+0.001),
#                  paired = FALSE,
#                  alternative = "two.sided")
#     resdf$pvalue[i] = res$p.value
#     
#     ggplot(data.frame(subset_xv = xv, ref_yv = yv, id = 1:length(xv)) %>% 
#              pivot_longer(cols = -id, names_to = "group", values_to = "norm_freq")) +
#       scale_x_discrete(name = "", labels = c(s, reference_subset)) +
#       scale_y_log10() +
#       scale_fill_manual(values = c("darkgrey", "firebrick1")) +
#       ggtitle(paste(g, ":", s, "vs", reference_subset, "\n log2FC =", round(resdf$log2FC[i], 3), "; pvalue =", round(res$p.value, 3))) +
#       geom_boxplot(aes(x = group, y = norm_freq, fill = group), alpha = 0.2) +
#       geom_point(aes(x = group, y = norm_freq)) +
#       geom_line(aes(x = group, y = norm_freq, group = id)) +
#       theme_bw() +
#       theme(legend.position = "none")
#     
#     ggsave(paste0(graphs.dir, comparison.dir, s, "_vs_", reference_subset, "_", g, ".pdf"), height = 10, width = 10, units = "cm")
#   }
# }
# 
# # Apply multiple hypothesis testing correction
# resdf$pvalue <- p.adjust(resdf$pvalue, method = "BH")
# 
# maxmlog10pvalue <- max(-log10(resdf$pvalue))
# maxlog2FC <- max(abs(resdf$log2FC))
# 
# resdf$dist0 <- sqrt(resdf$log2FC^2 + (-log10(resdf$pvalue) * maxlog2FC / maxmlog10pvalue)^2)
# resdf <- resdf %>% arrange(-dist0)
# 
# comparison.results <- rbind(comparison.results, resdf)
# write.csv(resdf, paste0(graphs.dir, comparison.dir, s, "_vs_", reference_subset, "_pvals.csv"))
# #### ####

#### Volcano + Heatmap plotting ####
# --- volcano 
resdf <- comparison.results

library(ggplot2)
library(ggrepel)
library(dplyr)

volcano <- ggplot(resdf) +
  ggtitle(paste(Group1, "vs", reference_subset, "Volcano")) +
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_continuous(limits = c(0, 4.5)) +
  geom_point(aes(x = log2FC, y = -log10(pvalue))) +
  geom_hline(yintercept = -log10(0.05), color = "firebrick1") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = -1, linetype = "dashed", color = "red") + 
  geom_text_repel(data = resdf %>% filter(pvalue < 0.05 & (log2FC <= -1 | log2FC >= 1)),
                  aes(x = log2FC, y = -log10(pvalue), label = target_gene),
                  nudge_x = ifelse(!is.na(resdf$log2FC) & resdf$log2FC < 0, -0.5, 0.5),
                  max.overlaps = 20) +
  theme(
    panel.background = element_rect(fill = "white", color = "black", size = 0.25),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_line(color = "#efefef"),
    panel.grid.minor = element_line(color = "#efefef"),
    axis.line = element_line(color = "black", size=0.25)
  )

print(volcano)


ggsave(
  filename = file.path(graphs.dir, paste0(Group1, "_vs_", Group2, "_Volcano_new.png")),
  plot = volcano,
  height = 10,
  width = 10,
  units = "cm"
)

### long heatmap code
#  Remove the column 'distr0' from comparison.results
comparison.results <- comparison.results %>% select(-dist0)

#export out the comparison results as a csv for use in other scripts
#change naming to be descriptive
write.csv(comparison.results, paste0(Group1, "_vs_", Group2, "_comparison_results.csv"), row.names = FALSE)

#taken from my script "log_abundance_volcano_and_heatmap"
#replacing the standard heatmap and volcano

library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)


heatmap.input <- read.csv(paste0(Group1, "_vs_", Group2, "_comparison_results.csv"))


#to format the result of the statistical test for a heatmap
# you can alter the p value threshold as you like


# heatmap specific code, to only show colours for significant p<0.01 values
heatmap.input[heatmap.input$pvalue>0.05,]$log2FC <- 0
#and the log2fc threshold - to match any log2fc threshold your volcano has, optional
heatmap.input[heatmap.input$log2FC > -1 & heatmap.input$log2FC < 1, ]$log2FC <- 0

heatmap.input <- heatmap.input %>% select(-pvalue)



#  Aggregating duplicate values using mean

heatmap.input <- heatmap.input %>%
  group_by(target_gene) %>%
  summarise(log2FC = mean(log2FC)) %>%
  ungroup()



# Calculate the average log2FC for each target_gene across all tissues
average_log2FC <- aggregate(log2FC ~ target_gene, data = heatmap.input, FUN = mean)

# Rank the target genes based on the average log2FC
# This will naturally order from most negative to most positive
rank_log2FC <- rank(average_log2FC$log2FC)

# Create a data frame with target_gene and rank_log2FC
ranked_genes <- data.frame(target_gene = average_log2FC$target_gene, 
                           rank_log2FC = rank_log2FC,
                           avg_log2FC = average_log2FC$log2FC)

heatmap.filtered <- heatmap.input %>%
  filter(target_gene %in% ranked_genes$target_gene)

heatmap.ordered <- heatmap.filtered %>%
  arrange(match(target_gene, ranked_genes$target_gene))


### option 1 - genes shown in order of rank_log2FC ###
# This will order from most negative to most positive log2FC
ranked_genes <- average_log2FC %>%
  arrange(log2FC) %>%
  mutate(rank_log2FC = row_number())


### option 2 - genes shown in alphabetical order ###
# better for when you have lots of heatmaps and want to compare results for a guide between them
# ranked_genes <- average_log2FC %>%
#   arrange(target_gene) %>%
#   mutate(rank_log2FC = row_number())

# all options rejoin here - reshape data into wide format: one row per gene, one column per tissue
heatmap.df <- as.data.frame(heatmap.ordered)
rownames(heatmap.df) <- heatmap.df$target_gene
heatmap.df <- heatmap.df %>% select(-target_gene)
heatmap.matrix <- as.matrix(heatmap.df)



# If you want to name the single column (e.g., "log2FC"), do this:
#colnames(heatmap.matrix) <- "log2FC"

# Define color palette and breaks
palette.length <- 50
heatmap.colors <- colorRampPalette(c("blue", "white", "red"))(palette.length)

heatmap.breaks <- c(
  seq(min(heatmap.matrix, na.rm = TRUE), 0, length.out = ceiling(palette.length / 2) + 1),
  seq(max(heatmap.matrix, na.rm = TRUE) / palette.length, max(heatmap.matrix, na.rm = TRUE), length.out = floor(palette.length / 2))
)

# Generate heatmap
# Generate heatmap object (not printed yet)


gene.impact.heatmap <- pheatmap(
  t(heatmap.matrix),
  color = heatmap.colors,
  breaks = heatmap.breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey",
  show_rownames = TRUE,
  show_colnames = TRUE,
  angle_col = 45,
  fontsize = 12,
  fontsize_row = 10,
  fontsize_col = 10,
  silent = TRUE  # Prevents auto-plotting
)

# Create title grob
heatmap.title <- textGrob(
  paste0(Group1, " vs ", reference_subset, " Heatmap (p=0.05)."),
  gp = gpar(fontsize = 14, fontface = "bold")
)

# Save with extra space for rotated labels
png(
  filename = paste0(graphs.dir, comparison.dir, Group1, "_vs_", reference_subset, "_heatmap.png"),
  height = 5,
  width = 24,
  units = "cm",
  res = 300
)

# Arrange title and heatmap with spacing
grid.arrange(
  heatmap.title,
  gene.impact.heatmap$gtable,
  ncol = 1,
  heights = c(1, 10)  # More space for heatmap
)

dev.off()



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




#### end of Volcano+Heatmap plotting ####







#### Non parametric code Option ####

#for adding the full functionality of allowing this code to also run the parametric tests
#I need to do tissue mins as normalized mins
# take the code from the old parametric compare t cell types which is in your validation folder
#later

# --- Now the non-parametric part

#I am assuming you have previously, in the masterscript, determined that non-parametric is more suitable

# Based on the analysis above, determine if you want to ranks and non-parametric test for Pancreas
# set this to TRUE or FALSE
# if you set it to FALSE then statistical tests will run on comparing percentage frequencies of the guides,
# if you set it to TRUE then the statistical tests will compare the guides non-parametrically, by ranks
USE_RANKS <- TRUE

# Translate ranks back to spleen values, i.e. assuming there is 
# the same distribution/TCR clonality in tissue and spleen - assuming this 
# (and thus removing the signal of excessive TCR expansion), 
# we want to measure the effect size in terms of logFC as in other comparisons
# it is an assumption, but should be reasonably accurate
# if you set it to FALSE, the code will output graphs with the x axis of change in ranks, not logfc
RANK_BASED_LOGFC <- TRUE

if ( USE_RANKS ) {
  pct_matrix_rank <- pct_matrix_corr %>% 
    group_by( group, mouse ) %>% 
    mutate( rank_raw = rank( norm_freq, ties.method = "average" ) ) %>% 
    mutate( mean_rank_raw = mean( rank_raw ) ) %>% 
    ungroup() %>%  
    
    # Normalize to deal with unequal number of mice per gene
    mutate( max_mean_rank_raw = max( mean_rank_raw ) ) %>% 
    mutate( rank = rank_raw * ( max_mean_rank_raw / mean_rank_raw ) )
}

#this is only an option if use_ranks is set to TRUE
# to check the rank-based log2fc 
if ( USE_RANKS ) {
  t_pct_matrix_rank <- tissue_pct_matrix %>% 
    group_by( tissue, mouse ) %>% 
    mutate( rank_raw = rank( freq, ties.method = "average" ) ) %>% 
    mutate( mean_rank_raw = mean( rank_raw ) ) %>% 
    ungroup() %>%  
    
    # Normalize to deal with unequal number of mice per gene
    mutate( max_mean_rank_raw = max( mean_rank_raw ) ) %>% 
    mutate( rank = rank_raw * ( max_mean_rank_raw / mean_rank_raw ) )
}




if ( CHECK_DATA_REVERT_TO_SPLEEN_DISTRIB <- TRUE ) {
  t_pct_matrix_rank %>% 
    group_by( tissue, target_gene ) %>% 
    mutate( med_freq = median( freq ),
            med_rank = median( rank ) ) %>%
    ungroup() %>% 
    ggplot( aes( x = med_rank, y = med_freq ) ) +
    geom_point() +
    geom_smooth( method = loess ) +
    facet_wrap( vars( tissue ) )
  
  #I think the N/As due to multiple sets, in freq, cause issue?
  #so I remove these rows, and the graph looks more as expected
  d_spleen_ranks <- t_pct_matrix_rank %>% 
    filter(tissue == "Spleen") %>%
    filter(!is.na(freq)) %>%
    group_by(tissue, target_gene) %>% 
    mutate( med_freq = median( freq ),
            med_rank = median( rank ) ) %>%
    ungroup()
  
  
  d_pancreas_ranks <- t_pct_matrix_rank %>% 
    filter( tissue == "Pancreas" ) %>%
    filter(!is.na(freq)) %>%
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

# --- now the statistical test can be run between the groups

## Simple Group1 vs Group2 freq comparison-- -- --- -- --- - ---
graphs.dir <- paste0("./", flowcode.seed, "_graphs/")

dir.create(graphs.dir, recursive = TRUE, showWarnings = TRUE)


if( !file.exists(graphs.dir) ){
  dir.create(graphs.dir)
}

comparison.dir <- paste0(Group1, "_vs_", Group2)
if( !file.exists( file.path(graphs.dir, comparison.dir) ) ){
  dir.create(file.path(graphs.dir, comparison.dir) )
}

reference_group <- Group2


# this performs the chosen statistical test and outputs boxplots and volcanos,
# it depends on if you have set USE_RANKS to TRUE or FALSE
# within the volcano code you can adjust p value significance level and log2fc thresholds

# VACLAV: The code for non-parametric approach is copied from below and modified to
# allow for USE_RANKS by using "flag"
# STILL FROM_VACLAV

#For correct graph title, set Tissue to your tissue of interest
# the tissue this data is for
Tissue <- "Pancreas"


if ( PLOT_BOXPLOTS_AND_VOLCANO_FOR_RANKS <- TRUE ) {
  
  
  comparison.results <- data.frame()
  
  if ( USE_RANKS ) { pct_matrix_corr <- pct_matrix_rank } 
  for (t in unique(pct_matrix_corr$group)) {
    
    if(t != reference_group){
      group.dir <- paste0(graphs.dir, comparison.dir, "/")
      
      if( !file.exists(group.dir)){
        dir.create(group.dir)
      }
      
      goi = pct_matrix_corr %>% filter(group == t)
      
      ref_goi = pct_matrix_corr %>% filter(group == reference_group)
      
      resdf = data.frame(target_gene = sort(unique(goi$target_gene)), 
                         log2FC = 0, pvalue_t = 1, pvalue = 1,
                         median_rank_diff = 0, pvalue_w = 1,
                         rank_based_logfc = 0 )
      
      for (i in 1:nrow(resdf)) {
        g = resdf$target_gene[i]
        xvdf = (goi %>% filter(target_gene == g) %>% arrange(mouse))
        yvdf = ref_goi %>% filter(target_gene == g)%>% filter(mouse %in% xvdf$mouse) %>% arrange(mouse)
        
        xv = (xvdf %>% filter(mouse %in% yvdf$mouse))$norm_freq
        yv = yvdf$norm_freq
        
        if (length(c(xv,yv)) >2){
          
          # Parametric testing
          resdf$log2FC[i] = mean(log2((xv + 0.001) / (yv + 0.001)))
          res_t <- t.test(x = log2(xv + 0.001), y = log2(yv + 0.001),
                          paired = TRUE, alternative = "two.sided")
          resdf$pvalue_t[i] = res_t$p.value
          
          ggplot(data.frame(goi = xv, ref = yv, id = 1:length(xv)) %>% 
                   pivot_longer(cols = -id, names_to = "group", values_to = "norm_freq")) +
            scale_x_discrete(name = "", labels = c(reference_group, t)) +
            scale_y_log10() +
            scale_fill_manual(values = c("darkgrey", "firebrick1")) +
            ggtitle(paste0(g, ": ", t, " vs ", reference_group, 
                           "\n log2FC = ", round(resdf$log2FC[i], 3), 
                           "; t-test p = ", round(res_t$p.value, 3))) +
            geom_boxplot(aes(x = group, y = norm_freq, fill = group), alpha = 0.2) +
            geom_point(aes(x = group, y = norm_freq)) +
            geom_line(aes(x = group, y = norm_freq, group = id)) +
            theme_bw() +
            theme(legend.position = "none")
          
          ggsave(paste0(group.dir, t, "_vs_", reference_group, "_", g, "_prop.pdf"), 
                 height = 10, width = 10, units = "cm")
          
          # Nonparametric testing
          if ( USE_RANKS ) {
            xv_r = (xvdf %>% filter(mouse %in% yvdf$mouse))$rank
            yv_r = yvdf$rank
            resdf$median_rank_diff[i] = median(xv_r - yv_r)
            res_w <- wilcox.test(x = xv_r + 0.001, y = yv_r + 0.001, 
                                 paired = TRUE, alternative = "two.sided")
            resdf$pvalue_w[i] = res_w$p.value
            
            ggplot(data.frame(goi = xv_r, spleen = yv_r, id = 1:length(xv_r)) %>% 
                     pivot_longer(cols = -id, names_to = "group", values_to = "rank")) +
              scale_x_discrete(name = "", labels = c(reference_group, t)) +
              # scale_y_log10() +
              scale_fill_manual(values = c("darkgrey", "firebrick1")) +
              ggtitle(paste0(g, ": ", t, " vs ", reference_group, 
                             "\n Rank diff = ", round(resdf$median_rank_diff[i], 3), 
                             "; Wilcox p = ", round(res_w$p.value, 3))) +
              geom_boxplot(aes(x = group, y = rank, fill = group), alpha = 0.2) +
              geom_point(aes(x = group, y = rank)) +
              geom_line(aes(x = group, y = rank, group = id)) +
              theme_bw() +
              theme(legend.position = "none")
            
            ggsave(paste0(group.dir, t, "_vs_", reference_group, "_", g, "_rank.pdf"), 
                   height = 10, width = 10, units = "cm")
            
          }
          if ( USE_RANKS & RANK_BASED_LOGFC ) {
            resdf$rank_based_logfc[i] = log2( 
              rank_to_value( median( xv_r ) ) / rank_to_value( median( yv_r ) ) )
            
          }
        }
      }
      
      
      
      resdf$group <- rep(t)
      
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
      write.csv(resdf, paste0(group.dir, t, "_vs_", reference_group, 
                              ifelse( USE_RANKS, "_ranks", "" ), ".csv"))
      
    }
    
  }
  

  
  if ( VG_VOLCANO_ON_RANKS_NORM <- TRUE ) {
    
    group_of_interest <- Group1
    reference_group <- Group2
  
    subset.results.vs.res <- read.csv(
      file.path(
        paste0(flowcode.seed, "_graphs"),
        paste0(Group1, "_vs_", Group2),
        paste0(Group1, "_vs_", reference_group, ifelse(USE_RANKS, "_ranks", ""), ".csv")
      ))
      
      
      
      # option 1 - normal plotting of everything with the log2fc threshold
    # or well it would be but as I haven't put in the parametric code fully, above, this
    # is actually just the rank difference.
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
      
      plot.figure <- paste0(group_of_interest, "_vs_", reference_group, Tissue,"_Volcano p=0.01_rank_diff.png" )
      
      resdf <- subset.results 
      resdf <- resdf %>% arrange(-dist0)
      
      # Store the ggplot object in a variable
      p <- ggplot(resdf) +
        ggtitle(paste(group_of_interest, "vs", reference_group, Tissue,"_Volcano p=0.01")) +
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
        plot.figure_2 <- paste0(Group1, "_vs_", Group2,"_", Tissue,
                                "_Volcano p=0.05_rank_logfc_labels.png" )
        p_rank_logfc <- ggplot(resdf) +
          ggtitle(paste(Group1, "vs", Group2, "_Volcano p=0.05")) +
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



#### End of the non parametric code ####













