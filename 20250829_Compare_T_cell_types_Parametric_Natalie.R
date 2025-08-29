#This is NOT the main analysis code. This is for comparing T cell types when parametric testing is suitable

# Comparing between T cell types
# uses pancreas guide frequencies normalized to spleen, to compare the important genes in T cell types to home to the pancreas
# you can easily change instances of pancreas to whatever your tissue of interest is
# plots volcanos of this
# and has the bar chart code for comparing T cell type counts and %s between tissues 

# this exact version was used for the validation study and runs the data as it is, no clone removal in cd8s, and is parametric
# is a masterscript that does all the between t cells comparisons, so you don't need 3 separate scripts
# not very efficient code as it is just cd8 vs tconv, then the treg vs tconv, then the cd8 vs treg sections. almost identical
# it assumes normally distributed data that is suitable for a paired t test

#last edited 29th August Natalie Ng - to make presentable

required.packages <- c("Rtsne", "ggplot2", "RColorBrewer", "dplyr", "emmeans",  
                       "EmbedSOM", "tidyr", "coda", "scattermore",
                       "data.table", "ggrepel", "pheatmap", "parallelly",
                       "purrr", "lsa", "stringr" )

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

flowcode.seed <- 20250709
flowcode.threads <- parallelly::availableCores() - 1
data.table::setDTthreads(flowcode.threads)

#### -- CD8 vs Tconv --- ####
#treg data is imported later, this is inefficient code but it works

# import CD8 perCell data --- --- --- --
perCellData <- fread(file = list.files(pattern = "CD8_corrected_FlowcodeDecoder_20250211_perCell_data.csv"))

non.procode.cells <- c( "no procode signal", "1 or 2 unexpected signals","3 or more unexpected signals")

perCellData$Procode_combination[perCellData$Procode_combination==""] <- "Untransduced"
perCellData$Id[perCellData$Id %in% non.procode.cells] <- "Other"

# add source and individual columns
perCellData <- perCellData %>% 
  separate(sample, remove = TRUE, sep = "[ _]", into = c("X1", "tissue", "X3", "X4", "X5", "X6",
                                                         "mouse", "X8", "X9" ) ) %>%
  select(-X1, -X3, -X4, -X5, -X6, -X8, -X9)

perCellData <- perCellData %>%
  mutate(subset = "CD8")


# optional: remove non-working guides here ###
library(stringr)

#Ids you want to remove
perCellData <- perCellData %>%
  filter(!str_detect(Id, "Hepacam|Celsr2"))

any(str_detect(perCellData$Id, "Hepacam|Celsr2"))

CD8_perCellData <- perCellData



# import Tconv perCell data--- --- --- --
perCellData <- fread(file = list.files(pattern = "Tconv.2_FlowcodeDecoder_20250211_perCell_data.csv"))

non.procode.cells <- c( "no procode signal", "1 or 2 unexpected signals","3 or more unexpected signals")

perCellData$Procode_combination[perCellData$Procode_combination==""] <- "Untransduced"
perCellData$Id[perCellData$Id %in% non.procode.cells] <- "Other"

# add source and individual columns
perCellData <- perCellData %>% 
  separate(sample, remove = TRUE, sep = "[ _]", into = c("X1", "tissue", "X3", "X4", "X5", "X6",
                                                         "mouse", "X8", "X9" ) ) %>%
  select(-X1, -X3, -X4, -X5, -X6, -X8, -X9)

perCellData <- perCellData %>%
  filter(!str_detect(Id, "Hepacam|Celsr2"))

any(str_detect(perCellData$Id, "Hepacam|Celsr2"))

perCellData <- perCellData %>%
  mutate(subset = "tconv")

tconv_perCellData <- perCellData

# Binding dataframes
perCellData <- bind_rows(tconv_perCellData, CD8_perCellData)

any(is.na(perCellData))

#removing untransduced cells becuase they aren't needed
perCellData <- perCellData %>% filter (Id != "Other")

# read in CD8 percent matrix file--- --- --- --

pct_matrix <- read.csv(file = list.files(pattern = "CD8_corrected_FlowcodeDecoder_20250211_pct_matrix.csv"))

## separate out sample names: this will need to be modified depending on your naming convention
pct_matrix <- pct_matrix %>% 
  separate(sample, remove = TRUE, sep = "[ _]", into = c("X1", "tissue", "X3", "X4", "X5",
                                                         "X6", "mouse", "X8", "X9", "X10", "X11" ) ) %>%
  select(-X1, -X3, -X4, -X5, -X6, -X8, -X9, -X10, -X11)

pct_matrix <- pct_matrix %>% select(-no.procode.signal, -X1.or.2.unexpected.signals, -X3.or.more.unexpected.signals ) 

# remove non-working guides here (change as needed)
pct_matrix <- pct_matrix %>% select( -Hepacam, -Celsr2)

# correct frequencies to be among transduced cells only
pct_matrix$total <- rowSums(pct_matrix[,-c(1,2)])

pct_matrix <- pct_matrix %>%
  mutate(subset = "CD8") 


## Important: remove any rows with total = 0 ###
pct_matrix <- pct_matrix %>% filter(total != 0)

pct_matrix <- pct_matrix %>% mutate(across(where(is.numeric)) *100/ total)
anyNA(pct_matrix)

pct_matrix <- pct_matrix %>% 
  pivot_longer(cols = -c(tissue, mouse, subset, total), names_to = "target_gene", values_to = "freq" )

CD8_pct_matrix <- pct_matrix

# read in tconv percent matrix file--- --- --- --

pct_matrix <- read.csv("Tconv.2_FlowcodeDecoder_20250211_pct_matrix.csv")

## separate out sample names: this will need to be modified depending on your naming convention
pct_matrix <- pct_matrix %>% 
  separate(sample, remove = TRUE, sep = "[ _]", into = c("X1", "tissue", "X3", "X4", "X5",
                                                         "X6", "mouse", "X8", "X9", "X10", "X11" ) ) %>%
  select(-X1, -X3, -X4, -X5, -X6, -X8, -X9, -X10, -X11)

pct_matrix <- pct_matrix %>% select(-no.procode.signal, -X1.or.2.unexpected.signals, -X3.or.more.unexpected.signals ) 

# remove non-working guides here (change as needed)
pct_matrix <- pct_matrix %>% select( -Hepacam, -Celsr2)

# correct frequencies to be among transduced cells only
pct_matrix$total <- rowSums(pct_matrix[,-c(1,2)])

pct_matrix <- pct_matrix %>%
  mutate(subset = "tconv") 

## Important: remove any rows with total = 0 ###
pct_matrix <- pct_matrix %>% filter(total != 0)

pct_matrix <- pct_matrix %>% mutate(across(where(is.numeric)) *100/ total)
anyNA(pct_matrix)

pct_matrix <- pct_matrix %>% 
  pivot_longer(cols = -c(tissue, mouse, total, subset), names_to = "target_gene", values_to = "freq" )

tconv_pct_matrix <- pct_matrix

### edited to include the normalization to spleen


CD8_spleen_pcts <- CD8_pct_matrix %>% filter (tissue == "Spleen")
CD8_pancreas_pcts <- CD8_pct_matrix %>% filter (tissue == "Pancreas")

CD8_spleen_pcts <- CD8_spleen_pcts %>%
  na.omit(data)

CD8_pancreas_pcts <- CD8_pancreas_pcts %>%
  na.omit(data)



# to normalise CD8s Merge pancreas and spleen data on mouse and target_gene
  norm_CD8_pcts <- CD8_pancreas_pcts %>%
  select(mouse, target_gene, freq_pancreas = freq) %>%
  inner_join(
    CD8_spleen_pcts %>%
      select(mouse, target_gene, freq_spleen = freq),
    by = c("mouse", "target_gene")
  )

norm_CD8_pcts <- norm_CD8_pcts %>%
  mutate(subset = "CD8") %>%
  mutate(norm_freq = freq_pancreas / freq_spleen)  # Compute normalized frequency
  
#normalising tconvs
tconv_spleen_pcts <- tconv_pct_matrix %>% filter (tissue == "Spleen")
tconv_pancreas_pcts <- tconv_pct_matrix %>% filter (tissue == "Pancreas")

#this could be causing mismatches, differeent numbers of values per mouse. try without. na's get removed later anyway
tconv_spleen_pcts <- tconv_spleen_pcts %>%
  na.omit(data)

tconv_pancreas_pcts <- tconv_pancreas_pcts %>%
  na.omit(data)



# to normalise tconvs Merge pancreas and spleen data on mouse and target_gene
norm_tconv_pcts <- tconv_pancreas_pcts %>%
  select(mouse, target_gene, freq_pancreas = freq) %>%
  inner_join(
    tconv_spleen_pcts %>%
      select(mouse, target_gene, freq_spleen = freq),
    by = c("mouse", "target_gene")
  )

norm_tconv_pcts <- norm_tconv_pcts %>%
  mutate(subset = "tconv")  %>%
  mutate(norm_freq = freq_pancreas / freq_spleen)  # Compute normalized frequency


#-- -- ---- --
#
#
#back to original code
# Binding dataframes
pct_matrix <- bind_rows(norm_tconv_pcts, norm_CD8_pcts)
pct_matrix<-pct_matrix %>%
  mutate(tissue = "Pancreas")

#remove NAs
pct_matrix <- pct_matrix %>%
  na.omit(data)

# calculate minimum detection frequency per sample from perCellData
#  normalise the pancreas min detection freq to the spleen for scale

tissue_counts <- perCellData %>%
  group_by( tissue, mouse, Id, subset ) %>%
  filter( Id != "Other") %>%
  summarize( count = n()) %>%
  mutate(pct = round(count/sum(count)*100, 2)) 

tissue_mins <- tissue_counts %>%
  group_by( tissue, mouse, subset ) %>%
  summarize( total = sum(count) ) %>%
  mutate( min_det = 100/total ) 

pancreas_mins <- tissue_mins %>% filter(tissue == "Pancreas")
spleen_mins <- tissue_mins %>% filter(tissue == "Spleen")

#joining the pancreas and spleen data
normalized_mins <- pancreas_mins %>%
  inner_join(spleen_mins, by = c("mouse", "subset"), suffix = c("_pancreas", "_spleen"))
 
#normalizing the pcts
normalized_mins <- normalized_mins %>%
  mutate(normalized_min_det = min_det_pancreas / min_det_spleen,
         normalized_half_min = normalized_min_det / 2)

normalized_mins <- normalized_mins %>%
  rename(tissue=tissue_pancreas)

normalized_mins <- normalized_mins %>% unite(col =  tissue, mouse, subset)

# use half-min frequency to replace zeros in percent matrix

#tissue_mins <- tissue_mins %>% filter(tissue=="Pancreas")

#tissue_mins$half_min <- tissue_mins$min_det/2
#tissue_mins <- tissue_mins %>% unite(col = sample, tissue, mouse, subset)

pct_matrix_corr <- pct_matrix
pct_matrix_corr <- pct_matrix_corr %>% unite(col =  tissue, mouse, subset)


for (s in unique(pct_matrix_corr$tissue)) {
  min_temp <- filter(normalized_mins, tissue == s)
  pct_temp <- filter(pct_matrix_corr, tissue == s)
  
  if (nrow(min_temp) == 1) {
    pct_temp$norm_freq[pct_temp$norm_freq == 0] <- min_temp$normalized_half_min
    pct_matrix_corr$norm_freq[pct_matrix_corr$tissue %in% unique(pct_temp$tissue)] <- pct_temp$norm_freq
  }
}

# Cell number cutoff: filter out sample with fewer than 20 cells

pct_matrix_corr <- left_join(pct_matrix_corr, normalized_mins, by = "tissue")

pct_matrix_corr <- pct_matrix_corr %>%
  filter(total_pancreas >= 20, total_spleen >= 20)


pct_matrix_corr <- pct_matrix_corr %>% separate(tissue, remove = TRUE, sep = "_", into = c( "mouse", "subset"))  %>%
  select(-total_pancreas, -total_spleen, -min_det_pancreas, -min_det_spleen, -normalized_half_min, -tissue_spleen, -normalized_min_det, -freq_pancreas, -freq_spleen)

anyNA(pct_matrix_corr)

#
pct_matrix_corr <- pct_matrix_corr %>% mutate(tissue="Pancreas")


## Simple CD8 vs Tconv freq comparison--- --- --- --- ---
# to which I have added the log and BH as I have to the base/other processing scripts
graphs.dir <- paste0("./", flowcode.seed, "_graphs/")

if( !file.exists(graphs.dir) ){
  dir.create(graphs.dir)
}

comparison.dir <- "CD8_vs_tconv/"
if( !file.exists( file.path(graphs.dir, comparison.dir) ) ){
  dir.create(file.path(graphs.dir, comparison.dir) )
}

# set up collection of p.vals and log2FC
comparison.results <- data.frame()

reference_subset <- "tconv"  # Specify your reference subset

### altered to not loop through tissues
for (t in unique(pct_matrix_corr$tissue)) {
  
  tissue_data <- pct_matrix_corr %>% filter(tissue == t)
  
  for (s in unique(tissue_data$subset)) {
    if (s != reference_subset) {
      subset_data <- tissue_data %>% filter(subset == s)
      ref_data <- tissue_data %>% filter(subset == reference_subset)
      
      resdf = data.frame(target_gene = sort(unique(subset_data$target_gene)), log2FC = 0, pvalue = 1 )
      
      for (i in 1:nrow(resdf)) {
        g = resdf$target_gene[i]
        
        subset_xvdf = subset_data %>% filter(target_gene == g) %>% arrange(mouse)
        ref_yvdf = ref_data %>% filter(target_gene == g) %>% filter(mouse %in% subset_xvdf$mouse) %>% arrange(mouse)
        
        xv = (subset_xvdf %>% filter(mouse %in% ref_yvdf$mouse))$norm_freq
        yv = ref_yvdf$norm_freq
        
        if (length(c(xv,yv)) > 2){
          resdf$log2FC[i] = mean(log2((xv + 0.001) / (yv + 0.001)))
          
          res = t.test(x = log2(xv+0.001),
                       y = log2(yv+0.001),
                       paired = TRUE,
                       alternative = "two.sided")
          resdf$pvalue[i] = res$p.value
          
          ggplot(data.frame(subset_xv = xv, ref_yv = yv, id = 1:length(xv)) %>% 
                   pivot_longer(cols = -id, names_to = "group", values_to = "norm_freq")) +
            scale_x_discrete(name = "", labels = c(s, reference_subset)) +
            scale_y_log10() +
            scale_fill_manual(values = c("darkgrey", "firebrick1")) +
            ggtitle(paste(g, ":", s, "vs", reference_subset, "in", t, "\n log2FC =", round(resdf$log2FC[i], 3), "; pvalue =", round(res$p.value, 3))) +
            geom_boxplot(aes(x = group, y = norm_freq, fill = group), alpha = 0.2) +
            geom_point(aes(x = group, y = norm_freq)) +
            geom_line(aes(x = group, y = norm_freq, group = id)) +
            theme_bw() +
            theme(legend.position = "none")
          
          ggsave(paste0(graphs.dir, comparison.dir, t, "_", s, "_vs_", reference_subset, "_", g, ".pdf"), height = 10, width = 10, units = "cm")
        }
      }
      
      resdf$tissue <- rep(t)
      resdf$subset <- rep(s)
      
      # Apply multiple hypothesis testing correction
      resdf$pvalue <- p.adjust(resdf$pvalue, method = "BH")
      
      maxmlog10pvalue <- max(-log10(resdf$pvalue))
      maxlog2FC <- max(abs(resdf$log2FC))
      
      resdf$dist0 <- sqrt(resdf$log2FC^2 + (-log10(resdf$pvalue) * maxlog2FC / maxmlog10pvalue)^2)
      resdf <- resdf %>% arrange(-dist0)
      
      comparison.results <- rbind(comparison.results, resdf)
      write.csv(resdf, paste0(graphs.dir, comparison.dir, t, "_", s, "_vs_", reference_subset, "_pvals.csv"))
    }
  }
}


# volcano plots -- - ---- --

# Set the reference subset
reference_subset <- "tconv"
#definitely consider changing x and y axis limits to suit your data

#remove unwanted guides if needed
genes_to_remove <- c("Itgb2l", "Cd22", "Ctnna2", "Itga10")

# Filter the dataframe
resdf <- resdf[!resdf$target_gene %in% genes_to_remove, ]

# --- volcano option 1

library(ggplot2)
library(ggrepel)
library(dplyr)

volcano <- ggplot(resdf) +
  ggtitle(paste(s, "vs", reference_subset, "in", t, "volcano")) +
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


ggsave(file.path(graphs.dir, comparison.dir, plot.figure), plot = volcano, height = 10, width = 10, units = "cm")

#### end ####





#### -- Treg vs Tconv ---####

# import treg perCell data --- --- --- --
perCellData <- fread(file = list.files(pattern = "Treg_VSVg_FlowcodeDecoder_20250306_perCell_data.csv"))

non.procode.cells <- c( "no procode signal", "1 or 2 unexpected signals","3 or more unexpected signals")

perCellData$Procode_combination[perCellData$Procode_combination==""] <- "Untransduced"
perCellData$Id[perCellData$Id %in% non.procode.cells] <- "Other"

# add source and individual columns
perCellData <- perCellData %>% 
  separate(sample, remove = TRUE, sep = "[ _]", into = c("X1", "tissue", "X3", "X4", "X5", "X6",
                                                         "mouse", "X8", "X9" ) ) %>%
  select(-X1, -X3, -X4, -X5, -X6, -X8, -X9)

perCellData <- perCellData %>%
  mutate(subset = "treg")


# optional: remove non-working guides here ###
library(stringr)

#Ids you want to remove
perCellData <- perCellData %>%
  filter(!str_detect(Id, "Hepacam|Celsr2"))

any(str_detect(perCellData$Id, "Hepacam|Celsr2"))

treg_perCellData <- perCellData

#for this comparison you need these
perCellData <- rbind(treg_perCellData, tconv_perCellData)

# read in treg percent matrix file--- --- --- --

pct_matrix <- read.csv(file = list.files(pattern = "Treg_VSVg_FlowcodeDecoder_20250306_pct_matrix.csv"))

## separate out sample names: this will need to be modified depending on your naming convention
pct_matrix <- pct_matrix %>% 
  separate(sample, remove = TRUE, sep = "[ _]", into = c("X1", "tissue", "X3", "X4", "X5",
                                                         "X6", "mouse", "X8", "X9", "X10", "X11" ) ) %>%
  select(-X1, -X3, -X4, -X5, -X6, -X8, -X9, -X10, -X11)

pct_matrix <- pct_matrix %>% select(-no.procode.signal, -X1.or.2.unexpected.signals, -X3.or.more.unexpected.signals ) 

# remove non-working guides here (change as needed)
pct_matrix <- pct_matrix %>% select( -Hepacam, -Celsr2)

# correct frequencies to be among transduced cells only
pct_matrix$total <- rowSums(pct_matrix[,-c(1,2)])

pct_matrix <- pct_matrix %>%
  mutate(subset = "treg") 


## Important: remove any rows with total = 0 ###
pct_matrix <- pct_matrix %>% filter(total != 0)

pct_matrix <- pct_matrix %>% mutate(across(where(is.numeric)) *100/ total)
anyNA(pct_matrix)

pct_matrix <- pct_matrix %>% 
  pivot_longer(cols = -c(tissue, mouse, subset, total), names_to = "target_gene", values_to = "freq" )

treg_pct_matrix <- pct_matrix

### edited to include the normalization to spleen


#normalising tregs
treg_spleen_pcts <- treg_pct_matrix %>% filter (tissue == "Spleen")
treg_pancreas_pcts <- treg_pct_matrix %>% filter (tissue == "Pancreas")

#this could be causing mismatches, differeent numbers of values per mouse. try without. na's get removed later anyway
treg_spleen_pcts <- treg_spleen_pcts %>%
  na.omit(data)

treg_pancreas_pcts <- treg_pancreas_pcts %>%
  na.omit(data)



# to normalise tregs Merge pancreas and spleen data on mouse and target_gene
norm_treg_pcts <- treg_pancreas_pcts %>%
  select(mouse, target_gene, freq_pancreas = freq) %>%
  inner_join(
    treg_spleen_pcts %>%
      select(mouse, target_gene, freq_spleen = freq),
    by = c("mouse", "target_gene")
  )

norm_treg_pcts <- norm_treg_pcts %>%
  mutate(subset = "treg")  %>%
  mutate(norm_freq = freq_pancreas / freq_spleen)  # Compute normalized frequency


#-- -- ---- --
#
#
#back to original code
# Binding dataframes
pct_matrix <- bind_rows(norm_treg_pcts, norm_tconv_pcts)
pct_matrix<-pct_matrix %>%
  mutate(tissue = "Pancreas")

#remove NAs
pct_matrix <- pct_matrix %>%
  na.omit(data)

# calculate minimum detection frequency per sample from perCellData
# statistically I should likewise normalise the min detection freq to the spleen for scale

tissue_counts <- perCellData %>%
  group_by( tissue, mouse, Id, subset ) %>%
  filter( Id != "Other") %>%
  summarize( count = n()) %>%
  mutate(pct = round(count/sum(count)*100, 2)) 

tissue_mins <- tissue_counts %>%
  group_by( tissue, mouse, subset ) %>%
  summarize( total = sum(count) ) %>%
  mutate( min_det = 100/total ) 

pancreas_mins <- tissue_mins %>% filter(tissue == "Pancreas")
spleen_mins <- tissue_mins %>% filter(tissue == "Spleen")

#joining the pancreas and spleen data
normalized_mins <- pancreas_mins %>%
  inner_join(spleen_mins, by = c("mouse", "subset"), suffix = c("_pancreas", "_spleen"))

#normalizing the pcts
normalized_mins <- normalized_mins %>%
  mutate(normalized_min_det = min_det_pancreas / min_det_spleen,
         normalized_half_min = normalized_min_det / 2)

normalized_mins <- normalized_mins %>%
  rename(tissue=tissue_pancreas)

normalized_mins <- normalized_mins %>% unite(col =  tissue, mouse, subset)

# use half-min frequency to replace zeros in percent matrix

#tissue_mins <- tissue_mins %>% filter(tissue=="Pancreas")

#tissue_mins$half_min <- tissue_mins$min_det/2
#tissue_mins <- tissue_mins %>% unite(col = sample, tissue, mouse, subset)

pct_matrix_corr <- pct_matrix
pct_matrix_corr <- pct_matrix_corr %>% unite(col =  tissue, mouse, subset)

#
for (s in unique(pct_matrix_corr$tissue)) {
  min_temp <- filter(normalized_mins, tissue == s)
  pct_temp <- filter(pct_matrix_corr, tissue == s)
  
  if (nrow(min_temp) == 1) {
    pct_temp$norm_freq[pct_temp$norm_freq == 0] <- min_temp$normalized_half_min
    pct_matrix_corr$norm_freq[pct_matrix_corr$tissue %in% unique(pct_temp$tissue)] <- pct_temp$norm_freq
  }
}

# Cell number cutoff: filter out sample with fewer than 20 cells
# sometimes the column is total, sometimes it is total.y and total.x or currently its total.pancreas

pct_matrix_corr <- left_join(pct_matrix_corr, normalized_mins, by = "tissue")

pct_matrix_corr <- pct_matrix_corr %>%
  filter(total_pancreas >= 20, total_spleen >= 20)


pct_matrix_corr <- pct_matrix_corr %>% separate(tissue, remove = TRUE, sep = "_", into = c( "mouse", "subset"))  %>%
  select(-total_pancreas, -total_spleen, -min_det_pancreas, -min_det_spleen, -normalized_half_min, -tissue_spleen, -normalized_min_det, -freq_pancreas, -freq_spleen)

anyNA(pct_matrix_corr)

#
pct_matrix_corr <- pct_matrix_corr %>% mutate(tissue="Pancreas")


## Simple Treg vs Tconv freq comparison--- --- --- --- ---
graphs.dir <- paste0("./", flowcode.seed, "_graphs/")

if( !file.exists(graphs.dir) ){
  dir.create(graphs.dir)
}

comparison.dir <- "treg_vs_tconv/"
if( !file.exists( file.path(graphs.dir, comparison.dir) ) ){
  dir.create(file.path(graphs.dir, comparison.dir) )
}

# set up collection of p.vals and log2FC
comparison.results <- data.frame()

reference_subset <- "tconv"  # Specify your reference subset

### altered to not loop through tissues
for (t in unique(pct_matrix_corr$tissue)) {
  
  tissue_data <- pct_matrix_corr %>% filter(tissue == t)
  
  for (s in unique(tissue_data$subset)) {
    if (s != reference_subset) {
      subset_data <- tissue_data %>% filter(subset == s)
      ref_data <- tissue_data %>% filter(subset == reference_subset)
      
      resdf = data.frame(target_gene = sort(unique(subset_data$target_gene)), log2FC = 0, pvalue = 1 )
      
      for (i in 1:nrow(resdf)) {
        g = resdf$target_gene[i]
        
        subset_xvdf = subset_data %>% filter(target_gene == g) %>% arrange(mouse)
        ref_yvdf = ref_data %>% filter(target_gene == g) %>% filter(mouse %in% subset_xvdf$mouse) %>% arrange(mouse)
        
        xv = (subset_xvdf %>% filter(mouse %in% ref_yvdf$mouse))$norm_freq
        yv = ref_yvdf$norm_freq
        
        if (length(c(xv,yv)) > 2){
          resdf$log2FC[i] = mean(log2((xv + 0.001) / (yv + 0.001)))
          
          res = t.test(x = log2(xv+0.001),
                       y = log2(yv+0.001),
                       paired = TRUE,
                       alternative = "two.sided")
          resdf$pvalue[i] = res$p.value
          
          ggplot(data.frame(subset_xv = xv, ref_yv = yv, id = 1:length(xv)) %>% 
                   pivot_longer(cols = -id, names_to = "group", values_to = "norm_freq")) +
            scale_x_discrete(name = "", labels = c(s, reference_subset)) +
            scale_y_log10() +
            scale_fill_manual(values = c("darkgrey", "firebrick1")) +
            ggtitle(paste(g, ":", s, "vs", reference_subset, "in", t, "\n log2FC =", round(resdf$log2FC[i], 3), "; pvalue =", round(res$p.value, 3))) +
            geom_boxplot(aes(x = group, y = norm_freq, fill = group), alpha = 0.2) +
            geom_point(aes(x = group, y = norm_freq)) +
            geom_line(aes(x = group, y = norm_freq, group = id)) +
            theme_bw() +
            theme(legend.position = "none")
          
          ggsave(paste0(graphs.dir, comparison.dir, t, "_", s, "_vs_", reference_subset, "_", g, ".pdf"), height = 10, width = 10, units = "cm")
        }
      }
      
      resdf$tissue <- rep(t)
      resdf$subset <- rep(s)
      
      # Apply multiple hypothesis testing correction
      resdf$pvalue <- p.adjust(resdf$pvalue, method = "BH")
      
      maxmlog10pvalue <- max(-log10(resdf$pvalue))
      maxlog2FC <- max(abs(resdf$log2FC))
      
      resdf$dist0 <- sqrt(resdf$log2FC^2 + (-log10(resdf$pvalue) * maxlog2FC / maxmlog10pvalue)^2)
      resdf <- resdf %>% arrange(-dist0)
      
      comparison.results <- rbind(comparison.results, resdf)
      write.csv(resdf, paste0(graphs.dir, comparison.dir, t, "_", s, "_vs_", reference_subset, "_pvals.csv"))
    }
  }
}

#a version of the pval code was removed. its in the old versions


# volcano plots -- - ---- --

# Set the reference subset
reference_subset <- "tconv"
#definitely consider changing x and y axis limits to suit your data

#remove unwanted guides
#remove unwanted guides
genes_to_remove <- c("Itgb2l", "Cd22", "Ctnna2", "Itga10")

# Filter the dataframe
resdf <- resdf[!resdf$target_gene %in% genes_to_remove, ]

# --- volcano option 1

library(ggplot2)
library(ggrepel)
library(dplyr)

volcano <- ggplot(resdf) +
  ggtitle(paste(s, "vs", reference_subset, "in", t, "volcano")) +
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


ggsave(file.path(graphs.dir, comparison.dir, plot.figure), plot = volcano, height = 10, width = 10, units = "cm")



#### end ####





#### CD8 vs Treg####
perCellData <- rbind(CD8_perCellData,treg_perCellData)

#back to original code
# Binding dataframes
pct_matrix <- bind_rows(norm_treg_pcts, norm_CD8_pcts)
pct_matrix<-pct_matrix %>%
  mutate(tissue = "Pancreas")

#remove NAs
pct_matrix <- pct_matrix %>%
  na.omit(data)

# calculate minimum detection frequency per sample from perCellData

tissue_counts <- perCellData %>%
  group_by( tissue, mouse, Id, subset ) %>%
  filter( Id != "Other") %>%
  summarize( count = n()) %>%
  mutate(pct = round(count/sum(count)*100, 2)) 

tissue_mins <- tissue_counts %>%
  group_by( tissue, mouse, subset ) %>%
  summarize( total = sum(count) ) %>%
  mutate( min_det = 100/total ) 

pancreas_mins <- tissue_mins %>% filter(tissue == "Pancreas")
spleen_mins <- tissue_mins %>% filter(tissue == "Spleen")

#joining the pancreas and spleen data
normalized_mins <- pancreas_mins %>%
  inner_join(spleen_mins, by = c("mouse", "subset"), suffix = c("_pancreas", "_spleen"))

#normalizing the pcts
normalized_mins <- normalized_mins %>%
  mutate(normalized_min_det = min_det_pancreas / min_det_spleen,
         normalized_half_min = normalized_min_det / 2)

normalized_mins <- normalized_mins %>%
  rename(tissue=tissue_pancreas)

normalized_mins <- normalized_mins %>% unite(col =  tissue, mouse, subset)

# use half-min frequency to replace zeros in percent matrix

#tissue_mins <- tissue_mins %>% filter(tissue=="Pancreas")

#tissue_mins$half_min <- tissue_mins$min_det/2
#tissue_mins <- tissue_mins %>% unite(col = sample, tissue, mouse, subset)

pct_matrix_corr <- pct_matrix
pct_matrix_corr <- pct_matrix_corr %>% unite(col =  tissue, mouse, subset)

#
for (s in unique(pct_matrix_corr$tissue)) {
  min_temp <- filter(normalized_mins, tissue == s)
  pct_temp <- filter(pct_matrix_corr, tissue == s)
  
  if (nrow(min_temp) == 1) {
    pct_temp$norm_freq[pct_temp$norm_freq == 0] <- min_temp$normalized_half_min
    pct_matrix_corr$norm_freq[pct_matrix_corr$tissue %in% unique(pct_temp$tissue)] <- pct_temp$norm_freq
  }
}

# Cell number cutoff: filter out sample with fewer than 20 cells
# sometimes the column is total, sometimes it is total.y and total.x or currently its total.pancreas

pct_matrix_corr <- left_join(pct_matrix_corr, normalized_mins, by = "tissue")

pct_matrix_corr <- pct_matrix_corr %>%
  filter(total_pancreas >= 20, total_spleen >= 20)


pct_matrix_corr <- pct_matrix_corr %>% separate(tissue, remove = TRUE, sep = "_", into = c( "mouse", "subset"))  %>%
  select(-total_pancreas, -total_spleen, -min_det_pancreas, -min_det_spleen, -normalized_half_min, -tissue_spleen, -normalized_min_det, -freq_pancreas, -freq_spleen)

anyNA(pct_matrix_corr)

#
pct_matrix_corr <- pct_matrix_corr %>% mutate(tissue="Pancreas")


## Simple cd8 vs treg freq comparison--- --- --- --- ---
graphs.dir <- paste0("./", flowcode.seed, "_graphs/")

if( !file.exists(graphs.dir) ){
  dir.create(graphs.dir)
}

comparison.dir <- "cd8_vs_treg/"
if( !file.exists( file.path(graphs.dir, comparison.dir) ) ){
  dir.create(file.path(graphs.dir, comparison.dir) )
}

# set up collection of p.vals and log2FC
comparison.results <- data.frame()

reference_subset <- "treg"  # Specify your reference subset

### altered to not loop through tissues
for (t in unique(pct_matrix_corr$tissue)) {
  
  tissue_data <- pct_matrix_corr %>% filter(tissue == t)
  
  for (s in unique(tissue_data$subset)) {
    if (s != reference_subset) {
      subset_data <- tissue_data %>% filter(subset == s)
      ref_data <- tissue_data %>% filter(subset == reference_subset)
      
      resdf = data.frame(target_gene = sort(unique(subset_data$target_gene)), log2FC = 0, pvalue = 1 )
      
      for (i in 1:nrow(resdf)) {
        g = resdf$target_gene[i]
        
        subset_xvdf = subset_data %>% filter(target_gene == g) %>% arrange(mouse)
        ref_yvdf = ref_data %>% filter(target_gene == g) %>% filter(mouse %in% subset_xvdf$mouse) %>% arrange(mouse)
        
        xv = (subset_xvdf %>% filter(mouse %in% ref_yvdf$mouse))$norm_freq
        yv = ref_yvdf$norm_freq
        
        if (length(c(xv,yv)) > 2){
          resdf$log2FC[i] = mean(log2((xv + 0.001) / (yv + 0.001)))
          
          res = t.test(x = log2(xv+0.001),
                       y = log2(yv+0.001),
                       paired = TRUE,
                       alternative = "two.sided")
          resdf$pvalue[i] = res$p.value
          
          ggplot(data.frame(subset_xv = xv, ref_yv = yv, id = 1:length(xv)) %>% 
                   pivot_longer(cols = -id, names_to = "group", values_to = "norm_freq")) +
            scale_x_discrete(name = "", labels = c(s, reference_subset)) +
            scale_y_log10() +
            scale_fill_manual(values = c("darkgrey", "firebrick1")) +
            ggtitle(paste(g, ":", s, "vs", reference_subset, "in", t, "\n log2FC =", round(resdf$log2FC[i], 3), "; pvalue =", round(res$p.value, 3))) +
            geom_boxplot(aes(x = group, y = norm_freq, fill = group), alpha = 0.2) +
            geom_point(aes(x = group, y = norm_freq)) +
            geom_line(aes(x = group, y = norm_freq, group = id)) +
            theme_bw() +
            theme(legend.position = "none")
          
          ggsave(paste0(graphs.dir, comparison.dir, t, "_", s, "_vs_", reference_subset, "_", g, ".pdf"), height = 10, width = 10, units = "cm")
        }
      }
      
      resdf$tissue <- rep(t)
      resdf$subset <- rep(s)
      
      # Apply multiple hypothesis testing correction
      resdf$pvalue <- p.adjust(resdf$pvalue, method = "BH")
      
      maxmlog10pvalue <- max(-log10(resdf$pvalue))
      maxlog2FC <- max(abs(resdf$log2FC))
      
      resdf$dist0 <- sqrt(resdf$log2FC^2 + (-log10(resdf$pvalue) * maxlog2FC / maxmlog10pvalue)^2)
      resdf <- resdf %>% arrange(-dist0)
      
      comparison.results <- rbind(comparison.results, resdf)
      write.csv(resdf, paste0(graphs.dir, comparison.dir, t, "_", s, "_vs_", reference_subset, "_pvals.csv"))
    }
  }
}

#a version of the pval code was removed. its in the old versions


# volcano plots -- - ---- --

# Set the reference subset
reference_subset <- "treg"
#definitely consider changing x and y axis limits to suit your data

#remove unwanted guides
#remove unwanted guides
genes_to_remove <- c("Itgb2l", "Cd22", "Ctnna2", "Itga10")

# Filter the dataframe
resdf <- resdf[!resdf$target_gene %in% genes_to_remove, ]

# --- volcano option 1

library(ggplot2)
library(ggrepel)
library(dplyr)

volcano <- ggplot(resdf) +
  ggtitle(paste(s, "vs", reference_subset, "in", t, "volcano")) +
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


ggsave(file.path(graphs.dir, comparison.dir, plot.figure), plot = volcano, height = 10, width = 10, units = "cm")




#### ####













# ---- Bar chart of the T cell types per tissue ----

#requires the perCell Data of the T cell types you want to compare already loaded

#  counts the cells
CD8 <- CD8_perCellData %>%
  group_by(tissue) %>%
  summarise(count = n())%>%
  mutate(Cell_type = "CD8") 

Tconv <- tconv_perCellData %>%
  group_by(tissue) %>%
  summarise(count = n())%>%
  mutate(Cell_type = "Tconv") 

Treg <- treg_perCellData %>%
  group_by(tissue) %>%
  summarise(count = n()) %>%
  mutate(Cell_type = "Treg") 

#--- Formatting to make the bar chart ---
result <- rbind (CD8, Tconv, Treg)

# Load necessary libraries
library(ggplot2)
library(reshape2)

result <- result %>%
  group_by(tissue, Cell_type)


result$tissue <- factor(result$tissue, levels = c("Spleen", "LN", "Blood", "Pancreas"))
result$Cell_type <- factor(result$Cell_type, levels = c("CD8", "Tconv", "Treg"))


# Plot the bar chart using ggplot2
plot <- ggplot(result, aes(x = tissue, y = count, fill = Cell_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("CD8" = "#0657B2", "Tconv" = "#22c3e7", "Treg" = "#50ECC4")) + # Customize colors
  theme_minimal() +
  theme(
    text = element_text(size = 16), # Customize font size
    axis.title = element_text(size = 16, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25, color = "black"),
    axis.text.y = element_text(size = 28, color = "black"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid.major.x = element_blank(), # Remove major vertical grid lines
    panel.grid.minor.x = element_blank(), # Remove minor vertical grid lines
    panel.grid.major.y = element_line(size = 0.5, linetype = 'solid', colour = "grey"), # Keep major horizontal grid lines
    panel.grid.minor.y = element_line(size = 0.25, linetype = 'solid', colour = "grey"), # Keep minor horizontal grid lines
    axis.line = element_line(size = 0.5, colour = "black") # Add axis lines
  ) +
  labs(
    title = "Cell Counts per T cell type",
    x = "Tissue",
    y = "Total Cell Count"
  )

# Print the plot
print(plot)

# save
ggsave("plot_cell_count.png", plot = plot, width = 10, height = 7, bg="white")


# --- now by proportion ---
result <- result %>%
  group_by(tissue) %>%
  mutate(tissue_sum = sum(count))

result <- result %>%
  mutate(tissue_pct = 100 * (count / tissue_sum))

plot <- ggplot(result, aes(x = tissue, y = tissue_pct, fill = Cell_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("CD8" = "#0657B2", "Tconv" = "#22c3e7", "Treg" = "#50ECC4")) + # Customize colors
  theme_minimal() +
  theme(
    text = element_text(size = 16), # Customize font size
    axis.title = element_text(size = 16, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25, color = "black"),
    axis.text.y = element_text(size = 28, color = "black"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid.major.x = element_blank(), # Remove major vertical grid lines
    panel.grid.minor.x = element_blank(), # Remove minor vertical grid lines
    panel.grid.major.y = element_line(size = 0.5, linetype = 'solid', colour = "grey"), # Keep major horizontal grid lines
    panel.grid.minor.y = element_line(size = 0.25, linetype = 'solid', colour = "grey"), # Keep minor horizontal grid lines
    axis.line = element_line(size = 0.5, colour = "black") # Add axis lines
  ) +
  labs(
    title = "Percentages per T cell type",
    x = "Tissue Group",
    y = "% of Tissue Sample"
  )

# Print the plot
print(plot)

# save
ggsave("plot_proportion.png", plot = plot, width = 10, height = 7, bg="white")




















# --- heatmap ----
#this part doesn't work very well 

# Prepare the data
heatmap.input <- comparison.results


# Define order of tissues based on groups
tissue_order <- c("Blood","LN","Pancreas","Spleen")

# Calculate the average log2FC for each target_gene across all tissues
average_log2FC_between_tissues <- aggregate(log2FC ~ target_gene, data = heatmap.input, FUN = mean)

# Rank the target genes based on the average log2FC
rank_log2FC <- rank(average_log2FC_between_tissues$log2FC)

# Create a data frame with target_gene and rank_log2FC
ranked_genes <- data.frame(target_gene = average_log2FC_between_tissues$target_gene, 
                           rank_log2FC = rank_log2FC,
                           avg_log2FC = average_log2FC_between_tissues$log2FC)

# Sort the ranked genes data frame based on avg_log2FC
ranked_genes <- ranked_genes[order(ranked_genes$avg_log2FC), ]

# Set log2FC to 0 for any p-value greater than 0.01
heatmap.input[heatmap.input$pvalue > 0.01, "log2FC"] <- 0

# Reshape data into wide format: one row per tissue, one column per target_gene
heatmap.input <- heatmap.input %>%
  select(tissue, target_gene, log2FC) %>%
  pivot_wider(names_from = target_gene, values_from = log2FC)

# Convert to matrix and set rownames as tissue
heatmap.rownames <- heatmap.input$tissue
heatmap.input <- heatmap.input %>% select(-tissue)
heatmap.input <- as.matrix(heatmap.input)
rownames(heatmap.input) <- heatmap.rownames

# Validate gene and tissue names for subsetting
valid_genes <- ranked_genes$target_gene[ranked_genes$target_gene %in% colnames(heatmap.input)]
valid_tissue_order <- tissue_order[tissue_order %in% rownames(heatmap.input)]

# Define color palette and breaks for the heatmap
palette.length <- 50
heatmap.colors <- colorRampPalette(c("blue", "white", "red"))(palette.length)
min_value <- min(heatmap.input, na.rm = TRUE)
max_value <- max(heatmap.input, na.rm = TRUE)
heatmap.breaks <- c(seq(min_value, 0, length.out = ceiling(palette.length / 2) + 1), 
                    seq(max_value / palette.length, max_value, length.out = floor(palette.length / 2)))

# Create and save the heatmap directly within pheatmap function
pheatmap(
  heatmap.input[valid_tissue_order, valid_genes],  # Subset using valid genes and tissues
  color = heatmap.colors,
  breaks = heatmap.breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  #border_color = "black",
  filename = paste0(graphs.dir, "Frequency_heatmap_p0.01.png"), # Save directly within pheatmap
  height = 10,
  width = 20,
  units = "cm"
)















