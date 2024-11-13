

#qPCR


setwd("C:/Users/15869/Desktop/CPP_LOCO_GF/qPCR/")
meta_ob1  <- read.table("/Users/15869/Desktop/CPP_LOCO_GF/qPCR/LOCO_GF_CPmg.txt", header=T, row.names=1)

GF_df <- meta_ob1[meta_ob1$Study=="GF",]

nrow(GF_df)

library(ggplot2)

GF_qPCR_plot <- ggplot(GF_df, aes(x = Group, y = log10(CPmg), color = factor(Group), fill = factor(Group))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, show.legend = FALSE, color = "black") +
  geom_point(position = position_jitterdodge(jitter.width = 0.6666, dodge.width = 0.7, seed = NULL), 
             size = 3, alpha = 0.75, stroke = 1, shape = 21, color = "black", show.legend = FALSE) +
  scale_shape_manual(values = c(22, 22, 22, 22, 22, 22, 22, 22)) +
  scale_color_manual(values = c("#FA0000", "#005DE6", "purple", "black", "magenta", "cyan", "pink", "orange", "brown")) +
  scale_fill_manual(values = c("#FA0000", "#005DE6", "purple", "black", "magenta", "cyan", "pink", "gray", "orange", "brown")) +
  theme_classic(base_size = 11) +
  theme(text = element_text(color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) +
  ylab("log(16S rRNA gene copies)") +
  scale_y_continuous(limits = c(0.0, 9),
                     labels = scales::number_format(accuracy = 0.1),
                     breaks = scales::pretty_breaks(n = 10)) +
  facet_wrap(~Time, scales = "free_x", ncol = 2, nrow = 1)

GF_qPCR_plot






###Sample size

setwd("C:/Users/15869/Desktop/CPP_LOCO_GF/16S/")

#####Germ-free Study

# Read the data from the file

meta_ob1 <- read.table("C:/Users/15869/Desktop/CPP_LOCO_GF/16S/CPP_LOCO_GF_meta.txt", header = TRUE, row.names = 1)

meta_ob1 <- meta_ob1[meta_ob1$Study== "GF",]
meta_ob1 <- meta_ob1[meta_ob1$Time== "End",]

unique(meta_ob1$Time)
unique(meta_ob1$Group)

nrow(meta_ob1)


# Load the dplyr library
library(dplyr)
# Group by Week, Type, and Treatment, then count the samples
result <- meta_ob1 %>%
  group_by(Study, Group, Sex, Sex_Group, Time) %>%
  summarise(count = n())

# Print the result
print(result, n=100)

#write.csv (result, file = "CPP_counts.csv")





###Alpha diversity

meta_ob1 <- read.table("C:/Users/15869/Desktop/CPP_LOCO_GF/16S/CPP_LOCO_GF_meta.txt", header = TRUE, row.names = 1)

meta_ob1 <- meta_ob1[meta_ob1$Study== "GF",]
meta_ob1 <- meta_ob1[meta_ob1$Time== "End",]


library(phyloseq)
library(metagMisc)
library(vegan)
library(ggplot2)

tax_ob1 <- import_mothur(mothur_constaxonomy = ("C:/Users/15869/Desktop/CPP_LOCO_GF/16S/CPP_LOCO_GF_taxa.txt"))

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("CPP_LOCO_GF_shared.csv",)
#fullMothur<-read.csv("mini_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtu))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

# Remove ASVs with zero counts across all samples
AF_mothur2 <- AF_mothur2[, colSums(AF_mothur2) > 0]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))


# Extract the taxonomy table
tax_table_df <- as.data.frame(tax_table(phyloseqobj.f))

# Remove taxa at Rank1 that are Archaea
tax_table_df <- tax_table_df[!(tax_table_df$Rank6 == "Remove"), ]

# Update the taxonomy table in the phyloseq object
tax_table(phyloseqobj.f) <- as.matrix(tax_table_df)

s_size <- min(sample_sums(phyloseqobj.f))

# Rename Firmicutes to Bacillota in the phyloseq object
tax_table(phyloseqobj.f)[tax_table(phyloseqobj.f) == "Firmicutes"] <- "Bacillota"

phyloseqobj.f2 <- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = s_size, rngseed = 1, replace = FALSE)

richness_out <- estimate_richness(phyloseqobj.f2, split = TRUE, measures = c("Chao1", "Shannon", "InvSimpson"))


otu_table <- otu_table(phyloseqobj.f2)
otu_df <- as.data.frame(otu_table)
calculate_good_coverage <- function(sample_counts) {
  total_reads <- sum(sample_counts)
  singletons <- sum(sample_counts == 1)
  coverage <- 1 - (singletons / total_reads)
  return(coverage)
}

Coverage_out <- as.data.frame(apply(otu_df, 1, calculate_good_coverage))

min(Coverage_out)


meta_ob2 <- cbind(meta_ob1$Sex, meta_ob1$Group, meta_ob1$Treatment, Coverage_out, richness_out$Chao1, richness_out$Shannon, richness_out$InvSimpson, meta_ob1$Sex_Group)
colnames(meta_ob2)[1] <- "Sex"
colnames(meta_ob2)[2] <- "Group"
colnames(meta_ob2)[3] <- "Treatment"
colnames(meta_ob2)[4] <- "Coverage"
colnames(meta_ob2)[5] <- "Chao1"
colnames(meta_ob2)[6] <- "Shannon"
colnames(meta_ob2)[7] <- "InvSimpson"
colnames(meta_ob2)[8] <- "Sex_Group"


library(ggplot2)

Coverage_plot <- ggplot(meta_ob2, aes(x=Sex_Group, y=Coverage, color=Sex_Group, fill=Sex_Group)) +
  geom_boxplot(outlier.shape = NA, alpha=0.5, show.legend = FALSE, color = "black") +
  geom_point(position=position_jitterdodge(jitter.width = 0.333, jitter.height = 0, dodge.width = 0.75, seed = NULL), 
             size=3, alpha=0.75, stroke=1, shape=21, color="black", show.legend = FALSE) +
  scale_shape_manual(values = c(22, 22, 22, 22, 22, 22, 22, 22)) +
  scale_color_manual(values=c("#FED2D2", "#88cafc","#FA0000", "#005DE6", "purple","black","magenta","cyan","pink","orange","brown")) +
  scale_fill_manual(values=c("#FED2D2", "#88cafc","#FA0000", "#005DE6", "purple","black","magenta","cyan","pink","gray","orange","brown")) +
  theme_classic(base_size = 11) +
  theme(text = element_text(color = "black"),   
        axis.title.x = element_blank(),    
        #axis.text.x = element_blank(),     
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))



chao1_plot <- ggplot(meta_ob2, aes(x=Sex_Group, y=Chao1, color=Sex_Group, fill=Sex_Group)) +
  geom_boxplot(outlier.shape = NA, alpha=0.5, show.legend = FALSE, color = "black") +
  geom_point(position=position_jitterdodge(jitter.width = 0.333, jitter.height = 0, dodge.width = 0.75, seed = NULL), 
             size=3, alpha=0.75, stroke=1, shape=21, color="black", show.legend = FALSE) +
  scale_shape_manual(values = c(22, 22, 22, 22, 22, 22, 22, 22)) +
  scale_color_manual(values=c("#FED2D2", "#88cafc","#FA0000", "#005DE6", "purple","black","magenta","cyan","pink","orange","brown")) +
  scale_fill_manual(values=c("#FED2D2", "#88cafc","#FA0000", "#005DE6", "purple","black","magenta","cyan","pink","gray","orange","brown")) +
  theme_classic(base_size = 11) +
  theme(text = element_text(color = "black"),   
        axis.title.x = element_blank(),    
        #axis.text.x = element_blank(),     
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))


shannon_plot <- ggplot(meta_ob2, aes(x=Sex_Group, y=Shannon, color=Sex_Group, fill=Sex_Group)) +
  geom_boxplot(outlier.shape = NA, alpha=0.5, show.legend = FALSE, color = "black") +
  geom_point(position=position_jitterdodge(jitter.width = 0.333, jitter.height = 0, dodge.width = 0.75, seed = NULL), 
             size=3, alpha=0.75, stroke=1, shape=21, color="black", show.legend = FALSE) +
  scale_shape_manual(values = c(22, 22, 22, 22, 22, 22, 22, 22)) +
  scale_color_manual(values=c("#FED2D2", "#88cafc","#FA0000", "#005DE6", "purple","black","magenta","cyan","pink","orange","brown")) +
  scale_fill_manual(values=c("#FED2D2", "#88cafc","#FA0000", "#005DE6", "purple","black","magenta","cyan","pink","gray","orange","brown")) +
  theme_classic(base_size = 11) +
  theme(text = element_text(color = "black"),   
        axis.title.x = element_blank(),    
        #axis.text.x = element_blank(),     
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))


invsimpson_plot <- ggplot(meta_ob2, aes(x=Sex_Group, y=InvSimpson, color=Sex_Group, fill=Sex_Group)) +
  geom_boxplot(outlier.shape = NA, alpha=0.5, show.legend = FALSE, color = "black") +
  geom_point(position=position_jitterdodge(jitter.width = 0.333, jitter.height = 0, dodge.width = 0.75, seed = NULL), 
             size=3, alpha=0.75, stroke=1, shape=21, color="black", show.legend = FALSE) +
  scale_color_manual(values=c("#FED2D2", "#88cafc","#FA0000", "#005DE6", "purple","black","magenta","cyan","pink","orange","brown")) +
  scale_fill_manual(values=c("#FED2D2", "#88cafc","#FA0000", "#005DE6", "purple","black","magenta","cyan","pink","gray","orange","brown")) +
  theme_classic(base_size = 11) +
  theme(text = element_text(color = "black"),   
        axis.title.x = element_blank(),    
        #axis.text.x = element_blank(),     
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black")) +
  ylab("Inverse Simpson")



chao1_plot
shannon_plot
invsimpson_plot




library(car)

library(emmeans)
options(max.print=1000000)
options(digits = 5)
emm_options(opt.digits = FALSE)

chao_model <- glm(Chao1 ~ Sex*Group, data = meta_ob2)
Anova(chao_model, type=3, test.statistic="F")
#Females with gretaer rcichness than males

shannon_model <- glm(Shannon ~ Sex*Group, data = meta_ob2)
Anova(shannon_model, type=3, test.statistic="F")
#Both Sex and Group results largely driven by high shannon values in female control mice

invsimpson_model <- glm(InvSimpson ~ Sex*Group, data = meta_ob2)
Anova(invsimpson_model, type=3, test.statistic="F")
#Both Sex and Group results largely driven by high shannon values in female control mice



library(emmeans)
options(max.print=1000000)
options(digits = 5)
emm_options(opt.digits = FALSE)


chao_emm <- emmeans(chao_model, ~ Sex*Group)
chao_tukey_out <- pairs(chao_emm, adjust = "tukey")
#write.csv(chao_tukey_out, file = "CPP_chao_tukey.csv")

shannon_emm <- emmeans(shannon_model, ~ Sex*Group)
shannon_tukey_out <- pairs(shannon_emm, adjust = "tukey")
#write.csv(shannon_tukey_out , file = "CPP_shannon_tukey.csv")

invsimpson_emm <- emmeans(invsimpson_model, ~ Sex*Group)
invsimpson_tukey_out <- pairs(invsimpson_emm, adjust = "tukey")
#write.csv(shannon_tukey_out , file = "CPP_shannon_tukey.csv")






#beta diversity


merged_phy_percent <- phyloseq_standardize_otu_abundance(phyloseqobj.f2, method = "total")
merged_phy_percent1 <- tax_glom(merged_phy_percent, taxrank="Rank7")


# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent1) {
  sd <- sample_data(merged_phy_percent1)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent1) {
  OTU <- otu_table(merged_phy_percent1)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


taxa_names(merged_phy_percent1)
tax_table(merged_phy_percent1)[,7]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,7]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
nrow(shared_phy_out1)

input_data  <- as.data.frame(shared_phy_out1)
input_metadata <- as.data.frame(meta_phy_out1)

nrow(input_metadata)



# Load necessary libraries
library(vegan)
library(ggplot2)
library(glue)
library(dplyr)
library(ggforce)


# Perform PCoA
dist_bray <- capscale(input_data~ 1, distance="bray", binary=TRUE, eig=TRUE)

# Extract proportion explained by the first two PCoA axes
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)

# Calculate scores for sites
data.scores <- scores(dist_bray, display = "sites", choices = c(1:3))
df <- data.frame(data.scores)

# Add Group and Sex columns to df dataframe
df$Sex_Group <- factor(input_metadata$Sex_Group)

# Labels for axes
labx1 <- glue("PCo 1 ({prop_explained[1]}%)")
laby1 <- glue("PCo 2 ({prop_explained[2]}%)")



# Extract scores for taxa (species)
taxa_scores <- scores(dist_bray, display = "species", choices = c(1, 2))
taxa_df <- as.data.frame(taxa_scores)

# Select top n taxa based on distance from origin (sqrt(MDS1^2 + MDS2^2))
taxa_df$distance <- sqrt(taxa_df$MDS1^2 + taxa_df$MDS2^2)
top_taxa <- taxa_df %>% 
  arrange(desc(distance)) %>% 
  head(3)
top_taxa$Taxa <- rownames(top_taxa)  # Add taxa names


# Create the plot

library(ggforce)

GF_Bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Sex_Group, fill = Sex_Group, shape = Sex_Group)) +
  scale_color_manual(values = c("#FED2D2", "#88cafc", "#FA0000", "#005DE6", "purple", "black", "magenta", "cyan", "pink", "orange", "brown")) +
  scale_fill_manual(values = c("#FED2D2", "#88cafc", "#FA0000", "#005DE6", "purple", "black", "magenta", "cyan", "pink", "gray", "orange", "brown")) +
  geom_point(aes(shape = Sex_Group), size = 3, alpha = 1, stroke = 1, color = "black") +
  scale_shape_manual(values = c(21, 21, 21, 21, 21, 21, 21, 21)) +  # Use shapes that support fill and stroke
  ggforce::geom_mark_ellipse(
    aes(fill = Sex_Group, color = Sex_Group),
    expand = unit(7.0, "mm"),   # Adjust the spacing between points and ellipse border
    label.fontsize = 10,      # Font size for the label inside the ellipse
    label.colour = "black",   # Color of the label text
    con.colour = "gray",      # Color of the connection line
    con.size = 0.5,           # Thickness of the connection line
    show.legend = TRUE,
    alpha= 0.5) +
  theme_classic() +
  labs(x = labx1, y = laby1) +
  theme(axis.text.x = element_text(color = "black"), 
        text = element_text(size = 11), 
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(size = 9, face = "plain"),
        axis.title.y = element_text(size = 9, face = "plain"),
        legend.title = element_blank(), 
        strip.text.x = element_blank(),
        legend.direction = 'vertical',
        legend.position = "bottom",
        legend.margin = margin(t = -10, r = 0, b = 0, l = -30),  # Shift legend slightly left with a negative left margin
        legend.text = element_text(size = 9)) +
  guides(color = guide_legend(ncol = 2),
         fill = guide_legend(ncol = 2),
         shape = guide_legend(ncol = 2)) +
  geom_text(data = top_taxa, aes(x = MDS1, y = MDS2, label = Taxa), color = "black", size = 3.5, fontface = "plain", inherit.aes = FALSE) +
  scale_x_continuous(limits = c(-0.9, 1.9), breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(limits = c(-1.9, 1.35), breaks = scales::pretty_breaks(n = 6))

GF_Bray_plot





library(vegan)
dist.bray <- vegdist(input_data, method="bray", binary=FALSE)
adonis2(formula = dist.bray  ~ Group*Sex, data = input_metadata, permutations=10000)


library(pairwiseAdonis)
x <- phyloseq::distance(merged_phy_percent1, method ="bray")
factors <- meta_phy_out1$Sex_Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="none"))








###Barplot

library(ggbreak)

create_bar_plot <- function(phyloseq_obj, tax_rank, num_taxa, group_factor, facet_factor = NULL, 
                            fill_colors = NULL, y_break1 = c(NULL, NULL), y_break2 = c(NULL, NULL), 
                            y_break3 = c(NULL, NULL), y_accuracy = NULL, n_breaks = NULL, 
                            y_limits = c(NULL, NULL), custom_taxa_list = NULL) {
  
  # Extract the taxonomy table
  tax_table_df <- as.data.frame(tax_table(phyloseq_obj))
  
  # Remove taxa at Rank1 that are Archaea
  tax_table_df <- tax_table_df[!(tax_table_df$Rank6 == "Remove"), ]
  
  # Update the taxonomy table in the phyloseq object
  tax_table(phyloseq_obj) <- as.matrix(tax_table_df)
  
  # Standardize OTU abundance
  merged_phy_percent <- phyloseq_standardize_otu_abundance(phyloseq_obj, method = "total")
  
  # Taxonomic glomming at the specified rank
  merged_phy_percent2 <- tax_glom(merged_phy_percent, taxrank = tax_rank)
  
  # Prune to the top 'num_taxa' taxa based on abundance
  merged_phy_percent2 <- prune_taxa(names(sort(taxa_sums(merged_phy_percent2), TRUE)[1:num_taxa]), merged_phy_percent2)
  
  # Melt the phyloseq object for ggplot
  physeq_df <- psmelt(merged_phy_percent2)
  
  # Filter for custom taxa list if provided
  if (!is.null(custom_taxa_list)) {
    physeq_df <- physeq_df %>% filter(!!sym(tax_rank) %in% custom_taxa_list)
  }
  
  # Calculate mean and standard error for each taxon within each sex/group
  summary_df <- physeq_df %>%
    group_by(!!sym(group_factor), !!sym(tax_rank)) %>%
    summarise(
      Mean_Abundance = mean(Abundance),
      SE_Abundance = sd(Abundance) / sqrt(n()),  # Standard Error
      .groups = 'drop'  # Avoid grouping after summarise
    )
  
  # Print the summary data frame for each taxon
  print(summary_df, n=100)
  
  # Sort summary_df by Mean Abundance in descending order
  summary_df <- summary_df %>%
    arrange(desc(Mean_Abundance)) %>%
    mutate(!!sym(tax_rank) := factor(!!sym(tax_rank), levels = unique(!!sym(tax_rank))))  # Set factor levels for the taxon rank
  
  # Define default colors if not provided
  if (is.null(fill_colors)) {
    fill_colors <- c("#FED2D2", "#88cafc", "#FA0000", "#005DE6", "purple", "black", "magenta", "cyan", "pink", "gray", "orange", "brown")
  }
  
  # Create a bar plot with mean abundance and error bars
  bar_plot <- ggplot(summary_df, aes_string(x = tax_rank, y = "Mean_Abundance", fill = group_factor)) +
    geom_bar(stat = "identity", position = "dodge", size = 0.25, alpha = 1, color = "black") +  # Add black border around bars
    geom_errorbar(aes(ymin = Mean_Abundance - SE_Abundance, ymax = Mean_Abundance + SE_Abundance), 
                  position = position_dodge(0.9), width = 0.25, color = "black", size = 0.3) +  # Set error bar color to black
    labs(y = "Relative\nAbundance", x = NULL) +  # Label the y-axis and remove x-axis label
    scale_color_manual(values = fill_colors) +  # Custom colors for fill
    scale_fill_manual(values = fill_colors) +  # Custom colors for fill
    theme_classic() +  # Add classic theme
    theme(
      legend.position = "bottom",
      legend.margin = margin(t = -10, r = 0, b = 0, l = -10),  # Shift legend slightly left with a negative left margin
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      text = element_text(size = 11),  # Set general text size
      axis.title.y = element_text(size = 10),
      legend.key.size = unit(0.375, 'cm'),
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = unit(c(0.0, 0.0, 0.0, 0.0), 'cm') # Set plot margins (top, right, bottom)
    ) +
    scale_y_continuous(
      limits = y_limits,  # Set y-axis limits
      labels = scales::number_format(accuracy = y_accuracy),
      breaks = scales::pretty_breaks(n = n_breaks),
      sec.axis = sec_axis(~., labels = NULL)  # Keeps the primary y-axis and removes secondary axis labels
    ) +
    guides(
      colour = guide_legend(ncol = 2),  # Make the legend two columns
      fill = guide_legend(ncol = 2),
      shape = guide_legend(ncol = 2)
    )
  
  # Remove the right y-axis and ticks
  bar_plot <- bar_plot + 
    theme(
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.line.y.right = element_blank()
    )
  
  # Apply scale_y_breaks for each y_break if at least one value is not NULL
  if (!is.null(y_break1[1]) && !is.null(y_break1[2])) {
    bar_plot <- bar_plot + scale_y_break(c(y_break1[1], y_break1[2]), scales = 1)
  }
  
  if (!is.null(y_break2[1]) && !is.null(y_break2[2])) {
    bar_plot <- bar_plot + scale_y_break(c(y_break2[1], y_break2[2]), scales = 1)
  }
  
  if (!is.null(y_break3[1]) && !is.null(y_break3[2])) {
    bar_plot <- bar_plot + scale_y_break(c(y_break3[1], y_break3[2]), scales = 1)
  }
  
  # Add faceting if a faceting factor is provided
  if (!is.null(facet_factor)) {
    bar_plot <- bar_plot + facet_wrap(as.formula(paste("~", facet_factor)))
  }
  
  # Return both the summary data frame and the bar plot object
  return(list(summary_df = summary_df, bar_plot = bar_plot))
}

# Example usage
custom_taxa <- c("Akkermansiaceae", "Ruminococcaceae")

GF_barplot_plot <- create_bar_plot(phyloseqobj.f2,
                                  tax_rank = "Rank5",
                                  num_taxa = 20,
                                  group_factor = "Sex_Group",
                                  facet_factor = NULL,
                                  fill_colors = c("#FED2D2", "#88cafc", "#FA0000", "#005DE6"),
                                  y_break1 = c(NULL, NULL),
                                  y_break2 = c(NULL, NULL),
                                  y_break3 = c(NULL, NULL),
                                  y_accuracy = 0.01,  # Set desired accuracy for y-axis
                                  n_breaks = 7,       # Set desired number of breaks for y-axis
                                  y_limits = c(0, 0.35),  # Set y-axis limits
                                  custom_taxa_list = custom_taxa)  # Custom taxa list

# Access the summary data frame and bar plot
summary_df <- GF_barplot_plot$summary_df
GF_barplot_object <- GF_barplot_plot$bar_plot

GF_barplot_object



# Filter for Bacilli and calculate the mean abundance
#mean_abundance_taxon <- summary_df %>%
  #filter(Rank4 == "Peptococcales") %>%
  #summarise(mean_abundance = mean(Mean_Abundance))


#abund_val <- mean_abundance_taxon * s_size
#abund_val







###Merge plots

library(ggpubr)

A_C_plots <- ggarrange(GF_qPCR_plot, chao1_plot, shannon_plot,
                       labels = c("A", "B", "C"),
                       common.legend = FALSE, legend = "bottom",
                       ncol = 3, nrow = 1)


D_F_plots <- ggarrange(invsimpson_plot, GF_Bray_plot, GF_barplot_object,
                       labels = c("D", "E", "F"),
                       common.legend = FALSE, legend = "bottom",
                       ncol = 3, nrow = 1)


GF_A_F_plots <- ggarrange(A_C_plots, D_F_plots,
                          common.legend = FALSE, legend = "bottom",
                          ncol = 1, nrow = 2)



width_in_inches <- 900 / 96
height_in_inches <- 525 / 96

ggsave("GF_A_F_plots.jpeg", plot = GF_A_F_plots, width = width_in_inches, height = height_in_inches, dpi = 600)





###Differential abundance analysis


###Maaslin2

merged_phy_percent1 <- tax_glom(phyloseqobj.f2 , taxrank="Rank5")

# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent1) {
  sd <- sample_data(merged_phy_percent1)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent1) {
  OTU <- otu_table(merged_phy_percent1)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


taxa_names(merged_phy_percent1)
tax_table(merged_phy_percent1)[,5]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,5]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
nrow(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)
row.names(input_metadata)



library(Maaslin2)

fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    min_abundance=500,
                    min_prevalence=0.34,
                    #min_variance=0.100,
                    normalization  = "TMM",
                    output         = "maaslin2_GF_Sex_Group_Rank5_500_34_110724",
                    fixed_effects  = c("Group", "Sex"),
                    reference      = c("Group,Conv", "Sex,F"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)





###Heatmap

library(ComplexHeatmap)
library(vegan)
library(RColorBrewer)


generate_taxa_heatmap <- function(phyloseq_obj, tax_rank = NULL, num_taxa = NULL, grouping_category = NULL, group_colors = NULL, column_order = NULL) {
  # Standardize OTU abundance
  merged_phy_percent <- phyloseq_standardize_otu_abundance(phyloseq_obj, method = "total")
  
  # Glom taxa by specified taxonomic rank
  merged_phy_percent1 <- tax_glom(merged_phy_percent, taxrank = tax_rank)
  
  # Prune the top 'num_taxa' based on abundance
  top_taxa_names <- names(sort(taxa_sums(merged_phy_percent1), TRUE)[1:num_taxa])
  merged_phy_percent2 <- prune_taxa(top_taxa_names, merged_phy_percent1)
  
  # Extract Rank labels for the pruned taxa
  rank_labels <- tax_table(merged_phy_percent1)[taxa_names(merged_phy_percent2), tax_rank]
  taxa_names(merged_phy_percent2) <- rank_labels
  
  # Extract metadata from the phyloseq object
  sample_data_df <- as.data.frame(sample_data(merged_phy_percent2))
  
  # Generate the heatmap_data from the phyloseq object and transpose it
  phyloseq_relative <- transform_sample_counts(merged_phy_percent2, function(x) x / sum(x))
  otu_table_matrix <- as.data.frame(otu_table(phyloseq_relative))
  heatmap_data <- as.matrix(otu_table_matrix)
  
  # Transpose the heatmap data
  heatmap_data <- t(heatmap_data)
  
  # Reorder columns if column_order is specified
  if (!is.null(column_order)) {
    if (is.character(column_order)) {
      column_order <- sample_data_df[[column_order]]
    }
    
    # Get unique groups in the column_order
    unique_groups <- unique(column_order)
    
    # Initialize an empty vector to store the new column order
    new_column_order <- c()
    
    # Loop through each group and perform clustering
    for (group in unique_groups) {
      group_indices <- which(column_order == group)
      group_data <- heatmap_data[, group_indices]
      
      # Perform clustering on the group data using Bray-Curtis distance and Ward's method
      dist_matrix <- vegdist(t(group_data), method = "bray")
      cluster_result <- hclust(dist_matrix, method = "ward.D2")
      
      # Get the order of columns based on clustering
      group_order <- group_indices[cluster_result$order]
      
      # Append the ordered columns to the new column order
      new_column_order <- c(new_column_order, group_order)
    }
    
    # Reorder the heatmap data and sample data based on the new column order
    heatmap_data <- heatmap_data[, new_column_order]
    sample_data_df <- sample_data_df[new_column_order, ]
  }
  
  # Create a mapping of colors for the grouping category
  if (is.null(group_colors)) {
    stop("Please provide a named vector of colors for the groups.")
  }
  
  # Ensure that all levels in the grouping category have corresponding colors
  unique_groups <- unique(sample_data_df[[grouping_category]])
  missing_colors <- setdiff(unique_groups, names(group_colors))
  if (length(missing_colors) > 0) {
    stop(paste("Missing colors for groups:", paste(missing_colors, collapse = ", ")))
  }
  

  # Group Annotation
  group_annotation <- HeatmapAnnotation(
    grouping_category = sample_data_df[[grouping_category]],
    col = list(grouping_category = group_colors),
    annotation_legend_param = list(grouping_category = list(direction = "vertical", nrow = 2)),
    gp = gpar(col = "black"),
    show_legend = TRUE)
  
  # Define the heatmap
  heatmap_plot <- Heatmap(
    heatmap_data,
    name = "Relative Abundance",  # Title for the heatmap legend
    col = colorRampPalette(brewer.pal(9, "Purples"))(100),  # Color scale for relative abundance
    cluster_rows = TRUE,  # Enable row clustering
    clustering_distance_rows = function(x) vegdist(log(x + 1), method = "bray"),  # Distance metric
    clustering_method_rows = "ward.D2",  # Clustering method for rows
    cluster_columns = FALSE,  # Disable column clustering
    show_row_names = TRUE,  # Show row names
    show_column_names = FALSE,  # Hide column names
    show_row_dend = TRUE,
    show_column_dend = FALSE,
    row_dend_width = unit(2.5, "mm"),
    top_annotation = group_annotation,  # Add group annotations
    heatmap_legend_param = list( title = "Relative Abundance",
                                 direction = "horizontal",
                                 title_gp = gpar(fontface = "plain", fontsize = 10), # Make the legend title plain and reduce font size
                                 labels_gp = gpar(fontsize = 8)), # Reduce font size of the legend scale values
    row_names_gp = gpar(fontsize = 10)) # Reduce font size of row names
  
  return(heatmap_plot)
}

custom_colors <- c("B_Group_Female" = "#FA0000", "B_Group_Male" = "#005DE6", "A_Group_Female" = "#FED2D2", "A_Group_Male" = "#88cafc")


GF_heatmap_object <- generate_taxa_heatmap(phyloseqobj.f2,
                                        tax_rank = "Rank7",
                                        num_taxa = 20,
                                        grouping_category = "Sex_Group",
                                        group_colors = custom_colors,
                                        column_order = "Sex_Group")


GF_heatmap_plot <- draw(GF_heatmap_object, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE)

GF_heatmap_grob <- grid.grabExpr(draw(GF_heatmap_plot , heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE))

GF_heatmap_plot




###Pathway plots

setwd("C:/Users/15869/Desktop/CPP_LOCO_GF/16S/GF_PWY/")

meta_ob1 <- read.table("C:/Users/15869/Desktop/CPP_LOCO_GF/16S/GF_PWY/GF_PWY_for_volcano.txt", header = TRUE, row.names = 1)

meta_ob1 <- meta_ob1[meta_ob1$metadata=="Group",]

data_df <- meta_ob1

# Define p_cutoff and F_cutoff
p_cutoff <- 10e-32
F_cutoff <- 5.0  # or any threshold you want to set for FC (fold change)

library(EnhancedVolcano)

# EnhancedVolcano plot
# Create a named vector of colors based on multiple conditions
custom_colors <- ifelse(data_df$FC >= F_cutoff & data_df$qval < p_cutoff, "red",
                        ifelse(data_df$FC <= -F_cutoff & data_df$qval < p_cutoff, "blue", "gray"))
names(custom_colors) <- data_df$PWY  # Ensure that the names correspond to the labels in your data

# EnhancedVolcano plot
GF_volcano_plot <- EnhancedVolcano(data_df,
                                   lab = "",
                                   x = 'FC',
                                   y = 'qval',
                                   colAlpha = 1,
                                   title = NULL,
                                   subtitle = NULL,
                                   pointSize = 2.5,
                                   labSize = 3,
                                   drawConnectors = TRUE,
                                   typeConnectors = 'closed',
                                   endsConnectors = 'last',
                                   legendPosition = 'none',
                                   pCutoff = p_cutoff,
                                   FCcutoff = F_cutoff,
                                   axisLabSize = 8,
                                   caption = NULL,
                                   max.overlaps = 15,
                                   colCustom = custom_colors,
                                   xlim = c(-7.5, 10)) +
  theme(plot.margin = unit(c(0, 0.2, 0, 1), 'lines'))


# Step 1: Filter pathways based on volcano plot conditions (FC >= F_cutoff & qval < p_cutoff for red, FC <= -F_cutoff & qval < p_cutoff for blue)
meta_filtered_selected <- data_df %>%
  filter((FC >= F_cutoff & qval < p_cutoff) | (FC <= -F_cutoff & qval < p_cutoff))

# Step 2: Select the top 10 most positive and 10 most negative pathways from meta_filtered_selected
top10_positive_selected <- meta_filtered_selected %>%
  filter(FC >= F_cutoff) %>%
  arrange(desc(FC)) %>%
  head(200)

top10_negative_selected <- meta_filtered_selected %>%
  filter(FC <= -F_cutoff) %>%
  arrange(FC) %>%
  head(200)

# Step 3: Combine the selected top 10 positive and negative pathways into one dataset
top_combined_selected <- bind_rows(top10_positive_selected, top10_negative_selected)

# Step 4: Create the bar plot with the combined dataset
GF_bar_plot_selected <- ggplot(top_combined_selected, aes(x = FC, y = reorder(Name, FC))) +
  geom_bar(stat = "identity", aes(fill = FC > 0), width = 0.7, color = "black") +  # Bars colored by positive/negative FC
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Vertical line at x = 0
  scale_fill_manual(values = c("blue", "red"), 
                    labels = c("Lower in ABX", "Higher in ABX")) +
  labs(fill = "Significance",
       x = "Log2 fold difference") +  # Add a title for the legend
  scale_x_continuous(limits = c(-7.5, 1.0), labels = scales::number_format(accuracy = 1), breaks = scales::pretty_breaks(n = 5)) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 8),
        legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 8, face = "plain"),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0))  # Remove extra margins


GF_volcano_plot
GF_bar_plot_selected


library(ggpubr)

# Merge the volcano plot and the bar plot
GF_PWY_plots_merged <- ggarrange(GF_volcano_plot, GF_bar_plot_selected,
                                 labels = c("A", "B"),
                                 common.legend = FALSE,
                                 ncol = 2, nrow = 1,
                                 widths = c(0.6666, 1.0))









###PWY Maaslin2


setwd("C:/Users/15869/Desktop/CPP_LOCO_GF/16S/GF_PWY")

meta_ob1 <- read.table("C:/Users/15869/Desktop/CPP_LOCO_GF/16S/GF_PWY/GF_PWY_meta.txt", header = TRUE, row.names = 1)

meta_ob1 <- meta_ob1[meta_ob1$Study== "GF",]
meta_ob1 <- meta_ob1[meta_ob1$Time== "End",]

nrow(meta_ob1)

library(phyloseq)
library(metagMisc)
library(vegan)
library(ggplot2)

tax_ob1 <- import_mothur(mothur_constaxonomy = ("C:/Users/15869/Desktop/CPP_LOCO_GF/16S/GF_PWY/GF_PWY_taxa.txt"))

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("C:/Users/15869/Desktop/CPP_LOCO_GF/16S/GF_PWY/GF_PWY_shared.csv",)
#fullMothur<-read.csv("mini_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtu))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
AF_mothur2[] <- lapply(AF_mothur2, function(x) as.numeric(as.character(x)))
sum(is.na(AF_mothur2))  # Check if there are any NA values
otu_tab <- otu_table(AF_mothur2, taxa_are_rows = FALSE)

phyloseqobj <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))


merged_phy_percent1 <- tax_glom(phyloseqobj , taxrank="Rank1")

# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent1) {
  sd <- sample_data(merged_phy_percent1)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent1) {
  OTU <- otu_table(merged_phy_percent1)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


taxa_names(merged_phy_percent1)
tax_table(merged_phy_percent1)[,1]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,1]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
nrow(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)
row.names(input_metadata)

library(Maaslin2)

fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=33,
                    min_prevalence=0.34,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "maaslin2_GF_Sex_Group_PWY",
                    fixed_effects  = c("Group", "Sex"),
                    reference      = c("Group,Conv", "Sex,F"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)



