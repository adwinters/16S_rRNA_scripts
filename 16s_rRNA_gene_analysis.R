








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


###Rank7
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




















###Maaslin2


###Rank6
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









###PWY plots

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



# Read the data from the file

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










































setwd("C:/Users/15869/Desktop/CPP_LOCO_GF/16S/GF_PWY/")



##PWY


meta_ob1 <- read.table("C:/Users/15869/Desktop/CPP_LOCO_GF/16S/GF_PWY/GF_PWY_meta.txt", header = TRUE, row.names = 1)

nrow(meta_ob1)

# Load the dplyr library
library(dplyr)
# Group by Week, Type, and Treatment, then count the samples
result <- meta_ob1 %>%
  group_by(Study, Sex, Group) %>%
  summarise(count = n())

# Print the result
print(result, n=100)

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

merged_phy_percent1  <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))

# Remove taxa with zero counts across all samples
merged_phy_percent1 <- prune_taxa(taxa_sums(merged_phy_percent1) > 0, merged_phy_percent1)



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

meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
nrow(shared_phy_out1)
ncol(shared_phy_out1)




library(eulerr)

custom_labels <- c("Conv" = "CONV", "GF" = "GF")

# Define the function with simplified "single" and "ninety" criteria
generate_venn <- function(data, metadata, category, criteria = c("single", "ninety"), fill_colors = c("#FED2D2", "#88cafc", "#FA0000", "#005DE6")) {
  # Extract the specified group column
  groups <- metadata[[category]]
  
  # Initialize lists for OTUs
  otu_list <- list()
  
  # Loop through unique groups and extract OTUs based on the specified criteria
  for (group in unique(groups)) {
    # Get the sample indices for the current group
    sample_indices <- which(groups == group)
    
    # Determine presence of pathways based on the criteria
    otus_present <- colSums(data[sample_indices, ] > 0)
    
    if (criteria == "single") {
      # Include pathways that have any presence (>0) in the group's samples
      otu_list[[group]] <- names(otus_present[otus_present > 0])
    } else if (criteria == "ninety") {
      # Include pathways present in at least 90% of the group's samples
      otu_list[[group]] <- names(otus_present[otus_present >= (length(sample_indices) * 0.9)])
    } else {
      stop("Invalid criteria specified. Choose 'single' or 'ninety'.")
    }
  }
  
  # Prepare data for eulerr
  euler_data <- euler(otu_list)
  
  # Generate the Euler diagram with proportional areas and labels
  venn_plot <- plot(euler_data,
                    fills = list(fill = fill_colors, alpha = 0.75),
                    quantities = list(cex = 0.85),
                    edges = "black",
                    legend = list(side = "bottom", nrow = 1, ncol = 2, cex = 0.85, labels = custom_labels))
  
  # Convert to ggplot object
  venn_ggplot <- as_ggplot(venn_plot)
  
  return(venn_ggplot)
}

# Example usage
GF_venn_plot <- generate_venn(shared_phy_out1, meta_phy_out1, "Group", "single", fill_colors = c("#005DE6", "#FA0000"))

# Adjust plot margins
GF_venn_without_Sex_ninety <- GF_venn_plot + theme(plot.margin = unit(c(0.5, 0.0, 0.25, 0.0), "cm")) # top, right, bottom, and left

GF_venn_without_Sex_ninety

























































# Load necessary libraries
library(phyloseq)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(RColorBrewer)

#The problem is here. After tax_glom by Rank2, the tax_names need to be renamed by Rank2.

# Prune to the top 'num_taxa' taxa based on abundance
merged_phy_percent <- phyloseq_standardize_otu_abundance(phyloseqobj.f2, method = "total")
merged_phy_percent1 <- tax_glom(merged_phy_percent, taxrank = "Rank2")
merged_phy_percent2 <- prune_taxa(names(sort(taxa_sums(merged_phy_percent1), TRUE)[1:5]), merged_phy_percent1)

taxa_names(merged_phy_percent2)



# Extract metadata from the phyloseq object
sample_data_df <- as.data.frame(sample_data(merged_phy_percent2))

# Generate the heatmap_data from the phyloseq object and transpose it
phyloseq_relative <- transform_sample_counts(merged_phy_percent2, function(x) x / sum(x))
otu_table_matrix <- as.data.frame(otu_table(phyloseq_relative))
heatmap_data <- as.matrix(otu_table_matrix)

# Transpose the heatmap data
heatmap_data <- t(heatmap_data)

# Check unique levels of Sex_Group
unique_sex_groups <- unique(sample_data_df$Sex_Group)
cat("Unique levels in Sex_Group:", unique_sex_groups, "\n")

# Create annotation for Sex_Group with updated colors
sex_group_colors <- c("GF_F" = "#FA0000", "GF_M" = "#005DE6", "Conv_F" = "#FED2D2", "Conv_M" = "#88cafc")

color <- colorRampPalette(brewer.pal(9, "Purples"))(100)

sex_group_annotation <- HeatmapAnnotation(
  Sex_Group = sample_data_df$Sex_Group,
  col = list(Sex_Group = sex_group_colors),
  which = "column"
)


bray_dist <- function(x) {
  vegdist(log(x+1), method = "bray")
}


# Generate the heatmap with annotations
heatmap_plot <- Heatmap(
  heatmap_data,
  name = "Relative Abundance",  # Title for the heatmap legend
  col = color,  # Color scale for relative abundance
  cluster_rows = TRUE,  # Enable row clustering
  clustering_distance_rows = bray_dist,  # Specify distance metric for rows
  clustering_method_rows = "ward.D2",  # Clustering method for rows
  cluster_columns = FALSE,  # Disable column clustering
  show_row_names = TRUE,  # Show row names
  show_column_names = FALSE,  # Show column names
  show_row_dend = FALSE,
  top_annotation = sex_group_annotation,  # Add sex group annotations
  heatmap_legend_param = list(  # Parameters for the heatmap legend
    title = "Relative Abundance",  # Title for the relative abundance legend
    direction = "horizontal",  # Set legend to horizontal
    legend_width = unit(6, "cm")  # Set width of the legend
  )
)

# Draw the heatmap with specified legend positioning
draw(
  heatmap_plot,
  heatmap_legend_side = "bottom",  # Position the heatmap legend at the bottom
  annotation_legend_side = "right",  # Position the annotations (including Sex_Group) on the right
  legend_grouping = "original"  # Maintain original grouping for the legend
)









# Convert row names to a column named "Taxa"
heatmap_df$Taxa <- rownames(heatmap_df)

# Now melt the dataframe
heatmap_long <- melt(heatmap_df, id.vars = "Taxa", variable.name = "Sample", value.name = "Relative_Abundance")

























# Load necessary libraries
library(ggalign)
library(reshape2)
library(viridis)  # For color scaling

# Assuming heatmap_long is already defined
# You previously melted your heatmap_df to long format:
# heatmap_long <- melt(heatmap_df, id.vars = "Taxa", variable.name = "Sample", value.name = "Relative_Abundance")

# Convert heatmap_long back to a wide format matrix for ggheatmap
heatmap_matrix <- acast(heatmap_long, Taxa ~ Sample, value.var = "Relative_Abundance")

# Create the heatmap using ggheatmap
heatmap_plot <- ggheatmap(heatmap_matrix) +
  scale_fill_viridis_c() +  # Set color scale
  hmanno("top") +  # Add top annotation
  align_dendro(aes(color = branch), k = 3) +  # Split into 3 groups (adjust k as needed)
  geom_point(aes(color = branch, y = y)) +  # Add points to the dendrogram
  scale_color_brewer(palette = "Dark2") +  # Set color palette for dendrogram
  labs(title = "Relative Abundance Heatmap", x = "Samples", y = "Taxa") +  # Add labels
  theme_minimal()  # Use a minimal theme

# Draw the aligned heatmap
print(heatmap_plot)

















#Heatmap


setwd("C:/Users/15869/Desktop/CPP_LOCO_GF/16S/")


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

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))


###Rank6
merged_phy_percent <- phyloseq_standardize_otu_abundance(phyloseqobj.f, method = "total")
merged_phy_percent <- tax_glom(merged_phy_percent, taxrank="Rank7")
merged_phy_percent1 <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:25]), merged_phy_percent)


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

row.names(input_data)
row.names(input_metadata)


F_data <- input_data
F_meta <- input_metadata



F_data_selected <- F_data[rownames(F_data) %in% rownames(F_meta), ]
F_data_selected <- F_data_selected[rownames(F_meta), ]
F_df_mat <- data.matrix(F_data_selected)

annotation_F_data <- F_meta[, c("Group", "Sex")]

library(vegan)
drows <- vegdist(F_df_mat, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)

library(pheatmap2)
# Adjusted color palette
color <- colorRampPalette(brewer.pal(9, "Purples"))(100) # Fixed palette size issue

library(ggplotify)
library(pheatmap)
library(RColorBrewer)

# Ensure the levels in annotation data match the colors
annotation_F_data$Group <- factor(annotation_F_data$Group, levels = c("Conv", "GF"))

# Update annotation_colors to match the actual levels
annotation_colors <- list(Group = c(Conv = "green", GF = "orange"),
                          Sex = c(F = "#5B84B1", M = "#FC766A"))

pheatmap2(t(F_df_mat),
          color = color,
          treeheight_row = 0,
          treeheight_col = 30,
          cluster_cols = FALSE,
          clustering_distance_rows = "bray",
          clustering_distance_cols = drows,
          annotation_col = annotation_F_data,
          annotation_colors = annotation_colors,
          fontsize_row = 10,
          fontsize_col = 10,
          angle_col = 45,
          gaps_col = c(12),
          border_color = NA)





































###Merged CPPand LOCO

setwd("C:/Users/15869/Desktop/CPP_LOCO_GF/16S/")


meta_ob1 <- read.table("C:/Users/15869/Desktop/CPP_LOCO_GF/16S/CPP_LOCO_GF_meta.txt", header = TRUE, row.names = 1)

meta_ob1 <- meta_ob1[meta_ob1$Study!="GF",]
meta_ob1 <- meta_ob1[meta_ob1$Time=="End",]
meta_ob1 <- meta_ob1[meta_ob1$Treatment!="B_SC" & meta_ob1$Treatment!="B_C",]


nrow(meta_ob1)


# Load the dplyr library
library(dplyr)
# Group by Week, Type, and Treatment, then count the samples
result <- meta_ob1 %>%
  group_by(Sex_Group, Group, Treatment, Sex, Time) %>%
  summarise(count = n())

# Print the result
print(result, n=100)

#write.csv (result, file = "CPP_counts.csv")










###Alpha diversity


meta_ob1 <- read.table("C:/Users/15869/Desktop/CPP_LOCO_GF/16S/CPP_LOCO_GF_meta.txt", header = TRUE, row.names = 1)

meta_ob1 <- meta_ob1[meta_ob1$Exclude_GF!="Yes",]
meta_ob1 <- meta_ob1[meta_ob1$Time=="End",]
meta_ob1 <- meta_ob1[meta_ob1$Treatment!="B_SC" & meta_ob1$Treatment!="B_C",]

nrow(meta_ob1)

library(phyloseq)
library(metagMisc)
library(vegan)
library(ggplot2)

tax_ob1 <- import_mothur(mothur_constaxonomy = ("C:/Users/15869/Desktop/CPP_LOCO_GF/16S/CPP_LOCO_GF_taxa.txt"))


GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("CPP_LOCO_GF_shared.csv",)
#fullMothur<-read.csv("mini_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtu))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))




# Extract the taxonomy table
tax_table_df <- as.data.frame(tax_table(phyloseqobj.f))

# Remove taxa at Rank1 that are Archaea
tax_table_df <- tax_table_df[!(tax_table_df$Rank6 == "Remove"), ]

# Update the taxonomy table in the phyloseq object
tax_table(phyloseqobj.f) <- as.matrix(tax_table_df)

phyloseqobj.f2 <-phyloseqobj.f

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


meta_ob2 <- cbind(meta_ob1$Sex, meta_ob1$Group, meta_ob1$Treatment, Coverage_out, richness_out$Chao1,
                  richness_out$Shannon, richness_out$InvSimpson, meta_ob1$Sex_Group)
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
  geom_boxplot(outlier.shape = NA, alpha=0.2, show.legend = FALSE) +
  geom_point(position=position_jitterdodge(jitter.width = 0.333, jitter.height = 0, dodge.width = 0.75, seed = NULL), 
             size=3, alpha=0.4, stroke=1, shape=21, color="black", show.legend = FALSE) +
  scale_shape_manual(values = c(22, 22, 22, 22, 22, 22, 22, 22)) +
  scale_color_manual(values=c("#FC766A", "#5B84B1", "green","orange","blue","black","magenta","cyan","pink","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#FC766A", "#5B84B1", "green","orange","blue","black","magenta","cyan","pink","gray","orange","brown","purple")) +
  theme_classic(base_size = 11) +
  theme(text = element_text(color = "black"),   
        axis.title.x = element_blank(),    
        axis.text.x = element_blank(),     
        axis.text.y = element_text(color = "black"))



chao1_plot <- ggplot(meta_ob2, aes(x=Sex_Group, y=Chao1, color=Sex_Group, fill=Sex_Group)) +
  geom_boxplot(outlier.shape = NA, alpha=0.2, show.legend = FALSE) +
  geom_point(position=position_jitterdodge(jitter.width = 0.333, jitter.height = 0, dodge.width = 0.75, seed = NULL), 
             size=3, alpha=0.4, stroke=1, shape=21, color="black", show.legend = FALSE) +
  scale_shape_manual(values = c(22, 22, 22, 22, 22, 22, 22, 22)) +
  scale_color_manual(values=c("#FC766A", "#5B84B1", "green","orange","blue","black","magenta","cyan","pink","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#FC766A", "#5B84B1", "green","orange","blue","black","magenta","cyan","pink","gray","orange","brown","purple")) +
  theme_classic(base_size = 11) +
  theme(text = element_text(color = "black"),   
        axis.title.x = element_blank(),    
        axis.text.x = element_blank(),     
        axis.text.y = element_text(color = "black"))


shannon_plot <- ggplot(meta_ob2, aes(x=Sex_Group, y=Shannon, color=Sex_Group, fill=Sex_Group)) +
  geom_boxplot(outlier.shape = NA, alpha=0.2, show.legend = FALSE) +
  geom_point(position=position_jitterdodge(jitter.width = 0.333, jitter.height = 0, dodge.width = 0.75, seed = NULL), 
             size=3, alpha=0.4, stroke=1, shape=21, color="black", show.legend = FALSE) +
  scale_shape_manual(values = c(22, 22, 22, 22, 22, 22, 22, 22)) +
  scale_color_manual(values=c("#FC766A", "#5B84B1", "green","orange","blue","black","magenta","cyan","pink","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#FC766A", "#5B84B1", "green","orange","blue","black","magenta","cyan","pink","gray","orange","brown","purple")) +
  theme_classic(base_size = 11) +
  theme(text = element_text(color = "black"),   
        axis.title.x = element_blank(),    
        axis.text.x = element_blank(),     
        axis.text.y = element_text(color = "black"))


invsimpson_plot <- ggplot(meta_ob2, aes(x=Sex_Group, y=InvSimpson, color=Sex_Group, fill=Sex_Group)) +
  geom_boxplot(outlier.shape = NA, alpha=0.2, show.legend = FALSE) +
  geom_point(position=position_jitterdodge(jitter.width = 0.333, jitter.height = 0, dodge.width = 0.75, seed = NULL), 
             size=3, alpha=0.4, stroke=1, shape=21, color="black", show.legend = FALSE) +
  scale_shape_manual(values = c(22, 22, 22, 22, 22, 22, 22, 22)) +
  scale_color_manual(values=c("#FC766A", "#5B84B1", "green","orange","blue","black","magenta","cyan","pink","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#FC766A", "#5B84B1", "green","orange","blue","black","magenta","cyan","pink","gray","orange","brown","purple")) +
  theme_classic(base_size = 11) +
  theme(text = element_text(color = "black"),   
        axis.title.x = element_blank(),    
        axis.text.x = element_blank(),     
        axis.text.y = element_text(color = "black"))


library(ggpubr)

alpha_merged <- ggarrange(chao1_plot, shannon_plot, invsimpson_plot ,
                          common.legend = TRUE, legend = "bottom",
                          ncol = 3, nrow = 1, align="hv")






library(car)

chao_model <- glm(Chao1 ~ Sex*Group, data = meta_ob2)
Anova(chao_model, type=3, test.statistic="F")

shannon_model <- glm(Shannon ~ Sex*Group, data = meta_ob2)
Anova(shannon_model, type=3, test.statistic="F")

invsimpson_model <- glm(InvSimpson ~ Sex*Group, data = meta_ob2)
Anova(invsimpson_model, type=3, test.statistic="F")


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
shannon_tukey_out <- pairs(invsimpson_emm, adjust = "tukey")
#write.csv(shannon_tukey_out , file = "CPP_shannon_tukey.csv")







#beta diversity



###Rank7
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



# Load necessary libraries
library(vegan)
library(ggplot2)
library(glue)
library(dplyr)

# Assuming Day1_meta and Day1_shared are pre-loaded data frames

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
df$Sex_Group <- factor(input_metadata$Sex_Group)
df$Group <- factor(input_metadata$Group)
df$Sex <- factor(input_metadata$Sex)

# Labels for axes
labx1 <- glue("PCo 1 ({prop_explained[1]}%)")
laby1 <- glue("PCo 2 ({prop_explained[2]}%)")


# Create the plot
Bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Sex_Group, fill = Sex_Group, shape = Sex_Group)) +
  scale_color_manual(values=c("#FC766A", "#5B84B1", "green","orange","blue","black","magenta","cyan","pink","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#FC766A", "#5B84B1", "green","orange","blue","black","magenta","cyan","pink","gray","orange","brown","purple")) +
  geom_point(aes(shape = Sex_Group), size = 3, alpha = 0.4, stroke = 1, color = "black") +
  scale_shape_manual(values = c(21, 21, 21, 21, 21, 21, 21, 21)) +  # Use shapes that support fill and stroke
  stat_ellipse(aes(group = Sex_Group), geom = "polygon", linetype = "dotted", level = 0.8, alpha = 0.2, size = 1, show.legend = FALSE) +
  theme_classic() +
  labs(x = labx1, y = laby1) +
  theme(axis.text.x = element_text(color = "black"), 
        text = element_text(size = 11), 
        axis.text.y = element_text(color = "black"),
        legend.title = element_blank(), 
        strip.text.x = element_blank())
        #legend.text = element_text(size = 11))



library(patchwork)

# Create the individual plots with labels
chao1_plot <- chao1_plot + labs(tag = "A")
shannon_plot <- shannon_plot + labs(tag = "B")
invsimpson_plot <- invsimpson_plot + labs(tag = "C")
Bray_plot <- Bray_plot + labs(tag = "D")

# Combine plots and add a label to the entire plot
combined_plot <- (chao1_plot | shannon_plot | invsimpson_plot) / Bray_plot +
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom') &
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 11)) 





library(vegan)
dist.bray <- vegdist(input_data, method="bray", binary=FALSE)
adonis2(formula = dist.bray  ~ Group*Sex, data = input_metadata, permutations=1000)


library(pairwiseAdonis)
x <- phyloseq::distance(merged_phy_percent1, method ="bray")
factors <- meta_phy_out1$Sex_Group
pairwise.adonis(x, factors, perm = 1000, p.adjust(method="none"))





###Maaslin2



###Rank6
merged_phy_percent1 <- tax_glom(phyloseqobj.f2 , taxrank="Rank6")


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
tax_table(merged_phy_percent1)[,6]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,6]


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
                    max_significance=0.05,
                    min_abundance=33,
                    min_prevalence=0.25,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "maaslin2_merged_LOCO_Sex_Group_Rank7",
                    fixed_effects  = c("Group", "Sex"),
                    reference      = c("Group,A_Control", "Sex,F"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)





#Heatmap


# Subset phyloseq object for each group
phy_A_Control <- subset_samples(phyloseqobj.f2, Group == "A_Control")
phy_B_ABX <- subset_samples(phyloseqobj.f2, Group == "B_ABX")

# Convert to relative abundance
phy_A_Control_percent <- phyloseq_standardize_otu_abundance(phy_A_Control, method = "total")
phy_B_ABX_percent <- phyloseq_standardize_otu_abundance(phy_B_ABX, method = "total")

# Aggregate by taxonomy rank
phy_A_Control_tax_glom <- tax_glom(phy_A_Control_percent, taxrank="Rank6")
phy_B_ABX_tax_glom <- tax_glom(phy_B_ABX_percent, taxrank="Rank6")

# Select top 20 ASVs based on abundance in each group
top_20_A_Control <- names(sort(taxa_sums(phy_A_Control_tax_glom), TRUE)[1:10])
top_20_B_ABX <- names(sort(taxa_sums(phy_B_ABX_tax_glom), TRUE)[1:10])

# Prune taxa to include only top 20 ASVs from each group
phy_A_Control_top20 <- prune_taxa(top_20_A_Control, phy_A_Control_tax_glom)
phy_B_ABX_top20 <- prune_taxa(top_20_B_ABX, phy_B_ABX_tax_glom)


# Combine the lists of top ASVs
combined_top_40 <- unique(c(top_20_A_Control, top_20_B_ABX))

# Prune the phyloseq object to include only the combined top ASVs
phy_combined_top_40 <- prune_taxa(combined_top_40, tax_glom(phyloseq_standardize_otu_abundance(phyloseqobj.f2, method = "total"), taxrank = "Rank6"))


# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(phy_combined_top_40) {
  sd <- sample_data(phy_combined_top_40)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(phy_combined_top_40) {
  OTU <- otu_table(phy_combined_top_40)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


taxa_names(phy_combined_top_40)
tax_table(phy_combined_top_40)[,6]
taxa_names(phy_combined_top_40) <- tax_table(phy_combined_top_40)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(phy_combined_top_40))
shared_phy_out1 <- as.data.frame(psotu2veg(phy_combined_top_40))

nrow(meta_phy_out1)
nrow(shared_phy_out1)

input_data  <- as.data.frame(shared_phy_out1)
input_metadata <- as.data.frame(meta_phy_out1)


F_data <- input_data
F_meta <- input_metadata


ncol(F_data)
#F_meta_selected <- F_meta[F_meta$Group=="IA",]
F_data_selected <-  F_data[rownames(F_data)%in%rownames(F_meta),]
#ncol(F_data_selected)

F_data_selected <- F_data_selected[rownames(F_meta),]
F_df_mat = data.matrix(F_data_selected)

annotation_F_data <- F_meta[,c("Group","Sex")]


library(vegan)
drows = vegdist(F_df_mat, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
#dcols = vegdist(F_df_mat, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)


library(pheatmap2)
# blue and green color palette
color <- colorRampPalette((c("white", "blue", "red")))(50)

library(ggplotify)


library(pheatmap)
library(RColorBrewer)

color <- colorRampPalette(brewer.pal(9, "Purples"))(100)

# Ensure the levels in annotation data match the colors
annotation_F_data$Group <- factor(annotation_F_data$Group, levels = c("A_Control", "B_ABX"))

# Update annotation_colors to match the actual levels
annotation_colors <- list(Group = c(A_Control = "green", B_ABX = "orange"),
                          Sex = c(F = "#FC766A", M = "#5B84B1"))

pheatmap2(t(F_df_mat),
          color = color,
          treeheight_row = 0,
          treeheight_col = 30,
          cluster_cols = F,
          clustering_distance_rows = "bray",
          clustering_distance_cols = dcols,
          annotation_col = annotation_F_data,
          annotation_colors = annotation_colors,
          fontsize_row = 10,
          fontsize_col = 10,
          angle_col = 45,
          #gaps_col = c(8, 16, 24),
          border_color = NA)
          
      
          
          
          
  
        
        
          
          
          #####LOCO
          
          # Read the data from the file
          
          meta_ob1 <- read.table("C:/Users/15869/Desktop/CPP_LOCO_GF/16S/CPP_LOCO_GF_meta.txt", header = TRUE, row.names = 1)
          
          meta_ob1 <- meta_ob1[meta_ob1$Study== "LOCO",]
          meta_ob1 <- meta_ob1[meta_ob1$Time== "End",]
          meta_ob1 <- meta_ob1[meta_ob1$Treatment== "A_S",]
          meta_ob1 <- meta_ob1[meta_ob1$Time!= "Other" & meta_ob1$Time!= "Blank" & meta_ob1$Time!= "Mock",]
          
          unique(meta_ob1$Time)
          unique(meta_ob1$Group)
          
          nrow(meta_ob1)
          
          
          # Load the dplyr library
          library(dplyr)
          # Group by Week, Type, and Treatment, then count the samples
          result <- meta_ob1 %>%
            group_by(Group, Treatment, Sex, Time) %>%
            summarise(count = n())
          
          # Print the result
          print(result, n=100)
          
          #write.csv (result, file = "CPP_counts.csv")
          
          
          
          
          
          ###Alpha diversity
          
          
          meta_ob1 <- read.table("C:/Users/15869/Desktop/CPP_LOCO_GF/16S/CPP_LOCO_GF_meta.txt", header = TRUE, row.names = 1)
          
          meta_ob1 <- meta_ob1[meta_ob1$Study== "LOCO",]
          meta_ob1 <- meta_ob1[meta_ob1$Time== "End",]
          meta_ob1 <- meta_ob1[meta_ob1$Treatment== "A_S",]
          meta_ob1 <- meta_ob1[meta_ob1$Time!= "Other" & meta_ob1$Time!= "Blank" & meta_ob1$Time!= "Mock",]
          
          
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
          
          #Make otu table and phyloseq object
          otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
          phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
          
          
          
          
          # Extract the taxonomy table
          tax_table_df <- as.data.frame(tax_table(phyloseqobj.f))
          
          # Remove taxa at Rank1 that are Archaea
          tax_table_df <- tax_table_df[!(tax_table_df$Rank6 == "Remove"), ]
          
          # Update the taxonomy table in the phyloseq object
          tax_table(phyloseqobj.f) <- as.matrix(tax_table_df)
          
          phyloseqobj.f2 <-phyloseqobj.f
          
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
            geom_boxplot(outlier.shape = NA, alpha=0.2, show.legend = FALSE) +
            geom_point(position=position_jitterdodge(jitter.width = 0.333, jitter.height = 0, dodge.width = 0.75, seed = NULL), 
                       size=3, alpha=0.4, stroke=1, shape=21, color="black", show.legend = FALSE) +
            scale_shape_manual(values = c(22, 22)) +
            scale_color_manual(values=c("#FC766A", "#5B84B1", "green","orange","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
            scale_fill_manual(values=c("#FC766A", "#5B84B1", "green","orange","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
            theme_classic(base_size = 11) +
            theme(text = element_text(color = "black"),   
                  axis.title.x = element_blank(),    
                  axis.text.x = element_blank(),     
                  axis.text.y = element_text(color = "black"))
          
          chao1_plot <- ggplot(meta_ob2, aes(x=Sex_Group, y=Chao1, color=Sex_Group, fill=Sex_Group)) +
            geom_boxplot(outlier.shape = NA, alpha=0.2, show.legend = FALSE) +
            geom_point(position=position_jitterdodge(jitter.width = 0.333, jitter.height = 0, dodge.width = 0.75, seed = NULL), 
                       size=3, alpha=0.4, stroke=1, shape=21, color="black", show.legend = FALSE) +
            scale_shape_manual(values = c(22, 22)) +
            scale_color_manual(values=c("#FC766A", "#5B84B1", "green","orange","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
            scale_fill_manual(values=c("#FC766A", "#5B84B1", "green","orange","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
            theme_classic(base_size = 11) +
            theme(text = element_text(color = "black"),   
                  axis.title.x = element_blank(),    
                  axis.text.x = element_blank(),     
                  axis.text.y = element_text(color = "black"))
          
          shannon_plot <- ggplot(meta_ob2, aes(x=Sex_Group, y=Shannon, color=Sex_Group, fill=Sex_Group)) +
            geom_boxplot(outlier.shape = NA, alpha=0.2, show.legend = FALSE) +
            geom_point(position=position_jitterdodge(jitter.width = 0.333, jitter.height = 0, dodge.width = 0.75, seed = NULL), 
                       size=3, alpha=0.4, stroke=1, shape=21, color="black", show.legend = FALSE) +
            scale_shape_manual(values = c(22, 22)) +
            scale_color_manual(values=c("#FC766A", "#5B84B1", "green","orange","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
            scale_fill_manual(values=c("#FC766A", "#5B84B1", "green","orange","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
            theme_classic(base_size = 11) +
            theme(text = element_text(color = "black"),   
                  axis.title.x = element_blank(),    
                  axis.text.x = element_blank(),     
                  axis.text.y = element_text(color = "black"))
          
          invsimpson_plot <- ggplot(meta_ob2, aes(x=Sex_Group, y=InvSimpson, color=Sex_Group, fill=Sex_Group)) +
            geom_boxplot(outlier.shape = NA, alpha=0.2, show.legend = FALSE) +
            geom_point(position=position_jitterdodge(jitter.width = 0.333, jitter.height = 0, dodge.width = 0.75, seed = NULL), 
                       size=3, alpha=0.4, stroke=1, shape=21, color="black", show.legend = FALSE) +
            scale_shape_manual(values = c(22, 22)) +
            scale_color_manual(values=c("#FC766A", "#5B84B1", "green","orange","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
            scale_fill_manual(values=c("#FC766A", "#5B84B1", "green","orange","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
            theme_classic(base_size = 11) +
            theme(text = element_text(color = "black"),   
                  axis.title.x = element_blank(),    
                  axis.text.x = element_blank(),     
                  axis.text.y = element_text(color = "black"))
          
          
          library(ggpubr)
          
          alpha_merged <- ggarrange(chao1_plot, shannon_plot, invsimpson_plot ,
                                    common.legend = TRUE, legend = "bottom",
                                    ncol = 1, nrow = 3, align="hv")
          
          
          
          
          
          
          
          library(car)
          
          chao_model <- glm(Chao1 ~ Sex*Group, data = meta_ob2)
          Anova(chao_model, type=3, test.statistic="F")
          
          shannon_model <- glm(Shannon ~ Sex*Group, data = meta_ob2)
          Anova(shannon_model, type=3, test.statistic="F")
          
          invsimpson_model <- glm(InvSimpson ~ Sex*Group, data = meta_ob2)
          Anova(invsimpson_model, type=3, test.statistic="F")
          
          
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
          shannon_tukey_out <- pairs(invsimpson_emm, adjust = "tukey")
          #write.csv(shannon_tukey_out , file = "CPP_shannon_tukey.csv")
          
          
          
          
          
          
          #beta diversity
          
          
          
          ###Rank7
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
          
          # Assuming Day1_meta and Day1_shared are pre-loaded data frames
          
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
          df$Group <- factor(input_metadata$Group)
          df$Sex <- factor(input_metadata$Sex)
          
          # Labels for axes
          labx1 <- glue("PCo 1 ({prop_explained[1]}%)")
          laby1 <- glue("PCo 2 ({prop_explained[2]}%)")
          
          # Create the plot
          Bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Sex_Group, fill = Sex_Group, shape = Sex_Group)) +
            scale_color_manual(values=c("#FC766A", "#5B84B1", "green","orange","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
            scale_fill_manual(values=c("#FC766A", "#5B84B1", "green","orange","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
            geom_point(aes(shape = Sex_Group), size = 3, alpha = 0.4, stroke = 1, color = "black") +
            scale_shape_manual(values = c(21, 21, 21, 21)) +  # Use shapes that support fill and stroke
            stat_ellipse(aes(group = Sex_Group), geom = "polygon", linetype = "dotted", level = 0.8, alpha = 0.2, size = 1, show.legend = FALSE) +
            theme_classic() +
            labs(x = labx1, y = laby1) +
            theme(axis.text.x = element_text(color = "black"), 
                  text = element_text(size = 11), 
                  axis.text.y = element_text(color = "black"),
                  legend.title = element_blank(), 
                  strip.text.x = element_blank())
          #legend.text = element_text(size = 11))
          
          
          
          
          library(patchwork)
          
          # Create the individual plots with labels
          chao1_plot <- chao1_plot + labs(tag = "A")
          shannon_plot <- shannon_plot + labs(tag = "B")
          invsimpson_plot <- invsimpson_plot + labs(tag = "C")
          Bray_plot <- Bray_plot + labs(tag = "D")
          
          # Combine plots and add a label to the entire plot
          combined_plot <- (chao1_plot | shannon_plot | invsimpson_plot) / Bray_plot +
            plot_layout(guides = 'collect') & 
            theme(legend.position = 'bottom') &
            plot_annotation(tag_levels = 'A') &
            theme(plot.tag = element_text(size = 11))
          
          
          
          
          
          
          library(vegan)
          dist.bray <- vegdist(input_data, method="bray", binary=FALSE)
          adonis2(formula = dist.bray  ~ Group*Sex, data = input_metadata, permutations=1000)
          
          
          library(pairwiseAdonis)
          
          x <- phyloseq::distance(merged_phy_percent1, method ="bray")
          factors <- meta_phy_out1$Sex_Group
          pairwise.adonis(x, factors, perm = 1000, p.adjust(method="none"))
          
          
          
          
          
        
          
          ###Maaslin2
          
          
          ###Rank6
          merged_phy_percent1 <- tax_glom(phyloseqobj.f2 , taxrank="Rank7")
          
          
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
          
          input_data  <- as.data.frame(t(shared_phy_out1))
          input_metadata <- as.data.frame(meta_phy_out1)
          
          row.names(input_data)
          row.names(input_metadata)
          
          
          unique(input_metadata$Study)
          unique(input_metadata$Sex)
          unique(input_metadata$Group)
          
          fit_data = Maaslin2(input_data     = input_data, 
                              input_metadata = input_metadata,
                              max_significance=0.05,
                              min_abundance=33,
                              min_prevalence=0.25,
                              #min_variance=20000,
                              normalization  = "TMM",
                              output         = "maaslin2_merged_LOCO_Sex_Group_Rank7",
                              fixed_effects  = c("Study", "Group", "Sex"),
                              reference      = c("Study,CPP", "Group,A_Control", "Sex,F"),
                              #random_effects  = "Mouse_ID",
                              analysis_method="NEGBIN",
                              transform = "NONE",
                              cores=12)
          
          
          
          
          
          
          
          
          #Heatmap
          
          
          # Subset phyloseq object for each group
          phy_A_Control <- subset_samples(phyloseqobj.f2, Group == "A_Control")
          phy_B_ABX <- subset_samples(phyloseqobj.f2, Group == "B_ABX")
          
          # Convert to relative abundance
          phy_A_Control_percent <- phyloseq_standardize_otu_abundance(phy_A_Control, method = "total")
          phy_B_ABX_percent <- phyloseq_standardize_otu_abundance(phy_B_ABX, method = "total")
          
          # Aggregate by taxonomy rank
          phy_A_Control_tax_glom <- tax_glom(phy_A_Control_percent, taxrank="Rank6")
          phy_B_ABX_tax_glom <- tax_glom(phy_B_ABX_percent, taxrank="Rank6")
          
          # Select top 20 ASVs based on abundance in each group
          top_20_A_Control <- names(sort(taxa_sums(phy_A_Control_tax_glom), TRUE)[1:10])
          top_20_B_ABX <- names(sort(taxa_sums(phy_B_ABX_tax_glom), TRUE)[1:10])
          
          # Prune taxa to include only top 20 ASVs from each group
          phy_A_Control_top20 <- prune_taxa(top_20_A_Control, phy_A_Control_tax_glom)
          phy_B_ABX_top20 <- prune_taxa(top_20_B_ABX, phy_B_ABX_tax_glom)
          
         
          # Combine the lists of top ASVs
          combined_top_40 <- unique(c(top_20_A_Control, top_20_B_ABX))
          
          # Prune the phyloseq object to include only the combined top ASVs
          phy_combined_top_40 <- prune_taxa(combined_top_40, tax_glom(phyloseq_standardize_otu_abundance(phyloseqobj.f2, method = "total"), taxrank = "Rank6"))
          
          
          # convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
          pssd2veg <- function(phy_combined_top_40) {
            sd <- sample_data(phy_combined_top_40)
            return(as(sd,"data.frame"))
          }
          
          # convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
          psotu2veg <- function(phy_combined_top_40) {
            OTU <- otu_table(phy_combined_top_40)
            if (taxa_are_rows(OTU)) {
              OTU <- t(OTU)
            }
            return(as(OTU, "matrix"))
          }
          
          
          taxa_names(phy_combined_top_40)
          tax_table(phy_combined_top_40)[,6]
          taxa_names(phy_combined_top_40) <- tax_table(phy_combined_top_40)[,6]
          
          
          meta_phy_out1 <- as.data.frame(pssd2veg(phy_combined_top_40))
          shared_phy_out1 <- as.data.frame(psotu2veg(phy_combined_top_40))
          
          nrow(meta_phy_out1)
          nrow(shared_phy_out1)
          
          input_data  <- as.data.frame(shared_phy_out1)
          input_metadata <- as.data.frame(meta_phy_out1)
          
          
          F_data <- input_data
          F_meta <- input_metadata
          
          
          ncol(F_data)
          #F_meta_selected <- F_meta[F_meta$Group=="IA",]
          F_data_selected <-  F_data[rownames(F_data)%in%rownames(F_meta),]
          #ncol(F_data_selected)
          
          F_data_selected <- F_data_selected[rownames(F_meta),]
          F_df_mat = data.matrix(F_data_selected)
          
          annotation_F_data <- F_meta[,c("Group","Sex")]
          
          
          library(vegan)
          drows = vegdist(F_df_mat, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
          #dcols = vegdist(F_df_mat, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
          
          
          library(pheatmap2)
          # blue and green color palette
          color <- colorRampPalette((c("white", "blue", "red")))(50)
          
          library(ggplotify)
          
          
          library(pheatmap)
          library(RColorBrewer)
          
          color <- colorRampPalette(brewer.pal(9, "Purples"))(100)
          
          # Ensure the levels in annotation data match the colors
          annotation_F_data$Group <- factor(annotation_F_data$Group, levels = c("A_Control", "B_ABX"))
          
          # Update annotation_colors to match the actual levels
          annotation_colors <- list(Group = c(A_Control = "green", B_ABX = "orange"),
                                    Sex = c(F = "#FC766A", M = "#5B84B1"))
          
          pheatmap2(t(F_df_mat),
                    color = color,
                    treeheight_row = 0,
                    treeheight_col = 30,
                    cluster_cols = F,
                    clustering_distance_rows = "bray",
                    clustering_distance_cols = dcols,
                    annotation_col = annotation_F_data,
                    annotation_colors = annotation_colors,
                    fontsize_row = 10,
                    fontsize_col = 10,
                    angle_col = 45,
                    gaps_col = c(8, 16, 24),
                    border_color = NA)
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
         
          
          
          
          
          
          
          
          
          
          
          
          
          
          library(ggplot2)
          library(reshape2)
          
          # Convert matrix to long format
          F_df_mat_long <- melt(t(F_df_mat))
          colnames(F_df_mat_long) <- c("Sample", "Variable", "Value")
          
          # Print the number of rows and columns
          print(paste("Number of rows in annotation_F_data:", nrow(annotation_F_data)))
          print(paste("Number of columns in F_df_mat:", length(colnames(F_df_mat))))
          
          # Ensure that annotation_F_data has the correct number of rows
          if (nrow(annotation_F_data) != length(colnames(F_df_mat))) {
            # If there are extra rows, subset annotation_F_data to match
            annotation_F_data <- annotation_F_data[1:length(colnames(F_df_mat)), , drop = FALSE]
            print("Subsetting annotation_F_data to match the columns of F_df_mat.")
            
            # Verify the subsetting
            print("New row names in annotation_F_data:")
            print(rownames(annotation_F_data))
          } else {
            print("The dimensions of annotation_F_data match the columns of F_df_mat.")
          }
          
          # Set row names of annotation_F_data to match the sample names from F_df_mat
          rownames(annotation_F_data) <- colnames(F_df_mat)
          
          # Check for any mismatches
          if (!all(colnames(F_df_mat) %in% rownames(annotation_F_data))) {
            missing_in_annotation <- setdiff(colnames(F_df_mat), rownames(annotation_F_data))
            missing_in_F_df_mat <- setdiff(rownames(annotation_F_data), colnames(F_df_mat))
            
            print("Samples in F_df_mat but not in annotation_F_data:")
            print(missing_in_annotation)
            
            print("Samples in annotation_F_data but not in F_df_mat:")
            print(missing_in_F_df_mat)
          } else {
            # Proceed with subsetting and merging if no mismatches
            annotation_F_data_subset <- annotation_F_data[colnames(F_df_mat), , drop = FALSE]
            annotation_F_data_for_merge <- data.frame(Sample = rownames(annotation_F_data_subset),
                                                      Group = annotation_F_data_subset$Group,
                                                      Sex = annotation_F_data_subset$Sex,
                                                      stringsAsFactors = FALSE)
            
            # Merge annotations with the long-format data
            F_df_mat_long <- merge(F_df_mat_long, annotation_F_data_for_merge, by = "Sample")
            
            # Print the merged data frame
            print("Merged data frame:")
            print(head(F_df_mat_long))
            
            # Create and print the heatmap with ggplot2
            heatmap_plot <- ggplot(F_df_mat_long, aes(x = Variable, y = Sample, fill = Value)) +
              geom_tile() +
              scale_fill_gradient(low = "white", high = "blue") +
              theme_minimal() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              labs(x = "Variable", y = "Sample", fill = "Value") +
              facet_grid(Group ~ Sex, scales = "free")  # Facet by Group and Sex
            
            # Add column annotations if needed
            heatmap_plot <- heatmap_plot +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5))  # Rotate column labels for better visibility
            
            print(heatmap_plot)
          }
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          


ncol(F_data)


F_data_selected <-  F_data[rownames(F_data)%in%rownames(F_meta),]
ncol(F_data_selected)

F_data_selected <- F_data_selected[rownames(F_meta),]

F_df_mat = data.matrix(F_data_selected)


# Define custom colors for each annotation
time_colors <- c("A_Start" = "gold", "End" = "darkorange")
treatment_colors <- c("A_Saline" = "skyblue", "Cocaine" = "darkblue")
group_colors <- c("A_Control" = "forestgreen", "ABX" = "purple")
sex_colors <- c("F" = "deeppink", "M" = "deepskyblue")

# Create the annotation colors list
annotation_colors <- list(
  Time = time_colors,
  Treatment = treatment_colors,
  Group = group_colors,
  Sex = sex_colors
)

# Load necessary libraries
library(vegan)
library(pheatmap2)
library(ggplotify)

# Compute distance matrix
drows <- vegdist(F_df_mat, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)

# Define the color palette for the heatmap
color <- colorRampPalette(c("white", "blue", "red"))(50)

# Generate the heatmap
pheatmap2(t(F_df_mat),
          color = color,
          treeheight_row = 0,
          cluster_cols = FALSE,
          annotation_col = annotation_F_data,
          annotation_colors = annotation_colors,  # Add custom annotation colors
          border_color = NA,
          clustering_distance_rows = "bray",
          show_colnames = FALSE)









###Heatmap




F_data <- input_data
F_meta <- input_metadata




ncol(F_data)
#F_meta_selected <- F_meta[F_meta$Group=="IA",]
F_data_selected <-  F_data[rownames(F_data)%in%rownames(F_meta),]
ncol(F_data_selected)

F_data_selected <- F_data_selected[rownames(F_meta),]

F_df_mat = data.matrix(F_data_selected)



#sample_Group <- as.data.frame(F_meta$Group)
#sample_Day <- as.data.frame(F_meta$DayPlot)
annotation_F_data <- F_meta[,c("Group","DayPlot")]
#row.names(annotation_F_data) <- colnames(F_data)






library(vegan)
drows = vegdist(F_df_mat, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
#dcols = vegdist(F_df_mat, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)


library(pheatmap2)
# blue and green color palette
color <- colorRampPalette((c("white", "blue", "red")))(50)

library(ggplotify)



pheatmap2(t(F_df_mat),
          color = color,
          treeheight_row = 0,
          cluster_cols = F,
          clustering_distance_rows = "bray",
          clustering_distance_cols = dcols,
          #scale = "row",
          annotation_col = annotation_F_data)

getwd()

setwd("/Users/15869/Desktop/Kuhn_Lab/CEF_16S/Run2/Merged_Sex_Results/heatmap/")
hm_plot_data <- as.data.frame(cbind(input_metadata$Group,input_data))

write.csv(hm_plot_data, file="hm_plot_data.csv")



























# Ensure Group is a factor with the correct levels
F_meta$Group <- factor(F_meta$Group, levels = c("A_Saline", "ABX"))
group_colors <- c("A_Control" = "#1F77B4", "ABX" = "#EA660C")

# Ensure Treatment is a factor with the correct levels
F_meta$Group <- factor(F_meta$Treatment, levels = c("A_Saline", "Cocaine"))
treatment_colors <- c("A_Saline" = "#1F77B4", "Cocaine" = "red")

# Ensure DayPlot is a factor with the correct levels
unique_days <- unique(F_meta$Time)
time_colors <- c("A_Start" = "#D89000", "End" = "#A3A500")

# Ensure Sex is a factor with the correct levels
F_meta$Sex <- factor(F_meta$Sex, levels = c("M", "F"))
sex_colors <- c("Female" = "#5B84B1", "Male" = "#FC766A")


# Combine them into a list
annotation_colors <- list(
  Group = group_colors,
  Treatment = treatment_colors ,
  Time = time_colors,
  Sex = sex_colors
)

# Print annotation_colors to verify
print(annotation_colors)

# Sample Group and Day (already in your code)
annotation_F_data <- F_meta[, c("Time", "Treatment", "Group", "Sex")]

# Library imports (already in your code)
library(vegan)
drows = vegdist(F_data, method = "bray", binary = FALSE, diag = FALSE, upper = FALSE, na.rm = FALSE)

library(pheatmap2)
color <- colorRampPalette((c("white", "blue", "red")))(50)



pheatmap2(t(F_df_mat),
          color = color,
          annotation_col = annotation_F_data,
          cluster_cols = F,)


pheatmap2(t(F_df_mat),
          color = color,
          treeheight_row = 0,
          cluster_cols = F,
          annotation_col = annotation_F_data,
          border_color = "NA",
          clustering_distance_rows = "bray",
          show_colnames = F,
          annotation_colors = annotation_colors)







# Generate the heatmap with the custom annotation colors
Female_heatmap <- pheatmap2(F_df_mat, 
                            cluster_cols = FALSE, 
                            clustering_method = "ward.D2", 
                            clustering_distance_rows = "bray", 
                            treeheight_row = 0, 
                            color = color, 
                            annotation_col = annotation_F_data, 
                            cluster_rows = TRUE, 
                            scale = "column", 
                            fontsize_col = 0.000001, 
                            fontsize_row = 9,
                            annotation_colors = annotation_colors)







































































F_data <- input_data
F_meta <- input_metadata




ncol(F_data)
#F_meta_selected <- F_meta[F_meta$Group=="IA",]
F_data_selected <-  F_data[rownames(F_data)%in%rownames(F_meta),]
ncol(F_data_selected)

F_data_selected <- F_data_selected[rownames(F_meta),]

F_df_mat = data.matrix(F_data_selected)



#sample_Group <- as.data.frame(F_meta$Group)
#sample_Day <- as.data.frame(F_meta$DayPlot)
annotation_F_data <- F_meta[,c("Sex", "Group","Time")]
#row.names(annotation_F_data) <- colnames(F_data)






library(vegan)
drows = vegdist(F_df_mat, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
#dcols = vegdist(F_df_mat, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)


library(pheatmap2)
# blue and green color palette
color <- colorRampPalette((c("white", "blue", "red")))(50)

library(ggplotify)



pheatmap2(t(F_df_mat),
          color = color,
          treeheight_row = 0,
          cluster_cols = F,
          clustering_distance_rows = "bray",
          clustering_distance_cols = dcols,
          #scale = "row",
          annotation_col = annotation_F_data,
          border_color = "NA")

getwd()

setwd("D://ABX_LOCO_2024/Merged_Sex_Results/heatmap/")
hm_plot_data <- as.data.frame(cbind(input_metadata$Group,input_data))

write.csv(hm_plot_data, file="hm_plot_data.csv")







#Female_heatmap <- pheatmap2(t(F_df_mat)), 
clustering_method = "ward.D2", 
clustering_distance_rows = "bray", 
treeheight_row = 0, 
color = color, 
#annotation_col = annotation_F_data,
cluster_cols = F,
scale = "row",
fontsize_row = 8,
fontsize_col = 5)
pheatmap2(t(F_df_mat))





























































ggplot(my_df, aes(x=as.factor(Day), y=mean, colour=Group)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
  geom_line() +
  geom_point() + 
 scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  ylab("Chao1") +
  theme_classic() +
  scale_y_continuous(limits = c(0,800), labels = scales::number_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = 7)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9),
        legend.title= element_blank())

















setwd("D://ABX_LOCO_2024/Merged_Sex_Results/")



library(lme4)
library(emmeans)
options(max.print=1000000)
options(digits = 5)
emm_options(opt.digits = FALSE)


#chao_mod <- lmer(as.numeric(meta_ob1$chao) ~ Group*as.factor(meta_ob1$Day)*meta_ob1$Sex + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=TRUE)
chao_mod <- lmer(as.numeric(meta_ob1$chao) ~ Group*as.factor(meta_ob1$Day) + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=TRUE)
chao_out <- emmeans(chao_mod, list(pairwise ~ Group*as.factor(Day)), adjust = "tukey")
chao_out <- as.data.frame(chao_out[2])
write.csv (chao_out, file = "chao_tukey.csv")
library(car)
Anova(chao_mod, type=3, test.statistic="F")

#shannon_mod <- lmer(as.numeric(meta_ob1$shannon) ~ Group*as.factor(meta_ob1$Day)*meta_ob1$Sex + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=TRUE)
shannon_mod <- lmer(as.numeric(meta_ob1$shannon) ~ Group*as.factor(meta_ob1$Day) + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=TRUE)
shannon_out <- emmeans(shannon_mod, list(pairwise ~ Group*as.factor(Day)), adjust = "tukey")
shannon_out <- as.data.frame(shannon_out[2])
write.csv (shannon_out, file = "shannon_tukey.csv")
library(car)
Anova(shannon_mod, type=3, test.statistic="F")

#invsimpson_mod <- lmer(as.numeric(meta_ob1$invsimpson) ~ Group*as.factor(meta_ob1$Day)*meta_ob1$Sex + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=TRUE)
invsimpson_mod <- lmer(as.numeric(meta_ob1$invsimpson) ~ Group*as.factor(meta_ob1$Day) + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=TRUE)
invsimpson_out <- emmeans(invsimpson_mod, list(pairwise ~ Group*as.factor(Day)), adjust = "tukey")
invsimpson_out <- as.data.frame(invsimpson_out[2])
write.csv (invsimpson_out, file = "invsimpson_tukey.csv")
library(car)
Anova(invsimpson_mod, type=3, test.statistic="F")


library(car)
Anova(mod1, type="II")


meta_ob1$chao <- as.numeric(meta_ob1$chao)




pd <- position_dodge(.1)  # Save the dodge spec because we use it repeatedly




ggplot(data = meta_ob1, mapping = aes(y = chao, x = as.factor(Day) ,group = Group, colour=Group, shape=Group)) +
   geom_point(size=5) 



ggplot(meta_ob1, aes(x=as.factor(Day), y=shannon, group=Group, color=Group)) + 
  geom_line() +
  geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
                position=position_dodge(0.05))











chao1_plot <- ggplot(data = meta_ob1,
       aes(x=as.factor(Day), y=chao, fill = factor(Group))) +
  geom_boxplot(fatten=NULL, outlier.shape = NA, position = position_dodge(width=0.85)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.333, dodge.width = 0.85), aes(fill = Group), pch = 21) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), position = position_dodge(width = 0.83), width = 0.74, size=1) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  #facet_wrap(~Sex, scales = "free_x", ncol = 2, nrow = 1) +
  ylab("Chao1") +
  theme_classic() +
  scale_y_continuous(limits = c(0,800), labels = scales::number_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = 7)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9),
        legend.title= element_blank())


  
shannon_plot <- ggplot(data = meta_ob1,
                     aes(x=as.factor(Day), y=shannon, fill = factor(Group))) +
  geom_boxplot(fatten=NULL, outlier.shape = NA, position = position_dodge(width=0.85)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.333, dodge.width = 0.85), aes(fill = Group), pch = 21) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), position = position_dodge(width = 0.83), width = 0.74, size=1) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  #facet_wrap(~Sex, scales = "free_x", ncol = 2, nrow = 1) +
  ylab("Shannon") +
  theme_classic() +
  scale_y_continuous(limits = c(0,5.5), labels = scales::number_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = 9)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9),
        strip.background = element_blank(),
        legend.title= element_blank(),
        strip.text.x = element_blank())




invsimpson_plot <- ggplot(data = meta_ob1,
                       aes(x=as.factor(Day), y=invsimpson, fill = factor(Group))) +
  geom_boxplot(fatten=NULL, outlier.shape = NA, position = position_dodge(width=0.85)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.333, dodge.width = 0.85), aes(fill = Group), pch = 21) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), position = position_dodge(width = 0.83), width = 0.74, size=1) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  #facet_wrap(~Sex, scales = "free_x", ncol = 2, nrow = 1) +
  ylab("Inverse Simpson") +
  theme_classic() +
  scale_y_continuous(limits = c(0,80), labels = scales::number_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = 7)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9),
        strip.background = element_blank(),
        legend.title= element_blank(),
        strip.text.x = element_blank())




library(ggpubr)

alpha_merged <- ggarrange(chao1_plot, shannon_plot, invsimpson_plot ,
                          common.legend = TRUE, legend = "right",
                          ncol = 1, nrow = 3, align="hv")



library(lme4)
library(emmeans)
options(max.print=1000000)
options(digits = 5)
emm_options(opt.digits = FALSE)



#chao
chao_mod <- lmer(chao ~ Group*as.factor(meta_ob1$Day)*Sex + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=FALSE)
chao_out <- emmeans(chao_mod, list(pairwise ~ Group*as.factor(Day)*Sex), adjust = "tukey")
chao_out <- as.data.frame(chao_out[2])
write.csv (chao_out, file = "chao_tukey.csv")

#shannon
shannon_mod <- lmer(shannon ~ Group*as.factor(meta_ob1$Day)*Sex + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=FALSE)
shannon_out <- emmeans(shannon_mod, list(pairwise ~ Group*as.factor(Day)*Sex), adjust = "tukey")
shannon_out <- as.data.frame(shannon_out[2])
write.csv (shannon_out, file = "shannon_tukey.csv")

#invsimpson
invsimpson_mod <- lmer(invsimpson ~ Group*as.factor(meta_ob1$Day)*Sex + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=FALSE)
invsimpson_out <- emmeans(invsimpson_mod, list(pairwise ~ Group*as.factor(Day)*Sex), adjust = "tukey")
invsimpson_out <- as.data.frame(invsimpson_out[2])





library(emmeans)
library(car)

#chao

chao_mod <- lmer(as.numeric(meta_ob1$chao) ~ Group*as.factor(meta_ob1$Day)*Sex + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=FALSE)
Anova(chao_mod, type="III")
chao_mod <- lmer(as.numeric(meta_ob1$chao) ~ Group*as.factor(meta_ob1$Day) + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=FALSE)
chao_out <- emmeans(chao_mod, list(pairwise ~ Group*as.factor(Day)), adjust = "tukey")
chao_out <- as.data.frame(chao_out[2])
write.csv (chao_out, file = "chao_merged_sex_tukey.csv")

#shannon
shannon_mod <- lmer(shannon ~ Group*as.factor(meta_ob1$Day)*Sex + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=FALSE)
Anova(shannon_mod, type="III")
shannon_mod <- lmer(shannon ~ Group*as.factor(meta_ob1$Day) + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=FALSE)
shannon_out <- emmeans(shannon_mod, list(pairwise ~ Group*as.factor(Day)), adjust = "tukey")
shannon_out <- as.data.frame(shannon_out[2])
write.csv (shannon_out, file = "shannon_merged_sex_tukey.csv")

#invsimpson
invsimpson_mod <- lmer(invsimpson ~ Group*as.factor(meta_ob1$Day)*Sex + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=FALSE)
Anova(invsimpson_mod, type="III")
invsimpson_mod <- lmer(invsimpson ~ Group*as.factor(meta_ob1$Day) + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=FALSE)
invsimpson_out <- emmeans(invsimpson_mod, list(pairwise ~ Group*as.factor(Day)), adjust = "tukey")
invsimpson_out <- as.data.frame(invsimpson_out[2])
write.csv (invsimpson_out, file = "invsimpson_merged_sex_tukey.csv")









#beta diversity

#devtools::install_github("vmikk/metagMisc")
library(phyloseq)
library(metagMisc)
library(vegan)
library(ggplot2)



setwd("D://ABX_LOCO_2024/")

meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)

tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa_138_ASV.txt"))

#tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/maaslin_selected_taxa.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Group!="GA",]


nrow(meta_ob1)


#Get sample names
GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("CEF_16S_shared.csv",)
#fullMothur<-read.csv("maaslin_selected_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
library("ape")
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))

#random_tree = rtree(ntaxa(phyloseqobj.f ), rooted=TRUE, tip.label=taxa_names(phyloseqobj.f ))

#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1), random_tree)
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))



phyloseqobj.f_percent <- phyloseq_standardize_otu_abundance(phyloseqobj.f, method = "total")

merged_phy_percent1 <- phyloseqobj.f_percent 

#merged_phy_percent1 <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:50]), merged_phy_percent)




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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)





#write.csv(input_data, "Genus_percent.csv")







Female_ps <- subset_samples(merged_phy_percent1, Sex=="Female")
Female_ps_merged = merge_samples(Female_ps, "DayPlot2")
Female_ps_merged_percent <- phyloseq_standardize_otu_abundance(Female_ps_merged, method = "total")
Female_ps_merged_percent_2 <- prune_taxa(names(sort(taxa_sums(Female_ps_merged_percent),TRUE)[1:3]), Female_ps_merged_percent)




Male_ps <- subset_samples(merged_phy_percent1, Sex=="Male")
Male_ps_merged = merge_samples(Male_ps, "DayPlot2") 
Male_ps_merged_percent <- phyloseq_standardize_otu_abundance(Male_ps_merged, method = "total")
Male_ps_merged_percent_2 <- prune_taxa(names(sort(taxa_sums(Male_ps_merged_percent),TRUE)[1:3]), Male_ps_merged_percent)



library(ggplot2)

F_FB_plot <- plot_bar(Female_ps_merged_percent_2, fill="Rank7") +
 theme(axis.text = element_text(size=2, color="black"),
        #strip.text.x = element_text(size=5, color="black"),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        axis.text.x= element_text(color="black", size=1),
        axis.text.y= element_text(color="black", size=5),
        text = element_text(size=5, face=NULL),
        legend.text=element_text(size=6),
        axis.title.x = element_blank(),
        legend.position="bottom", legend.box = "vertical")


M_FB_plot <- plot_bar(Male_ps_merged_percent_2, fill="Rank2") +
  theme(axis.text = element_text(size=2, color="black"),
        #strip.text.x = element_text(size=5, color="black"),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        axis.text.x= element_text(color="black", size=1),
        axis.text.y= element_text(color="black", size=5),
        text = element_text(size=5, face=NULL),
        legend.text=element_text(size=6),
        axis.title.x = element_blank(),
        legend.position="bottom", legend.box = "vertical")

  
#F_FB_plot_data <- ggplot_build(F_FB_plot)$plot$data
#write.csv(F_FB_plot_data, file="F_FB_plot_data.csv")


#M_FB_plot_data <- ggplot_build(M_FB_plot)$plot$data
#write.csv(M_FB_plot_data, file="M_FB_plot_data.csv")





FB_df <- read.table("D://ABX_LOCO_2024/FB_plot_data.txt", header=T, row.names=1)


ggplot(FB_df, aes(fill=Rank2, y=Abundance, x=Day)) + 
  geom_bar(position="stack", stat="identity") +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  facet_grid(Sex~Group, space = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.45)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9))
#legend.justification=c(1,1), legend.position=c(1,1),legend.title=element_blank())


  
  
  



plot_data_df <- read.table("D://ABX_LOCO_2024/plot_data.txt", header=T, row.names=1)


library(dplyr)
#meta_ob1 %>% count(Sex, Group, Day, Cohort)

means <- plot_data_df %>%
  group_by(Sex, Group, DayPlot, Rank6) %>%
  summarise_at(vars(Abundance), list(name = mean)) %>%
  print(n=100)
write.csv(means, file="Rank6_Abundance_means.csv")


df <- read.table("D://ABX_LOCO_2024/Rank6_Abundance_means.txt", header=T, row.names=1)






merged_phy_percent2 <- prune_taxa(names(sort(taxa_sums(merged_phy_percent1),TRUE)[1:4]), merged_phy_percent1)
#sampleOrder = unique(sample_names(SI_pruned))
#taxaOrder = rev(unique(taxa_names(SI_pruned)))




Female_ps <- subset_samples(merged_phy_percent2, Sex=="Female")


Male_ps <- subset_samples(merged_phy_percent2, Sex=="Male")






Female_ps_C_plot <- plot_bar(Female_ps, fill="Rank2") +
  facet_wrap(~DayPlot2, scales = "free_x", nrow=1) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","gray27","magenta","pink","cyan","gray","orange","brown","purple", "#006003", "#FF0600")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","gray27","magenta","pink","cyan","gray","orange","brown","purple", "#006003",  "#FF0600")) +
  theme_minimal() +
  theme(legend.position="bottom", legend.box = "bottom",
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


Male_ps_C_plot <- plot_bar(Male_ps, fill="Rank2") +
  facet_wrap(~DayPlot2, scales = "free_x", nrow=1) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","gray27","magenta","pink","cyan","gray","orange","brown","purple", "#006003", "#FF0600")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","gray27","magenta","pink","cyan","gray","orange","brown","purple", "#006003",  "#FF0600")) +
  theme_minimal() +
  theme(legend.position="bottom", legend.box = "bottom",
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


library(ggpubr)

Control_Cohort_plot2 <- ggarrange(Female_ps_C_plot,
                                  Male_ps_C_plot,
                                  common.legend = TRUE, legend = "bottom",
                                  ncol = 1, nrow = 2, align="hv")
















###Day5 plot

merged_phy_percent2 <- subset_samples(merged_phy_percent1, DayPlot=="C5")
merged_phy_percent3 <- phyloseq_standardize_otu_abundance(merged_phy_percent2, method = "total")
Female_ps_merged_percent_4 <- prune_taxa(names(sort(taxa_sums(merged_phy_percent3),TRUE)[1:10]), merged_phy_percent3)



F_Day5_plot <- plot_bar(Female_ps_merged_percent_4, fill="Rank6") +
  theme(axis.text = element_text(size=2, color="black"),
        #strip.text.x = element_text(size=5, color="black"),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        axis.text.x= element_text(color="black", size=1),
        axis.text.y= element_text(color="black", size=5),
        text = element_text(size=5, face=NULL),
        legend.text=element_text(size=6),
        axis.title.x = element_blank(),
        legend.position="bottom", legend.box = "vertical") +
  facet_wrap(~Sex_Group, scales = "free_x", ncol = 3, nrow = 2)
  



library(pairwiseAdonis)

merged_phy_percent2 <- subset_samples(merged_phy_percent1, DayPlot=="C5" & Group=="GA")
x <- phyloseq::distance(merged_phy_percent2, method ="bray")
Day5_meta <- meta_ob1[meta_ob1$DayPlot=="C5" & meta_ob1$Group=="GA",]
factors <- Day5_meta$Sex_Group
pairwise.adonis(x, factors, perm = 1000000, p.adjust(method="none"))
















###FB_ratio

FB_df <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)

FB_df$FB <- FB_df$Firmicutes/FB_df$Bacteroidetes


ggplot(data = FB_df,
       aes(x=as.factor(Day), y=log10(FB+0.0000001), fill = factor(Group))) +
  geom_boxplot(fatten=NULL, outlier.shape = NA, position = position_dodge(width=0.85)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.333, dodge.width = 0.85), aes(fill = Group), pch = 21) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), position = position_dodge(width = 0.83), width = 0.74, size=1) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  facet_wrap(~Sex, scales = "free_x", ncol = 2, nrow = 1) +
  ylim(0,100) +
  ylab("Shannon") +
  theme_classic() +
  scale_y_continuous(limits = c(0,5.5), labels = scales::number_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = 9)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9),
        strip.background = element_blank(),
        legend.title= element_blank(),
        strip.text.x = element_blank())





















#dist.bray <- vegdist(shared_phy_out1, method="bray", binary=FALSE)
#adonis2(formula = dist.bray  ~ Mouse_ID*Group*Day*Sex, data = meta_phy_out1, permutations=1000)







#Female_Controls_PCoA


Female_ps <- subset_samples(merged_phy_percent1, Sex=="Female")
Female_C_ps <- subset_samples(Female_ps,Group=="C")


#Female_Day1
pssd2veg <- function(Female_C_ps) {
  sd <- sample_data(Female_C_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Female_C_ps) {
  OTU <- otu_table(Female_C_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Female_C_meta <- as.data.frame(pssd2veg(Female_C_ps))
Female_C_shared <- as.data.frame(psotu2veg(Female_C_ps))


dist_bray <- capscale(Female_C_shared~1, distance="bray", binary=FALSE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Female_C_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Female_C_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Female_C_meta$CD, fill=Female_C_meta$CD)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  #scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Female_C_meta$CD), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))




library(pairwiseAdonis)

packageVersion("pairwiseAdonis")

x <- phyloseq::distance(Female_C_ps, method ="bray") 
factors <- Female_C_meta$CD
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="none"))



dist_bray1 <- vegdist(Female_C_shared, method="bray", binary=FALSE)
adonis2(formula = dist_bray1 ~ CD, data = Female_C_meta, permutations=1000)





#Male_Controls_PCoA


Male_ps <- subset_samples(merged_phy_percent1, Sex=="Male")
Male_C_ps <- subset_samples(Male_ps,Group=="C")


#Male_Day1
pssd2veg <- function(Male_C_ps) {
  sd <- sample_data(Male_C_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Male_C_ps) {
  OTU <- otu_table(Male_C_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Male_C_meta <- as.data.frame(pssd2veg(Male_C_ps))
Male_C_shared <- as.data.frame(psotu2veg(Male_C_ps))


dist_bray <- capscale(Male_C_shared~1, distance="bray", binary=FALSE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Male_C_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Male_C_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Male_C_meta$CD, fill=Male_C_meta$CD)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  #scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Male_C_meta$CD), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))




library(pairwiseAdonis)



x <- phyloseq::distance(Male_C_ps, method ="bray") 
factors <- Male_C_meta$CD
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="none"))






dist_bray1 <- vegdist(Male_C_shared, method="bray", binary=FALSE)
adonis2(formula = dist_bray1 ~ CD, data = Male_C_meta, permutations=1000)







#Control Cohort barplot
Female_ps <- subset_samples(merged_phy_percent1, Sex=="Female")
Female_C_ps_200 <- subset_samples(Female_ps,Group=="C")

Female_ps_C_plot <- plot_bar(Female_C_ps_200, fill="Rank6") +
  facet_wrap(~CD, scales = "free_x", nrow=1) +
    theme(axis.text = element_text(size=2, color="black"),
        #strip.text.x = element_text(size=5, color="black"),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        axis.text.x= element_text(color="black", size=1),
        axis.text.y= element_text(color="black", size=5),
        text = element_text(size=5, face=NULL),
        legend.text=element_text(size=6),
        axis.title.x = element_blank(),
        legend.position="bottom", legend.box = "vertical")




#Control Cohort barplot
Male_ps <- subset_samples(merged_phy_percent1, Sex=="Male")
Male_C_ps_200 <- subset_samples(Male_ps,Group=="C")

Male_ps_C_plot <- plot_bar(Male_C_ps_200, fill="Rank6") +
  facet_wrap(~CD, scales = "free_x", nrow=1) +
  theme(axis.text = element_text(size=2, color="black"),
        #strip.text.x = element_text(size=5, color="black"),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        axis.text.x= element_text(color="black", size=1),
        axis.text.y= element_text(color="black", size=5),
        text = element_text(size=5, face=NULL),
        legend.text=element_text(size=6),
        axis.title.x = element_blank(),
        legend.position="bottom", legend.box = "vertical")


Control_Cohort_plot <- ggarrange(Female_ps_C_plot ,
                                 Female_ps_C_plot ,
                                 common.legend = TRUE, legend = "bottom",
                                 ncol = 1, nrow = 2, align="hv")



library(dplyr)




F_ps <- merge_samples(Female_ps, Female_C_meta$CD, fun=mean)

M_ps <- merge_samples(Male_ps, Male_C_meta$CD, fun=mean)


Female_ps_C_plot <- plot_bar(Female_ps, fill="Rank6") +
  facet_wrap(~CD, scales = "free_x", nrow=1) +
  theme(axis.text = element_text(size=2, color="black"),
        #strip.text.x = element_text(size=5, color="black"),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        axis.text.x= element_text(color="black", size=5),
        axis.text.y= element_text(color="black", size=5),
        text = element_text(size=5, face=NULL),
        legend.text=element_text(size=6),
        axis.title.x = element_blank(),
        legend.position="bottom", legend.box = "right")

Male_ps_C_plot <- plot_bar(Male_ps, fill="Rank6") +
  facet_wrap(~CD, scales = "free_x", nrow=1) +
  theme(axis.text = element_text(size=2, color="black"),
        #strip.text.x = element_text(size=5, color="black"),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        axis.text.x= element_text(color="black", size=5),
        axis.text.y= element_text(color="black", size=5),
        text = element_text(size=5, face=NULL),
        legend.text=element_text(size=6),
        axis.title.x = element_blank(),
        legend.position="bottom", legend.box = "right")


library(ggpubr)

Control_Cohort_plot2 <- ggarrange(Female_ps_C_plot,
                                  Male_ps_C_plot,
                                 common.legend = TRUE, legend = "right",
                                 ncol = 1, nrow = 2, align="hv")









nrow(meta_phy_out1)



###Split male and female
Female_ps <- subset_samples(merged_phy_percent1, Sex=="Female")
Male_ps <- subset_samples(merged_phy_percent1, Sex=="Male")

###Split days
Female_Day1_ps <- subset_samples(Female_ps, Day=="1")
Female_Day2_ps <- subset_samples(Female_ps, Day=="2")
Female_Day5_ps <- subset_samples(Female_ps, Day=="5")
Female_Day8_ps <- subset_samples(Female_ps, Day=="8")
Female_Day12_ps <- subset_samples(Female_ps, Day=="12")
Female_Day65_ps <- subset_samples(Female_ps, Day=="65")
Female_Day100_ps <- subset_samples(Female_ps, Day=="100")

Male_Day1_ps <- subset_samples(Male_ps, Day=="1")
Male_Day2_ps <- subset_samples(Male_ps, Day=="2")
Male_Day5_ps <- subset_samples(Male_ps, Day=="5")
Male_Day8_ps <- subset_samples(Male_ps, Day=="8")
Male_Day12_ps <- subset_samples(Male_ps, Day=="12")
Male_Day65_ps <- subset_samples(Male_ps, Day=="65")
Male_Day100_ps <- subset_samples(Male_ps, Day=="100")



###Get top 16 taxa from each day
Female_Day1_ps_16 <- prune_taxa(names(sort(taxa_sums(Female_Day1_ps),TRUE)[1:16]), Female_Day1_ps)
Female_Day2_ps_16 <- prune_taxa(names(sort(taxa_sums(Female_Day2_ps),TRUE)[1:16]), Female_Day2_ps)
Female_Day5_ps_16 <- prune_taxa(names(sort(taxa_sums(Female_Day5_ps),TRUE)[1:16]), Female_Day5_ps)
Female_Day8_ps_16 <- prune_taxa(names(sort(taxa_sums(Female_Day8_ps),TRUE)[1:16]), Female_Day8_ps)
Female_Day12_ps_16 <- prune_taxa(names(sort(taxa_sums(Female_Day12_ps),TRUE)[1:16]), Female_Day12_ps)
Female_Day65_ps_16 <- prune_taxa(names(sort(taxa_sums(Female_Day65_ps),TRUE)[1:16]), Female_Day65_ps)
Female_Day100_ps_16 <- prune_taxa(names(sort(taxa_sums(Female_Day100_ps),TRUE)[1:16]), Female_Day100_ps)



Male_Day1_ps_16 <- prune_taxa(names(sort(taxa_sums(Male_Day1_ps),TRUE)[1:16]), Male_Day1_ps)
Male_Day2_ps_16 <- prune_taxa(names(sort(taxa_sums(Male_Day2_ps),TRUE)[1:16]), Male_Day2_ps)
Male_Day5_ps_16 <- prune_taxa(names(sort(taxa_sums(Male_Day5_ps),TRUE)[1:16]), Male_Day5_ps)
Male_Day8_ps_16 <- prune_taxa(names(sort(taxa_sums(Male_Day8_ps),TRUE)[1:16]), Male_Day8_ps)
Male_Day12_ps_16 <- prune_taxa(names(sort(taxa_sums(Male_Day12_ps),TRUE)[1:16]), Male_Day12_ps)
Male_Day65_ps_16 <- prune_taxa(names(sort(taxa_sums(Male_Day65_ps),TRUE)[1:16]), Male_Day65_ps)
Male_Day100_ps_16 <- prune_taxa(names(sort(taxa_sums(Male_Day100_ps),TRUE)[1:16]), Male_Day100_ps)



#Split groups
Female_Day1_C_ps_16 <- subset_samples(Female_Day1_ps_16, Group=="C")
Female_Day2_C_ps_16 <- subset_samples(Female_Day2_ps_16, Group=="C")
Female_Day5_C_ps_16 <- subset_samples(Female_Day5_ps_16, Group=="C")
Female_Day8_C_ps_16 <- subset_samples(Female_Day8_ps_16, Group=="C")
Female_Day12_C_ps_16 <- subset_samples(Female_Day12_ps_16, Group=="C")
Female_Day65_C_ps_16 <- subset_samples(Female_Day65_ps_16, Group=="C")
Female_Day100_C_ps_16 <- subset_samples(Female_Day100_ps_16, Group=="C")

Female_Day1_GA_ps_16 <- subset_samples(Female_Day1_ps_16, Group=="GA")
Female_Day2_GA_ps_16 <- subset_samples(Female_Day2_ps_16, Group=="GA")
Female_Day5_GA_ps_16 <- subset_samples(Female_Day5_ps_16, Group=="GA")
Female_Day8_GA_ps_16 <- subset_samples(Female_Day8_ps_16, Group=="GA")
Female_Day12_GA_ps_16 <- subset_samples(Female_Day12_ps_16, Group=="GA")
Female_Day65_GA_ps_16 <- subset_samples(Female_Day65_ps_16, Group=="GA")
Female_Day100_GA_ps_16 <- subset_samples(Female_Day100_ps_16, Group=="GA")

Female_Day1_IA_ps_16 <- subset_samples(Female_Day1_ps_16, Group=="IA")
Female_Day2_IA_ps_16 <- subset_samples(Female_Day2_ps_16, Group=="IA")
Female_Day5_IA_ps_16 <- subset_samples(Female_Day5_ps_16, Group=="IA")
Female_Day8_IA_ps_16 <- subset_samples(Female_Day8_ps_16, Group=="IA")
Female_Day12_IA_ps_16 <- subset_samples(Female_Day12_ps_16, Group=="IA")
Female_Day65_IA_ps_16 <- subset_samples(Female_Day65_ps_16, Group=="IA")
Female_Day100_IA_ps_16 <- subset_samples(Female_Day100_ps_16, Group=="IA")




Male_Day1_C_ps_16 <- subset_samples(Male_Day1_ps_16, Group=="C")
Male_Day2_C_ps_16 <- subset_samples(Male_Day2_ps_16, Group=="C")
Male_Day5_C_ps_16 <- subset_samples(Male_Day5_ps_16, Group=="C")
Male_Day8_C_ps_16 <- subset_samples(Male_Day8_ps_16, Group=="C")
Male_Day12_C_ps_16 <- subset_samples(Male_Day12_ps_16, Group=="C")
Male_Day65_C_ps_16 <- subset_samples(Male_Day65_ps_16, Group=="C")
Male_Day100_C_ps_16 <- subset_samples(Male_Day100_ps_16, Group=="C")

Male_Day1_GA_ps_16 <- subset_samples(Male_Day1_ps_16, Group=="GA")
Male_Day2_GA_ps_16 <- subset_samples(Male_Day2_ps_16, Group=="GA")
Male_Day5_GA_ps_16 <- subset_samples(Male_Day5_ps_16, Group=="GA")
Male_Day8_GA_ps_16 <- subset_samples(Male_Day8_ps_16, Group=="GA")
Male_Day12_GA_ps_16 <- subset_samples(Male_Day12_ps_16, Group=="GA")
Male_Day65_GA_ps_16 <- subset_samples(Male_Day65_ps_16, Group=="GA")
Male_Day100_GA_ps_16 <- subset_samples(Male_Day100_ps_16, Group=="GA")

Male_Day1_IA_ps_16 <- subset_samples(Male_Day1_ps_16, Group=="IA")
Male_Day2_IA_ps_16 <- subset_samples(Male_Day2_ps_16, Group=="IA")
Male_Day5_IA_ps_16 <- subset_samples(Male_Day5_ps_16, Group=="IA")
Male_Day8_IA_ps_16 <- subset_samples(Male_Day8_ps_16, Group=="IA")
Male_Day12_IA_ps_16 <- subset_samples(Male_Day12_ps_16, Group=="IA")
Male_Day65_IA_ps_16 <- subset_samples(Male_Day65_ps_16, Group=="IA")
Male_Day100_IA_ps_16 <- subset_samples(Male_Day100_ps_16, Group=="IA")






###Merge groups
Female_merged_ps <- merge_phyloseq(Female_Day1_ps_16, Female_Day2_ps_16, Female_Day5_ps_16, Female_Day8_ps_16, Female_Day12_ps_16, Female_Day65_ps_16, Female_Day100_ps_16) 

Male_merged_ps <- merge_phyloseq(Male_Day1_ps_16, Male_Day2_ps_16, Male_Day5_ps_16, Male_Day8_ps_16, Male_Day12_ps_16, Male_Day65_ps_16, Male_Day100_ps_16) 


F_gen_list <- taxa_names(Female_merged_ps)
str(F_gen_list)

M_gen_list <- taxa_names(Male_merged_ps)
str(M_gen_list)



###Get selected taxa from original phyloseq object
F_gen_of_int <- subset_taxa(merged_phy_percent1, Rank6 %in% F_gen_list[1:35])
M_gen_of_int <- subset_taxa(merged_phy_percent1, Rank6 %in% M_gen_list[1:34])


Female_selected_35_ps <- subset_samples(F_gen_of_int, Sex=="Female")
Male_selected_34_ps <- subset_samples(M_gen_of_int, Sex=="Male")










###Separate groups again containing all selected taxa for clustering



###Split days

Female_35_Day1_ps <- subset_samples(Female_selected_35_ps, Day=="1")
Female_35_Day2_ps <- subset_samples(Female_selected_35_ps, Day=="2")
Female_35_Day5_ps <- subset_samples(Female_selected_35_ps, Day=="5")
Female_35_Day8_ps <- subset_samples(Female_selected_35_ps, Day=="8")
Female_35_Day12_ps <- subset_samples(Female_selected_35_ps, Day=="12")
Female_35_Day65_ps <- subset_samples(Female_selected_35_ps, Day=="65")
Female_35_Day100_ps <- subset_samples(Female_selected_35_ps, Day=="100")

Male_34_Day1_ps <- subset_samples(Male_selected_34_ps, Day=="1")
Male_34_Day2_ps <- subset_samples(Male_selected_34_ps, Day=="2")
Male_34_Day5_ps <- subset_samples(Male_selected_34_ps, Day=="5")
Male_34_Day8_ps <- subset_samples(Male_selected_34_ps, Day=="8")
Male_34_Day12_ps <- subset_samples(Male_selected_34_ps, Day=="12")
Male_34_Day65_ps <- subset_samples(Male_selected_34_ps, Day=="65")
Male_34_Day100_ps <- subset_samples(Male_selected_34_ps, Day=="100")





#Split groups

Female_35_Day1_C_ps_16 <- subset_samples(Female_35_Day1_ps, Group=="C")
Female_35_Day2_C_ps_16 <- subset_samples(Female_35_Day2_ps, Group=="C")
Female_35_Day5_C_ps_16 <- subset_samples(Female_35_Day5_ps, Group=="C")
Female_35_Day8_C_ps_16 <- subset_samples(Female_35_Day8_ps, Group=="C")
Female_35_Day12_C_ps_16 <- subset_samples(Female_35_Day12_ps, Group=="C")
Female_35_Day65_C_ps_16 <- subset_samples(Female_35_Day65_ps, Group=="C")
Female_35_Day100_C_ps_16 <- subset_samples(Female_35_Day100_ps, Group=="C")

Female_35_Day1_GA_ps_16 <- subset_samples(Female_35_Day1_ps, Group=="GA")
Female_35_Day2_GA_ps_16 <- subset_samples(Female_35_Day2_ps, Group=="GA")
Female_35_Day5_GA_ps_16 <- subset_samples(Female_35_Day5_ps, Group=="GA")
Female_35_Day8_GA_ps_16 <- subset_samples(Female_35_Day8_ps, Group=="GA")
Female_35_Day12_GA_ps_16 <- subset_samples(Female_35_Day12_ps, Group=="GA")
Female_35_Day65_GA_ps_16 <- subset_samples(Female_35_Day65_ps, Group=="GA")
Female_35_Day100_GA_ps_16 <- subset_samples(Female_35_Day100_ps, Group=="GA")

Female_35_Day1_IA_ps_16 <- subset_samples(Female_35_Day1_ps, Group=="IA")
Female_35_Day2_IA_ps_16 <- subset_samples(Female_35_Day2_ps, Group=="IA")
Female_35_Day5_IA_ps_16 <- subset_samples(Female_35_Day5_ps, Group=="IA")
Female_35_Day8_IA_ps_16 <- subset_samples(Female_35_Day8_ps, Group=="IA")
Female_35_Day12_IA_ps_16 <- subset_samples(Female_35_Day12_ps, Group=="IA")
Female_35_Day65_IA_ps_16 <- subset_samples(Female_35_Day65_ps, Group=="IA")
Female_35_Day100_IA_ps_16 <- subset_samples(Female_35_Day100_ps, Group=="IA")



Male_34_Day1_C_ps_16 <- subset_samples(Male_34_Day1_ps, Group=="C")
Male_34_Day2_C_ps_16 <- subset_samples(Male_34_Day2_ps, Group=="C")
Male_34_Day5_C_ps_16 <- subset_samples(Male_34_Day5_ps, Group=="C")
Male_34_Day8_C_ps_16 <- subset_samples(Male_34_Day8_ps, Group=="C")
Male_34_Day12_C_ps_16 <- subset_samples(Male_34_Day12_ps, Group=="C")
Male_34_Day65_C_ps_16 <- subset_samples(Male_34_Day65_ps, Group=="C")
Male_34_Day100_C_ps_16 <- subset_samples(Male_34_Day100_ps, Group=="C")

Male_34_Day1_GA_ps_16 <- subset_samples(Male_34_Day1_ps, Group=="GA")
Male_34_Day2_GA_ps_16 <- subset_samples(Male_34_Day2_ps, Group=="GA")
Male_34_Day5_GA_ps_16 <- subset_samples(Male_34_Day5_ps, Group=="GA")
Male_34_Day8_GA_ps_16 <- subset_samples(Male_34_Day8_ps, Group=="GA")
Male_34_Day12_GA_ps_16 <- subset_samples(Male_34_Day12_ps, Group=="GA")
Male_34_Day65_GA_ps_16 <- subset_samples(Male_34_Day65_ps, Group=="GA")
Male_34_Day100_GA_ps_16 <- subset_samples(Male_34_Day100_ps, Group=="GA")

Male_34_Day1_IA_ps_16 <- subset_samples(Male_34_Day1_ps, Group=="IA")
Male_34_Day2_IA_ps_16 <- subset_samples(Male_34_Day2_ps, Group=="IA")
Male_34_Day5_IA_ps_16 <- subset_samples(Male_34_Day5_ps, Group=="IA")
Male_34_Day8_IA_ps_16 <- subset_samples(Male_34_Day8_ps, Group=="IA")
Male_34_Day12_IA_ps_16 <- subset_samples(Male_34_Day12_ps, Group=="IA")
Male_34_Day65_IA_ps_16 <- subset_samples(Male_34_Day65_ps, Group=="IA")
Male_34_Day100_IA_ps_16 <- subset_samples(Male_34_Day100_ps, Group=="IA")








###Cluster each group individually using all 34 and 31 taxa for each sample

library(RFLPtools)
library(dendextend)
distance = distance(Female_35_Day1_C_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day1_C_ps_16.txt", "hello", k =12, append = FALSE, dec = ",")

distance = distance(Female_35_Day2_C_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day2_C_ps_16.txt", "hello", k =6, append = FALSE, dec = ",")

distance = distance(Female_35_Day5_C_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day5_C_ps_16.txt", "hello", k =12, append = FALSE, dec = ",")

distance = distance(Female_35_Day8_C_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day8_C_ps_16.txt", "hello", k =6, append = FALSE, dec = ",")

distance = distance(Female_35_Day12_C_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day12_C_ps_16.txt", "hello", k =6, append = FALSE, dec = ",")

distance = distance(Female_35_Day65_C_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day12_C_ps_65.txt", "hello", k =6, append = FALSE, dec = ",")

distance = distance(Female_35_Day100_C_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day100_C_ps_65.txt", "hello", k =6, append = FALSE, dec = ",")






distance = distance(Female_35_Day1_GA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day1_GA_ps_16.txt", "hello", k =14, append = FALSE, dec = ",")

distance = distance(Female_35_Day2_GA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day2_GA_ps_16.txt", "hello", k =7, append = FALSE, dec = ",")

distance = distance(Female_35_Day5_GA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day5_GA_ps_16.txt", "hello", k =14, append = FALSE, dec = ",")

distance = distance(Female_35_Day8_GA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day8_GA_ps_16.txt", "hello", k =7, append = FALSE, dec = ",")

distance = distance(Female_35_Day12_GA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day12_GA_ps_16.txt", "hello", k =7, append = FALSE, dec = ",")

distance = distance(Female_35_Day65_GA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day12_GA_ps_65.txt", "hello", k =7, append = FALSE, dec = ",")

distance = distance(Female_35_Day100_GA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day100_GA_ps_65.txt", "hello", k =7, append = FALSE, dec = ",")








distance = distance(Female_35_Day1_IA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day1_IA_ps_16.txt", "hello", k =14, append = FALSE, dec = ",")

distance = distance(Female_35_Day2_IA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day2_IA_ps_16.txt", "hello", k =7, append = FALSE, dec = ",")

distance = distance(Female_35_Day5_IA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day5_IA_ps_16.txt", "hello", k =14, append = FALSE, dec = ",")

distance = distance(Female_35_Day8_IA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day8_IA_ps_16.txt", "hello", k =7, append = FALSE, dec = ",")

distance = distance(Female_35_Day12_IA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day12_IA_ps_16.txt", "hello", k =7, append = FALSE, dec = ",")

distance = distance(Female_35_Day65_IA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day12_IA_ps_65.txt", "hello", k =7, append = FALSE, dec = ",")

distance = distance(Female_35_Day100_IA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_35_Day100_IA_ps_65.txt", "hello", k =7, append = FALSE, dec = ",")












distance = distance(Male_34_Day1_C_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day1_C_ps_16.txt", "hello", k =12, append = FALSE, dec = ",")

distance = distance(Male_34_Day2_C_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day2_C_ps_16.txt", "hello", k =6, append = FALSE, dec = ",")

distance = distance(Male_34_Day5_C_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day5_C_ps_16.txt", "hello", k =12, append = FALSE, dec = ",")

distance = distance(Male_34_Day8_C_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day8_C_ps_16.txt", "hello", k =6, append = FALSE, dec = ",")

distance = distance(Male_34_Day12_C_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day12_C_ps_16.txt", "hello", k =6, append = FALSE, dec = ",")

distance = distance(Male_34_Day65_C_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day12_C_ps_65.txt", "hello", k =6, append = FALSE, dec = ",")

distance = distance(Male_34_Day100_C_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day100_C_ps_65.txt", "hello", k =6, append = FALSE, dec = ",")






distance = distance(Male_34_Day1_GA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day1_GA_ps_16.txt", "hello", k =14, append = FALSE, dec = ",")

distance = distance(Male_34_Day2_GA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day2_GA_ps_16.txt", "hello", k =7, append = FALSE, dec = ",")

distance = distance(Male_34_Day5_GA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day5_GA_ps_16.txt", "hello", k =14, append = FALSE, dec = ",")

distance = distance(Male_34_Day8_GA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day8_GA_ps_16.txt", "hello", k =7, append = FALSE, dec = ",")

distance = distance(Male_34_Day12_GA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day12_GA_ps_16.txt", "hello", k =7, append = FALSE, dec = ",")

distance = distance(Male_34_Day65_GA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day12_GA_ps_65.txt", "hello", k =7, append = FALSE, dec = ",")

distance = distance(Male_34_Day100_GA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day100_GA_ps_65.txt", "hello", k =7, append = FALSE, dec = ",")








distance = distance(Male_34_Day1_IA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day1_IA_ps_16.txt", "hello", k =14, append = FALSE, dec = ",")

distance = distance(Male_34_Day2_IA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day2_IA_ps_16.txt", "hello", k =7, append = FALSE, dec = ",")

distance = distance(Male_34_Day5_IA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day5_IA_ps_16.txt", "hello", k =14, append = FALSE, dec = ",")

distance = distance(Male_34_Day8_IA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day8_IA_ps_16.txt", "hello", k =7, append = FALSE, dec = ",")

distance = distance(Male_34_Day12_IA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day12_IA_ps_16.txt", "hello", k =7, append = FALSE, dec = ",")

distance = distance(Male_34_Day65_IA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day12_IA_ps_65.txt", "hello", k =7, append = FALSE, dec = ",")

distance = distance(Male_34_Day100_IA_ps_16, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Male_34_Day100_IA_ps_65.txt", "hello", k =7, append = FALSE, dec = ",")






# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(Female_selected_35_ps) {
  sd <- sample_data(Female_selected_35_ps)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(Female_selected_35_ps) {
  OTU <- otu_table(Female_selected_35_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


taxa_names(Female_selected_35_ps)
tax_table(Female_selected_35_ps)[,6]
taxa_names(Female_selected_35_ps) <- tax_table(Female_selected_35_ps)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(Female_selected_35_ps))
shared_phy_out1 <- as.data.frame(psotu2veg(Female_selected_35_ps))





write.table(shared_phy_out1, file='D://ABX_LOCO_2024/Female_selected_35_shared.txt', sep='\t', col.names = NA)
write.table(meta_phy_out1, file='D://ABX_LOCO_2024/Female_selected_35_meta.txt', sep='\t', col.names = NA)













# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(Male_selected_34_ps) {
  sd <- sample_data(Male_selected_34_ps)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(Male_selected_34_ps) {
  OTU <- otu_table(Male_selected_34_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


taxa_names(Male_selected_34_ps)
tax_table(Male_selected_34_ps)[,6]
taxa_names(Male_selected_34_ps) <- tax_table(Male_selected_34_ps)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(Male_selected_34_ps))
shared_phy_out1 <- as.data.frame(psotu2veg(Male_selected_34_ps))





write.table(shared_phy_out1, file='D://ABX_LOCO_2024/Male_selected_34_shared.txt', sep='\t', col.names = NA)
write.table(meta_phy_out1, file='D://ABX_LOCO_2024/Male_selected_34_meta.txt', sep='\t', col.names = NA)







###Generate heatmaps with selected taxa




#input_data  <- as.data.frame(t(shared_phy_out1))
#input_metadata <- as.data.frame(meta_phy_out1)


#write.table(input_data, file='D://ABX_LOCO_2024/input_data.txt', sep='\t', col.names = NA)
#write.table(input_metadata, file='D://ABX_LOCO_2024/input_metadata.txt', sep='\t', col.names = NA)


F_data <- read.table("D://ABX_LOCO_2024/input_data.txt", header=T, row.names="gene")
F_meta <- read.table("D://ABX_LOCO_2024/input_metadata.txt", header=T, row.names=1)



ncol(F_data)
#F_meta_selected <- F_meta[F_meta$Group=="IA",]
F_data_selected <-  F_data[,colnames(F_data)%in%rownames(F_meta)]
ncol(F_data_selected)

F_data_selected <- F_data_selected[,rownames(F_meta)]

F_df_mat = data.matrix(F_data_selected)



#sample_Group <- as.data.frame(F_meta$Group)
#sample_Day <- as.data.frame(F_meta$DayPlot)
annotation_F_data <- F_meta[,c("Group","DayPlot")]
#row.names(annotation_F_data) <- colnames(F_data)


library(vegan)
drows = vegdist(F_data, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
#dcols =vegdist(t(F_data), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)


library(pheatmap2)
# blue and green color palette
color <- colorRampPalette((c("blue", "white", "red")))(50)

library(ggplotify)
Female_heatmap <- pheatmap2(F_df_mat, cluster_cols = FALSE, clustering_method = "ward.D2", clustering_distance_rows = "bray", treeheight_row=0, color = color, annotation_col = annotation_F_data, cluster_rows = TRUE, scale = "column", fontsize_col=0.000001, fontsize_row = 9)



str(Female_heatmap)









F_data <- read.table("D://ABX_LOCO_2024/Female_selected_35_shared.txt", header=T, row.names="gene")
F_meta <- read.table("D://ABX_LOCO_2024/Female_selected_35_meta.txt", header=T, row.names=1)

M_data <- read.table("D://ABX_LOCO_2024/Male_selected_34_shared.txt", header=T, row.names="gene")
M_meta <- read.table("D://ABX_LOCO_2024/Male_selected_34_meta.txt", header=T, row.names=1)





ncol(F_data)
F_meta_selected <- F_meta[F_meta$Group=="GA",]
F_data_selected <-  F_data[setdiff(names(F_data), rownames(F_meta_selected))]
ncol(F_data_selected)
nrow(F_meta_selected)
F_meta <- F_meta[F_meta$Group!="GA",]
nrow(F_meta)

ncol(M_data)
M_meta_selected <- M_meta[M_meta$Group=="GA",]
nrow(M_meta_selected)
M_data_selected <-  M_data[setdiff(names(M_data), rownames(M_meta_selected ))]
ncol(M_data_selected)
M_meta <- M_meta[M_meta$Group!="GA",]
nrow(M_meta)


ncol(M_data_selected)
nrow(M_data_selected)



F_data <- F_data_selected

M_data <- M_data_selected




F_df_mat = data.matrix(F_data)


sample_Group <- as.data.frame(F_meta$Group)
sample_Day <- as.data.frame(F_meta$DayPlot)
annotation_F_data <- as.data.frame(cbind(sample_Group, sample_Day))
row.names(annotation_F_data) <- colnames(F_data)


library(vegan)
drows = vegdist(F_data, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
#dcols =vegdist(t(F_data), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)


library(pheatmap2)
# blue and green color palette
color <- colorRampPalette((c("blue", "white", "red")))(50)

library(ggplotify)
Female_heatmap <- pheatmap2(F_df_mat, cluster_cols = FALSE, clustering_method = "ward.D2", clustering_distance_rows = "bray", treeheight_row=0, color = color, annotation_col = annotation_F_data, cluster_rows = TRUE, scale = "column", fontsize_col=0.000001, fontsize_row = 9)






M_df_mat = data.matrix(M_data)


sample_Group <- as.data.frame(M_meta$Group)
sample_Day <- as.data.frame(M_meta$DayPlot)
annotation_M_data <- as.data.frame(cbind(sample_Group, sample_Day))
row.names(annotation_M_data) <- colnames(M_data)


library(vegan)
drows = vegdist(M_data, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
#dcols =vegdist(t(M_data), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)


library(pheatmap2)
# blue and green color palette
color <- colorRampPalette((c("blue", "white", "red")))(50)

library(ggplotify)
Male_heatmap <- pheatmap2(M_df_mat, cluster_cols = FALSE, clustering_method = "ward.D2", clustering_distance_rows = "bray", treeheight_row=0, color = color, annotation_col = annotation_M_data, cluster_rows = TRUE, scale = "column", fontsize_col=0.000001, fontsize_row = 9)














M_df_mat = data.matrix(M_data)


sample_Group <- as.data.frame(M_meta$Group)
sample_Day <- as.data.frame(M_meta$DayPlot)
annotation_M_data <- as.data.frame(cbind(sample_Group, sample_Day))
row.names(annotation_M_data) <- colnames(M_data)


library(vegan)
drows = vegdist(M_data, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
#dcols =vegdist(t(M_data), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)


library(pheatmap2)
# blue and green color palette
color <- colorRampPalette((c("blue", "white", "red")))(50)

Male_heatmap <- pheatmap2(M_df_mat, cluster_cols = FALSE, clustering_method = "ward.D2", clustering_distance_rows = "bray", treeheight_row=0, color = color, annotation_col = annotation_M_data, cluster_rows = TRUE, scale = "column", fontsize_col=0.000001, fontsize_row = 9)




















"#808000"
"#ff4500"
"#c71585"
"#00ff00"
"#00ffff"
"#0000ff"
"#1e90ff"



"#191970"
"#006400"
"#ff0000"
"#ffd700"
"#00ffff"
"#ff00ff"
"#ffb6c1"









library(ggplot2)
library(ggplotify)
library(pheatmap)
library(patchwork)















library(RFLPtools)
library(dendextend)
distance = distance(Female_Day1_C_ps_15, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_Day1_C_ps_15.txt", "hello", k =12, append = FALSE, dec = ",")

distance = distance(Female_Day2_C_ps_15, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_Day2_C_ps_15.txt", "hello", k =6, append = FALSE, dec = ",")

distance = distance(Female_Day5_C_ps_15, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_Day5_C_ps_15.txt", "hello", k =12, append = FALSE, dec = ",")

distance = distance(Female_Day8_C_ps_15, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_Day8_C_ps_15.txt", "hello", k =6, append = FALSE, dec = ",")

distance = distance(Female_Day12_C_ps_15, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_Day12_C_ps_15.txt", "hello", k =6, append = FALSE, dec = ",")

distance = distance(Female_Day65_C_ps_15, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_Day12_C_ps_65.txt", "hello", k =6, append = FALSE, dec = ",")

distance = distance(Female_Day100_C_ps_15, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_Day100_C_ps_65.txt", "hello", k =6, append = FALSE, dec = ",")






distance = distance(Female_Day1_GA_ps_15, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_Day1_GA_ps_15.txt", "hello", k =14, append = FALSE, dec = ",")

distance = distance(Female_Day2_GA_ps_15, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_Day2_GA_ps_15.txt", "hello", k =7, append = FALSE, dec = ",")

distance = distance(Female_Day5_GA_ps_15, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_Day5_GA_ps_15.txt", "hello", k =14, append = FALSE, dec = ",")

distance = distance(Female_Day8_GA_ps_15, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_Day8_GA_ps_15.txt", "hello", k =7, append = FALSE, dec = ",")

distance = distance(Female_Day12_GA_ps_15, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_Day12_GA_ps_15.txt", "hello", k =7, append = FALSE, dec = ",")

distance = distance(Female_Day65_GA_ps_15, method="bray", binary=FALSE)
hcluster = hclust(distance, method ="ward.D2")
write.hclust(hcluster, "Female_Day12_GA_ps_65.txt")












library(dendextend)
library(circlize)
library(magrittr)

distance = distance(merged_phy_percent1, method="bray", binary=FALSE)
#distance = dist(phyloseqobj.f, method ="bray")    

hcluster = hclust(distance, method ="ward.D2")
dend <- as.dendrogram(hcluster)

#plot(dend)


cols <- c("#3F94C1", "#084287", "#5A9D5A", "#F6D832", "#E41A1C", "#CF9C76", "#7fff00", "#084287", "#8DD3C7", "#E7298A", "#FFAD12", "#D0CD66", "#45A939", "#91569A", "#D58EC4", "#B6742A","#999999")

dend <- color_branches(dend, k = 3, col = cols)
dend %<>% set("labels_col", value = cols, k= 3)
dend %<>% set("labels_cex", 1)
dend %<>% set("branches_lwd", 1)

plot(dend)

#circlize_dendrogram(dend)



library(RFLPtools)
write.hclust(hcluster, "clusters3.txt", "hello", k = 3, append = FALSE, dec = ",")



















###PCoA


Day1_ps_no_GA <- subset_samples(merged_phy_percent1, Group!="GA")

Day1_ps <- subset_samples(Day1_ps_no_GA, Day=="1")
Day2_ps <- subset_samples(Day1_ps_no_GA, Day=="2")
Day5_ps <- subset_samples(Day1_ps_no_GA, Day=="5")
Day8_ps <- subset_samples(Day1_ps_no_GA, Day=="8")
Day12_ps <- subset_samples(Day1_ps_no_GA, Day=="12")
Day65_ps <- subset_samples(Day1_ps_no_GA, Day=="65")
Day100_ps <- subset_samples(Day1_ps_no_GA, Day=="100")






#Day1
pssd2veg <- function(Day1_ps) {
  sd <- sample_data(Day1_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Day1_ps) {
  OTU <- otu_table(Day1_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Day1_meta <- as.data.frame(pssd2veg(Day1_ps))
Day1_shared <- as.data.frame(psotu2veg(Day1_ps))


#Day2
pssd2veg <- function(Day2_ps) {
  sd <- sample_data(Day2_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Day2_ps) {
  OTU <- otu_table(Day2_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Day2_meta <- as.data.frame(pssd2veg(Day2_ps))
Day2_shared <- as.data.frame(psotu2veg(Day2_ps))

#Day5
pssd2veg <- function(Day5_ps) {
  sd <- sample_data(Day5_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Day5_ps) {
  OTU <- otu_table(Day5_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Day5_meta <- as.data.frame(pssd2veg(Day5_ps))
Day5_shared <- as.data.frame(psotu2veg(Day5_ps))

#Day8
pssd2veg <- function(Day8_ps) {
  sd <- sample_data(Day8_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Day8_ps) {
  OTU <- otu_table(Day8_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Day8_meta <- as.data.frame(pssd2veg(Day8_ps))
Day8_shared <- as.data.frame(psotu2veg(Day8_ps))

#Day12
pssd2veg <- function(Day12_ps) {
  sd <- sample_data(Day12_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Day12_ps) {
  OTU <- otu_table(Day12_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Day12_meta <- as.data.frame(pssd2veg(Day12_ps))
Day12_shared <- as.data.frame(psotu2veg(Day12_ps))


#Day65
pssd2veg <- function(Day65_ps) {
  sd <- sample_data(Day65_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Day65_ps) {
  OTU <- otu_table(Day65_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Day65_meta <- as.data.frame(pssd2veg(Day65_ps))
Day65_shared <- as.data.frame(psotu2veg(Day65_ps))


#Day100
pssd2veg <- function(Day100_ps) {
  sd <- sample_data(Day100_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Day100_ps) {
  OTU <- otu_table(Day100_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Day100_meta <- as.data.frame(pssd2veg(Day100_ps))
Day100_shared <- as.data.frame(psotu2veg(Day100_ps))



library(dplyr)


#Day1_PCoA
dist_bray <- capscale(Day1_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Day1_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Day1_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Day1_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 2, 21, 22)) +
  geom_point(mapping = aes(colour=Day1_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.0, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))


#Day2_PCoA
dist_bray <- capscale(Day2_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Day2_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Day2_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Day2_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 2, 21, 22)) +
  geom_point(mapping = aes(colour=Day2_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.0, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))


#Day5_PCoA
dist_bray <- capscale(Day5_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Day5_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Day5_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Day5_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 2, 21, 22)) +
  geom_point(mapping = aes(colour=Day5_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.0, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))


#Day8_PCoA
dist_bray <- capscale(Day8_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Day8_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Day8_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Day8_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 2, 21, 22)) +
  geom_point(mapping = aes(colour=Day8_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.0, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))


#Day12_PCoA
dist_bray <- capscale(Day12_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Day12_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Day12_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Day12_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 2, 21, 22)) +
  geom_point(mapping = aes(colour=Day12_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.0, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))


#Day65_PCoA
dist_bray <- capscale(Day65_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Day65_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Day65_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Day65_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 2, 21, 22)) +
  geom_point(mapping = aes(colour=Day65_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.0, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))


#Day100_PCoA
dist_bray <- capscale(Day100_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Day100_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Day100_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Day100_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 2, 21, 22)) +
  geom_point(mapping = aes(colour=Day100_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.0, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))






library(ggpubr)

bray_merged <- ggarrange(Day1_bray_plot,
                                Day2_bray_plot,
                                Day5_bray_plot,
                                Day8_bray_plot,
                                Day12_bray_plot,
                                Day65_bray_plot,
                                Day100_bray_plot,
                                common.legend = TRUE, legend = "none",
                                ncol = 1, nrow = 7, align="hv")





dist_bray1 <- vegdist(Day1_shared, method="bray", binary=FALSE)
adonis2(formula = dist_bray1 ~ Sex + Group, data = Day1_meta, permutations=10000)

dist_bray2 <- vegdist(Day2_shared, method="bray", binary=FALSE)
adonis2(formula = dist_bray2 ~ Sex + Group, data = Day2_meta, permutations=10000)

dist_bray5 <- vegdist(Day5_shared, method="bray", binary=FALSE)
adonis2(formula = dist_bray5 ~ Sex + Group, data = Day5_meta, permutations=10000)

dist_bray8 <- vegdist(Day8_shared, method="bray", binary=FALSE)
adonis2(formula = dist_bray8 ~ Sex + Group, data = Day8_meta, permutations=10000)

dist_bray12 <- vegdist(Day12_shared, method="bray", binary=FALSE)
adonis2(formula = dist_bray12 ~ Sex + Group, data = Day12_meta, permutations=10000)

dist_bray65 <- vegdist(Day65_shared, method="bray", binary=FALSE)
adonis2(formula = dist_bray65 ~ Sex + Group, data = Day65_meta, permutations=1000)

dist_bray100 <- vegdist(Day100_shared, method="bray", binary=FALSE)
adonis2(formula = dist_bray100 ~ Sex + Group, data = Day100_meta, permutations=1000000)



























#Female_Day1
pssd2veg <- function(Female_Day1_ps) {
  sd <- sample_data(Female_Day1_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Female_Day1_ps) {
  OTU <- otu_table(Female_Day1_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Female_Day1_meta <- as.data.frame(pssd2veg(Female_Day1_ps))
Female_Day1_shared <- as.data.frame(psotu2veg(Female_Day1_ps))

#Female_Day2
pssd2veg <- function(Female_Day2_ps) {
  sd <- sample_data(Female_Day2_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Female_Day2_ps) {
  OTU <- otu_table(Female_Day2_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Female_Day2_meta <- as.data.frame(pssd2veg(Female_Day2_ps))
Female_Day2_shared <- as.data.frame(psotu2veg(Female_Day2_ps))

#Female_Day5
pssd2veg <- function(Female_Day5_ps) {
  sd <- sample_data(Female_Day5_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Female_Day5_ps) {
  OTU <- otu_table(Female_Day5_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Female_Day5_meta <- as.data.frame(pssd2veg(Female_Day5_ps))
Female_Day5_shared <- as.data.frame(psotu2veg(Female_Day5_ps))

#Female_Day8
pssd2veg <- function(Female_Day8_ps) {
  sd <- sample_data(Female_Day8_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Female_Day8_ps) {
  OTU <- otu_table(Female_Day8_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Female_Day8_meta <- as.data.frame(pssd2veg(Female_Day8_ps))
Female_Day8_shared <- as.data.frame(psotu2veg(Female_Day8_ps))

#Female_Day12
pssd2veg <- function(Female_Day12_ps) {
  sd <- sample_data(Female_Day12_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Female_Day12_ps) {
  OTU <- otu_table(Female_Day12_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Female_Day12_meta <- as.data.frame(pssd2veg(Female_Day12_ps))
Female_Day12_shared <- as.data.frame(psotu2veg(Female_Day12_ps))

#Female_Day65
pssd2veg <- function(Female_Day65_ps) {
  sd <- sample_data(Female_Day65_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Female_Day65_ps) {
  OTU <- otu_table(Female_Day65_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Female_Day65_meta <- as.data.frame(pssd2veg(Female_Day65_ps))
Female_Day65_shared <- as.data.frame(psotu2veg(Female_Day65_ps))

#Female_Day100
pssd2veg <- function(Female_Day100_ps) {
  sd <- sample_data(Female_Day100_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Female_Day100_ps) {
  OTU <- otu_table(Female_Day100_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Female_Day100_meta <- as.data.frame(pssd2veg(Female_Day100_ps))
Female_Day100_shared <- as.data.frame(psotu2veg(Female_Day100_ps))





#Male_Day1
pssd2veg <- function(Male_Day1_ps) {
  sd <- sample_data(Male_Day1_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Male_Day1_ps) {
  OTU <- otu_table(Male_Day1_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Male_Day1_meta <- as.data.frame(pssd2veg(Male_Day1_ps))
Male_Day1_shared <- as.data.frame(psotu2veg(Male_Day1_ps))

#Male_Day2
pssd2veg <- function(Male_Day2_ps) {
  sd <- sample_data(Male_Day2_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Male_Day2_ps) {
  OTU <- otu_table(Male_Day2_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Male_Day2_meta <- as.data.frame(pssd2veg(Male_Day2_ps))
Male_Day2_shared <- as.data.frame(psotu2veg(Male_Day2_ps))

#Male_Day5
pssd2veg <- function(Male_Day5_ps) {
  sd <- sample_data(Male_Day5_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Male_Day5_ps) {
  OTU <- otu_table(Male_Day5_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Male_Day5_meta <- as.data.frame(pssd2veg(Male_Day5_ps))
Male_Day5_shared <- as.data.frame(psotu2veg(Male_Day5_ps))

#Male_Day8
pssd2veg <- function(Male_Day8_ps) {
  sd <- sample_data(Male_Day8_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Male_Day8_ps) {
  OTU <- otu_table(Male_Day8_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Male_Day8_meta <- as.data.frame(pssd2veg(Male_Day8_ps))
Male_Day8_shared <- as.data.frame(psotu2veg(Male_Day8_ps))

#Male_Day12
pssd2veg <- function(Male_Day12_ps) {
  sd <- sample_data(Male_Day12_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Male_Day12_ps) {
  OTU <- otu_table(Male_Day12_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Male_Day12_meta <- as.data.frame(pssd2veg(Male_Day12_ps))
Male_Day12_shared <- as.data.frame(psotu2veg(Male_Day12_ps))

#Male_Day65
pssd2veg <- function(Male_Day65_ps) {
  sd <- sample_data(Male_Day65_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Male_Day65_ps) {
  OTU <- otu_table(Male_Day65_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Male_Day65_meta <- as.data.frame(pssd2veg(Male_Day65_ps))
Male_Day65_shared <- as.data.frame(psotu2veg(Male_Day65_ps))

#Male_Day100
pssd2veg <- function(Male_Day100_ps) {
  sd <- sample_data(Male_Day100_ps)
  return(as(sd,"data.frame"))}
psotu2veg <- function(Male_Day100_ps) {
  OTU <- otu_table(Male_Day100_ps)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
Male_Day100_meta <- as.data.frame(pssd2veg(Male_Day100_ps))
Male_Day100_shared <- as.data.frame(psotu2veg(Male_Day100_ps))



library(dplyr)
library(ggrepel)



#Mantel test
veg_bray <- vegdist(Female_Day100_shared, distance="bray", binary=TRUE)
veg_bray <- as.matrix(vegdist(Female_Day100_shared, distance="bray", binary=FALSE))
mantel(veg_bray, veg_bray, method = "spearman", permutations =1000)











###Merge sex PCoA


###bray


#Female_Day1_PCoA
dist_bray <- capscale(Female_Day1_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Female_Day1_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Female_Day1_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Female_Day1_meta$Group, fill=Female_Day1_meta$Group, shape=Female_Day1_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Female_Day1_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))


#Female_Day2_PCoA
dist_bray <- capscale(Female_Day2_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Female_Day2_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Female_Day2_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Female_Day2_meta$Group, fill=Female_Day2_meta$Group, shape=Female_Day2_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Female_Day2_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

#Female_Day5_PCoA
dist_bray <- capscale(Female_Day5_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Female_Day5_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Female_Day5_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Female_Day5_meta$Group, fill=Female_Day5_meta$Group, shape=Female_Day5_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Female_Day5_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

#Female_Day8_PCoA
dist_bray <- capscale(Female_Day8_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Female_Day8_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Female_Day8_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Female_Day8_meta$Group, fill=Female_Day8_meta$Group, shape=Female_Day8_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Female_Day8_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

#Female_Day12_PCoA
dist_bray <- capscale(Female_Day12_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Female_Day12_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Female_Day12_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Female_Day12_meta$Group, fill=Female_Day12_meta$Group, shape=Female_Day12_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Female_Day12_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

#Female_Day65_PCoA
dist_bray <- capscale(Female_Day65_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Female_Day65_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Female_Day65_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Female_Day65_meta$Group, fill=Female_Day65_meta$Group, shape=Female_Day65_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Female_Day65_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

#Female_Day100_PCoA
dist_bray <- capscale(Female_Day100_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Female_Day100_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Female_Day100_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Female_Day100_meta$Group, fill=Female_Day100_meta$Group, shape=Female_Day100_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Female_Day100_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))








#Male_Day1_PCoA
dist_bray <- capscale(Male_Day1_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Male_Day1_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Male_Day1_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Male_Day1_meta$Group, fill=Male_Day1_meta$Group, shape=Male_Day1_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Male_Day1_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

#Male_Day2_PCoA
dist_bray <- capscale(Male_Day2_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Male_Day2_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Male_Day2_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Male_Day2_meta$Group, fill=Male_Day2_meta$Group, shape=Male_Day2_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Male_Day2_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

#Male_Day5_PCoA
dist_bray <- capscale(Male_Day5_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Male_Day5_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Male_Day5_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Male_Day5_meta$Group, fill=Male_Day5_meta$Group, shape=Male_Day5_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Male_Day5_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

#Male_Day8_PCoA
dist_bray <- capscale(Male_Day8_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Male_Day8_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Male_Day8_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Male_Day8_meta$Group, fill=Male_Day8_meta$Group, shape=Male_Day8_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Male_Day8_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))
4
#Male_Day12_PCoA
dist_bray <- capscale(Male_Day12_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Male_Day12_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Male_Day12_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Male_Day12_meta$Group, fill=Male_Day12_meta$Group, shape=Male_Day12_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Male_Day12_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

#Male_Day65_PCoA
dist_bray <- capscale(Male_Day65_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Male_Day65_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Male_Day65_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Male_Day65_meta$Group, fill=Male_Day65_meta$Group, shape=Male_Day65_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Male_Day65_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

#Male_Day100_PCoA
dist_bray <- capscale(Male_Day100_shared~1, distance="bray", binary=TRUE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Male_Day100_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Male_Day100_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Male_Day100_meta$Group, fill=Male_Day100_meta$Group, shape=Male_Day100_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Male_Day100_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))





library(ggpubr)

Female_bray_merged <- ggarrange(Female_Day1_bray_plot,
                                   Female_Day2_bray_plot,
                                   Female_Day5_bray_plot,
                                   Female_Day8_bray_plot,
                                   Female_Day12_bray_plot,
                                   Female_Day65_bray_plot,
                                   Female_Day100_bray_plot,
                                   common.legend = TRUE, legend = "none",
                                   ncol = 1, nrow = 7, align="hv")



Male_bray_merged <- ggarrange(Male_Day1_bray_plot,
                                 Male_Day2_bray_plot,
                                 Male_Day5_bray_plot,
                                 Male_Day8_bray_plot,
                                 Male_Day12_bray_plot,
                                 Male_Day65_bray_plot,
                                 Male_Day100_bray_plot,
                                 common.legend = TRUE, legend = "none",
                                 ncol = 1, nrow = 7, align="hv")



bray_merged <- ggarrange(Female_bray_merged ,
                            Male_bray_merged,
                            common.legend = TRUE, legend = "right",
                            ncol = 2, nrow = 1, align="hv")

















setwd("D://ABX_LOCO_2024/")



library(tidyverse)
library(vegan)

library(tidyverse)
library(vegan)

#days_wanted <- c(4:10, 141:150)

#eary_late_df <-
  
eary_late_df <- read_tsv("CEF_16S.shared") %>%
  select(Group, starts_with("ASV")) %>%
  pivot_longer(-Group) %>%
  separate(Group, into=c("animal", "day"), sep="X",
           remove=FALSE, convert=TRUE) %>%
  #filter(day %in% days_wanted) %>%
  group_by(Group) %>%
  mutate(N = sum(value)) %>%
  ungroup() %>%
  filter(N >= 8363) %>%
  select(-N, -animal, -day) %>%
  pivot_wider(names_from="name", values_from="value", values_fill=0) %>%
  column_to_rownames("Group")



bray <- avgdist(eary_late_df, dmethod="bray", sample=8363) %>%
  as.matrix() %>%
  as_tibble(rownames = "A") %>%
  pivot_longer(-A, names_to="B", values_to="distances")



#write.csv(bray, file="pairwise_bray_distances.csv")


meta_ob1 <- read.table("D://ABX_LOCO_2024/pairwise_bray_distances.txt", header=T, row.names=1)

#meta_ob1  <- meta_ob1[meta_ob1$Sex=="Male",]


library(ggplot2)


#8363 BC
ggplot(data = meta_ob1,
       aes(x=as.factor(Day), y=BC2, fill = factor(vs))) +
  geom_boxplot(fatten=NULL, outlier.shape = NA, position = position_dodge(width=0.85), alpha=1) +
  #geom_point(position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.85), aes(fill = vs), pch = 21, size=1) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), position = position_dodge(width = 0.83), width = 0.74, size=1) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  facet_wrap(~Sex, scales = "free_x", ncol = 1, nrow = 2) +
  #ylab("16S rRNA gene copies per mg of feces relative to same day control sample (%)") +
  theme_classic() +
  theme(axis.title.y = element_text(size=9)) +
  #scale_y_continuous(limits = c(0.975,1), labels = scales::number_format(accuracy = 0.001), breaks = scales::pretty_breaks(n = 4)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.45)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9),
       legend.title=element_blank())






ggplot(meta_ob1, aes(Group, Group2)) +
  geom_tile(aes(fill = BC2)) + 
  #geom_text(aes(label = round(BC, 1))) +
  scale_fill_gradient(low = "white", high = "red") +
  facet_wrap(~Day, scales = "free_x", ncol = 7)
  



























#p = Female_Day1_ps, m = "bray", s = "Sample", d = "Group"


require("phyloseq")
require("tidyverse")
require("dplyr")

wu.m = phyloseq::distance(Female_Day1_ps, "bray") %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(cols = c(everything(), -Var1), names_to = "Var2", values_to = "value") %>%
  filter(as.character(Var1) != as.character(Var2))
sd = data.frame(sample_data(Female_Day1_ps)) %>%
  select("Sample", "Group") %>%
  mutate_if(is.factor, as.character)
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")
colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

Female_Day1_bray_plot =  ggplot(data = wu.sd,
                      aes(x=Type2, y=value)) +
  geom_boxplot(aes(color = NULL)) +
  geom_jitter(width=0.2) +
scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  theme_classic() +
  theme(axis.title.y = element_text(size=9)) +
  #scale_y_continuous(limits = c(0.975,1), labels = scales::number_format(accuracy = 0.001), breaks = scales::pretty_breaks(n = 4)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.45)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9))













wu.m = phyloseq::distance(Female_Day2_ps, "bray") %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(cols = c(everything(), -Var1), names_to = "Var2", values_to = "value") %>%
  filter(as.character(Var1) != as.character(Var2))
sd = data.frame(sample_data(Female_Day2_ps)) %>%
  select("Sample", "Group") %>%
  mutate_if(is.factor, as.character)
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")
colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

Female_Day2_bray_plot =  ggplot(data = wu.sd,
                                aes(x=Type2, y=value)) +
  geom_boxplot(aes(color = NULL)) +
  geom_jitter(width=0.2) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  theme_classic() +
  theme(axis.title.y = element_text(size=9)) +
  #scale_y_continuous(limits = c(0.975,1), labels = scales::number_format(accuracy = 0.001), breaks = scales::pretty_breaks(n = 4)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.45)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9))


wu.m = phyloseq::distance(Female_Day5_ps, "bray") %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(cols = c(everything(), -Var1), names_to = "Var2", values_to = "value") %>%
  filter(as.character(Var1) != as.character(Var2))
sd = data.frame(sample_data(Female_Day5_ps)) %>%
  select("Sample", "Group") %>%
  mutate_if(is.factor, as.character)
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")
colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

Female_Day5_bray_plot =  ggplot(data = wu.sd,
                                aes(x=Type2, y=value)) +
  geom_boxplot(aes(color = NULL)) +
  geom_jitter(width=0.2) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  theme_classic() +
  theme(axis.title.y = element_text(size=9)) +
  #scale_y_continuous(limits = c(0.975,1), labels = scales::number_format(accuracy = 0.001), breaks = scales::pretty_breaks(n = 4)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.45)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9))


wu.m = phyloseq::distance(Female_Day8_ps, "bray") %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(cols = c(everything(), -Var1), names_to = "Var2", values_to = "value") %>%
  filter(as.character(Var1) != as.character(Var2))
sd = data.frame(sample_data(Female_Day8_ps)) %>%
  select("Sample", "Group") %>%
  mutate_if(is.factor, as.character)
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")
colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

Female_Day8_bray_plot =  ggplot(data = wu.sd,
                                aes(x=Type2, y=value)) +
  geom_boxplot(aes(color = NULL)) +
  geom_jitter(width=0.2) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  theme_classic() +
  theme(axis.title.y = element_text(size=9)) +
  #scale_y_continuous(limits = c(0.975,1), labels = scales::number_format(accuracy = 0.001), breaks = scales::pretty_breaks(n = 4)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.45)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9))

wu.m = phyloseq::distance(Female_Day12_ps, "bray") %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(cols = c(everything(), -Var1), names_to = "Var2", values_to = "value") %>%
  filter(as.character(Var1) != as.character(Var2))
sd = data.frame(sample_data(Female_Day12_ps)) %>%
  select("Sample", "Group") %>%
  mutate_if(is.factor, as.character)
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")
colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

Female_Day12_bray_plot =  ggplot(data = wu.sd,
                                 aes(x=Type2, y=value)) +
  geom_boxplot(aes(color = NULL)) +
  geom_jitter(width=0.2) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  theme_classic() +
  theme(axis.title.y = element_text(size=9)) +
  #scale_y_continuous(limits = c(0.975,1), labels = scales::number_format(accuracy = 0.001), breaks = scales::pretty_breaks(n = 4)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.45)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9))


wu.m = phyloseq::distance(Female_Day65_ps, "bray") %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(cols = c(everything(), -Var1), names_to = "Var2", values_to = "value") %>%
  filter(as.character(Var1) != as.character(Var2))
sd = data.frame(sample_data(Female_Day65_ps)) %>%
  select("Sample", "Group") %>%
  mutate_if(is.factor, as.character)
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")
colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

Female_Day65_bray_plot =  ggplot(data = wu.sd,
                                 aes(x=Type2, y=value)) +
  geom_boxplot(aes(color = NULL)) +
  geom_jitter(width=0.2) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  theme_classic() +
  theme(axis.title.y = element_text(size=9)) +
  #scale_y_continuous(limits = c(0.975,1), labels = scales::number_format(accuracy = 0.001), breaks = scales::pretty_breaks(n = 4)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.45)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9))


wu.m = phyloseq::distance(Female_Day100_ps, "bray") %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(cols = c(everything(), -Var1), names_to = "Var2", values_to = "value") %>%
  filter(as.character(Var1) != as.character(Var2))
sd = data.frame(sample_data(Female_Day100_ps)) %>%
  select("Sample", "Group") %>%
  mutate_if(is.factor, as.character)
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")
colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

Female_Day100_bray_plot =  ggplot(data = wu.sd,
                                  aes(x=Type2, y=value)) +
  geom_boxplot(aes(color = NULL)) +
  geom_jitter(width=0.2) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  theme_classic() +
  theme(axis.title.y = element_text(size=9)) +
  #scale_y_continuous(limits = c(0.975,1), labels = scales::number_format(accuracy = 0.001), breaks = scales::pretty_breaks(n = 4)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.45)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9))


library(ggpubr)
Female_bray_dist_merged <- ggarrange(Female_Day1_bray_plot,
                               Female_Day2_bray_plot,
                               Female_Day5_bray_plot,
                               Female_Day8_bray_plot,
                               Female_Day12_bray_plot,
                               Female_Day65_bray_plot,
                               Female_Day100_bray_plot,
                               common.legend = TRUE, legend = "none",
                               ncol = 7, nrow = 1, align="hv")



wu.m = phyloseq::distance(Male_Day1_ps, "bray") %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(cols = c(everything(), -Var1), names_to = "Var2", values_to = "value") %>%
  filter(as.character(Var1) != as.character(Var2))
sd = data.frame(sample_data(Male_Day1_ps)) %>%
  select("Sample", "Group") %>%
  mutate_if(is.factor, as.character)
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")
colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

Male_Day1_bray_plot =  ggplot(data = wu.sd,
                              aes(x=Type2, y=value)) +
  geom_boxplot(aes(color = NULL)) +
  geom_jitter(width=0.2) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  theme_classic() +
  theme(axis.title.y = element_text(size=9)) +
  #scale_y_continuous(limits = c(0.975,1), labels = scales::number_format(accuracy = 0.001), breaks = scales::pretty_breaks(n = 4)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.45)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9))


wu.m = phyloseq::distance(Male_Day2_ps, "bray") %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(cols = c(everything(), -Var1), names_to = "Var2", values_to = "value") %>%
  filter(as.character(Var1) != as.character(Var2))
sd = data.frame(sample_data(Male_Day2_ps)) %>%
  select("Sample", "Group") %>%
  mutate_if(is.factor, as.character)
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")
colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

Male_Day2_bray_plot =  ggplot(data = wu.sd,
                              aes(x=Type2, y=value)) +
  geom_boxplot(aes(color = NULL)) +
  geom_jitter(width=0.2) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  theme_classic() +
  theme(axis.title.y = element_text(size=9)) +
  #scale_y_continuous(limits = c(0.975,1), labels = scales::number_format(accuracy = 0.001), breaks = scales::pretty_breaks(n = 4)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.45)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9))


wu.m = phyloseq::distance(Male_Day5_ps, "bray") %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(cols = c(everything(), -Var1), names_to = "Var2", values_to = "value") %>%
  filter(as.character(Var1) != as.character(Var2))
sd = data.frame(sample_data(Male_Day5_ps)) %>%
  select("Sample", "Group") %>%
  mutate_if(is.factor, as.character)
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")
colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

Male_Day5_bray_plot =  ggplot(data = wu.sd,
                              aes(x=Type2, y=value)) +
  geom_boxplot(aes(color = NULL)) +
  geom_jitter(width=0.2) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  theme_classic() +
  theme(axis.title.y = element_text(size=9)) +
  #scale_y_continuous(limits = c(0.975,1), labels = scales::number_format(accuracy = 0.001), breaks = scales::pretty_breaks(n = 4)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.45)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9))


wu.m = phyloseq::distance(Male_Day8_ps, "bray") %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(cols = c(everything(), -Var1), names_to = "Var2", values_to = "value") %>%
  filter(as.character(Var1) != as.character(Var2))
sd = data.frame(sample_data(Male_Day8_ps)) %>%
  select("Sample", "Group") %>%
  mutate_if(is.factor, as.character)
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")
colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

Male_Day8_bray_plot =  ggplot(data = wu.sd,
                              aes(x=Type2, y=value)) +
  geom_boxplot(aes(color = NULL)) +
  geom_jitter(width=0.2) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  theme_classic() +
  theme(axis.title.y = element_text(size=9)) +
  #scale_y_continuous(limits = c(0.975,1), labels = scales::number_format(accuracy = 0.001), breaks = scales::pretty_breaks(n = 4)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.45)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9))

wu.m = phyloseq::distance(Male_Day12_ps, "bray") %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(cols = c(everything(), -Var1), names_to = "Var2", values_to = "value") %>%
  filter(as.character(Var1) != as.character(Var2))
sd = data.frame(sample_data(Male_Day12_ps)) %>%
  select("Sample", "Group") %>%
  mutate_if(is.factor, as.character)
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")
colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

Male_Day12_bray_plot =  ggplot(data = wu.sd,
                               aes(x=Type2, y=value)) +
  geom_boxplot(aes(color = NULL)) +
  geom_jitter(width=0.2) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  theme_classic() +
  theme(axis.title.y = element_text(size=9)) +
  #scale_y_continuous(limits = c(0.975,1), labels = scales::number_format(accuracy = 0.001), breaks = scales::pretty_breaks(n = 4)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.45)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9))


wu.m = phyloseq::distance(Male_Day65_ps, "bray") %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(cols = c(everything(), -Var1), names_to = "Var2", values_to = "value") %>%
  filter(as.character(Var1) != as.character(Var2))
sd = data.frame(sample_data(Male_Day65_ps)) %>%
  select("Sample", "Group") %>%
  mutate_if(is.factor, as.character)
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")
colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

Male_Day65_bray_plot =  ggplot(data = wu.sd,
                               aes(x=Type2, y=value)) +
  geom_boxplot(aes(color = NULL)) +
  geom_jitter(width=0.2) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  theme_classic() +
  theme(axis.title.y = element_text(size=9)) +
  #scale_y_continuous(limits = c(0.975,1), labels = scales::number_format(accuracy = 0.001), breaks = scales::pretty_breaks(n = 4)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.45)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9))


wu.m = phyloseq::distance(Male_Day100_ps, "bray") %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(cols = c(everything(), -Var1), names_to = "Var2", values_to = "value") %>%
  filter(as.character(Var1) != as.character(Var2))
sd = data.frame(sample_data(Male_Day100_ps)) %>%
  select("Sample", "Group") %>%
  mutate_if(is.factor, as.character)
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")
colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

Male_Day100_bray_plot =  ggplot(data = wu.sd,
                                aes(x=Type2, y=value)) +
  geom_boxplot(aes(color = NULL)) +
  geom_jitter(width=0.2) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  theme_classic() +
  theme(axis.title.y = element_text(size=9)) +
  #scale_y_continuous(limits = c(0.975,1), labels = scales::number_format(accuracy = 0.001), breaks = scales::pretty_breaks(n = 4)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.45)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9))


library(ggpubr)
Male_bray_dist_merged <- ggarrange(Male_Day1_bray_plot,
                                   Male_Day2_bray_plot,
                                   Male_Day5_bray_plot,
                                   Male_Day8_bray_plot,
                                   Male_Day12_bray_plot,
                                   Male_Day65_bray_plot,
                                   Male_Day100_bray_plot,
                                   common.legend = TRUE, legend = "none",
                                   ncol = 7, nrow = 1, align="hv")




merged_bray_dist_plot <- ggarrange(Male_bray_dist_merged,
                                   Male_bray_dist_merged,
                                   common.legend = TRUE, legend = "none",
                                   ncol = 1, nrow = 2, align="hv")



dev.off()



























library(ggplot2)
library(tidyverse)
dat <- matrix(rnorm(100, 3, 1), ncol = 10)
## the matrix needs names
names(dat) <- paste("X", 1:10)

## convert to tibble, add row identifier, and shape "long"
dat2 <-
  mat %>%
  as_tibble() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
  mutate(
    Var1 = factor(Var1, levels = 1:300),
    Var2 = factor(gsub("V", "", Var2), levels = 1:300)
  )
#> Warning: The `x` argument of `as_tibble.matrix()` must have unique column names if
#> `.name_repair` is omitted as of tibble 2.0.0.
#> ℹ Using compatibility `.name_repair`.

ggplot(dat2, aes(Var1, Var2)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = round(value, 1))) +
  scale_fill_gradient(low = "white", high = "red")






















###bray


#Female_Day1_PCoA
dist_bray <- capscale(Female_Day1_shared~1, distance="bray", binary=FALSE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Female_Day1_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Female_Day1_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Female_Day1_meta$Group, fill=Female_Day1_meta$Group, shape=Female_Day1_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Female_Day1_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))


#Female_Day2_PCoA
dist_bray <- capscale(Female_Day2_shared~1, distance="bray", binary=FALSE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Female_Day2_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Female_Day2_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Female_Day2_meta$Group, fill=Female_Day2_meta$Group, shape=Female_Day2_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Female_Day2_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

#Female_Day5_PCoA
dist_bray <- capscale(Female_Day5_shared~1, distance="bray", binary=FALSE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Female_Day5_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Female_Day5_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Female_Day5_meta$Group, fill=Female_Day5_meta$Group, shape=Female_Day5_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Female_Day5_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

#Female_Day8_PCoA
dist_bray <- capscale(Female_Day8_shared~1, distance="bray", binary=FALSE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Female_Day8_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Female_Day8_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Female_Day8_meta$Group, fill=Female_Day8_meta$Group, shape=Female_Day8_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Female_Day8_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

#Female_Day12_PCoA
dist_bray <- capscale(Female_Day12_shared~1, distance="bray", binary=FALSE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Female_Day12_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Female_Day12_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Female_Day12_meta$Group, fill=Female_Day12_meta$Group, shape=Female_Day12_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Female_Day12_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

#Female_Day65_PCoA
dist_bray <- capscale(Female_Day65_shared~1, distance="bray", binary=FALSE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Female_Day65_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Female_Day65_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Female_Day65_meta$Group, fill=Female_Day65_meta$Group, shape=Female_Day65_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Female_Day65_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

#Female_Day100_PCoA
dist_bray <- capscale(Female_Day100_shared~1, distance="bray", binary=FALSE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Female_Day100_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Female_Day100_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Female_Day100_meta$Group, fill=Female_Day100_meta$Group, shape=Female_Day100_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Female_Day100_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))








#Male_Day1_PCoA
dist_bray <- capscale(Male_Day1_shared~1, distance="bray", binary=FALSE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Male_Day1_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Male_Day1_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Male_Day1_meta$Group, fill=Male_Day1_meta$Group, shape=Male_Day1_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Male_Day1_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

#Male_Day2_PCoA
dist_bray <- capscale(Male_Day2_shared~1, distance="bray", binary=FALSE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Male_Day2_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Male_Day2_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Male_Day2_meta$Group, fill=Male_Day2_meta$Group, shape=Male_Day2_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Male_Day2_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

#Male_Day5_PCoA
dist_bray <- capscale(Male_Day5_shared~1, distance="bray", binary=FALSE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Male_Day5_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Male_Day5_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Male_Day5_meta$Group, fill=Male_Day5_meta$Group, shape=Male_Day5_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Male_Day5_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

#Male_Day8_PCoA
dist_bray <- capscale(Male_Day8_shared~1, distance="bray", binary=FALSE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Male_Day8_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Male_Day8_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Male_Day8_meta$Group, fill=Male_Day8_meta$Group, shape=Male_Day8_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Male_Day8_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))


#Male_Day12_PCoA
dist_bray <- capscale(Male_Day12_shared~1, distance="bray", binary=FALSE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Male_Day12_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Male_Day12_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Male_Day12_meta$Group, fill=Male_Day12_meta$Group, shape=Male_Day12_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Male_Day12_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

#Male_Day65_PCoA
dist_bray <- capscale(Male_Day65_shared~1, distance="bray", binary=FALSE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Male_Day65_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Male_Day65_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Male_Day65_meta$Group, fill=Male_Day65_meta$Group, shape=Male_Day65_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Male_Day65_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

#Male_Day100_PCoA
dist_bray <- capscale(Male_Day100_shared~1, distance="bray", binary=FALSE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = Male_Day100_shared)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Male_Day100_bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = Male_Day100_meta$Group, fill=Male_Day100_meta$Group, shape=Male_Day100_meta$Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = Male_Day100_meta$Group), size=1.5, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=5, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))





library(ggpubr)

Female_bray_merged <- ggarrange(Female_Day1_bray_plot,
                                Female_Day2_bray_plot,
                                Female_Day5_bray_plot,
                                Female_Day8_bray_plot,
                                Female_Day12_bray_plot,
                                Female_Day65_bray_plot,
                                Female_Day100_bray_plot,
                                common.legend = TRUE, legend = "right",
                                ncol = 1, nrow = 7, align="hv")



Male_bray_merged <- ggarrange(Male_Day1_bray_plot,
                              Male_Day2_bray_plot,
                              Male_Day5_bray_plot,
                              Male_Day8_bray_plot,
                              Male_Day12_bray_plot,
                              Male_Day65_bray_plot,
                              Male_Day100_bray_plot,
                              common.legend = TRUE, legend = "none",
                              ncol = 1, nrow = 7, align="hv")



bray_merged <- ggarrange(Female_bray_merged ,
                         Male_bray_merged,
                         common.legend = TRUE, legend = "right",
                         ncol = 2, nrow = 1, align="hv")













library(pairwiseAdonis)



x <- phyloseq::distance(merged_phy_percent1, method ="bray") 
factors <- meta_phy_out1$Cohort
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="none"))

dist_bray1 <- vegdist(shared_phy_out1, method="bray", binary=FALSE)
adonis2(formula = dist_bray1 ~ Cohort*Sex*Day, data = meta_phy_out1, permutations=1000)





x <- phyloseq::distance(Female_Day1_ps, method ="bray") 
factors <- Female_Day1_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Female_Day2_ps, method ="bray") 
factors <- Female_Day2_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Female_Day5_ps, method ="bray") 
factors <- Female_Day5_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Female_Day8_ps, method ="bray") 
factors <- Female_Day8_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Female_Day12_ps, method ="bray") 
factors <- Female_Day12_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Female_Day65_ps, method ="bray") 
factors <- Female_Day2_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Female_Day100_ps, method ="bray") 
factors <- Female_Day2_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))





x <- phyloseq::distance(Male_Day1_ps, method ="bray", binary=FALSE) 
factors <- Male_Day1_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Male_Day2_ps, method ="bray") 
factors <- Male_Day2_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Male_Day5_ps, method ="bray") 
factors <- Male_Day5_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Male_Day8_ps, method ="bray") 
factors <- Male_Day8_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Male_Day12_ps, method ="bray") 
factors <- Male_Day12_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Male_Day65_ps, method ="bray") 
factors <- Male_Day2_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Male_Day100_ps, method ="bray") 
factors <- Male_Day2_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))





x <- phyloseq::distance(Female_Day1_ps, method ="bray") 
factors <- Female_Day1_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Female_Day2_ps, method ="bray") 
factors <- Female_Day2_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Female_Day5_ps, method ="bray") 
factors <- Female_Day5_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Female_Day8_ps, method ="bray") 
factors <- Female_Day8_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Female_Day12_ps, method ="bray") 
factors <- Female_Day12_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Female_Day65_ps, method ="bray") 
factors <- Female_Day2_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Female_Day100_ps, method ="bray") 
factors <- Female_Day2_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))








x <- phyloseq::distance(Male_Day1_ps, method ="bray") 
factors <- Male_Day1_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Male_Day2_ps, method ="bray") 
factors <- Male_Day2_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Male_Day5_ps, method ="bray") 
factors <- Male_Day5_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Male_Day8_ps, method ="bray") 
factors <- Male_Day8_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Male_Day12_ps, method ="bray") 
factors <- Male_Day12_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Male_Day65_ps, method ="bray") 
factors <- Male_Day2_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))

x <- phyloseq::distance(Male_Day100_ps, method ="bray") 
factors <- Male_Day2_meta$Group
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))










dev.off()


Female_ps <- subset_samples(merged_phy_percent1, Sex=="Female")
Female_ps_C <- subset_samples(Female_ps, Group=="C")
Female_ps_GA <- subset_samples(Female_ps, Group=="GA")
Female_ps_IA <- subset_samples(Female_ps, Group=="IA")


Male_ps <- subset_samples(merged_phy_percent1, Sex=="Male")
Male_ps_C <- subset_samples(Male_ps, Group=="C")
Male_ps_GA <- subset_samples(Male_ps, Group=="GA")
Male_ps_IA <- subset_samples(Male_ps, Group=="IA")



#merged_phy_percent2 <- phyloseq_standardize_otu_abundance(merged_phy_percent , method = "total")


#SI_pruned <- prune_taxa(names(sort(taxa_sums(merged_phy_percent2),TRUE)[1:5]), merged_phy_percent2)
#sampleOrder = unique(sample_names(SI_pruned))
#taxaOrder = rev(unique(taxa_names(SI_pruned)))


Female_ps_C_plot <- plot_bar(Female_ps_C, fill="Rank1") +
  facet_wrap(~DayPlot2, scales = "free_x", nrow=1) +
  theme(axis.text = element_text(size=2, color="black"),
        #strip.text.x = element_text(size=5, color="black"),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        axis.text.x= element_text(color="black", size=1),
        axis.text.y= element_text(color="black", size=5),
        text = element_text(size=5, face=NULL),
        legend.text=element_text(size=6),
        axis.title.x = element_blank(),
        legend.position="right", legend.box = "vertical")
        
Female_ps_GA_plot <- plot_bar(Female_ps_GA, fill="Rank6") +
  facet_wrap(~DayPlot2, scales = "free_x", nrow=1) +
  theme(axis.text = element_text(size=2, color="black"),
        #strip.text.x = element_text(size=5, color="black"),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        axis.text.x= element_text(color="black", size=1),
        axis.text.y= element_text(color="black", size=5),
        text = element_text(size=5, face=NULL),
        legend.text=element_text(size=6),
        axis.title.x = element_blank(),
        legend.position="right", legend.box = "vertical")

Female_ps_IA_plot <- plot_bar(Female_ps_IA, fill="Rank6") +
  facet_wrap(~DayPlot2, scales = "free_x", nrow=1) +
  theme(axis.text = element_text(size=2, color="black"),
        #strip.text.x = element_text(size=5, color="black"),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        axis.text.x= element_text(color="black", size=1),
        axis.text.y= element_text(color="black", size=5),
        text = element_text(size=5, face=NULL),
        legend.text=element_text(size=6),
        axis.title.x = element_blank(),
        legend.position="right", legend.box = "vertical")

library(ggpubr)

Female_taxa_plot_merged <- ggarrange(Female_ps_C_plot, Female_ps_GA_plot, Female_ps_IA_plot,
                         common.legend = TRUE, legend = "bottom",
                         ncol = 1, nrow = 3, align="hv")










Male_ps_C_plot <- plot_bar(Male_ps_C, fill="Rank6") +
  facet_wrap(~DayPlot2, scales = "free_x", nrow=1) +
  theme(axis.text = element_text(size=2, color="black"),
        #strip.text.x = element_text(size=5, color="black"),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        axis.text.x= element_text(color="black", size=1),
        axis.text.y= element_text(color="black", size=5),
        text = element_text(size=5, face=NULL),
        legend.text=element_text(size=6),
        axis.title.x = element_blank(),
        legend.position="right", legend.box = "vertical")

Male_ps_GA_plot <- plot_bar(Male_ps_GA, fill="Rank6") +
  facet_wrap(~DayPlot2, scales = "free_x", nrow=1) +
  theme(axis.text = element_text(size=2, color="black"),
        #strip.text.x = element_text(size=5, color="black"),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        axis.text.x= element_text(color="black", size=1),
        axis.text.y= element_text(color="black", size=5),
        text = element_text(size=5, face=NULL),
        legend.text=element_text(size=6),
        axis.title.x = element_blank(),
        legend.position="right", legend.box = "vertical")

Male_ps_IA_plot <- plot_bar(Male_ps_IA, fill="Rank6") +
  facet_wrap(~DayPlot2, scales = "free_x", nrow=1) +
  theme(axis.text = element_text(size=2, color="black"),
        #strip.text.x = element_text(size=5, color="black"),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        axis.text.x= element_text(color="black", size=1),
        axis.text.y= element_text(color="black", size=5),
        text = element_text(size=5, face=NULL),
        legend.text=element_text(size=6),
        axis.title.x = element_blank(),
        legend.position="right", legend.box = "vertical")



library(ggpubr)

Male_taxa_plot_merged <- ggarrange(Male_ps_C_plot, Male_ps_GA_plot, Male_ps_IA_plot,
                                   common.legend = TRUE, legend = "bottom",
                                   ncol = 1, nrow = 3, align="hv")































































#Differential abundance



library(phyloseq)
library(metagMisc)


setwd("D://ABX_LOCO_2024/")






setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/ABX_LOCO_meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/ABX_LOCO_taxa.txt"))


#meta_ob1 <- meta_ob1[meta_ob1$Outliers=="No",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)
#fullMothur<-read.csv("PWY_shared.csv",)
fullMothur<-read.csv("ABX_LOCO_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtu))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))


phyloseqobj.f <- tax_glom(phyloseqobj.f, taxrank="Rank6")

merged_phy_percent1 <- phyloseqobj.f 


#colnames(tax_table(phyloseqobj.f)) <- c("kingdom", "phylum", "class", "order", "family",  "genus")



#merged_phy_percent1 <- phyloseq_standardize_otu_abundance(phyloseqobj.f , method = "total")


#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
#merged_phy_percent1 <- tax_glom(pruned_ps, taxrank="Rank6")

rank_names(merged_phy_percent1)

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


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)



#ps <- psmelt(merged_phy_percent1)



ps <- plot_bar(merged_phy_percent1  , fill="Rank1") + facet_wrap(~Sex, scales = "free_x") +
  theme(legend.position = "bottom")

plot_data <- ggplot_build(ps)$plot$data
write.csv(plot_data, file="shared_phy_out1.csv")












#beta diversity

#devtools::install_github("vmikk/metagMisc")
library(phyloseq)
library(metagMisc)
library(vegan)
library(ggplot2)









###Maaslin2 Merged Sex Rank1

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa_138.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Group!="GA" & meta_ob1$DayPlot=="B2",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names

fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank2")

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
tax_table(merged_phy_percent1)[,2]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,2]

meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)

setwd("D://ABX_LOCO_2024/Merged_Sex_Results/")

library(Maaslin2)

fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    min_prevalence=0.5,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "test_phyla_Group_Sex_B2", 
                    fixed_effects  = c("Group", "Sex"),
                    reference      = c("Group,C", "Sex,Female"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)





setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa_138.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Group!="GA" & meta_ob1$DayPlot=="C5",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names

fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank2")

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
tax_table(merged_phy_percent1)[,2]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,2]

meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)

setwd("D://ABX_LOCO_2024/Merged_Sex_Results/")

library(Maaslin2)

fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    min_prevalence=0.5,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "phyla_Group_Sex_C5", 
                    fixed_effects  = c("Group", "Sex"),
                    reference      = c("Group,C", "Sex,Female"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)








setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa_138.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Group!="GA" & meta_ob1$DayPlot=="D8",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names

fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank2")

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
tax_table(merged_phy_percent1)[,2]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,2]

meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)

setwd("D://ABX_LOCO_2024/Merged_Sex_Results/")

library(Maaslin2)

fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    min_prevalence=0.5,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "phyla_Group_Sex_D8", 
                    fixed_effects  = c("Group", "Sex"),
                    reference      = c("Group,C", "Sex,Female"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)








setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa_138.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Group!="GA" & meta_ob1$DayPlot=="E12",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names

fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank2")

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
tax_table(merged_phy_percent1)[,2]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,2]

meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)

setwd("D://ABX_LOCO_2024/Merged_Sex_Results/")

library(Maaslin2)

fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    min_prevalence=0.5,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "phyla_Group_Sex_E12", 
                    fixed_effects  = c("Group", "Sex"),
                    reference      = c("Group,C", "Sex,Female"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)








setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa_138.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Group!="GA" & meta_ob1$DayPlot=="G65",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names

fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank2")

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
tax_table(merged_phy_percent1)[,2]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,2]

meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)

setwd("D://ABX_LOCO_2024/Merged_Sex_Results/")

library(Maaslin2)

fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    min_prevalence=0.5,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "phyla_Group_Sex_G65", 
                    fixed_effects  = c("Group", "Sex"),
                    reference      = c("Group,C", "Sex,Female"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)








setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa_138.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Group!="GA" & meta_ob1$DayPlot=="H100",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names

fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank2")

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
tax_table(merged_phy_percent1)[,2]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,2]

meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)

setwd("D://ABX_LOCO_2024/Merged_Sex_Results/")

library(Maaslin2)

fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    min_prevalence=0.5,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "phyla_Group_Sex_H100", 
                    fixed_effects  = c("Group", "Sex"),
                    reference      = c("Group,C", "Sex,Female"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)




































###Maaslin2 Merged Sex Rank6

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa_138_Rank6.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Group!="GA" & meta_ob1$DayPlot=="B2",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names

fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank2")

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
tax_table(merged_phy_percent1)[,2]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,2]

meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)

setwd("D://ABX_LOCO_2024/Merged_Sex_Results/")

library(Maaslin2)

fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    min_prevalence=0.5,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Rank6_Group_Sex_B2", 
                    fixed_effects  = c("Group", "Sex"),
                    reference      = c("Group,C", "Sex,Female"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)





setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa_138_Rank6.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Group!="GA" & meta_ob1$DayPlot=="C5",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names

fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank2")

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
tax_table(merged_phy_percent1)[,2]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,2]

meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)

setwd("D://ABX_LOCO_2024/Merged_Sex_Results/")

library(Maaslin2)

fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    min_prevalence=0.5,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Rank6_Group_Sex_C5", 
                    fixed_effects  = c("Group", "Sex"),
                    reference      = c("Group,C", "Sex,Female"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)








setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa_138_Rank6.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Group!="GA" & meta_ob1$DayPlot=="D8",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names

fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank2")

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
tax_table(merged_phy_percent1)[,2]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,2]

meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)

setwd("D://ABX_LOCO_2024/Merged_Sex_Results/")

library(Maaslin2)

fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    min_prevalence=0.5,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Rank6_Group_Sex_D8", 
                    fixed_effects  = c("Group", "Sex"),
                    reference      = c("Group,C", "Sex,Female"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)








setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa_138_Rank6.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Group!="GA" & meta_ob1$DayPlot=="E12",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names

fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank2")

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
tax_table(merged_phy_percent1)[,2]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,2]

meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)

setwd("D://ABX_LOCO_2024/Merged_Sex_Results/")

library(Maaslin2)

fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    min_prevalence=0.5,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Rank6_Group_Sex_E12", 
                    fixed_effects  = c("Group", "Sex"),
                    reference      = c("Group,C", "Sex,Female"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)








setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa_138_Rank6.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Group!="GA" & meta_ob1$DayPlot=="G65",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names

fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank2")

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
tax_table(merged_phy_percent1)[,2]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,2]

meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)

setwd("D://ABX_LOCO_2024/Merged_Sex_Results/")

library(Maaslin2)

fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    min_prevalence=0.5,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Rank6_Group_Sex_G65", 
                    fixed_effects  = c("Group", "Sex"),
                    reference      = c("Group,C", "Sex,Female"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)








setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa_138_Rank6.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Group!="GA" & meta_ob1$DayPlot=="H100",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names

fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank2")

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
tax_table(merged_phy_percent1)[,2]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,2]

meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)

setwd("D://ABX_LOCO_2024/Merged_Sex_Results/")

library(Maaslin2)

fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    min_prevalence=0.5,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Rank6_Group_Sex_H100", 
                    fixed_effects  = c("Group", "Sex"),
                    reference      = c("Group,C", "Sex,Female"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)


























###Male controls



setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa.txt"))



meta_ob1 <- meta_ob1[meta_ob1$Sex=="Male" & meta_ob1$Control_Group!="Other",]



nrow(meta_ob1)



GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)
#fullMothur<-read.csv("CEF_PWY_shared.csv",)
fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

#merged_phy_percent1 <- phyloseqobj.f 


#merged_phy_percent1 <- phyloseq_standardize_otu_abundance(phyloseqobj.f, method = "total")


merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank6")


#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
#merged_phy_percent1 <- tax_glom(pruned_ps, taxrank="Rank6")


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
tax_table(merged_phy_percent1)[,6]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)






library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Male_Controls_NEGBIN", 
                    fixed_effects  = c("DayPlot"),
                    reference      = c("DayPlot,A1"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)










###Female Plot controls

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa.txt"))


meta_ob1 <- meta_ob1[meta_ob1$Sex=="Female" & meta_ob1$Group=="C",]


nrow(meta_ob1)


GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))


merged_phy_percent1 <- phyloseq_standardize_otu_abundance(phyloseqobj.f, method = "total")


merged_phy_percent1 <- tax_glom(merged_phy_percent1, taxrank="Rank6")



#merged_phy_percent1 <- prune_taxa(names(sort(taxa_sums(merged_phy_percent1),TRUE)[1:25]), merged_phy_percent1)


#merged_phy_percent1 <- merge_samples(merged_phy_percent1, meta_ob1$Control_Group, fun=mean)


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
tax_table(merged_phy_percent1)[,6]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)


ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent1),TRUE)[1:25]), merged_phy_percent1 )
sampleOrder = unique(sample_names(ps))
taxaOrder = rev(unique(taxa_names(ps)))








#write.csv(shared_phy_out1, file="top_25.csv")




Female_Control_bar_plot <- plot_bar(ps, fill="Rank6") + facet_wrap(~Control_Group, scales = "free_x", ncol = 7, nrow = 1) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = 4)) +
  theme_minimal() +
  theme(legend.position="bottom", legend.box = "bottom",
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size=10, face=NULL))



x <- phyloseq::distance(merged_phy_percent1, method ="bray") 
factors <- meta_phy_out1$Control_Group
pairwise.adonis(x, factors, perm = 1000, p.adjust(method="FDR"))




dist_bray <- capscale(shared_phy_out1~1, distance="bray", binary=FALSE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = shared_phy_out1)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Female_C_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = input_metadata$Control_Group, fill=input_metadata$Control_Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  #scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = input_metadata$Control_Group), size=3, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=8, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines')) +
  guides(colour = guide_legend(nrow = 1))
  #annotate("text", colour = "red", x =-0.777329209, y =0.221109825, label = "Microbacteriaceae") +
  #annotate("text", colour = "red", x = -0.553126116, y = -0.083112784, label = "Akkermansia") +
  #annotate("text", colour = "red", x = 0.565300189, y = -0.42666947, label = "Lachnospiraceae_UCG-010") +
  #annotate("text", colour = "red", x = 0.489226656, y = 	0.343370012, label = "Enterococcaceae") +
  #annotate("text", colour = "red", x = 0.095690727, y = 	-0.78777082, label = "Prevotellaceae") +
  #annotate("text", colour = "red", x = 0.318129332, y = 	0.305375509, label = "Clostridium sensu stricto 1")










###Male Plot controls

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa.txt"))



meta_ob1 <- meta_ob1[meta_ob1$Sex=="Male" & meta_ob1$Group=="C",]

nrow(meta_ob1)


GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))


merged_phy_percent1 <- phyloseq_standardize_otu_abundance(phyloseqobj.f, method = "total")


#merged_phy_percent1 <- tax_glom(merged_phy_percent1, taxrank="Rank6")


#merged_phy_percent1 <- merge_samples(merged_phy_percent1, meta_ob1$Control_Group, fun=mean)


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
tax_table(merged_phy_percent1)[,6]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)


ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent1),TRUE)[1:25]), merged_phy_percent1 )
sampleOrder = unique(sample_names(ps))
taxaOrder = rev(unique(taxa_names(ps)))


Male_Control_bar_plot <- plot_bar(ps, fill="Rank6") + facet_wrap(~Control_Group, scales = "free_x", ncol = 7, nrow = 1) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = 4)) +
  theme_minimal() +
  theme(legend.position="bottom", legend.box = "bottom",
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size=10, face=NULL))




x <- phyloseq::distance(merged_phy_percent1, method ="bray") 
factors <- meta_phy_out1$Control_Group
pairwise.adonis(x, factors, perm = 1000, p.adjust(method="FDR"))




dist_bray <- capscale(shared_phy_out1~1, distance="bray", binary=FALSE, eig=TRUE)
prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)
data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)
df <- data.frame(data.scores)
library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = shared_phy_out1)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)
Male_C_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = input_metadata$Control_Group, fill=input_metadata$Control_Group)) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  #scale_shape_manual(values = c(1, 0, 2)) +
  geom_point(mapping = aes(colour = input_metadata$Control_Group), size=3, alpha=1.0, stroke=0.5) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=8, face=NULL), axis.text.y = element_text(color="black"),
        legend.title= element_blank(),
        strip.text.x = element_blank()) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines')) +
  guides(colour = guide_legend(nrow = 1))















library(ggpubr)

abundance_plot_merged <- ggarrange(Female_Control_bar_plot,
                                   Male_Control_bar_plot,
                                   common.legend = TRUE, legend = "bottom",
                                   ncol = 2, nrow = 1, align="hv")


PCoA_plot_merged <- ggarrange(Female_C_plot,
                              Male_C_plot,
                                   common.legend = TRUE, legend = "top",
                                   ncol = 2, nrow = 1, align="hv")


abundance_PCoA_plot_merged <- ggarrange(PCoA_plot_merged,
                                        abundance_plot_merged,
                                        common.legend = TRUE, legend = "top",
                                        ncol = 1, nrow = 2, align="hv")



###Rank6 Controls Sex adonis


setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa.txt"))


meta_ob1 <- meta_ob1[meta_ob1$Control_Group!="Other",]


nrow(meta_ob1)


GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))


merged_phy_percent1 <- phyloseq_standardize_otu_abundance(phyloseqobj.f, method = "total")


merged_phy_percent1 <- tax_glom(merged_phy_percent1, taxrank="Rank6")


#merged_phy_percent1 <- merge_samples(merged_phy_percent1, meta_ob1$Control_Group, fun=mean)


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
tax_table(merged_phy_percent1)[,6]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)







library(pairwiseAdonis)


dist.bray <- vegdist(shared_phy_out1, method="bray", binary=FALSE)
adonis2(formula = dist.bray  ~ Sex + Control_Group, data = meta_phy_out1, permutations=1000)


x <- phyloseq::distance(merged_phy_percent1, method ="bray") 
factors <- meta_phy_out1$DayPlot
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))











###Rank6 maaslin




###Female Day A1

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa.txt"))



#meta_ob1 <- meta_ob1[meta_ob1$Sex=="Female" & meta_ob1$DayPlot=="A1",]
meta_ob1 <- meta_ob1[meta_ob1$DayPlot=="C5",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)
#fullMothur<-read.csv("CEF_PWY_shared.csv",)
fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

#merged_phy_percent1 <- phyloseqobj.f 


#merged_phy_percent1 <- phyloseq_standardize_otu_abundance(phyloseqobj.f, method = "total")


merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank6")

#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
#merged_phy_percent1 <- tax_glom(pruned_ps, taxrank="Rank6")


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
tax_table(merged_phy_percent1)[,6]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)






library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Female_NEGBIN_A1", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)







###Female Day B2

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Female" & meta_ob1$DayPlot=="B2",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)
#fullMothur<-read.csv("CEF_PWY_shared.csv",)
fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

#merged_phy_percent1 <- phyloseqobj.f 

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank6")

#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
#merged_phy_percent1 <- tax_glom(pruned_ps, taxrank="Rank6")



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
tax_table(merged_phy_percent1)[,6]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Female_NEGBIN_B2", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)







###Female Day C5

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Female" & meta_ob1$DayPlot=="C5",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)
#fullMothur<-read.csv("PWY_shared.csv",)
fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

#merged_phy_percent1 <- phyloseqobj.f 

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank6")

#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
#merged_phy_percent1 <- tax_glom(pruned_ps, taxrank="Rank6")


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
tax_table(merged_phy_percent1)[,6]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Female_NEGBIN_C5", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)







###Female Day D8

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Female" & meta_ob1$DayPlot=="D8",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)
#fullMothur<-read.csv("CEF_PWY_shared.csv",)
fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

#merged_phy_percent1 <- phyloseqobj.f 

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank6")

#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
#merged_phy_percent1 <- tax_glom(pruned_ps, taxrank="Rank6")


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
tax_table(merged_phy_percent1)[,6]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Female_NEGBIN_D8", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)







###Female Day E12

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Female" & meta_ob1$DayPlot=="E12",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)
#fullMothur<-read.csv("CEF_PWY_shared.csv",)
fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

#merged_phy_percent1 <- phyloseqobj.f 

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank6")

#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
#merged_phy_percent1 <- tax_glom(pruned_ps, taxrank="Rank6")


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
tax_table(merged_phy_percent1)[,6]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Female_NEGBIN_E12", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)






###Female Day G65

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Female" & meta_ob1$DayPlot=="G65",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)
#fullMothur<-read.csv("CEF_PWY_shared.csv",)
fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

#merged_phy_percent1 <- phyloseqobj.f 

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank6")

#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
#merged_phy_percent1 <- tax_glom(pruned_ps, taxrank="Rank6")


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
tax_table(merged_phy_percent1)[,6]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Female_NEGBIN_G65", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)






###Female Day H100

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Female" & meta_ob1$DayPlot=="H100",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)
#fullMothur<-read.csv("CEF_PWY_shared.csv",)
fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

#merged_phy_percent1 <- phyloseqobj.f 

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank6")

#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
#merged_phy_percent1 <- tax_glom(pruned_ps, taxrank="Rank6")


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
tax_table(merged_phy_percent1)[,6]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Female_NEGBIN_H100", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)













###Male Day A1

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Male" & meta_ob1$DayPlot=="A1",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)
#fullMothur<-read.csv("CEF_PWY_shared.csv",)
fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

#merged_phy_percent1 <- phyloseqobj.f 

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank6")

#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
#merged_phy_percent1 <- tax_glom(pruned_ps, taxrank="Rank6")


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
tax_table(merged_phy_percent1)[,6]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Male_NEGBIN_A1", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)







###Male Day B2

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Male" & meta_ob1$DayPlot=="B2",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)
#fullMothur<-read.csv("CEF_PWY_shared.csv",)
fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

#merged_phy_percent1 <- phyloseqobj.f 

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank6")

#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
#merged_phy_percent1 <- tax_glom(pruned_ps, taxrank="Rank6")


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
tax_table(merged_phy_percent1)[,6]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)



library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Male_NEGBIN_B2", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)







###Male Day C5

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Male" & meta_ob1$DayPlot=="C5",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)
#fullMothur<-read.csv("CEF_PWY_shared.csv",)
fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

#merged_phy_percent1 <- phyloseqobj.f 

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank6")

#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
#merged_phy_percent1 <- tax_glom(pruned_ps, taxrank="Rank6")


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
tax_table(merged_phy_percent1)[,6]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Male_NEGBIN_C5", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)







###Male Day D8

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Male" & meta_ob1$DayPlot=="D8",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)
#fullMothur<-read.csv("CEF_PWY_shared.csv",)
fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

#merged_phy_percent1 <- phyloseqobj.f 

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank6")

#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
#merged_phy_percent1 <- tax_glom(pruned_ps, taxrank="Rank6")


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
tax_table(merged_phy_percent1)[,6]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Male_NEGBIN_D8", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)







###Male Day E12

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Male" & meta_ob1$DayPlot=="E12",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)
#fullMothur<-read.csv("CEF_PWY_shared.csv",)
fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

#merged_phy_percent1 <- phyloseqobj.f 

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank6")

#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
#merged_phy_percent1 <- tax_glom(pruned_ps, taxrank="Rank6")


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
tax_table(merged_phy_percent1)[,6]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)



library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Male_NEGBIN_E12", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)






###Male Day G65

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Male" & meta_ob1$DayPlot=="G65",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)
#fullMothur<-read.csv("CEF_PWY_shared.csv",)
fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

#merged_phy_percent1 <- phyloseqobj.f 

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank6")

#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
#merged_phy_percent1 <- tax_glom(pruned_ps, taxrank="Rank6")


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
tax_table(merged_phy_percent1)[,6]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Male_NEGBIN_G65", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)






###Male Day H100

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Male" & meta_ob1$DayPlot=="H100",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)
#fullMothur<-read.csv("CEF_PWY_shared.csv",)
fullMothur<-read.csv("CEF_16S_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

#merged_phy_percent1 <- phyloseqobj.f 

merged_phy_percent1 <- tax_glom(phyloseqobj.f, taxrank="Rank6")

#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
#merged_phy_percent1 <- tax_glom(pruned_ps, taxrank="Rank6")


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
tax_table(merged_phy_percent1)[,6]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)






library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Male_NEGBIN_H100", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)
















library("ggpicrust2")



setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)


meta_ob1 <- meta_ob1[meta_ob1$Sex=="Female",]


nrow(meta_ob1)


GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("pred_metagenome_unstrat.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 





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


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1)) 


row.names(meta_phy_out1) == colnames(t(shared_phy_out1))

write.table(meta_phy_out1, file='D://ABX_LOCO_2024/F_8_meta.txt', sep='\t', col.names = NA)
write.table(t(shared_phy_out1), file='D://ABX_LOCO_2024/F_8_pred_metagenome_unstrat.tsv', sep='\t', col.names = NA)



library(tidyverse)

metadata <-
  read_delim(
    "D://ABX_LOCO_2024/F_8_meta.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )
group <- "Enviroment"
daa_results_list <-
  ggpicrust2(
    file = "D://ABX_LOCO_2024/F_8_pred_metagenome_unstrat.tsv",
    metadata = metadata,
    group = "Group",
    pathway = "EC",
    daa_method = "Maaslin2",
    order = "Group",
    ko_to_kegg = FALSE,
    x_lab = "description",
    p.adjust = "BH",
    select = NULL,
    reference ="Group,C")













###PWY maaslin



#beta diversity

#devtools::install_github("vmikk/metagMisc")
library(phyloseq)
library(metagMisc)
library(vegan)
library(ggplot2)




###Day A1

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Group!="GA" & meta_ob1$DayPlot=="A1",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 


#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")


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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)



library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Merged_Sex_Day1", 
                    fixed_effects  = c("Sex", "Group"),
                    reference      = c("Sex,Female", "Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)









###Day B2

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Group!="GA" & meta_ob1$DayPlot=="B2",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 


#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")


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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)



library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Merged_Sex_Day2", 
                    fixed_effects  = c("Sex", "Group"),
                    reference      = c("Sex,Female", "Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)













###Day C5

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Group!="GA" & meta_ob1$DayPlot=="C5",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 


#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")


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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)



library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Merged_Sex_Day5", 
                    fixed_effects  = c("Sex", "Group"),
                    reference      = c("Sex,Female", "Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)









###Day D8

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Group!="GA" & meta_ob1$DayPlot=="D8",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 


#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")


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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)



library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Merged_Sex_Day8", 
                    fixed_effects  = c("Sex", "Group"),
                    reference      = c("Sex,Female", "Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)











###Day E12

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Group!="GA" & meta_ob1$DayPlot=="E12",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 


#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")


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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)



library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Merged_Sex_Day12", 
                    fixed_effects  = c("Sex", "Group"),
                    reference      = c("Sex,Female", "Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)










###Day G65

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Group!="GA" & meta_ob1$DayPlot=="G65",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 


#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")


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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)



library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Merged_Sex_Day65", 
                    fixed_effects  = c("Sex", "Group"),
                    reference      = c("Sex,Female", "Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)









###Day H100

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Group!="GA" & meta_ob1$DayPlot=="H100",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 


#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")


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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)



library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Merged_Sex_Day100", 
                    fixed_effects  = c("Sex", "Group"),
                    reference      = c("Sex,Female", "Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)





###PWY plots



###PWY plots

setwd("D://ABX_LOCO_2024/")

data_df <- read.table("D://ABX_LOCO_2024/Merged_Sex_Results/PWY_output/PWY_merged_Sex.txt", header=T, row.names=1, check.names = FALSE)

library(dplyr)

#F_data_df <- subset(data_df, abs(FC) > 2.0 & qval<0.01 & Sex=="Female")
#M_data_df <- subset(data_df, abs(FC) > 2.0 & qval<0.01 & Sex=="Male")


data_df <- subset(data_df, N.not.0 > 16)



###Split days
Day2_PWY <- subset(data_df, Day=="2")
Day5_PWY <- subset(data_df, Day=="5")
Day8_PWY <- subset(data_df, Day=="8")
Day12_PWY <- subset(data_df, Day=="12")
Day65_PWY <- subset(data_df, Day=="65")
Day100_PWY <- subset(data_df, Day=="100")


library(EnhancedVolcano)









keyvals.shape2 <- Female_Day2_PWY$shape
names(keyvals.shape2)[keyvals.shape2 == 1] <- 'GA'
#names(keyvals.shape2)[keyvals.shape2 == 2] <- 'IA'

keyvals.shape5 <- Female_Day5_PWY$shape
names(keyvals.shape5)[keyvals.shape5 == 1] <- 'GA'
#names(keyvals.shape5)[keyvals.shape5 == 2] <- 'IA'

keyvals.shape8 <- Female_Day8_PWY$shape
names(keyvals.shape8)[keyvals.shape8 == 1] <- 'GA'
#names(keyvals.shape8)[keyvals.shape8 == 2] <- 'IA'

keyvals.shape12 <- Female_Day12_PWY$shape
names(keyvals.shape12)[keyvals.shape12 == 1] <- 'GA'
#names(keyvals.shape12)[keyvals.shape12 == 2] <- 'IA'

keyvals.shape65 <- Female_Day65_PWY$shape
names(keyvals.shape65)[keyvals.shape65 == 1] <- 'GA'
#names(keyvals.shape65)[keyvals.shape65 == 2] <- 'IA'

keyvals.shape100 <- Female_Day100_PWY$shape
names(keyvals.shape100)[keyvals.shape100 == 1] <- 'GA'
#names(keyvals.shape100)[keyvals.shape100 == 2] <- 'IA'







Day2_PWY_plot <- EnhancedVolcano(Day2_PWY,
                                 lab=Day2_PWY$feature,x='FC',y='qval',
                                 colAlpha = 1,
                                 title = NULL,
                                 subtitle = NULL,
                                 pointSize = 3.0,
                                 labSize = 3,
                                 drawConnectors = TRUE,
                                 typeConnectors='closed',
                                 endsConnectors='last',
                                 legendPosition = 'blank',
                                 pCutoff = 0.000001,
                                 FCcutoff = 2.0,
                                 axisLabSize=4,
                                 caption=NULL,
                                 #shapeCustom = F_keyvals.shape2,
                                 max.overlaps = 100) +
  #coord_cartesian(xlim=c(-6,8), ylim=c(0,155)) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))


Day5_PWY_plot <- EnhancedVolcano(Day5_PWY,
                                 lab=Day5_PWY$feature,x='FC',y='qval',
                                 colAlpha = 1,
                                 title = NULL,
                                 subtitle = NULL,
                                 pointSize = 3.0,
                                 labSize = 3,
                                 drawConnectors = TRUE,
                                 typeConnectors='closed',
                                 endsConnectors='last',
                                 legendPosition = 'blank',
                                 pCutoff = 0.000001,
                                 FCcutoff = 2.0,
                                 axisLabSize=4,
                                 caption=NULL,
                                 #shapeCustom = F_keyvals.shape2,
                                 max.overlaps = 100) +
  #coord_cartesian(xlim=c(-6,8), ylim=c(0,170)) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))




Day8_PWY_plot <- EnhancedVolcano(Day8_PWY,
                                 lab=Day8_PWY$feature,x='FC',y='qval',
                                 colAlpha = 1,
                                 title = NULL,
                                 subtitle = NULL,
                                 pointSize = 3.0,
                                 labSize = 3,
                                 drawConnectors = TRUE,
                                 typeConnectors='closed',
                                 endsConnectors='last',
                                 legendPosition = 'blank',
                                 pCutoff = 0.000001,
                                 FCcutoff = 2.0,
                                 axisLabSize=4,
                                 caption=NULL,
                                 #shapeCustom = F_keyvals.shape2,
                                 max.overlaps = 100) +
  #coord_cartesian(xlim=c(-6,8), ylim=c(0,155)) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))



Day12_PWY_plot <- EnhancedVolcano(Day12_PWY,
                                 lab=Day12_PWY$feature,x='FC',y='qval',
                                 colAlpha = 1,
                                 title = NULL,
                                 subtitle = NULL,
                                 pointSize = 3.0,
                                 labSize = 3,
                                 drawConnectors = TRUE,
                                 typeConnectors='closed',
                                 endsConnectors='last',
                                 legendPosition = 'blank',
                                 pCutoff = 0.000001,
                                 FCcutoff = 2.0,
                                 axisLabSize=4,
                                 caption=NULL,
                                 #shapeCustom = F_keyvals.shape2,
                                 max.overlaps = 100) +
  #coord_cartesian(xlim=c(-6,8), ylim=c(0,155)) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))





Day65_PWY_plot <- EnhancedVolcano(Day65_PWY,
                                  lab=Day65_PWY$feature,x='FC',y='qval',
                                  colAlpha = 1,
                                  title = NULL,
                                  subtitle = NULL,
                                  pointSize = 3.0,
                                  labSize = 3,
                                  drawConnectors = TRUE,
                                  typeConnectors='closed',
                                  endsConnectors='last',
                                  legendPosition = 'blank',
                                  pCutoff = 0.000001,
                                  FCcutoff = 2.0,
                                  axisLabSize=4,
                                  caption=NULL,
                                  #shapeCustom = F_keyvals.shape2,
                                  max.overlaps = 100) +
  #coord_cartesian(xlim=c(-6,8), ylim=c(0,155)) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))




Day100_PWY_plot <- EnhancedVolcano(Day100_PWY,
                                  lab=Day100_PWY$feature,x='FC',y='qval',
                                  colAlpha = 1,
                                  title = NULL,
                                  subtitle = NULL,
                                  pointSize = 3.0,
                                  labSize = 3,
                                  drawConnectors = TRUE,
                                  typeConnectors='closed',
                                  endsConnectors='last',
                                  legendPosition = 'blank',
                                  pCutoff = 0.000001,
                                  FCcutoff = 2.0,
                                  axisLabSize=4,
                                  caption=NULL,
                                  #shapeCustom = F_keyvals.shape2,
                                  max.overlaps = 100) +
  #coord_cartesian(xlim=c(-6,8), ylim=c(0,155)) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))



library(ggpubr)

PWY_plot_merged <- ggarrange(Day2_PWY_plot,
                         Day5_PWY_plot,
                         Day8_PWY_plot,
                         Day12_PWY_plot,
                         Day65_PWY_plot,
                         Day100_PWY_plot,
                         common.legend = FALSE, legend = "none",
                         ncol = 3, nrow = 2, align="hv")























###Female Day A1

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))

#meta_ob1 <- meta_ob1[meta_ob1$Sex=="Female" & meta_ob1$DayPlot=="A1",]
#meta_ob1 <- meta_ob1[meta_ob1$Group=="IA" & meta_ob1$DayPlot=="C5",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 







#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")



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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)





library(Maaslin2)




fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Sex_IA_PWY_C5", 
                    fixed_effects  = c("Sex"),
                    reference      = c("Sex,Female"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)


















fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Female_NEGBIN_PWY_A1", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)







###Female Day B2

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Female" & meta_ob1$DayPlot=="B2",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 



#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")



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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Female_NEGBIN_PWY_B2", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)







###Female Day C5

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Female" & meta_ob1$DayPlot=="C5",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 



#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")


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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Female_NEGBIN_PWY_C5", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)







###Female Day D8

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Female" & meta_ob1$DayPlot=="D8",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 



#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")


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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Female_NEGBIN_PWY_D8", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)







###Female Day E12

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Female" & meta_ob1$DayPlot=="E12",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 



#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")


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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Female_NEGBIN_PWY_E12", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)






###Female Day G65

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Female" & meta_ob1$DayPlot=="G65",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 



#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")


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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Female_NEGBIN_PWY_G65", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)






###Female Day H100

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Female" & meta_ob1$DayPlot=="H100",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 



#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")


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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Female_NEGBIN_PWY_H100", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)













###Male Day A1

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Male" & meta_ob1$DayPlot=="A1",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 



#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")


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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Male_NEGBIN_PWY_A1", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)







###Male Day B2

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Male" & meta_ob1$DayPlot=="B2",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 



#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")


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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)



library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Male_NEGBIN_PWY_B2", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)







###Male Day C5

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Male" & meta_ob1$DayPlot=="C5",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 



#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")


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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Male_NEGBIN_PWY_C5", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)







###Male Day D8

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Male" & meta_ob1$DayPlot=="D8",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 



#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")


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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Male_NEGBIN_PWY_D8", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)







###Male Day E12

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Male" & meta_ob1$DayPlot=="E12",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 



#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")


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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)



library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Male_NEGBIN_PWY_E12", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)






###Male Day G65

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Male" & meta_ob1$DayPlot=="G65",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 



#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")


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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Male_NEGBIN_PWY_G65", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)






###Male Day H100

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/PWY_taxa.txt"))
meta_ob1 <- meta_ob1[meta_ob1$Sex=="Male" & meta_ob1$DayPlot=="H100",]

nrow(meta_ob1)

GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.csv",)

fullMothur<-read.csv("PWY_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f 



#pruned_ps <- prune_taxa(names(sort(taxa_sums(merged_phy_percent),TRUE)[1:200]), merged_phy_percent)
merged_phy_percent1 <- tax_glom(merged_phy_percent1 , taxrank="Rank1")


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
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)






library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Male_NEGBIN_PWY_H100", 
                    fixed_effects  = c("Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)














###Volcano plots




library(EnhancedVolcano)


data_dif <- read.table("D://ABX_LOCO_2024/PWY_FC_Female.txt", header=T, row.names=1, check.names = FALSE)


#Female_data_dif <- data_dif[data_dif $Sex=="Female",]
#Female_data_dif2 <- Female_data_dif[,c("feature","LF", "qval")]




EnhancedVolcano(data_dif,
                lab = rownames(data_dif),
                x = 'LF', #log2(exp(coef))
                y = 'qval')


EnhancedVolcano(data_dif,
                lab = rownames(data_dif),
                x = 'LF',
                y = 'qval',
                #pCutoff = 0.1,
                title = NULL,
                subtitle = NULL,
                #FCcutoff = 25,
                pointSize = 3.0,
                labSize = 3,
                col=c('gray', 'blue', 'green', 'red3'),
                colAlpha = 1,
                drawConnectors = TRUE,
                typeConnectors='closed',
                endsConnectors='last',
                #legendPosition = 'blank',
                max.overlaps=15,
                #legendLabSize = 10,
                #legendIconSize = 3.0,
                #legendLabels=c('Not sig.','Log (base 2) FC','', 'qvalue & Log (base 2) FC'))
                pCutoff = 0.1,
                FCcutoff = 0.0) + 
  ggplot2::coord_cartesian(xlim=c(-2, 2)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-8,10, 1)) +
  ggplot2::coord_cartesian(ylim=c(0, 8)) +
  ggplot2::scale_y_continuous(
    breaks=seq(0,175, 25)) +
  theme(axis.text=element_text(colour="black"))


EnhancedVolcano(data_dif,
                lab = rownames(data_dif),
                x = 'LF',
                y = 'qval',
                title = 'N061011 versus N61311',
                #pCutoff = 10e-16,
                #FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 6.0,
                col=c('black', 'blue', 'green', 'red3'),
                colAlpha = 1)



p1 <- EnhancedVolcano(data_dif,
                      lab = rownames(data_dif),
                      x = "LF",
                      y = "qval",
                      #pCutoff = 0.1,
                      max.overlaps = 40,
                      FCcutoff = 0,
                      #ylim = c(0, -log10(10e-12)),
                      #pointSize = c(ifelse(data_dif$coef>2, 8, 1)),
                      pointSize = 3.5,
                      labSize = 6.0,
                      #shape = c(6, 6, 19, 16),
                      title = "Saline vs Cocaine",
                      subtitle = "Differential abundance",
                      #caption = bquote(~Log[2]~ "fold change cutoff, 2; p-value cutoff, 10e-4"),
                      legendPosition = "right",
                      legendLabSize = 14,
                      col = c("grey30", "forestgreen", "royalblue", "red2"),
                      colAlpha = 0.9,
                      #drawConnectors = TRUE,
                      #hline = c(10e-8),
                      widthConnectors = 0.5)

p1


















































































Female_df <- read.table("D://ABX_LOCO_2024/Female_plot_data_selected_250.txt", header=T, row.names=1)

#df <- meta_ob1[meta_ob1$Sex=="Female",]



Female_abundance_plot <- ggplot(data = Female_df,
       aes(x=as.factor(Rank6), y=Abundance, fill = factor(Group))) +
  geom_boxplot(fatten=NULL, outlier.shape = NA, position = position_dodge(width=0.85), alpha=1) +
  #geom_point(position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.85), aes(fill = Group), pch = 21, size=1) +
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), position = position_dodge(width = 0.83), width = 0.74, linewidth=1) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  facet_wrap(~DayPlot, scales = "free_y", ncol = 1, nrow = 7) +
  #ylab("16S rRNA gene copies per mg of feces relative to same day control sample (%)") +
  theme_classic() +
  theme(axis.title.y = element_text(size=8)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), breaks = scales::pretty_breaks(n = 4)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=8),
        axis.text.x = element_text(color="black", size=8),
        axis.text.y = element_text(color="black", size=8),
        legend.title=element_blank()) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )








Male_df <- read.table("D://ABX_LOCO_2024/Male_plot_data_selected_250.txt", header=T, row.names=1)

#df <- meta_ob1[meta_ob1$Sex=="Female",]



Male_abundance_plot <- ggplot(data = Male_df,
                                aes(x=as.factor(Rank6), y=Abundance, fill = factor(Group))) +
  geom_boxplot(fatten=NULL, outlier.shape = NA, position = position_dodge(width=0.85), alpha=1) +
  #geom_point(position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.85), aes(fill = vs), pch = 21, size=1) +
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), position = position_dodge(width = 0.83), width = 0.74, linewidth=1) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  facet_wrap(~DayPlot, scales = "free_y", ncol = 1, nrow = 7) +
  #ylab("16S rRNA gene copies per mg of feces relative to same day control sample (%)") +
  theme_classic() +
  theme(axis.title.y = element_text(size=8)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), breaks = scales::pretty_breaks(n = 4)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=8),
        axis.text.x = element_text(color="black", size=8),
        axis.text.y = element_text(color="black", size=8),
        legend.title=element_blank()) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )







library(ggpubr)

abundance_plot_merged <- ggarrange(Female_abundance_plot,
                         Male_abundance_plot,
                         common.legend = TRUE, legend = "right",
                         ncol = 2, nrow = 1, align="hv")




dev.off()







#input_metadata$GADay = (input_metadata$Group == "GA") *
  #input_metadata$Day
#input_metadata$IADay = (input_metadata$Group == "IA") *
  #input_metadata$Day



#colnames(input_data) == row.names(input_metadata)


library(Maaslin2)




fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "Female_NEGBIN_Day1", 
                    fixed_effects  = c("Group", "CPmg"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)







fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "Female_NEGBIN", 
                    fixed_effects  = c("DayPlot", "Group", "CPmg"),
                    reference      = c("DayPlot,A1", "Group,C"),
                    random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)



fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    #min_prevalence=0.033,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "Male_PWY_with_CPmg_out_NEGBIN", 
                    fixed_effects  = c("Group", "DayPlot", "CPmg"),
                    reference      = c("Group,C", "DayPlot,A1"),
                    random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)



fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    #min_prevalence=0.033,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "Male_with_CPmg_out_NEGBIN", 
                    fixed_effects  = c("Group", "DayPlot", "CPmg"),
                    reference      = c("Group,C", "DayPlot,A1"),
                    random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)


















###GAM


df <- as.data.frame(cbind(meta_phy_out1$Day, meta_phy_out1$Group, meta_phy_out1$Sex, meta_phy_out1$CPmg))
colnames(df)[1] <- "Day"
colnames(df)[2] <- "Group"
colnames(df)[3] <- "Sex"
colnames(df)[4] <- "CPmg"
df$CPmg <- as.numeric(df$CPmg)
df$Group <- as.factor(df$Group)
df$Sex <- as.factor(df$Sex)






df <- read.table("D://ABX_LOCO_2024/PWY_for_GAM.txt", row.names=1, header=T)
df$CPmg <- as.numeric(df$CPmg)
df$Group <- as.factor(df$Group)
df$Sex <- as.factor(df$Sex)



df <- subset(df, Day=="8")



ggplot(df, aes(Lactobacillus, Peptidoglycan_Recycling_I, colour = df$Group, fill=df$Group, shape=df$Group)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  scale_y_continuous() +
  facet_wrap(~Sex, scales = "free_x", ncol = 2, nrow = 1)







###PWY vs. genus maaslin

setwd("D://ABX_LOCO_2024/")
meta_ob1 <- read.table("D://ABX_LOCO_2024/PWY_for_GAM.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("D://ABX_LOCO_2024/taxa.txt"))




meta_ob1 <- meta_ob1[meta_ob1$Day=="8" & meta_ob1$Sex=="Male",]

nrow(meta_ob1)




#Get sample names
GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("CEF_16S_shared.csv",)
#fullMothur<-read.csv("maaslin_selected_shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
library("ape")
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))

#random_tree = rtree(ntaxa(phyloseqobj.f ), rooted=TRUE, tip.label=taxa_names(phyloseqobj.f ))

#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1), random_tree)
#phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

phyloseqobj.f <- tax_glom(phyloseqobj.f, taxrank="Rank6")

merged_phy_percent1 <- phyloseqobj.f


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
tax_table(merged_phy_percent1)[,6]
taxa_names(merged_phy_percent1) <- tax_table(merged_phy_percent1)[,6]


meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

nrow(meta_phy_out1)
ncol(shared_phy_out1)

input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)





library(Maaslin2)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.01,
                    #min_abundance=0.0,
                    #min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "TMM",
                    output         = "Peptidoglycan_Recycling_I_Male_Day_8", 
                    fixed_effects  = c("Peptidoglycan_Recycling_I", "Group"),
                    reference      = c("Group,C"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)



#Butanediol_Biosynthesis
#Peptidoglycan_Recycling_I

















library(mgcv)

mod1 <- gam(Lactobacillus ~ s(Day, bs = "cc", k = 7), method = "REML", data=df)
plot(mod1)


mod2 <- gam(Lactobacillus ~ s(Day, by = df$Group, k = 7), method = "REML", data=df)
plot(mod2)


mod3 <- gam(Lactobacillus ~ s(Day, by = df$Group, k = 7), data=df)
plot(mod3)



dev.off()

anova(mod3)



ggplot(df, aes(Day, Lactobacillus, colour = df$Group, fill=df$Group, shape=df$Group)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 7)) +
  scale_y_continuous()



ggplot(df, aes(Day, Peptidoglycan_Recycling_I, colour = df$Group, fill=df$Group, shape=df$Group)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 7)) +
  scale_y_continuous() +
  facet_wrap(~Sex, scales = "free_x", ncol = 2, nrow = 1)
  






library(ggplot2)

ggplot(aes(Day, Lactobacillus)) +
  geom_point() +
  geom_abline(intercept = coef(lm.fit)[1], slope = coef(lm.fit)[2])





#"Bacilli" =ASV1
#"Bacteroidia" = ASV5
#"Clostridia" = ASV9
#"Gammaproteobacteria" = ASV14




#write.table(df, file='/Users/15869/Desktop/Kuhn_Lab/CEF_16S/df.txt', sep='\t')
#df <- read.table("/Users/15869/Desktop/Kuhn_Lab/CEF_16S/df.txt", header=T)


library(mgcv)


ggplot(df, aes(x, Bacilli, colour = Group, fill=Group, shape=Group)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3))



  



dev.off()

alpha_plot <- plot_richness(merged_phy_percent1, x ="DayPlot", color="Group" , measures=c("Chao1", "Shannon", "Simpson")) +
facet_wrap(Sex~variable, scales = "free_y")

#alpha_plot_data <- ggplot_build(alpha_plot)$plot$data
#write.csv(alpha_plot_data, file="alpha_plot_data.csv")




#Alpha diversity





library(lmerTest)
library(car)
library(caret)
library(lmtest)




library(lme4)
model1 <- lmer(invsimpson ~ Group*as.factor(meta_ob1$Day)*Sex + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=FALSE)

model1 <- lm(invsimpson ~ Group*as.factor(meta_ob1$Day)*Sex, data = meta_ob1)


library(emmeans)
options(max.print=1000000)
options(digits = 5)
emm_options(opt.digits = FALSE)
emmeans(model1, list(pairwise ~ Group*as.factor(Day)*Sex), adjust = "tukey")










#bray PCoA!
dist_bray <- capscale(shared_phy_out1~1, distance="bray", binary=TRUE, eig=TRUE)

prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)

data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)

df <- data.frame(data.scores)

library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))

library(dplyr)
library(ggrepel)
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = shared_phy_out1)

#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)


bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = meta_phy_out1$DayPlot, fill=meta_phy_out1$DayPlot, shape=meta_phy_out1$DayPlot)) +
  scale_shape_manual(values = c(15, 16, 17, 18, 19, 25)) +
  geom_point(mapping = aes(colour = meta_phy_out1$DayPlot), size=5, alpha=1.0) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  #facet_wrap(~meta_ob1$Sex, scales = "free_x", ncol = 2) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=10, face=NULL), axis.text.y = element_text(color="black"))
  







#Bray PCoA
dist_bray <- capscale(shared_phy_out1~1, distance="bray", binary=TRUE, eig=TRUE)

prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)

data.scores = scores(dist_bray, display = "sites", choices = c(1:3))
#summary(data.scores)

df <- data.frame(data.scores)

library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))

library(dplyr)
library(ggrepel)
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = shared_phy_out1)

#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)


bray_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = meta_phy_out1$Cluster5, fill=meta_phy_out1$Cluster5, shape=meta_phy_out1$Cluster5)) +
  scale_color_manual(values=c("#084287", "#5A9D5A", "#F6D832", "#E41A1C", "#CF9C76", "#7fff00", "#084287", "#3F94C3", "#8DD3C7", "#E7298A", "#FFAD12", "#D0CD66", "#45A939", "#91569A", "#D58EC4", "#B6742A", "#999999")) +
  scale_fill_manual(values=c("#084287", "#5A9D5A", "#F6D832", "#E41A1C", "#CF9C76", "#7fff00", "#084287", "#3F94C3", "#8DD3C7", "#E7298A", "#FFAD12", "#D0CD66", "#45A939", "#91569A", "#D58EC4", "#B6742A", "#999999")) +
    scale_shape_manual(values = c(15, 16, 17, 18, 19, 25)) +
  geom_point(mapping = aes(colour = meta_phy_out1$Cluster5), size=5, alpha=1.0) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  theme_classic() +
  labs(x=labx1, y=laby1) +
  #facet_wrap(~meta_ob1$Sex, scales = "free_x", ncol = 2) +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=10, face=NULL), axis.text.y = element_text(color="black"))



library(ggpubr)

beta_merged <- ggarrange(bray_plot, bray_plot,
                          common.legend = TRUE, legend = "right",
                          ncol = 1, nrow = 2, align="hv")




library(pairwiseAdonis)


dist.bray <- vegdist(shared_phy_out1, method="bray", binary=FALSE)
adonis2(formula = dist.bray  ~ DayPlot, data = meta_phy_out1, permutations=1000)



adonis2(formula = x ~ Mouse_ID, data = meta_phy_out1, permutations = 100000)
adonis2(formula = x ~ DayPlot, data = meta_phy_out1, permutations = 100000)

x <- phyloseq::distance(merged_phy_percent1, method ="bray") 
factors <- meta_phy_out1$DayPlot
pairwise.adonis(x, factors, perm = 10000, p.adjust(method="FDR"))
#pairwise.adonis(x, factors, p.adjust(method="holm"))





merged_phy_percent1 <- phyloseqobj.f

meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))

















#Calculate bray distance matrix on dataset using vegan
dist.jacc <- vegdist(shared_phy_out1, method="bray", binary=TRUE)
#Generate plottable bray distance matrix object
veg_jacc <- capscale(dist.jacc~1, distance="bray", binary=TRUE)
summary(eigenvals(veg_jacc))
#Calculate Bray-Curtis distance matrix on dataset using vegan
dist.bray <- vegdist(shared_phy_out1, method="bray", binary=FALSE)
#Generate plottable Bray-Curtis distance matrix object
veg_bray <- capscale(dist.bray ~1, distance="bray", binary=FALSE)
summary(eigenvals(veg_bray))


cols = c("navyblue", "red3", "green", "magenta", "yellow","blue","black","pink","cyan","gray","orange","brown","purple")



library(dplyr)



# Convert metadata to numeric for plotting
data_meta = meta_ob1 %>%
  mutate(Sex = case_when(Sex == "Female" ~ 1,
                         Sex == "Male" ~ 2))

data_meta = meta_ob1 %>%
  mutate(Group = case_when(Group == "C" ~ 1,
                           Group == "GA" ~ 2,
                           Group == "IA" ~ 3))


data_meta = meta_ob1 %>%
  mutate(DayPlot = case_when(DayPlot == "A1" ~ 1,
                         DayPlot == "B2" ~ 2,
                         DayPlot == "C5" ~ 3,
                         DayPlot == "D8" ~ 4,
                         DayPlot == "E12" ~ 5,
                         DayPlot == "F100" ~ 6))
                         
                         
                    


plot(veg_jacc$Ybar[,c(1,2)], xlab = "PC1 (XX.X%)", ylab = "PC2 (XX.X%)", col = cols[data_meta$Group], pch = c(21,21)[data_meta$Group], lwd=2, cex=2, xaxs='r', yaxs='r', main="bray")
with(data_meta, ordiellipse(veg_bray$Ybar[,c(1,2)], Group, col= cols, lty=1, lwd=2, label = FALSE, cex=2.0, alpha = .5))
#with(data_meta , ordihull(veg_bray$Ybar[,c(1,2)], Sex, col= "black", lty=1, lwd=0.1, label = FALSE))
legend("topright", legend = c("C", "GA", "IA"), bty = "n", col = cols, pch = 21, cex=2.0)


plot(veg_bray$Ybar[,c(1,2)], xlab = "PC1 (XX.X%)", ylab = "PC2 (XX.X%)", col = cols[data_meta$Group], pch = c(21,21)[data_meta$Group], lwd=2, cex=2, xaxs='r', yaxs='r', main="bray")
with(data_meta, ordiellipse(veg_bray$Ybar[,c(1,2)], Group, col= cols, lty=1, lwd=2, label = FALSE, cex=2.0, alpha = .5))
#with(data_meta , ordihull(veg_bray$Ybar[,c(1,2)], Sex, col= "black", lty=1, lwd=0.1, label = FALSE))
legend("topright", legend = c("C", "GA", "IA"), bty = "n", col = cols, pch = 21, cex=2.0)








adonis2(formula = dist.jacc  ~ Mouse_ID, data = meta_phy_out1, permutations=1000)
adonis2(formula = dist.bray  ~ DayPlot, data = meta_phy_out1, permutations=1000)







pdf(file="bray_PCoA.pdf", width = 6.4, height = 5.5)
veg_bray <- capscale(shared_phy_out1~1, distance="bray", binary=FALSE)
summary(eigenvals(veg_bray))
ef <- envfit(veg_bray, meta_phy_out1, permu = 100, na.rm=TRUE)
ef
plot(veg_bray, type="n", xlab = "PC1 (XX.X%)", ylab = "PC2 (XX.X%)")

cols=c("dodgerblue", "chartreuse4", "firebrick3", "palevioletred1", "gold", "slateblue1", "black", "yellow4", "chartreuse", "darkgreen", "gray67", "coral3", "darkolivegreen3", "darkmagenta", "turquoise2", "red", "mediumvioletred", "magenta", "blue")

#plot(ef, col = "black", p.max = 0.05)
points(veg_bray, col = cols[meta_phy_out1$Group], pch = 19, cex=2.0)
legend("topright", legend = unique(meta_phy_out1$Group),
       bty = "n", col = cols, pch = 19, cex=1.5)
scrs <- scores(veg_bray, display = c("sites", "species"),)
with(scrs, text(species, labels = rownames(species),
                col = "black", cex = 1.0))
#with(scrs, text(sites, labels = rownames(sites),
#col = "red", cex = 0.8))
with(meta_phy_out1, ordiellipse(veg_bray, Group, col=cols, lty=5, lwd=2, kind = "se", conf = 0.95, label = TRUE, cex=1.0))
dev.off()



pdf(file="bray_PCoA.pdf", width = 6.4, height = 5.5)
veg_bray <- capscale(shared_phy_out1~1, distance="bray", binary=FALSE)
summary(eigenvals(veg_bray))
ef <- envfit(veg_bray, meta_phy_out1, permu = 100, na.rm=TRUE)
ef
plot(veg_bray, type="n", xlab = "PC1 (XX.X%)", ylab = "PC2 (XX.X%)")

cols=c("dodgerblue", "chartreuse4", "firebrick3", "palevioletred1", "gold", "slateblue1", "black", "yellow4", "chartreuse", "darkgreen", "gray67", "coral3", "darkolivegreen3", "darkmagenta", "turquoise2", "red", "mediumvioletred", "magenta", "blue")

#plot(ef, col = "black", p.max = 0.05)
points(veg_bray, col = cols[meta_phy_out1$Group], pch = 19, cex=2.0)
legend("topright", legend = unique(meta_phy_out1$Group),
       bty = "n", col = cols, pch = 19, cex=1.5)
scrs <- scores(veg_bray, display = c("sites", "species"),)
with(scrs, text(species, labels = rownames(species),
                col = "black", cex = 1.0))
#with(scrs, text(sites, labels = rownames(sites),
#col = "red", cex = 0.8))
with(meta_phy_out1, ordiellipse(veg_bray, Group, col=cols, lty=5, lwd=2, kind = "se", conf = 0.95, label = TRUE, cex=1.0))
dev.off()







dev.off()


merged_phy_percent2 <- tax_glom(merged_phy_percent1, taxrank="Rank3")

#merged_phy_percent2 <- phyloseq_standardize_otu_abundance(merged_phy_percent , method = "total")


SI_pruned <- prune_taxa(names(sort(taxa_sums(merged_phy_percent2),TRUE)[1:5]), merged_phy_percent2)
sampleOrder = unique(sample_names(SI_pruned))
taxaOrder = rev(unique(taxa_names(SI_pruned)))


bar_plot <- plot_bar(merged_phy_percent1  , fill="Rank1") + facet_wrap(~Sex, scales = "free_x") +
  theme(legend.position = "bottom")
#scale_x_discrete(labels = bar_plot$data$Mouse_ID)










input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)




#colnames(input_data) == row.names(input_metadata)


library(Maaslin2)

fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    #min_prevalence=0.033,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "out_NEGBIN", 
                    fixed_effects  = c("Group", "Sex", "DayPlot"),
                    reference      = c("Group,C", "Sex,Female", "DayPlot,A1"),
                    random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)



fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    #min_prevalence=0.033,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "out_LC_NEGBIN", 
                    fixed_effects  = c("Group", "Sex", "DayPlot", "LC"),
                    reference      = c("Group,C", "Sex,Female", "DayPlot,A1"),
                    random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)






library(dendextend)
library(circlize)
library(magrittr)

distance = distance(merged_phy_percent1, method="bray", binary=FALSE)
#distance = dist(phyloseqobj.f, method ="bray")    

hcluster = hclust(distance, method ="ward.D2")
dend <- as.dendrogram(hcluster)

#plot(dend)


cols <- c("#3F94C1", "#084287", "#5A9D5A", "#F6D832", "#E41A1C", "#CF9C76", "#7fff00", "#084287", "#8DD3C7", "#E7298A", "#FFAD12", "#D0CD66", "#45A939", "#91569A", "#D58EC4", "#B6742A","#999999")

dend <- color_branches(dend, k = 3, col = cols)
dend %<>% set("labels_col", value = cols, k= 3)
dend %<>% set("labels_cex", 1)
dend %<>% set("branches_lwd", 1)

plot(dend)

#circlize_dendrogram(dend)

dev.off()



library(RFLPtools)
write.hclust(hcluster, "clusters3.txt", "hello", k = 3, append = FALSE, dec = ",")















setwd("D://ABX_LOCO_2024/")

meta_ob1 <- read.table("D://ABX_LOCO_2024/meta.txt", header=T, row.names=1)



#meta_ob1 <- meta_ob1[meta_ob1$Cohort=="C2",]


dev.off()

library(ggalluvial)


meta_ob1_Female <- meta_ob1[meta_ob1$Sex=="Female",]

meta_ob1_Male <- meta_ob1[meta_ob1$Sex=="Male",]



Female_plot <- ggplot(meta_ob1_Female, aes(x = `DayPlot`, stratum = GenusCluster4, alluvium = Mouse_ID, fill = GenusCluster4, label = GenusCluster4)) +
  scale_fill_manual(values=c("#084287", "#5A9D5A", "#CF9C76", "#E41A1C",  "#8DD3C7", "#F6D832", "#CF9C76", "#7fff00", "#084287", "#3F94C3", "#8DD3C7", "#E7298A", "#FFAD12", "#D0CD66", "#45A939", "#91569A", "#D58EC4", "#B6742A", "#999999")) +
  #scale_fill_brewer(type = "qual", palette = "Spectral") +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray") +
  geom_stratum() + theme(legend.position = "right") +
  facet_wrap(~meta_ob1_Female$Group, ncol = 3)
  


Male_plot <- ggplot(meta_ob1_Male, aes(x = `DayPlot`, stratum = GenusCluster4, alluvium = Mouse_ID, fill = GenusCluster4, label = GenusCluster4)) +
  scale_fill_manual(values=c("#084287", "#5A9D5A", "#CF9C76", "#E41A1C",  "#8DD3C7", "#F6D832", "#CF9C76", "#7fff00", "#084287", "#3F94C3", "#8DD3C7", "#E7298A", "#FFAD12", "#D0CD66", "#45A939", "#91569A", "#D58EC4", "#B6742A", "#999999")) +
  #scale_fill_brewer(type = "qual", palette = "Spectral") +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray") +
  geom_stratum() + theme(legend.position = "right") +
  facet_wrap(~meta_ob1_Female$Group, ncol = 3)












Female_plot <- ggplot(meta_ob1_Female,
       aes(x = DayPlot, stratum = GenusCluster3, alluvium = Mouse_ID,
           fill = GenusCluster3, label = GenusCluster3)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_stratum() +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  theme(legend.position = "bottom") +
  ggtitle("Female") +
  facet_wrap(~meta_ob1_Female$Group, ncol = 3)


Male_plot <- ggplot(meta_ob1_Male,
       aes(x = DayPlot, stratum = GenusCluster3, alluvium = Mouse_ID,
           fill = GenusCluster3, label = GenusCluster3)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_stratum() +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  theme(legend.position = "bottom") +
  ggtitle("Male") +
  facet_wrap(~meta_ob1_Female$Group, ncol = 3)






library(ggpubr)

merged_plot <- ggarrange(Female_plot, Male_plot,
                         common.legend = FALSE, legend = "right",
                         ncol = 1, nrow = 2, align="hv")








ggplot(as.data.frame(meta_ob1_Female),
       aes(y = Freq,
           axis1 = DayPlot, axis2 = Mouse_ID)) +
  geom_flow(aes(fill = GenusCluster4), width = .4, curve_type = "quintic") +
  geom_stratum(width = .4) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Class", "Sex")) +
  facet_wrap(~ Sex, scales = "fixed")












###PWY plots

setwd("D://ABX_LOCO_2024/")

data_df <- read.table("D://ABX_LOCO_2024/PWY_FC.txt", header=T, row.names=1, check.names = FALSE)

library(dplyr)

#F_data_df <- subset(data_df, abs(FC) > 2.0 & qval<0.01 & Sex=="Female")
#M_data_df <- subset(data_df, abs(FC) > 2.0 & qval<0.01 & Sex=="Male")


F_data_df <- subset(data_df, Sex=="Female" & N.not.0 > 16)
M_data_df <- subset(data_df, Sex=="Male" & N.not.0 > 16)


###Split days
Female_Day2_PWY <- subset(F_data_df, Day=="B2")
Female_Day5_PWY <- subset(F_data_df, Day=="C5")
Female_Day8_PWY <- subset(F_data_df, Day=="D8")
Female_Day12_PWY <- subset(F_data_df, Day=="E12")
Female_Day65_PWY <- subset(F_data_df, Day=="G65")
Female_Day100_PWY <- subset(F_data_df, Day=="H100")

Male_Day2_PWY <- subset(M_data_df, Day=="B2")
Male_Day5_PWY <- subset(M_data_df, Day=="C5")
Male_Day8_PWY <- subset(M_data_df, Day=="D8")
Male_Day12_PWY <- subset(M_data_df, Day=="E12")
Male_Day65_PWY <- subset(M_data_df, Day=="G65")
Male_Day100_PWY <- subset(M_data_df, Day=="H100")


library(EnhancedVolcano)





F_keyvals.shape2 <- Female_Day2_PWY$shape
names(F_keyvals.shape2)[F_keyvals.shape2 == 1] <- 'GA'
names(F_keyvals.shape2)[F_keyvals.shape2 == 2] <- 'IA'

F_keyvals.shape5 <- Female_Day5_PWY$shape
names(F_keyvals.shape5)[F_keyvals.shape5 == 1] <- 'GA'
names(F_keyvals.shape5)[F_keyvals.shape5 == 2] <- 'IA'

F_keyvals.shape8 <- Female_Day8_PWY$shape
names(F_keyvals.shape8)[F_keyvals.shape8 == 1] <- 'GA'
names(F_keyvals.shape8)[F_keyvals.shape8 == 2] <- 'IA'

F_keyvals.shape12 <- Female_Day12_PWY$shape
names(F_keyvals.shape12)[F_keyvals.shape12 == 1] <- 'GA'
names(F_keyvals.shape12)[F_keyvals.shape12 == 2] <- 'IA'

F_keyvals.shape65 <- Female_Day65_PWY$shape
names(F_keyvals.shape65)[F_keyvals.shape65 == 1] <- 'GA'
names(F_keyvals.shape65)[F_keyvals.shape65 == 2] <- 'IA'

F_keyvals.shape100 <- Female_Day100_PWY$shape
names(F_keyvals.shape100)[F_keyvals.shape100 == 1] <- 'GA'
names(F_keyvals.shape100)[F_keyvals.shape100 == 2] <- 'IA'




Female_Day2_PWY_plot <- EnhancedVolcano(Female_Day2_PWY,
                                        lab=Female_Day2_PWY$feature,x='FC',y='qval',
                                        colAlpha = 1,
                                        title = NULL,
                                        subtitle = NULL,
                                        pointSize = 3.0,
                                        labSize = 3,
                                        drawConnectors = TRUE,
                                        typeConnectors='closed',
                                        endsConnectors='last',
                                        legendPosition = 'blank',
                                        pCutoff = 0.0000001,
                                        FCcutoff = 2.0,
                                        axisLabSize=4,
                                        caption=NULL,
                                        shapeCustom = F_keyvals.shape2,
                                        max.overlaps = 100) +
  coord_cartesian(xlim=c(-6,8), ylim=c(0,155)) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))


Female_Day5_PWY_plot <- EnhancedVolcano(Female_Day5_PWY,
                                        lab=Female_Day5_PWY$feature,x='FC',y='qval',
                                        colAlpha = 1,
                                        title = NULL,
                                        subtitle = NULL,
                                        pointSize = 3.0,
                                        labSize = 3,
                                        drawConnectors = TRUE,
                                        typeConnectors='closed',
                                        endsConnectors='last',
                                        legendPosition = 'blank',
                                        pCutoff = 0.0000001,
                                        FCcutoff = 2.0,
                                        axisLabSize=4,
                                        caption=NULL,
                                        shapeCustom = F_keyvals.shape5,
                                        max.overlaps = 100) +
  coord_cartesian(xlim=c(-6,8), ylim=c(0,155)) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))



Female_Day8_PWY_plot <- EnhancedVolcano(Female_Day8_PWY,
                                        lab=Female_Day8_PWY$feature,x='FC',y='qval',
                                        colAlpha = 1,
                                        title = NULL,
                                        subtitle = NULL,
                                        pointSize = 3.0,
                                        labSize = 3,
                                        drawConnectors = TRUE,
                                        typeConnectors='closed',
                                        endsConnectors='last',
                                        legendPosition = 'blank',
                                        pCutoff = 0.0000001,
                                        FCcutoff = 2.0,
                                        axisLabSize=4,
                                        caption=NULL,
                                        shapeCustom = F_keyvals.shape8,
                                        max.overlaps = 100) +
  coord_cartesian(xlim=c(-6,8), ylim=c(0,155)) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

Female_Day12_PWY_plot <- EnhancedVolcano(Female_Day12_PWY,
                                         lab=Female_Day12_PWY$feature,x='FC',y='qval',
                                         colAlpha = 1,
                                         title = NULL,
                                         subtitle = NULL,
                                         pointSize = 3.0,
                                         labSize = 3,
                                         drawConnectors = TRUE,
                                         typeConnectors='closed',
                                         endsConnectors='last',
                                         legendPosition = 'blank',
                                         pCutoff = 0.0000001,
                                         FCcutoff = 2.0,
                                         axisLabSize=4,
                                         caption=NULL,
                                         shapeCustom = F_keyvals.shape12,
                                         max.overlaps = 100) +
  coord_cartesian(xlim=c(-6,8), ylim=c(0,155)) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))


Female_Day65_PWY_plot <- EnhancedVolcano(Female_Day65_PWY,
                                         lab=Female_Day65_PWY$feature,x='FC',y='qval',
                                         colAlpha = 1,
                                         title = NULL,
                                         subtitle = NULL,
                                         pointSize = 3.0,
                                         labSize = 3,
                                         drawConnectors = TRUE,
                                         typeConnectors='closed',
                                         endsConnectors='last',
                                         legendPosition = 'blank',
                                         pCutoff = 0.0000001,
                                         FCcutoff = 2.0,
                                         axisLabSize=4,
                                         caption=NULL,
                                         shapeCustom = F_keyvals.shape65,
                                         max.overlaps = 100) +
  coord_cartesian(xlim=c(-6,8), ylim=c(0,155)) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))


Female_Day100_PWY_plot <- EnhancedVolcano(Female_Day100_PWY,
                                          lab=Female_Day100_PWY$feature,x='FC',y='qval',
                                          colAlpha = 1,
                                          title = NULL,
                                          subtitle = NULL,
                                          pointSize = 3.0,
                                          labSize = 3,
                                          drawConnectors = TRUE,
                                          typeConnectors='closed',
                                          endsConnectors='last',
                                          legendPosition = 'blank',
                                          pCutoff = 0.0000001,
                                          FCcutoff = 2.0,
                                          axisLabSize=4,
                                          caption=NULL,
                                          shapeCustom = F_keyvals.shape100,
                                          max.overlaps = 100) +
  coord_cartesian(xlim=c(-6,8), ylim=c(0,155)) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))







library(ggpubr)
Female_PWY_plot_merged <- ggarrange(Female_Day2_PWY_plot,
                             Female_Day5_PWY_plot,
                             Female_Day8_PWY_plot,
                             Female_Day12_PWY_plot,
                             Female_Day65_PWY_plot,
                             Female_Day100_PWY_plot,
                            common.legend = TRUE, legend = "none",
                            ncol = 3, nrow = 2, align="hv")




ggsave("Female_PWY_plot_merged.jpeg", Female_PWY_plot_merged, width=10.5, height=8, units="in", limitsize = FALSE)


ggsave("Female_PWY_plot_merged.pdf", Female_PWY_plot_merged, width=10.5, height=8, units="in", limitsize = FALSE)













M_keyvals.shape2 <- Male_Day2_PWY$shape
names(M_keyvals.shape2)[M_keyvals.shape2 == 1] <- 'GA'
names(M_keyvals.shape2)[M_keyvals.shape2 == 2] <- 'IA'

M_keyvals.shape5 <- Male_Day5_PWY$shape
names(M_keyvals.shape5)[M_keyvals.shape5 == 1] <- 'GA'
names(M_keyvals.shape5)[M_keyvals.shape5 == 2] <- 'IA'

M_keyvals.shape8 <- Male_Day8_PWY$shape
names(M_keyvals.shape8)[M_keyvals.shape8 == 1] <- 'GA'
names(M_keyvals.shape8)[M_keyvals.shape8 == 2] <- 'IA'

M_keyvals.shape12 <- Male_Day12_PWY$shape
names(M_keyvals.shape12)[M_keyvals.shape12 == 1] <- 'GA'
names(M_keyvals.shape12)[M_keyvals.shape12 == 2] <- 'IA'

M_keyvals.shape65 <- Male_Day65_PWY$shape
names(M_keyvals.shape65)[M_keyvals.shape65 == 1] <- 'GA'
names(M_keyvals.shape65)[M_keyvals.shape65 == 2] <- 'IA'

M_keyvals.shape100 <- Male_Day100_PWY$shape
names(M_keyvals.shape100)[M_keyvals.shape100 == 1] <- 'GA'
names(M_keyvals.shape100)[M_keyvals.shape100 == 2] <- 'IA'




Male_Day2_PWY_plot <- EnhancedVolcano(Male_Day2_PWY,
                                      lab=Male_Day2_PWY$feature,x='FC',y='qval',
                                      colAlpha = 1,
                                      title = NULL,
                                      subtitle = NULL,
                                      pointSize = 3.0,
                                      labSize = 3,
                                      drawConnectors = TRUE,
                                      typeConnectors='closed',
                                      endsConnectors='last',
                                      legendPosition = 'blank',
                                      pCutoff = 0.0000001,
                                      FCcutoff = 2.0,
                                      axisLabSize=4,
                                      caption=NULL,
                                      shapeCustom = M_keyvals.shape2,
                                      max.overlaps = 100) +
  coord_cartesian(xlim=c(-8,10), ylim=c(0,125)) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))


Male_Day5_PWY_plot <- EnhancedVolcano(Male_Day5_PWY,
                                      lab=Male_Day5_PWY$feature,x='FC',y='qval',
                                      colAlpha = 1,
                                      title = NULL,
                                      subtitle = NULL,
                                      pointSize = 3.0,
                                      labSize = 3,
                                      drawConnectors = TRUE,
                                      typeConnectors='closed',
                                      endsConnectors='last',
                                      legendPosition = 'blank',
                                      pCutoff = 0.0000001,
                                      FCcutoff = 2.0,
                                      axisLabSize=4,
                                      caption=NULL,
                                      shapeCustom = M_keyvals.shape5,
                                      max.overlaps = 100) +
  coord_cartesian(xlim=c(-8,10), ylim=c(0,125)) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))


Male_Day8_PWY_plot <- EnhancedVolcano(Male_Day8_PWY,
                                      lab=Male_Day8_PWY$feature,x='FC',y='qval',
                                      colAlpha = 1,
                                      title = NULL,
                                      subtitle = NULL,
                                      pointSize = 3.0,
                                      labSize = 3,
                                      drawConnectors = TRUE,
                                      typeConnectors='closed',
                                      endsConnectors='last',
                                      legendPosition = 'blank',
                                      pCutoff = 0.0000001,
                                      FCcutoff = 2.0,
                                      axisLabSize=4,
                                      caption=NULL,
                                      shapeCustom = M_keyvals.shape8,
                                      max.overlaps = 100) +
  coord_cartesian(xlim=c(-8,10), ylim=c(0,125)) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))


Male_Day12_PWY_plot <- EnhancedVolcano(Male_Day12_PWY,
                                       lab=Male_Day12_PWY$feature,x='FC',y='qval',
                                       colAlpha = 1,
                                       title = NULL,
                                       subtitle = NULL,
                                       pointSize = 3.0,
                                       labSize = 3,
                                       drawConnectors = TRUE,
                                       typeConnectors='closed',
                                       endsConnectors='last',
                                       legendPosition = 'blank',
                                       pCutoff = 0.0000001,
                                       FCcutoff = 2.0,
                                       axisLabSize=4,
                                       caption=NULL,
                                       shapeCustom = M_keyvals.shape12,
                                       max.overlaps = 100) +
  coord_cartesian(xlim=c(-8,10), ylim=c(0,125)) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))


Male_Day65_PWY_plot <- EnhancedVolcano(Male_Day65_PWY,
                                       lab=Male_Day65_PWY$feature,x='FC',y='qval',
                                       colAlpha = 1,
                                       title = NULL,
                                       subtitle = NULL,
                                       pointSize = 3.0,
                                       labSize = 3,
                                       drawConnectors = TRUE,
                                       typeConnectors='closed',
                                       endsConnectors='last',
                                       legendPosition = 'blank',
                                       pCutoff = 0.0000001,
                                       FCcutoff = 2.0,
                                       axisLabSize=4,
                                       caption=NULL,
                                       shapeCustom = M_keyvals.shape65,
                                       max.overlaps = 100) +
  coord_cartesian(xlim=c(-8,10), ylim=c(0,125)) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))
  



Male_Day100_PWY_plot <- EnhancedVolcano(Male_Day100_PWY,
                                        lab=Male_Day100_PWY$feature,x='FC',y='qval',
                                        colAlpha = 1,
                                        title = NULL,
                                        subtitle = NULL,
                                        pointSize = 3.0,
                                        labSize = 3,
                                        drawConnectors = TRUE,
                                        typeConnectors='closed',
                                        endsConnectors='last',
                                        legendPosition = 'blank',
                                        pCutoff = 0.0000001,
                                        FCcutoff = 2.0,
                                        axisLabSize=4,
                                        caption=NULL,
                                        shapeCustom = M_keyvals.shape100,
                                        max.overlaps = 100) +
  coord_cartesian(xlim=c(-8,10), ylim=c(0,125)) +
  theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))











library(ggpubr)
Male_PWY_plot_merged <- ggarrange(Male_Day2_PWY_plot,
                                  Male_Day5_PWY_plot,
                                  Male_Day8_PWY_plot,
                                  Male_Day12_PWY_plot,
                                  Male_Day65_PWY_plot,
                                  Male_Day100_PWY_plot,
                                  common.legend = TRUE, legend = "none",
                                  ncol = 3, nrow = 2, align="hv")




ggsave("Male_PWY_plot_merged.jpeg", Male_PWY_plot_merged , width=10.5, height=8, units="in", limitsize = FALSE)

ggsave("Male_PWY_plot_merged.pdf", Male_PWY_plot_merged , width=10.5, height=8, units="in", limitsize = FALSE)













ggplot2::coord_cartesian(xlim=c(-2, 2)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-8,10, 1)) +
  ggplot2::coord_cartesian(ylim=c(0, 8)) +
  ggplot2::scale_y_continuous(
    breaks=seq(0,175, 25)) +
  theme(axis.text=element_text(colour="black"))







library(dplyr)
#meta_ob1 %>% count(Sex, Group, Day, Cohort)

means_FC <- data_df  %>%
  group_by(feature, Sex, Day) %>%
  summarise_at(vars(FC), list(name = mean)) %>%
  print(n=100)
write.csv(means_FC, file="means_FC.csv")

means_qval <- data_df  %>%
  group_by(feature, Sex, Day) %>%
  summarise_at(vars(qval), list(name = mean)) %>%
  print(n=100)
write.csv(means_qval, file="means_qval.csv")

means_df <- cbind(means_FC[1:4], means_qval[4])
means_df <- as.data.frame(means_df)
colnames(means_df)[4] <- "FC"
colnames(means_df)[5] <- "qval"


F_data_df <- subset(means_df, abs(FC) > 2.0 & Sex=="Female")
M_data_df <- subset(means_df, abs(FC) > 2.0 & Sex=="Male")







library(ggplot2)






require("ggrepel")
set.seed(1)


F_PWY_df_plot <- ggplot(F_data_df, aes(x = FC, y = -log10(qval), color = Group, shape=Group))+
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(16, 17)) +
  geom_point(size = 2, alpha = 1, na.rm = T) +
  theme_bw() + 
  geom_hline(yintercept = 2.0, colour="gray71", linetype="dashed") + 
  geom_vline(xintercept = 2.0, colour="gray71", linetype="dashed") + 
  geom_vline(xintercept = -2.0, colour="gray71", linetype="dashed") + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.45)) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(color="black", size=11),
        axis.text.y = element_text(color="black", size=11)) +
  scale_y_continuous(trans = "log1p",
                     labels = scales::number_format(accuracy = 0.1),
                     breaks = scales::pretty_breaks(n = 27)) +
  facet_wrap(~Day, scales = "fixed", ncol = 3, nrow = 2) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none") +
  #geom_text_repel(aes(label = feature),
                  #size=4, force_pull=10, box.padding = 0.75, max.overlaps = 100, xlim = c(-5, 5), max.iter=1000000000)
  
geom_text_repel(aes(label = feature),
                size=4, 
box.padding = 0.5, 
point.padding = 1.6,
segment.color = 'gray23',
segment.size = 0.5,
arrow = arrow(length = unit(0.01, 'npc')),
force = 10,
max.iter = 1000000000,
max.overlaps = 100,
min.segment.length = unit(0, 'lines'))



dev.off()

ggsave("F_PWY_plot2.pdf", F_PWY_df_plot, width=3.6, height=5, units="in", scale=5)



M_PWY_df_plot <- ggplot(M_data_df, aes(x = FC, y = -log10(qval), color = Group))+
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_shape_manual(values = c(16, 17)) +
  geom_point(size = 2, alpha = 1, na.rm = T) +
  theme_bw() + 
  geom_hline(yintercept = 2.0, colour="gray71", linetype="dashed") + 
  geom_vline(xintercept = 2.0, colour="gray71", linetype="dashed") + 
  geom_vline(xintercept = -2.0, colour="gray71", linetype="dashed") + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.45)) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(color="black", size=11),
        axis.text.y = element_text(color="black", size=11)) +
  scale_y_continuous(trans = "log1p",
                     labels = scales::number_format(accuracy = 0.1),
                     breaks = scales::pretty_breaks(n = 27)) +
  facet_wrap(~Day, scales = "fixed", ncol = 3, nrow = 2) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none") +
  #geom_text_repel(aes(label = feature),
  #size=4, force_pull=10, box.padding = 0.75, max.overlaps = 100, xlim = c(-5, 5), max.iter=1000000000)
  
  geom_text_repel(aes(label = feature),
                  size=4, 
                  box.padding = 0.5, 
                  point.padding = 1.6,
                  segment.color = 'gray23',
                  segment.size = 0.5,
                  arrow = arrow(length = unit(0.01, 'npc')),
                  force = 10,
                  max.iter = 1000000000,
                  max.overlaps = 100,
                  min.segment.length = unit(0, 'lines'))

#ggsave("M_PWY_plot2.pdf", M_PWY_df_plot, width=3.6, height=5, units="in", scale=5)




data_df <- read.table("D://ABX_LOCO_2024/PWY_FC.txt", header=T, row.names=1, check.names = FALSE)

data_df <- subset(data_df, abs(FC) > 2.0)






limits <- aes(ymax = `FC` + stderr2, ymin = `FC` - stderr2)


F_FC_plot <- ggplot(data = F_data_df,
       aes(x = Order, y = `FC`, width=0.8, group = factor(Group))) +
  geom_bar(stat = "identity",color="black", size=0.5, alpha=0.75,
           aes(fill = factor(Group)),
           position = position_dodge2(preserve = "single")) +
  #ylim(-2.0, 8.0) +
  #ylim(-2.0, 8.0) +
  ylab("Log fold change (natural log)") +
  geom_errorbar(limits, width = 0.2, size=0.5, position = position_dodge(width = 0.9)) +
  coord_flip() +
  theme(axis.title.x=element_blank(),
        text = element_text(size=12), axis.text.y = element_text(color="black")) +
  theme(panel.background = element_rect(fill = "white"),
        legend.key  = element_rect(fill = "white"),
        axis.line.x = element_line(colour = "black", size = 0.5),
        axis.line.y = element_line(colour = "black", size = 0.5)) +
  geom_hline(yintercept=0, linetype="dotted", color = "black", size=1) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y = element_blank()) +
  facet_grid(Sex~Day, space = "free") +
  theme(axis.text.y = element_text(color="black", size=10)) +
  #theme(plot.margin = unit(c(0.5,0.5,0.75,2), "lines")) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        legend.position="right", legend.box = "vertical",
        legend.title= element_blank(),
        axis.text.x = element_text(size=8)) 
  





limits <- aes(ymax = `FC` + stderr2, ymin = `FC` - stderr2)


M_FC_plot <- ggplot(data = M_data_df,
                    aes(x = Order, y = `FC`, width=0.8, group = factor(Group))) +
  geom_bar(stat = "identity",color="black", size=0.5, alpha=0.75,
           aes(fill = factor(Group)),
           position = position_dodge2(preserve = "single")) +
  #ylim(-2.0, 8.0) +
  #ylim(-2.0, 8.0) +
  ylab("Log fold change (natural log)") +
  geom_errorbar(limits, width = 0.2, size=0.5, position = position_dodge(width = 0.9)) +
  coord_flip() +
  theme(axis.title.x=element_blank(),
        text = element_text(size=12), axis.text.y = element_text(color="black")) +
  theme(panel.background = element_rect(fill = "white"),
        legend.key  = element_rect(fill = "white"),
        axis.line.x = element_line(colour = "black", size = 0.5),
        axis.line.y = element_line(colour = "black", size = 0.5)) +
  geom_hline(yintercept=0, linetype="dotted", color = "black", size=1) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y = element_blank()) +
  facet_grid(Sex~Day, space = "free") +
  theme(axis.text.y = element_text(color="black", size=10)) +
  #theme(plot.margin = unit(c(0.5,0.5,0.75,2), "lines")) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        legend.position="right", legend.box = "vertical",
        legend.title= element_blank(),
        axis.text.x = element_text(size=8)) 



library(ggpubr)
FC_merged_plot <- ggarrange(F_FC_plot, M_FC_plot,
                         common.legend = TRUE, legend = "right",
                         ncol = 1, nrow = 2, align="hv")








F_plot <- ggplot(F_data_df, aes(x=feature, y=FC, fill=Group), show_guide=FALSE) + 
  geom_bar(position="dodge",stat="identity") + 
  coord_flip() +
  facet_wrap(~Day, scales = "free_y", ncol = 2, nrow = 3) +
  theme(axis.text.y = element_text(color="black", size=12))

ggsave("F_FC_plot.pdf", F_plot, width=7.2, height=5, units="in", scale=3)


M_plot <- ggplot(M_data_df, aes(x=feature, y=FC, fill=Group), show_guide=FALSE) + 
  geom_bar(position="dodge",stat="identity") + 
  coord_flip() +
  facet_wrap(~Day, scales = "free_y", ncol = 2, nrow = 3) +
  theme(axis.text.y = element_text(color="black", size=12))

ggsave("M_FC_plot.pdf", M_plot, width=7.2, height=5, units="in", scale=3)







#ggsave("FC_plot.pdf", FC_plot, width=7.2, height=5, units="in", scale=2)



dev.off()
























limits <- aes(ymax = `FC` + stderr, ymin = `FC` - stderr)


plot1 <- ggplot(data = F_data_df ,
  aes(x = feature, y = `FC`, width=0.8, group = factor(Group))) +
  geom_bar(stat = "identity",color="black", size=0.5, alpha=0.75,
  aes(fill = factor(Group)),
  position = position_dodge2(preserve = "single")) +
  #ylim(-2.0, 8.0) +
  ylab("Log fold change (natural log)") +
  geom_errorbar(limits, width = 0.2, size=0.5, position = position_dodge(width = 0.9)) +
  coord_flip() +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, color="black"), text = element_text(size=11), axis.text.y = element_text(color="black")) +
  theme(panel.background = element_rect(fill = "white"),
        legend.key  = element_rect(fill = "white"),
        axis.line.x = element_line(colour = "black", size = 0.5),
        axis.line.y = element_line(colour = "black", size = 0.5)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  geom_hline(yintercept=0, linetype="dotted", color = "black", size=1) +
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=0, hjust=0.5)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y = element_blank()) +
    facet_wrap(~Day, scales = "free_y", ncol = 2, nrow = 3) +
  theme(axis.text.y = element_text(color="black", size=12))

ggsave("F.pdf", plot1, width=5, height=7.2, units="in", scale=3)






limits <- aes(ymax = `FC` + stderr, ymin = `FC` - stderr)


plot2 <- ggplot(data = M_data_df ,
  aes(x = feature, y = `FC`, width=0.8, group = factor(Group))) +
  geom_bar(stat = "identity",color="black", size=0.5, alpha=0.75,
  aes(fill = factor(Group)),
  position = position_dodge2(preserve = "single")) +
  #ylim(-2.0, 8.0) +
  #ylim(-2.0, 8.0) +
  ylab("Log fold change (natural log)") +
  geom_errorbar(limits, width = 0.2, size=0.5, position = position_dodge(width = 0.9)) +
  coord_flip() +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, color="black"), text = element_text(size=11), axis.text.y = element_text(color="black")) +
  theme(panel.background = element_rect(fill = "white"),
        legend.key  = element_rect(fill = "white"),
        axis.line.x = element_line(colour = "black", size = 0.5),
        axis.line.y = element_line(colour = "black", size = 0.5)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  geom_hline(yintercept=0, linetype="dotted", color = "black", size=1) +
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=0, hjust=0.5)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y = element_blank()) +
    facet_wrap(~Day, scales = "free_y", ncol = 2, nrow = 3) +
  theme(axis.text.y = element_text(color="black", size=12))

ggsave("M.pdf", plot2, width=5, height=7.2, units="in", scale=3)









data_df <- read.table("D://ABX_LOCO_2024/PWY_FC.txt", header=T, row.names=1, check.names = FALSE)

data_df <- subset(data_df, abs(FC) > 2.0)



F_data_df <- subset(data_df, abs(FC) > 2.0 & Sex=="Female")

M_data_df <- subset(data_df, abs(FC) > 2.0 & Sex=="Male")


limits <- aes(ymax = `FC` + stderr2, ymin = `FC` - stderr2)


plot3 <- ggplot(data = data_df,
                aes(x = feature, y = `FC`, width=0.8, group = factor(Group))) +
  geom_bar(stat = "identity",color="black", size=0.5, alpha=0.75,
           aes(fill = factor(Group)),
           position = position_dodge2(preserve = "single")) +
  #ylim(-2.0, 8.0) +
  #ylim(-2.0, 8.0) +
  ylab("Log fold change (natural log)") +
  geom_errorbar(limits, width = 0.2, size=0.5, position = position_dodge(width = 0.9)) +
  #coord_flip() +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, color="black"), text = element_text(size=11), axis.text.y = element_text(color="black")) +
  theme(panel.background = element_rect(fill = "white"),
        legend.key  = element_rect(fill = "white"),
        axis.line.x = element_line(colour = "black", size = 0.5),
        axis.line.y = element_line(colour = "black", size = 0.5)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  geom_hline(yintercept=0, linetype="dotted", color = "black", size=1) +
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=0, hjust=0.5)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y = element_blank()) +
  facet_grid(Day~Sex, space = "free") +
  theme(axis.text.y = element_text(color="black", size=12))

ggsave("F_M_3.pdf", plot3, width=5, height=7.2, units="in", scale=3)

























  #ylim(-2.0, 8.0) +
  ylab("Log fold change (natural log)") +
  geom_errorbar(limits, width = 0.2, size=1, position = position_dodge(width = 0.9)) +
  coord_flip() +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, color="black"), text = element_text(size=11), axis.text.y = element_text(color="black")) +
  theme(panel.background = element_rect(fill = "white"),
        legend.key  = element_rect(fill = "white"),
        axis.line.x = element_line(colour = "black", size = 0.5),
        axis.line.y = element_line(colour = "black", size = 0.5)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  geom_hline(yintercept=0, linetype="dotted", color = "black", size=1) +
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=0, hjust=0.5)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y = element_blank()) +
  facet_wrap(~Day, scales = "free_y", ncol = 2, nrow = 3) +
  theme(axis.text.y = element_text(color="black", size=12))

ggsave("M.pdf", plot1, width=7.2, height=5, units="in", scale=3)




















F_p <- ggplot(F_PWY_df_Day5, aes(x=LF, y=-log10(qval), col=factor(Group))) + 
  geom_point(aes()) + 
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  geom_vline(xintercept = 0.0) +
  #theme_classic(base_size = 10) +
  #scale_y_continuous(limits = c(0.0,200), labels = scales::number_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = 4)) +
  #facet_wrap(~as.factor(Day), ncol = 7) +
  theme(axis.title.x=element_blank(),
        legend.title= element_blank())
#strip.text.x = element_blank())

F_p + geom_text_repel(aes(label = feature),
                      size=2.0, box.padding = 0.75, max.overlaps = 10)












































M_p <- ggplot(data_dif_M, aes(x=LF, y=-log10(qval), col=factor(Group))) + 
  geom_point(aes()) + 
  scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  geom_vline(xintercept = 0.0) +
  #theme_classic(base_size = 10) +
  #scale_y_continuous(limits = c(0.0,200), labels = scales::number_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = 4)) +
  facet_wrap(~as.factor(Day), ncol = 7) +
  theme(axis.title.x=element_blank(),
        legend.title= element_blank())
#strip.text.x = element_blank())

M_p + geom_text_repel(aes(label = feature),
                      size=2.0, box.padding = 0.75, max.overlaps = 10)






#geom_point(aes(shape=factor(NULL))) +






geom_point(position = position_jitterdodge(jitter.width = 0.333, dodge.width = 0.85), aes(fill = Group), pch = 21) +
 scale_color_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  




data_dif$Significant <- ifelse(data_dif$qval < 0.05, "FDR < 0.05", "Not Sig")
ggplot(data_dif, aes(x = LF, y = log10(qval))) +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("red", "grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  geom_text_repel(
    data = subset(data_dif, qval < 0.05),
    aes(label = rownames(data_dif)),
    size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )














data_df <- read.table("D://ABX_LOCO_2024/PWY_FC_2.txt", header=T, row.names=1, check.names = FALSE)

library(dplyr)


F_PWY_df_Day1 <- subset(data_df, abs(FC) > 2.0 & Sex=="Female" & Day=="A1")
F_PWY_df_Day2 <- subset(data_df, abs(FC) > 2.0 & Sex=="Female" & Day=="B2")
F_PWY_df_Day5 <- subset(data_df, abs(FC) > 2.0 & Sex=="Female" & Day=="C5")
F_PWY_df_Day8 <- subset(data_df, abs(FC) > 2.0 & Sex=="Female" & Day=="D8")
F_PWY_df_Day12 <- subset(data_df, abs(FC) > 2.0 & Sex=="Female" & Day=="E12")
F_PWY_df_Day65 <- subset(data_df, abs(FC) > 2.0 & Sex=="Female" & Day=="G65")
F_PWY_df_Day100 <- subset(data_df, abs(FC) > 2.0 & Sex=="Female" & Day=="H100")





library(EnhancedVolcano)




EnhancedVolcano(F_PWY_df_Day5,
                lab = rownames(F_PWY_df_Day2),
                x = 'FC',
                y = 'qval',
                #pCutoff = 0.1,
                title = NULL,
                subtitle = NULL,
                pointSize = 3.0,
                labSize = 3,
                col=c('gray', 'blue', 'green', 'red3'),
                colAlpha = 1,
                drawConnectors = TRUE,
                typeConnectors='closed',
                endsConnectors='last',
                legendPosition = 'blank',
                max.overlaps=15,
                pCutoff = 0.05,
                FCcutoff = 0.0)









  ggplot2::coord_cartesian(xlim=c(-2, 2)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-8,10, 1)) +
  ggplot2::coord_cartesian(ylim=c(0, 8)) +
  ggplot2::scale_y_continuous(
    breaks=seq(0,175, 25)) +
  theme(axis.text=element_text(colour="black"))
















EnhancedVolcano(data_dif,
                lab = rownames(data_dif),
                x = 'LF', #log2(exp(coef))
                y = 'qval')






EnhancedVolcano(data_dif,
                lab = rownames(data_dif),
                x = 'LF',
                y = 'qval',
                title = 'N061011 versus N61311',
                #pCutoff = 10e-16,
                #FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 6.0,
                col=c('black', 'blue', 'green', 'red3'),
                colAlpha = 1)



p1 <- EnhancedVolcano(data_dif,
                      lab = rownames(data_dif),
                      x = "LF",
                      y = "qval",
                      #pCutoff = 0.1,
                      max.overlaps = 40,
                      FCcutoff = 0,
                      #ylim = c(0, -log10(10e-12)),
                      #pointSize = c(ifelse(data_dif$coef>2, 8, 1)),
                      pointSize = 3.5,
                      labSize = 6.0,
                      #shape = c(6, 6, 19, 16),
                      title = "Saline vs Cocaine",
                      subtitle = "Differential abundance",
                      #caption = bquote(~Log[2]~ "fold change cutoff, 2; p-value cutoff, 10e-4"),
                      legendPosition = "right",
                      legendLabSize = 14,
                      col = c("grey30", "forestgreen", "royalblue", "red2"),
                      colAlpha = 0.9,
                      #drawConnectors = TRUE,
                      #hline = c(10e-8),
                      widthConnectors = 0.5)

p1








p1 <- EnhancedVolcano(data_dif,
                      lab = rownames(data_dif),
                      x = "LF",
                      y = "qval",
                      pCutoff = 0.1,
                      #FCcutoff = 2,
                      #ylim = c(0, -log10(10e-12)),
                      pointSize = c(ifelse(data_dif$LF>2, 8, 1)),
                      labSize = 3.0,
                      #shape = c(6, 6, 19, 16),
                      title = "Saline vs Cocaine in Control Mice",
                      subtitle = "Differential abundance",
                      #caption = bquote(~Log[2]~ "fold change cutoff, 2; p-value cutoff, 10e-4"),
                      legendPosition = "right",
                      legendLabSize = 14,
                      colAlpha = 0.9,
                      #colGradient = c('red3', 'royalblue'),
                      #drawConnectors = TRUE,
                      hline = c(10e-8),
                      widthConnectors = 0.5)

p1






EnhancedVolcano(data_dif,
                lab = rownames(data_dif),
                x = 'LF',
                y = 'pval',
                title = 'N061011 versus N61311',
                pCutoff = 0.1,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)


























p1 <- EnhancedVolcano(data_dif,
                      lab = rownames(data_dif),
                      x = "LF",
                      y = "qval",
                      pCutoff = 0.1,
                      #FCcutoff = 2,
                      #ylim = c(0, -log10(10e-12)),
                      pointSize = c(ifelse(data_dif$coef>2, 8, 1)),
                      labSize = 3.0,
                      #shape = c(6, 6, 19, 16),
                      title = "Saline vs Cocaine in Control Mice",
                      subtitle = "Differential abundance",
                      #caption = bquote(~Log[2]~ "fold change cutoff, 2; p-value cutoff, 10e-4"),
                      legendPosition = "right",
                      legendLabSize = 14,
                      colAlpha = 0.9,
                      #colGradient = c('red3', 'royalblue'),
                      drawConnectors = TRUE,
                      hline = c(10e-8),
                      widthConnectors = 0.5)

p1




library(EnhancedVolcano)

EnhancedVolcano(data_dif,
                lab = rownames(data_dif),
                x = 'coef',
                y = 'pval',
                pCutoff = 0.1,
                #FCcutoff = 25,
                pointSize = 3.0,
                labSize = 6.0)



EnhancedVolcano(data_dif,
                lab = rownames(data_dif),
                x = 'LF',
                y = 'pval',
                title = 'N061011 versus N61311',
                #pCutoff = 0.1,
                #FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)
























#xlim=c(0, 2000)

dev.off()
warnings()

library(rstatix)
#value <- as.data.frame(log10(meta_ob1$CPmg))
#identify_outliers(value)
value <- as.data.frame(meta_ob1$LC)
identify_outliers(value)


model <- aov(formula = LC ~ Cohort, data = meta_ob1)
TukeyHSD(model, conf.level=.95)

#library(car)
Anova(model, type="II")



model <- aov(formula = log10(CPmg) ~ Group + Sex, data = meta_ob1)
TukeyHSD(model, conf.level=.95)

library(car)
Anova(model, type="II")




meta_ob1 <- meta_ob1[meta_ob1$Unique!="IA2_F_1",] #3264524583
#meta_ob1 <- meta_ob1[meta_ob1$Unique!="GA4_M_1",] #3086845492


nrow(meta_ob1)


library(outliers)
test <- grubbs.test(log10(meta_ob1$CPmg))
test

library(outliers)
test <- grubbs.test(log10(meta_ob1$Diff))
test


library(lme4)
library(lmerTest)
library(car)
library(caret)
library(lmtest)
#model <- lmer(LC ~ meta_ob1$Cohort + (1|meta_ob1$Mouse_ID), data = meta_ob1)
#model <- lmer(LC ~ meta_ob1$Group*as.factor(meta_ob1$Day)*Sex + (1|meta_ob1$Mouse_ID), data = meta_ob1)
#model <- lmer(LC ~ meta_ob1$Group*as.factor(meta_ob1$Day)*Sex*Cohort + (1|meta_ob1$Mouse_ID), data = meta_ob1)


model1 <- lmer(LC ~ Group*as.factor(meta_ob1$Day)*Sex + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=FALSE)

#model2 <- lmer(CPmg ~ Group:as.factor(meta_ob1$Day):Sex:Cohort + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=FALSE)
#model3 <- lmer(CPmg ~ Group:as.factor(meta_ob1$Day):Sex + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=FALSE)

#anova(model1, model2, model3)


library(hnp)
d <- function(obj) resid(obj, type="pearson")
s <- function(n, obj) simulate(obj)[[1]]
f <- function(y.) refit(model1, y.)

hnp(model1, newclass=TRUE, diagfun=d, simfun=s, fitfun=f)






library(fitdistrplus)
descdist(meta_ob1$LC, discrete = FALSE)

normal_dist <- fitdist(meta_ob1$LC, "norm")
plot(normal_dist)




dev.off()

library(dplyr)
#meta_ob1 %>% count(Sex, Group, Day, Cohort)

means <- meta_ob1 %>%
group_by(Sex, Group, Day) %>%
summarise_at(vars(CPmg), list(name = mean)) %>%
print(n=100)
write.csv(means, file="CPmg_means.csv")


distBCMod <- caret::BoxCoxTrans(meta_ob1$LC)
print(distBCMod)
meta_ob1 <- cbind(meta_ob1, LC_new=predict(distBCMod, meta_ob1$LC)) # append the transformed variable to cars
head(meta_ob2)
model4 <- lmer(LC_new ~ Group:as.factor(meta_ob1$Day):Sex + Cohort + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=FALSE)

anova(model1, model2, model3, model4)




plot(model1)

hist((resid(model1) - mean(resid(model1))) / sd(resid(model1)), freq = FALSE); curve(dnorm, add = TRUE)

plot(ranef(model1))

lattice::dotplot(ranef(model1, condVar = TRUE))



dev.off()



Anova(model1, type=3)
summary(model1)
#vcov(model)




library(emmeans)
options(max.print=1000000)
options(digits = 5)
emm_options(opt.digits = FALSE)
emmeans(model1, list(pairwise ~ Group*as.factor(Day)*Sex), adjust = "tukey")


library(ggpubr)
ggdensity(meta_ob1$LC, 
          main = "Density plot of LC",
          xlab = "LC")



ggqqplot(meta_ob1$CPmg)
         
qqPlot(meta_ob1$LC)
       
shapiro.test(meta_ob1$CPmg)


#meta_ob1 <- na.omit(meta_ob1)

#meta_ob1 <- meta_ob1[meta_ob1$Day!="2" & meta_ob1$Day!="5"  & meta_ob1$Day!="12" & meta_ob1$Day!="51"  & meta_ob1$Day!="65" & meta_ob1$Sex=="M" ,]
#meta_ob1 <- meta_ob1[meta_ob1$Day!="1" & meta_ob1$Day!="8"  & meta_ob1$Day!="12" & meta_ob1$Day!="51"  & meta_ob1$Day!="65" & meta_ob1$Sex=="F" ,]

meta_ob1 <- meta_ob1[meta_ob1$Unique!="GA4_M_1" & meta_ob1$Unique!="IA2_F_1",]

#meta_ob1 <- meta_ob1[meta_ob1$Day=="2" & meta_ob1$Sex=="F",]






p <- ggplot(DATA, (aes(x=DBH, y=Epiphyte, color=Location, shape=Location))) + theme_classic() +ylab("Epiphyte Abundance") + xlab("DBH (cm)")
p + geom_point() + geom_smooth(method=lm, se=FALSE, fullrange=TRUE)








library(ggplot2)


min.mean.se.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


copies_plot_1 <- ggplot(meta_ob1, aes(x=as.factor(Day), y=LC, color=Group)) +
  #geom_boxplot(outlier.shape = NA, alpha=0, fatten=NA) +
  #stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
              #width = 0.75, size = 1, linetype = "solid")+
  stat_summary(fun.data = min.mean.se.max, geom = "boxplot") +
  geom_point(position=position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=1, alpha=0.7, stroke=1.25) +
  scale_shape_manual(values = c(1, 2, 0)) +
  scale_color_manual(values=c("navyblue", "red3", "darkgreen","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("navyblue", "red3", "darekgreen","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  facet_wrap(~Sex, scales = "free_x", ncol = 10) +
  ylab("Log 16S rRNA gene copies per mg of feces (N=360)") +
  theme_classic() +
  theme(axis.title.y = element_text(size=9)) 


#plot_data <- ggplot_build(copies_plot_1)$plot$data

#write.csv(plot_data, file="plot_data.csv")









ggplot(meta_ob1, aes(x = as.factor(Day), y = LC, fill = Group)) +
geom_boxplot(position = position_dodge(width = 0.7)) +
geom_point(position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.7), aes(fill = Group), pch = 21) +
facet_wrap(~Sex, scales = "free_x", ncol = 10) +
ylab("Log 16S rRNA gene copies per mg of feces (N=351)") +
theme_classic() +
theme(axis.title.y = element_text(size=9)) 







library(ggplot2)


meta_ob1 <- meta_ob1[meta_ob1$ID!="H2O" & meta_ob1$Global_LC_Extreme=="No" & meta_ob1$Day!="8" & meta_ob1$Day!="12" & meta_ob1$Day!="51" & meta_ob1$Day!="65",]


ggplot(meta_ob1, (aes(x=Day, y=LC, color=Group))) +
 theme_classic() +
ylab("Log 16S rRNA gene copies per mg of feces (N=193)") + xlab("Days") +
 geom_point(aes(fill = Group), pch = 1, size=5) +
 geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  facet_wrap(~Sex, scales = "free_x", ncol = 1)
  

library(lme4)
library(lmerTest)
library(car)
model1 <- lmer(LC ~ Group*Sex +(1|meta_ob1$Mouse_ID), data = meta_ob1, REML=FALSE)
Anova(model1)
summary(model)


options(max.print=1000000)
options(digits = 5)
emm_options(opt.digits = FALSE)
emmeans(ancova_model, list(pairwise ~ Group*Sex), adjust = "tukey")




library(multcomp)
postHocs <- glht(ancova_model, linfct = mcp(Group = "Tukey"))
summary(postHocs)


dev.off()









setwd("/Users/15869/Desktop/Kuhn_Lab/CEF_qPCR/")

meta_ob1  <- read.table("/Users/15869/Desktop/Kuhn_Lab/CEF_qPCR/CEF_qPCR.txt", header=T, row.names=1)

nrow(meta_ob1)

copies_plot_2 <- ggplot(meta_ob1, aes(x=as.factor(Day), y=log10(normCP3ul), color=Group, shape=Sex)) +
  #geom_boxplot(outlier.shape = NA, alpha=0, fatten=NA) +
  #stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
  #width = 0.75, size = 1, linetype = "solid")+
  stat_summary(fun.data = min.mean.se.max, geom = "boxplot") +
  geom_point(position=position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=1, alpha=0.7, stroke=1.25) +
  scale_shape_manual(values = c(1, 2)) +
  scale_color_manual(values=c("navyblue", "red3", "darkgreen","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("navyblue", "red3", "darekgreen","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  #facet_wrap(~Sex, scales = "free_x", ncol = 5) +
  ylab("log10(16S rRNA gene copies per reaction) N=240") +
  theme_classic() +
  theme(axis.title.y = element_text(size=9))



library(ggpubr)
merged_plot <- ggarrange(copies_plot_1, copies_plot_2,
                          common.legend = TRUE, legend = "right",
                          ncol = 1, nrow = 2, align="hv")












setwd("/Users/15869/Desktop/Kuhn_Lab/CEF_qPCR/")

meta_ob1  <- read.table("/Users/15869/Desktop/Kuhn_Lab/CEF_qPCR/CEF_qPCR_normalized.txt", header=T, row.names=1)

#meta_ob1 <- na.omit(meta_ob1)

meta_ob1 <- meta_ob1[meta_ob1$Day!="2" & meta_ob1$Day!="5"  & meta_ob1$Day!="12" & meta_ob1$Day!="51"  & meta_ob1$Day!="65" ,]

nrow(meta_ob1)




#model <- aov(formula = log10(CPmg) ~ ID, data = meta_ob1)
model <- aov(formula = log10(CPmg) ~ ID, data = meta_ob1)
TukeyHSD(model, conf.level=.95)

library(car)
Anova(model, type="II")





#model <- aov(formula = log10(CPmg) ~ ID, data = meta_ob1)
model <- aov(formula = log10(normCP3ul) ~ as.factor(Day) + Group + Sex, data = meta_ob1)
TukeyHSD(model, conf.level=.95)

library(car)
Anova(model, type="II")



df <- meta_ob1 


df <- as.data.frame(cbind(df$Day, input_metadata$CPP))
colnames(df)[1] <- "PWY138"







input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)





#colnames(input_data) == row.names(input_metadata)





fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    min_prevalence=0.1,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "M_Clusters_PWY_CPP", 
                    fixed_effects  = c("ET", "CPP"),
                    reference      = c("ET,C1"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)









library(multcomp)
model <- aov(CPP ~ log10(CPmg), data = meta_ob1)
postHocs <- glht(ancova_model, linfct = mcp(ID4 = "Tukey"))
summary(postHocs)







library(car)
ancova_model <- aov(log10(CPmg) ~ ID2, data = meta_ob1)
Anova(ancova_model, type="II")




female_chao_mod <- aov(log10(CPmg) ~ ID2, data = meta_ob1)
TukeyHSD(female_chao_mod, conf.level=.95)




str(CPP_slopes)


names(CPP_slopes)

CPP_slopes$test
CPP_slopes_out <- rbind(CPP_slopes$test)


write.csv(CPP_slopes_out, file="CPP_slopes.csv")



























copies_plot <- ggplot(meta_ob1, aes(x=ID, y=log10(CPmg), color=Group, shape=NULL)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  #stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
  #width = 0.75, size = 1, linetype = "solid")+
  geom_point(position=position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.7, stroke=1.25) +
  scale_shape_manual(values = c(1, 2)) +
  scale_color_manual(values=c("navyblue", "red3", "darkgreen","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("navyblue", "red3", "darekgreen","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  facet_wrap(~ID2, scales = "free_x", ncol = 6) +
  ylab("log10(16S rRNA gene copies per mg of feces)") +
  theme_classic()



dev.off()


df <- meta_ob1[meta_ob1$ID2== "IA_M",]

pairwise.t.test(log10(df$CPmg), df$ID, p.adjust.method = "none", paired = TRUE)





out <- meta_ob1 %>% 
  group_by(ID)  %>% 
  summarise(as_tibble(rbind(summary(CPmg))))

summary_stats <- select(out, ID, Mean)

write.table(summary_stats, file='/Users/zebra/Desktop/CEF_qPCR/summary_stats.txt', sep='\t')








meta_ob1 <- meta_ob1[meta_ob1$Day!= "2" & meta_ob1$Day!= "5" & meta_ob1$Sex== "F",]





library(lme4)
library(lmerTest)
library(car)
model <- lmer (log10(meta_ob1$CPmg) ~ meta_ob1$Group + meta_ob1$Day + (1|Mouse_ID), data = meta_ob1)
Anova(model)
summary(model)


library(emmeans)
emmeans(model, list(pairwise ~ Group), adjust = "tukey")





meta_ob1 <- meta_ob1[meta_ob1$Day!= "2" & meta_ob1$Day!= "5",]


model <- aov(log10(meta_ob1$CPmg)~ Group + Sex, data = meta_ob1)
summary(model)



meta_ob1





library(agricolae)


meta_ob1 <- meta_ob1[meta_ob1$Day== "51" & meta_ob1$Sex== "F",]




plant.lm <- lm(log10(meta_ob1$CPmg) ~ Group, data = meta_ob1)
plant.av <- aov(plant.lm)
summary(plant.av)

tukey.test <- TukeyHSD(plant.av)
tukey.test


plot(tukey.test)






summary_out <- tapply(log10(meta_ob1$CPmg), meta_ob1$ID, summary)  


library(tidyverse)
out <- meta_ob1 %>% 
  group_by(ID) %>% 
  summarize(vals = list(summary(log10(meta_ob1$CPmg)))) %>% 
  unnest_wider(vals)

out <- as.data.frame(out)


    
    







    
    
    
ggplot(meta_ob1 , aes(x=ID, y=CPmg)) +
  #geom_point(aes(colour = meta_phy_out1$Sex)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_jitter(aes(color = NULL, shape=NULL), width = 0.1, size = 3) +
  scale_color_manual(values=c("navyblue", "red3", "green", "magenta", "yellow","blue","black","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("navyblue", "red3", "green","magenta", "yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  theme_bw() +
  theme(axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  facet_wrap(~Group, scales = "free_x", ncol = 2)
#facet_wrap(~factor(meta_phy_out1$ID3, c('F_Control', 'F_ABX', 'M_Control', 'M_ABX')))














plot_data <- ggplot_build(copies1_plot)$plot$data



pairwise.t.test(plot_data$copies, plot_data$Group, p.adjust.method = "none", paired = FALSE)



C_SS <- plot_data[plot_data$Group=="C",]
GA_SS <- plot_data[plot_data$Group=="GA",]
IA_SS <- plot_data[plot_data$Group=="IA",]


shapiro.test(C_SS$copies)
shapiro.test(GA_SS$copies)
shapiro.test(IA_SS$copies)





#devtools::install_github("vmikk/metagMisc")
library(phyloseq)
library(metagMisc)
library(vegan)
library(ggplot2)






setwd("/Users/zebra/Desktop/GF/")

meta_ob1  <- read.table("/Users/zebra/Desktop/GF/meta.txt", header=T, row.names=1)
#tax_ob1 <- import_mothur(mothur_constaxonomy = ("/Users/zebra/Desktop/GF/taxa.txt"))
tax_ob1 <- import_mothur(mothur_constaxonomy = ("/Users/zebra/Desktop/GF/PWY_taxa.txt"))



#meta_ob1  <- subset(meta_ob1, Type!="Other" ,)
meta_ob1  <- subset(meta_ob1, Type=="GF_27" ,)



unique(meta_ob1$Type)



nrow(meta_ob1)






#Alpha diversity


chao1_plot <- ggplot(meta_ob1, aes(x=Type, y=meta_ob1$chao, color=Type, shape=NULL)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.333, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.7, stroke=1.25) +
  scale_shape_manual(values = c(22, 22)) +
  scale_color_manual(values=c("navyblue", "red3", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("navyblue", "red3", "red", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  #facet_wrap(~meta_ob1$Sex, scales = "free_x", nrow = 1) +
  theme_classic()
 
shannon_plot <- ggplot(meta_ob1, aes(x=Type, y=meta_ob1$shannon, color=Type, shape=NULL)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.333, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.7, stroke=1.25) +
  scale_shape_manual(values = c(21, 21)) +
  scale_color_manual(values=c("navyblue", "red3", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("navyblue", "red3", "red", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  #facet_wrap(~meta_ob1$Sex, scales = "free_x", nrow = 1) +
  theme_classic()

invsimpson_plot <- ggplot(meta_ob1, aes(x=Type, y=meta_ob1$invsimpson, color=Type, shape=NULL)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.333, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.7, stroke=1.25) +
  scale_shape_manual(values = c(21, 21)) +
  scale_color_manual(values=c("navyblue", "red3", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("navyblue", "red3", "red", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  #facet_wrap(~meta_ob1$Sex, scales = "free_x", nrow = 1) +
  theme_classic()



library(ggpubr)
alpha_merged <- ggarrange(CPP_plot, chao1_plot, shannon_plot, invsimpson_plot ,
                       common.legend = TRUE, legend = "right",
                       ncol = 1, nrow = 4, align="hv")


alpha_merged <- ggarrange(chao1_plot, shannon_plot, invsimpson_plot ,
                          common.legend = TRUE, legend = "right",
                          ncol = 1, nrow = 3, align="hv")


plot_data <- ggplot_build(invsimpson_plot)$plot$data

#write.csv(plot_data ,file="plot_data.csv")


as.vector()





pairwise.t.test(plot_data$chao, plot_data$Type, p.adjust.method = "none", paired = FALSE)

pairwise.t.test(plot_data$shannon, plot_data$Type, p.adjust.method = "none", paired = FALSE)

pairwise.t.test(plot_data$invsimpson, plot_data$Type, p.adjust.method = "none", paired = FALSE)















ttest_out_CPP <- pairwise.t.test(plot_data$CPP, plot_data$ID4, p.adjust.method = "none", paired = FALSE)
ttest_out_CPP <- rbind(ttest_out_CPP$p.value, ttest_out_CPP$p.adjust.method, ttest_out_CPP$method, ttest_out_CPP$data.name)





M <- as.matrix(plot_data$ID4)

do.call(expand.grid,split(M,rep(1:nrow(M),ncol(M))))

M1<-matrix(rpois(10,2),ncol=5)
do.call(expand.grid,split(M1,rep(1:nrow(M1),ncol(M1))))




V1 <- unique(plot_data$ID4)
V2 <- unique(plot_data$ID4)

out <- expand.grid(V1, V2)

corr[ttest_out_CPP$p.value < 0.05]

write.csv(out, file="out.csv")



ANCOVA_meta <- read.table("/Users/zebra/Desktop/GF/CPP_ANCOVA.txt", header=T, row.names=1)

library(car)
ancova_model <- aov(ANCOVA_meta ~ ID4 + SS, data = data)
Anova(ancova_model, type="III")


p1 <- ggplot(ANCOVA_meta, aes(SS, SC, colour = ID4)) + geom_point(size = 3) + theme(legend.position="top")


ggplot(data = ANCOVA_meta, mapping = aes(y = CPP, x = ID4, group = ID4, colour=ID4, shape=ID4)) +
  geom_line(linewidth=1.2) +
  geom_point(size=5) 




library(car)
ancova_model <- aov(CPP ~ ID4, data = meta_ob1)
Anova(ancova_model, type="III")

str(ancova_model)

meta_ob1$ID4 <- as.factor(meta_ob1$ID4)

library(multcomp)
ancova_model <- aov(CPP ~ ID4, data = meta_ob1)
postHocs <- glht(ancova_model, linfct = mcp(ID4 = "Tukey"))
CPP_slopes <- summary(postHocs)

str(CPP_slopes)


names(CPP_slopes)

CPP_slopes$test
CPP_slopes_out <- rbind(CPP_slopes$test)


write.csv(CPP_slopes_out, file="CPP_slopes.csv")




wilcox_out_CPP <- pairwise.wilcox.test(plot_data$CPP, plot_data$ID4, p.adjust.method = "none", paired = FALSE)
wilcox_out_CPP <- rbind(wilcox_out_CPP$p.value, wilcox_out_CPP$p.adjust.method, wilcox_out_CPP$method, wilcox_out_CPP$data.name)

wilcox_out_chao <- pairwise.wilcox.test(plot_data$chao, plot_data$ID4, p.adjust.method = "none", paired = FALSE)
wilcox_out_chao <- rbind(wilcox_out_chao$p.value, wilcox_out_chao$p.adjust.method, wilcox_out_chao$method, wilcox_out_chao$data.name)

wilcox_out_shannon <- pairwise.wilcox.test(plot_data$shannon, plot_data$ID4, p.adjust.method = "BH", paired = FALSE)
wilcox_out_shannon <- rbind(wilcox_out_shannon$p.value, wilcox_out_shannon$p.adjust.method, wilcox_out_shannon$method, wilcox_out_shannon$data.name)

wilcox_out_invsimpson <- pairwise.wilcox.test(plot_data$invsimpson, plot_data$ID4, p.adjust.method = "BH", paired = FALSE)
wilcox_out_invsimpson <- rbind(wilcox_out_invsimpson$p.value, wilcox_out_invsimpson$p.adjust.method, wilcox_out_invsimpson$method, wilcox_out_invsimpson$data.name)


wilcox_out_merged <- rbind(wilcox_out_CPP, wilcox_out_chao, wilcox_out_shannon, wilcox_out_invsimpson)



write.csv(wilcox_out_CPP, file="wilcox_out_CPP.csv")


F_A_Control_A_SS <- plot_data[plot_data$ID2=="F_A_Control_A_SS",]

F_A_Control_B_SC <- plot_data[plot_data$ID2=="F_A_Control_B_SC",]

F_B_ABX_A_SS <- plot_data[plot_data$ID2=="F_B_ABX_A_SS",]

F_B_ABX_B_SC <- plot_data[plot_data$ID2=="F_B_ABX_B_SC",]

M_A_Control_A_SS <- plot_data[plot_data$ID2=="M_A_Control_A_SS",]

M_A_Control_B_SC <- plot_data[plot_data$ID2=="M_A_Control_B_SC",]

M_B_ABX_A_SS <- plot_data[plot_data$ID2=="M_B_ABX_A_SS",]

M_B_ABX_B_SC <- plot_data[plot_data$ID2=="M_B_ABX_B_SC",]



shapiro.test(F_A_Control_A_SS$CPP)
shapiro.test(F_A_Control_B_SC$CPP)
shapiro.test(F_B_ABX_A_SS$CPP)
shapiro.test(F_B_ABX_B_SC$CPP)
shapiro.test(M_A_Control_A_SS$CPP)
shapiro.test(M_A_Control_B_SC$CPP)
shapiro.test(M_B_ABX_A_SS$CPP)
shapiro.test(M_B_ABX_B_SC$CPP)

shapiro.test(F_A_Control_A_SS$chao)
shapiro.test(F_A_Control_B_SC$chao)
shapiro.test(F_B_ABX_A_SS$chao)
shapiro.test(F_B_ABX_B_SC$chao)
shapiro.test(M_A_Control_A_SS$chao)
shapiro.test(M_A_Control_B_SC$chao)
shapiro.test(M_B_ABX_A_SS$chao)
shapiro.test(M_B_ABX_B_SC$chao)

shapiro.test(F_A_Control_A_SS$shannon)
shapiro.test(F_A_Control_B_SC$shannon)
shapiro.test(F_B_ABX_A_SS$shannon)
shapiro.test(F_B_ABX_B_SC$shannon)
shapiro.test(M_A_Control_A_SS$shannon)
shapiro.test(M_A_Control_B_SC$shannon)
shapiro.test(M_B_ABX_A_SS$shannon)
shapiro.test(M_B_ABX_B_SC$shannon)

shapiro.test(F_A_Control_A_SS$invsimpson)
shapiro.test(F_A_Control_B_SC$invsimpson)
shapiro.test(F_B_ABX_A_SS$invsimpson)
shapiro.test(F_B_ABX_B_SC$invsimpson)
shapiro.test(M_A_Control_A_SS$invsimpson)
shapiro.test(M_A_Control_B_SC$invsimpson)
shapiro.test(M_B_ABX_A_SS$invsimpson)
shapiro.test(M_B_ABX_B_SC$invsimpson)



female_data <- rbind(F_B_ABX_A_SS, F_B_ABX_B_SC)

female_data <- rbind(F_A_Control_A_SS, F_A_Control_B_SC, F_B_ABX_A_SS, F_B_ABX_B_SC)

male_data <- rbind(M_A_Control_A_SS, M_A_Control_B_SC, M_B_ABX_A_SS, M_B_ABX_B_SC)


female_chao_mod <- aov(formula = chao ~ ID2, data = female_data)
TukeyHSD(female_chao_mod, conf.level=.95)

female_shannon_mod <- aov(formula = shannon ~ ID2, data = female_data)
TukeyHSD(female_shannon_mod, conf.level=.95)

female_invsimpson_mod <- aov(formula = invsimpson ~ ID2, data = female_data)
TukeyHSD(female_invsimpson_mod, conf.level=.95)



male_chao_mod <- aov(formula = chao ~ ID2, data = male_data)
TukeyHSD(male_chao_mod, conf.level=.95)

male_shannon_mod <- aov(formula = shannon ~ ID2, data = male_data)
TukeyHSD(male_shannon_mod, conf.level=.95)

male_invsimpson_mod <- aov(formula = invsimpson ~ ID2, data = male_data)
TukeyHSD(male_invsimpson_mod, conf.level=.95)







library(multcomp)



post.hoc <- glht(male_chao_mod, linfct = mcp(ID2 = 'Tukey'))






summary(glht(mod1, linfct = mcp(ID2 = "Tukey")), test = adjusted("BH"))


post.hoc <- glht(model, linfct = mcp(reatment = 'Tukey'))


model <- lm(chao ~ Type, data = plot_data)


















#Beta diversity





#write.csv(shared_phy_out1, file ="F_ABX_picrust_shared.csv")



dev.off()



meta_ob1  <- read.table("/Users/zebra/Desktop/GF/meta.12523.txt", header=T, row.names=1)
#meta_ob1 <- meta_ob1[meta_ob1$Sample_ID!= "wk1" & meta_ob1$Sample_ID!= "wk2" & meta_ob1$Type2=="test" & meta_ob1$Type!="B_ABX",]
meta_ob1 <- meta_ob1[meta_ob1$Sample_ID!= "wk1" & meta_ob1$Sample_ID!= "wk2" & meta_ob1$Type2=="test",]
meta_ob1 <- meta_ob1[meta_ob1$ID3=="F_A_Control",]




#Get sample names
GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("shared.12523.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))



merged_phy_percent1 <- phyloseq_standardize_otu_abundance(phyloseqobj.f, method = "total")


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


meta_phy_out1 <- pssd2veg(merged_phy_percent1)
shared_phy_out1 <- psotu2veg(merged_phy_percent1)




#dist_bray1 <- capscale(shared_phy_out1~1, distance="bray", binary=TRUE, eig=TRUE)
dist_bray1 <- capscale(shared_phy_out1~1, distance="bray", binary=FALSE, eig=TRUE)


prop_explained <- summary(eigenvals(dist_bray1))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)


data.scores = scores(dist_bray1, display = "sites", choices = c(1:3))
#summary(data.scores)


df <- data.frame(data.scores)

library(glue)
labx1 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby1 <- c(glue("PCo 2 ({prop_explained [2]}%)"))


library(dplyr)
library(ggrepel)
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = shared_phy_out1)

#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)







M_B_ABX_out <- ggplot(df, aes(x = MDS1, y = MDS2, colour = meta_phy_out1$ID5, fill=meta_phy_out1$ID5)) +
  #scale_shape_manual(values = c(21, 21)) +
  #scale_shape_manual(values = c(16, 17)) +
  geom_point(mapping = aes(colour = meta_phy_out1$ID5), size=5, alpha=1.0) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  scale_color_manual(values=c("navyblue", "red3", "green", "magenta", "yellow","blue","black","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("navyblue", "red3", "green","magenta", "yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  #scale_color_manual(values=c("#B81A1A", "#283DA8", "#56D43A", "#5A5C59")) +
  #scale_fill_manual(values=c("pink", "dodgerblue", "#56D43A", "#5A5C59")) +
  #coord_fixed()+
  theme_classic() +
  labs(x=labx1, y=laby1) +
  #adds Weighted average score points for ASVs
  ##geom_point(data = WAscores_df, colour = "black",aes(x = MDS1, y = MDS2), inherit.aes = FALSE) +
  ##geom_text(data = WAscores_df, inherit.aes = FALSE, aes(x=MDS1, y= MDS2), colour="black", label = rownames(WAscores_df))+
  #ggtitle("Small intestine") +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=10, face=NULL), axis.text.y = element_text(color="black")) +
  #legend.background = element_rect(fill="white", color="black")
  #theme(legend.position = c(0.9, 0.4)) +
  #annotate("text", x=0.55, y=0.950, label = "Test: R2 = 0.055, p = 0.696", size=4) +
  #annotate("text", x=0.55, y=0.850, label = "Weight: R2 = 0.081, p = 0.005", size=4) +
  #theme(legend.title=element_blank())
  #theme(legend.position="botTypetom", legend.box = "horizontal") +
  #facet_wrap(~meta_ob1$ID3, scales = "free_x", ncol = 2) +
  theme(legend.title=element_blank())

#(file='M_B_ABX.pdf', plot = M_B_ABX_out, width=8, height=8) 






library(ggpubr)
female_plot <- ggarrange(F_A_Control_out, F_B_ABX_out, M_A_Control_out, M_B_ABX_out, labels = c("Female_Control", "Female_Abx", "Male_Control", "Male_Abx"),
                       common.legend = TRUE, legend = "bottom",
                       ncol = 2, nrow = 2, align="hv")









dist_bray1 <- vegdist(shared_phy_out1, method="bray", binary=TRUE)
#dist_bray1 <- vegdist(shared_phy_out1, method="bray", binary=FALSE)
adonis2(formula = dist_bray1 ~ Type + Sex, data = meta_phy_out1, permutations=1000)
adonis2(formula = dist_bray1 ~ Sample_ID, data = meta_phy_out1, permutations=1000, strata=meta_phy_out1$Mouse_ID)
adonis2(formula = dist_bray1 ~ Sex*Type, data = meta_phy_out1, permutations=1000)
adonis2(formula = dist_bray1 ~ Sex, data = meta_phy_out1, permutations=1000, strata=meta_phy_out1$Mouse_ID)


#adonis2(formula = dist_bray1 ~ Type, data = meta_phy_out1, permutations=10000)

#adonis2(formula = dist_bray1 ~ Type, data = meta_phy_out1, permutations=100000)


adonis2(formula = dist_bray1 ~ Mouse_ID, data = meta_phy_out1, permutations=1000)
adonis2(formula = dist_bray1 ~ Sex, data = meta_phy_out1, permutations=1000)
adonis2(formula = dist_bray1 ~ Type, data = meta_phy_out1, permutations=1000)
adonis2(formula = dist_bray1 ~ Type, data = meta_phy_out1, permutations=1000)
adonis2(formula = dist_bray1 ~ Sample_ID, data = meta_phy_out1, permutations=1000)

adonis2(formula = dist_bray1 ~ Sample_ID + Type, data = meta_phy_out1, permutations=1000, strata=meta_phy_out1$Mouse_ID)




adonis2(formula = dist_bray1 ~ Type, data = meta_phy_out1, permutations=100000)





#adonis2(formula = dist_bray1 ~ Sex, data = meta_phy_out1, permutations=1000)



#adonis2(formula = dist_bray1 ~ Sex*Type, data = meta_phy_out1, permutations=1000)




library(pairwiseAdonis)
pairwise.adonis2(dist_bray1 ~ID, data = meta_phy_out1, perm=10000)



pairwise.adonis2(dist_bray1 ~ ID5, data = meta_phy_out1, strata = 'Mouse_ID')






library(Maaslin2)



#setwd("/Users/zebra/Desktop/GF/")
setwd("/Users/zebra/Desktop/GF/Sex_differences")


meta_ob1  <- read.table("/Users/zebra/Desktop/GF/meta.12523.txt", header=T, row.names=1)
meta_ob1 <- meta_ob1[meta_ob1$Sample_ID!= "wk1" & meta_ob1$Sample_ID!= "wk2" & meta_ob1$Type2=="test",]
#meta_ob1 <- meta_ob1[meta_ob1$ID3=="F_A_Control",]
#meta_ob1 <- meta_ob1[meta_ob1$Type=="A_Control",]
#meta_ob1 <- meta_ob1[meta_ob1$Type!="B_ABX",]


#meta_ob1 <- meta_ob1[meta_ob1$Sample_ID!= "wk1" & meta_ob1$Sample_ID!= "wk2" & meta_ob1$Type2=="test",]
meta_ob1 <- meta_ob1[meta_ob1$Sample_ID!= "wk1" & meta_ob1$Sample_ID!= "wk2" & meta_ob1$Type2=="test",]
#meta_ob1 <- meta_ob1[meta_ob1$Sample_ID!= "wk1" & meta_ob1$Sample_ID!= "wk2" & meta_ob1$Type2=="test" & meta_ob1$Mouse_ID!="M_A8" & meta_ob1$Sex== "F" ,]

meta_ob1 <- meta_ob1[meta_ob1$Mouse_ID!="M_A8",]
meta_ob1 <- meta_ob1[meta_ob1$Mouse_ID!="M_A15",]
meta_ob1 <- meta_ob1[meta_ob1$Mouse_ID!="F_C10",]
meta_ob1 <- meta_ob1[meta_ob1$Mouse_ID!="F_A4",]


#meta_ob1 <- meta_ob1[meta_ob1$Type=="B_ABX" & meta_ob1$Sex=="F",]

tax_ob1 <- import_mothur(mothur_constaxonomy = ("/Users/zebra/Desktop/GF/taxa.txt"))




meta_ob1 <- meta_ob1[meta_ob1$Sample_ID!= "wk1" & meta_ob1$Sample_ID!= "wk2" & meta_ob1$Sex== "F" & meta_ob1$Type2=="test" & meta_ob1$Type=="B_ABX",]




nrow(meta_ob1)


#Get sample names
GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
#fullMothur<-read.csv("shared.12523.csv",)
fullMothur<-read.csv("F_ABX_unstrat.tsv.csv",)

AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtus))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1))

merged_phy_percent1 <- phyloseqobj.f

#merged_phy_percent1 <- phyloseq_standardize_otu_abundance(phyloseqobj.f, method = "total")


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


meta_phy_out1 <- pssd2veg(merged_phy_percent2)
shared_phy_out1 <- psotu2veg(merged_phy_percent2)



input_data  <- as.data.frame(t(shared_phy_out1))
input_metadata <- as.data.frame(meta_phy_out1)

row.names(input_data)





#colnames(input_data) == row.names(input_metadata)









fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "GF_ASV1_PWY", 
                    fixed_effects  = c("ASV1"),
                    #reference      = c("Type,CS_27", "Sex,F"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)






fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "GF_taxa_Rank2_Type_Sex_NEGBIN", 
                    fixed_effects  = c("Type", "Sex"),
                    reference      = c("Type,CS_27", "Sex,F"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)





fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "GF_taxa_Rank5_Type_Sex_NEGBIN", 
                    fixed_effects  = c("Type", "Sex"),
                    reference      = c("Type,CS_27", "Sex,F"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)





fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "GF_taxa_Type_Sex_NEGBIN", 
                    fixed_effects  = c("Type", "Sex"),
                    reference      = c("Type,CS_27", "Sex,F"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)




fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "GF_PWY_Type_Sex_NEGBIN", 
                    fixed_effects  = c("Type", "Sex"),
                    reference      = c("Type,CS_27", "Sex,F"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)




fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "CPP_PWY_58_ABX_samples_wk2_and_wk3_NEGBIN", 
                    fixed_effects  = c("Sex", "CPP"),
                    reference      = c("Sex,F"),
                    random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)











fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "CPP_PWY_31_pretest_Control_samples_NEGBIN", 
                    fixed_effects  = c("Sex", "CPP"),
                    reference      = c("Sex,F"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)















fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "CPP60_PWY_60_pretest_samples_NEGBIN", 
                    fixed_effects  = c("Type", "Sex", "CPP60"),
                    reference      = c("Type,A_Control", "Sex,F", "CPP60,No"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)




fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "CPP60_PWY_180_samples_NEGBIN", 
                    fixed_effects  = c("Type", "Sex", "Sample_ID", "CPP60", "Type"),
                    reference      = c("Type,A_Control", "Sex,F", "Sample_ID,C_wk3", "CPP60,No", "Type,A_SS"),
                    random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)





fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "CPP_60_pretest_samples", 
                    fixed_effects  = c("Type", "Sex", "CPP"),
                    reference      = c("Type,A_Control", "Sex,F"),
                    #random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)










fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "maaslin_F_ABX_SC_pre_post", 
                    fixed_effects  = c("ID5"),
                    reference      = c("ID5,C_wk3_B_SC"),
                    random_effects  = "Mouse_ID",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)








fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    min_prevalence=0.15,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "maaslin_F_ABX_Pathway_tax_out", 
                    fixed_effects  = c("Type", "CPP"),
                    reference      = c("Type,A_SS"),
                    #random_effects  = "Type",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)

fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    min_prevalence=0.15,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "maaslin_F_ABX_Pathway_Type_out", 
                    fixed_effects  = c("Type"),
                    reference      = c("Type,A_SS"),
                    #random_effects  = "Type",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    min_prevalence=0.15,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "maaslin_F_ABX_Pathway_CPP_out", 
                    fixed_effects  = c("CPP"),
                    #reference      = c("Type,A_SS"),
                    #random_effects  = "Type",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)



fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "maaslin_F_ABX_glom_tax_Type_out", 
                    fixed_effects  = c("Type"),
                    reference      = c("Type,A_SS"),
                    #random_effects  = "Type",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)

fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "maaslin_F_ABX_glom_tax_CPP_out", 
                    fixed_effects  = c("CPP", "Type"),
                    #reference      = c("Type,A_SS"),
                    #random_effects  = "Type",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)


fit_data = Maaslin2(input_data     = input_data, 
                    input_metadata = input_metadata,
                    max_significance=0.1,
                    #min_abundance=0.0,
                    min_prevalence=0.0,
                    #min_variance=20000,
                    normalization  = "NONE",
                    output         = "maaslin_Sex_Controls_tax_out", 
                    fixed_effects  = "Sex",
                    reference      = "Sex,A_SS,F",
                    #random_effects  = "Type",
                    analysis_method="NEGBIN",
                    transform = "NONE",
                    cores=12)







library(pROC)

library(vip)


meta_subset<- as.data.frame(pssd2veg(merged_phy_percent1))
eset <- as.data.frame(psotu2veg(merged_phy_percent1))





#make this example reproducible
set.seed(1)

#use 70% of dataset as training set and 30% as test set
sample <- sample(c(TRUE, FALSE), nrow(meta_subset), replace=TRUE, prob=c(0.7,0.3))
train  <- meta_subset[sample, ]
test   <- meta_subset[!sample, ]



#eset2 <- eset[, c("ASV6", "ASV2", "ASV267")]



##fit model with all data and feat importance
#uses values in eset as predictors and Group as predictor of PTL_PPROM with sample size set to the minimum number of samples belonging to a group (Control/PTL_PPROM)
#also assesses importance of predictors
modb <- randomForest::randomForest(x= data.frame(eset),y=factor(meta_subset$CPP60), sampsize=c(min(table(meta_subset$CPP60)),min(table(meta_subset$CPP60))),ntree = 1300,importance=T)
mod1 <- randomForest::randomForest(x= data.frame(eset),y=factor(meta_subset$CPP60), sampsize=c(min(table(meta_subset$CPP60)),min(table(meta_subset$CPP60))),ntree = 2000,importance=T)


#converts model to text
df<- vip::vi(mod1)



#takes the first 50 rows from df (i.e., 50 most important ASV features) and plots random forest
g <- df %>% dplyr::slice_head(n=50) %>% ggplot2::ggplot(aes(x=forcats::fct_reorder(Variable,Importance),y=Importance)) +
  ggplot2::geom_col(color="black", fill=rgb(red=51,green=51,blue=153,max=255))+
  ggplot2::coord_flip()+
  ggplot2::labs(x="") + 
  #ggplot2::ggtitle(paste(outcome, ifelse(names(outcomes[i])=="34", "<34", ""),summary_stat,GA_samples_max)) + 
  ggplot2::theme_light()

g









outcs <-NULL
predps <- NULL
perfs <- NULL

    #gets probabilities for test data based on random forest training model
    out=predict(modb,eset,type="prob")
    #gets the 90th percent quantile of "D" in the "C" subset of predictions
    cut=quantile(out[meta_subset$CPP60=="No","Yes"],0.9)
    #gets the "D" probabilities assoicated with PTL_PPROM
    pred=predict(modb,eset,type="prob")[,"Yes"]
    #if prediction in D is greater than 90% quantile 1 otherwise 0
    Yh=ifelse(pred>cut,1,0)
    #converts test group to numeric
    outcs=c(outcs,as.numeric(meta_subset$CPP60=="Yes"))
    #predictions added
    predps=c(predps,pred)
    #gets area under the curve from ROC plot
    perfs=c(perfs,pROC::auc(pROC::roc(response=as.numeric(meta_subset$CPP60=="Yes"),predictor=pred)))
  

  #stores ROC into RC
  #RC=pROC::roc(response=outcs,predictor=predps)
  plot(pROC::roc(response=outcs,predictor=predps))
  return(list(RC, aucs))















shared_phy_out1_percent <- shared_phy_out1 / rowSums(shared_phy_out1) * 100


df <- as.data.frame(shared_phy_out1_percent)
meta <- as.data.frame(meta_phy_out1)

nrow(df)
nrow(meta)




SI_abund <- as.data.frame(cbind(df$ASV206, df$ASV7, df$ASV250, df$ASV89, df$ASV6, df$ASV11, df$ASV49, df$ASV1, df$ASV137, df$ASV184))
colnames(SI_abund)[1] <- "Pseudomonadales"
colnames(SI_abund)[2] <- "Verrucomicrobiales"
colnames(SI_abund)[3] <- "Bacteroidia_NA"
colnames(SI_abund)[4] <- "Erysipelotrichales"
colnames(SI_abund)[5] <- "Bacteroidales"
colnames(SI_abund)[6] <- "Clostridiales"
colnames(SI_abund)[7] <- "Coriobacteriales"
colnames(SI_abund)[8] <- "Enterobacteriales"
colnames(SI_abund)[9] <- "Clostridia_NA"
colnames(SI_abund)[10] <- "Gastranaerophilales"
#colnames(SI_abund)[11] <- "Sex"






p1 <- ggplot(SI_abund ,aes(x=meta$Sex, y=SI_abund$Pseudomonadales, color=meta$Sex, fill=meta$Sex, shape=meta$Sex)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  ylim(0, 0.6)

p2 <- ggplot(SI_abund ,aes(x=meta$Sex, y=SI_abund$Verrucomicrobiales, color=meta$Sex, fill=meta$Sex, shape=meta$Sex)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  ylim(0, 15)

p3 <- ggplot(SI_abund ,aes(x=meta$Sex, y=SI_abund$Bacteroidia_NA, color=meta$Sex, fill=meta$Sex, shape=meta$Sex)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  ylim(0, 0.25)

p4 <- ggplot(SI_abund ,aes(x=meta$Sex, y=SI_abund$Erysipelotrichales, color=meta$Sex, fill=meta$Sex, shape=meta$Sex)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  ylim(0, 0.8)

p5 <- ggplot(SI_abund ,aes(x=meta$Sex, y=SI_abund$Bacteroidales, color=meta$Sex, fill=meta$Sex, shape=meta$Sex)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  ylim(0, 45)

p6 <- ggplot(SI_abund ,aes(x=meta$Sex, y=SI_abund$Clostridiales, color=meta$Sex, fill=meta$Sex, shape=meta$Sex)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  ylim(0, 100)

p7 <- ggplot(SI_abund ,aes(x=meta$Sex, y=SI_abund$Coriobacteriales, color=meta$Sex, fill=meta$Sex, shape=meta$Sex)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  ylim(0, 2.0)

p8 <- ggplot(SI_abund ,aes(x=meta$Sex, y=SI_abund$Enterobacteriales, color=meta$Sex, fill=meta$Sex, shape=meta$Sex)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  ylim(0, 0.125)

p9 <- ggplot(SI_abund ,aes(x=meta$Sex, y=SI_abund$Clostridia_NA, color=meta$Sex, fill=meta$Sex, shape=meta$Sex)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  ylim(0, 0.85)

p10 <- ggplot(SI_abund ,aes(x=meta$Sex, y=SI_abund$Gastranaerophilales, color=meta$Sex, fill=meta$Sex, shape=meta$Sex)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  ylim(0, 1.0)




library(ggpubr)
plot_merged_1 <- ggarrange(p5, p2, p4, p1, p3,
                           labels = c("Bacteroidales", "Verrucomicrobiales", "Erysipelotrichales", "Pseudomonadales", "Bacteroidia_NA"),
                           common.legend =TRUE, legend = "bottom",
                           ncol = 1, nrow = 5, align="hv")

plot_merged_2 <- ggarrange(p6, p7, p8, p10, p9,
                           labels = c("Clostridiales", "Coriobacteriales", "Enterobacteriales", "Gastranaerophilales", "Clostridia_NA"),
                           common.legend =TRUE, legend = "bottom",
                           ncol = 1, nrow = 5, align="hv")

plot_merged_3 <- ggarrange(plot_merged_1, plot_merged_2,
                           common.legend =TRUE, legend = "none",
                           ncol = 2, nrow = 1, align="hv")






ggplot(SI_abund ,aes(x=Sex, y=Pseudomonadales, color=Sex, fill=Sex)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75), size=4, alpha=0.5, stroke=1.25) +
  #scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("navyblue", "red3", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("navyblue", "red3", "red", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  theme_classic() +
  #theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  ylim(-1, 3.5)
#theme(axis.title.y=element_blank())



















p1 <- ggplot(SI_abund ,aes(x=Sex, y=Pseudomonadales, color=Sex, fill=Sex)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75), size=4, alpha=0.5, stroke=1.25) +
  #scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("navyblue", "red3", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("navyblue", "red3", "red", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  ylim(0, 3.5)





































#Bacteroides_thetaiotaomicron_G5
#Bacteroides_uniformis_G10
#Duncaniella_sp_B8_G8
#Parabacteroides_distasonis_G25
#Duncaniella_dubosii_G7
#Bifidobacterium_pseudolongum_G11




SI_abund <- as.data.frame(cbind(df$CATECHOL_ORTHO_CLEAVAGE_PWY, df$PWY_5417, df$PWY_6182, df$PWY_5431, df$PWY_6185, df$PWY_6629, meta$CPP, meta$Type))
colnames(SI_abund)[1] <- "CATECHOL.ORTHO.CLEAVAGE.PWY"
colnames(SI_abund)[2] <- "PWY.5417"
colnames(SI_abund)[3] <- "PWY.6182"
colnames(SI_abund)[4] <- "PWY.5431"
colnames(SI_abund)[5] <- "PWY.6185"
colnames(SI_abund)[6] <- "PWY.6629"
colnames(SI_abund)[7] <- "CPP"
colnames(SI_abund)[8] <- "Type"



#write.csv(SI_abund ,file="F_ABX_plot_data.csv")


SI_abund <- read.table("/Users/zebra/Desktop/GF/F_ABX_plot_data.txt", header=T, row.names=1)



ggplot(SI_abund, aes(x = Type, y = CATECHOL.ORTHO.CLEAVAGE.PWY)) + 
  geom_point() + 
  scale_y_continuous() + 
  geom_jitter()




p1 <- ggplot(SI_abund ,aes(x=Type, y=log10(CATECHOL.ORTHO.CLEAVAGE.PWY), color=Type, fill=Type)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75), size=4, alpha=0.5, stroke=1.25) +
  #scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("navyblue", "red3", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("navyblue", "red3", "red", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  ylim(0, 3.5)
  #theme(axis.title.y=element_blank())

  p3 <- ggplot(SI_abund ,aes(x=Type, y=log10(PWY.5417), color=Type, fill=Type)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75), size=4, alpha=0.5, stroke=1.25) +
  #scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("navyblue", "red3", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("navyblue", "red3", "red", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
    ylim(0, 3.5)
  #theme(axis.title.y=element_blank())
  
  
  p5 <- ggplot(SI_abund ,aes(x=Type, y=log10(PWY.6182), color=Type, fill=Type)) +
    geom_boxplot(outlier.shape = NA, alpha=0) +
    geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75), size=4, alpha=0.5, stroke=1.25) +
    #scale_shape_manual(values = c(21, 21, 1, 1)) +
    scale_color_manual(values=c("navyblue", "red3", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
    scale_fill_manual(values=c("navyblue", "red3", "red", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
    theme_classic() +
    theme(axis.text.x=element_blank()) +
    theme(axis.text.y=element_text(size=10)) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank()) +
    ylim(0, 3.5)
  #theme(axis.title.y=element_blank())
  
  p7 <- ggplot(SI_abund ,aes(x=Type, y=log10(PWY.5431), color=Type, fill=Type)) +
    geom_boxplot(outlier.shape = NA, alpha=0) +
    geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75), size=4, alpha=0.5, stroke=1.25) +
    #scale_shape_manual(values = c(21, 21, 1, 1)) +
    scale_color_manual(values=c("navyblue", "red3", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
    scale_fill_manual(values=c("navyblue", "red3", "red", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
    theme_classic() +
    theme(axis.text.x=element_blank()) +
    theme(axis.text.y=element_text(size=10)) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank()) +
    ylim(0, 3.5)
  #theme(axis.title.y=element_blank())
  
  p9 <- ggplot(SI_abund ,aes(x=Type, y=log10(PWY.6185), color=Type, fill=Type)) +
    geom_boxplot(outlier.shape = NA, alpha=0) +
    geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75), size=4, alpha=0.5, stroke=1.25) +
    #scale_shape_manual(values = c(21, 21, 1, 1)) +
    scale_color_manual(values=c("navyblue", "red3", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
    scale_fill_manual(values=c("navyblue", "red3", "red", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
    theme_classic() +
    theme(axis.text.x=element_blank()) +
    theme(axis.text.y=element_text(size=10)) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank()) +
    ylim(0, 3.5)
  #theme(axis.title.y=element_blank())
  
  
  p11 <- ggplot(SI_abund ,aes(x=Type, y=log10(PWY.6629+000.1), color=Type, fill=Type)) +
    geom_boxplot(outlier.shape = NA, alpha=0) +
    geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75), size=4, alpha=0.5, stroke=1.25) +
    #scale_shape_manual(values = c(21, 21, 1, 1)) +
    scale_color_manual(values=c("navyblue", "red3", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
    scale_fill_manual(values=c("navyblue", "red3", "red", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
    theme_classic() +
    #theme(axis.text.x=element_blank()) +
    theme(axis.text.y=element_text(size=10)) +
    theme(axis.ticks.x=element_blank(),
          axis.title.x=element_blank()) +
    ylim(-1, 3.5)
  #theme(axis.title.y=element_blank())
  
  
  
  library(ggpmisc)

p2 <- ggplot(SI_abund, aes(CPP, log10(CATECHOL.ORTHO.CLEAVAGE.PWY))) +
  stat_poly_line() +
  stat_poly_eq() +
  geom_point() +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
theme(axis.title.y=element_blank())



p4 <- ggplot(SI_abund, aes(CPP, log10(PWY.5417))) +
  stat_poly_line() +
  stat_poly_eq() +
  geom_point() +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())

p6 <- ggplot(SI_abund, aes(CPP, log10(PWY.6182))) +
  stat_poly_line() +
  stat_poly_eq() +
  geom_point() +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())
  
p8 <- ggplot(SI_abund, aes(CPP, log10(PWY.5431))) +
  stat_poly_line() +
  stat_poly_eq() +
  geom_point() +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())

p10 <- ggplot(SI_abund, aes(CPP, log10(PWY.6185))) +
  stat_poly_line() +
  stat_poly_eq() +
  geom_point() +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())

p12 <- ggplot(SI_abund, aes(CPP, log10(PWY.6629+0.0001))) +
  stat_poly_line() +
  stat_poly_eq() +
  geom_point() +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())







library(ggpubr)
plot_merged_1 <- ggarrange(p3, p7, p9, p11,
                         common.legend =TRUE, legend = "right",
                         ncol = 1, nrow = 4, align="hv")

plot_merged_2 <- ggarrange(p4, p8, p10, p12,
                         common.legend =TRUE, legend = "right",
                         ncol = 1, nrow = 4, align="hv")

plot_merged_3 <- ggarrange(plot_merged_1, plot_merged_2,
                           common.legend=FALSE, legend = "right",
                           ncol = 2, nrow = 1, align="hv")







dev.off()


merged_phy_percent2 <- tax_glom(merged_phy_percent1, taxrank="Rank2")

#merged_phy_percent2 <- phyloseq_standardize_otu_abundance(merged_phy_percent , method = "total")


SI_pruned <- prune_taxa(names(sort(taxa_sums(merged_phy_percent1),TRUE)[1:200]), merged_phy_percent1)
sampleOrder = unique(sample_names(SI_pruned))
taxaOrder = rev(unique(taxa_names(SI_pruned)))


bar_plot <- plot_bar(merged_phy_percent2, fill="Rank2") + facet_wrap(~Type, scales = "free_x", ncol = 2) +
  theme(legend.position = "bottom")
  #scale_x_discrete(labels = bar_plot$data$Mouse_ID)


  
  
  




plot_heatmap(SI_pruned, "PCoA", "bray", sample.order = "Sex", taxa.label = "Rank1", low="white", high="dodgerblue3", na.value="azure1") +
  facet_grid(~Type, scales="free_x", switch="both") +
  #theme(axis.text.x=element_blank()) +
  #theme(axis.text.x=element_text(size=4)) +
  theme(axis.text.y=element_text(size=8)) +
  #theme(axis.text.x=element_blank(),
  #axis.ticks.x=element_blank(),
  #axis.title.x=element_blank()) +

#theme(strip.text.x = element_blank()) +
#theme(axis.title.y=element_blank())












plot_heatmap(SI_pruned , trans = NULL, method="PCOA", low="white", high="dodgerblue3", na.value="azure1", title=NULL) +
  facet_wrap(~Type, scales = "free_x", ncol = 2)
theme(axis.text.x=element_blank()) +
  #theme(axis.text.x=element_text(size=4)) +
  theme(axis.text.y=element_text(size=8)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  facet_wrap(~Sample_ID, scales = "free_x", nrow = 5)
#theme(strip.text.x = element_blank()) +
#theme(axis.title.y=element_blank())

  
  
  


bar_plot <- plot_bar(SI_pruned, fill="Rank6") + facet_wrap(~Sample_ID, scales = "free_x", nrow = 5) 


phyla_plot_data <- ggplot_build(bar_plot)$plot$data

dev.off()

#write.csv(phyla_plot_data ,file="phyla_plot_data.csv")




SI_pruned <- prune_taxa(names(sort(taxa_sums(merged_phy_percent1),TRUE)[1:5]), merged_phy_percent1)
sampleOrder = unique(sample_names(SI_pruned))
taxaOrder = rev(unique(taxa_names(SI_pruned)))







  
  
  plot_heatmap(SI_pruned, "PCoA", "bray", taxa.label = "Rank1", low="white", high="dodgerblue3", na.value="azure1") +
    facet_grid(~Type, scales="free_x", switch="both") +
    #theme(axis.text.x=element_blank()) +
    #theme(axis.text.x=element_text(size=4)) +
    theme(axis.text.y=element_text(size=8)) +
    #theme(axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    #axis.title.x=element_blank()) +
    facet_wrap(~Sample_ID, scales = "free_x", nrow = 5)
  #theme(strip.text.x = element_blank()) +
  #theme(axis.title.y=element_blank())
  
  dev.off()


  
  
library(ggpubr)
plot_merged <- ggarrange(bar_plot,heatmap,
                       common.legend = FALSE, legend = "right",
                       ncol = 1, nrow = 2, align="hv")




library(dendextend)
library(circlize)
library(magrittr)

distance = distance(SI_pruned, method="bray", binary=FALSE)
#distance = dist(phyloseqobj.f, method ="bray")    

hcluster = hclust(distance, method ="ward.D2")
dend <- as.dendrogram(hcluster)

#plot(dend)


cols <- c("#e6194B", "#f58231", "#ffe119", "#3cb44b", "#000075", "#9A6324", "#911eb4")

dend <- color_branches(dend, k = 40, col = cols)
dend %<>% set("labels_col", value = cols, k= 40)
dend %<>% set("labels_cex", 1)
dend %<>% set("branches_lwd", 1)

circlize_dendrogram(dend)

dev.off()



library(RFLPtools)
write.hclust(hcluster, "clusters.txt", "hello", k = 40, append = FALSE, dec = ",")














meta_phy_out1$wk3p


merged_phy_percent2 <- tax_glom(merged_phy_percent1, taxrank="Rank2")



SI_pruned <- prune_taxa(names(sort(taxa_sums(merged_phy_percent1),TRUE)[1:35]), merged_phy_percent1)
sampleOrder = unique(sample_names(SI_pruned))
taxaOrder = rev(unique(taxa_names(SI_pruned)))


bar_plot <- plot_bar(merged_phy_percent2 , fill="Rank2") + facet_wrap(~ID3, scales = "free_x", nrow = 1) 


plot_heatmap(SI_pruned , trans = NULL, sample.order = sampleOrder, taxa.order = taxaOrder, low="white", high="dodgerblue3", na.value="azure1",
             title=NULL) +
  theme(axis.text.x=element_blank()) +
  #theme(axis.text.x=element_text(size=4)) +
  theme(axis.text.y=element_text(size=8)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  facet_wrap(~ID4, scales = "free_x", nrow = 1)










plot_data <- ggplot_build(bar_plot)$plot$data

write.csv(plot_data ,file="plot_data2.csv")


SI_HM <- plot_heatmap(SI_pruned , trans = NULL, sample.order = sampleOrder, taxa.order = taxaOrder, low="white", high="dodgerblue3", na.value="azure1",
                      title=NULL) +
  theme(axis.text.x=element_blank()) +
  #theme(axis.text.x=element_text(size=4)) +
  theme(axis.text.y=element_text(size=8)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  facet_wrap(~wk3p, scales = "free_x", nrow = 1)
  #theme(strip.text.x = element_blank()) +
  #theme(axis.title.y=element_blank())


  plot_heatmap(SI_pruned, sample.order = sampleOrder, taxa.order = taxaOrder)

  
  heatmap <- plot_heatmap(SI_pruned, method = "PCoA", distance = "bray", "Mouse_ID", "Rank8", trans = NULL, sample.order="Mouse_ID", low="azure1", high="dodgerblue3", na.value="white") +
  facet_wrap(~wk3p, scales = "free_x", nrow = 1)
  
  
SI_bar <- plot_bar(SI_pruned , fill="Rank4") + facet_wrap(~ID3, scales = "free_x", nrow = 1) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
  theme(legend.position="right")
  theme(strip.text.x = element_blank()) +
  theme(legend.position="right")






library(ggpubr)
SI_merged <- ggarrange(SI_plot, SI_HM, SI_bar_BD,
                              common.legend = FALSE, legend = "right",
                              ncol = 3, nrow = 1, align="hv")















###Top 10 taxa



merged_phy_percent1 <- phyloseq_standardize_otu_abundance(phyloseqobj.f, method = "total")

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


meta_phy_out1 <- pssd2veg(merged_phy_percent1)
shared_phy_out1 <- psotu2veg(merged_phy_percent1)

df <- as.data.frame(shared_phy_out1)
meta <- as.data.frame(meta_phy_out1)

nrow(df)
nrow(meta)

colnames(df)






phyla_abund_CPP <- as.data.frame(cbind(df$ASV2, df$ASV3, df$ASV6, df$ASV7, df$ASV41, df$ASV49, df$ASV171, df$ASV184, df$ASV284, df$ASV332, df$ASV573, df$ASV1157, meta$CPP, meta$ID4, meta$Sex, meta$Type, meta$Type, meta$Mouse_ID))
colnames(phyla_abund_CPP)[1] <- "ASV2"
colnames(phyla_abund_CPP)[2] <- "ASV3"
colnames(phyla_abund_CPP)[3] <- "ASV6"
colnames(phyla_abund_CPP)[4] <- "ASV7"
colnames(phyla_abund_CPP)[5] <- "ASV41"
colnames(phyla_abund_CPP)[6] <- "ASV49"
colnames(phyla_abund_CPP)[7] <- "ASV171"
colnames(phyla_abund_CPP)[8] <- "ASV184"
colnames(phyla_abund_CPP)[9] <- "ASV284"
colnames(phyla_abund_CPP)[10] <- "ASV332"
colnames(phyla_abund_CPP)[11] <- "ASV573"
colnames(phyla_abund_CPP)[12] <- "ASV1157"
colnames(phyla_abund_CPP)[13] <- "CPP"
colnames(phyla_abund_CPP)[14] <- "ID4"
colnames(phyla_abund_CPP)[15] <- "Sex"
colnames(phyla_abund_CPP)[16] <- "Type"
colnames(phyla_abund_CPP)[17] <- "Type"
colnames(phyla_abund_CPP)[18] <- "Mouse_ID"








#write.csv(phyla_abund_CPP,file="phyla_abund_CPP.csv")

phyla_abund_CPP  <- read.table("/Users/zebra/Desktop/GF/phyla_abund_CPP.txt", header=T, row.names=1)

#phyla_abund_CPP <- phyla_abund_CPP[phyla_abund_CPP$Sex== "F" & phyla_abund_CPP$Type2=="test" & phyla_abund_CPP$Type=="B_ABX",]

phyla_abund_CPP <- phyla_abund_CPP[phyla_abund_CPP$Sex== "F",]

meta_ob1 <- meta_ob1[meta_ob1$Mouse_ID!="M_A8",]
meta_ob1 <- meta_ob1[meta_ob1$Mouse_ID!="M_A15",]
meta_ob1 <- meta_ob1[meta_ob1$Mouse_ID!="F_C10",]
meta_ob1 <- meta_ob1[meta_ob1$Mouse_ID!="F_A4",]




meta_phy_out1 <- as.data.frame(pssd2veg(merged_phy_percent1))
shared_phy_out1 <- as.data.frame(psotu2veg(merged_phy_percent1))


library(ggplot2)
library(ggpmisc)











p1 <- ggplot(shared_phy_out1, aes(x=meta_phy_out1$CPP, y=log10(shared_phy_out1$PWY120+1))) +
  geom_point(aes(colour = meta_phy_out1$Sex)) +
  geom_smooth(method = "lm") +
  stat_poly_line() +
  stat_poly_eq() +
  theme_bw() +
  theme(axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  ylab("Protocatechuate degradation II_PWY120")

p2 <- ggplot(shared_phy_out1, aes(x=meta_phy_out1$CPP, y=log10(shared_phy_out1$PWY139+1))) +
  geom_point(aes(colour = meta_phy_out1$Sex)) +
  geom_smooth(method = "lm") +
  stat_poly_line() +
  stat_poly_eq() +
  theme_bw() +
  theme(axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  ylab("L-histidine degradation II_PWY139")

p3 <- ggplot(shared_phy_out1, aes(x=meta_phy_out1$CPP, y=log10(shared_phy_out1$PWY162+1))) +
  geom_point(aes(colour = meta_phy_out1$Sex)) +
  geom_smooth(method = "lm") +
  stat_poly_line() +
  stat_poly_eq() +
  theme_bw() +
  theme(axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  ylab("catechol degradation III_PWY162")

p4 <- ggplot(shared_phy_out1, aes(x=meta_phy_out1$CPP, y=log10(shared_phy_out1$PWY225+1))) +
  geom_point(aes(colour = meta_phy_out1$Sex)) +
  geom_smooth(method = "lm") +
  stat_poly_line() +
  stat_poly_eq() +
  theme_bw() +
  theme(axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  ylab("4-methylcatechol degradation_PWY225")

p5 <- ggplot(shared_phy_out1, aes(x=meta_phy_out1$CPP, y=log10(shared_phy_out1$PWY18+1))) +
  geom_point(aes(colour = meta_phy_out1$Sex)) +
  geom_smooth(method = "lm") +
  stat_poly_line() +
  stat_poly_eq() +
  theme_bw() +
  theme(axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  ylab("Catechol degradation to beta ketoadipate_PWY18")

p6 <- ggplot(shared_phy_out1, aes(x=meta_phy_out1$CPP, y=log10(shared_phy_out1$PWY382+1))) +
  geom_point(aes(colour = meta_phy_out1$Sex)) +
  geom_smooth(method = "lm") +
  stat_poly_line() +
  stat_poly_eq() +
  theme_bw() +
  theme(axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  ylab("Teichoic acid (poly-glycerol) biosynthesis_PWY382")


  




library(ggpubr)
CCP_plot <- ggarrange(p1, p2, p3, p4, p5, p6,
                        common.legend = TRUE, legend = "bottom",
                        ncol = 3, nrow = 2, align="hv")










ggplot(shared_phy_out1, aes(x=meta_phy_out1$CPP, y=shared_phy_out1$PWY139, colour = factor(meta_phy_out1$Sex))) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  stat_poly_line() +
  stat_poly_eq() +



library(ggpmisc)

p2 <- ggplot(phyla_abund_CPP, aes(CPP, ASV2)) +
  stat_poly_line() +
  stat_poly_eq() +
  geom_point() +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  facet_wrap(~ID4, scales = "free_y", nrow = 2)



















SI_abund <- as.data.frame(cbind(df$ASV5, df$ASV30, df$ASV171, df$ASV34, df$ASV68, meta))
colnames(SI_abund)[1] <- "ASV5"
colnames(SI_abund)[2] <- "ASV30"
colnames(SI_abund)[3] <- "ASV171"
colnames(SI_abund)[4] <- "ASV34"
colnames(SI_abund)[5] <- "ASV68"
colnames(SI_abund)[6] <- "G11"





G1_p <- ggplot(SI_abund ,aes(x=CPP, y=df$ASV34, color=ID4, fill=ID4, shape=ID4)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  #scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("navyblue", "red3", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  scale_fill_manual(values=c("navyblue", "red3", "red", "red", "green","yellow","blue","black","magenta","pink","cyan","gray","orange","brown","purple")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


ASV5_plot <- ggplot(SI_abund , aes(CPP, ASV68)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ID4, scales = "free_x", nrow = 1)



ASV30_plot <- ggplot(SI_abund , aes(CPP, ASV30)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ID4, scales = "free_x", nrow = 1)











library(tidyverse)
library(ggrepel)


# A short function for outputting the tables
knitr_table <- function(x) {
  x %>% 
    knitr::kable(format = "html", digits = Inf, 
                 format.args = list(big.mark = ",")) %>%
    kableExtra::kable_styling(font_size = 15)
}


# Import data
data <- read_tsv("https://raw.githubusercontent.com/sdgamboa/misc_datasets/master/L0_vs_L20.tsv")
dim(data)

head(data) %>% 
  knitr_table()

p1 <- ggplot(data, aes(log10(coef+0.0000001), FDR)) + # -log10 conversion  
  geom_point(size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))
p1




data <- data %>% 
  mutate(
    Expression = case_when(logFC >= log(2) & FDR <= 0.05 ~ "Up-regulated",
                           logFC <= -log(2) & FDR <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
head(data)

p2 <- ggplot(data, aes(logFC, -log(FDR,10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 
p2


data <- data %>% 
  mutate(
    Significance = case_when(
      abs(logFC) >= log(2) & FDR <= 0.05 & FDR > 0.01 ~ "FDR 0.05", 
      abs(logFC) >= log(2) & FDR <= 0.01 & FDR > 0.001 ~ "FDR 0.01",
      abs(logFC) >= log(2) & FDR <= 0.001 ~ "FDR 0.001", 
      TRUE ~ "Unchanged")
  )
head(data) %>% 
  knitr_table()


p3 <- ggplot(data, aes(logFC, -log(FDR,10))) +
  geom_point(aes(color = Significance), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

p3



top <- 10
top_genes <- bind_rows(
  data %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(FDR, desc(abs(logFC))) %>% 
    head(top),
  data %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(FDR, desc(abs(logFC))) %>% 
    head(top)
)



p3 <-  p3 +
  geom_label_repel(data = top_genes,
                   mapping = aes(logFC, -log(FDR,10), label = Genes),
                   size = 2)
p3















data_dif <- read.table("/Users/zebra/Desktop/GF/maaslin_CPP_Both_tax_out/all_results.txt", header=T, row.names=1)








p1 <- EnhancedVolcano(data_dif,
                      lab = rownames(data_dif),
                      x = "coef",
                      y = "qval",
                      #pCutoff = 0.1,
                      max.overlaps = 40,
                      FCcutoff = 0,
                      #ylim = c(0, -log10(10e-12)),
                      #pointSize = c(ifelse(data_dif$coef>2, 8, 1)),
                      pointSize = 3.5,
                      labSize = 6.0,
                      #shape = c(6, 6, 19, 16),
                      title = "Saline vs Cocaine",
                      subtitle = "Differential abundance",
                      #caption = bquote(~Log[2]~ "fold change cutoff, 2; p-value cutoff, 10e-4"),
                      legendPosition = "right",
                      legendLabSize = 14,
                      col = c("grey30", "forestgreen", "royalblue", "red2"),
                      colAlpha = 0.9,
                      drawConnectors = TRUE,
                      #hline = c(10e-8),
                      widthConnectors = 0.5)

p1












p1 <- EnhancedVolcano(data_dif,
                      lab = rownames(data_dif),
                      x = "coef",
                      y = "qval",
                      pCutoff = 0.1,
                      #FCcutoff = 2,
                      #ylim = c(0, -log10(10e-12)),
                      pointSize = c(ifelse(data_dif$coef>2, 8, 1)),
                      labSize = 3.0,
                      #shape = c(6, 6, 19, 16),
                      title = "Saline vs Cocaine in Control Mice",
                      subtitle = "Differential abundance",
                      #caption = bquote(~Log[2]~ "fold change cutoff, 2; p-value cutoff, 10e-4"),
                      legendPosition = "right",
                      legendLabSize = 14,
                      colAlpha = 0.9,
                      #colGradient = c('red3', 'royalblue'),
                      drawConnectors = TRUE,
                      hline = c(10e-8),
                      widthConnectors = 0.5)

p1


  

library(EnhancedVolcano)

EnhancedVolcano(data_dif,
                lab = rownames(data_dif),
                x = 'coef',
                y = 'pval',
pCutoff = 0.1,
#FCcutoff = 25,
pointSize = 3.0,
labSize = 6.0)
  
  

EnhancedVolcano(data_dif,
                lab = rownames(data_dif),
                x = 'coef',
                y = 'pval',
                title = 'N061011 versus N61311',
                #pCutoff = 0.1,
                #FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)
  
























  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 







ggplot(data=data_dif, aes(x=coeff, y=-log10(pval))) +
  geom_point( size=1 ,aes(color=as.factor(thershold))) +
  xlim(c(-10, 10)) + ylim(c(0, 6)) +
  xlab("log2 fold change") + 
  ylab("-log10 p-value")  + 
  annotate("label", x =c(-8,5), y = 4.75, label = c("400","120"), col=c("red","steelblue"))+
  annotate("text", x =c(-8,5), y = 5, label = c("50 FC>4","8FC <-4"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")

#SI_10 <- prune_taxa(names(sort(taxa_sums(merged_phy_percent1),TRUE)[1:10]), merged_phy_percent1)





SI_10 <- shared_phy_out1

#merged_phy_counts <- merge_phyloseq(shared_phy_out1)

#SI_10 <- phyloseq_standardize_otu_abundance(merged_phy_counts, method = "total")



# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(SI_10) {
  sd <- sample_data(SI_10)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(SI_10) {
  OTU <- otu_table(SI_10)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}







SI_10_meta <- pssd2veg(SI_10)
SI_10_shared <- psotu2veg(SI_10)

df <- as.data.frame(SI_10_shared)
meta <- as.data.frame(SI_10_meta)

nrow(df)
nrow(meta)




#Bacteroides_thetaiotaomicron_G5
#Bacteroides_uniformis_G10
#Duncaniella_sp_B8_G8
#Parabacteroides_distasonis_G25
#Duncaniella_dubosii_G7
#Bifidobacterium_pseudolongum_G11




SI_abund <- as.data.frame(cbind(df$G5, df$G10, df$G8, df$G25, df$G7, df$G11, meta))
colnames(SI_abund)[1] <- "G5"
colnames(SI_abund)[2] <- "G10"
colnames(SI_abund)[3] <- "G8"
colnames(SI_abund)[4] <- "G25"
colnames(SI_abund)[5] <- "G7"
colnames(SI_abund)[6] <- "G11"



write.csv(SI_abund,file="plot_data.csv")




G1_p <- ggplot(SI_abund ,aes(x=Order, y=df$G1, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G2_p <- ggplot(SI_abund ,aes(x=Order, y=df$G2, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())



G3_p <- ggplot(SI_abund ,aes(x=Order, y=df$G3, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G4_p <- ggplot(SI_abund ,aes(x=Order, y=df$G4, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())

G5_p <- ggplot(SI_abund ,aes(x=Order, y=df$G5, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G6_p <- ggplot(SI_abund ,aes(x=Order, y=df$G6, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G7_p <- ggplot(SI_abund ,aes(x=Order, y=df$G7, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=5, alpha=0.5) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G8_p <- ggplot(SI_abund ,aes(x=Order, y=df$G8, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G11_p <- ggplot(SI_abund ,aes(x=Order, y=df$G11, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G33_p <- ggplot(SI_abund ,aes(x=Order, y=df$G33, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())




library(ggpubr)
SI_barplot <- ggarrange(G1_p, G2_p, G3_p, G4_p, G5_p, G6_p, G7_p, G8_p, G11_p, G33_p, labels = c("Muribaculum_G1", 	"Akkermansia_G2", 	"Lactobacillus_G3", 	"Muribaculum_G4", 	"Bacteroides_G5", 	"Lactobacillus_G6", 	"Duncaniella_G7", 	"Duncaniella_G8", 	"Bifidobacterium_G11", "Candidatus_Arthromitus_G33"),
          common.legend = TRUE, legend = "right",
          ncol = 2, nrow = 5, align="hv")


SI_plot_data <- G33_p$data

#write.csv(SI_plot_data,file="SI_plot_data.csv")








SI_weight_df <- data.frame(SI_plot_data$G2, SI_plot_data$G7, SI_plot_data$G8, SI_plot_data$Weight, SI_plot_data$Test, SI_plot_data$Order2)

colnames(SI_weight_df)[1] <- "G2"
colnames(SI_weight_df)[2] <- "G7"
colnames(SI_weight_df)[3] <- "G8"
colnames(SI_weight_df)[4] <- "Weight"
colnames(SI_weight_df)[5] <- "Test"
colnames(SI_weight_df)[6] <- "Order2"




G2_plot <- ggplot(SI_weight_df, aes(Weight, G2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~Order2, scales = "free_x", nrow = 1)



G7_plot <- ggplot(SI_weight_df, aes(Weight, G7)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~Order2, scales = "free_x", nrow = 1)


G8_plot <- ggplot(SI_weight_df, aes(Weight, G8)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~Order2, scales = "free_x", nrow = 1)






SI_Weight <- ggarrange(G2_plot, G7_plot, G8_plot, labels = c("Akkermansia_G2", "Duncaniella_G7", 	"Duncaniella_G8"),
                       common.legend = TRUE, legend = "right",
                       ncol = 1, nrow = 3, align="hv")


















###Cecum

setwd("/Users/zebra/Desktop/GF/")

meta_ob1  <- read.table("/Users/zebra/Desktop/GF/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("/Users/zebra/Desktop/GF/taxa.txt"))


meta_ob1 <-  meta_ob1[meta_ob1$Origin=='Neonatal' & meta_ob1$Order=="B_Ce", ]
#meta_ob1 <- meta_ob1[meta_ob1$Origin=='Neonatal' & meta_ob1$Order=="B_Ce" & meta_ob1$Order2!="C" & meta_ob1$Order2!="D",]
#meta_ob1 <- meta_ob1[meta_ob1$Origin=='Neonatal' & meta_ob1$Order=="B_Ce" & meta_ob1$Order2!="A" & meta_ob1$Order2!="B",]


nrow(meta_ob1)


write.csv(meta_ob1,file="seqs_to_upload.csv")


#any(duplicated(colnames(meta_ob1)))


#Get sample names
GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtu))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))


#Rarefy
#minimum <- min(sample_sums(phyloseqobj.f))
out2 <- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1278327, rngseed = 1, replace = FALSE)

#Count to percent
phy1 <- phyloseq_standardize_otu_abundance(out2, method = "total")


merged_phy_percent2 <- phyloseq_standardize_otu_abundance(phy1, method = "total")


#Ce_phy <- merged_phy_percent2 


# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent2) {
  sd <- sample_data(merged_phy_percent2)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent2) {
  OTU <- otu_table(merged_phy_percent2)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


meta_phy_out2 <- pssd2veg(merged_phy_percent2)
shared_phy_out2 <- psotu2veg(merged_phy_percent2)



dist_bray2 <- capscale(shared_phy_out2~1, distance="bray", binary=FALSE, eig=TRUE)

prop_explained <- summary(eigenvals(dist_bray2))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)


data.scores = scores(dist_bray2, display = "sites", choices = c(1:10))
#summary(data.scores)


df <- data.frame(data.scores)

library(glue)
labx2 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby2 <- c(glue("PCo 2 ({prop_explained [2]}%)"))


library(dplyr)
library(ggrepel)
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = shared_phy_out2)

#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)






Ce_AB_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = meta_phy_out2$Order2, fill=meta_phy_out2$Order2, shape=meta_phy_out2$Order2)) +
  #scale_shape_manual(values = c(22, 22)) +
  scale_shape_manual(values = c(0, 0)) +
  geom_point(mapping = aes(colour = meta_phy_out2$Order2), size=5, alpha=1.0) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  #scale_color_manual(values=c("#B81A1A", "#283DA8", "#56D43A", "#5A5C59")) +
  #scale_fill_manual(values=c("pink", "dodgerblue", "#56D43A", "#5A5C59")) +
  #coord_fixed()+
  theme_classic() +
  labs(x=labx2, y=laby2) +
  #adds Weighted average score points for ASVs
  ##geom_point(data = WAscores_df, colour = "black",aes(x = MDS1, y = MDS2), inherit.aes = FALSE) +
  ##geom_text(data = WAscores_df, inherit.aes = FALSE, aes(x=MDS1, y= MDS2), colour="black", label = rownames(WAscores_df))+
  #ggtitle("Cecum") +
  #theme(legend.position="bottom", legend.box = "vertical") +
  #theme(axis.text.x = element_text(color="black"), text = element_text(size=10, face=NULL), axis.text.y = element_text(color="black")) +
  #legend.background = element_rect(fill="white", color="black")
  #theme(legend.position = c(0.9, 0.4)) + 
  #annotate("text", x=0.55, y=0.950, label = "Test: R2 = 0.155, p = 0.0002", size=4) +
  #annotate("text", x=0.55, y=0.850, label = "Weight: R2 = 0.029, p = 0.200", size=4) +
#theme(legend.title=element_blank())
#theme(legend.position="bottom", legend.box = "horizontal")
  theme(legend.title=element_blank())

meta_phy_out2 <- pssd2veg(merged_phy_percent2)
shared_phy_out2 <- psotu2veg(merged_phy_percent2)


dist_bray2 <- vegdist(shared_phy_out2, method="bray", binary=FALSE)
#adonis2(formula = dist_bray2 ~ dam_ID, data = meta_phy_out2, permutations=10000)
adonis2(formula = dist_bray2 ~ Test + Weight, data = meta_phy_out2, permutations=10000, strata=meta_phy_out2$dam_ID)
adonis2(formula = dist_bray2 ~ Test + Weight, data = meta_phy_out2, permutations=10000)

#anosim(shared_phy_out2, meta_phy_out2$Test, distance = "bray", permutations = 9999)
#anosim(shared_phy_out2, meta_phy_out2$Test, distance = "bray", permutations = 9999)






Ce_pruned <- prune_taxa(names(sort(taxa_sums(merged_phy_percent2),TRUE)[1:20]), merged_phy_percent2)
sampleOrder = unique(sample_names(Ce_pruned))
taxaOrder = rev(unique(taxa_names(Ce_pruned)))

Ce_HM <- plot_heatmap(Ce_pruned , trans = NULL, sample.order = sampleOrder, taxa.order = taxaOrder, low="white", high="dodgerblue3", na.value="azure1",
                      title=NULL) +
  theme(axis.text.x=element_blank()) +
  #theme(axis.text.x=element_text(size=4)) +
  theme(axis.text.y=element_text(size=8)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  facet_wrap(~Order2, scales = "free_x", nrow = 1) +
  theme(strip.text.x = element_blank()) +
  theme(axis.title.y=element_blank())




Ce_bar <- plot_bar(SI_pruned , fill="Rank5") + facet_wrap(~Order2, scales = "free_x", nrow = 4) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
  theme(strip.text.x = element_blank())






library(ggpubr)
Ce_merged <- ggarrange(Ce_bar, Ce_HM, Ce_bar_BD,
                       common.legend = FALSE, legend = "right",
                       ncol = 3, nrow = 1, align="hv")












###Top 10 taxa


Ce_10 <- prune_taxa(names(sort(taxa_sums(merged_phy_percent2),TRUE)[1:10]), merged_phy_percent2)

# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(Ce_10) {
  sd <- sample_data(Ce_10)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(Ce_10) {
  OTU <- otu_table(Ce_10)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


Ce_10_meta <- pssd2veg(Ce_10)
Ce_10_shared <- psotu2veg(Ce_10)

df <- as.data.frame(Ce_10_shared)
meta <- as.data.frame(Ce_10_meta)

colnames(Ce_10_shared)


Se_abund <- as.data.frame(cbind(df$G1, df$G2, df$G3, df$G4, df$G5, df$G6, df$G7, df$G8, df$G9, df$G12, meta))
colnames(Se_abund)[1] <- "G1"
colnames(Se_abund)[2] <- "G2"
colnames(Se_abund)[3] <- "G3"
colnames(Se_abund)[4] <- "G4"
colnames(Se_abund)[5] <- "G5"
colnames(Se_abund)[6] <- "G6"
colnames(Se_abund)[7] <- "G7"
colnames(Se_abund)[8] <- "G8"
colnames(Se_abund)[9] <- "G9"
colnames(Se_abund)[10] <- "G12"





G1_p <- ggplot(Se_abund ,aes(x=Order, y=df$G1, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G2_p <- ggplot(Se_abund ,aes(x=Order, y=df$G2, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())



G3_p <- ggplot(Se_abund ,aes(x=Order, y=df$G3, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G4_p <- ggplot(Se_abund ,aes(x=Order, y=df$G4, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())

G5_p <- ggplot(Se_abund ,aes(x=Order, y=df$G5, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G6_p <- ggplot(Se_abund ,aes(x=Order, y=df$G6, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G7_p <- ggplot(Se_abund ,aes(x=Order, y=df$G7, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=5, alpha=0.5) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G8_p <- ggplot(Se_abund ,aes(x=Order, y=df$G8, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G9_p <- ggplot(Se_abund ,aes(x=Order, y=df$G9, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G12_p <- ggplot(Se_abund ,aes(x=Order, y=df$G12, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())




Ce_plot_data <- G12_p$data

#write.csv(Ce_plot_data,file="Ce_plot_data.csv")








Ce_weight_df <- data.frame(Ce_plot_data$G2, Ce_plot_data$G7, Ce_plot_data$G8, Ce_plot_data$Weight, Ce_plot_data$Test, Ce_plot_data$Order2)

colnames(Ce_weight_df)[1] <- "G2"
colnames(Ce_weight_df)[2] <- "G7"
colnames(Ce_weight_df)[3] <- "G8"
colnames(Ce_weight_df)[4] <- "Weight"
colnames(Ce_weight_df)[5] <- "Test"
colnames(Ce_weight_df)[6] <- "Order2"




G2_plot <- ggplot(Ce_weight_df, aes(Weight, G2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~Order2, scales = "free_x", nrow = 1)



G7_plot <- ggplot(Ce_weight_df, aes(Weight, G7)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~Order2, scales = "free_x", nrow = 1)


G8_plot <- ggplot(Ce_weight_df, aes(Weight, G8)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~Order2, scales = "free_x", nrow = 1)






Ce_Weight <- ggarrange(G2_plot, G7_plot, G8_plot, labels = c("Akkermansia_G2", "Duncaniella_G7", 	"Duncaniella_G8"),
                       common.legend = TRUE, legend = "right",
                       ncol = 1, nrow = 3, align="hv")











library(ggpubr)
Ce_barplot <- ggarrange(SI_Weight, Ce_Weight, Ce_Weight, labels = c("Akkermansia_G2", "Duncaniella_G7", 	"Duncaniella_G8"),
                        common.legend = TRUE, legend = "right",
                        ncol = 2, nrow = 5, align="hv")







library(ggpubr)
Ce_barplot <- ggarrange(G1_p, G2_p, G3_p, G4_p, G5_p, G6_p, G7_p, G8_p, G9_p, G12_p, labels = c("Muribaculum_G1", 	"Akkermansia_G2", 	"Lactobacillus_G3", 	"Muribaculum_G4", 	"Bacteroides_G5", 	"Lactobacillus_G6", 	"Duncaniella_G7", 	"Duncaniella_G8", "Acutalibacter_G9", "Lachnoclostridium_G12"),
                        common.legend = TRUE, legend = "right",
                        ncol = 2, nrow = 5, align="hv")



























###Colon

setwd("/Users/zebra/Desktop/GF/")

meta_ob1  <- read.table("/Users/zebra/Desktop/GF/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("/Users/zebra/Desktop/GF/taxa.txt"))


meta_ob1 <- meta_ob1[meta_ob1$Origin=='Neonatal' & meta_ob1$Order=="C_Co",]
#meta_ob1 <- meta_ob1[meta_ob1$Origin=='Neonatal' & meta_ob1$Order=="C_Co" & meta_ob1$Order2!="C" & meta_ob1$Order2!="D",]
#meta_ob1 <- meta_ob1[meta_ob1$Origin=='Neonatal' & meta_ob1$Order=="C_Co" & meta_ob1$Order2!="A" & meta_ob1$Order2!="B",]


nrow(meta_ob1)


#Get sample names
GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtu))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))


#Rarefy
#minimum <- min(sample_sums(phyloseqobj.f))
out3 <- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1278327, rngseed = 1, replace = FALSE)

#Count to percent
phy1 <- phyloseq_standardize_otu_abundance(out3, method = "total")



merged_phy_percent3 <- phyloseq_standardize_otu_abundance(phy1, method = "total")

#Co_phy <- merged_phy_percent3


# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent3) {
  sd <- sample_data(merged_phy_percent3)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent3) {
  OTU <- otu_table(merged_phy_percent3)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


meta_phy_out3 <- pssd2veg(merged_phy_percent3)
shared_phy_out3 <- psotu2veg(merged_phy_percent3)



dist_bray3 <- capscale(shared_phy_out3~1, distance="bray", binary=FALSE, eig=TRUE)

prop_explained <- summary(eigenvals(dist_bray3))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)


data.scores = scores(dist_bray3, display = "sites", choices = c(1:10))
#summary(data.scores)


df <- data.frame(data.scores)

library(glue)
labx3 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby3 <- c(glue("PCo 2 ({prop_explained [2]}%)"))


library(dplyr)
library(ggrepel)
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = shared_phy_out3)

#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)






Co_AB_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = meta_phy_out3$Order2, fill=meta_phy_out3$Order2, shape=meta_phy_out3$Order2)) +
  #scale_shape_manual(values = c(24, 24)) +
  scale_shape_manual(values = c(2, 2)) +
  geom_point(mapping = aes(colour = meta_phy_out3$Order2), size=5, alpha=1.0) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.1, show.legend=FALSE) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  #scale_color_manual(values=c("#B81A1A", "#283DA8", "#56D43A", "#5A5C59")) +
  #scale_fill_manual(values=c("pink", "dodgerblue", "#56D43A", "#5A5C59")) +
  #coord_fixed()+
  theme_classic() +
  labs(x=labx3, y=laby3) +
  #adds Weighted average score points for ASVs
  ##geom_point(data = WAscores_df, colour = "black",aes(x = MDS1, y = MDS2), inherit.aes = FALSE) +
  ##geom_text(data = WAscores_df, inherit.aes = FALSE, aes(x=MDS1, y= MDS2), colour="black", label = rownames(WAscores_df))+
  #ggtitle("Colon") +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=10, face=NULL), axis.text.y = element_text(color="black")) +
  #legend.background = element_rect(fill="white", color="black")
  #theme(legend.position = c(0.9, 0.4)) +
  #annotate("text", x=0.55, y=0.950, label = "Test: R2 = 0.198, p = 0.0003", size=4) +
  #annotate("text", x=0.55, y=0.850, label = "Weight: R2 = 0.020, p = 0.603", size=4) +
 #theme(legend.title=element_blank())
 #theme(legend.position="bottom", legend.box = "horizontal")
  theme(legend.title=element_blank())
  


meta_phy_out3 <- pssd2veg(merged_phy_percent3)
shared_phy_out3 <- psotu2veg(merged_phy_percent3)



dist_bray3 <- vegdist(shared_phy_out3, method="bray", binary=FALSE)
#adonis2(formula = dist_bray3 ~ dam_ID, data = meta_phy_out3, permutations=10000)
adonis2(formula = dist_bray3 ~ Test + Weight, data = meta_phy_out3, permutations=10000, strata=meta_phy_out3$dam_ID)
adonis2(formula = dist_bray3 ~ Test + Weight, data = meta_phy_out3, permutations=10000)



Co_pruned <- prune_taxa(names(sort(taxa_sums(merged_phy_percent3),TRUE)[1:20]), merged_phy_percent3)
sampleOrder = unique(sample_names(Co_pruned))
taxaOrder = rev(unique(taxa_names(Co_pruned)))


Co_HM <- plot_heatmap(Co_pruned , trans = NULL, sample.order = sampleOrder, taxa.order = taxaOrder, low="white", high="dodgerblue3", na.value="azure1",
                      title=NULL) +
  theme(axis.text.x=element_blank()) +
  #theme(axis.text.x=element_text(size=4)) +
  theme(axis.text.y=element_text(size=8)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  facet_wrap(~Order2, scales = "free_x", nrow = 1)



  #theme(strip.text.x = element_blank()) +
  #theme(axis.title.y=element_blank())



  
Co_bar <- plot_bar(SI_pruned , fill="Rank5") + facet_wrap(~Order2, scales = "free_x", nrow = 4) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
  theme(strip.text.x = element_blank())

  
  
  
  
  
  library(ggpubr)
  
  
  Ce_Co_PCoA_merged <- ggarrange(Ce_AB_plot, Ce_CD_plot, Co_AB_plot, Co_CD_plot,
                         common.legend = TRUE, legend = "right",
                         ncol = 1, nrow = 4)
  
  
  Ce_Co_bar <- ggarrange(Ce_bar, Co_bar,
                                 common.legend = TRUE, legend = "bottom",
                                 ncol = 1, nrow = 2)
  
  
  
  
  
  
  
  
  
  Co_merged <- ggarrange(Co_plot, Co_HM, Co_bar_BD,
                         common.legend = FALSE, legend = "right",
                         ncol = 3, nrow = 1, align="hv")
  
  
  
  Co_Ce_Bar_merged <- ggarrange(Ce_bar_BD, Co_bar_BD,
                         common.legend = TRUE, legend = "right",
                         ncol = 1, nrow = 2, align="hv")
  
  
  
  
  
  

  Ce_Co_PCoA_merged <- ggarrange(Ce_plot, Co_plot,
                               common.legend = TRUE, legend = "right",
                               ncol = 1, nrow = 2, align="hv")
  
  
  
  

  Ce_Co_HM_merged <- ggarrange(Ce_HM, Co_HM,
                         common.legend = TRUE, legend = "right",
                         ncol = 1, nrow = 2, align="hv")
  
  
  Ce_Co_bar_merged <- ggarrange(Ce_bar_BD, Co_bar_BD,
                                 common.legend = TRUE, legend = "right",
                                 ncol = 1, nrow = 2, align="hv")
  
  
  
  
  



###Top 10 taxa


Co_10 <- prune_taxa(names(sort(taxa_sums(merged_phy_percent3),TRUE)[1:10]), merged_phy_percent3)

# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(Co_10) {
  sd <- sample_data(Co_10)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(Co_10) {
  OTU <- otu_table(Co_10)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


Co_10_meta <- pssd2veg(Co_10)
Co_10_shared <- psotu2veg(Co_10)

df <- as.data.frame(Co_10_shared)
meta <- as.data.frame(Co_10_meta)

colnames(Co_10_shared)


Se_abund <- as.data.frame(cbind(df$G1, df$G2, df$G3, df$G4, df$G5, df$G6, df$G7, df$G8, df$G9, df$G11, meta))
colnames(Se_abund)[1] <- "G1"
colnames(Se_abund)[2] <- "G2"
colnames(Se_abund)[3] <- "G3"
colnames(Se_abund)[4] <- "G4"
colnames(Se_abund)[5] <- "G5"
colnames(Se_abund)[6] <- "G6"
colnames(Se_abund)[7] <- "G7"
colnames(Se_abund)[8] <- "G8"
colnames(Se_abund)[9] <- "G9"
colnames(Se_abund)[10] <- "G11"





G1_p <- ggplot(Se_abund ,aes(x=Order, y=df$G1, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G2_p <- ggplot(Se_abund ,aes(x=Order, y=df$G2, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())



G3_p <- ggplot(Se_abund ,aes(x=Order, y=df$G3, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G4_p <- ggplot(Se_abund ,aes(x=Order, y=df$G4, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())

G5_p <- ggplot(Se_abund ,aes(x=Order, y=df$G5, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G6_p <- ggplot(Se_abund ,aes(x=Order, y=df$G6, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G7_p <- ggplot(Se_abund ,aes(x=Order, y=df$G7, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=5, alpha=0.5) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G8_p <- ggplot(Se_abund ,aes(x=Order, y=df$G8, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G9_p <- ggplot(Se_abund ,aes(x=Order, y=df$G9, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())


G11_p <- ggplot(Se_abund ,aes(x=Order, y=df$G11, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=4, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())







Co_plot_data <- G11_p$data

#write.csv(Co_plot_data,file="Co_plot_data.csv")




Co_weight_df <- data.frame(Co_plot_data$G2, Co_plot_data$G7, Co_plot_data$G8, Co_plot_data$Weight, Co_plot_data$Test, Co_plot_data$Order2)

colnames(Co_weight_df)[1] <- "G2"
colnames(Co_weight_df)[2] <- "G7"
colnames(Co_weight_df)[3] <- "G8"
colnames(Co_weight_df)[4] <- "Weight"
colnames(Co_weight_df)[5] <- "Test"
colnames(Co_weight_df)[6] <- "Order2"




G2_plot <- ggplot(Co_weight_df, aes(Weight, G2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~Order2, scales = "free_x", nrow = 1)



G7_plot <- ggplot(Co_weight_df, aes(Weight, G7)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~Order2, scales = "free_x", nrow = 1)


G8_plot <- ggplot(Co_weight_df, aes(Weight, G8)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~Order2, scales = "free_x", nrow = 1)






Co_Weight <- ggarrange(G2_plot, G7_plot, G8_plot, labels = c("Akkermansia_G2", "Duncaniella_G7", 	"Duncaniella_G8"),
                       common.legend = TRUE, legend = "right",
                       ncol = 1, nrow = 3, align="hv")









library(ggpubr)
Weight_plot <- ggarrange(SI_Weight, Ce_Weight, Co_Weight, labels = c("SI", "Ce", "Co"),
                        common.legend = TRUE, legend = "right",
                        ncol = 1, nrow = 3, align="hv")








Colon_pruned_BD <- prune_taxa(names(sort(taxa_sums(percent_Colon_BD),TRUE)[1:30]), percent_Cecum_BD)
Ce_plot_BD <- plot_heatmap(Cecum_pruned_BD , method = "PCoA", distance = "bray", "Injection", "Rank1",
                           trans = NULL, sample.order="Injection", low="azure1", high="dodgerblue3", na.value="white",
                           title="Neonatal Cecum B vs D") +
  theme(axis.text.x=element_text(size=4)) +
  theme(axis.text.y=element_text(size=8)) +
  theme(axis.title.x=element_blank())












###Alpha diversity


# Load necessary libraries
library(ggplot2)

# Read the data from the file
meta_ob1 <- read.table("/Users/zebra/Desktop/GW_2024/GW_meta.txt", header = TRUE, row.names = 1)



chao_plot <- ggplot(data = meta_ob1,
                    aes(x=as.factor(Week), y=chao, fill = factor(Treatment))) +
  geom_boxplot(fatten=NULL, outlier.shape = NA, position = position_dodge(width=0.85)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.333, dodge.width = 0.85), aes(fill =Treatment), pch = 21) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), position = position_dodge(width = 0.83), width = 0.74, size=1) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "#53B400", "#A58AFF", "#FFD700", "#00bf7d","purple", "blue","black","magenta","pink","cyan","gray","orange","brown")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "#53B400", "#A58AFF", "#FFD700", "#00bf7d","purple", "black","magenta","pink","cyan","gray","orange","brown","purple")) +
  ylab("Chao1") +
  theme_classic() +
  #scale_y_continuous(limits = c(0,800), labels = scales::number_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = 7)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9),
        legend.title= element_blank()) +
  facet_wrap(~Type, scales = "free_x", ncol = 2, nrow = 1)


shannon_plot <- ggplot(data = meta_ob1,
                       aes(x=as.factor(Week), y=shannon, fill = factor(Treatment))) +
  geom_boxplot(fatten=NULL, outlier.shape = NA, position = position_dodge(width=0.85)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.333, dodge.width = 0.85), aes(fill =Treatment), pch = 21) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), position = position_dodge(width = 0.83), width = 0.74, size=1) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "#53B400", "#A58AFF", "#FFD700", "#00bf7d","purple", "blue","black","magenta","pink","cyan","gray","orange","brown")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "#53B400", "#A58AFF", "#FFD700", "#00bf7d","purple", "black","magenta","pink","cyan","gray","orange","brown","purple")) +
  ylab("Shannon") +
  theme_classic() +
  #scale_y_continuous(limits = c(0,800), labels = scales::number_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = 7)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9),
        legend.title= element_blank()) +
  facet_wrap(~Type, scales = "free_x", ncol = 2, nrow = 1)



invsimpson_plot <- ggplot(data = meta_ob1,
                          aes(x=as.factor(Week), y=invsimpson, fill = factor(Treatment))) +
  geom_boxplot(fatten=NULL, outlier.shape = NA, position = position_dodge(width=0.85)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.333, dodge.width = 0.85), aes(fill =Treatment), pch = 21) +
  stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), position = position_dodge(width = 0.83), width = 0.74, size=1) +
  scale_color_manual(values=c("#5B84B1", "#FC766A", "#53B400", "#A58AFF", "#FFD700", "#00bf7d","purple", "blue","black","magenta","pink","cyan","gray","orange","brown")) +
  scale_fill_manual(values=c("#5B84B1", "#FC766A", "#53B400", "#A58AFF", "#FFD700", "#00bf7d","purple", "black","magenta","pink","cyan","gray","orange","brown","purple")) +
  ylab("Inverse Simpson") +
  theme_classic() +
  #scale_y_continuous(limits = c(0,800), labels = scales::number_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = 7)) +
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(color="black", size=9),
        axis.text.x = element_text(color="black", size=8.5),
        axis.text.y = element_text(color="black", size=9),
        legend.title= element_blank()) +
  facet_wrap(~Type, scales = "free_x", ncol = 2, nrow = 1)







library(ggpubr)

alpha_merged <- ggarrange(chao_plot, shannon_plot, invsimpson_plot ,
                          common.legend = TRUE, legend = "right",
                          ncol = 1, nrow = 3, align="hv")













library(lme4)
library(emmeans)
library(car)

#chao
chao_mod <- lmer(as.numeric(meta_ob1$chao) ~ Treatment*Week*Type + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=FALSE)
Anova(chao_mod, type="III")
chao_out <- emmeans(chao_mod, list(pairwise ~ Treatment*Week*Type), adjust = "tukey")
chao_out <- as.data.frame(chao_out[2])
write.csv (chao_out, file = "chao_tukey.csv")


#shannon
shannon_mod <- lmer(as.numeric(meta_ob1$shannon) ~ Treatment*Week*Type + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=FALSE)
Anova(shannon_mod, type="III")
shannon_out <- emmeans(shannon_mod, list(pairwise ~ Treatment*Week*Type), adjust = "tukey")
shannon_out <- as.data.frame(shannon_out[2])
write.csv (shannon_out, file = "shannon_tukey.csv")


#invsimpson
invsimpson_mod <- lmer(as.numeric(meta_ob1$invsimpson) ~ Treatment*Week*Type + (1|meta_ob1$Mouse_ID), data = meta_ob1, REML=FALSE)
Anova(invsimpson_mod, type="III")
invsimpson_out <- emmeans(invsimpson_mod, list(pairwise ~ Treatment*Week*Type), adjust = "tukey")
invsimpson_out <- as.data.frame(invsimpson_out[2])
write.csv (invsimpson_out, file = "invsimpson_tukey.csv")


































































#######










setwd("/Users/zebra/Desktop/GF/")

meta_ob1  <- read.table("/Users/zebra/Desktop/GF/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("/Users/zebra/Desktop/GF/taxa.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Origin=='Neonatal' & meta_ob1$Order!="C_Co",]

nrow(meta_ob1)


#Get sample names
GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtu))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))


#Rarefy
#minimum <- min(sample_sums(phyloseqobj.f))
out1 <- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1278327, rngseed = 1, replace = FALSE)
#out1<- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1000, rngseed = 1, replace = FALSE)

#Count to percent
phy1 <- phyloseq_standardize_otu_abundance(out1, method = "total")


merged_phy_percent2 <- phyloseq_standardize_otu_abundance(phy1, method = "total")



# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent2) {
  sd <- sample_data(merged_phy_percent2)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent2) {
  OTU <- otu_table(merged_phy_percent2)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


meta_phy_out2 <- pssd2veg(merged_phy_percent2)
shared_phy_out2 <- psotu2veg(merged_phy_percent2)



dist_bray2 <- capscale(shared_phy_out2~1, distance="bray", binary=FALS, eig=TRUE)

prop_explained <- summary(eigenvals(dist_bray2))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)


data.scores = scores(dist_bray2, display = "sites", choices = c(1:10))
#summary(data.scores)


df <- data.frame(data.scores)

library(glue)
labx2 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby2 <- c(glue("PCo 2 ({prop_explained [2]}%)"))


library(dplyr)
library(ggrepel)
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = shared_phy_out2)

#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)


Co_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = meta_phy_out2$Test, fill=meta_phy_out2$Test)) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.25, show.legend=FALSE, shape=21) +
  geom_point(mapping = aes(colour = meta_phy_out2$Test), size=5, alpha=1.0) +
  scale_color_manual(values=c("#B81A1A", "#283DA8", "#56D43A", "#5A5C59")) +
  scale_fill_manual(values=c("pink", "dodgerblue", "#56D43A", "#5A5C59")) +
  #coord_fixed()+
  theme_classic() +
  labs(x=labx2, y=laby2) +
  #adds Weighted average score points for ASVs
  ##geom_point(data = WAscores_df, colour = "black",aes(x = MDS1, y = MDS2), inherit.aes = FALSE) +
  ##geom_text(data = WAscores_df, inherit.aes = FALSE, aes(x=MDS1, y= MDS2), colour="black", label = rownames(WAscores_df))+
  ggtitle("Colon") +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=12, face="bold"), axis.text.y = element_text(color="black")) +
  #legend.background = element_rect(fill="white", color="black")
  #theme(legend.position = c(0.9, 0.4)) +
  annotate("text", x=0.55, y=1.0, label = "R2 = 0.152, p < 0.0001", size=4)
 #theme(legend.title=element_blank())
#theme(legend.position="bottom", legend.box = "horizontal")


meta_phy_out2 <- pssd2veg(merged_phy_percent2)
shared_phy_out2 <- psotu2veg(merged_phy_percent2)


#dist_bray2 <- vegdist(shared_phy_out2, method="bray", binary=FALSE)
#adonis2(formula = dist_bray2 ~ Test, data = meta_phy_out2, permutations=10000)








###Colon

setwd("/Users/zebra/Desktop/GF/")

meta_ob1  <- read.table("/Users/zebra/Desktop/GF/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("/Users/zebra/Desktop/GF/taxa.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Origin=='Neonatal' & meta_ob1$Order!="A_SI",]

nrow(meta_ob1)


#Get sample names
GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtu))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))


#Rarefy
#minimum <- min(sample_sums(phyloseqobj.f))
out1 <- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1278327, rngseed = 1, replace = FALSE)
#out1<- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1000, rngseed = 1, replace = FALSE)

#Count to percent
phy1 <- phyloseq_standardize_otu_abundance(out1, method = "total")


merged_phy_percent3 <- phyloseq_standardize_otu_abundance(phy1, method = "total")



# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent3) {
  sd <- sample_data(merged_phy_percent3)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent3) {
  OTU <- otu_table(merged_phy_percent3)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


meta_phy_out3 <- pssd2veg(merged_phy_percent3)
shared_phy_out3 <- psotu2veg(merged_phy_percent3)



dist_bray3 <- capscale(shared_phy_out3~1, distance="bray", binary=FALS, eig=TRUE)

prop_explained <- summary(eigenvals(dist_bray3))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)


data.scores = scores(dist_bray3, display = "sites", choices = c(1:10))
#summary(data.scores)


df <- data.frame(data.scores)

library(glue)
labx3 <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby3 <- c(glue("PCo 2 ({prop_explained [2]}%)"))


library(dplyr)
library(ggrepel)
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = shared_phy_out3)

#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)


Colon_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = meta_phy_out3$Test, fill=meta_phy_out3$Test)) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.25, show.legend=FALSE) +
  geom_point(mapping = aes(colour = meta_phy_out3$Test), size=5, shape = 24) +
  scale_color_manual(values=c("#B81A1A", "#283DA8", "#56D43A", "#5A5C59")) +
  scale_fill_manual(values=c("pink", "dodgerblue", "#56D43A", "#5A5C59")) +
  #coord_fixed()+
  theme_classic() +
  labs(x=labx3, y=laby3) +
  #adds Weighted average score points for ASVs
  ##geom_point(data = WAscores_df, colour = "black",aes(x = MDS1, y = MDS2), inherit.aes = FALSE) +
  ##geom_text(data = WAscores_df, inherit.aes = FALSE, aes(x=MDS1, y= MDS2), colour="black", label = rownames(WAscores_df))+
  ggtitle("Colon") +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=12, face="bold"), axis.text.y = element_text(color="black")) +
  #legend.background = element_rect(fill="white", color="black")
  #theme(legend.position = c(0.9, 0.4)) +
  annotate("text", x=0.55, y=1.0, label = "R2 = 0.152, p < 0.0001", size=4)
#theme(legend.title=element_blank())
#theme(legend.position="bottom", legend.box = "horizontal")


meta_phy_out3 <- pssd2veg(merged_phy_percent3)
shared_phy_out3 <- psotu2veg(merged_phy_percent3)


#dist_bray3 <- vegdist(shared_phy_out3, method="bray", binary=FALSE)
#adonis2(formula = dist_bray3 ~ Test, data = meta_phy_out3, permutations=10000)








library(ggpubr)
PCoA_plot_AC <- ggarrange(All_plot, Cecum_plot, Colon_plot,
                   common.legend = TRUE, legend = "bottom",
                   ncol = 1, nrow = 3, align="hv")

PCoA_plot_BD <- ggarrange(All_plot, Cecum_plot, Colon_plot, labels = c("1", "2", "3"),
                          common.legend = TRUE, legend = "bottom",
                          ncol = 1, nrow = 3, align="hv")




merged_PCoA <- ggarrange(PCoA_plot_BD, PCoA_plot_AC, labels = c("BD", "AC"),
                          common.legend = TRUE, legend = "bottom",
                          ncol = 1, nrow = 3, align="hv")




















###All samples for global





merged_phy_counts <- merge_phyloseq(out1, out2, out3)

merged_phy_percent <- phyloseq_standardize_otu_abundance(merged_phy_counts, method = "total")




# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent) {
  sd <- sample_data(merged_phy_percent)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent) {
  OTU <- otu_table(merged_phy_percent)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


meta_phy_out_All <- pssd2veg(merged_phy_percent)
shared_phy_out_All <- psotu2veg(merged_phy_percent)


nrow(meta_phy_out_All)



merged_phy_percent_CD_1 <- subset_samples(merged_phy_percent , Order2!="A")
merged_phy_percent_AB <- subset_samples(merged_phy_percent_CD_1 , Order2!="B")




# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent_AB) {
  sd <- sample_data(merged_phy_percent_AB)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent_AB) {
  OTU <- otu_table(merged_phy_percent_AB)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


meta_phy_out_All_AB <- pssd2veg(merged_phy_percent_AB)
shared_phy_out_All_AB <- psotu2veg(merged_phy_percent_AB)


nrow(meta_phy_out_All_AB)





dist_bray <- vegdist(shared_phy_out_All_AB, method="bray", binary=FALSE)
#adonis2(formula = dist_bray ~ Test + Type + Weight, data = meta_phy_out_All_AB, permutations=10000, strata=meta_phy_out_All_AB$dam_ID)


adonis2(formula = dist_bray ~ Test + Type, data = meta_phy_out_All_AB, permutations=10000)
















merged_phy_counts <- merge_phyloseq(out3)

merged_phy_percent <- phyloseq_standardize_otu_abundance(merged_phy_counts, method = "total")




# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent) {
  sd <- sample_data(merged_phy_percent)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent) {
  OTU <- otu_table(merged_phy_percent)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


meta_phy_out_All <- pssd2veg(merged_phy_percent)
shared_phy_out_All <- psotu2veg(merged_phy_percent)


nrow(meta_phy_out_All)



merged_phy_percent_AB_1 <- subset_samples(merged_phy_percent , Order2!="A")
merged_phy_percent_CD <- subset_samples(merged_phy_percent_AB_1, Order2!="B")






# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent_CD) {
  sd <- sample_data(merged_phy_percent_CD)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent_CD) {
  OTU <- otu_table(merged_phy_percent_CD)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


meta_phy_out_All_CD <- pssd2veg(merged_phy_percent_CD)
shared_phy_out_All_CD <- psotu2veg(merged_phy_percent_CD)


nrow(meta_phy_out_All_CD)





dist_bray <- vegdist(shared_phy_out_All_CD , method="bray", binary=FALSE)
adonis2(formula = dist_bray ~ Test + Type + Weight, data = meta_phy_out_All_CD, permutations=10000, strata=meta_phy_out_All_CD$dam_ID)


adonis2(formula = dist_bray ~ Test + Weight, data = meta_phy_out_All_CD, permutations=10000, strata=meta_phy_out_All_CD$dam_ID)












#merged_phy_counts <- merge_samples(merged_phy_counts, meta_phy_out_All$Order2, fun=mean)









Ce_Co_pruned <- prune_taxa(names(sort(taxa_sums(merged_phy_counts ),TRUE)[1:20]), merged_phy_counts)
sampleOrder = unique(sample_names(Ce_Co_pruned))
taxaOrder = rev(unique(taxa_names(Ce_Co_pruned)))


Ce_Co_HM <- plot_heatmap(Ce_Co_pruned , trans = NULL, sample.order = sampleOrder, taxa.order = taxaOrder, low="white", high="dodgerblue3", na.value="azure1",
                      title=NULL) +
                      facet_wrap(~Order3, scales = "free_x", ncol = 2, dir="h") +
  theme(axis.text.x=element_blank()) +
  #theme(axis.text.x=element_text(size=4)) +
  theme(axis.text.y=element_text(size=8)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
#theme(strip.text = element_text(size = rel(0.1))) +
theme(strip.text.x = element_blank()) +
theme(axis.title.y=element_blank()) +
  theme(strip.background =element_rect(fill="white"))



plot_heatmap(Ce_Co_pruned , method = "PCoA", distance = "bray", "Injection", "Rank8",
             trans = NULL, sample.order="Injection", low="azure1", high="dodgerblue3", na.value="white")









Ce_Co_pruned <- merge_samples(merged_phy_counts, meta_phy_out_All$Order3, fun=sum)


merged_phy_counts <- phyloseq_standardize_otu_abundance(merged_phy_counts, method = "total")




prev_data <- psmelt(merged_phy_counts)

write.csv(as.matrix(prev_data), "prev_data2.csv")

#merged_phy_counts <- merge_samples(merged_phy_counts, meta_phy_out_All$Order2, fun=mean)

plot_out <- plot_bar(Ce_Co_pruned , fill="Rank2")








plot_out <- plot_bar(Ce_Co_pruned , fill="Rank2") + facet_wrap(~Order3, scales = "free_x", nrow = 6) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
  theme(strip.text = element_text(size = rel(0.1))) +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_blank()) +
  #theme(axis.text.y=element_text(size=8)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  theme(legend.position="bottom")





plot_data <- ggplot_build(plot_out)$plot$data

write.csv(as.matrix(plot_data), "barplot_data_taxa.csv")






  theme(legend.position="bottom")
theme(strip.text.x = element_blank()) +
  theme(legend.position="bottom")







Ce_Co_HM <- plot_heatmap(SI_pruned , trans = NULL, sample.order = sampleOrder, taxa.order = taxaOrder, low="white", high="dodgerblue3", na.value="azure1")
  theme(axis.text.x=element_blank()) +
  #theme(axis.text.x=element_text(size=4)) +
  theme(axis.text.y=element_text(size=8)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  facet_wrap(~Order2, scales = "free_x", nrow = 1)
#theme(strip.text.x = element_blank()) +
#theme(axis.title.y=element_blank())











merged_phy_counts <- merge_phyloseq(out2, out3)
#merged_phy_counts <- merge_phyloseq(out1)
#merged_phy_counts <- merge_phyloseq(out2)
#merged_phy_counts <- merge_phyloseq(out3)

merged_phy_counts <- subset_samples(merged_phy_counts , Order2!="C")
merged_phy_counts <- subset_samples(merged_phy_counts , Order2!="D")



# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_counts) {
  sd <- sample_data(merged_phy_counts)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_counts) {
  OTU <- otu_table(merged_phy_counts)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


meta_phy_out_All <- pssd2veg(merged_phy_counts)
shared_phy_out_All <- psotu2veg(merged_phy_counts)



dist_bray <- capscale(shared_phy_out_All~1, distance="bray", binary=FALSE, eig=TRUE)

prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)


data.scores = scores(dist_bray, display = "sites", choices = c(1:10))
#summary(data.scores)


df <- data.frame(data.scores)

library(glue)
labx <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby <- c(glue("PCo 2 ({prop_explained [2]}%)"))


library(dplyr)
library(ggrepel)
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = shared_phy_out_All)

#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)


All_plot <- ggplot(df, aes(x = MDS1, y = MDS2, colour = meta_phy_out_All$Order2, fill=meta_phy_out_All$Order2, shape=meta_phy_out_All$Order2)) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  geom_point(mapping = aes(colour = meta_phy_out_All$Order2), size=5, alpha=1.0) +
  stat_ellipse(geom="polygon", level=0.8, alpha=0.25, show.legend=FALSE) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  #scale_color_manual(values=c("#B81A1A", "#283DA8", "#56D43A", "#5A5C59")) +
  #scale_fill_manual(values=c("pink", "dodgerblue", "#56D43A", "#5A5C59")) +
  #coord_fixed()+
  theme_classic() +
  labs(x=labx, y=laby) +
  #adds Weighted average score points for ASVs
  ##geom_point(data = WAscores_df, colour = "black",aes(x = MDS1, y = MDS2), inherit.aes = FALSE) +
  ##geom_text(data = WAscores_df, inherit.aes = FALSE, aes(x=MDS1, y= MDS2), colour="black", label = rownames(WAscores_df))+
  #ggtitle("Small intestine") +
  theme(legend.position="bottom", legend.box = "vertical") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=10, face=NULL), axis.text.y = element_text(color="black")) +
  #legend.background = element_rect(fill="white", color="black")
  #theme(legend.position = c(0.9, 0.4)) +
  annotate("text", x=0.55, y=0.950, label = "Test: R2 = 0.055, p = 0.696", size=4) +
  annotate("text", x=0.55, y=0.850, label = "Weight: R2 = 0.081, p = 0.005", size=4) +
  #theme(legend.title=element_blank())
  #theme(legend.position="bottom", legend.box = "horizontal") +
  theme(legend.title=element_blank())




dist_bray <- vegdist(shared_phy_out_All, method="bray", binary=FALSE)
adonis2(formula = dist_bray ~ Test + Weight, data = meta_phy_out_All, permutations=10000, strata=meta_phy_out_All$dam_ID)





###Top 25 taxa



merged_phy_counts <- merge_phyloseq(out3)
#merged_phy_counts <- merge_phyloseq(out3)

merged_phy_counts <- subset_samples(merged_phy_counts , Order2!="C")
merged_phy_counts <- subset_samples(merged_phy_counts , Order2!="D")



All_25 <- prune_taxa(names(sort(taxa_sums(merged_phy_counts),TRUE)[1:25]), merged_phy_counts)



# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(All_25) {
  sd <- sample_data(All_25)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(All_25) {
  OTU <- otu_table(All_25)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


All_25_meta <- pssd2veg(All_25)
All_25_shared <- psotu2veg(All_25)


rownames(All_25_meta)  == rownames(All_25_shared)


ASV_table <- as.data.frame(All_25_shared)

meta <- as.data.frame(All_25_meta)


rownames(ASV_table) == rownames(meta)

meta$Test <- as.factor(meta$Test)


nrow(ASV_table)
nrow(meta)

levels(meta$Test)




PBS_Iso_model_out <- MW_model_function(meta=meta, ASV_table = ASV_table, Group_col = "Test", reference = "PBS_Iso")
#PBS_anti6_model_out <- MW_model_function(meta=meta, ASV_table = ASV_table, Group_col = "Test", reference = "PBS_anti6")




write.csv(PBS_Iso_model_out,file="Ce_AB_holm_MW.csv")
write.csv(PBS_Iso_model_out,file="Co_AB_holm_MW.csv")

write.csv(PBS_anti6_model_out,file="Ce_CD_holm_MW.csv")
write.csv(PBS_anti6_model_out,file="Co_CD_holm_MW.csv")






























###selected




All_25_shared <- All_25_shared[, c("G5",	"G10",	"G8",	"G7",	"G1",	"G3",	"G20") ]

#All_25_shared <- All_25_shared[, c("G11", "G2",  "G19") ]

write.csv(All_25_shared,file="Ce_Co_top_25_shared.csv")
write.csv(All_25_meta,file="Ce_Co_top_25_meta.csv")


ASV_table <- as.data.frame(All_25_shared)
#ASV_table <- ASV_table + 0.0000000001 

meta <- as.data.frame(All_25_meta)

#meta$Test <- factor(meta$Test, levels = c("PBS_Iso", "IL1_Iso"))
meta$Test <- factor(meta$Test, levels = c("PBS_Iso", "IL1_Iso"))


meta$Test <- as.factor(meta$Test)
levels(meta$Test)


nrow(meta)


#############################

rownames(ASV_table) == rownames(meta)



#note that metadata and ASV table must have same rownames

library(glmmTMB)
#sets seed for reproducibility
set.seed(1)
#instantiates res which holds the results in the outer loop
res=NULL
#outer loop
#for all groups except the first (i.e. the control group)
for(group in levels(meta[[Group_col]])[levels(meta[[Group_col]]) != reference]){
  
  #subsets metadata on rows by `Group` column to `Control`s and `group`
  #meta_subset=meta[meta$Test%in%c("PBS_Iso",group),]
  meta_subset=meta[meta$Test%in%c("PBS_Iso",group),]
  
  
  #`tres` instantiated. Temporary results
  tres=NULL
  #for each ASV in `ASV_table` with ASVs on columns
  for (ASV in colnames(ASV_table)){
    ##NB binary
    #copies over one ASV column to subsetted metadata and names it column `Y`
    meta_subset$Y=ASV_table[rownames(meta_subset),ASV]
    #fits a generalized linear mixed model using template model builder for negative binomial model with offset by the log of the total number of reads in a sample
    #may not work for non-mixed effects
    #modnb=glmmTMB::glmmTMB(Y~Group+GA2+BMI60+Nulliparity+Preterm2+(1|Main_Index), data = meta_subset, offset=log(nR),family=glmmTMB::nbinom2)
    #modnb=glmmTMB::glmmTMB(Y~Origin + Test + (1|Mouse_ID), data = meta_subset, family=glmmTMB::nbinom2)
    #modnb=glmmTMB::glmmTMB(Y ~ Test + Type + Weight + (1|dam_ID) + (1|pup_ID), data = meta_subset, family=glmmTMB::nbinom2)
    #modnb=glmmTMB::glmmTMB(Y ~ Test + Type + Weight + (1|dam_ID) + (1|pup_ID), data = meta_subset, family=glmmTMB::nbinom2)
    modnb=glmmTMB::glmmTMB(Y ~ Test + Weight + (1|dam_ID), data = meta_subset, family=glmmTMB::nbinom2)
        #gets the estimate and P value
    nbb=summary(modnb)$coef$cond[2,c("Estimate","Pr(>|z|)")]
    #stores the temporary results into `tres`
    tres=rbind(tres,nbb)
    #lets you know how far along the loop is by printing the ASV
    print(ASV)
  }
  #correction for false discovery rate
  tres=cbind(tres,p.adjust(tres[,2],"fdr"))
  #renames the column names of `tres`
  colnames(tres)<-c(paste("Coef_",group,"vsControl",sep=""),paste("P_",group,"vsControl",sep=""),paste("Q_",group,"vsControl",sep=""))
  #stores looped through ASV model results into `res`
  res=cbind(res,tres)
}
#adds a Taxonomy column to the results
res=data.frame(Taxonomy=colnames(ASV_table),res)
#reorders results by p value in the `paste("P_",group,"vsControl",sep="")` column
#rename to group of interest
#res=res[order(res$P_PTL_PPROMvsControl),]
write.csv(res,file="new.csv")





library(car)
mod1 <- glmmTMB(Y ~ Test + Type + Weight + (1|dam_ID), data = meta_subset, family=glmmTMB::nbinom2)
Anova(mod1, type="III", test="Chisq")

diagnose(mod1)















###aov



#note that metadata and ASV table must have same rownames

library(car)
#sets seed for reproducibility
set.seed(1)
#instantiates res which holds the results in the outer loop
res=NULL
#outer loop
#for all groups except the first (i.e. the control group)
for(group in levels(meta$Test)[-1]){
  
  #subsets metadata on rows by `Group` column to `Control`s and `group`

  meta_subset=meta[meta$Test%in%c("PBS_Iso",group),]
  
  
  #`tres` instantiated. Temporary results
  tres=NULL
  #for each ASV in `ASV_table` with ASVs on columns
  for (ASV in colnames(ASV_table)){
    ##NB binary
    #copies over one ASV column to subsetted metadata and names it column `Y`
    meta_subset$Y=ASV_table[rownames(meta_subset),ASV]

    #modnb=glmmTMB::glmmTMB(Y ~ Test + Weight + (1|dam_ID), data = meta_subset, family=glmmTMB::nbinom2)
    modnb=aov(Y ~ Test + Weight, data = meta_subset)
    #gets the estimate and P value
    
    #nbb=summary(modnb)$coef$cond[2,c("Estimate","Pr(>|z|)")]
    
    nbb=summary(modnb)[[1]][["Pr(>F)"]][1]
    
    #stores the temporary results into `tres`
    tres=rbind(tres,nbb)
    #lets you know how far along the loop is by printing the ASV
    print(ASV)
  }
  #correction for false discovery rate
  tres=cbind(tres,p.adjust(tres[,2],"fdr"))
  #renames the column names of `tres`
  colnames(tres)<-c(paste("Coef_",group,"vsControl",sep=""),paste("P_",group,"vsControl",sep=""),paste("Q_",group,"vsControl",sep=""))
  #stores looped through ASV model results into `res`
  res=cbind(res,tres)
}
#adds a Taxonomy column to the results
res=data.frame(Taxonomy=colnames(ASV_table),res)
#reorders results by p value in the `paste("P_",group,"vsControl",sep="")` column
#rename to group of interest
#res=res[order(res$P_PTL_PPROMvsControl),]
write.csv(res,file="aov.csv")








Ce_Co_top_25_df <- read.table("/Users/zebra/Desktop/GF/Ce_Co_top_25.txt", header=T, row.names=1)
Ce_Co_top_25_df <- Ce_Co_top_25_df[Ce_Co_top_25_df$Origin=='Neonatal' & Ce_Co_top_25_df$Order2!="C" & Ce_Co_top_25_df$Order2!="D",]




summary(G1_mod)[[1]][["Pr(>F)"]][1]




coef(summary.lm())


G1_mod <- aov(data=Ce_Co_top_25_df,G1~ Test + Weight)

summary(G1_mod)

str(G1_mod)

G1_mod <- aov(data=Ce_Co_top_25_df,G1~ Test + Weight)
G2_mod <- aov(data=Ce_Co_top_25_df,G2~ Test + Weight)
G3_mod <- aov(data=Ce_Co_top_25_df,G3~ Test + Weight)
G4_mod <- aov(data=Ce_Co_top_25_df,G4~ Test + Weight)
G5_mod <- aov(data=Ce_Co_top_25_df,G5~ Test + Weight)
G6_mod <- aov(data=Ce_Co_top_25_df,G6~ Test + Weight)
G7_mod <- aov(data=Ce_Co_top_25_df,G7~ Test + Weight)
G8_mod <- aov(data=Ce_Co_top_25_df,G8~ Test + Weight)
G9_mod <- aov(data=Ce_Co_top_25_df,G9~ Test + Weight)
G10_mod <- aov(data=Ce_Co_top_25_df,G10~ Test + Weight)
G11_mod <- aov(data=Ce_Co_top_25_df,G11~ Test + Weight)
G12_mod <- aov(data=Ce_Co_top_25_df,G12~ Test + Weight)
G13_mod <- aov(data=Ce_Co_top_25_df,G13~ Test + Weight)
G14_mod <- aov(data=Ce_Co_top_25_df,G14~ Test + Weight)
G15_mod <- aov(data=Ce_Co_top_25_df,G15~ Test + Weight)
G16_mod <- aov(data=Ce_Co_top_25_df,G16~ Test + Weight)
G18_mod <- aov(data=Ce_Co_top_25_df,G18~ Test + Weight)
G19_mod <- aov(data=Ce_Co_top_25_df,G19~ Test + Weight)
G20_mod <- aov(data=Ce_Co_top_25_df,G20~ Test + Weight)
G21_mod <- aov(data=Ce_Co_top_25_df,G21~ Test + Weight)
G22_mod <- aov(data=Ce_Co_top_25_df,G22~ Test + Weight)
G26_mod <- aov(data=Ce_Co_top_25_df,G26~ Test + Weight)
G29_mod <- aov(data=Ce_Co_top_25_df,G29~ Test + Weight)
G30_mod <- aov(data=Ce_Co_top_25_df,G30~ Test + Weight)
G32_mod <- aov(data=Ce_Co_top_25_df,G32~ Test + Weight)


library(car)

G1_summary <- Anova(G1_mod, type="III", test="Chisq")
G2_summary <- Anova(G2_mod, type="III", test="Chisq")
G3_summary <- Anova(G3_mod, type="III", test="Chisq")
G4_summary <- Anova(G4_mod, type="III", test="Chisq")
G5_summary <- Anova(G5_mod, type="III", test="Chisq")
G6_summary <- Anova(G6_mod, type="III", test="Chisq")
G7_summary <- Anova(G7_mod, type="III", test="Chisq")
G8_summary <- Anova(G8_mod, type="III", test="Chisq")
G9_summary <- Anova(G9_mod, type="III", test="Chisq")
G10_summary <- Anova(G10_mod, type="III", test="Chisq")
G11_summary <- Anova(G11_mod, type="III", test="Chisq")
G12_summary <- Anova(G12_mod, type="III", test="Chisq")
G13_summary <- Anova(G13_mod, type="III", test="Chisq")
G14_summary <- Anova(G14_mod, type="III", test="Chisq")
G15_summary <- Anova(G15_mod, type="III", test="Chisq")
G16_summary <- Anova(G16_mod, type="III", test="Chisq")
G18_summary <- Anova(G18_mod, type="III", test="Chisq")
G19_summary <- Anova(G19_mod, type="III", test="Chisq")
G20_summary <- Anova(G20_mod, type="III", test="Chisq")
G21_summary <- Anova(G21_mod, type="III", test="Chisq")
G22_summary <- Anova(G22_mod, type="III", test="Chisq")
G26_summary <- Anova(G26_mod, type="III", test="Chisq")
G29_summary <- Anova(G29_mod, type="III", test="Chisq")
G30_summary <- Anova(G30_mod, type="III", test="Chisq")
G32_summary <- Anova(G32_mod, type="III", test="Chisq")



G1_p_value <- G1_summary[2,4]
G2_p_value <- G2_summary[2,4]
G3_p_value <- G3_summary[2,4]
G4_p_value <- G4_summary[2,4]
G5_p_value <- G5_summary[2,4]
G6_p_value <- G6_summary[2,4]
G7_p_value <- G7_summary[2,4]
G8_p_value <- G8_summary[2,4]
G9_p_value <- G9_summary[2,4]
G10_p_value <- G10_summary[2,4]
G11_p_value <- G11_summary[2,4]
G12_p_value <- G12_summary[2,4]
G13_p_value <- G13_summary[2,4]
G14_p_value <- G14_summary[2,4]
G15_p_value <- G15_summary[2,4]
G16_p_value <- G16_summary[2,4]
G18_p_value <- G18_summary[2,4]
G19_p_value <- G19_summary[2,4]
G20_p_value <- G20_summary[2,4]
G21_p_value <- G21_summary[2,4]
G22_p_value <- G22_summary[2,4]
G26_p_value <- G26_summary[2,4]
G29_p_value <- G29_summary[2,4]
G30_p_value <- G30_summary[2,4]
G32_p_value <- G32_summary[2,4]





















merged_phy_counts <- merge_phyloseq(out1, out2, out3)


#Count to percent
merged_phy_percent <- phyloseq_standardize_otu_abundance(merged_phy_counts, method = "total")




# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent) {
  sd <- sample_data(merged_phy_percent)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent) {
  OTU <- otu_table(merged_phy_percent)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


meta <- pssd2veg(merged_phy_percent)
df <- psotu2veg(merged_phy_percent)


meta <- as.data.frame(meta)
df <- as.data.frame(df)

nrow(df)
nrow(meta)

colnames(df)


All_abund <- as.data.frame(cbind(df$G25, meta))
colnames(All_abund)[1] <- "G25"






G25_plot <- ggplot(All_abund ,aes(x=Order, y=G25, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x=element_blank())
  

G10_plot <- ggplot(All_abund ,aes(x=Order, y=df$G10, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())

G8_plot <- ggplot(All_abund ,aes(x=Order, y=df$G8, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x=element_blank())

G7_plot <- ggplot(All_abund ,aes(x=Order, y=df$G7, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())



G1_plot <- ggplot(All_abund ,aes(x=Order, y=df$G1, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())

G3_plot <- ggplot(All_abund ,aes(x=Order, y=df$G3, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())

G20_plot <- ggplot(All_abund ,aes(x=Order, y=df$G10, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())

G11_plot <- ggplot(All_abund ,aes(x=Order, y=df$G11, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())

G2_plot <- ggplot(All_abund ,aes(x=Order, y=df$G2, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())

G19_plot <- ggplot(All_abund ,aes(x=Order, y=df$G19, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x=element_blank())





library(ggpubr)
IL1_iso_barplot <- ggarrange(G10_plot, G20_plot, G1_plot, G3_plot, G7_plot, G8_plot, G5_plot, labels = c("Bacteroides_uniformis_G10", "Lactobacillus_animalis_G20", "Muribaculum_intestinale_G1", "Lactobacillus_murinus_G3", "Duncaniella_dubosii_G7", "Duncaniella_sp_B8_G8", "Bacteroides_thetaiotaomicron_G5"),
                         common.legend = TRUE, legend = "right",
                         ncol = 2, nrow = 4, align="hv")

IL1_anti6_barplot <- ggarrange(G2_plot, G11_plot, G19_plot, labels = c("Akkermansia_muciniphila_G2", "Bifidobacterium_pseudolongum_G11", "Bacteroides_xylanisolvens_G19"),
                         common.legend = TRUE, legend = "right",
                         ncol = 1, nrow = 3, align="hv")








plot_data <- G25_plot$data

write.csv(plot_data ,file="G25_plot_data.csv")






gpvag_phy <- subset_samples(gpt2, Type == 'vaginal')
gprect_phy <- subset_samples(gpt2, Type == 'rectal')
gppen_phy <- subset_samples(gpt2, Type == 'penile')
control_phy<-subset_samples(gpt2, Type=='Control')

meta_gp$Gp_dpc<-paste(meta_gp$Gp,meta_gp$DPC,sep="_")

df_vag<-phyloseq_to_df(gpvag_phy, addtax = T, addtot = F, addmaxrank = F,
                       sorting = "abundance")
df_vag2<-df_vag[,-c(1,2,3,4,5,6)]
rownames(df_vag2)<-df_vag2$Rank6
df_vag2<-df_vag2[,-c(1)]
df_vag2<-as.data.frame(t(df_vag2))
df_vag2$samp<-rownames(df_vag2)

df_vag_long <- df_vag2 %>%
  pivot_longer(!samp, names_to = "taxa", values_to = "abund")

df_vag_long2<- merge(df_vag_long, meta_gp,
                     by.x="samp", by.y = 'row.names', all = FALSE)

p <-ggplot(df_vag_long2, aes(x = Gp_dpc)) + geom_bar(aes(weight=abund, fill = taxa))+
  facet_grid(.~ Gp, scale="free")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p












###Alpha diversity


merged_phy_counts <- merge_phyloseq(out1, out2, out3)
#merged_phy_counts <- merge_phyloseq(out1)
#merged_phy_counts <- merge_phyloseq(out2)
#merged_phy_counts <- merge_phyloseq(out3)


#merged_phy_counts <- phyloseq_standardize_otu_abundance(merged_phy_counts, method = "total")




meta_phy_out_All <- pssd2veg(merged_phy_counts)
shared_phy_out_All <- psotu2veg(merged_phy_counts)

#View(meta_phy_out_All)

SI_phy <- subset_samples(merged_phy_counts, Order=="A_SI")
Ce_phy <- subset_samples(merged_phy_counts, Order=="B_Ce")
Co_phy <- subset_samples(merged_phy_counts, Order=="C_Co")



SI_alpha_plot <- plot_richness(SI_phy, x="Order", "Order3", measures=c("Chao1", "Shannon", "InvSimpson"))
Ce_alpha_plot <- plot_richness(Ce_phy, x="Order", "Order3", measures=c("Chao1", "Shannon", "InvSimpson"))
Co_alpha_plot <- plot_richness(Co_phy, x="Order", "Order3", measures=c("Chao1", "Shannon", "InvSimpson"))



SI_alpha_plot_data <- ggplot_build(SI_alpha_plot)$plot$data
Ce_alpha_plot_data <- ggplot_build(Ce_alpha_plot)$plot$data
Co_alpha_plot_data <- ggplot_build(Co_alpha_plot)$plot$data



#Normality

SI_Chao1 <- SI_alpha_plot_data[SI_alpha_plot_data$Test=="PBS_Iso" & SI_alpha_plot_data$variable=="Chao1",]
shapiro.test(SI_Chao1$value)

SI_Shannon <- SI_alpha_plot_data[SI_alpha_plot_data$Test=="PBS_Iso" & SI_alpha_plot_data$variable=="Shannon",]
shapiro.test(SI_Shannon$value)

SI_InvSimpson <- SI_alpha_plot_data[SI_alpha_plot_data$Test=="PBS_Iso" & SI_alpha_plot_data$variable=="InvSimpson",]
shapiro.test(SI_InvSimpson$value)



Ce_Chao1 <- Ce_alpha_plot_data[Ce_alpha_plot_data$Test=="PBS_Iso" & Ce_alpha_plot_data$variable=="Chao1",]
shapiro.test(Ce_Chao1$value)

Ce_Shannon <- Ce_alpha_plot_data[Ce_alpha_plot_data$Test=="PBS_Iso" & Ce_alpha_plot_data$variable=="Shannon",]
shapiro.test(Ce_Shannon$value)


Ce_InvSimpson <- Ce_alpha_plot_data[Ce_alpha_plot_data$Test=="PBS_Iso" & Ce_alpha_plot_data$variable=="InvSimpson",]
shapiro.test(Ce_InvSimpson$value)



Co_Chao1 <- Co_alpha_plot_data[Co_alpha_plot_data$Test=="PBS_Iso" & Co_alpha_plot_data$variable=="Chao1",]
shapiro.test(Co_Chao1$value)

Co_Shannon <- Co_alpha_plot_data[Co_alpha_plot_data$Test=="PBS_Iso" & Co_alpha_plot_data$variable=="Shannon",]
shapiro.test(Co_Shannon$value)


Co_InvSimpson <- Co_alpha_plot_data[Co_alpha_plot_data$Test=="PBS_Iso" & Co_alpha_plot_data$variable=="InvSimpson",]
shapiro.test(Co_InvSimpson$value)






SI_alpha_AB_Chao1 <- SI_alpha_plot_data[SI_alpha_plot_data$Order2!="C" & SI_alpha_plot_data$Order2!="D" & SI_alpha_plot_data$variable=="Chao1",]
SI_alpha_AB_Shannon <- SI_alpha_plot_data[SI_alpha_plot_data$Order2!="C" & SI_alpha_plot_data$Order2!="D" & SI_alpha_plot_data$variable=="Shannon",]
SI_alpha_AB_InvSimpson <- SI_alpha_plot_data[SI_alpha_plot_data$Order2!="C" & SI_alpha_plot_data$Order2!="D" & SI_alpha_plot_data$variable=="InvSimpson",]

Ce_alpha_AB_Chao1 <- Ce_alpha_plot_data[Ce_alpha_plot_data$Order2!="C" & Ce_alpha_plot_data$Order2!="D" & Ce_alpha_plot_data$variable=="Chao1",]
Ce_alpha_AB_Shannon <- Ce_alpha_plot_data[Ce_alpha_plot_data$Order2!="C" & Ce_alpha_plot_data$Order2!="D" & Ce_alpha_plot_data$variable=="Shannon",]
Ce_alpha_AB_InvSimpson <- Ce_alpha_plot_data[Ce_alpha_plot_data$Order2!="C" & Ce_alpha_plot_data$Order2!="D" & Ce_alpha_plot_data$variable=="InvSimpson",]

Co_alpha_AB_Chao1 <- Co_alpha_plot_data[Co_alpha_plot_data$Order2!="C" & Co_alpha_plot_data$Order2!="D" & Co_alpha_plot_data$variable=="Chao1",]
Co_alpha_AB_Shannon <- Co_alpha_plot_data[Co_alpha_plot_data$Order2!="C" & Co_alpha_plot_data$Order2!="D" & Co_alpha_plot_data$variable=="Shannon",]
Co_alpha_AB_InvSimpson <- Co_alpha_plot_data[Co_alpha_plot_data$Order2!="C" & Co_alpha_plot_data$Order2!="D" & Co_alpha_plot_data$variable=="InvSimpson",]



SI_alpha_CD_Chao1 <- SI_alpha_plot_data[SI_alpha_plot_data$Order2!="A" & SI_alpha_plot_data$Order2!="B" & SI_alpha_plot_data$variable=="Chao1",]
SI_alpha_CD_Shannon <- SI_alpha_plot_data[SI_alpha_plot_data$Order2!="A" & SI_alpha_plot_data$Order2!="B" & SI_alpha_plot_data$variable=="Shannon",]
SI_alpha_CD_InvSimpson <- SI_alpha_plot_data[SI_alpha_plot_data$Order2!="A" & SI_alpha_plot_data$Order2!="B" & SI_alpha_plot_data$variable=="InvSimpson",]

Ce_alpha_CD_Chao1 <- Ce_alpha_plot_data[Ce_alpha_plot_data$Order2!="A" & Ce_alpha_plot_data$Order2!="B" & Ce_alpha_plot_data$variable=="Chao1",]
Ce_alpha_CD_Shannon <- Ce_alpha_plot_data[Ce_alpha_plot_data$Order2!="A" & Ce_alpha_plot_data$Order2!="B" & Ce_alpha_plot_data$variable=="Shannon",]
Ce_alpha_CD_InvSimpson <- Ce_alpha_plot_data[Ce_alpha_plot_data$Order2!="A" & Ce_alpha_plot_data$Order2!="B" & Ce_alpha_plot_data$variable=="InvSimpson",]

Co_alpha_CD_Chao1 <- Co_alpha_plot_data[Co_alpha_plot_data$Order2!="A" & Co_alpha_plot_data$Order2!="B" & Co_alpha_plot_data$variable=="Chao1",]
Co_alpha_CD_Shannon <- Co_alpha_plot_data[Co_alpha_plot_data$Order2!="A" & Co_alpha_plot_data$Order2!="B" & Co_alpha_plot_data$variable=="Shannon",]
Co_alpha_CD_InvSimpson <- Co_alpha_plot_data[Co_alpha_plot_data$Order2!="A" & Co_alpha_plot_data$Order2!="B" & Co_alpha_plot_data$variable=="InvSimpson",]




mod1 <- wilcox.test(value ~ Test, data=SI_alpha_AB_Chao1, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)

str(mod1)


#Perform the Mann-Whitney U test
#SI_AB
wilcox.test(value ~ Test, data=SI_alpha_AB_Chao1, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(value ~ Test, data=SI_alpha_AB_Shannon, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(value ~ Test, data=SI_alpha_AB_InvSimpson, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)

#Ce_AB
wilcox.test(value ~ Test, data=Ce_alpha_AB_Chao1, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(value ~ Test, data=Ce_alpha_AB_Shannon, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(value ~ Test, data=Ce_alpha_AB_InvSimpson, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)

#Co_AB
wilcox.test(value ~ Test, data=Co_alpha_AB_Chao1, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(value ~ Test, data=Co_alpha_AB_Shannon, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(value ~ Test, data=Co_alpha_AB_InvSimpson, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)


#SI_CD
wilcox.test(value ~ Test, data=SI_alpha_CD_Chao1, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(value ~ Test, data=SI_alpha_CD_Shannon, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(value ~ Test, data=SI_alpha_CD_InvSimpson, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)

#Ce_CD
wilcox.test(value ~ Test, data=Ce_alpha_CD_Chao1, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(value ~ Test, data=Ce_alpha_CD_Shannon, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(value ~ Test, data=Ce_alpha_CD_InvSimpson, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)

#Co_CD
wilcox.test(value ~ Test, data=Co_alpha_CD_Chao1, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(value ~ Test, data=Co_alpha_CD_Shannon, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(value ~ Test, data=Co_alpha_CD_InvSimpson, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)







friedman.test(formula = SI_alpha_AB_Chao1$value ~ SI_alpha_AB_Chao1$Test | SI_alpha_AB_Chao1$Weight, data = SI_alpha_AB_Chao1)

friedman.test(formula = SI_alpha_AB_Chao1$value ~ SI_alpha_AB_Chao1$Test | SI_alpha_AB_Chao1$Weight, data = SI_alpha_AB_Chao1)

friedman.test(formula = disease ~ locations | plant_var, data = df_long)




library(stats)
friedman.test(Ce_Chao1$value ~ Ce_Chao1$Test + Ce_Chao1$Weight )

friedman.test(formula = Test ~ Weight, data = Ce_Chao1)


friedman.test(Ce_Chao1$value = Ce_Chao1$Test +Weight, groups = df_long$locations, blocks = df_long$plant_var)









Co_Chao1 <- Co_alpha_plot_data[Co_alpha_plot_data$Test=="PBS_Iso" & Co_alpha_plot_data$variable=="Chao1",]
shapiro.test(Co_Chao1$value)

Co_Shannon <- Co_alpha_plot_data[Co_alpha_plot_data$Test=="PBS_Iso" & Co_alpha_plot_data$variable=="Shannon",]
shapiro.test(Co_Shannon$value)

#Not normal PBS_Iso
Co_InvSimpson <- Co_alpha_plot_data[Co_alpha_plot_data$Test=="PBS_Iso" & Co_alpha_plot_data$variable=="InvSimpson",]
shapiro.test(Co_InvSimpson$value)



















Co_Chao1 <- Co_alpha_plot_data[Co_alpha_plot_data$Test=="PBS_Iso" & Co_alpha_plot_data$variable=="Chao1",]
shapiro.test(Co_Chao1$value)

Co_Shannon <- Co_alpha_plot_data[Co_alpha_plot_data$Test=="PBS_Iso" & Co_alpha_plot_data$variable=="Shannon",]
shapiro.test(Co_Shannon$value)

Co_InvSimpson <- Co_alpha_plot_data[Co_alpha_plot_data$Test=="PBS_Iso" & Co_alpha_plot_data$variable=="InvSimpson",]
shapiro.test(Co_InvSimpson$value)






SI_alpha_plot <- ggplot(SI_alpha_plot_data ,aes(x=Order2, y=value, color=Order2, fill=Order3, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  facet_wrap(~variable, scales = "free_y", nrow = 1)
  

Ce_alpha_plot <- ggplot(Ce_alpha_plot_data ,aes(x=Order2, y=value, color=Order2, fill=Order3, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  facet_wrap(~variable, scales = "free_y", nrow = 1)


Co_alpha_plot <- ggplot(Co_alpha_plot_data ,aes(x=Order2, y=value, color=Order2, fill=Order3, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  facet_wrap(~variable, scales = "free_y", nrow = 1)



library(ggpubr)
merged_alpha_plots <- ggarrange(SI_alpha_plot, Ce_alpha_plot, Co_alpha_plot, labels = c("SI", "Ce", "Co"),
                        common.legend = TRUE, legend = "right",
                        ncol = 1, nrow = 3, align="hv")





#View(Ce_Co_alpha_plot_data)


Ce_Co_alpha_plot_data <- rbind(Ce_alpha_plot_data)

Ce_Co_alpha_plot_data  <- Ce_Co_alpha_plot_data[Ce_Co_alpha_plot_data$Order2!='C' & Ce_Co_alpha_plot_data$Order2!='D',]
#Ce_Co_alpha_plot_data  <- Ce_Co_alpha_plot_data[Ce_Co_alpha_plot_data$Order=='B_Ce' ,]



Ce_Co_Chao1_data <- (Ce_Co_alpha_plot_data[Ce_Co_alpha_plot_data$variable=='Chao1' , ])
Ce_Co_Shannon_data <- (Ce_Co_alpha_plot_data[Ce_Co_alpha_plot_data$variable=='Shannon' , ])
Ce_Co_InvSimpson_data <- (Ce_Co_alpha_plot_data[Ce_Co_alpha_plot_data$variable=='InvSimpson' , ])




library(vegan)
packageVersion("vegan")


library(lme4)
library(car)




Chao1_mod <- lmer(value ~ Order2 + Weight + (1|dam_ID), data = Ce_Co_Chao1_data, REML=TRUE)
Anova(Chao1_mod, type="II")

Chao1_mod <- aov(value ~ Order2 + Weight, data = Ce_Co_Chao1_data)
Anova(Chao1_mod, type="II")

Chao1_mod <- aov(value ~ Order2, data = Ce_Co_Chao1_data)
Anova(Chao1_mod, type="II")




Shannon_mod <- lmer(value ~ Order2 + Weight + (1|dam_ID), data = Ce_Co_Shannon_data, REML=TRUE)
Anova(Shannon_mod, type="II")

Shannon_mod <- aov(value ~ Order2 + Weight, data = Ce_Co_Shannon_data)
Anova(Shannon_mod, type="II")

Shannon_mod <- aov(value ~ Order2, data = Ce_Co_Shannon_data)
Anova(Shannon_mod, type="II")




Shannon_mod <- lmer(value ~ Order2 + Weight + (1|dam_ID), data = Ce_Co_InvSimpson_data, REML=TRUE)
Anova(InvSimpson_mod, type="II")

InvSimpson_mod <- aov(value ~ Order2 + Weight, data = Ce_Co_InvSimpson_data)
Anova(InvSimpson_mod, type="II")

InvSimpson_mod <- aov(value ~ Order2, data = Ce_Co_InvSimpson_data)
Anova(InvSimpson_mod, type="II")






Chao1_mod <- lmer(value ~ Order2 + Order + Weight + (1|dam_ID) + (1|pup_ID), data = Ce_Co_Chao1_data, REML=TRUE)
Anova(Chao1_mod, type="II")

Shannon_mod <- lmer(value ~ Order2 + Order + Weight + (1|dam_ID) + (1|pup_ID), data = Ce_Co_Shannon_data, REML=TRUE)
Anova(Shannon_mod, type="II")

InvSimpson_mod <- lmer(value ~ Order2 + Order + Weight + (1|dam_ID) + (1|pup_ID), data = Ce_Co_InvSimpson_data, REML=TRUE)
Anova(InvSimpson_mod, type="II")









Chao1_mod <- aov(value ~ Order2 + Order + Weight, data = Ce_Co_Chao1_data)
Anova(Chao1_mod, type="II")

Shannon_mod <- aov(value ~ Order2 + Order + Weight, data = Ce_Co_Shannon_data)
Anova(Shannon_mod, type="II")

InvSimpson_mod <- aov(value ~ Order2 + Order + Weight, data = Ce_Co_InvSimpson_data)
Anova(InvSimpson_mod, type="II")







Chao1_mod <- aov(value ~ Order2 + Weight, data = Ce_Co_Chao1_data)
Anova(Chao1_mod, type="II")

Shannon_mod <- aov(value ~ Order2 + Weight, data = Ce_Co_Shannon_data)
Anova(Shannon_mod, type="II")

InvSimpson_mod <- aov(value ~ Order2 + Weight, data = Ce_Co_InvSimpson_data)
Anova(InvSimpson_mod, type="II")











library(lmtest)

mod1 <- lmer(value ~ Test + Type + Weight + (1|dam_ID) + (1|pup_ID), data = Ce_Co_InvSimpson_data, REML=FALSE)
Anova(InvSimpson_mod, type="II")

mod2 <- lmer(value ~ Test + Type + Weight + (1|dam_ID) + (1|pup_ID), data = Ce_Co_InvSimpson_data, REML=TRUE)
Anova(InvSimpson_mod, type="II")

lrtest(mod1, mod2)






























###Heatmaps and barplots

























###Weight








merged_phy_percent <- merge_phyloseq(out1, out2, out3)


# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent) {
  sd <- sample_data(merged_phy_percent)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent) {
  OTU <- otu_table(merged_phy_percent)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


meta_phy_out_All <- pssd2veg(merged_phy_percent)
shared_phy_out_All <- psotu2veg(merged_phy_percent)

shared_phy_out_All <- as.data.frame(shared_phy_out_All)

colnames(shared_phy_out_All)


All_abund <- as.data.frame(cbind(shared_phy_out_All$G2, shared_phy_out_All$G7, shared_phy_out_All$G8, meta_phy_out_All))
colnames(All_abund)[1] <- "G2"
colnames(All_abund)[2] <- "G7"
colnames(All_abund)[3] <- "G8"




eq <- function(x,y) {
  m <- lm(y ~ x)
  as.character(
    as.expression(
      substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                 list(a = format(coef(m)[1], digits = 4),
                      b = format(coef(m)[2], digits = 4),
                      r2 = format(summary(m)$r.squared, digits = 3)))
    )
  )
}





G2_p <- ggplot(All_abund, aes(Weight, G2)) +
        geom_point() +
        geom_smooth(method = "lm") +
        stat_regline_equation(label.x = 8, label.y = 750000, aes(label = ..rr.label..), size=8)
  

G7_p <- ggplot(All_abund, aes(Weight, G7)) +
        geom_point() +
        geom_smooth(method = "lm") +
        stat_regline_equation(label.x = 8, label.y = 110000, aes(label = ..rr.label..), size=8)
      


G8_p <- ggplot(All_abund, aes(Weight, G8)) +
        geom_point() +
        geom_smooth(method = "lm") +
        stat_regline_equation(label.x = 8, label.y = 125000, aes(label = ..rr.label..), size=8)






library(ggpubr)
All_weight_plot <- ggarrange(G2_p, G7_p, G8_p, labels = c("Akkermansia_G2", 	"Duncaniella_G7", "Duncaniella_G8"),
                        common.legend = TRUE, legend = "right",
                        ncol = 3, nrow = 1, align="hv")



weight_mod1 <- lmer(G8 ~ Weight + Test + Type + (1|dam_ID) + (1|pup_ID), data = All_abund, REML=TRUE)
Anova(weight_mod1, type="II")



Co_alpha_plot <- ggplot(All_abund, aes(x=Order2, y=Weight, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic()

























All_plot_data <- G2_p$data

#rite.csv(All_plot_data,file="All_plot_data.csv")




meta <- pssd2veg(merged_phy_counts)
df <- psotu2veg(merged_phy_counts)


meta <- as.data.frame(meta)
df <- as.data.frame(df)


SI_alpha_plot <- ggplot(SI_alpha_plot_data ,aes(x=Order2, y=weight, color=Order2, fill=Order2, shape=Order2)) +
  geom_boxplot(outlier.shape = NA, alpha=0) +
  geom_point(position=position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.5, stroke=1.25) +
  scale_shape_manual(values = c(21, 21, 1, 1)) +
  scale_color_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  scale_fill_manual(values=c("#0000FF", "#FF0000", "#0000FF", "#FF0000")) +
  theme_classic() +
  facet_wrap(~variable, scales = "free_y", nrow = 1)


colnames(df)







All_weight_df <- data.frame(All_plot_data$G2, All_plot_data$G7, All_plot_data$G8, All_plot_data$Weight)

colnames(All_weight_df)[1] <- "G2"
colnames(All_weight_df)[2] <- "G7"
colnames(All_weight_df)[3] <- "G8"
colnames(All_weight_df)[4] <- "Weight"



G2_plot <- ggplot(All_weight_df, aes(Weight, G2)) +
  geom_point() +
  geom_smooth(method = "lm")

mod1 <- aov(data=G2_df, G2_df$G2 ~ G2_df$Weight)
summary(mod1)



G7_plot <- ggplot(All_weight_df, aes(Weight, G7)) +
  geom_point() +
  geom_smooth(method = "lm")

mod1 <- aov(data=G2_df, G2_df$G2 ~ G2_df$Weight)
summary(mod1)



G8_plot <- ggplot(All_weight_df, aes(Weight, G8)) +
  geom_point() +
  geom_smooth(method = "lm")

mod1 <- aov(data=G2_df, G2_df$G2 ~ G2_df$Weight)
summary(mod1)


All_Weight <- ggarrange(G2_plot, G7_plot, G8_plot, labels = c("Akkermansia_G2", "Duncaniella_G7", 	"Duncaniella_G8"),
                        common.legend = TRUE, legend = "right",
                        ncol = 3, nrow = 1, align="hv")










All_weight_df <- data.frame(All_plot_data$G2, All_plot_data$G7, All_plot_data$G8, All_plot_data$Weight)

colnames(All_weight_df)[1] <- "G2"
colnames(All_weight_df)[2] <- "G7"
colnames(All_weight_df)[3] <- "G8"
colnames(All_weight_df)[4] <- "Weight"



G2_plot <- ggplot(All_weight_df, aes(Weight, G2)) +
  geom_point() +
  geom_smooth(method = "lm")

mod1 <- aov(data=G2_df, G2_df$G2 ~ G2_df$Weight)
summary(mod1)



G7_plot <- ggplot(All_weight_df, aes(Weight, G7)) +
  geom_point() +
  geom_smooth(method = "lm")

mod1 <- aov(data=G2_df, G2_df$G2 ~ G2_df$Weight)
summary(mod1)



G8_plot <- ggplot(All_weight_df, aes(Weight, G8)) +
  geom_point() +
  geom_smooth(method = "lm")

mod1 <- aov(data=G2_df, G2_df$G2 ~ G2_df$Weight)
summary(mod1)


All_Weight <- ggarrange(G2_plot, G7_plot, G8_plot, labels = c("Akkermansia_G2", "Duncaniella_G7", 	"Duncaniella_G8"),
                        common.legend = TRUE, legend = "right",
                        ncol = 1, nrow = 3, align="hv")












































































###Heatmaps and barplots




















#devtools::install_github("vmikk/metagMisc")
library(phyloseq)
library(metagMisc)
library(vegan)
library(ggplot2)






#B vs D



meta_ob1  <- read.table("/Users/zebra/Desktop/GF/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("/Users/zebra/Desktop/GF/taxa.txt"))


meta_ob_BD <- meta_ob1[meta_ob1$Origin=='Neonatal' & meta_ob1$Type!="Small_gut" & meta_ob1$Group!="A" & meta_ob1$Group!="C",]

nrow(meta_ob_BD)



#Get sample names
GroupSamples<-list(rownames(meta_ob_BD))

#Subset based on sample names
fullMothur<-read.csv("shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtu))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob_BD), tax_table(tax_ob1))


#Rarefy
#minimum <- min(sample_sums(phyloseqobj.f))
out1 <- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1278327, rngseed = 1, replace = FALSE)
#out1<- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1000, rngseed = 1, replace = FALSE)

#Count to percent
phy1 <- phyloseq_standardize_otu_abundance(out1, method = "total")


merged_phy_percent_for_plots <- phyloseq_standardize_otu_abundance(phy1, method = "total")



percent_Cecum_BD <- subset_samples(merged_phy_percent_for_plots , Type=="Cecum")
Cecum_pruned_BD <- prune_taxa(names(sort(taxa_sums(percent_Cecum_BD),TRUE)[1:30]), percent_Cecum_BD)
Ce_plot_BD <- plot_heatmap(Cecum_pruned_BD , method = "PCoA", distance = "bray", "Injection", "Rank1",
                           trans = NULL, sample.order="Injection", low="azure1", high="dodgerblue3", na.value="white",
                           title="Neonatal Cecum B vs D") +
  theme(axis.text.x=element_text(size=4)) +
  theme(axis.text.y=element_text(size=8)) +
  theme(axis.title.x=element_blank())


percent_Colon_BD <- subset_samples(merged_phy_percent_for_plots , Type=="Colon")
Colon_pruned_BD <- prune_taxa(names(sort(taxa_sums(percent_Cecum_BD),TRUE)[1:30]), percent_Colon_BD)
Co_plot_BD <- plot_heatmap(Colon_pruned_BD , method = "PCoA", distance = "bray", "Injection", "Rank1",
                           trans = NULL, sample.order="Injection", low="azure1", high="dodgerblue3", na.value="white",
                           title="Neonatal Colon B vs D") +
  theme(axis.text.x=element_text(size=4)) +
  theme(axis.text.y=element_text(size=8)) +
  theme(axis.title.x=element_blank())



Cecum_bar_BD <- plot_bar(Cecum_pruned_BD , fill="Rank2") + facet_wrap(~Test, scales = "free_x", nrow = 1) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())

Colon_bar_BD <- plot_bar(Colon_pruned_BD , fill="Rank2") + facet_wrap(~Test, scales = "free_x", nrow = 1) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())




#A vs C




meta_ob1  <- read.table("/Users/zebra/Desktop/GF/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("/Users/zebra/Desktop/GF/taxa.txt"))


meta_ob_AC <- meta_ob1[meta_ob1$Origin=='Neonatal' & meta_ob1$Type!="Small_gut" & meta_ob1$Group!="B" & meta_ob1$Group!="D",]

nrow(meta_ob_AC)



#Get sample names
GroupSamples<-list(rownames(meta_ob_AC))

#Subset based on sample names
fullMothur<-read.csv("shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtu))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob_AC), tax_table(tax_ob1))


#Rarefy
#minimum <- min(sample_sums(phyloseqobj.f))
out1 <- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1278327, rngseed = 1, replace = FALSE)
#out1<- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1000, rngseed = 1, replace = FALSE)

#Count to percent
phy1 <- phyloseq_standardize_otu_abundance(out1, method = "total")


merged_phy_percent_for_plots <- phyloseq_standardize_otu_abundance(phy1, method = "total")



percent_Cecum_AC <- subset_samples(merged_phy_percent_for_plots , Type=="Cecum")
Cecum_pruned_AC <- prune_taxa(names(sort(taxa_sums(percent_Cecum_AC),TRUE)[1:30]), percent_Cecum_AC)
Ce_plot_AC <- plot_heatmap(Cecum_pruned_AC , method = "PCoA", distance = "bray", "Injection", "Rank1",
                           trans = NULL, sample.order="Injection", low="azure1", high="dodgerblue3", na.value="white",
                           title="Neonatal Cecum A vs C") +
  theme(axis.text.x=element_text(size=4)) +
  theme(axis.text.y=element_text(size=8)) +
  theme(axis.title.x=element_blank())


percent_Colon_AC <- subset_samples(merged_phy_percent_for_plots , Type=="Colon")
Colon_pruned_AC <- prune_taxa(names(sort(taxa_sums(percent_Cecum_AC),TRUE)[1:30]), percent_Colon_AC)
Co_plot_AC <- plot_heatmap(Colon_pruned_AC , method = "PCoA", distance = "bray", "Injection", "Rank1",
                           trans = NULL, sample.order="Injection", low="azure1", high="dodgerblue3", na.value="white",
                           title="Neonatal Colon A vs C") +
  theme(axis.text.x=element_text(size=4)) +
  theme(axis.text.y=element_text(size=8)) +
  theme(axis.title.x=element_blank())


Cecum_bar_AC <- plot_bar(Cecum_pruned_AC , fill="Rank2") + facet_wrap(~Test, scales = "free_x", nrow = 1) +
theme(axis.title.x=element_blank(), axis.text.x=element_blank())

Colon_bar_AC <- plot_bar(Colon_pruned_AC , fill="Rank2") + facet_wrap(~Test, scales = "free_x", nrow = 1) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())






library(ggpubr)
Heatmap_plot <- ggarrange(Ce_plot_BD, Co_plot_BD, Ce_plot_AC, Co_plot_AC, labels = c("Ce_BD", "Co_BD", "Ce_AC", "Co_AC"),
                       common.legend = TRUE, legend = "right",
                       ncol = 2, nrow = 2, align="hv")


bar_plots <- ggarrange(Cecum_bar_BD, Colon_bar_BD, Cecum_bar_AC, Colon_bar_AC, labels = c("Cecum", "Colon", "Cecum", "Colon"),
                       common.legend = TRUE, legend = "right",
                       ncol = 2, nrow = 2, align="hv")



























































































#Co-occurrences of strains



setwd("/Users/zebra/Desktop/GF/")

meta_ob1  <- read.table("/Users/zebra/Desktop/GF/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("/Users/zebra/Desktop/GF/taxa.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Origin!='2Neonatal' ,]

nrow(meta_ob1)


#Get sample names
GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtu))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))


#Rarefy
#minimum <- min(sample_sums(phyloseqobj.f))
#out1 <- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1278327, rngseed = 1, replace = FALSE)

#Count to percent
phy1 <- phyloseq_standardize_otu_abundance(phyloseqobj.f , method = "total")


merged_phy_percent <- phyloseq_standardize_otu_abundance(phy1, method = "total")









taxa_df <- subset_taxa(merged_phy_percent2, Rank2=="Duncaniella")
pruned <- prune_samples(sample_sums(taxa_df) > 0, taxa_df)


plot_heatmap(pruned, method = "PCoA", distance = "bray", "Sample_ID", "Rank3", trans = NULL, sample.order="Group", low="azure1", high="dodgerblue3", na.value="white")

































###Neonatal Cecum PCoA B_D


setwd("/Users/zebra/Desktop/GF/")

meta_ob1  <- read.table("/Users/zebra/Desktop/GF/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("/Users/zebra/Desktop/GF/taxa.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Type=="Cecum" & meta_ob1$Origin=='Neonatal' ,]

nrow(meta_ob1)


#Get sample names
GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtu))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))


#Rarefy
#minimum <- min(sample_sums(phyloseqobj.f))
out1 <- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1278327, rngseed = 1, replace = FALSE)
#out1<- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1000, rngseed = 1, replace = FALSE)

#Count to percent
phy1 <- phyloseq_standardize_otu_abundance(out1, method = "total")


merged_phy_percent <- phyloseq_standardize_otu_abundance(phy1, method = "total")



# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent) {
  sd <- sample_data(merged_phy_percent)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent) {
  OTU <- otu_table(merged_phy_percent)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


meta_phy_out <- pssd2veg(merged_phy_percent)
shared_phy_out <- psotu2veg(merged_phy_percent)


#dist_bray <- vegdist(shared_phy_out, method="bray", binary=FALSE)
#adonis2(formula = dist_bray ~ Type + Test, data = meta_phy_out, permutations=10000, strata=meta_phy_out$pup_ID)



dist_bray <- capscale(shared_phy_out~1, distance="bray", binary=FALS, eig=TRUE)

prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)


data.scores = scores(dist_bray, display = "sites", choices = c(1:10))
#summary(data.scores)


df <- data.frame(data.scores)

library(glue)
labx <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby <- c(glue("PCo 2 ({prop_explained [2]}%)"))


library(dplyr)
library(ggrepel)
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = shared_phy_out)

#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)

#figure out which ASVs you want to plot
#ASV_vector <- c("G42", "G201")

#can include code here to select specific ASVs
#WAscores_df <- WAscores_df %>% filter(rownames(.) %in% ASV_vector)

#, shape=meta_phy_out$Type

ggplot(df, aes(x = MDS1, y = MDS2, colour = meta_phy_out$Test, fill=meta_phy_out$Test)) +
  stat_ellipse(geom="polygon", level=0.75, alpha=0.25, show.legend=FALSE) +
  geom_point(mapping = aes(colour = meta_phy_out$Test), size=5, alpha=1.0) +
  scale_color_manual(values=c("#B81A1A", "#283DA8", "#56D43A", "#5A5C59")) +
  scale_fill_manual(values=c("pink", "dodgerblue", "#56D43A", "#5A5C59")) +
  #coord_fixed()+
  theme_classic() +
  labs(x=labx, y=laby) +
  #adds Weighted average score points for ASVs
  geom_point(data = WAscores_df, colour = "black",aes(x = MDS1, y = MDS2), inherit.aes = FALSE) +
  geom_text(data = WAscores_df, inherit.aes = FALSE, aes(x=MDS1, y= MDS2), colour="black", label = rownames(WAscores_df))+
  #xlab("PC1 (22.5%)") +
  #ylab("PC2 (13.8%)") +
  ggtitle("Neonatal Cecum Bacterial Community Structure") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=12, face="bold"), axis.text.y = element_text(color="black"))
#legend.background = element_rect(fill="white", color="black")
#theme(legend.position = c(0.9, 0.4)) +
#annotate("text", x=0.3, y=0.7, label = "F = 1.71, R2 = 0.02, p = 0.6021", size=4) +
#annotate("text", x=0.3, y=0.56, label = "Far/BLK: F = 0.918, R2 = 0.051, p = 0.4614 ", size=4) +
#annotate("text", x=0.3, y=0.42, label = "Near/Far (Mouse): F = 1.015, R2 = 0.4616, p = 0.409", size=4) +
#annotate("text", x=0.3, y=0.28, label = "Near/Far (Site): F = 1.106, R2 = 0.084, p = 0.245", size=4) +
#theme(legend.title=element_blank()) +
#theme(legend.position="bottom", legend.box = "horizontal")





###Neonatal Cecum adonis


merged_phy_percent_A <- subset_samples(merged_phy_percent, Group=="A")
merged_phy_percent_B <- subset_samples(merged_phy_percent, Group=="B")
merged_phy_percent_A_B <- merge_phyloseq(merged_phy_percent_A, merged_phy_percent_B)

# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent_A_B) {
  sd <- sample_data(merged_phy_percent_A_B)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent_A_B) {
  OTU <- otu_table(merged_phy_percent_A_B)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

meta_phy_out <- pssd2veg(merged_phy_percent_A_B)
shared_phy_out <- psotu2veg(merged_phy_percent_A_B)

dist_bray <- vegdist(shared_phy_out, method="bray", binary=FALSE)
adonis2(formula = dist_bray ~ Group, data = meta_phy_out, permutations=10000)




merged_phy_percent_A <- subset_samples(merged_phy_percent, Group=="A")
merged_phy_percent_C <- subset_samples(merged_phy_percent, Group=="C")
merged_phy_percent_A_C <- merge_phyloseq(merged_phy_percent_A, merged_phy_percent_C)

# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent_A_C) {
  sd <- sample_data(merged_phy_percent_A_C)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent_A_C) {
  OTU <- otu_table(merged_phy_percent_A_C)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

meta_phy_out <- pssd2veg(merged_phy_percent_A_C)
shared_phy_out <- psotu2veg(merged_phy_percent_A_C)

dist_bray <- vegdist(shared_phy_out, method="bray", binary=FALSE)
adonis2(formula = dist_bray ~ Group, data = meta_phy_out, permutations=10000)





merged_phy_percent_A <- subset_samples(merged_phy_percent, Group=="A")
merged_phy_percent_D <- subset_samples(merged_phy_percent, Group=="D")
merged_phy_percent_A_D <- merge_phyloseq(merged_phy_percent_A, merged_phy_percent_D)

# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent_A_D) {
  sd <- sample_data(merged_phy_percent_A_D)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent_A_D) {
  OTU <- otu_table(merged_phy_percent_A_D)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

meta_phy_out <- pssd2veg(merged_phy_percent_A_D)
shared_phy_out <- psotu2veg(merged_phy_percent_A_D)

dist_bray <- vegdist(shared_phy_out, method="bray", binary=FALSE)
adonis2(formula = dist_bray ~ Group, data = meta_phy_out, permutations=10000)







merged_phy_percent_B <- subset_samples(merged_phy_percent, Group=="B")
merged_phy_percent_C <- subset_samples(merged_phy_percent, Group=="C")
merged_phy_percent_B_C <- merge_phyloseq(merged_phy_percent_B, merged_phy_percent_C)

# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent_B_C) {
  sd <- sample_data(merged_phy_percent_B_C)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent_B_C) {
  OTU <- otu_table(merged_phy_percent_B_C)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

meta_phy_out <- pssd2veg(merged_phy_percent_B_C)
shared_phy_out <- psotu2veg(merged_phy_percent_B_C)

dist_bray <- vegdist(shared_phy_out, method="bray", binary=FALSE)
adonis2(formula = dist_bray ~ Group, data = meta_phy_out, permutations=10000)





merged_phy_percent_B <- subset_samples(merged_phy_percent, Group=="B")
merged_phy_percent_D <- subset_samples(merged_phy_percent, Group=="D")
merged_phy_percent_B_D <- merge_phyloseq(merged_phy_percent_B, merged_phy_percent_D)

# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent_B_D) {
  sd <- sample_data(merged_phy_percent_B_D)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent_B_D) {
  OTU <- otu_table(merged_phy_percent_B_D)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

meta_phy_out <- pssd2veg(merged_phy_percent_B_D)
shared_phy_out <- psotu2veg(merged_phy_percent_B_D)

dist_bray <- vegdist(shared_phy_out, method="bray", binary=FALSE)
adonis2(formula = dist_bray ~ Group, data = meta_phy_out, permutations=10000)





merged_phy_percent_C <- subset_samples(merged_phy_percent, Group=="C")
merged_phy_percent_D <- subset_samples(merged_phy_percent, Group=="D")
merged_phy_percent_C_D <- merge_phyloseq(merged_phy_percent_C, merged_phy_percent_D)

# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent_C_D) {
  sd <- sample_data(merged_phy_percent_C_D)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent_C_D) {
  OTU <- otu_table(merged_phy_percent_C_D)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

meta_phy_out <- pssd2veg(merged_phy_percent_C_D)
shared_phy_out <- psotu2veg(merged_phy_percent_C_D)

dist_bray <- vegdist(shared_phy_out, method="bray", binary=FALSE)
adonis2(formula = dist_bray ~ Group, data = meta_phy_out, permutations=10000)








merged_phy_percent_A_C_prune <- prune_taxa(names(sort(taxa_sums(merged_phy_percent_A_C),TRUE)[1:295]), merged_phy_percent_A_C)
plot_heatmap(merged_phy_percent_A_C_prune, method = "PCoA", distance = "bray", "Test", "Rank1", trans = NULL, sample.order="Test", low="azure1", high="dodgerblue3", na.value="white")


merged_phy_percent_B_D_prune <- prune_taxa(names(sort(taxa_sums(merged_phy_percent_B_D),TRUE)[1:295]), merged_phy_percent_B_D)
plot_heatmap(merged_phy_percent_B_D_prune, method = "PCoA", distance = "bray", "Test", "Rank1", trans = NULL, sample.order="Test", low="azure1", high="dodgerblue3", na.value="white")



plot_bar(merged_phy_percent_A_C_prune, fill="Rank2") + facet_wrap(~Test, scales = "free_x", nrow = 2)

plot_bar(merged_phy_percent_B_D_prune, fill="Rank2") + facet_wrap(~Test, scales = "free_x", nrow = 2)




taxa_df <- subset_taxa(merged_phy_percent, Rank2=="Enterobacter")
pruned <- prune_samples(sample_sums(taxa_df) > 0, taxa_df)


plot_heatmap(pruned, method = "PCoA", distance = "bray", "Sample_ID", "Rank1", trans = NULL, sample.order="Group", low="azure1", high="dodgerblue3", na.value="white")



#devtools::install_github("vmikk/metagMisc")
library(phyloseq)
library(metagMisc)
library(vegan)
library(ggplot2)




###Colon PCoA


setwd("/Users/zebra/Desktop/GF/")

meta_ob1  <- read.table("/Users/zebra/Desktop/GF/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("/Users/zebra/Desktop/GF/taxa.txt"))

meta_ob1 <- meta_ob1[meta_ob1$Type=="Colon" & meta_ob1$Origin=='Neonatal' ,]

nrow(meta_ob1)


#Get sample names
GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtu))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))


#Rarefy
#minimum <- min(sample_sums(phyloseqobj.f))
out1 <- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1278327, rngseed = 1, replace = FALSE)
#out1<- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1000, rngseed = 1, replace = FALSE)

#Count to percent
phy1 <- phyloseq_standardize_otu_abundance(out1, method = "total")


merged_phy_percent <- phyloseq_standardize_otu_abundance(phy1, method = "total")



# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent) {
  sd <- sample_data(merged_phy_percent)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent) {
  OTU <- otu_table(merged_phy_percent)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


meta_phy_out <- pssd2veg(merged_phy_percent)
shared_phy_out <- psotu2veg(merged_phy_percent)


#dist_bray <- vegdist(shared_phy_out, method="bray", binary=FALSE)
#adonis2(formula = dist_bray ~ Type + Test, data = meta_phy_out, permutations=10000, strata=meta_phy_out$pup_ID)



dist_bray <- capscale(shared_phy_out~1, distance="bray", binary=FALS, eig=TRUE)

prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)


data.scores = scores(dist_bray, display = "sites", choices = c(1:10))
#summary(data.scores)


df <- data.frame(data.scores)

library(glue)
labx <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby <- c(glue("PCo 2 ({prop_explained [2]}%)"))


library(dplyr)
library(ggrepel)
#calculates Weighted average scores
WAscores <- vegan::wascores(x = data.scores, w = shared_phy_out)
#converts to dataframe and takes first two dimensions
WAscores_df <- as.data.frame(WAscores) %>% dplyr::select(1:2)

#figure out which ASVs you want to plot
ASV_vector <- c("G42", "G201")
#can include code here to select specific ASVs
WAscores_df <- WAscores_df %>% filter(rownames(.) %in% ASV_vector)

#, shape=meta_phy_out$Type

ggplot(df, aes(x = MDS1, y = MDS2, colour = meta_phy_out$Test, fill=meta_phy_out$Test)) +
  stat_ellipse(geom="polygon", level=0.75, alpha=0.25, show.legend=FALSE) +
  geom_point(mapping = aes(colour = meta_phy_out$Test), size=5, alpha=1.0) +
  scale_color_manual(values=c("#B81A1A", "#283DA8", "#56D43A", "#5A5C59")) +
  scale_fill_manual(values=c("pink", "dodgerblue", "#56D43A", "#5A5C59")) +
  #coord_fixed()+
  theme_classic() +
  labs(x=labx, y=laby) +
  #adds Weighted average score points for ASVs
  geom_point(data = WAscores_df, colour = "black",aes(x = MDS1, y = MDS2), inherit.aes = FALSE) +
  geom_text(data = WAscores_df, inherit.aes = FALSE, aes(x=MDS1, y= MDS2), colour="black", label = rownames(WAscores_df))+
  #xlab("PC1 (22.5%)") +
  #ylab("PC2 (13.8%)") +
  ggtitle("Neonatal Colon Bacterial Community Structure") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=12, face="bold"), axis.text.y = element_text(color="black"))
#legend.background = element_rect(fill="white", color="black")
#theme(legend.position = c(0.9, 0.4)) +
#annotate("text", x=0.3, y=0.7, label = "F = 1.71, R2 = 0.02, p = 0.6021", size=4) +
#annotate("text", x=0.3, y=0.56, label = "Far/BLK: F = 0.918, R2 = 0.051, p = 0.4614 ", size=4) +
#annotate("text", x=0.3, y=0.42, label = "Near/Far (Mouse): F = 1.015, R2 = 0.4616, p = 0.409", size=4) +
#annotate("text", x=0.3, y=0.28, label = "Near/Far (Site): F = 1.106, R2 = 0.084, p = 0.245", size=4) +
#theme(legend.title=element_blank()) +
#theme(legend.position="bottom", legend.box = "horizontal")









###Neonatal Cecum adonis


merged_phy_percent_A <- subset_samples(merged_phy_percent, Group=="A")
merged_phy_percent_B <- subset_samples(merged_phy_percent, Group=="B")
merged_phy_percent_A_B <- merge_phyloseq(merged_phy_percent_A, merged_phy_percent_B)

# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent_A_B) {
  sd <- sample_data(merged_phy_percent_A_B)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent_A_B) {
  OTU <- otu_table(merged_phy_percent_A_B)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

meta_phy_out <- pssd2veg(merged_phy_percent_A_B)
shared_phy_out <- psotu2veg(merged_phy_percent_A_B)

dist_bray <- vegdist(shared_phy_out, method="bray", binary=FALSE)
adonis2(formula = dist_bray ~ Group, data = meta_phy_out, permutations=10000)




merged_phy_percent_A <- subset_samples(merged_phy_percent, Group=="A")
merged_phy_percent_C <- subset_samples(merged_phy_percent, Group=="C")
merged_phy_percent_A_C <- merge_phyloseq(merged_phy_percent_A, merged_phy_percent_C)

# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent_A_C) {
  sd <- sample_data(merged_phy_percent_A_C)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent_A_C) {
  OTU <- otu_table(merged_phy_percent_A_C)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

meta_phy_out <- pssd2veg(merged_phy_percent_A_C)
shared_phy_out <- psotu2veg(merged_phy_percent_A_C)

dist_bray <- vegdist(shared_phy_out, method="bray", binary=FALSE)
adonis2(formula = dist_bray ~ Group, data = meta_phy_out, permutations=10000)





merged_phy_percent_A <- subset_samples(merged_phy_percent, Group=="A")
merged_phy_percent_D <- subset_samples(merged_phy_percent, Group=="D")
merged_phy_percent_A_D <- merge_phyloseq(merged_phy_percent_A, merged_phy_percent_D)

# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent_A_D) {
  sd <- sample_data(merged_phy_percent_A_D)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent_A_D) {
  OTU <- otu_table(merged_phy_percent_A_D)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

meta_phy_out <- pssd2veg(merged_phy_percent_A_D)
shared_phy_out <- psotu2veg(merged_phy_percent_A_D)

dist_bray <- vegdist(shared_phy_out, method="bray", binary=FALSE)
adonis2(formula = dist_bray ~ Group, data = meta_phy_out, permutations=10000)







merged_phy_percent_B <- subset_samples(merged_phy_percent, Group=="B")
merged_phy_percent_C <- subset_samples(merged_phy_percent, Group=="C")
merged_phy_percent_B_C <- merge_phyloseq(merged_phy_percent_B, merged_phy_percent_C)

# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent_B_C) {
  sd <- sample_data(merged_phy_percent_B_C)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent_B_C) {
  OTU <- otu_table(merged_phy_percent_B_C)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

meta_phy_out <- pssd2veg(merged_phy_percent_B_C)
shared_phy_out <- psotu2veg(merged_phy_percent_B_C)

dist_bray <- vegdist(shared_phy_out, method="bray", binary=FALSE)
adonis2(formula = dist_bray ~ Group, data = meta_phy_out, permutations=10000)





merged_phy_percent_B <- subset_samples(merged_phy_percent, Group=="B")
merged_phy_percent_D <- subset_samples(merged_phy_percent, Group=="D")
merged_phy_percent_B_D <- merge_phyloseq(merged_phy_percent_B, merged_phy_percent_D)

# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent_B_D) {
  sd <- sample_data(merged_phy_percent_B_D)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent_B_D) {
  OTU <- otu_table(merged_phy_percent_B_D)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

meta_phy_out <- pssd2veg(merged_phy_percent_B_D)
shared_phy_out <- psotu2veg(merged_phy_percent_B_D)

dist_bray <- vegdist(shared_phy_out, method="bray", binary=FALSE)
adonis2(formula = dist_bray ~ Group, data = meta_phy_out, permutations=10000)





merged_phy_percent_C <- subset_samples(merged_phy_percent, Group=="C")
merged_phy_percent_D <- subset_samples(merged_phy_percent, Group=="D")
merged_phy_percent_C_D <- merge_phyloseq(merged_phy_percent_C, merged_phy_percent_D)

# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent_C_D) {
  sd <- sample_data(merged_phy_percent_C_D)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent_C_D) {
  OTU <- otu_table(merged_phy_percent_C_D)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

meta_phy_out <- pssd2veg(merged_phy_percent_C_D)
shared_phy_out <- psotu2veg(merged_phy_percent_C_D)

dist_bray <- vegdist(shared_phy_out, method="bray", binary=FALSE)
adonis2(formula = dist_bray ~ Group, data = meta_phy_out, permutations=10000)







merged_phy_percent_A_C_prune <- prune_taxa(names(sort(taxa_sums(merged_phy_percent_A_C),TRUE)[1:40]), merged_phy_percent_A_C)
plot_heatmap(merged_phy_percent_A_C_prune, method = "PCoA", distance = "bray", "Test", "Rank1", trans = NULL, sample.order="Test", low="azure1", high="dodgerblue3", na.value="white")


merged_phy_percent_B_D_prune <- prune_taxa(names(sort(taxa_sums(merged_phy_percent_B_D),TRUE)[1:40]), merged_phy_percent_B_D)
plot_heatmap(merged_phy_percent_B_D_prune, method = "PCoA", distance = "bray", "Test", "Rank1", trans = NULL, sample.order="Test", low="azure1", high="dodgerblue3", na.value="white")




















#Extract for glm



###Cecum and colon


setwd("/Users/zebra/Desktop/GF/")

meta_ob1  <- read.table("/Users/zebra/Desktop/GF/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("/Users/zebra/Desktop/GF/taxa.txt"))


meta_ob1 <- meta_ob1[meta_ob1$Type=="Cecum",]
meta_ob1 <- meta_ob1[meta_ob1$Group!="A" & meta_ob1$Group!="C" & meta_ob1$Origin=='Neonatal' ,]
#meta_ob1 <- meta_ob1[meta_ob1$Origin=='Neonatal' ,]
#meta_ob1 <- meta_ob1[meta_ob1$Group!="A" & meta_ob1$Group!="C" ,]



nrow(meta_ob1)

#Get sample names
GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtu))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))


#Rarefy
#minimum <- min(sample_sums(phyloseqobj.f))
out1 <- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1278327, rngseed = 1, replace = FALSE)
#out1<- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1000, rngseed = 1, replace = FALSE)

#Count to percent
#phy1 <- phyloseq_standardize_otu_abundance(out1, method = "total")





###Small_gut


setwd("/Users/zebra/Desktop/GF/")

meta_ob2  <- read.table("/Users/zebra/Desktop/GF/meta.txt", header=T, row.names=1)
tax_ob2 <- import_mothur(mothur_constaxonomy = ("/Users/zebra/Desktop/GF/taxa.txt"))


meta_ob2 <- meta_ob2[meta_ob2$Type=="Colon",]
meta_ob2 <- meta_ob2[meta_ob2$Group!="A" & meta_ob2$Group!="C" & meta_ob2$Origin=='Neonatal' ,]
#meta_ob1 <- meta_ob1[meta_ob1$Group!="A" & meta_ob1$Group!="C" ,]



nrow(meta_ob2)

#Get sample names
GroupSamples<-list(rownames(meta_ob2))

nrow(meta_ob2)

#Subset based on sample names
fullMothur<-read.csv("shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtu))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob2), tax_table(tax_ob2))


#Rarefy
#minimum <- min(sample_sums(phyloseqobj.f))
#out2<- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 50347, rngseed = 1, replace = FALSE)
out2 <- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1278327, rngseed = 1, replace = FALSE)


#Count to percent
#phy2 <- phyloseq_standardize_otu_abundance(out2, method = "total")





#Rarefactiuons curves of count data
merged_phy_counts <- merge_phyloseq(out1, out2)

#merged_phy_counts <- phyloseq_standardize_otu_abundance(merged_phy_counts, method = "total")
#rarecurve(t(otu_table(merged_phy_counts)), step=50, cex=0.5)

#merged_phy_counts <- subset_samples(merged_phy_counts, Type=="Cecum")

#merged_phy_counts <- prune_taxa(names(sort(taxa_sums(merged_phy_counts),TRUE)[1:10]), merged_phy_counts)

#GPr  = transform_sample_counts(merged_phy_counts, function(x) x / sum(x) )
#merged_phy_counts = filter_taxa(GPr, function(x) mean(x) > 0.000001, TRUE); merged_phy_counts
           


# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_counts) {
  sd <- sample_data(merged_phy_counts)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_counts) {
  OTU <- otu_table(merged_phy_counts)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}



meta <- pssd2veg(merged_phy_counts)
meta$Test <- factor(meta$Test, levels = c("PBS_Iso", "IL1_Iso"))
ASV_table <- psotu2veg(merged_phy_counts)


levels(meta$Test)




nrow(meta)
nrow(ASV_table)


package_version("glmmTMB")

#note that metadata and ASV table must have same rownames

library(glmmTMB)
#sets seed for reproducibility
set.seed(1)
#instantiates res which holds the results in the outer loop
res=NULL
#outer loop
#for all groups except the first (i.e. the control group)
for(group in levels(meta$Test)[-1]){
  #subsets metadata on rows by `Group` column to `Control`s and `group`
  meta_subset=meta[meta$Test%in%c("PBS_Iso",group),]
  #`tres` instantiated. Temporary results
  tres=NULL
  #for each ASV in `ASV_table` with ASVs on columns
  for (ASV in colnames(ASV_table)){
    ##NB binary
    #copies over one ASV column to subsetted metadata and names it column `Y`
    meta_subset$Y=ASV_table[rownames(meta_subset),ASV]
    #fits a generalized linear mixed model using template model builder for negative binomial model with offset by the log of the total number of reads in a sample
    #may not work for non-mixed effects
    #modnb=glmmTMB::glmmTMB(Y~Group+GA2+BMI60+Nulliparity+Preterm2+(1|Main_Index), data = meta_subset, offset=log(nR),family=glmmTMB::nbinom2)
    #modnb=glmmTMB::glmmTMB(Y~Origin + Test + (1|Mouse_ID), data = meta_subset, family=glmmTMB::nbinom2)
    modnb=glmmTMB::glmmTMB(Y~Type + Test + (1|dam_ID), data = meta_subset, family=glmmTMB::nbinom2)
    #gets the estimate and P value
    nbb=summary(modnb)$coef$cond[2,c("Estimate","Pr(>|z|)")]
    #stores the temporary results into `tres`
    tres=rbind(tres,nbb)
    #lets you know how far along the loop is by printing the ASV
    print(ASV)
  }
  #correction for false discovery rate
  tres=cbind(tres,p.adjust(tres[,2],"fdr"))
  #renames the column names of `tres`
  colnames(tres)<-c(paste("Coef_",group,"vsControl",sep=""),paste("P_",group,"vsControl",sep=""),paste("Q_",group,"vsControl",sep=""))
  #stores looped through ASV model results into `res`
  res=cbind(res,tres)
}
#adds a Taxonomy column to the results
res=data.frame(Taxonomy=colnames(ASV_table),res)
#reorders results by p value in the `paste("P_",group,"vsControl",sep="")` column
#rename to group of interest
#res=res[order(res$P_PTL_PPROMvsControl),]
write.csv(res,file="ASV_model_results.csv")













#Rarefactiuons curves of count data
merged_phy_counts <- merge_phyloseq(out1, out2)
#rarecurve(t(otu_table(merged_phy_counts)), step=50, cex=0.5)

#merged_phy_counts_Cecum <- subset_samples(merged_phy_counts, Type=="Cecum")

merged_phy_percent <- phyloseq_standardize_otu_abundance(merged_phy_counts, method = "total")




# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent) {
  sd <- sample_data(merged_phy_percent)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent) {
  OTU <- otu_table(merged_phy_percent)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}




meta <- pssd2veg(merged_phy_percent)
ASV_table_percent <- psotu2veg(merged_phy_percent)


#nrow(meta)
#nrow(ASV_table_percent)


df <- as.data.frame(ASV_table_percent)
meta <- as.data.frame(meta)





abund <- as.data.frame(cbind(df$G134, df$G12, df$G73, df$G76, df$G13, df$G45, df$G15, df$G30, df$G9, df$G69, df$G96, df$G39, df$G32, df$G112, df$G71, df$G122, df$G100, df$G78, df$G72, df$G5, meta))
colnames(abund)[1] <- "Lachnoclostridium_sp_YL32"
colnames(abund)[2] <- "Flintibacter_sp_KGMB00164"
colnames(abund)[3] <- "Anaerobutyricum_hallii"
colnames(abund)[4] <- "Enterocloster_bolteae"
colnames(abund)[5] <- "Enterococcus_faecalis"
colnames(abund)[6] <- "Clostridium_scindens"
colnames(abund)[7] <- "Clostridium_innocuum"
colnames(abund)[8] <- "Acutalibacter_muris"
colnames(abund)[9] <- "Intestinimonas_butyriciproducens"
colnames(abund)[10] <- "Anaerostipes_hadrus"
colnames(abund)[11] <- "Dysosmobacter_welbionis"
colnames(abund)[12] <- "Flavonifractor_plautii"
colnames(abund)[13] <- "Ruthenibacterium_lactatiformans"
colnames(abund)[14] <- "Oscillibacter_sp_NSJ_62"
colnames(abund)[15] <- "Lachnospira_eligens"
colnames(abund)[16] <- "Turicibacter_sp_H121"
colnames(abund)[17] <- "Curtobacterium_flaccumfaciens"

abund$df$G12

G12_p <- ggplot(abund,aes(x=Order, y=df$G12, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()

G73_p <- ggplot(abund,aes(x=Order, y=df$G121, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()

G73_p <- ggplot(abund,aes(x=Order, y=df$G147, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()

G73_p <- ggplot(abund,aes(x=Order, y=df$G132, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()

G73_p <- ggplot(abund,aes(x=Order, y=df$G112, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()







G12_p <- ggplot(abund,aes(x=Order, y=df$G12, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()

G73_p <- ggplot(abund,aes(x=Order, y=df$G73, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()

G76_p <- ggplot(abund,aes(x=Order, y=df$G76, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()

G13_p <- ggplot(abund,aes(x=Order, y=df$G13, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()

G45_p <- ggplot(abund,aes(x=Order, y=df$G45, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()

G15_p <- ggplot(abund,aes(x=Order, y=df$G15, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()

G30_p <- ggplot(abund,aes(x=Order, y=df$G30, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()

G9_p <- ggplot(abund,aes(x=Order, y=df$G9, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()

G69_p <- ggplot(abund,aes(x=Order, y=df$G69, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()

G96_p <- ggplot(abund,aes(x=Order, y=df$G96, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()

G39_p <- ggplot(abund,aes(x=Order, y=df$G39, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()

G32_p <- ggplot(abund,aes(x=Order, y=df$G32, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()

G112_p <- ggplot(abund,aes(x=Order, y=df$G112, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()

G71_p <- ggplot(abund,aes(x=Order, y=df$G71, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()

G122_p <- ggplot(abund,aes(x=Order, y=df$G122, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()

#G100_p <- ggplot(abund,aes(x=Order, y=df$G100, color=Test)) +
  #geom_boxplot(outlier.shape = NA) +
  #geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  #scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  #scale_y_continuous(labels=function(x)x*100) +
  #theme_classic()

G78_p <- ggplot(abund,aes(x=Order, y=df$G78, color=Test)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0, dodge.width = 0.75, seed = NULL), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  scale_y_continuous(labels=function(x)x*100) +
  theme_classic()






library(ggpubr)
ggarrange(G12_p, G73_p, G76_p, G13_p, G45_p, G15_p, G30_p, G9_p, G69_p, G96_p, G39_p, G32_p, G112_p, G71_p, G122_p, G78_p, labels = c("Lachnoclostridium_sp_YL32", "Flintibacter_sp_KGMB00164", "Anaerobutyricum_hallii", "Enterocloster_bolteae", "Enterococcus_faecalis", "Clostridium_scindens", "Clostridium_innocuum", "Acutalibacter_muris", "Intestinimonas_butyriciproducens", "Anaerostipes_hadrus", "Dysosmobacter_welbionis", "Flavonifractor_plautii", "Ruthenibacterium_lactatiformans", "Oscillibacter_sp_NSJ_62", "Lachnospira_eligens", "Curtobacterium_flaccumfaciens"),
          common.legend = TRUE, legend = "right",
          ncol = 3, nrow = 6, align="hv")











subset <- subset_samples(merged_phy_counts, Type!="2Cecum" & Origin=="Neonatal")

merged_phy_percent <- phyloseq_standardize_otu_abundance(subset, method = "total")




# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy_percent) {
  sd <- sample_data(merged_phy_percent)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy_percent) {
  OTU <- otu_table(merged_phy_percent)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}




meta_phy_out <- pssd2veg(merged_phy_percent)
shared_phy_out <- psotu2veg(merged_phy_percent)





dist_bray <- vegdist(shared_phy_out, method="bray", binary=FALSE)
adonis2(formula = dist_bray ~ Type + Test, data = meta_phy_out, permutations=10000, strata=meta_phy_out$pup_ID)
#adonis2(formula = dist_bray ~ dam_ID + Type + Test, data = meta_phy_out, permutations=10000)









dist_bray <- capscale(shared_phy_out~1, distance="bray", binary=FALS, eig=TRUE)

prop_explained <- summary(eigenvals(dist_bray))
prop_explained <- prop_explained[2, 1:2]
prop_explained <- format(round((100 * prop_explained), digits=1), nsmall=1, trim=TRUE)



data.scores = scores(dist_bray, display = "sites", choices = c(1:10))
#summary(data.scores)


df <- data.frame(data.scores)

library(glue)
labx <- c(glue("PCo 1 ({prop_explained [1]}%)"))
laby <- c(glue("PCo 2 ({prop_explained [2]}%)"))


#, shape=meta_phy_out$Type

ggplot(df, aes(x = MDS1, y = MDS2, colour = meta_phy_out$Test, fill=meta_phy_out$Test, shape=meta_phy_out$Type)) +
  stat_ellipse(geom="polygon", level=0.75, alpha=0.25, show.legend=FALSE) +
  geom_point(mapping = aes(colour = meta_phy_out$Test), size=5, alpha=1.0) +
  scale_color_manual(values=c("#B81A1A", "#283DA8", "#56D43A", "#FFFFFF")) +
  scale_fill_manual(values=c("pink", "dodgerblue", "#56D43A", "#FFFFFF")) +
  #coord_fixed()+
  theme_classic() +
  labs(x=labx, y=laby) +
  #xlab("PC1 (22.5%)") +
  #ylab("PC2 (13.8%)") +
  ggtitle("Bacterial Community Structure") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=12, face="bold"), axis.text.y = element_text(color="black"))
#legend.background = element_rect(fill="white", color="black")
#theme(legend.position = c(0.9, 0.4)) +
#annotate("text", x=0.3, y=0.7, label = "F = 1.71, R2 = 0.02, p = 0.6021", size=4) +
#annotate("text", x=0.3, y=0.56, label = "Far/BLK: F = 0.918, R2 = 0.051, p = 0.4614 ", size=4) +
#annotate("text", x=0.3, y=0.42, label = "Near/Far (Mouse): F = 1.015, R2 = 0.4616, p = 0.409", size=4) +
#annotate("text", x=0.3, y=0.28, label = "Near/Far (Site): F = 1.106, R2 = 0.084, p = 0.245", size=4) +
#theme(legend.title=element_blank()) +
#theme(legend.position="bottom", legend.box = "horizontal")


















#rarefaction curves
#remotes::install_github("mahendra-mariadassou/phyloseq-extended", ref = "dev")
library(phyloseq.extended)

p1 <- ggrare(out1, step = 5000, color = "Type", label = NULL, se = FALSE, parallel=TRUE)
p1 <- p1 + facet_wrap(~Type, ncol = 2)
p1

p2 <- ggrare(out2, step = 2000, color = "Type", label = NULL, se = FALSE, parallel=TRUE)
p2 <- p2 + facet_wrap(~Type, ncol = 1)
p2

dev.off()


#Merge phyloseq objects with relative abundances          
merged_phy <- merge_phyloseq(phy1, phy2)



library(microbiome)
meta <- function(x) {
  df <- as(sample_data(x), "data.frame")
  rownames(df) <- sample_names(x)
  df
}

merged_meta_ob <- meta(merged_phy)





# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(merged_phy) {
  sd <- sample_data(merged_phy)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(merged_phy) {
  OTU <- otu_table(merged_phy)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}






meta_phy_out <- pssd2veg(merged_phy)
shared_phy_out <- psotu2veg(merged_phy)




dist_bray <- vegdist(shared_phy_out, method="bray", binary=FALSE)
adonis2(formula = dist_bray ~ Type + Test, data = meta_phy_out, permutations=10000)




dist_bray <- capscale(shared_phy_out~1, distance="bray", binary=FALSE)

summary(eigenvals(dist_bray))




data.scores = scores(dist_bray, display = "sites", choices = c(1:10))
#summary(data.scores)


df <- data.frame(data.scores)


#, shape=meta_phy_out$Type

ggplot(df, aes(x = MDS1, y = MDS2, colour = meta_phy_out$Order)) +
  geom_point(mapping = aes(colour = meta_phy_out$Order, shape=meta_phy_out$Test), size=5, alpha=1.0) +
  #scale_color_manual(values=c("#B81A1A", "#283DA8", "#56D43A", "#FFFFFF")) +
  #coord_fixed()+
  theme_classic() +
  xlab("PC1 (22.5%)") +
  ylab("PC2 (13.8%)") +
  ggtitle("Bacterial Community Structure") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=12, face="bold"), axis.text.y = element_text(color="black"))
#theme(legend.position = c(0.9, 0.4)) +
#annotate("text", x=0.3, y=0.7, label = "F = 1.71, R2 = 0.02, p = 0.6021", size=4) +
#annotate("text", x=0.3, y=0.56, label = "Far/BLK: F = 0.918, R2 = 0.051, p = 0.4614 ", size=4) +
#annotate("text", x=0.3, y=0.42, label = "Near/Far (Mouse): F = 1.015, R2 = 0.4616, p = 0.409", size=4) +
#annotate("text", x=0.3, y=0.28, label = "Near/Far (Site): F = 1.106, R2 = 0.084, p = 0.245", size=4) +
#theme(legend.title=element_blank()) +
#theme(legend.position="bottom", legend.box = "horizontal")








###Maternal

Small_gut_phy <- subset_samples(merged_phy, Type == 'Small_gut' & Origin=='Maternal')
Cecum_phy <- subset_samples(merged_phy, Type == 'Cecum' & Origin=='Maternal')
Colon_phy <- subset_samples(merged_phy, Type == 'Colon' & Origin=='Maternal')




#Small gut


# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(Small_gut_phy) {
  sd <- sample_data(Small_gut_phy)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(Small_gut_phy) {
  OTU <- otu_table(Small_gut_phy)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

Small_gut_phy_meta <- pssd2veg(Small_gut_phy)
Small_gut_phy_shared <- psotu2veg(Small_gut_phy)


Small_gut_dist_bray <- vegdist(Small_gut_phy_shared, method="bray", binary=FALSE)
adonis2(formula = Small_gut_dist_bray ~ Mouse_ID + Group, data = Small_gut_phy_meta)

library(pairwiseAdonis)
pairwise.adonis2(Small_gut_dist_bray ~ Group, data = Small_gut_phy_meta, perm=10000)




#Cecum

# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(Cecum_phy) {
  sd <- sample_data(Cecum_phy)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(Cecum_phy) {
  OTU <- otu_table(Cecum_phy)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

Cecum_phy_meta <- pssd2veg(Cecum_phy)
Cecum_phy_shared <- psotu2veg(Cecum_phy)


Cecum_dist_bray <- vegdist(Cecum_phy_shared, method="bray", binary=FALSE)
adonis(formula = Cecum_dist_bray ~ Mouse_ID + Group, data = Cecum_phy_meta)

library(pairwiseAdonis)
pairwise.adonis2(Cecum_dist_bray ~ Group, data = Cecum_phy_meta, perm=10000)





#Colon


# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(Colon_phy) {
  sd <- sample_data(Colon_phy)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(Colon_phy) {
  OTU <- otu_table(Colon_phy)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

Colon_phy_meta <- pssd2veg(Colon_phy)
Colon_phy_shared <- psotu2veg(Colon_phy)

library(vegan)
Colon_dist_bray <- vegdist(Colon_phy_shared, method="bray", binary=FALSE)
adonis2(formula = Colon_dist_bray ~ Mouse_ID + Group, data = Colon_phy_meta)

library(pairwiseAdonis)
pairwise.adonis2(Colon_dist_bray ~ Group, data = Colon_phy_meta, perm=10000)















###Neonatal

Small_gut_phy <- subset_samples(merged_phy, Type == 'Small_gut' & Origin=='Neonatal')
Cecum_phy <- subset_samples(merged_phy, Type == 'Cecum' & Origin=='Neonatal')
Colon_phy <- subset_samples(merged_phy, Type == 'Colon' & Origin=='Neonatal')




#Small gut


# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(Small_gut_phy) {
  sd <- sample_data(Small_gut_phy)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(Small_gut_phy) {
  OTU <- otu_table(Small_gut_phy)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

Small_gut_phy_meta <- pssd2veg(Small_gut_phy)
Small_gut_phy_shared <- psotu2veg(Small_gut_phy)


Small_gut_dist_bray <- vegdist(Small_gut_phy_shared, method="bray", binary=FALSE)
adonis(formula = Small_gut_dist_bray ~ Mouse_ID + Group, data = Small_gut_phy_meta)

library(pairwiseAdonis)
pairwise.adonis2(Small_gut_dist_bray ~ Group, data = Small_gut_phy_meta, perm=10000)




#Cecum

# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(Cecum_phy) {
  sd <- sample_data(Cecum_phy)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(Cecum_phy) {
  OTU <- otu_table(Cecum_phy)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

Cecum_phy_meta <- pssd2veg(Cecum_phy)
Cecum_phy_shared <- psotu2veg(Cecum_phy)


Cecum_dist_bray <- vegdist(Cecum_phy_shared, method="bray", binary=FALSE)
adonis(formula = Cecum_dist_bray ~ Mouse_ID + Group, data = Cecum_phy_meta)

library(pairwiseAdonis)
pairwise.adonis2(Cecum_dist_bray ~ Group, data = Cecum_phy_meta, perm=10000)





#Colon


# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(Colon_phy) {
  sd <- sample_data(Colon_phy)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(Colon_phy) {
  OTU <- otu_table(Colon_phy)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

Colon_phy_meta <- pssd2veg(Colon_phy)
Colon_phy_shared <- psotu2veg(Colon_phy)


Colon_dist_bray <- vegdist(Colon_phy_shared, method="bray", binary=FALSE)
adonis(formula = Colon_dist_bray ~ Mouse_ID + Group, data = Colon_phy_meta)

library(pairwiseAdonis)
pairwise.adonis2(Colon_dist_bray ~ Group, data = Colon_phy_meta, perm=10000)












###Neonatal Colon


setwd("/Users/zebra/Desktop/GF/")

meta_ob1  <- read.table("/Users/zebra/Desktop/GF/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("/Users/zebra/Desktop/GF/taxa.txt"))
C

meta_ob1 <- meta_ob1[meta_ob1$Origin=="Neonatal" & meta_ob1$Type=="Colon" & meta_ob1$Group!="A" & meta_ob1$Group!="C"  ,]



nrow(meta_ob1)

#Get sample names
GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtu))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))


#Rarefy
#minimum <- min(sample_sums(phyloseqobj.f))
out1<- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1278327, rngseed = 1, replace = FALSE)
#out1<- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1000, rngseed = 1, replace = FALSE)

#Count to percent
phy_neo_colon <- phyloseq_standardize_otu_abundance(out1, method = "total")






#20 most abundant Bacteria tax_ob1a across all samples
phy_neo_colon2 <- prune_taxa(names(sort(taxa_sums(phy_neo_colon ),TRUE)[1:40]), phy_neo_colon )






#plot_heatmap(Jameyt, method = "NMDS", distance = "bray", "ID", "Rank5", trans = NULL, sample.order="ID2")
plot_heatmap(phy_neo_colon2, method = "PCoA", distance = "bray", "Test", "Rank1", trans = NULL, sample.order="Test", low="azure1", high="dodgerblue3", na.value="white")


#plot_bar(phy_neo_colon2, fill="Rank1") + facet_wrap(~Test, scales = "free_x", nrow = 4)






















###Plot neonatal cecum and colon




setwd("/Users/zebra/Desktop/GF/")

meta_ob1  <- read.table("/Users/zebra/Desktop/GF/meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constaxonomy = ("/Users/zebra/Desktop/GF/taxa.txt"))


meta_ob1 <- meta_ob1[meta_ob1$Origin=="Neonatal" & meta_ob1$Group!="C" & meta_ob1$Group!="D" ,]

#& meta_ob1$Type=="Colon" 

nrow(meta_ob1)

#Get sample names
GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtu))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, taxa_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_table(tax_ob1))


#Rarefy
#minimum <- min(sample_sums(phyloseqobj.f))
out1<- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1278327, rngseed = 1, replace = FALSE)
#out1<- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1000, rngseed = 1, replace = FALSE)

#Count to percent
phy_neo_colon <- phyloseq_standardize_otu_abundance(out1, method = "total")



# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(phy_neo_colon) {
  sd <- sample_data(phy_neo_colon)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(phy_neo_colon) {
  OTU <- otu_table(phy_neo_colon)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

meta_phy_out <- pssd2veg(phy_neo_colon)
shared_phy_out <- psotu2veg(phy_neo_colon)




df <- as.data.frame(shared_phy_out)
meta <- as.data.frame(meta_phy_out)





abund <- as.data.frame(cbind(df$G2, df$G5, df$G11, meta))
colnames(abund)[1] <- "Akkermansia_muciniphila"
colnames(abund)[2] <- "Bacteroides_thetaiotaomicron"
colnames(abund)[3] <- "Bifidobacterium_pseudolongum"



G2_p <- ggplot(abund,aes(x=Order, y=Akkermansia_muciniphila, color=Test)) +
  geom_boxplot() +
  geom_point(position=position_jitterdodge(), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  theme_classic()


G5_p <- ggplot(abund,aes(x=Order, y=Bacteroides_thetaiotaomicron, color=Test)) +
  geom_boxplot() +
  geom_point(position=position_jitterdodge(), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  theme_classic()


G11_p <- ggplot(abund,aes(x=Order, y=Bifidobacterium_pseudolongum, color=Test)) +
  geom_boxplot() +
  geom_point(position=position_jitterdodge(), size=3, alpha=0.8) +
  scale_color_manual(values=c("red3", "mediumblue", "darkgreen")) +
  theme_classic()





library(ggpubr)
ggarrange(G2_p, G5_p, G11_p, labels = c("Akkermansia_muciniphila", "Bacteroides_thetaiotaomicron", "Bifidobacterium_pseudolongum"), 
          common.legend = TRUE, legend = "right",
          ncol = 1, nrow = 3, align="hv")







df <- as.data.frame(shared_phy_out)
meta <- as.data.frame(meta_phy_out)



dist_bray <- capscale(df~1, distance="bray", binary=FALSE)

summary(eigenvals(dist_bray))




data.scores = scores(dist_bray, display = "sites", choices = c(1:10))
#summary(data.scores)


df <- data.frame(data.scores)


#, shape=meta_phy_out$Type

p2 <-  ggplot(df, aes(x = MDS1, y = MDS2, colour = meta$Test)) +
  geom_point(mapping = aes(colour = meta$Test), size=5, alpha=0.6) +
  scale_color_manual(values=c("azure3", "mediumblue", "darkgreen")) +
  #coord_fixed()+
  theme_classic() +
  #xlab("PC1 (49.5%)") +
  #ylab("PC2 (18.1%)") +
  #ggtitle("Bray") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=12, face="bold"), axis.text.y = element_text(color="black"))
#theme(legend.position = c(0.9, 0.4)) +
#annotate("text", x=0.3, y=0.7, label = "F = 1.71, R2 = 0.02, p = 0.6021", size=4) +
#annotate("text", x=0.3, y=0.56, label = "Far/BLK: F = 0.918, R2 = 0.051, p = 0.4614 ", size=4) +
#annotate("text", x=0.3, y=0.42, label = "Near/Far (Mouse): F = 1.015, R2 = 0.4616, p = 0.409", size=4) +
#annotate("text", x=0.3, y=0.28, label = "Near/Far (Site): F = 1.106, R2 = 0.084, p = 0.245", size=4) +
#theme(legend.title=element_blank()) +
#theme(legend.position="bottom", legend.box = "horizontal")











lm <- lmer(mergedf$ASV1 ~ df$Test + (1|mergedf$Mouse), data = mergedf, REML=TRUE)
Anova(lm, type="III", test="Chisq")



library(lme4)
library(car)


#lm <- lmer(df$G2 ~ meta$Test, data = df, REML=TRUE)

lm <- aov(df$G2 ~ meta$Type + meta$Test, data=df)
Anova(lm, type="III", test="Chisq")

lm <- aov(df$G5 ~ meta$Type + meta$Test, data=df)
Anova(lm, type="III", test="Chisq")

lm <- aov(df$G11 ~ meta$Type + meta$Test, data=df)
Anova(lm, type="III", test="Chisq")









































# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(Cecum_phy) {
  sd <- sample_data(Cecum_phy)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(Cecum_phy) {
  OTU <- otu_table(Cecum_phy)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

Cecum_phy_meta <- pssd2veg(Cecum_phy)
Cecum_phy_shared <- psotu2veg(Cecum_phy)


Cecum_dist_bray <- vegdist(Cecum_phy_shared, method="bray", binary=FALSE)
adonis(formula = Cecum_dist_bray ~ Mouse_ID + Group, data = Cecum_phy_meta)

pairwise.adonis2(Cecum_dist_bray ~ Group, data = Cecum_phy_meta)





# convert the sample_data() within a merged_phyloseq object to a vegan compatible data object
pssd2veg <- function(Colon_phy) {
  sd <- sample_data(Colon_phy)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a merged_phyloseq object to a vegan compatible data object
psotu2veg <- function(Colon_phy) {
  OTU <- otu_table(Colon_phy)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

Colon_phy_meta <- pssd2veg(Colon_phy)
Colon_phy_shared <- psotu2veg(Colon_phy)



Colon_dist_bray <- vegdist(Colon_phy_shared, method="bray", binary=FALSE)
adonis(formula = Colon_dist_bray ~ Mouse_ID + Group, data = Colon_phy_meta)

pairwise.adonis2(Colon_dist_bray ~ Group, data = Colon_phy_meta)












Small_gut_meta <- sample_data(Small_gut_phy) 
Cecum_meta <- sample_data(Cecum_phy)
Colon_meta <- sample_data(Colon_phy) 







Small_gut_dist <- phyloseq::distance(Small_gut_phy , method = "bray")
pairwise.adonis(Small_gut_dist, sample_data(Small_gut_phy)$Group, perm=10000)

Cecum_dist <- phyloseq::distance(Cecum_phy , method = "bray")
pairwise.adonis(Cecum_dist, sample_data(Cecum_phy)$Group, perm=10000)

Colon_dist <- phyloseq::distance(Colon_phy , method = "bray")
pairwise.adonis(Colon_dist, sample_data(Colon_phy)$Group, perm=10000)








#adonis
adonis(distance(Small_gut_phy, method="bray") ~ Mouse_ID+Origin+Type+Group+Injection+Type, data = merged_meta_ob, permutations=999)

#adonis
adonis(distance(merged_phy, method="bray") ~ Type, data = merged_meta_ob, permutations=999, strata=merged_meta_ob$Mouse_ID)

#adonis
adonis(distance(merged_phy, method="bray") ~ Type+Group+Origin+Injection+Type, data = merged_meta_ob, strata=merged_meta_ob$Mouse_ID, permutations=999)



























dist.uf <- phyloseq::distance(merged_phy , method = "bray")
pairwise.adonis(dist.uf, sample_data(merged_phy)$Group, perm=10000)



data_df <- sample_data(merged_phy) 



dist_bray <- capscale(shared_phy_out~1, distance="bray", binary=FALSE)





pairwise.adonis(dist.uf ~ meta_phy_out$Mouse_ID, strata = 'Mouse_ID', data = data_df)


nrow(dist.uf)
nrow(meta_phy_out)












dist_bray <- capscale(shared_phy_out~1, distance="bray", binary=FALSE)

dist_bray <- vegdist(shared_phy_out, method="bray", binary=FALSE)








#write.csv(as.matrix(meta), "meta.csv")
#write.csv(as.matrix(shared), "shared.csv")


dist_bray <- capscale(dist_bray~1, distance="bray", binary=FALSE)

#summary(eigenvals(dist_bray))

data.scores = scores(dist_bray, display = "sites", choices = c(1:10))
#summary(data.scores)


df <- data.frame(data.scores)




p2 <-  ggplot(df, aes(x = MDS1, y = MDS2, colour = meta$Mouse_ID, shape=meta$Origin)) +
  geom_point(mapping = aes(colour = meta$Mouse_ID), size=5) +
  scale_color_manual(values=c("azure4", "#333399", "#A50021", "chartreuse4", "goldenrod1")) +
  #coord_fixed()+
  theme_classic() +
  #xlab("PC1 (49.5%)") +
  #ylab("PC2 (18.1%)") +
  #ggtitle("Bray") +
  theme(axis.text.x = element_text(color="black"), text = element_text(size=12, face="bold"), axis.text.y = element_text(color="black"))
#theme(legend.position = c(0.9, 0.4)) +
#annotate("text", x=0.3, y=0.7, label = "F = 1.71, R2 = 0.02, p = 0.6021", size=4) +
#annotate("text", x=0.3, y=0.56, label = "Far/BLK: F = 0.918, R2 = 0.051, p = 0.4614 ", size=4) +
#annotate("text", x=0.3, y=0.42, label = "Near/Far (Mouse): F = 1.015, R2 = 0.4616, p = 0.409", size=4) +
#annotate("text", x=0.3, y=0.28, label = "Near/Far (Site): F = 1.106, R2 = 0.084, p = 0.245", size=4) +
#theme(legend.title=element_blank()) +
#theme(legend.position="bottom", legend.box = "horizontal")








#20 most abundant Bacteria tax_ob1a across all samples
gpt1 <- prune_taxa(names(sort(taxa_sums(phy1),TRUE)[1:20]), phy1)


#plot_heatmap(gpt1, method = "PCoA", distance = "bray", "ID", "Rank1", trans = NULL, sample.order="ID", low="azure1", high="dodgerblue3", na.value="white")










#20 most abundant Bacteria tax_ob1a across all samples
gpt2 <- prune_taxa(names(sort(taxa_sums(phy2),TRUE)[1:30]), phy2)



#plot_heatmap(gpt2, method = "PCoA", distance = "bray", "ID", "Rank1", trans = NULL, sample.order="ID", low="azure1", high="dodgerblue3", na.value="white")








#20 most abundant Bacteria tax_ob1a across all samples
gpt3 <- prune_taxa(names(sort(taxa_sums(merged_phy),TRUE)[1:30]), merged_phy)


#plot_heatmap(gpt, method = "NMDS", distance = "bray", "ID", "Rank5", trans = NULL, sample.order="ID2")
plot_heatmap(gpt3, method = "PCoA", distance = "bray", "ID", "Rank1", trans = NULL, sample.order="ID", low="azure1", high="dodgerblue3", na.value="white")








#adonis
adonis(distance(phy, method="bray") ~ ID, data = meta_ob1, permutations=999)

#Plot distribution of ASVs across all samples
plot(sort(tax_ob1a_sums(phy), TRUE), type="h", ylim=c(0, 3), xlim=c(0,30))
#plot_bar(phyloseqobj.f, fill="Rank2") + facet_wrap(~ID, scales = "free_x", nrow = 4)

#20 most abundant Bacteria tax_ob1a across all samples
gpt <- prune_tax_ob1a(names(sort(tax_ob1a_sums(phy),TRUE)[1:30]), phy)

#Barplot of selected ASVs with tax_ob1onomy
plot_ob <- plot_bar(gpt, fill="Rank1") + facet_wrap(~ID, scales = "free_x", nrow = 4)

#plot_heatmap(gpt, method = "NMDS", distance = "bray", "ID", "Rank5", trans = NULL, sample.order="ID2")
plot_heatmap(gpt, method = "PCoA", distance = "bray", "ID", "Rank1", trans = NULL, sample.order="ID", low="azure1", high="dodgerblue3", na.value="white")



dev.off()





















#AF_Near to BLKs

setwd("/Users/awinters/Desktop/")

meta_ob1  <- read.table("/Users/awinters/Desktop/All_meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constax_ob1onomy = ("/Users/awinters/Desktop/All_tax_ob1.txt"))

#Sort into controls and age ranges
meta_ob1 <-meta_ob1[meta_ob1$Type1=="AF" & meta_ob1$Group!="AF_F",]

#Get sample names
GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtu))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, tax_ob1a_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_ob1_table(tax_ob1))

#Rarefy
#minimum <- min(sample_sums(phyloseqobj.f))
out<- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 2546, rngseed = 1, replace = FALSE)

#Count to percent
phy <- phyloseq_standardize_otu_abundance(out, method = "total")

#adonis
adonis(distance(phy, method="bray") ~ ID, data = meta_ob1, permutations=999)

#Plot distribution of ASVs across all samples
plot(sort(tax_ob1a_sums(phy), TRUE), type="h", ylim=c(0, 3), xlim=c(0,100))

#20 most abundant Bacteria tax_ob1a across all samples
gpt <- prune_tax_ob1a(names(sort(tax_ob1a_sums(phy),TRUE)[1:20]), phy)

#Barplot of selected ASVs with tax_ob1onomy
plot_bar(gpt, fill="Rank6") + facet_wrap(~Group, scales = "free_x", nrow = 2)

#plot_heatmap(gpt, method = "NMDS", distance = "bray", "ID", "Rank5", trans = NULL, sample.order="ID2")
plot_heatmap(gpt, method = "PCoA", distance = "bray", "ID", "Rank6", trans = NULL, sample.order="Type2", low="azure1", high="dodgerblue3", na.value="white")












#AF_Far to BLKs

setwd("/Users/awinters/Desktop/")

meta_ob1  <- read.table("/Users/awinters/Desktop/All_meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constax_ob1onomy = ("/Users/awinters/Desktop/All_tax_ob1.txt"))

#Sort into controls and age ranges
meta_ob1 <-meta_ob1[meta_ob1$Type1=="AF" & meta_ob1$Group!="AF_N",]

#Get sample names
GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtu))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, tax_ob1a_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_ob1_table(tax_ob1))

#Rarefy
#minimum <- min(sample_sums(phyloseqobj.f))
out<- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 2546, rngseed = 1, replace = FALSE)

#Count to percent
phy <- phyloseq_standardize_otu_abundance(out, method = "total")

#adonis
adonis(distance(phy, method="bray") ~ ID, data = meta_ob1, permutations=999)

#Plot distribution of ASVs across all samples
plot(sort(tax_ob1a_sums(phy), TRUE), type="h", ylim=c(0, 3), xlim=c(0,100))

#20 most abundant Bacteria tax_ob1a across all samples
gpt <- prune_tax_ob1a(names(sort(tax_ob1a_sums(phy),TRUE)[1:20]), phy)

#Barplot of selected ASVs with tax_ob1onomy
plot_bar(gpt, fill="Rank6") + facet_wrap(~Group, scales = "free_x", nrow = 2)

#plot_heatmap(gpt, method = "NMDS", distance = "bray", "ID", "Rank5", trans = NULL, sample.order="ID2")
plot_heatmap(gpt, method = "PCoA", distance = "bray", "ID", "Rank6", trans = NULL, sample.order="Type2", low="azure1", high="dodgerblue3", na.value="white")











#AF_Near to AF_Far

setwd("/Users/awinters/Desktop/")

meta_ob1  <- read.table("/Users/awinters/Desktop/All_meta.txt", header=T, row.names=1)
tax_ob1 <- import_mothur(mothur_constax_ob1onomy = ("/Users/awinters/Desktop/All_tax_ob1.txt"))

#Sort into controls and age ranges
meta_ob1 <-meta_ob1[meta_ob1$Type1=="AF" & meta_ob1$Type2!="Control",]
#meta_ob1 <-meta_ob1[meta_ob1$Type1=="AF" & meta_ob1$Type2!="Control" & meta_ob1$Low=="No",]

#Get sample names
GroupSamples<-list(rownames(meta_ob1))

#Subset based on sample names
fullMothur<-read.csv("shared.csv",)
AF_mothur<-fullMothur[fullMothur$Group %in% c(unlist(GroupSamples)),]

#Reformat files to go into otu tables
AF_mothur<-subset(AF_mothur, select = -c(label, numOtu))
AF_mothur2 <- AF_mothur[,-1]
rownames(AF_mothur2) <- AF_mothur[,1]

#Make otu table and phyloseq object
otu_tab<-otu_table(AF_mothur2, tax_ob1a_are_rows = FALSE)
phyloseqobj.f <- phyloseq(otu_table(otu_tab), sample_data(meta_ob1), tax_ob1_table(tax_ob1))

#Decontaminating phyloseq object
#Bring in contaminants in a single column in a csv file
contam_list<-read.csv("AF_contam_list.csv",)
#Take your otu table out of your phyloseq object
otu_pract=as(otu_table(phyloseqobj.f),"matrix")
#Select only non contaminant otus from your extracted otu table
decontam_otu_tab<-otu_pract[, !(colnames(otu_pract) %in% unlist(contam_list))]
#Make these into a list
decon_list<-colnames(decontam_otu_tab)
#Prune phyloseq object using your list of good otus
phyloseqobj.f=prune_tax_ob1a(unlist(decon_list),phyloseqobj.f)

#Rarefy
#minimum <- min(sample_sums(phyloseqobj.f))
out<- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1366, rngseed = 1, replace = FALSE)
#out<- rarefy_even_depth(physeq = phyloseqobj.f, sample.size = 1947, rngseed = 1, replace = FALSE)

#Count to percent
phy <- phyloseq_standardize_otu_abundance(out, method = "total")

#adonis
adonis(distance(phy, method="bray") ~ ID, data = meta_ob1, permutations=999, strata=meta_ob1$Mouse)

#Plot distribution of ASVs across all samples
plot(sort(tax_ob1a_sums(phy), TRUE), type="h", ylim=c(0, 3), xlim=c(0,100))

#20 most abundant Bacteria tax_ob1a across all samples
gpt <- prune_tax_ob1a(names(sort(tax_ob1a_sums(phy),TRUE)[1:100]), phy)

#Barplot of selected ASVs with tax_ob1onomy
plot_bar(gpt, fill="Rank6") + facet_wrap(~Group, scales = "free_x", nrow = 2)

#plot_heatmap(gpt, method = "NMDS", distance = "bray", "ID", "Rank5", trans = NULL, sample.order="ID2")
plot_heatmap(gpt, method = "PCoA", distance = "bray", "ID", "Rank6", trans = NULL, sample.order="Type2", low="azure1", high="black", na.value="white")









