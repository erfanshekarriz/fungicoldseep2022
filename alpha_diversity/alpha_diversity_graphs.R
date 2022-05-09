library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(patchwork)
library(rstudioapi)
library(vegan)
library(phyloseq)
library(ggpubr)
library(ggpmisc)



setwd(dirname(getActiveDocumentContext()$path))       # Set working directory to source file location
getwd()                                               # Check updated working directory


### ALPHA DIVERSITY ###

#16S
physeqA <- readRDS("/Users/erfanshekarriz/Desktop/Hongbin LIU Lab/haima_24_soil_samples/bacteriaA_16s/phyloseq_objects/physeq.BacteriaA.ALLasv.rds")
sample_sums(physeqA)
max(sample_sums(physeqA))
min(sample_sums(physeqA))
physeqA <- prune_taxa(taxa_sums(physeqA) > 0, physeqA)

# rarecurve(t(otu_table(physeqA)), step=200, cex=0.5)
physeqA.rar <- rarefy_even_depth(physeqA, 
                                 sample.size = min(sample_sums(physeqA)),
                                 rngseed = 193, 
                                 replace = TRUE, 
                                 trimOTUs = TRUE, 
                                 verbose = TRUE)

plot_richness(physeqA.rar, 
              x="end_depth", 
              color = "layer",
              shape = "site") 

richness.plot.A <- plot_richness(physeqA.rar, 
                                 x="end_depth", 
                                 color = "layer",
                                 shape = "site",
                                 measures = c("Chao1", "Shannon","Simpson")) 
richness.plot.A

plot.data.A <- richness.plot.A$data


#ITS
physeqB <- readRDS("/Users/erfanshekarriz/Desktop/Hongbin LIU Lab/haima_24_soil_samples/fungi_ITS/phyloseq_objects/physeq.ITS.ALLasv.rds")
sample_sums(physeqB)
max(sample_sums(physeqB))
min(sample_sums(physeqB))
physeqB <- prune_taxa(taxa_sums(physeqB) > 0, physeqB)

# rarecurve(t(otu_table(physeqB)), step=200, cex=0.5)
physeqB.rar <- rarefy_even_depth(physeqB, 
                                 sample.size = min(sample_sums(physeqB)),
                                 rngseed = 193, 
                                 replace = TRUE, 
                                 trimOTUs = TRUE, 
                                 verbose = TRUE)


plot_richness(physeqB.rar, 
              x="end_depth", 
              color = "layer",
              shape = "site") 

richness.plot.B <- plot_richness(physeqB.rar, 
                                 x="end_depth", 
                                 color = "layer",
                                 shape = "site",
                                 measures = c("Chao1", "Shannon","Simpson")) 
richness.plot.B

plot.data.B<- richness.plot.B$data



#### plot 
plot.data.B %>% filter(variable =="Chao1") %>%
  ggplot(aes(x= end_depth, y=value)) + 
  geom_point(size=2, aes(shape = site, color=site)) + 
  geom_smooth(method="glm",se=TRUE,
              formula=y ~ x + I(x^2), 
              size =0.4, 
              linetype="dashed") + 
  stat_poly_eq(formula =  y ~ x + I(x^2), 
               aes(label = paste(..rr.label.., ..p.value.label.., sep = "*`,`~")), 
               parse = TRUE,
               label.x.npc = "left",
               vstep = 0.05, 
               size = 3) + # sets vertical spacing
  # stat_regline_equation(aes(label = ..rr.label.., color=Site),
  #                        formula = y ~ x + I(x^2), position = "identity") +
  # stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "*`,`~")),
  #          label.y= Inf, label.x = Inf, vjust = 1, hjust = 1.1, size = 3)+
  # scale_color_grey()+ 
  theme_bw() + 
  labs(title = "ITS Shannon Index Alpha Diversity ") +
  xlab("Sampling Depth") + 
  ylab("Chao1 Index Value")



ggsave("ITS.Chao1.Alpha.Diversity.png" , 
       width = 12,
       height = 10,
       units = "cm",
       dpi = 1000 )


plot.data.A %>% filter(variable =="Chao1") %>%
  filter(samples != "ROV3.8.0") %>%
  ggplot(aes(x= end_depth, y=value)) + 
  geom_point(size=2, aes(shape = site, color=site)) + 
  geom_smooth(method="lm",se=TRUE,
              formula=y ~ x, 
              size =0.4, 
              linetype="dashed") + 
  stat_poly_eq(formula =  y ~ x, 
               aes(label = paste(..rr.label.., ..p.value.label.., sep = "*`,`~")), 
               parse = TRUE,
               label.x.npc = "left",
               vstep = 0.05, 
               size = 3) + # sets vertical spacing
  theme_bw() + 
  labs(title = "16S Shannon Index Alpha Diversity ") +
  xlab("Sampling Depth") + 
  ylab("Chao1 Index Value")



ggsave("16S.Shannon.Alpha.Diversity.png" , 
       width = 12,
       height = 10,
       units = "cm",
       dpi = 300 )

