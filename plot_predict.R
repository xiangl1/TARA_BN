## plot prediction value
setwd("~/Dissertation/tara_BN/tara_bac/")

library(ggplot2)
library(tidyr)
library(plyr)
library(RColorBrewer)
library(dplyr)
library(scales)
library(funrar)
library(data.table)

pre_data <- read.csv('test_bn/bac_boots_no_meta.csv')

# long format
long_df <- gather(data = pre_data, key = Phylum_class , value = Abundance,-1)

# get plot colour
colourCount <- length(pre_data[,-1])
getPalette <- colorRampPalette(brewer.pal(12,"Paired"),bias = 2)

# make stacked barplot object
phylum.heatmap <- ggplot(data=long_df, aes(x=sample_name, y = Abundance, fill=Phylum_class)) + 
  geom_bar(stat = 'identity') +
  geom_tile() +
  xlab(label = "Sampling Time") +
  ylab(label = "Relative Abundance (%)") + 
  #scale_x_discrete(labels = function(x) strftime(strptime(x,'%Y-%m-%d'), '%b-%d'))+
  scale_y_continuous(expand = c(0,0),limits = c(0,100)) +
  #facet_grid(Level~Year, switch = "x", scales = "free_x", space = "free_x")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(strip.text.x = element_text(size = 18)) + 
  theme(strip.text.y = element_text(size = 18)) +
  theme(axis.title.x=element_text(size = 18),
        axis.text.x=element_text(angle = 45, size = 12,hjust = 1))+
  #axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_text(size = 18),
        axis.text.y=element_text(size = 12)
        #axis.ticks.y=element_blank()
  ) +
  theme(panel.spacing = unit(0.5, "lines")) +
  
  # match colour with lake mendota
  scale_fill_manual(name ='Phylum_class',values = getPalette(colourCount))+
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5,size = 20)) +
  ggtitle("Tara Bacterial Relative Abundance" )

# write plot to tiff
tiff("test_figure",units="in", width=30, height=10,res = 300)
print(phylum.heatmap)
dev.off()