# Tara bacterial abundance BN compare different modeling
setwd("~/Dissertation/tara_BN/tara_bac/test_bn/")

library(ggplot2)
library(tidyr)
library(plyr)
library(RColorBrewer)
library(dplyr)
library(scales)
library(funrar)
library(data.table)

bc_df <- read.csv('../pre_results/cv_bc.csv')
# long format
long_df <- gather(data = bc_df, key = BC)

# add group
#group <- rep(c("No expert knowledge","Continous only","With Discrete","Discrete Extra"),c(9,72,72,72))
group <- rep(c("No expert knowledge","Continous only","With Discrete"),c(9,72,72))
#group <- rep(c("No expert knowledge","Continous only"),c(9,72))
long_df <- cbind(long_df,group)
  
bc_boxplot <- ggplot(data=long_df, aes(x = BC, y = value,fill = group)) + 
  geom_boxplot() +
  #facet_wrap(~group,scales = "free_x")+
  xlab(label = "Model name") +
  ylab(label = "Bray-Curtis similarity") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(strip.text.x = element_text(size = 18)) + 
  theme(axis.title.x=element_text(size = 18),
        axis.text.x=element_text(angle = 45, size = 12,hjust = 1))+
  #axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_text(size = 18),
        axis.text.y=element_text(size = 12)
        #axis.ticks.y=element_blank()
  ) +
  theme(panel.spacing = unit(0.5, "lines")) +
  ggtitle("BC Similarity of Different Modeling Strategies" ) +
  theme(plot.title = element_text(lineheight=.8, size= 18,face="bold",hjust=0.5))

# write plot to jpg
jpeg("cv_cg_bc_figure.jpg",units="in", width=15, height=10,res = 300)
print(bc_boxplot)
dev.off()

#---------------------------------------------------------------------------------------
# plot mse_df
mse_df <- read.csv('../pre_results/cg_mse.csv')
long_mse <- gather(data = mse_df, key = Phylum_class,value = mse,-1)

mse_lineplot <- ggplot(data=long_mse, aes(x = Phylum_class, y = mse,color = Models,group = Models)) + 
  geom_line() +
  geom_point() +
  #scale_color_manual(values = brewer.pal(12,"Paired"))+
  xlab(label = "Organism") +
  ylab(label = "MSE") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(strip.text.x = element_text(size = 18)) + 
  theme(axis.title.x=element_text(size = 18),
        axis.text.x=element_text(angle = 45,size = 12,hjust = 1))+
  #axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_text(size = 18),
        axis.text.y=element_text(size = 12)
        #axis.ticks.y=element_blank()
  ) + 
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"in")) +
  ggtitle("MSE of Different Modeling Strategies (With discrete)" ) +
  theme(plot.title = element_text(lineheight=.8, size= 18,face="bold",hjust=0.5))
  
# write plot to jpg
jpeg("../pre_results/mse_cg.jpg", units="in", width=20, height=10,res = 300)
print(mse_lineplot)
dev.off()
