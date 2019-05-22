# Tara bacterial abundance BN test (with metadata)
setwd("~/Dissertation/tara_BN/tara_bac/test_bn/")

library(bnlearn)
library(Rgraphviz)
library(tibble)
library(Metrics)

# input data
main_df <- read.csv('bac_BN_with_meta.csv')
complete_df <- na.omit(main_df)

# drop catagorical environmental variables
bacterial_data <- complete_df[,c(-1,-3,-4,-5,-6,-7,-8,-9,-10)]

# extract column names of different variable type
bac_name <- colnames(bacterial_data[c(14:ncol(bacterial_data))])
geo_meta <- colnames(bacterial_data[1])
phys_meta <- colnames(bacterial_data[c(2:5)])
nutri_meta <- colnames(bacterial_data[c(6:13)])

# log10 transform relative abundance
bacterial_data[c(14:ncol(bacterial_data))] <- log10(bacterial_data[c(14:ncol(bacterial_data))]+1)

#----------------------------------------------------------
# split data to train and test set at 10-fold
smp_size <- floor(0.9 * nrow(bacterial_data))
set.seed(123)
train_ind <- sample(seq_len(nrow(bacterial_data)),size = smp_size)

train_bac <- bacterial_data[train_ind,]
test_bac <- bacterial_data[-train_ind,]

#----------------------------------------------------------
## Case 1: no expert information
# 1.1 bootstrap approach, hill-climbing with BIC score on training data
boot_st <- boot.strength(train_bac, R = 200, algorithm = "hc")
ave_dag <- averaged.network(boot_st,threshold = 0.7)
#set.arc(ave_dag,"Other","Cyanobacteria_Chloroplast")
#set.arc(ave_dag,"Planctomycetes_Planctomycetacia","Verrucomicrobia_Spartobacteria")

# parameter learning MLE with training data
fit <- bn.fit(ave_dag,train_bac)
#----------------------------------------------------------

# make predicted file of test data
pre_df <- data.frame(matrix(ncol = length(bac_name),nrow = nrow(test_bac)))

for (i in 1:length(bac_name)) {
  pre_df[,i] <- predict(fit,node=bac_name[i],data=test_bac)
}

# transfer back to relative abundance, recalculate percentage and extract taxa data from test data
pre_df <- 10^pre_df -1
pre_df[pre_df < 0] <- 0
pre_df <- as.data.frame(t(apply(pre_df,1,function(x) {round(x*100/sum(x),3)})))
obser_df <- 10^test_bac[,c(14:ncol(test_bac))] -1 

#----------------------------------------------------------
# evaluate predict errors
# calculate Bray-Curtis dissimilarity
BC_boot <- c()
for (i in 1:nrow(pre_df)) {
  pred <- pre_df[i,]
  obser <- obser_df[i,]
  BC_boot[i] <- 2*sum(pmin(pred,obser))/(sum(pred)+sum(obser))
}

# calculate predict errors (mse for continues)
mse_boot <- c()
for (i in 1:ncol(pre_df)) {
  mse_boot[i] <- round(mse(pre_df[,i],obser_df[,i]),3)
}
# write files
#write.csv(output,'test_out.csv',row.names=FALSE)
#write.csv(bc_df,"bc_test.csv",row.names=FALSE)

#----------------------------------------------------------
# plot BN

tiff("bn_boots",units="in",width =15,height = 10,res = 300)
bn_boots <- strength.plot(ave_dag,boot_st,shape = "ellipse",layout = "fdp",render=FALSE)

# modify network
# modify edges
bac_a <- arcs(subgraph(ave_dag, bac_name))
bac_a <- as.character(interaction(bac_a[,"from"], bac_a[,"to"], sep = "~"))

geo_a <- arcs(subgraph(ave_dag, geo_meta))
geo_a <- as.character(interaction(geo_a[,"from"], geo_a[,"to"], sep = "~"))

phys_a <- arcs(subgraph(ave_dag, phys_meta))
phys_a <- as.character(interaction(phys_a[,"from"], phys_a[,"to"], sep = "~"))

nutri_a <- arcs(subgraph(ave_dag, nutri_meta))
nutri_a <- as.character(interaction(nutri_a[,"from"], nutri_a[,"to"], sep = "~"))

edgeRenderInfo(bn_boots)$col = "grey"
edgeRenderInfo(bn_boots)$col[geo_a] = "yellow3"
edgeRenderInfo(bn_boots)$col[phys_a] = "lightblue"
edgeRenderInfo(bn_boots)$col[nutri_a] = "slateblue1"
edgeRenderInfo(bn_boots)$col[bac_a] = "limegreen"

# modify node
nodeRenderInfo(bn_boots)$fill[geo_meta] <- "yellow3"
nodeRenderInfo(bn_boots)$fill[phys_meta] <- "lightblue"
nodeRenderInfo(bn_boots)$fill[nutri_meta] <- "slateblue1"
nodeRenderInfo(bn_boots)$fill[bac_name] <- "limegreen"
nodeRenderInfo(bn_boots)$fontsize <- 18
nodeRenderInfo(bn_boots)$rWidth[bac_name] <- 225

renderGraph(bn_boots)
dev.off()