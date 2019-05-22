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
## Case 1: no expert information
## 1.2 cross-validation approach: 10-fold, 10 times
cv_bn <- bn.cv(bacterial_data,bn="hc",method = "hold-out", k=10, m=9,fit = "mle",runs = 10)

#-------------------------------------------------------------------------
# make average structure
arc_str <- custom.strength(cv_bn)
ave_dag_cv <- averaged.network(arc_str,threshold = 0.7)

#-------------------------------------------------------------------------
# evaluate bary-curtis similarity and MSE
BC_cv <- matrix(nrow=10,ncol=9)
mse_cv <- matrix(nrow=10,ncol=23)

for (k in 1:10) {
  pre_df_list <- list()
  obser_df_list <- list()
  # make predicted file of test data
  for (i in 1:10) {
    test_bac <- bacterial_data[cv_bn[[k]][[i]]$test,]
    pre_df_list[[i]] <- as.data.frame(matrix(ncol = length(bac_name),nrow = nrow(test_bac)))
    obser_df_list[[i]] <- 10^test_bac[,c(14:ncol(test_bac))]-1
  
    # make predict dataframes
    for (j in 1:length(bac_name)) {
      pre_df_list[[i]][,j] <- predict(cv_bn[[k]][[i]]$fitted,node=bac_name[j],data=test_bac)
    }  
    # transfer back to relative abundance, recalculate percentage and extract taxa data from test data
    pre_df_list[[i]] <-10^pre_df_list[[i]] -1
    pre_df_list[[i]][pre_df_list[[i]] < 0 ] <- 0
    pre_df_list[[i]] <- as.data.frame(t(apply(pre_df_list[[i]],1,function(x) {round(x*100/sum(x),3)})))
  }

  # calculate Bray-Curtis dissimilarity
  for (i in 1:10) {
    for (j in 1:nrow(pre_df_list[[i]])) {
      pred <- pre_df_list[[i]][j,]
      obser <- obser_df_list[[i]][j,]
      BC_cv[k,j] <- 2*sum(pmin(pred,obser))/(sum(pred)+sum(obser))
    }
  }
  # calculate predict errors (mse for continues)
  for (i in 1:10) {
    for (j in 1:ncol(pre_df_list[[i]])) {
      mse_cv[k,j] <- round(mse(pre_df_list[[i]][,j],obser_df_list[[i]][,j]),3)
    }
  }
}

# get final average BC and MSE
bc_cv_ave <- colMeans(BC_cv)
mse_cv_ave <- colMeans(mse_cv)

#-------------------------------------------------------------------------
# plot BN

tiff("bn_cv",units="in",width =20,height = 15,res = 300)
bn_boots <- strength.plot(ave_dag_cv,arc_str,shape = "ellipse",layout = "fdp",render=FALSE)

# modify network
# modify edges
bac_a <- arcs(subgraph(ave_dag_cv, bac_name))
bac_a <- as.character(interaction(bac_a[,"from"], bac_a[,"to"], sep = "~"))

geo_a <- arcs(subgraph(ave_dag_cv, geo_meta))
geo_a <- as.character(interaction(geo_a[,"from"], geo_a[,"to"], sep = "~"))

phys_a <- arcs(subgraph(ave_dag_cv, phys_meta))
phys_a <- as.character(interaction(phys_a[,"from"], phys_a[,"to"], sep = "~"))

nutri_a <- arcs(subgraph(ave_dag_cv, nutri_meta))
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