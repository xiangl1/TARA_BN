# Tara bacterial abundance BN test (continuous data only)
setwd("~/Dissertation/tara_BN/tara_bac/test_bn/")

library(bnlearn)
library(Rgraphviz)
library(tibble)
library(Metrics)

# input data
main_df <- read.csv('bac_BN_with_meta.csv')
complete_df <- na.omit(main_df)

# drop catagorical environmental variables
bacterial_data <- complete_df[,c(-1,-4,-5,-6,-7,-8,-9,-10)]

# extract column names of different variable type
bac_name <- colnames(bacterial_data[c(16:ncol(bacterial_data))])
geo_meta <- colnames(bacterial_data[c(1,2)])
phys_meta <- colnames(bacterial_data[c(3:7)])
nutri_meta <- colnames(bacterial_data[c(8:15)])

# log10 transform relative abundance
bacterial_data[c(16:ncol(bacterial_data))] <- log10(bacterial_data[c(16:ncol(bacterial_data))]+1)

#----------------------------------------------------------
## Case 1: no expert information
#cv_bn <- bn.cv(bacterial_data,bn="hc",fit = "mle",runs = 10)

#----------------------------------------------------------
## Case 2: with expert information
# 1. Latitude and longitude have no parents, 
#    Temperature parents: latitude, longitude and PAR.day
#    PAR.day parents: Latitude and longitude
# 2. whitelist: latitude -> temperature,longitude -> temperature, PAR.day -> Temperature 
#               latitude -> PAR.day, longitude -> PAR.day,
#               PAR.day -> chlorophyll_a, PAR.day-> Cyanobacteria,
# 3. no arcs between N, P, Fe and Si

#-------------------------------------------------------------------------
# make whitelist
wl <- matrix(c("PAR.day","Temperature",
               "PAR.day","Cyanobacteria_Chloroplast",
               "PAR.day","Cyanobacteria_Cyanobacteria"), 
             ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))

#-------------------------------------------------------------------------
# make blacklist
# block out arcs to Latitude and longitude
bl_latitude_from <- colnames(bacterial_data)[-1]
bl_latitude_to <-rep("Latitude",length(colnames(bacterial_data)[-1]))

bl_longi_from <- colnames(bacterial_data)[-2]
bl_longi_to <-rep("Longitude",length(colnames(bacterial_data)[-2]))

# block out arcs to temperature except latitude,longitude and PAR
bl_temp_from <- colnames(bacterial_data)[-c(1,2,4,6)]
bl_temp_to <-rep("Temperature",length(colnames(bacterial_data)[-c(1,2,4,6)]))

# block out arcs to PAR.day except Latitude and longitude
bl_par_from <- colnames(bacterial_data)[-c(1,2,6)]
bl_par_to <-rep("PAR.day",length(colnames(bacterial_data)[-c(1,2,6)]))

# block out arcs between N,P,Fe and Si
# block out N~P,N~Fe, and N~Si 
bl_n_from <- c(rep("Nitrite.NO2..",3),rep("Nitrate.NO3..",3),rep("Nitrate.and.Nitrite",3),
               rep("X.NH4..",3),rep(c("Phosphate","Fe.tot","Silicate"),4))
bl_n_to <- c(rep(c("Phosphate","Fe.tot","Silicate"),4),rep("Nitrite.NO2..",3),
             rep("Nitrate.NO3..",3),rep("Nitrate.and.Nitrite",3),rep("X.NH4..",3))

# block out P~Fe, P~Si and Fe~Si
bl_other_from <- c("Phosphate","Phosphate","Fe.tot","Fe.tot","Silicate","Silicate")
bl_other_to <- c("Fe.tot","Silicate","Silicate","Phosphate","Fe.tot","Phosphate")

# make final blacklist
bl <- cbind(c(bl_latitude_from,bl_longi_from,bl_temp_from,bl_par_from,bl_n_from,bl_other_from),
            c(bl_latitude_to,bl_longi_to,bl_temp_to,bl_par_to,bl_n_to,bl_other_to))
colnames(bl) <- c("from","to")
#-------------------------------------------------------------------------
# do cross-validation with maxparents = inf
cv_bn <- bn.cv(bacterial_data,bn="hc",
               algorithm.args = list(whitelist=wl,blacklist=bl,maxp=7),
               fit = "mle",runs = 10)

#-------------------------------------------------------------------------
# make average structure
arc_str <- custom.strength(cv_bn,cpdag=FALSE)
ave_dag_cv <- averaged.network(arc_str,threshold = 0.5)

# evaluate bary-curtis similarity and MSE
BC_cv_ave <- matrix(nrow=10,ncol=9)
mse_cv_ave <- matrix(nrow=10,ncol=23)

for (k in 1:10) {
  pre_df_list <- list()
  pre_log_list <- list()
  obser_log_list <- list()
  obser_df_list <- list()
  BC_cv <- matrix(nrow=10,ncol=9)
  mse_cv <- matrix(nrow=10,ncol=23)
  
  # make predicted file of test data
  for (i in 1:10) {
    test_bac <- bacterial_data[cv_bn[[k]][[i]]$test,]
    pre_log_list[[i]] <- as.data.frame(matrix(ncol = length(bac_name),nrow = nrow(test_bac)))
    pre_df_list[[i]] <- as.data.frame(matrix(ncol = length(bac_name),nrow = nrow(test_bac)))
    obser_df_list[[i]] <- 10^test_bac[,c(16:ncol(test_bac))]-1
    obser_log_list[[i]]<- test_bac[,c(16:ncol(test_bac))]  
    
    # make predict dataframes
    for (j in 1:length(bac_name)) {
      pre_log_list[[i]][,j] <- predict(cv_bn[[k]][[i]]$fitted,node=bac_name[j],data=test_bac,method = "bayes-lw")
    }  
    # save for mse calculation
    
    # transfer back to relative abundance for BC calculation, recalculate percentage and extract taxa data from test data
    pre_df_list[[i]] <-10^pre_log_list[[i]] -1
    pre_df_list[[i]][pre_df_list[[i]] < 0 ] <- 0
    pre_df_list[[i]] <- as.data.frame(t(apply(pre_df_list[[i]],1,function(x) {round(x*100/sum(x),3)})))
  }

  # calculate Bray-Curtis dissimilarity
  for (a in 1:10) {
    for (b in 1:nrow(pre_df_list[[a]])) {
      pred <- pre_df_list[[a]][b,]
      obser <- obser_df_list[[a]][b,]
      BC_cv[a,b] <- 2*sum(pmin(pred,obser))/(sum(pred)+sum(obser))
    }
  }
  # calculate predict errors (mse for continues)
  for (a in 1:10) {
    for (b in 1:ncol(pre_log_list[[a]])) {
      mse_cv[a,b] <- mse(pre_log_list[[a]][,b],obser_log_list[[a]][,b])
    }
  }
  BC_cv_ave[k,] <- colMeans(BC_cv,na.rm=TRUE)
  mse_cv_ave[k,] <- colMeans(mse_cv,na.rm=TRUE)
}

# get final average BC and MSE
bc_cv_final <- round(rowMeans(BC_cv_ave),3)
mse_cv_final <- round(colMeans(mse_cv_ave),3)

# write files
capture.output(ave_dag_cv, file = "../pre_results/p7_expert_cv/p7_expert_bn_cv_summary.txt")
capture.output(ave_dag_cv$nodes, file = "../pre_results/p7_expert_cv/p7_expert_bn_cv_nodes.txt")
write.csv(ave_dag_cv$arcs,"../pre_results/p7_expert_cv/p7_expert_bn_cv_arcs.csv")
write.csv(as.data.frame(bc_cv_final),"../pre_results/p7_expert_cv/p7_expert_cv_bc.csv",row.names = FALSE)
mse_cv_final <- t(as.data.frame(mse_cv_final))
colnames(mse_cv_final) <- bac_name
write.csv(mse_cv_final,"../pre_results/p7_expert_cv/p7_expert_cv_mse.csv",row.names = FALSE)

#-------------------------------------------------------------------------
# plot BN

jpeg("../pre_results/p7_expert_cv/bn_cv_p7_expert.jpg",units="in",width =30,height = 20,res = 300)
bn_boots <- strength.plot(ave_dag_cv,arc_str,shape = "ellipse",layout = "fdp",render=FALSE)

# modify network
# modify edges
bac_a <- arcs(subgraph(ave_dag_cv, bac_name))
bac_a <- as.character(interaction(bac_a[,"from"], bac_a[,"to"], sep = "~"))

geo_bac_a <- arcs(subgraph(ave_dag_cv, c(geo_meta,bac_name)))
geo_bac_a <- as.character(interaction(geo_bac_a[,"from"], geo_bac_a[,"to"], sep = "~"))

phys_bac_a <- arcs(subgraph(ave_dag_cv, c(phys_meta,bac_name)))
phys_bac_a <- as.character(interaction(phys_bac_a[,"from"], phys_bac_a[,"to"], sep = "~"))

nutri_bac_a <- arcs(subgraph(ave_dag_cv, c(nutri_meta,bac_name)))
nutri_bac_a <- as.character(interaction(nutri_bac_a[,"from"], nutri_bac_a[,"to"], sep = "~"))

geo_a <- arcs(subgraph(ave_dag_cv, geo_meta))
geo_a <- as.character(interaction(geo_a[,"from"], geo_a[,"to"], sep = "~"))

phys_a <- arcs(subgraph(ave_dag_cv,phys_meta))
phys_a <- as.character(interaction(phys_a[,"from"], phys_a[,"to"], sep = "~"))

nutri_a <- arcs(subgraph(ave_dag_cv, nutri_meta))
nutri_a <- as.character(interaction(nutri_a[,"from"], nutri_a[,"to"], sep = "~"))

edgeRenderInfo(bn_boots)$col = "grey"
edgeRenderInfo(bn_boots)$col[geo_bac_a] = "red"
edgeRenderInfo(bn_boots)$col[phys_bac_a] = "red"
edgeRenderInfo(bn_boots)$col[nutri_bac_a] = "red"
edgeRenderInfo(bn_boots)$col[bac_a] = "darkgreen"

edgeRenderInfo(bn_boots)$col[geo_a] = "grey"
edgeRenderInfo(bn_boots)$col[phys_a] = "grey"
edgeRenderInfo(bn_boots)$col[nutri_a] = "grey"

# modify node
nodeRenderInfo(bn_boots)$fill[geo_meta] <- "yellow3"
nodeRenderInfo(bn_boots)$fill[phys_meta] <- "lightblue"
nodeRenderInfo(bn_boots)$fill[nutri_meta] <- "slateblue1"
nodeRenderInfo(bn_boots)$fill[bac_name] <- "limegreen"
nodeRenderInfo(bn_boots)$fontsize <- 18
nodeRenderInfo(bn_boots)$rWidth[bac_name] <- 225

renderGraph(bn_boots)
dev.off()
rm(list=ls())