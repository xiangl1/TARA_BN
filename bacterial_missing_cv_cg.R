# Tara bacterial abundance BN test (incomplete data)
setwd("~/Dissertation/tara_BN/tara_bac/test_bn/")

library(bnlearn)
library(Rgraphviz)
library(tibble)
library(Metrics)

# input data
main_df <- read.csv('bac_BN_with_meta.csv')

# drop catagorical environmental variables
bacterial_data <- main_df[,c(-1,-5,-6,-7,-10)]

# remove SO extrme environment
#bacterial_data <- bacterial_data[!bacterial_data$Env.feature == "MIX",]
#bacterial_data <- bacterial_data[!bacterial_data$OS.region == "SO",]
#bacterial_data <- bacterial_data[!bacterial_data$Season == "summer",]
#bacterial_data<- droplevels(bacterial_data)

# extract column names of different variable type
bac_name <- colnames(bacterial_data[c(19:ncol(bacterial_data))])
geo_meta <- colnames(bacterial_data[c(1:5)])
phys_meta <- colnames(bacterial_data[c(6:10)])
nutri_meta <- colnames(bacterial_data[c(11:18)])

# log10 transform relative abundance
bacterial_data[bac_name] <- log10(bacterial_data[bac_name]+1)

#----------------------------------------------------------
## Case 1: with expert information and incomplete data
# 1. Categorical variables have no continuious parents:
#    Latitude have no parents
#    Longitude have no partents
#    OS.region parents: Latitude,Longitude
#    Env.feature have no parents
#    Season have no parents
#    Temperature have parents:latitude, longitude,PAR.day, Env.feature
#    PAR.day have parents:Latitude, longitude,Env.feature and season
#    Oxygen have parent: Env.feature

# 2. whitelist: latitude->temperature,Longitude -> temperature, PAR.day -> Temperature,
#               Env.feature->temperature,
#               latitude->PAR.day, longitude->PAR.day,Env.feature->PAR.day, Season->PAR.day
#               latitude-> Os.region, Longitude -> Os.region,
#               Env.feature -> Oxygen,
#               PAR.day-> chlorophyll_a and PAR.day->Cyanobacteria
# 3. no arcs between N, P, Fe and Si, no other arcs to categorical variables

#-------------------------------------------------------------------------
# make whitelist
wl <- matrix(c("Env.feature","Temperature","PAR.day","Temperature",
               "Env.feature","PAR.day","Env.feature","Oxygen",
               "PAR.day","Cyanobacteria_Chloroplast",
               "PAR.day","Cyanobacteria_Cyanobacteria"), 
             ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))

#-------------------------------------------------------------------------
# make blacklist
# block out all arcs to Latitude and longitude
bl_latitude_from <- colnames(bacterial_data)[-1]
bl_latitude_to <-rep("Latitude",length(colnames(bacterial_data)[-1]))

bl_longi_from <- colnames(bacterial_data)[-2]
bl_longi_to <-rep("Longitude",length(colnames(bacterial_data)[-2]))

# block out all arcs to Env.feature
bl_env_from <- colnames(bacterial_data)[-3]
bl_env_to <-rep("Env.feature",length(colnames(bacterial_data)[-3]))

# block out all arcs to OS.region except latitude and longitude
bl_os_from <- colnames(bacterial_data)[-c(1,2,4)]
bl_os_to <-rep("OS.region",length(colnames(bacterial_data)[-c(1,2,4)]))

# block out all arcs to Season
bl_season_from <- colnames(bacterial_data)[-5]
bl_season_to <-rep("Season",length(colnames(bacterial_data)[-5]))

# block out arcs to temperature except latitude, longitude, PAR, Season,and Env.feature
bl_temp_from <- colnames(bacterial_data)[-c(1,2,3,9)]
bl_temp_to <-rep("Temperature",length(colnames(bacterial_data)[-c(1,2,3,9)]))

# block out arcs to PAR.day except Latitude, Longitude, Env.feature and Season
bl_par_from <- colnames(bacterial_data)[-c(1,2,3,5)]
bl_par_to <-rep("PAR.day",length(colnames(bacterial_data)[-c(1,2,3,5)]))

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
bl <- cbind(c(bl_latitude_from,bl_longi_from,bl_env_from,bl_os_from,bl_season_from,
              bl_temp_from,bl_par_from,bl_n_from,bl_other_from),
            c(bl_latitude_to,bl_longi_to,bl_env_from,bl_os_from,bl_season_from,
              bl_temp_to,bl_par_to,bl_n_to,bl_other_to))
colnames(bl) <- c("from","to")
#-------------------------------------------------------------------------
# do 10-fold 10 times cross-validation manully because bn.cv not supprot structure EM

# creat matrix to store bc and mse
BC_cv_ave <- matrix(nrow=10,ncol=14)
mse_cv_ave <- matrix(nrow=10,ncol=23)

run_list <- rep(list(list()), 10) 

#-------------------------------------------------------------------------
# do cross-validation
for (run in 1:10) {
  rm(list=setdiff(ls(), c("bacterial_data","BC_cv_ave","run_list","mse_cv_ave","wl","bl","run",
                          "bac_name","geo_meta","phys_meta","nutri_meta")))
  pre_df_list <- list()
  pre_log_list <- list()
  obser_log_list <- list()
  obser_df_list <- list()
  BC_cv <- matrix(nrow=10,ncol=14)
  mse_cv <- matrix(nrow=10,ncol=23)
  
  # shuffle the data to get unbiased splits
  n <- nrow(bacterial_data)
  k <- 10
  kcv <- split(sample(n), seq_len(k))
  
  run_list[[run]] <- rep(list(list()), 10)
  
  # make predicted file of test data
  for (i in 1:10) {
    # get test and training data
    test_ind <- kcv[[i]]
    train_bac <- bacterial_data[-test_ind,]
    test_bac <- bacterial_data[test_ind,]
    
    # do structural EM learning from missing data
    bn_dag <- structural.em(train_bac,maximize = "hc",
                            maximize.args = list(whitelist = wl, blacklist=bl,maxp=5),
                            fit = "mle", max.iter = 10,return.all = TRUE)
    fit <- bn_dag$fitted
    
    # save bn class
    run_list[[run]][[i]] <- bn_dag$dag
    
    # impute missing value in test data
    if (anyNA(test_bac)) {
      impute_test <- impute(fit,test_bac)
    } else {
      impute_test <- test_bac
    }
    
    # store predict files 
    pre_log_list[[i]] <- as.data.frame(matrix(ncol = length(bac_name),nrow = nrow(test_bac)))
    pre_df_list[[i]] <- as.data.frame(matrix(ncol = length(bac_name),nrow = nrow(test_bac)))
    obser_df_list[[i]] <- 10^test_bac[,c(19:ncol(test_bac))]-1
    obser_log_list[[i]]<- test_bac[,c(19:ncol(test_bac))] 

    # make predict dataframes
    for (j in 1:length(bac_name)) {
      pre_log_list[[i]][,j] <- predict(fit,node=bac_name[j],data=impute_test,method = "bayes-lw")
    }  
    
    # transfer back to relative abundance, recalculate percentage and extract taxa data from test data
    pre_df_list[[i]] <-10^pre_log_list[[i]] -1
    pre_df_list[[i]][pre_df_list[[i]] < 0 ] <- 0
    pre_df_list[[i]] <- as.data.frame(t(apply(pre_df_list[[i]],1,function(x) round(x*100/sum(x),3))))
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
  BC_cv_ave[run,] <- colMeans(BC_cv,na.rm=TRUE)
  mse_cv_ave[run] <- colMeans(mse_cv,na.rm=TRUE)
}

# get final average BC and MSE
bc_cv_final <- round(rowMeans(BC_cv_ave),3)
mse_cv_final <- round(colMeans(mse_cv_ave),3)

#-------------------------------------------------------------------------
# make average structure
arclist <- list()

for (i in 1:10) {
  for (j in 1:10)
    arclist[[length(arclist) + 1]] = arcs(run_list[[i]][[j]])
}

nodes = unique(unlist(arclist))

arc_str <- custom.strength(arclist, nodes = nodes,cpdag = FALSE)
ave_dag_cv <- averaged.network(arc_str,threshold = 0.5)

#-------------------------------------------------------------------------
# write files
capture.output(ave_dag_cv, file = "../pre_results/p1_missing_cg/p1_missing_bn_cg_summary.txt")
capture.output(ave_dag_cv$nodes, file = "../pre_results/p1_missing_cg/p1_missing_bn_cg_nodes.txt")
write.csv(ave_dag_cv$arcs,"../pre_results/p1_missing_cg/p1_missing_bn_cg_arcs.csv")
write.csv(as.data.frame(bc_cv_final),"../pre_results/p1_missing_cg/p1_missing_cg_bc.csv",row.names = FALSE)
mse_cv_final <- t(as.data.frame(mse_cv_final))
colnames(mse_cv_final) <- bac_name
write.csv(mse_cv_final,"../pre_results/p1_missing_cg/p1_missing_cg_mse.csv",row.names = FALSE)

#-------------------------------------------------------------------------
# plot BN

jpeg("../pre_results/p1_missing_cg/p1_bn_missing_cg.jpg",units="in",width =30,height = 20,res = 300)
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