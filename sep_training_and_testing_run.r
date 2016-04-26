####################################################################
#### Load Functions
####################################################################
#packages
require(rgdal)
require(raster)
require(maptools)
require(randomForest)
#setwd("~/Dropbox/Caileigh_stuff/YFDP")
setwd("C:/Users/shootc/Dropbox/Caileigh_stuff/YFDP/")

#Setup Functions
setup.predictors <- function(file_list = cgs_under5m_outlier_fileList, 
                             mask_raster = example_raster_2m, 
                             end_of_colname = "under_5m"){
  #Clip the LiDAR and Topo predictors to the Species Matrix
  Predictors.matrix <- Reduce("cbind", lapply(lapply(lapply(file_list, raster), crop, y = raster(mask_raster)), values))
  fileNameSplit <- sapply(strsplit(file_list, split="[_/]"), function(x) x[-1])
  colnames(Predictors.matrix) <-  sapply(fileNameSplit, 
                                         function(x) paste(paste(x[-c(1:5, length(x))], collapse="_"), end_of_colname, collapse="_"))
  print("Predictor Matrix has been created")
  return(Predictors.matrix)
}

CreateMatchRaster <- function(data_in, cellsize, proj4)
{
  width  <- extent(data_in)@xmax - extent(data_in)@xmin
  height <- extent(data_in)@ymax - extent(data_in)@ymin
  if((width %% cellsize) < (cellsize / 2))
    remainder_w <- -1 * (width %% cellsize)
  else
    remainder_w <- cellsize - (width %% cellsize)
  if((height %% cellsize) < (cellsize / 2))
    remainder_h <- -1 * (height %% cellsize)
  else
    remainder_h <- cellsize - (height %% cellsize)
  ras_match <- raster(nrows=((height + remainder_h) / cellsize),
                      ncols=((width + remainder_w) / cellsize),
                      xmn=extent(data_in)@xmin, xmx=(extent(data_in)@xmax + remainder_w),
                      ymn=extent(data_in)@ymin, ymx=(extent(data_in)@ymax + remainder_h), crs=proj4)
  return(ras_match)
}


species.response <- function(example_raster=example_raster_2m, 
                             species.shp=all.shrub.shp.file,
                             plot_boundary_shapefile=plot.boundary.shp.file){ 
  #set the resamp_resolution to the same res as your other rasters if you don't want to resample
  #set it to a smaller resolution if you do want to resample
  species.response.matrix <- NULL
  
  #rasterize the species shapefile
  plot_boundary <- rasterize(plot_boundary_shapefile, example_raster, 1)
  in.plot <- !is.na(values(plot_boundary))
  
  for(spp in unique(species.shp$Species)){
    current.spp = species.shp[species.shp$Species == spp,]
    print(paste("starting the rasterization and resampling for ", spp))
    spp_presence = rasterize(current.spp,example_raster)
    print(paste("ending rasterization for ", spp))
    
    v = values(spp_presence)
    v[!in.plot] <- NA
    v[is.na(v) & in.plot] <- 0
    v[in.plot & v != 0] <- 1
    values(spp_presence) = v
    species.response.matrix <- cbind(species.response.matrix, v)
  }
  colnames(species.response.matrix) = unique(species.shp$Species)
  return(species.response.matrix)
}

#This function takes the sum of all the rows within the species matrix to create a presence/absence matrix
#Rowsums > 1 are set to 1 outside of the function later on
pres.abs.response <- function(species.matrix){
  pres.abs <- rowSums(species.matrix)
  return(pres.abs)
} 

### This creates the growth form Matrix, it's no longer needed
# growth.form.response <- function(species.matrix, 
#                                  growthForm1.spp=c("SYMO", "ROBR", "CECO","CEPA", "RIRO","RUPA"), 
#                                  growthForm2.spp=c("VAUL","LEDA","RINE","CEIN","CHSE","ARPA"), 
#                                  growthForm3.spp=c("COCOC","RHOC","SARA","COSE","ABCO"), 
#                                  pres.abs = T){
#   spp <- spp.matrix
#   form1 <- rowSums(spp[, growthForm1.spp])
#   form2 <- rowSums(spp[, growthForm2.spp])
#   form3 <- rowSums(spp[, growthForm3.spp])
#   growth.form.matrix <- cbind(form1, form2, form3)
#   if(pres.abs)
#   {
#     return(growth.form.matrix)
#   }
#   
#   growth.forms.numbered <- t(t(growth.form.matrix) * 1:3)
#   all.growth.forms <- as.matrix(rowSums(growth.forms.numbered))
#   return(all.growth.forms) #This multiplies the 3 rows by 1, 2, and 3 in order
# }

make_classError_matrix <- function(NClasses, actual, predicted, ncols = 2,
                                   matrix_rownames = c("Abs", "Pres"), 
                                   matrix_colnames = c("Abs", "Pres")){
  # N classes * actual + predicted 
  # 2 classes for pres/abs, 4 for growth forms, etc.
  values <- NClasses * actual + predicted
  confusion <- matrix(sapply(0:3, function(n) length(which(values == n))), ncol= ncols)
  rownames(confusion) <- matrix_rownames
  colnames(confusion) <- matrix_colnames
  return(confusion)
}

class_error_presAbs <-  function(confusion){
  class_error_absence <- confusion[1, 2] / sum(confusion[1, ])
  class_error_presence <- confusion[2, 1] / sum(confusion[2, ])
  overall_error <- (confusion[1, 2] + confusion[2, 1]) / sum(confusion)
  
  class_error_all <- matrix(data=rbind(class_error_absence, class_error_presence, overall_error), ncol=3)
  colnames(class_error_all) <- c("Class Error Absence", "Class Error Presence", "Overall Error")
  return(class_error_all)
}

verifyParsimonious <- function(predictors, 
                               responses, 
                               rf, 
                               subset, #This ensures that we verify parsimonious predictors on the same parts of the data set that we trained on.
                               npredictors #This is the number of predictors that it should output (top 10 most predictive)
){
  require(randomForest)
  response. = responses[subset] #subset the response set
  imp <- rf$importance #Get a list of most important variable from the rf model
  #Take each of the importance values,divide them by the sum of all the values, then multiply by 100
  #in order to get the %overall importance, then round to the 1s place, and put them in decreasing order
  imp <- round((imp[,1]/sum(imp[,1]))*100,1) 
  imp <- imp[order(imp, decreasing=T)]
  
  # Print the RF run that we're doing, and print the importance list. It may be good to eliminate the importance list soon.
  print(rf)
  #print(imp)
  
  #Get a list of predictor names
  everything.predictors.names <- names(imp) 
  
  
  matrix.names = matrix(data=0, nrow=length(everything.predictors.names), ncol=1)
  rownames(matrix.names) = c(everything.predictors.names)
  names.list = matrix(data=c(everything.predictors.names), nrow=length(everything.predictors.names), ncol=1)
  
  rsq. = 0
  #add in initial rsq of everything
  parsimonious.predictors.names = c(everything.predictors.names[1]) 
  for (jj in 2:npredictors) {
    for (ii in 1:length(matrix.names)) {
      
      tmp.names = c(parsimonious.predictors.names, names.list[ii])
      
      predictors. = predictors[subset, c(tmp.names)]
      
      rf = randomForest(x=predictors., y=as.factor(response.), importance=T)
      # pct.corr = sum(diag(nrow(rf$confusion)) * rf$confusion[, -ncol(rf$confusion)]) / sum(rf$confusion[, -ncol(rf$confusion)])
      pct.corr = 1-rf$err.rate[500]
      matrix.names[ii] = pct.corr#rf$rsq[500]; head(matrix.names)
      #print(paste(jj, ii, round(matrix.names[ii],digits=3), names.list[ii]))
    } # for (ii in 2:length(matrix.names))
    
    length(matrix.names)
    max(matrix.names)
    index. = which(matrix.names==max(matrix.names))
    parsimonious.predictors.names = c(parsimonious.predictors.names, names.list[index.])
    # print(paste(index., names.list[index.], round(matrix.names[index.],digits=3)))
    rsq. = c(rsq., round(matrix.names[index.],digits=3))
    # print(matrix.names)
    # print(names.list)
    matrix.names = matrix.names[c(-index.)]
    names.list = names.list[c(-index.)]
    # print(rsq.)
    #     print("")
    #     print("")
  } # for (jj in 1:10)
  
  #   print(parsimonious.predictors.names)
  #   print(rsq.)
  parsimonious.matrix <- as.matrix(cbind(parsimonious.predictors.names, rsq.))
  print(parsimonious.matrix)
  return(parsimonious.matrix)
} # function(predictors, responses, rf, subset)

####################################################################
#### Load in Data/Prep Data
####################################################################
# variables we need to set up
cgs_no_outlier_fileList         <- c((list.files(path = "YFDP_Metrics_Ground_Points_Retained/Metrics_0to150m/Metrics_2METERS", pattern= ".img", full.names=TRUE)), 
                                     (list.files(path = "YFDP_Metrics_Ground_Points_Retained/Metrics_0to150m/StrataMetrics_2METERS",pattern= ".img", full.names=TRUE)),
                                     (list.files(path = "YFDP_Metrics_Ground_Points_Retained/Metrics_0to150m/StrataCoverMetrics_2METERS",pattern= ".img", full.names=TRUE)),
                                     (list.files(path = "YFDP_Metrics_Ground_Points_Retained/Metrics_0to150m/FineTopoMetrics_2METERS",pattern= ".img", full.names=TRUE)),
                                     (list.files(path = "YFDP_Metrics_Ground_Points_Retained/Metrics_0to150m/Normalized_FineTPIMetrics_2METERS",pattern= ".img", full.names=TRUE))) #Here you create a list of all the files 

cgs_under5m_outlier_fileList     <- c((list.files(path = "YFDP_Metrics_Ground_Points_Retained/Metrics_0to5m/Metrics_2METERS", pattern= ".img", full.names=TRUE)),
                                      (list.files(path = "YFDP_Metrics_Ground_Points_Retained/Metrics_0to5m/StrataCoverMetrics_2METERS",pattern= ".img", full.names=TRUE)),
                                      (list.files(path = "YFDP_Metrics_Ground_Points_Retained/Metrics_0to5m/FineTopoMetrics_2METERS",pattern= ".img", full.names=TRUE)),
                                      (list.files(path = "YFDP_Metrics_Ground_Points_Retained/Metrics_0to5m/Normalized_FineTPIMetrics_2METERS",pattern= ".img", full.names=TRUE)))

all.shrub.shp.file     <- readOGR(dsn ="All_YFDP_NecessaryFiles_UTM10N", layer ="YFDP_Shrub_UTM10N")
plot.boundary.shp.file <- readOGR(dsn ="All_YFDP_NecessaryFiles_UTM10N", layer ="YFDP_Boundary_UTM10N")

example_raster_2m <- crop(x=raster("YFDP_Metrics_Ground_Points_Retained/Metrics_0to5m/Metrics_2METERS/1st_cnt_above_mean_2METERS.img"), 
                          y=plot.boundary.shp.file)
# YFDP.dtm               <- crop(x = raster("LiDAR_data/YFDP/be_YFDP_DEM_100m_buffer.asc"), y = plot.boundary.shp.file)

#Run Functions
all_predictors <- cbind(setup.predictors(cgs_under5m_outlier_fileList, 
                                         mask_raster = example_raster_2m,
                                         end_of_colname = ""),
                        setup.predictors(cgs_no_outlier_fileList,
                                         mask_raster = example_raster_2m,
                                         end_of_colname = ""))

# Set all NA values in LiDAR to 0
all_predictors[is.na(all_predictors)] <- 0
#write.csv(all_predictors, file="all_predictors_NAs_included.csv")

spp.matrix <- species.response(example_raster=example_raster_2m, 
                               species.shp=all.shrub.shp.file,
                               plot_boundary_shapefile=plot.boundary.shp.file)

#This is what causes the speckles in the Rasters - these are the spots with more than one shrub. 
all.complete <- complete.cases(cbind(all_predictors, spp.matrix))
all_predictors_complete  <- all_predictors[all.complete, ]
spp.matrix_complete          <- spp.matrix[all.complete, ]

presAbs.response.cgs <- as.matrix(pres.abs.response(species.matrix=spp.matrix_complete))
presAbs.response.cgs[presAbs.response.cgs > 1] <- 1

topo_complete <- all_predictors_complete[ ,c(which(grepl("FineTopoMetrics", colnames(all_predictors_complete))))]
LiDAR_complete <- all_predictors_complete[ ,-c(which(grepl("FineTopoMetrics", colnames(all_predictors_complete))))]

#Break up testing and training data
plot_boundary <- rasterize( plot.boundary.shp.file, example_raster_2m, 1)
test_raster_AA <- plot_boundary
values(test_raster_AA) <- c(rep(c(1,1,0,1,0,1), each = 67, length.out = 402*41),
                            rep(c(0,1,1,1,0,1), each = 67, length.out = 402*42),
                            rep(c(0,1,0,1,1,1), each = 67, length.out = 402*41),
                            rep(c(0,1,0,1,0,1), each = 67, length.out = 402*42))
plot(test_raster_AA)

v <- as.matrix(values(test_raster_AA))
v <- v[all.complete, ]

#       training_subset <- which(as.matrix(values(test_raster_AA == 1))) #2/3 of data
#       testing_subset <- which(as.matrix(values(test_raster_AA == 0))) #1/3 of data
#       under5m_metrics_topoandLiDAR_complete <- all_predictors_complete[ ,c(which(grepl("0to5m", colnames(all_predictors_complete))))] 
#       upto150m_metrics_topoandLiDAR_complete <- all_predictors_complete[ ,c(which(grepl("0to150m", colnames(all_predictors_complete))))]

######################################################################
#Run Analysis
######################################################################
## Create the run sub
# pres    <- which(presAbs.response.cgs == 1)
# set.seed(100)
# rand.not.pres <- sample(which(presAbs.response.cgs == 0), length(pres))
# set.seed(100)
# pres    <- sample(pres, as.integer(0.67 * length(pres)))

## Run Sub 2
v_pres    <- which(presAbs.response.cgs == 1 & v == 1)
v_rand.not.pres <- sample(which(presAbs.response.cgs == 0 & v == 1), length(v_pres))
v_run.sub       <- c(v_pres, v_rand.not.pres)

######################################
#1
######################################
###All Shrub presence/absence with topo and LiDAR
#Set up results table
v_all.presAbs_topoANDLiDAR.results = matrix(data=NA, ncol=3, nrow=1)
rownames(v_all.presAbs_topoANDLiDAR.results) = colnames(presAbs.response.cgs)
colnames(v_all.presAbs_topoANDLiDAR.results) = c("OOB", "absent", "present")

#Run Random Forest
v_presAbs_topoANDLiDAR_rf <- randomForest(y=as.factor(presAbs.response.cgs[v_run.sub, ]), 
                                        x=(all_predictors_complete[v_run.sub, ]))

v_all.presAbs_topoANDLiDAR.results[1,"OOB"] = round(1-v_presAbs_topoANDLiDAR_rf$err.rate[500,1], digits=3)
v_all.presAbs_topoANDLiDAR.results[1,"absent"] = round(1-v_presAbs_topoANDLiDAR_rf$err.rate[500,2], digits=3)
v_all.presAbs_topoANDLiDAR.results[1,"present"] = round(1-v_presAbs_topoANDLiDAR_rf$err.rate[500,3], digits=3);
v_all.presAbs_topoANDLiDAR.results
#save(presAbs_topoANDLiDAR_rf, file="presAbs_topoANDLiDAR_rf.rda")

#write.csv(all.presAbs_topoANDLiDAR.results, file="all.presAbs_topoANDLiDAR.results.csv")


######################################
#2
######################################
###All Shrub presence/absence with LiDAR
#Set up results table
v_all.presAbs_LiDAR.results = matrix(data=NA, ncol=3, nrow=1)
rownames(v_all.presAbs_LiDAR.results) = colnames(presAbs.response.cgs)
colnames(v_all.presAbs_LiDAR.results) = c("OOB", "absent", "present")

#Run Random Forest
v_presAbs_LiDAR_rf <- randomForest(y=as.factor(presAbs.response.cgs[v_run.sub, ]), x=(LiDAR_complete[v_run.sub,]))

v_all.presAbs_LiDAR.results[1,"OOB"] = round(1-v_presAbs_LiDAR_rf$err.rate[500,1], digits=3)
v_all.presAbs_LiDAR.results[1,"absent"] = round(1-v_presAbs_LiDAR_rf$err.rate[500,2], digits=3)
v_all.presAbs_LiDAR.results[1,"present"] = round(1-v_presAbs_LiDAR_rf$err.rate[500,3], digits=3);
v_all.presAbs_LiDAR.results
#save(presAbs_LiDAR_rf, file="presAbs_LiDAR_rf.rda")

#write.csv(all.presAbs_LiDAR.results, file="all.presAbs_LiDAR.results.csv")

######################################
#3
######################################
###All Shrub presence/absence with topo
#Set up results table
v_all.presAbs_topo.results = matrix(data=NA, ncol=3, nrow=1)
rownames(v_all.presAbs_topo.results) = colnames(presAbs.response.cgs)
colnames(v_all.presAbs_topo.results) = c("OOB", "absent", "present")

#Run Random Forest
v_presAbs_topo_rf <- randomForest(y=as.factor(presAbs.response.cgs[v_run.sub, ]), x=(topo_complete[v_run.sub,]))

v_all.presAbs_topo.results[1,"OOB"] = round(1-v_presAbs_topo_rf$err.rate[500,1], digits=3)
v_all.presAbs_topo.results[1,"absent"] = round(1-v_presAbs_topo_rf$err.rate[500,2], digits=3)
v_all.presAbs_topo.results[1,"present"] = round(1-v_presAbs_topo_rf$err.rate[500,3], digits=3);
v_all.presAbs_topo.results
#save(presAbs_topo_rf, file="presAbs_topo_rf.rda")


########################################## 
#### Run presence/absence predictions ####
##########################################
#### All Shrub Presence/Absence

v_topo_presAbs_predicted <- as.integer(as.character(predict(object=v_presAbs_topo_rf, topo_complete)))
v_topo_presAbs_err_matrix <- make_classError_matrix(NClasses=2, actual=presAbs.response.cgs[v == 0], predicted=v_topo_presAbs_predicted[v == 0])
v_topo_class_error_presAbs <- class_error_presAbs(v_topo_presAbs_err_matrix)

v_LiDAR_presAbs_predicted <- as.integer(as.character(predict(object=v_presAbs_LiDAR_rf, LiDAR_complete)))
v_lidar_presAbs_err_matrix <- make_classError_matrix(NClasses=2, actual = presAbs.response.cgs[v == 0], predicted = v_LiDAR_presAbs_predicted[v == 0])
v_lidar_class_error_presAbs <- class_error_presAbs(v_lidar_presAbs_err_matrix)

v_topo_and_LiDAR_PresAbs_predicted <- as.integer(as.character(predict(object=v_presAbs_topoANDLiDAR_rf, all_predictors_complete)))
v_lidarANDtopo_presAbs_err_matrix <- make_classError_matrix(NClasses=2, actual = presAbs.response.cgs[v == 0], predicted = v_topo_and_LiDAR_PresAbs_predicted[v == 0])
v_lidarANDtopo_class_error_presAbs <- class_error_presAbs(v_lidarANDtopo_presAbs_err_matrix)

#### Locations of shrubs
actual.shrub.locations <- matrix(data = NA, nrow = length(all.complete))

actual.shrub.locations[all.complete] <- presAbs.response.cgs

actual.shrub.locations_raster <- raster(actual.shrub.locations, template = example_raster_2m)
plot(actual.shrub.locations_raster, main = "Actual Shrub Locations", col =c("#d8b365", "#5ab4ac"))

v_only_present_shrub_locations <- actual.shrub.locations_raster
v_only_present_shrub_locations[v_only_present_shrub_locations == 0] <- NA
v_only_present_shrub_locations[values(test_raster_AA) == 1] <- NA
plot(v_only_present_shrub_locations, col = "RED")

##################################
#### Topo predicted locations ####
##################################
v_topo.predicted.allshrub.locations <- matrix(data = NA, nrow = length(all.complete))

v_topo.predicted.allshrub.locations[all.complete] <- v_topo_presAbs_predicted
v_topo.predicted.allshrub.locations[values(test_raster_AA) == 1] <- NA

v_topo.predicted.allshrub.locations_raster <- raster(v_topo.predicted.allshrub.locations, template = example_raster_2m)
#writeRaster(topo.predicted.allshrub.locations_raster + (2*actual.shrub.locations_raster), "topo.predicted.allshrub.locations_raster.tif")


#jpeg("topo.predicted.allshrub.locations_raster.jpg", width = 10, height = 8, units = "in", res = 200)
plot(v_topo.predicted.allshrub.locations_raster, main = "Topo All Shrub Predicted Locations", col =c("#D3D3D3", "#4daf4a"), legend = FALSE)
plot(v_only_present_shrub_locations, col="#e41a1c", add= TRUE, legend = FALSE)
legend("topright", legend = c("Predicted Shrub Absence", "Predicted Shrub Presence", "Actual Shrub Locations"), fill = c("#D3D3D3", "#4daf4a", "#e41a1c"))
#dev.off()

#Topo
jpeg("topo_AllShrubPredicted.jpg", width = 12, height = 10, units = "in", res = 400)
plot(v_topo.predicted.allshrub.locations_raster + (2*actual.shrub.locations_raster), 
     main = "Topography Predicted Locations", legend = F, 
     col = c("#377eb8", "#e41a1c","black","#4daf4a"))
legend("bottomleft", legend = c("Correct Absence",
                                "Error of Omission",
                                "Error of Commission", 
                                "Correct Presence"), 
       fill = c("#377eb8","black", "#e41a1c","#4daf4a"),
       bty="n")
dev.off()

###################################
#### LiDAR predicted locations ####
###################################
v_LiDAR.predicted.allshrub.locations <- matrix(data = NA, nrow = length(all.complete))

v_LiDAR.predicted.allshrub.locations[all.complete] <- v_LiDAR_presAbs_predicted
v_LiDAR.predicted.allshrub.locations[values(test_raster_AA) == 1] <- NA

v_LiDAR.predicted.allshrub.locations_raster <- raster(v_LiDAR.predicted.allshrub.locations, template = example_raster_2m)
#writeRaster(LiDAR.predicted.allshrub.locations_raster + (2*actual.shrub.locations_raster), "LiDAR.predicted.allshrub.locations_raster.tif")

#jpeg("LiDAR.predicted.allshrub.locations_raster.jpg", width = 10, height = 8, units = "in", res = 200)
plot(v_LiDAR.predicted.allshrub.locations_raster, main = "LiDAR All Shrub Predicted Locations",  col =c("#D3D3D3", "#4daf4a"), legend = FALSE)
plot(v_only_present_shrub_locations, col="#e41a1c", add= TRUE, legend = FALSE)
legend("topright", legend = c("Predicted Shrub Absence", "Predicted Shrub Presence", "Actual Shrub Locations"), fill = c("#D3D3D3", "#4daf4a", "#e41a1c"))
#dev.off()

#LiDAR
jpeg("LiDAR_AllShrubPredicted.jpg", width = 12, height = 10, units = "in", res = 400)
plot(v_LiDAR.predicted.allshrub.locations_raster + (2*actual.shrub.locations_raster), 
     main = "LiDAR Predicted Locations", legend = F, 
     col = c("#377eb8", "#e41a1c","black","#4daf4a"))
legend("bottomleft", legend = c("Correct Absence",
                                "Error of Omission",
                                "Error of Commission", 
                                "Correct Presence"), 
       fill = c("#377eb8","black", "#e41a1c","#4daf4a"),
       bty="n")
dev.off()

############################################
#### Topo and LiDAR predicted locations ####
############################################
v_topo_and_LiDAR.predicted.allshrub.locations <- matrix(data = NA, nrow = length(all.complete))

v_topo_and_LiDAR.predicted.allshrub.locations[all.complete] <- v_topo_and_LiDAR_PresAbs_predicted
v_topo_and_LiDAR.predicted.allshrub.locations[values(test_raster_AA) == 1] <- NA

v_topo_and_LiDAR.predicted.allshrub.locations_raster <- raster(v_topo_and_LiDAR.predicted.allshrub.locations, template = example_raster_2m)
#writeRaster(topo_and_LiDAR.predicted.allshrub.locations_raster + (2*actual.shrub.locations_raster), "topo_and_LiDAR.predicted.allshrub.locations_raster.tif")

#jpeg("topo_and_LiDAR.predicted.allshrub.locations_raster.jpg", width = 10, height = 8, units = "in", res = 200)
plot(v_topo_and_LiDAR.predicted.allshrub.locations_raster, main = "Topo and LiDAR All Shrub Predicted Locations", col =c("#D3D3D3", "#4daf4a"), legend = FALSE)
plot(v_only_present_shrub_locations, col="#e41a1c", add= TRUE, legend = FALSE)
legend("topright", legend = c("Predicted Shrub Absence", "Predicted Shrub Presence", "Actual Shrub Locations"), fill = c("#D3D3D3", "#4daf4a", "#e41a1c"))
#dev.off()

#Topo and LiDAR
jpeg("topo_and_LiDAR_AllShrubPredicted.jpg", width = 12, height = 10, units = "in", res = 400)
plot(v_topo_and_LiDAR.predicted.allshrub.locations_raster + (2*actual.shrub.locations_raster), 
     main = "Topography and LiDAR Predicted Locations", legend = F, 
     col = c("#377eb8", "#e41a1c","black","#4daf4a"))
legend("bottomleft", legend = c("Correct Absence",
                                "Error of Omission",
                                "Error of Commission", 
                                "Correct Presence"), 
       fill = c("#377eb8","black", "#e41a1c","#4daf4a"),
       bty="n")
dev.off()


# ########################################################
# #### Add in CWD into predicted locations
# #######################################################
# in.plot <- !is.na(values(plot_boundary))
# 
# cwd <- readOGR(dsn = "All_YFDP_NecessaryFiles_UTM10N", layer = "YFDP_CWD_UTM10N")
# cwd_raster <- rasterize(cwd, example_raster_2m, field="Length", fun="min")
# cwd_raster[cwd_raster >= 0] <- 1
# cwd_raster[is.na(cwd_raster)] <- 0
# cwd_raster[!in.plot] <- NA
# plot(is.na(cwd_raster, "NA Values in CWD Raster"))
# plot(cwd_raster, main = "CWD Raster")
# cwd_presAbs <- as.matrix(values(cwd_raster))
# cwd_presAbs <- cwd_presAbs[all.complete]
# 
# presAbs.response.withCWD <- presAbs.response.cgs + cwd_presAbs
# presAbs.response.withCWD[presAbs.response.withCWD > 1] <- 1
# 
# topo_presAbs_err_matrix_cwdANDshrubs <- make_classError_matrix(NClasses=2, actual=presAbs.response.withCWD, predicted=topo_presAbs_predicted)
# topo_class_error_presAbs_cwdANDshrubs <- class_error_presAbs(topo_presAbs_err_matrix)
# 
# lidar_presAbs_err_matrix_cwdANDshrubs <- make_classError_matrix(NClasses=2, actual = presAbs.response.withCWD, predicted = LiDAR_presAbs_predicted)
# lidar_class_error_presAbs_cwdANDshrubs <- class_error_presAbs(lidar_presAbs_err_matrix)
# 
# lidarANDtopo_presAbs_err_matrix_cwdANDshrubs <- make_classError_matrix(NClasses=2, actual = presAbs.response.withCWD, predicted = topo_and_LiDAR_PresAbs_predicted)
# lidarANDtopo_class_error_presAbs_cwdANDshrubs <- class_error_presAbs(lidarANDtopo_presAbs_err_matrix)
# 
# actual.shrubANDcwd.locations <- matrix(data = NA, nrow = length(all.complete))
# actual.shrubANDcwd.locations[all.complete] <- presAbs.response.withCWD
# actual.shrubANDcwd.locations_raster <- raster(actual.shrubANDcwd.locations, template = example_raster_2m)
# 
# 
# #Topo
# jpeg("topo_AllShrubPredicted_withCWD.jpg", width = 12, height = 10, units = "in", res = 400)
# plot(topo.predicted.allshrub.locations_raster + (2*actual.shrubANDcwd.locations_raster), 
#      main = "Topography Predicted Locations of shrubs and CWD", legend = F, 
#      col = c("#377eb8", "#e41a1c","black","#4daf4a"))
# legend("bottomleft", legend = c("Correct Absence",
#                                 "Error of Omission",
#                                 "Error of Commission", 
#                                 "Correct Presence"), 
#        fill = c("#377eb8","black", "#e41a1c","#4daf4a"),
#        bty="n")
# dev.off()
# 
# #LiDAR
# jpeg("LiDAR_AllShrubPredicted_withCWD.jpg", width = 12, height = 10, units = "in", res = 400)
# plot(LiDAR.predicted.allshrub.locations_raster + (2*actual.shrubANDcwd.locations_raster), 
#      main = "LiDAR Predicted Locations of shrubs and CWD", legend = F, 
#      col = c("#377eb8", "#e41a1c","black","#4daf4a"))
# legend("bottomleft", legend = c("Correct Absence",
#                                 "Error of Omission",
#                                 "Error of Commission", 
#                                 "Correct Presence"), 
#        fill = c("#377eb8","black", "#e41a1c","#4daf4a"),
#        bty="n")
# dev.off()
# #Topo and LiDAR
# jpeg("topo_and_LiDAR_AllShrubPredicted_withCWD.jpg", width = 12, height = 10, units = "in", res = 400)
# plot(topo_and_LiDAR.predicted.allshrub.locations_raster + (2*actual.shrubANDcwd.locations_raster), 
#      main = "Topography and LiDAR Predicted Locations of shrubs and CWD", legend = F, 
#      col = c("#377eb8", "#e41a1c","black","#4daf4a"))
# legend("bottomleft", legend = c("Correct Absence",
#                                 "Error of Omission",
#                                 "Error of Commission", 
#                                 "Correct Presence"), 
#        fill = c("#377eb8","black", "#e41a1c","#4daf4a"),
#        bty="n")
# dev.off()
# 
# ########################################################
# #### Add in trees into predicted locations
# #######################################################
# trees <- readOGR(dsn ="All_YFDP_NecessaryFiles_UTM10N", layer ="YFDP_Tree_UTM10N")
# trees[trees$DBH > 5] <- NA 
# trees_raster <- rasterize(trees, example_raster_2m, field="DBH", fun="min")
# trees_raster[trees_raster >= 0] <- 1
# trees_raster[is.na(trees_raster)] <- 0
# trees_raster[!in.plot] <- NA
# plot(is.na(trees_raster), main - "NA Values in Trees Raster")
# plot(trees_raster, main = "Tree Raster")
# trees_presAbs <- as.matrix(values(trees_raster))
# trees_presAbs <- trees_presAbs[all.complete]
# 
# presAbs.response.withtrees <- presAbs.response.cgs + trees_presAbs
# presAbs.response.withtrees[presAbs.response.withtrees > 1] <- 1
# 
# topo_presAbs_err_matrix_treesANDshrubs <- make_classError_matrix(NClasses=2, actual=presAbs.response.withtrees, predicted=topo_presAbs_predicted)
# topo_class_error_presAbs_treesANDshrubs <- class_error_presAbs(topo_presAbs_err_matrix)
# 
# lidar_presAbs_err_matrix_treesANDshrubs <- make_classError_matrix(NClasses=2, actual = presAbs.response.withtrees, predicted = LiDAR_presAbs_predicted)
# lidar_class_error_presAbs_treesANDshrubs <- class_error_presAbs(lidar_presAbs_err_matrix)
# 
# lidarANDtopo_presAbs_err_matrix_treesANDshrubs <- make_classError_matrix(NClasses=2, actual = presAbs.response.withtrees, predicted = topo_and_LiDAR_PresAbs_predicted)
# lidarANDtopo_class_error_presAbs_treesANDshrubs <- class_error_presAbs(lidarANDtopo_presAbs_err_matrix)
# 
# actual.shrubANDtrees.locations <- matrix(data = NA, nrow = length(all.complete))
# actual.shrubANDtrees.locations[all.complete] <- presAbs.response.withtrees
# actual.shrubANDtrees.locations_raster <- raster(actual.shrubANDtrees.locations, template = example_raster_2m)
# 
# #Topo
# jpeg("topo_AllShrubPredicted_withtrees.jpg", width = 12, height = 10, units = "in", res = 400)
# plot(topo.predicted.allshrub.locations_raster + (2*actual.shrubANDtrees.locations_raster), 
#      main = "Topography Predicted Locations of Shrubs and Trees", legend = F, 
#      col = c("#377eb8", "#e41a1c","black","#4daf4a"))
# legend("bottomleft", legend = c("Correct Absence",
#                                 "Error of Omission",
#                                 "Error of Commission", 
#                                 "Correct Presence"), 
#        fill = c("#377eb8","black", "#e41a1c","#4daf4a"),
#        bty="n")
# dev.off()
# 
# #LiDAR
# jpeg("LiDAR_AllShrubPredicted_withtrees.jpg", width = 12, height = 10, units = "in", res = 400)
# plot(LiDAR.predicted.allshrub.locations_raster + (2*actual.shrubANDtrees.locations_raster), 
#      main = "LiDAR Predicted Locations of Shrubs and Trees", legend = F, 
#      col = c("#377eb8", "#e41a1c","black","#4daf4a"))
# legend("bottomleft", legend = c("Correct Absence",
#                                 "Error of Omission",
#                                 "Error of Commission", 
#                                 "Correct Presence"), 
#        fill = c("#377eb8","black", "#e41a1c","#4daf4a"),
#        bty="n")
# dev.off()
# #Topo and LiDAR
# jpeg("topo_and_LiDAR_AllShrubPredicted_withtrees.jpg", width = 12, height = 10, units = "in", res = 400)
# plot(topo_and_LiDAR.predicted.allshrub.locations_raster + (2*actual.shrubANDtrees.locations_raster), 
#      main = "Topography and LiDAR Predicted Locations of Shrubs and Trees", legend = F, 
#      col = c("#377eb8", "#e41a1c","black","#4daf4a"))
# legend("bottomleft", legend = c("Correct Absence",
#                                 "Error of Omission",
#                                 "Error of Commission", 
#                                 "Correct Presence"), 
#        fill = c("#377eb8","black", "#e41a1c","#4daf4a"),
#        bty="n")
# dev.off()
# 
# ########################################################
# #### Add in trees, shrubs, and cwd into predicted locations
# #######################################################
# presAbs.response.withtreesANDcwd <- presAbs.response.cgs + trees_presAbs + cwd_presAbs
# presAbs.response.withtreesANDcwd[presAbs.response.withtreesANDcwd > 1] <- 1
# 
# topo_presAbs_err_matrix_treesANDshrubsANDcwd <- make_classError_matrix(NClasses=2, actual=presAbs.response.withtreesANDcwd, predicted=topo_presAbs_predicted)
# topo_class_error_presAbs_treesANDshrubsANDcwd <- class_error_presAbs(topo_presAbs_err_matrix)
# 
# lidar_presAbs_err_matrix_treesANDshrubsANDcwd <- make_classError_matrix(NClasses=2, actual = presAbs.response.withtreesANDcwd, predicted = LiDAR_presAbs_predicted)
# lidar_class_error_presAbs_treesANDshrubsANDcwd <- class_error_presAbs(lidar_presAbs_err_matrix)
# 
# lidarANDtopo_presAbs_err_matrix_treesANDshrubsANDcwd <- make_classError_matrix(NClasses=2, actual = presAbs.response.withtreesANDcwd, predicted = topo_and_LiDAR_PresAbs_predicted)
# lidarANDtopo_class_error_presAbs_treesANDshrubsANDcwd <- class_error_presAbs(lidarANDtopo_presAbs_err_matrix)
# 
# actual.shrubANDtreesANDcwd.locations <- matrix(data = NA, nrow = length(all.complete))
# actual.shrubANDtreesANDcwd.locations[all.complete] <- presAbs.response.withtreesANDcwd
# actual.shrubANDtreesANDcwd.locations_raster <- raster(actual.shrubANDtreesANDcwd.locations, template = example_raster_2m)
# 
# #Topo
# jpeg("topo_AllShrubPredicted_withTreesANDCWD.jpg", width = 12, height = 10, units = "in", res = 400)
# plot(topo.predicted.allshrub.locations_raster + (2*actual.shrubANDtreesANDcwd.locations_raster), 
#      main = "Topography Predicted Locations of Shrubs and Trees and CWD", legend = F, 
#      col = c("#377eb8", "#e41a1c","black","#4daf4a"))
# legend("bottomleft", legend = c("Correct Absence",
#                                 "Error of Omission",
#                                 "Error of Commission", 
#                                 "Correct Presence"), 
#        fill = c("#377eb8","black", "#e41a1c","#4daf4a"),
#        bty="n")
# dev.off()
# 
# #LiDAR
# jpeg("LiDAR_AllShrubPredicted_withTreesANDCWD.jpg", width = 12, height = 10, units = "in", res = 400)
# plot(LiDAR.predicted.allshrub.locations_raster + (2*actual.shrubANDtreesANDcwd.locations_raster), 
#      main = "LiDAR Predicted Locations of Shrubs and Trees and CWD", legend = F, 
#      col = c("#377eb8", "#e41a1c","black","#4daf4a"))
# legend("bottomleft", legend = c("Correct Absence",
#                                 "Error of Omission",
#                                 "Error of Commission", 
#                                 "Correct Presence"), 
#        fill = c("#377eb8","black", "#e41a1c","#4daf4a"),
#        bty="n")
# dev.off()
# #Topo and LiDAR
# jpeg("topo_and_LiDAR_AllShrubPredicted_withTreesANDCWD.jpg", width = 12, height = 10, units = "in", res = 400)
# plot(topo_and_LiDAR.predicted.allshrub.locations_raster + (2*actual.shrubANDtreesANDcwd.locations_raster), 
#      main = "Topography and LiDAR Predicted Locations of Shrubs and Trees and CWD", legend = F, 
#      col = c("#377eb8", "#e41a1c","black","#4daf4a"))
# legend("bottomleft", legend = c("Correct Absence",
#                                 "Error of Omission",
#                                 "Error of Commission", 
#                                 "Correct Presence"), 
#        fill = c("#377eb8","black", "#e41a1c","#4daf4a"),
#        bty="n")
# dev.off()

#########################################################
#### Verify Parsimonious Predictors
#########################################################
## Test sub
t_pres    <- which(presAbs.response.cgs == 1 & v == 0)
t_rand.not.pres <- sample(which(presAbs.response.cgs == 0 & v == 0), length(t_pres))
t_test.sub <- c(t_pres, t_rand.not.pres)

output <- file("v_new_gridmetrics_parsimonious_sep_training_and_testing.txt")
sink(output, append=TRUE)
sink(output, append=TRUE, type="message")

# This will echo all input and not truncate 150+ character lines...
################################
# Topo and LiDAR
# Commenting out because this one was ok
t_parsimonious_topoANDLiDAR_presAbs <- verifyParsimonious(predictors=all_predictors_complete, 
                                                        responses=presAbs.response.cgs, 
                                                        rf=v_presAbs_topoANDLiDAR_rf, 
                                                        subset=t_test.sub,
                                                        npredictors=10)
gc()
# Topo   
t_parsimonious_topo_presAbs <- verifyParsimonious(predictors=topo_complete, 
                                                responses=presAbs.response.cgs, 
                                                rf=v_presAbs_topo_rf, 
                                                subset=t_test.sub,
                                                npredictors=10)
gc()
# LiDAR
t_parsimonious_LiDAR_presAbs <- verifyParsimonious(predictors= LiDAR_complete, 
                                                 responses=presAbs.response.cgs, 
                                                 rf=v_presAbs_LiDAR_rf, 
                                                 subset=t_test.sub,
                                                 npredictors=10)

v_parsimonious_topoANDLiDAR_presAbs <- verifyParsimonious(predictors=all_predictors_complete, 
                                                          responses=presAbs.response.cgs, 
                                                          rf=v_presAbs_topoANDLiDAR_rf, 
                                                          subset=v_run.sub,
                                                          npredictors=10)
gc()
# Topo
v_parsimonious_topo_presAbs <- verifyParsimonious(predictors=topo_complete, 
                                                  responses=presAbs.response.cgs, 
                                                  rf=v_presAbs_topo_rf, 
                                                  subset=v_run.sub,
                                                  npredictors=10)
gc()
# LiDAR
v_parsimonious_LiDAR_presAbs <- verifyParsimonious(predictors= LiDAR_complete, 
                                                   responses=presAbs.response.cgs, 
                                                   rf=v_presAbs_LiDAR_rf, 
                                                   subset=v_run.sub,
                                                   npredictors=10)

# Restore output to console
sink() 
sink(type="message")
