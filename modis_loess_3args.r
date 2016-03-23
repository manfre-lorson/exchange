#
# loess regession on Modis ndvi with Quality
# calc on parts of the data
# read e.g. fiste 10 lines of each image over the full timeseries
# do the regession and write the result in a new hdf file

##--- librarys:
library(rhdf5)
library(rgdal)
library(raster)

source("O:/GIS_data/_r_script/RSTools.R")

strtTi <- Sys.time()
args = commandArgs(trailingOnly=TRUE)
##--------------------------------------------------------------------

##--- settings
tile <- args[1]
#tile <- 'h11v10'
root <- 'R:/Florian_Wolf/' #"D:/Florian_Wolf/hdf2h5/" #working path
files <- list.files(paste(root,"hdf5_",tile,sep="")) #hdf list in folder
files <- unique(substr(files,1,44))
op <- paste(root,"hdf5_",tile, "_QA", "/",sep="") #output path
##--------------------------------------------------------------------


##---------------------------------------------------------------------
##Function definitions

##-----------------------------------------
##def the loess functions on single pixel
vector.loess <- function(xi, QAi, n=11, b=0.5, s=0.1, max.na=5){
  
  if (max(na.consecutive(xi), na.rm = T) >= max.na) {return(rep(NA,length(xi)))}
  else
  {
    ti <-1:length(xi)
    wi <- 1/((b*QAi)+1)
    span = n/length(xi)
    loess.fit1 <- loess(xi ~ ti,
                        weight = wi,
                        span = span, 
                        control=loess.control(surface='direct'))
    
    loess.pred1 <- predict(loess.fit1, newdata=ti)
    
    diff <- xi - loess.pred1
    wi.2 <- wi/(1+(abs(diff)/(s*sd(diff, na.rm=T))))  
    wi.2[which(diff>0)] <- wi[which(diff>0)]
    
    loess.fit2 <- loess(xi ~ ti,
                        weight = wi.2,
                        span = span, 
                        control=loess.control(surface='direct')) 
    loess.pred2 <- as.integer(predict(loess.fit2, newdata=ti))
    #loess.pred2[loess.pred2 > 9999] <- xi[loess.pred2 > 9999]
    return(loess.pred2)
    #return(loess.pred2)
  }
}

##function to apply the loess function on the matrix
matrix.loess <- function(x, QA, n=11, b=0.5, s=0.1){
  t(sapply(1:NROW(x),       
           FUN=function(ii){
             tryCatch(vector.loess(x[ii,], QA[ii,], n=n, b=b, s=s),
                      error=function(err) err=as.numeric(rep(NA, length(x[ii,]))))
           }))
}


##get Quality as Number from the qa-layer
qaWeightFast <- function(x){
  bitwShiftR(x,2)%%16
}

## Filters
modland_qa.fun <- function(x){
  x%%4 #  modland_qa <- modland_qa.fun(qa) ; ndvi[modland_qa>1] <- NA
}

aerosol_qa.fun <- function(x){
  bitwShiftR(x,6)%%4  # aerosol_qa <- areosol_qa.fun(qa) ; ndvi[areosol_qa>2] <- NA
}

adjacent_cloud.fun <- function(x){
  bitwShiftR(x,8)%%2  # adjacent_cloud <- adjacent_cloud.fun(qa) ; ndvi[adjacent_cloud>0] <- NA
}

mixed_clouds.fun <- function(x){
  bitwShiftR(x,10)%%2 # mixed_clouds <- mixed_clouds.fun(qa); ndvi[mixed_clouds>0] <- NA
}

land_water.fun <- function(x){
  bitwShiftR(x,11)%%8 # land_water <- land_water.fun(qa); ndvi[land_water!=1] <- NA
}

possible_shadow.fun <- function(x){
  bitwShiftR(x,15)%%2 # possible_shadows <- possible_shadows.fun(qa); ndvi[possible_shadows>0] <- NA
}



# Find number of consecutive NA
na.consecutive <-function(x) is.na(x)*sequence(rle(is.na(x))$lengths)




lineWrite <- function(output,l,op,lineStart, lineStopp, name){
  h5write(as.matrix(output[ ,l]),file=paste(op,opNameList[l],sep=''),name=name,index=list(NULL,seq(lineStart,lineStopp)))
  H5close()
}

## end Definitions -------------------
##----------------------------------------------------------



opNameList <- list.files(paste(root,'hdf5_',tile,"_QA", sep=''))

com.args1 <- as.integer(args[2])
com.args2 <- as.integer(args[3])
#com.args1 <- 2221
#com.args2 <- 2221
#p=2221
for (p in seq(com.args1,com.args2,1)){ 
  strtTi <- Sys.time()
  lineStart <- p
  lineStopp <- p
  
  ## create the ndvi and qa matrix [1:4800,1:362]
  for (i in 1:length(files)){
    file <- files[i]
    h5 <- paste(root,"hdf5_",tile,"/",file, sep="" )
    meta <- h5ls(h5)
    i.ndvi <- h5read(h5 ,paste(meta[4,1],meta[6,2], sep="/"), index=list(NULL,seq(lineStart,lineStopp)))
    i.qa <- h5read(h5, paste(meta[4,1],meta[8,2], sep="/"), index= list(NULL,seq(lineStart,lineStopp)))
    if (i == 1) {ndvi <- i.ndvi
    qa <- i.qa} 
    else {
      ndvi <- cbind(ndvi, i.ndvi)
      qa <- cbind(qa, i.qa)      
    }
    
  }
  
  ndvi.orig <- ndvi
  qa.orig <- qa
  
  #filter_1
  qa.f1 <- matrix(qaWeightFast(qa.orig), ncol=ncol(qa))
  landwater <- matrix(land_water.fun(qa), ncol=ncol(qa))
  ndvi.f1 <- ndvi.orig
  ndvi.f1[qa.f1 > 12] <- NA
  ndvi.f1[landwater!= 1] <- NA
  ndvi.f1[ndvi.f1== -3000] <-NA
  ndvi.f1[ndvi.f1<0]<- 0

  ## do the regession
  output <- matrix.loess(ndvi.f1, qa.f1)
  output[is.na(output)]<- -3000
  #output[output>9999] <- ndvi.f1[output > 9999]
  
  
  
  
  ## write the calc data to the hdf
  for (l in seq(1,dim(output)[2])){
    #for (f in seq(1,3)){
    er <- tryCatch(lineWrite(l=l,output=output,op=op,lineStart = lineStart, lineStopp = lineStopp, name='qaModis'),error=function(e) e, warning=function(w) w)
    H5close()
    if(is(er,"warning")) {
      lineWrite(l=l,output=output,op=op,lineStart = lineStart, lineStopp = lineStopp, name='qaModis' )
      H5close()
      #lineWrite(l=l,output=output,lineStart = lineStart, lineStopp = lineStopp)
      #lineWrite(l=l,output=output.f1[,l],op=op,lineStart = lineStart, lineStopp = lineStopp, name=filters[[f]] )
      print('there was a warning')}
  }#}
  print(paste("line:",p,";  duration:", Sys.time()-strtTi))
}

print("done")
