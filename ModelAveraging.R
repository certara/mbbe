set.seed(1)
setwd("c:/fda/mbbe/modelaveraging")
rm(list = ls())
source("funcs.R")
library(xpose)
library(dplyr)  
library(stringr)
library(tidyr)
library(xml2) 
library(reshape2)
library(plyr)
set.seed(1) 
nmodels <- 5
numParallel <- 16 # must be multiple of nmodels
samp.size <- 100 #  # of samples for bootstrap
home.dir <- "c:/fda/mbbe/modelaveraging" 
model.source <- "u:/fda/mbbe/modelaveraging"
# copy source .mod files to working dir
for(this.model in 1:nmodels){
  sim.dir = file.path(home.dir,paste0("model",this.model))
  if(dir.exists(sim.dir)){
    unlink(sim.dir,recursive = TRUE, force = TRUE)
  }
  dir.create(sim.dir) 
  source.dir <- file.path(model.source,paste0("model",this.model))
  modfile <- file.path(model.source,paste0("model",this.model),list.files(source.dir))
  file.copy(modfile,file.path(sim.dir,paste0("bs",this.model,".mod")))
}
# get data file name
con <- file(file.path(home.dir,"model1","bs1.mod"), "r")
suppressWarnings(control <- readLines(con,encoding = "UTF-8"))
data.line <- str_trim(control[grep("\\$DATA",control)])
close(con)
data.line <- str_trim(str_replace(data.line,"\\$DATA",""),side ="both")
    # if file name is quoted, just put out part in in quotes, otherwise get first white space
any.quotes= grep("^\"",data.line)
if(length(any.quotes)>0){
    # find 2nd 
    pos <- gregexpr(pattern ='\"',data.line)
    data.file <- substr(data.line, 1,pos[[1]][2]) 
    rest.of.data.line <- substr(data.line,pos[[1]][2],nchar(data.line))
  }else{
    # find first white space
    pos <- gregexpr(pattern ="\\s",data.line)
    data.file  <- str_trim(substr(data.line, 1,pos[[1]][1]),side ="both")
    rest.of.data.line = substr(data.line,pos[[1]][1],nchar(data.line))
  } 
  # get datafile
  # do my model, possibly different data files in different models
org.data = read.csv(data.file,header=TRUE)
  # get IDs
cols <- colnames(org.data)
num.data.items <- length(cols)
  #cannot have repeat IDs!!!!
idID <- match(cols,"ID")
IDcol <- which(idID %in% 1) 
IDs <-  org.data[IDcol] %>% distinct(ID)
nsubs <- dim(IDs)[1]
  # create BS data sets
for (this.samp in 1:samp.size){
    who <- sample(IDs$ID,nsubs,replace=TRUE)
    this.data <- data.frame(matrix(-99,nrow=0,ncol=num.data.items))
    colnames(this.data) <- cols
    for(this.rep in 1:nsubs){ # need to be sure not to have adjacent with same ID, just just number IDs sequentially
      next_sub <- org.data %>% filter(ID==who[this.rep]) %>% mutate(ID=this.rep)
      this.data <- rbind(this.data,next_sub )
    } 
    # write data, will be in parent directory of run directory
    datafile.name <- paste0("data_rep",this.samp,".csv")
    write.csv(this.data,datafile.name,quote = FALSE,row.names = FALSE)
    # replace $DATA for each control, each model
     
    for(this.model in 1:nmodels){
      con <- file(file.path(home.dir,paste0("model",this.model),paste0("bs",this.model,".mod")), "r")
      suppressWarnings(control <- readLines(con,encoding = "UTF-8"))
      close(con) 
      data.line =  grep("\\$DATA",control)
      # get rest of data line
      # !!! note, will need to read entire $DATA block, maybe on more than one line!!!! 
      newdata.line = paste0("$DATA ",file.path("..","..",paste0("data_rep",this.samp,".csv ", rest.of.data.line)))
      control[data.line] =  newdata.line
      # remove covariance
      cov.line =  grep("\\$COV",control)
      if(length(cov.line)>0) control[cov.line] = ";; no covariance step"
      # and any tables
      control = str_replace(control,"$TABLE",";$TABLE")
      dir.create(file.path(home.dir,paste0("model",this.model),this.samp))
                 
      con <- file(file.path(home.dir,paste0("model",this.model),this.samp,paste0("bsrep",this.samp,".mod")), "w")
      writeLines(control,con)
      close(con)
    }
  }
# and run bootstrap
library(future)
library(promises) 
plan(multisession,workers=numParallel)   
for(this.model in 1:nmodels){
    print(paste("Start",this.model,", Time  =",Sys.time()))
    for(this.samp in 1:samp.size){ 
    future({  
      setwd(file.path(home.dir,paste0("model",this.model),this.samp))
      command = paste0("nmfe74 bsrep",this.samp,".mod bsrep",this.samp,".lst -nmexec=",paste0("bsrep",this.model,"_",this.samp,".exe"))
      shell(command,wait = TRUE)
    })
 }
}

plan(sequential)
setwd(home.dir)
# need to wait until all are done??? xml note written??
BICS <- data.frame(matrix(-99999,nrow = samp.size,ncol=nmodels+1))
colnames(BICS) <- c(paste0("Model",seq(1:nmodels)),"Best")
# BIC=k*ln(n) -2LL where is is the number of estimated parameters and n is the number of data
nparms <- data.frame(ntheta=as.integer()) 
parameters = list()    
this.model = 1
this.samp = 1 
for(this.model in 1:nmodels){ 
      # readxml  
    num_successful = 1 # only save parameters if model finished, doesn't need to converge, just finish
    for(this.samp in 1:samp.size){
      # first just get the number of parameters and observations for this model
       xml_file <- file.path(home.dir,paste0("model",this.model),this.samp,paste0("bsrep",this.samp,".xml"))
       data <- read_xml(xml_file, encoding = "ASCII") 
       problem_node <- xml_find_all(data, "//nm:problem_information")
       contents <- xml_contents(problem_node)
       text <- xml_text(contents) 
       text = as.list(unlist(strsplit(text, "\n"))) 
       nobsline = grep("TOT. NO. OF OBS RECS:",text)
       nobs <- as.integer(str_replace(text[nobsline],"TOT. NO. OF OBS RECS:","")) 
       info_node <- xml_find_all(data,"//nm:problem_options")
       info_contents <- xml_attrs(info_node)
       ntheta <- as.numeric(info_contents[[1]]['nthetat']) 
       if(this.samp == 1){
         nparms = rbind(nparms,data.frame(ntheta=ntheta)) 
         parameters.this.model =  data.frame(matrix(-9999,nrow=samp.size,ncol=ntheta) ) 
         colnames(parameters.this.model) = paste0("THETA",seq(1:ntheta))}
       #TOT. NO. OF OBS RECS:     6018
       #0LENGTH OF THETA:  27 
       # and OFV
       estim.node  <- xml_find_all(data,"//nm:estimation")
       estim_contents <- xml_attrs(estim.node)
         #nm:final_objective_function
       OFV.node <-  xml_find_all(data,"//nm:final_objective_function")
       OFV_contents <- xml_contents(OFV.node)
       OFV <- as.numeric(xml_text(OFV_contents))
       if(length(OFV) > 0) {
         BICS[this.samp,this.model] <- ntheta*log(nobs) + OFV
       } else {
         BICS[this.samp,this.model] = 999999999
       } 
      theta_node <- xml_find_all(data,"//nm:theta")
      theta_children <- xml_children(theta_node)
      theta_contents <- xml_contents(theta_node)
      theta <- as.numeric(xml_text(theta_contents))
        # length of theta will be 1 if a crash
      if(length(theta) > 1){
          parameters.this.model[num_successful,]  <- theta
          num_successful = num_successful+ 1
        } 
    } 
    parameters <- append(parameters,list(parameters.this.model))
}
# and select best model  
for(this.samp in 1:samp.size){
  BICS$Best[this.samp] <- which.min(BICS[this.samp,1:nmodels])
}
# and the simulation
# create the control files
base.models <- vector(mode='list', length = 0) # all but $THETA, $SIM, $TABLE
final.models <- vector(mode='list', length =0) 
this.model = 1
for(this.model in 1:nmodels){
  con <- file(file.path(paste0("model",this.model),paste0("bs",this.model,".mod")), "r")
  suppressWarnings(control <- readLines(con,encoding = "UTF-8"))
  control <- control[1:grep(';;;;.*Start EST', control)]
  close(con)
  base.models <- append(base.models, list(control))
   # select model base on BEST
    # add $SIM $TABLE and $THETAs  
}
model.indices <- rep(0,nmodels) # which parameter set to use when this model is selected, roll over if not enough samples,
                                # shouldn't happen, as we shouldn't use all parameter sets for any model
this.samp = 1
for(this.samp in 1:samp.size){
  which.model <- BICS$Best[this.samp]
  # use BS parameters, when you run out (as some BS samples fail), just start over, so need different random seed in $SIM
  # model is first index, e.g., parameters[[2]]$THETA1[1] is THETA[1] for model 2, first sample
  model.indices[which.model] <- model.indices[which.model] + 1
  if(model.indices[which.model] > dim(parameters[[1]])[1]) model.indices[which.model] <- 1
  full.control <- base.models[which.model]
  # need to get sequence in here, calculated in $PK
  seed <- round(runif(1,0,10000),0)
  full.control <- c(control,paste("$SIM ONLYSIM (",seed,")"),"$THETA") 
  for(this.parm in 1:nparms$ntheta[this.model]){
    full.control <- c(full.control,paste0(parameters[[which.model]][model.indices[which.model],][this.parm],"  ;; THETA",this.parm))
  }  
  full.control <- c(full.control,"$TABLE ID TIME GROUP OCC SEQUENCE DV EVID NOPRINT NOAPPEND FILE=OUT.DAT ONEHEADER")
  sim.dir <- file.path(home.dir,paste0("Sim",this.samp))
  if(dir.exists(sim.dir)){
    unlink(sim.dir,recursive = TRUE, force = TRUE)
    }
  dir.create(sim.dir) 
  con <- file(file.path(sim.dir,"sim1.mod"), "w")
  writeLines(full.control, con)
  close(con)
  final.models <- append(final.models,list(full.control))
  } 

# simulations 
library(future)
library(promises)
library(PKNCA)
plan(multisession,workers=4) 
AUCs  <- vector("list", nsamps) # list of data frames with AUCs for each sample, best model only
Cmaxs <- vector("list", nsamps) # list of data frames with CMAXs for each sample, best model only
source("NCA.r")
for(i in 1:samp.size){ 
  future({  
    setwd(file.path(home.dir,paste0("sim",i)))
    shell("nmfe74 sim1.mod sim1.lst",wait = TRUE)
  })
}
source(file.path(home.dir,"NCA.R"))
i = for(i in 1:nsamps){ 
  future({  
    getNCA(home.dir,i)
  })
}
print("OK")

plan(sequential)
 