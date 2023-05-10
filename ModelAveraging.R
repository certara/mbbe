# the term "model" refers the the source model, the models to be averaged. The term sample refers to the bootstrap sample of the data set
# and the model run against that bootstrap sample of the data set, and the set of Monte Carlo simulations based on the best model from the 
# models run against that bootstrap sample of the data set
rm(list = ls())
crash_value <- 999999
nmodels <- 3 # how many models in the model averaging?
ngroups <- 4 # number of treatment groups, i.e, 4 for ABBA
Reference.groups <- c(1,2) # which of the group is reference
Test.groups <- c(3,4)  # and which is test
numParallel <- 16 
samp.size <- 20 #  # of samples for bootstrap and Monte Carlo simulation
home.dir <- "c:/fda/mbbe/modelaveraging" # where to run the models,
model.source <- "u:/fda/mbbe/modelaveraging" # model control files will be in model.source\model1, model.source\model2 etc. 
# must be only 1 file, with .mod extension in the folder
nmfe.path <- "c:/nm744/util/nmfe74.bat"
delta.parms <- 0.1 # 10% different between parameters for identifiability testing
use_check_identifiable <- TRUE # whether to check
NCA.end.time <- 72 # end time for AUClast in NCA
rndseed <- 1 # random seed
use.simulation.data <- TRUE
simulation.data.path <- "U:\\fda\\mbbe\\Modelaveraging\\data_sim.csv"
library(dplyr)  
library(stringr) 
library(xml2)  
library(future)
library(PKNCA) 
library(DescTools)
library(PKNCA)
library(installr)
set.seed(rndseed)  
setwd(home.dir) 
#' get.block
#' Get a block (e.g., $DATA) from a control file with multiple lines
#' return block as a single character string
#'
#' @param stem Block title, e.g., $DATA to look for
#' @param control control file text 
#'
#' @examples
#' get.block("$DATA",c("$PROB test","$INPUT ...)) 
get.block <- function(stem,control){  
  
  tryCatch({ 
    nlines <- length(control)
    for(this.line in 1:nlines){
      # remove comments
      comment.pos <- gregexpr(";",control[this.line])[[1]][1]
      if(comment.pos > 0){
        control[this.line] <- str_trim(substr(control[this.line],1,comment.pos-1))
      }else{
        control[this.line] <- str_trim(control[this.line])
      }
    }
    start.line <- grep(paste0("^\\",stem),control)
    next.line <- grep("^\\$",control[start.line+1:length(control)])[1]
    if(is.empty(next.line)){# $stem is last
      last.line <- length(control)
    }else{
      last.line <- start.line+next.line
    }
    block <- ""
    for(line in start.line:(last.line-1)){
      block <- paste(block,control[line])
    } 
    return(block)
  },error = function(err) {  
    return(FALSE)
  }
  )
}
#' check.requirements
#' 1. Is there exactly 1 .mod file in each source folder
#' 2. Does that .mod file contain ;;;;.*Start EST
#' 3. is the sum of Reference.groups and Test.groups == ngroups?
#' 4. Any duplicates in Reference and Test groups?
#' 5. check if data file is present
#' 6.  check if nmfe path is correct
#' 7. check for saddle_reset if requested
#' 8. check for repeat IDs in data set 
#' 9. if use.simulation.data, see if data available
#' 
#' @param source.dir source  model files
#' @param nmodels number of models
#' @param ngroups number of groups in data set
#' @param Reference.groups Which groups are reference
#' @param Test.groups Which groups are Test
#' @param nmfe.path nmfe??.bat
#' @param use_check_identifiable - logical, whether to check for identifability, will check if SADDLE_RESET is in $EST
#' @param use.simulation.data logical, if the simulation will be done with a diferent data set than the bootstrap
#' @param simulation.data.path, if use.simulation.data, this is the path to the data file
#' @examples
#' check.requirements("c:/models",5,"c:/nm744/util",TRUE)) 
check.requirements <- function(source.dir,nmodels,ngroups,Reference.groups,Test.groups,nmfe.path,use_check_identifiable,use.simulation.data,simulation.data.path=NULL){
  
  msg <- list()
  result = tryCatch({ 
    if(!file.exists(nmfe.path)){ # if DOS path, convert to R/linux
      return(list(rval=FALSE,msg=paste("Cannot find nmfe?? at",nmfe.path,", exiting")))
    }
    # check number   in Reference and Test groups 
    if(sum(length(Reference.groups),length(Test.groups)) != ngroups){
      
      return(list(rval=FALSE,msg=paste("number of Reference groups ",length(Reference.groups),"+ Test groups ",length(Test.groups),
                                       "doesn't equal the number of groups", ngroups,", exiting")))
    }
    # no duplicated in Reference and Test groups 
    if(anyDuplicated(c(Reference.groups,Test.groups)) > 0) {
      return(list(rval=FALSE,msg=paste("There are duplicated group numbers between Reference and Test group, exiting")))
    }
    if(anyDuplicated(Reference.groups) > 0) {
      return(list(rval=FALSE,msg=paste("There are duplicated group numbers in the Reference group, exiting")))
    }
    if(anyDuplicated(Test.groups) > 0) {
      return(list(rval=FALSE,msg=paste("There are duplicated group numbers in the test group, exiting")))
    }
  },error = function(err) {  
    
    msg = list(rval=FALSE,msg=paste("Error in finding nfme??.bat,",nmfe.path,", exiting")) 
  }
  )
  for(this.model in 1:nmodels){
    result <- tryCatch({
      for(this.model in 1:nmodels){
        if(!file.exists(file.path(source.dir,paste0("model",this.model)))){
          return(list(rval=FALSE,msg=paste("Cannot find",file.path(source.dir,paste0("model",this.model)),", exiting")))
        }else{ 
          model.dir <- file.path(source.dir,paste0("model",this.model))
          modfiles <- list.files(model.dir,pattern = ".mod")
          if(length(modfiles)!=1){
            msg=append(list,list(rval=FALSE,msg=paste("There must be exactly 1 .mod file in",model.dir,", exiting")))
          } 
          # open file
          con <- file(file.path(model.dir,modfiles), "r")
          suppressWarnings(control <- readLines(con,encoding = "UTF-8"))
          #control <- readLines(con,encoding = "UTF-8")
          close(con)
          data.line <- get.block("$DATA",control)
          #data.line <- str_trim(control[grep("\\$DATA",control)]) 
          data.line <- str_trim(str_replace(data.line,"\\$DATA",""),side ="both")
          any.quotes = grep("^\"",data.line)
          if(length(any.quotes)>0){
            # find 2nd
            pos <- gregexpr(pattern ='\"',data.line)
            data.file <- substr(data.line, 1,pos[[1]][2])
          }else{
            # find first white space
            pos <- gregexpr(pattern ="\\s",data.line)
            data.file <- str_trim(substr(data.line, 1,pos[[1]][1]),side ="both")
          }
          if(!file.exists(data.file)){
            msg <- append(msg,list(rval=FALSE,msg=paste0("Cannot find ",data.file,", exiting")))
          }else{
            # repeat IDs??
            
            data <- read.csv(data.file)
            newIDs <- data %>% select(ID) %>% 
              mutate(newID=if_else(ID==lag(ID),0,1)) %>% 
              mutate(newID=if_else(is.na(newID),1,newID)) %>%
              filter(newID==1)  
            num_newIDs <- dim(newIDs)[1] 
            num_oldIDs <- dim(data %>% distinct(ID))[1]
            if(num_newIDs != num_oldIDs){
              msg = append(msg,list(rval=FALSE,msg=paste0("There appears to be repeat IDs in",data.file,
                                                          "\n, for bootstrap sampling IDs must not repeat, exiting")))
            }  
          }
          
          # $EST
          # must be on one line, need to fix this
          # contains ;;;; Start EST?
          contains_start = grepl(";;;;.*Start\\s+EST",control, ignore.case = TRUE)
          if(!any(contains_start)){
            msg=append(msg,list(rval=FALSE,msg=paste0("The control file ",file.path(model.dir,modfiles)," does not contain \";;;; Start EST\", required before the $EST and the $THETA records and after $OMEGA and $SIGMA")))
          }
          if(use_check_identifiable){
            EST.line <- get.block("$EST",control) 
            if(!grepl("SADDLE_RESET\\s*=\\s*1", EST.line, fixed = FALSE)){
              msg=append(msg,list(rval=FALSE,msg=paste0("Identifiability check requested, but SADDLE_RESET not set to 1 in ",file.path(model.dir,modfiles,", exiting"))))
            }      
          }
          if(use.simulation.data){
            if(!file.exists(simulation.data.path)){
              msg <- append(msg,list(rval=FALSE,msg=paste0("Cannot find simulation data file ",simulation.data.path,", exiting")))
            } 
          }
        }
      }
    },
    error = function(err) {  
      msg = append(msg,list(rval=FALSE,msg=paste("Error in number of model files in",model.dir,", exiting")))
    }
    )
    if(length(msg)>0){
      return(msg)
    }else{
      return(list(rval=TRUE,msg=paste("passes requirements check")))
    }
  }
}
#' split_path
#' split a path into parents
#' 
#' @param  path, string, path name to be split
#' @param mustWork logical
#' @examples
split_path <- function(path, mustWork = FALSE) {
  
  #' split_path("c:/modelaveragaing/runmodel/model1/1)
  output <- c(strsplit(dirname(normalizePath(path,mustWork = TRUE)),
                       "/|\\\\")[[1]], basename(path))
  return(output)
} 
#' copy.model.files
#' copy NONMEM model files from a source to a run directory
#' 
#' @param  model.source folder name were subfolders (model1..modelnmodels) with model files can be found. There must be exactly one .mod file (with any stem) in the folder
#' @param home.dir Folder where models are to be run
#' @param nmodels how many models are there
#' @examples
#' copy.model.files("c:/modelaveraging","c:/modelaveraging/run",4)
copy.model.files <- function(model.source,home.dir,nmodels){
  
  if(!file.exists(home.dir)){
    # split home.dir, check each parent
    normpath = SplitPath(home.dir)$normpath
    parents = split_path(normpath)
    nparents = length(parents)
    cur.path = parents[1]
    for(this.parent in 2:nparents){
      cur.path = file.path(cur.path,parents[this.parent])
      if(!file.exists(cur.path)){
        dir.create(cur.path)
      }
    }
  }
  for(this.model in 1:nmodels){
    sim.dir <- file.path(home.dir,paste0("model",this.model))
    if(dir.exists(sim.dir)){
      unlink(sim.dir,recursive = TRUE, force = TRUE)
    }
    dir.create(sim.dir) 
    source.dir <- file.path(model.source,paste0("model",this.model))
    modfile <- file.path(model.source,paste0("model",this.model),list.files(source.dir))
    file.copy(modfile,file.path(sim.dir,paste0("bs",this.model,".mod")))
  }
}

#' create boostrap samples of NONMEM data set, placed in home.dir, file names = data_repN.csv
#' 
#' @param home.dir Folder where models are to be run
#' @param samp.size how many bootstrap samples, default is 100
#' @param nmodels how many models are there
#' @examples
#' sample.data("c:/modelaveraging",100,4)
sample.data <- function(home.dir,nmodels,samp.size){
  
  con <- file(file.path(home.dir,"model1","bs1.mod"), "r")
  suppressWarnings(control <- readLines(con,encoding = "UTF-8"))
  close(con)
  data.line <- str_trim(control[grep("\\$DATA",control)]) 
  data.line <- str_trim(str_replace(data.line,"\\$DATA",""),side ="both")
  # if file name is quoted, just put out part in in quotes, otherwise get first white space
  any.quotes = grep("^\"",data.line)
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
  #  possibly different data files in different models? not supported at this time
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
      
      con <- file(file.path(home.dir,paste0("model",this.model),this.samp,paste0("bsSamp",this.samp,".mod")), "w")
      writeLines(control,con)
      close(con)
    }
  }
}

#' run the bootstrap models/samples
#' 
#' @param nmfe.path path to nmfe??.bat
#' @param home.dir Folder where models are to be run 
#' @param nmodels how many models are there
#' @param samp.size how many samples are there
#' @param numParallel how many to run in parallel, default = availableCores()
#' @examples
#' run.bootstrap("c:/nmfe744/util/nmfe74.bat", "c:/modelaveraging",8)
run.bootstrap <- function(nmfe.path,home.dir,nmodels,samp.size,numParallel=availableCores()){
  
  rval = tryCatch({
    plan(multisession,workers=numParallel)   
    this.model = 1
    this.samp = 1
    for(this.model in 1:nmodels){
      print(paste0("Starting first bootstrap sample model ",this.model,", Time = ",Sys.time()))
      for(this.samp in 1:samp.size){ 
        future({  
          setwd(file.path(home.dir,paste0("model",this.model),this.samp))
          command = paste0(nmfe.path," bsSamp",this.samp,".mod bsSamp",this.samp,".lst -nmexec=",paste0("bsSamp",this.model,"_",this.samp,".exe"))
          shell(command,wait = TRUE)
          delete_files(getwd())
        })
      }
      print(paste0("Staring last bootstrap sample model ",this.model,", Time = ",Sys.time()))
    }
    plan(sequential) 
    return(TRUE)
  },
  error = function(cond) {
    return(FALSE)
  }
  ) 
}

#' read the xml file for each boostrap/model
#'  return list of BICs (with the BEST column) and parameters for all models/samples
#'  
#' @param home.dir Folder where models are to be run 
#' @param nmodels how many models are there
#' @param samp.size how many samples are there 
#' @param delta.parms criteria for identifability test
#' @examples
#' get.parameters("c:/modelaveraging",4,100,0.1)
get.parameters <- function(home.dir,nmodels,samp.size,delta.parms){ 
  
  BICS <- data.frame(matrix(crash_value,nrow = samp.size,ncol=nmodels+1))
  colnames(BICS) <- c(paste0("Model",seq(1:nmodels)),"Best")
  # BIC=k*ln(n) -2LL where k is is the number of estimated parameters and n is the number of data
  nparms <- data.frame(ntheta=as.integer()) 
  
  # only need parameters for best model, but need to read xml file to get the number of observtions and parameter, so may
  # as well get THETA etc now for all
  this.model <- this.samp <- 1
  nfailed.ident <- 0 # number of samples in this model that fail the identifabilit test
  # do by sample first only save parameters for the best model
  selected.parameters = list() # parameters for selected model, 
  this.samp = 1
  for(this.samp in 1:samp.size){ 
    num_successful <- 1 # only save parameters if model finished, doesn't need to converge, just finish, this is the number of successful samples for this model
    parameters.this.sample <- list()
    this.model = 1
    for(this.model in 1:nmodels){ 
      # if using saddle_reset, test for identifiability
      if(!use_check_identifiable){
        identifiable_ok <- TRUE
      }else{
        identifiable_ok<- check_identifiable(home.dir,this.model,this.samp,delta.parms )
      }
      if(identifiable_ok){
        # first just get the number of parameters and observations for this model 
        xml_file <- file.path(home.dir,paste0("model",this.model),this.samp,paste0("bsSamp",this.samp,".xml"))
        if(!file.exists(xml_file)){
          BICS[this.samp,this.model] <- crash_value
          next
        }
        
        data <- read_xml(xml_file, encoding = "ASCII") 
        problem_node <- xml_find_all(data, "//nm:problem_information")
        contents <- xml_contents(problem_node)
        text <- xml_text(contents) 
        text <- as.list(unlist(strsplit(text, "\n"))) 
        nobsline <- grep("TOT. NO. OF OBS RECS:",text)
        nobs <- as.integer(str_replace(text[nobsline],"TOT. NO. OF OBS RECS:","")) 
        info_node <- xml_find_all(data,"//nm:problem_options")
        info_contents <- xml_attrs(info_node)
        ntheta <- as.numeric(info_contents[[1]]['nthetat']) 
        if(this.samp == 1){ # have to do this here, don't know how many thetas until here
          nparms <- rbind(nparms,data.frame(ntheta=ntheta)) 
          parameters.this.model <- data.frame(matrix(-9999,ncol=ntheta) ) 
          colnames(parameters.this.model) = paste0("THETA",seq(1:ntheta))
        }
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
          BICS[this.samp,this.model] = crash_value
        } 
        theta_node <- xml_find_all(data,"//nm:theta")
        theta_children <- xml_children(theta_node)
        theta_contents <- xml_contents(theta_node)
        theta <- as.numeric(xml_text(theta_contents))
        # length of theta will be 1 if a crash
        if(length(theta) > 1){
          parameters.this.model <- theta
          num_successful <- num_successful+ 1
        }   
      }else{ # fails identifiable test
        BICS[this.samp,this.model] <- crash_value
        nfailed.ident <- nfailed.ident + 1
      }
      # different models have differnt # of pameters
      
      # myList[[length(myList)+1]] <- list(sample(1:3))
      parameters.this.sample[[this.model]] = parameters.this.model
    } 
    # and select best model  
    best <- which.min(BICS[this.samp,1:nmodels])
    BICS$Best[this.samp] <- which.min(best)
    # copy selected model parameters 
    # myList[[length(myList)+1]] <- list(sample(1:3))
    selected.parameters[[this.samp]] = unlist(parameters.this.sample[best])
    
  }
  # write out results
  write.csv(BICS,file.path(home.dir,"BICS.csv"),quote=FALSE,row.names = FALSE)
  rval=list(BICS = BICS,parameters = selected.parameters)
  return(rval)
}

#' for each model get the control file used for the bootstrap
#' return a list of the models  
#' @param nmodels how many models are there 
#' @examples
#' get.base.model(4)
get.base.model <- function(nmodels){  
  
  
  # need error trapping for no ;;;;;.*Start EST 
  base.models <- vector(mode='list', length = 0) # all but $THETA, $SIM, $TABLE 
  
  for(this.model in 1:nmodels){
    con <- file(file.path(paste0("model",this.model),paste0("bs",this.model,".mod")), "r")
    suppressWarnings(control <- readLines(con,encoding = "UTF-8"))
    control <- control[1:grep(';;;;.*Start EST', control)]
    close(con)
    base.models <- append(base.models, list(control)) 
  }
  return(base.models)
  
  
}

#' Edits the best based model for each sample, replaces the original parameters with the bootstrap parameters
#' and adds the $SIM and $TABLE
#'  
#' @param home.dir Folder where models are to be run 
#' @param parms list that include BICs and parameters for each sample/model
#' @param base.models list of the text in each model used for the bootstrap
#' @param samp.size how many samples are there 
#' @examples
#' write.sim.controls("c:/modelaveraging",parms,base.models,100)
write.sim.controls <- function(home.dir,parms,base.models,samp.size,use.simulation.data,simulation.data.path=NULL){ 
  
  nmodels <- length(base.models)
  final.models <- vector(mode='list', length = samp.size) 
  model.indices <- rep(0,nmodels) # which parameter set to use when this model is selected, roll over if not enough samples,
  # shouldn't happen, as we shouldn't use all parameter sets for any model
  this.samp = 1
  for(this.samp in 1:samp.size){
    which.model <- parms$BICS$Best[this.samp]
    # use BS parameters, when you run out (as some BS samples fail), just start over, so need different random seed in $SIM
    # model is first index, e.g., parameters[[2]]$THETA1[1] is THETA[1] for model 2, first sample
    model.indices[which.model] <- model.indices[which.model] + 1
    # if not enough model parameters (because some model crashed?), recycle. But, shouldn't need to
    # as model that crashes should have crash_value BIC
    if(model.indices[which.model] > length(parms$parameters[[1]])[1]) model.indices[which.model] <- 1
    full.control <- base.models[which.model][[1]]
    # need to get sequence in here, calculated in $PK
    seed <- round(runif(1,0,10000),0)
    full.control <- c(full.control,paste("$SIM ONLYSIM (",seed,")"),"$THETA") 
    ntheta <- length(parms$parameters[[which.model]])
    this.parm = 1
    if(use.simulation.data){
      # get $DATA
      start.line <- grep("^\\$DATA",full.control)
      next.line <- grep("^\\$",full.control[start.line+1:length(full.control)])[1]
      if(is.empty(next.line)){ # $DATA is last
        last.line <- length(full.control[[1]])
      }else{
        last.line <- start.line+next.line
      }
      # replace with simulation data set
      full.control <-  c(full.control[1:(start.line-1)],
                         paste("$DATA",simulation.data.path,"IGNORE=@"),
                         full.control[(last.line):length(full.control)]) 
    }
    for(this.parm in 1:ntheta){
      full.control <- c(full.control,paste0(parms$parameters[[which.model]][this.parm],"  ;; THETA(",this.parm,")"))
    }  
    full.control <- c(full.control,"$TABLE ID TIME GROUP OCC SEQ DV EVID NOPRINT NOAPPEND FILE=OUT.DAT ONEHEADER")
    sim.dir <- file.path(home.dir,paste0("Sim",this.samp))
    if(dir.exists(sim.dir)){
      unlink(sim.dir,recursive = TRUE, force = TRUE)
    }
    dir.create(sim.dir) 
    con <- file(file.path(sim.dir,paste0("sim",this.samp,".mod")), "w")
    writeLines(unlist(full.control), con)
    close(con)
    final.models <- append(final.models,list(full.control))
  } 
  return(final.models) 
}

#' make array of all PIDs, check when all return Null 
#' Stem for bootstrap is bsSamp,  
#' Note this assume at least one run has started, could be aproblem if all runs
#' are still running nmtran/compiler
#' No time out on runs??
#' @param nmodels number of base models
#' @param samp.size number of bootstrap sample for each model
#' @examples
#' wait.for_bs(4,100) 
# names of runs will be stem(modelnum)_(samp_num).exe
wait.for.bs = function(nmodels, samp.size){
  
  PID <- vector("numeric",length=nmodels*samp.size)
  this.run <- 0
  for(this.model in 1:nmodels){ 
    for(this.samp in 1:samp.size){
      this.run <- this.run + 1
      PID[this.run] <- max(0,get_pid(paste0("bsSamp",this.model,"_",this.samp,".exe")))
    }
  }
  # loop until all are 0
  while(any(PID>0)){
    Sys.sleep(5)
    this.run <- 0
    for(this.model in 1:nmodels){ 
      for(this.samp in 1:samp.size){
        this.run <- this.run +1
        if(PID[this.run]>0){
          PID[this.run] <- max(0,get_pid(paste0("bsSamp",this.model,"_",this.samp,".exe")))
        }
      }
    }
  }
}

#' make array of all PIDs, check when all return Null 
#' Stem for simulation is sim
#' Note this assume at least one run has started, could be aproblem if all runs
#' are still running nmtran/compiler
#' No time out on runs??
#' @param samp.size number of bootstrap sample for each model
#' @examples
#' wait.for.sim(100) 
# names of runs will be sim(samp_num).exe
wait.for.sim = function(samp.size){
  
  
  PID <- vector("numeric",length=nmodels*samp.size)
  this.run <- 0 
  for(this.samp in 1:samp.size){
    this.run <- this.run + 1
    PID[this.run] <- max(0,get_pid(paste0("sim",this.samp,".exe")))
  } 
  # loop until all are 0
  while(any(PID>0)){
    Sys.sleep(5)
    this.run <- 0 
    for(this.samp in 1:samp.size){
      this.run <- this.run +1
      if(PID[this.run]>0){
        PID[this.run] <- max(0,get_pid(paste0("sim",this.samp,".exe"))) # get_pid returns NULL is no such process
        
      }
    }
  }
}
#' run all simulations
#' @param nmfe.path path to nmfe??.bat
#' @param home.dir home directory
#' @param numParallel number of models to run in parallel, default = avaiableCores()
#' @param samp.size number of bootstrap sample for each model
#' @examples
#' run.simulations("c:/nmfe75/util/nmfe75.bat","c:/mbbe",8,100)
run.simulations <- function(nmfe.path,home.dir,samp.size,numParallel=availableCores()){
  plan(multisession,workers=numParallel)  
  for(this.samp in 1:samp.size){ 
    future({  
      setwd(file.path(home.dir,paste0("sim",this.samp)))
      command <- paste0(nmfe.path," sim",this.samp,".mod sim",this.samp,".lst -nmexec=sim",this.samp,".exe") 
      shell(command,wait = TRUE)
      delete_files(getwd())
    }) 
  }
  
}

#' call getNCA for each sample in 1:samp.size
#' and do NCA (Cmax, AUCin,AUClst, from 0 to end.time)
#' if check.identifability, do that with delta.parms as criteria
#' @param home.dir, string home directory  
#' @param nGroups, integer how many groups, e.g., 4 for ABBA
#' @param Reference.groups, arrays which groups are reference  
#' @param Test.groups, array, which groups are  test
#' @param NCA.end.time, numeric, end time for AUClast and AUCinf
#' @param samp.size,numeric integer
#' @param delta.parms, numeric absolute difference in parameters criteria for failig identifability
#' @examples
#' calc.NCA("c:/runmodels",4,c(1,2),c(3,4)),72,100,0.1)
calc.NCA <- function(home.dir,ngroups,Reference.groups,Test.groups,NCA.end.time,samp.size){
  
  this.samp = 1
  for(this.samp in 1:samp.size){ 
    future({  
      # writes to file, future can't return a value?
      getNCA(home.dir,this.samp,ngroups,Reference.groups,Test.groups,NCA.end.time)
    })
  }
}
#' read .lst file, find out where saddle reset occurs then get parameter before reset from .ext file (Pre.parms)
#' compare with final parameters (Post.parms), see if any differ by delta.parms
#' open .lst 
#' @param home.dir home directory
#' @param this.model integer
#' @param this.sample integer
#' @param delta.parms absolute difference in parameters criteria for failig identifability
#' @examples
#' check_identifiable("c:/runmodels",1,1,0.1)
check_identifiable <- function(home.dir,this.model,this.sample,delta.parms){
  lstfile <- file.path(home.dir,paste0("model",this.model),this.sample,paste0("bsSamp",this.sample,".lst"))
  if(!file.exists(lstfile)){
    return(FALSE)  
  }else{
    tryCatch({
      con <- file(lstfile, "r")
      suppressWarnings(output <- readLines(con,encoding = "UTF-8"))
      close(con) 
      reset.line <- grep("^0SADDLE POINT RESET",output)
      
      # get previous "0ITERATION NO.: "
      first.output <- output[1:reset.line]
      first.output <- first.output[grep("^0ITERATION NO.:",first.output)]
      first.output <- first.output[length(first.output)] 
      reset.iteration <- as.integer(substr(first.output,16,24))
      last.output <- output[reset.line:length(output)]
      last.output <- last.output[grep("^0ITERATION NO.:",last.output)]
      last.output <- last.output[1] 
      last.iteration <-  as.integer(substr(last.output,16,24))
      # read parameters from .ext
      ext <- read.table(file.path(home.dir,paste0("model",this.model),this.sample,paste0("bsSamp",this.sample,".ext")),header=TRUE,skip=1)
      Pre.parms <- ext %>% filter(ITERATION == reset.iteration) 
      Post.parms <- ext %>% filter(ITERATION == last.iteration) 
      nparms <- length(Pre.parms)
      Passes.identifiability <- TRUE
      for(this.parm in 2:nparms){
        if(Pre.parms[this.parm] != 0){
          difference <- abs((Pre.parms[this.parm]-Post.parms[this.parm])/Pre.parms[this.parm])
          if(difference > delta.parms){ 
            Passes.identifiability <- FALSE
            break
          }
        }
      }
      return(Passes.identifiability)
    },
    error = function(cond) {
      return(FALSE)
    }
    )
  }
}

#' read $TABLE output from simulation (file name out.dat)
#' and do NCA (Cmax, AUCin,AUClst, from 0 to end.time)
#' if check.identifability, do that with delta.parms as criteria
#' @param home.dir home directory 
#' @param this.sample integer
#' @param NumGroups integer how many groups, e.g., 4 for ABBA
#' @param Reference.groups list of arrays, which groups are reference  
#' @param Test.groups list of arrays, which groups are  test
#' @param NCA.end.time end time for AUClast and AUCinf 
#' @examples
#' getNCA("c:/runmodels",1,4,c(1,2),c(3,4)),72,0.1)
getNCA = function(home.dir,this.sample,NumGroups,Reference.groups,Test.groups,NCA.end.time){
  group.NCA.results=  All.NCA.results = data.frame(ID=as.integer(),treatment = as.integer(),period=as.integer(),
                                                   sequence=as.integer(),Cmax=as.numeric(),AUCinf=as.numeric(),AUClast=as.numeric())
  data <- read.table(file.path(home.dir,paste0("sim",this.sample),"out.dat"),skip=1,header=TRUE) 
  data <- data %>% filter(EVID==0) %>% 
    filter(DV>0)   
  
  for(this.group in 1:NumGroups){
    group.data <- data %>% filter(GROUP==this.group) 
    period <- group.data %>% group_by(ID) %>% 
      distinct(ID,.keep_all=TRUE) %>% 
      select(ID,OCC,SEQ) %>%  arrange(ID)
    
    # insert conc=0 at time = 0
    zero.time <- group.data %>% distinct(ID,.keep_all = TRUE)
    zero.time$TIME <- 0
    zero.time$DV <- 0
    group.data <- rbind(group.data,zero.time) %>% 
      arrange(ID,TIME)
    conc_obj <- PKNCAconc(group.data, DV~TIME|ID)
    data_obj <- PKNCAdata(data.conc=conc_obj,
                          intervals=data.frame(start=0,
                                               end=NCA.end.time, 
                                               aucinf.obs=TRUE,
                                               auclast=TRUE,
                                               cmax=TRUE))
    results_obj <- pk.nca(data_obj)$result
    AUCinf <- results_obj %>% filter(PPTESTCD=="aucinf.obs") %>% 
      select(ID,PPORRES)
    AUClast <- results_obj %>% filter(PPTESTCD=="auclast") %>% 
      select(ID,PPORRES)
    CMAX <- results_obj %>% filter(PPTESTCD=="cmax") %>% 
      select(ID,PPORRES)
    if(this.group %in% Reference.groups){
      treatment <- "Reference"
    }else{
      treatment <- "Test"
    }
    group.NCA.results <- data.frame(ID=AUCinf$ID,treatment = treatment,period=period$OCC,
                                    sequence=period$SEQ,Cmax=CMAX$PPORRES,AUCinf=AUCinf$PPORRES,AUClast=AUClast$PPORRES)
    All.NCA.results <- rbind(All.NCA.results,group.NCA.results)
  }
  output.file <- file.path(home.dir,paste0("sim",this.sample),paste0("NCAresults",this.sample,".csv"))
  if(file.exists(output.file)){
    file.remove(output.file)
  }
  write.csv(All.NCA.results,file=output.file,quote=FALSE,row.names = FALSE) 
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
run.mbbe <- function(){ 
msg <- check.requirements(model.source,nmodels,ngroups,Reference.groups,Test.groups,nmfe.path,use_check_identifiable,simulation.data.path)
if(msg$rval){
  # return error code from each of these??
  copy.model.files(model.source,home.dir,nmodels)
  sample.data(home.dir,nmodels,samp.size)
  if(!run.bootstrap(nmfe.path,home.dir,nmodels,samp.size)){
    print("Failed bootstrap")
  }else{
    # need to wait until all are done, this returns when all are started.
    Sys.sleep(40) # in case all model are still compiling
    wait.for.bs(nmodels, samp.size)
    setwd(home.dir)
    parms <- get.parameters(home.dir,nmodels,samp.size,delta.parms)
    base.models <- get.base.model(nmodels) # get all nmodels base model
    final.models <- write.sim.controls(home.dir,parms,base.models,samp.size,use.simulation.data,simulation.data.path ) # don't really do anything with final models, already written to disc
    run.simulations(nmfe.path,home.dir,samp.size)
    Sys.sleep(40) # in case all model are still compiling, no executable yet
    wait.for.sim(samp.size)
    calc.NCA(home.dir,ngroups,Reference.groups,Test.groups,NCA.end.time,samp.size)
    # read NCA output and do stats
    BEsuccess <- data.frame(Cmax.success = as.logical(),AUCinf.success = as.logical(),AUClast.success = as.logical())
    for(this.model in 1:nmodels){
      for(this.samp in 1:samp.size){
      }
    }
  }
}else{
  print(msg$msg)
}
plan(sequential)
}
