rm(list = ls())
samp_size <- 10
crash_value <- 10000
residual_prop_error <- 0.2
residual_add_error <- 1
# source('u:/fda/tacrolimus/search_vcl/func.r')
defaultW <- getOption("warn")
options(warn = -1)
suppressMessages(library(DescTools))
suppressMessages(library(dplyr))
home.dir <- getwd()
folder <- tail(strsplit(home.dir, "/")[[1]], 2)
file_name_root <- paste0("nm_", folder[1], "_", folder[2])
output_file <- file.path(home.dir, paste0(file_name_root, ".lst"))
control_file <- file.path(home.dir, paste0(file_name_root, ".mod"))
xml_file <- file.path(home.dir, paste0(file_name_root, ".xml"))
nbins <- ceiling((log(samp_size)) * 2)
set.seed(1)
suppressMessages(library(transport))
Get_parms <- function(xml_file) {
  count <- 0
  while (!file.exists(xml_file) & count < 20) {
    count <- count + 1
    Sys.sleep(0.5)
  }
  if (!file.exists(xml_file)) {
    return(-1)
  }
  suppressMessages(library(xml2))
  data <- read_xml(xml_file, encoding = "ASCII")
  covariance_node <- xml_find_all(data, "//nm:covariance")
  if (length(covariance_node) == 0) {
    return(-2)
  }
  covariance_children <- xml_children(covariance_node)
  info_node <- xml_find_all(data, "//nm:problem_options")
  info_contents <- xml_attrs(info_node)
  ntheta <- as.numeric(info_contents[[1]]["nthetat"])
  theta <- vector("numeric", length = ntheta)
  theta_node <- xml_find_all(data, "//nm:theta")
  theta_children <- xml_children(theta_node)
  theta_contents <- xml_contents(theta_node)
  theta <- as.numeric(xml_text(theta_contents))
  # all omegas are 1
  covar_matrix <- matrix(NA, ncol = ntheta, nrow = ntheta)
  for (i in 1:ntheta) {
    node <- covariance_children[i]
    names <- xml_attrs(node)
    contents <- xml_contents(node)
    text <- xml_text(contents)
    for (n in 1:i) {
      covar_matrix[i, n] <- covar_matrix[n, i] <- as.numeric(text[n])
    }
  }
  return(list(THETA = theta, SIGMA = covar_matrix))
}


CalcNCAVec <- function(data) {
  # first sort, just to make is all lines up, will need for cases when Tmax == 24, then nothing in logauc data
  # create global ID
  tryCatch(
    {
      NCAparms <- data %>%
        group_by(ID, GROUP, SIM) %>%
        filter(TIME <= 24) %>%
        mutate(CMAX = max(DV), LastTime = lag(TIME), LastConc = lag(DV)) %>%
        rowwise() %>%
        mutate(pAUC = ifelse(anyNA(c(LastTime, LastConc, TIME, DV)), 0, ifelse(LastConc <= DV | DV == 0, (TIME -
          LastTime) * (DV + LastConc) / 2, (DV - LastConc) / (log(DV / LastConc)) * (TIME - LastTime)))) %>%
        group_by(ID, GROUP, SIM) %>%
        summarise(AUC = sum(pAUC), CMAX = CMAX, .groups = "keep")
    },
    error = function(e) {
      return("Error in CalcNCAVec")
    }
  )
  # remove any that are NA
  NCAparms <- NCAparms %>%
    filter(!is.na(CMAX), !is.na(AUC)) %>%
    distinct(ID, SIM, .keep_all = TRUE)
  return(NCAparms)
}



Get_TrueNCA <- function(org.data.file.path) {
  tryCatch(
    {
      # READ IN ORGINAL DATA:

      org.data <- read.table(org.data.file.path, skip = 1, header = 1)
      org.data$SIM <- 1
      reference1.data <- org.data %>%
        filter(GROUP == 1)
      reference2.data <- org.data %>%
        filter(GROUP == 2)
      test1.data <- org.data %>%
        filter(GROUP == 3)
      test2.data <- org.data %>%
        filter(GROUP == 4)

      # get distribution of 'real' reference and test data tacrolimus ABBA design, will match geomean of Cmax and AUC
      # =IF(G2='R1',1,IF(G2='R2',2,IF(G2='T1',3,IF(G2='T2',4)))) rename REALOBS to DV
      colnames(reference1.data) <- colnames(reference2.data) <- colnames(test1.data) <- colnames(test2.data) <- c(
        "ID",
        "TIME", "EVID", "DV", "TRT", "GROUP", "SIM"
      )
      true.Reference1.parms <- CalcNCAVec(reference1.data) # %>% dplyr::select(ID,TIME,DV))
      true.Reference2.parms <- CalcNCAVec(reference2.data) # %>% dplyr::select(ID,TIME,DV))
      true.test1.parms <- CalcNCAVec(test1.data) # %>% dplyr::select(ID,TIME,DV))
      true.test2.parms <- CalcNCAVec(test2.data) # %>% dplyr::select(ID,TIME,DV))
      # bootstrap true values bootstrap should be 'best case', with perfect precision of parameter, as this is just a
      # single sample from the distribution any results with parameter uncertainty should be worse (higher variance)
      # but mean may (will) be different as well

      NSubs <- length(true.Reference1.parms$CMAX)
      Ref1.BSParms <- data.frame(matrix(NA, nrow = samp_size, ncol = 3))
      colnames(Ref1.BSParms) <- c("Sample", "Cmax", "AUC")
      Ref2.BSParms <- data.frame(matrix(NA, nrow = samp_size, ncol = 3))
      colnames(Ref2.BSParms) <- c("Sample", "Cmax", "AUC")
      test1.BSParms <- data.frame(matrix(NA, nrow = samp_size, ncol = 3))
      colnames(test1.BSParms) <- c("Sample", "Cmax", "AUC")
      test2.BSParms <- data.frame(matrix(NA, nrow = samp_size, ncol = 3))
      colnames(test2.BSParms) <- c("Sample", "Cmax", "AUC")

      for (i in 1:samp_size) {
        samp <- true.Reference1.parms[sample(nrow(true.Reference1.parms), NSubs, replace = TRUE), ]
        Ref1.BSParms$Cmax[i] <- as.numeric(exp(mean(log(samp$CMAX))))
        Ref1.BSParms$AUC[i] <- as.numeric(exp(mean(log(samp$AUC))))
        Ref1.BSParms$Sample[i] <- i

        samp <- true.Reference2.parms[sample(nrow(true.Reference2.parms), NSubs, replace = TRUE), ]
        Ref2.BSParms$Cmax[i] <- as.numeric(exp(mean(log(samp$CMAX))))
        Ref2.BSParms$AUC[i] <- as.numeric(exp(mean(log(samp$AUC))))
        Ref2.BSParms$Sample[i] <- i

        samp <- true.test1.parms[sample(nrow(true.test1.parms), NSubs, replace = TRUE), ]
        test1.BSParms$Cmax[i] <- as.numeric(exp(mean(log(samp$CMAX))))
        test1.BSParms$AUC[i] <- as.numeric(exp(mean(log(samp$AUC))))
        test1.BSParms$Sample[i] <- i

        samp <- true.test2.parms[sample(nrow(true.test2.parms), NSubs, replace = TRUE), ]
        test2.BSParms$Cmax[i] <- as.numeric(exp(mean(log(samp$CMAX))))
        test2.BSParms$AUC[i] <- as.numeric(exp(mean(log(samp$AUC))))
        test2.BSParms$Sample[i] <- i
      }
    },
    error = function(e) {
      return("Error in Get_TrueNCA ")
    }
  )
  return(list(Ref1 = Ref1.BSParms, Ref2 = Ref2.BSParms, test1 = test1.BSParms, test2 = test2.BSParms))
}

unlinking_dir <- function(filePath) {
  tryCatch(
    {
      count <- 0
      while (dir.exists(filePath)) {
        unlink(filePath, recursive = TRUE, force = TRUE)
        count <- count + 1
        if (count > 20) {
          stop("Cannot delete ", filePath)
        }
        Sys.sleep(0.5)
      }
    },
    error = function(e) {
      return("Error in unlinking_dir")
    }
  )

  return(TRUE)
}

Get_SimNCA <- function(xml_file) {
  suppressMessages(library(matrixcalc))
  suppressMessages(library(MASS))
  suppressMessages(library(stringr))
  suppressMessages(library(stringi))
  unlinking_dir(file.path(getwd(), "sim"))
  parms <- Get_parms(xml_file)
  if (!is.list(parms)) {
    if (parms == -1) {
      return(-1)
    }
    if (parms == -2) {
      return(-2)
    }
  }
  # set residual error to 20%
  # need to match THETAs in template.txt
  parms$THETA[8] <- log(residual_prop_error^2) # 20% residual CV, DEFINED IN SEARCH_SPPC.R
  parms$THETA[9] <- log(residual_add_error^2) # small
  sim.parms <- data.frame(mvrnorm(samp_size, parms$THETA, parms$SIGMA))
  sim.parms <- unname(sim.parms)
  con <- file(control_file, "r")
  org.control <- readLines(con, encoding = "UTF-8")
  close(con)
  # get code without THETA All.but.THETA <- org.control[1:grep('^;;;;.*Start THETA', org.control)] # only need part
  # before this, code without sub
  seeds <- as.integer(runif(samp_size, 0, 10) * 10000)
  base.code <- org.control[1:grep(";;;;.*Start subs", org.control)] # only need part before this,
  # code with $EST, replace with $SIM and $Table out?.dat
  Initial.code <- org.control[1:grep(";;.*Start EST", org.control)]
  sim.control <- c(
    Initial.code, paste0(
      "$TABLE  ID TIME EVID REALOBS TRT GROUP NOAPPEND FILE=ORG.DAT NOPRINT  ONEHEADER",
      "\n$SIM ONLYSIM (", seeds[1], ") \n", "$TABLE ID TIME EVID DV TRT GROUP NOAPPEND FILE=OUT1.DAT NOPRINT ONEHEADER"
    ),
    "\n$THETA", unlist(sim.parms[1, ])
  )

  for (i in 2:samp_size) {
    THETA.Block <- unlist(c("$THETA", sim.parms[i, ]))
    sim.control <- c(sim.control, base.code, THETA.Block, paste0(
      "$SIM ONLYSIM (", seeds[i], ") \n", "$TABLE ID TIME EVID DV TRT GROUP NOAPPEND FILE=OUT",
      i, ".DAT NOPRINT ONEHEADER"
    ))
  }

  dir.create("sim")
  # file.copy(file.path(home.dir,data.file),file.path(home.dir,'sim',data.file))
  con <- file(file.path(home.dir, "sim", "sim1.mod"), "w")
  writeLines(sim.control, con)
  close(con)

  nmoutput <- processx::run("nmfe74.bat", args = c("sim1.mod", "sim1.lst"), wd = file.path(home.dir, "sim"))

  raw_data <- data.frame()

  tryCatch(
    {
      for (this.sim in 1:samp_size) {
        file_name <- file.path(home.dir, "sim", paste0("out", this.sim, ".dat"))
        count <- 0
        while (!file.exists(file_name) & count < 10) {
          count <- count + 1
          Sys.sleep(0.2)
        }
        if (file.exists(file_name)) {
          data <- read.table(file_name, skip = 1, header = TRUE)
          data <- data %>%
            filter(DV >= 0)
          data$SIM <- this.sim
          raw_data <- rbind(raw_data, data)
        } else {
          return(-3)
        }
      }
      NCAparms <- CalcNCAVec(raw_data)

      # get geo means
      NCAparms <- NCAparms %>%
        group_by(SIM, GROUP) %>%
        summarise(GMeanCmax = exp(mean(log(CMAX), na.rm = TRUE)), GMeanAUC = exp(mean(log(AUC), na.rm = TRUE)), .groups = "keep")

      colnames(NCAparms) <- c("Sample", "Group", "Cmax", "AUC")
    },
    error = function(e) {
      return(-3)
    }
  )

  return(NCAparms)
  #
}

######### Sim
setwd(home.dir)
simparms <- Get_SimNCA(xml_file)
if (!is.list(simparms)) {
  if (simparms == -1) {
    c(crash_value, paste(xml_file, "not found"))
  }
  if (simparms == -2) {
    c(crash_value, paste("Covariance results not available in", xml_file))
  }
  if (simparms == -3) {
    c(crash_value, paste("Simulation failed for", home.dir))
  }
} else {
  simparms$TRT <- "Simulated"
  simparms <- simparms %>%
    group_by(Group)
  Sim.Cmax.means <- simparms %>%
    summarise(Cmax = exp(mean(log(Cmax)))) %>%
    mutate(Group = if_else(Group < 3, "Reference", "Test")) %>%
    group_by(Group) %>%
    summarise(Cmax = mean(Cmax))

  Sim.AUC.means <- simparms %>%
    summarise(AUC = exp(mean(log(AUC)))) %>%
    mutate(Group = if_else(Group < 3, "Reference", "Test")) %>%
    group_by(Group) %>%
    summarise(AUC = mean(AUC))
  Sim.Cmax.cvs <- simparms %>%
    summarise(CmaxCV = sd(log(Cmax))) %>%
    mutate(Group = if_else(Group < 3, "Reference", "Test")) %>%
    group_by(Group) %>%
    summarise(Cmax = mean(CmaxCV))
  Sim.AUC.cvs <- simparms %>%
    summarise(AUCCV = sd(log(AUC))) %>%
    mutate(Group = if_else(Group < 3, "Reference", "Test")) %>%
    group_by(Group) %>%
    summarise(AUC = mean(AUCCV))
  # ratios, group 1 and 2 are reference, 3 and 4 are test

  # mean for reference and test this return bootstrap samples of original data parameters
  bsparms <- Get_TrueNCA(file.path(home.dir, "sim", "org.dat"))
  bs.Ref.means <- rbind(bsparms$Ref1, bsparms$Ref2) %>%
    summarise(AUC = exp(mean(log(AUC))), Cmax = exp(mean(log(Cmax))))
  bs.Test.means <- rbind(bsparms$test1, bsparms$test2) %>%
    summarise(AUC = exp(mean(log(AUC))), Cmax = exp(mean(log(Cmax))))

  # ratios

  Sim.Cmax.Ratio <- mean(Sim.Cmax.means[2, ]$Cmax / Sim.Cmax.means[1, ]$Cmax)
  Sim.AUC.Ratio <- mean(Sim.AUC.means[2, ]$AUC / Sim.AUC.means[1, ]$AUC)

  bs.Cmax.Ratio <- bs.Test.means$Cmax / bs.Ref.means$Cmax
  bs.AUC.Ratio <- bs.Test.means$AUC / bs.Ref.means$AUC

  # penalties, ref, test, Cmax, auc, mean and cv
  Cmax.ratio.penalty <- abs(0.25 * crash_value * (bs.Cmax.Ratio - Sim.Cmax.Ratio) / bs.Cmax.Ratio)

  AUC.ratio.penalty <- abs(0.25 * crash_value * (bs.AUC.Ratio - Sim.AUC.Ratio) / bs.AUC.Ratio)

  sim.Ref.mean.Cmax <- Sim.Cmax.means[1, ]$Cmax
  bs.Ref.mean.Cmax <- bs.Ref.means$Cmax


  sim.Ref.mean.AUC <- Sim.AUC.means[1, ]$AUC
  bs.Ref.mean.AUC <- bs.Ref.means$AUC

  Ref.Cmax.mean.penalty <- abs(0.25 * crash_value * (bs.Ref.mean.Cmax - sim.Ref.mean.Cmax) / bs.Ref.mean.Cmax)

  sim.Test.mean.Cmax <- Sim.Cmax.means[2, ]$Cmax
  bs.Test.mean.Cmax <- bs.Test.means$Cmax
  Test.Cmax.mean.penalty <- abs(0.25 * crash_value * (bs.Test.mean.Cmax - sim.Test.mean.Cmax) / bs.Test.mean.Cmax)

  ## AUC
  sim.Ref.mean.AUC <- Sim.AUC.means[1, ]$AUC
  bs.Ref.mean.AUC <- bs.Ref.means$AUC
  Ref.AUC.mean.penalty <- abs(0.25 * crash_value * (bs.Ref.mean.AUC - sim.Ref.mean.AUC) / bs.Ref.mean.AUC)

  sim.Test.mean.AUC <- Sim.AUC.means[2, ]$AUC
  bs.Test.mean.AUC <- bs.Test.means$AUC
  Test.AUC.mean.penalty <- abs(0.25 * crash_value * (bs.Test.mean.AUC - sim.Test.mean.AUC) / bs.Test.mean.AUC)

  # cvs, only for simulation, reference is a constant
  Cmax.cv.penalty <- sum(Sim.Cmax.cvs$Cmax) * crash_value / 10
  AUC.cv.penalty <- sum(Sim.AUC.cvs$AUC) * crash_value / 10

  penalty <- Cmax.ratio.penalty + AUC.ratio.penalty + Ref.Cmax.mean.penalty + Test.Cmax.mean.penalty + Ref.AUC.mean.penalty +
    Test.AUC.mean.penalty + Cmax.cv.penalty + AUC.cv.penalty
  penalty <- min(crash_value, penalty)
  lstfile <- file(output_file, "a")
  writeLines(paste("\nReference Sim Cmax mean =", as.character(round(sim.Ref.mean.Cmax, 2))), lstfile)
  writeLines(paste("Reference Observed Cmax mean =", as.character(round(bs.Ref.mean.Cmax, 2))), lstfile)
  writeLines(paste("Test Sim Cmax mean =", as.character(round(sim.Test.mean.Cmax, 2))), lstfile)
  writeLines(paste("Test Observed Cmax mean =", as.character(round(bs.Test.mean.Cmax, 2))), lstfile)

  writeLines(paste("Reference Sim AUC mean =", as.character(round(sim.Ref.mean.AUC, 1))), lstfile)
  writeLines(paste("Reference Observed AUC mean =", as.character(round(bs.Ref.mean.AUC, 1))), lstfile)
  writeLines(paste("Test Sim AUC mean =", as.character(round(sim.Test.mean.AUC, 1))), lstfile)
  writeLines(paste("Test Observed AUC mean =", as.character(round(bs.Test.mean.AUC, 1))), lstfile)


  writeLines(paste("Reference Sim Cmax cv =", as.character(round(Sim.Cmax.cvs[1, ]$Cmax, 3))), lstfile)
  writeLines(paste("Text Sim Cmax cv =", as.character(round(Sim.Cmax.cvs[2, ]$Cmax, 3))), lstfile)
  writeLines(paste("Reference Sim AUC cv =", as.character(round(Sim.AUC.cvs[1, ]$AUC, 3))), lstfile)
  writeLines(paste("Text Sim AUC cv =", as.character(round(Sim.AUC.cvs[2, ]$AUC, 3))), lstfile)

  writeLines(paste0("Difference in Reference Cmax = ", as.character(round(
    100 * (bs.Ref.mean.Cmax - sim.Ref.mean.Cmax) / bs.Ref.mean.Cmax,
    2
  )), "%"), lstfile)
  writeLines(paste0("Difference in Reference AUC = ", as.character(round(
    100 * (bs.Ref.mean.AUC - sim.Ref.mean.AUC) / bs.Ref.mean.AUC,
    2
  )), "%"), lstfile)
  writeLines(paste0("Difference in Test Cmax =", as.character(round(
    100 * (bs.Test.mean.Cmax - sim.Test.mean.Cmax) / bs.Test.mean.Cmax,
    2
  )), "%"), lstfile)
  writeLines(paste0("Difference in Test AUC =", as.character(round(
    100 * (bs.Test.mean.AUC - sim.Test.mean.AUC) / bs.Test.mean.AUC,
    2
  )), "%"), lstfile)
  writeLines(paste0("Cmax Test/Reference for Real data = ", as.character(round(bs.Cmax.Ratio, 3))), lstfile)
  writeLines(paste0("AUC Test/Reference for Real data = ", as.character(round(bs.AUC.Ratio, 3))), lstfile)
  writeLines(paste0("Cmax Test/Reference for Simulated data = ", as.character(round(Sim.Cmax.Ratio, 3))), lstfile)
  writeLines(paste0("AUC Test/Reference for Simulated data = ", as.character(round(Sim.AUC.Ratio, 3))), lstfile)

  writeLines(paste("Ratio of AUC Ratio Sim/Real  = ", as.character(round(Sim.AUC.Ratio/bs.AUC.Ratio, 3))), lstfile)
  writeLines(paste("Ratio of Cmax Ratio Sim/Real  = ", as.character(round(Sim.Cmax.Ratio/bs.Cmax.Ratio, 3))), lstfile)
  close(lstfile)
  if (file.exists(file.path(home.dir, "sim"))) {
    unlinking_dir("sim")
  }
  c(penalty, "See above, written by R code") # text written by R, don't return any other text
}
