sampleIterations <- 25000; warmupIterations <- 25000; nChains <- 1
for(i in 1:1){ #hacky way of loading everything
  
  library(rethinking)
  keratinIndex <- 1
  keratin <- c("removeCorrectly", "removeIncorrectly")[keratinIndex]
  dataSizeIndex <- 2
  dataSize <- c("full", "default", "sparse", "defaultfixed")[dataSizeIndex]
  addNewPatients <- T
  EDA <- T
  removeP9 <- F
  if(removeP9){
    dataSize <- paste0(dataSize, "_removeP9")
  }
  
  #finding files
  if(Sys.info()['sysname'] == "Darwin") {
    setwd("/Volumes/macOS/Users/nikolai/heart_transplant/")
  } else {
    setwd("C:\\Users\\Nikolai\\Desktop\\Kate Dissertation\\CSVs")
  }
  
  files <- paste0("spectralCountData/", list.files("spectralCountData"))
  names <- sapply(1:length(files), function (x) strsplit(files[x], split = ".csv")[[1]])
  
  #loading and wrangling files
  columnNames <- c("spectral_count", "patient", "protein", "timepoint", "timepoint_date", "date_of_transplant", "date_of_rejection")
  
  ######################################################################
  # loading in the new data
  d_new <- read.csv("newData.csv", header = T)
  n_timepoints <- 19
  p_new <- data.frame(matrix(NA, nrow = n_timepoints*nrow(d_new), ncol = length(columnNames))); colnames(p_new) <- columnNames
  p_new$spectral_count <- c(sapply(7:25, function(x) d_new[,x]))
  p_new$patient <- c(sapply(7:25, function(x) rep(substr(colnames(d_new)[x],1,2), length(d_new[,x]))))
  p_new$protein <- rep(trimws(d_new$Identified.Proteins..1109.), n_timepoints)
  p_new$timepoint <- c(sapply(7:25, function(x) rep(as.numeric(substr(colnames(d_new)[x],5,5))+1, length(d_new[,x]))))
  timepoint_dates <- read.csv("timepoint_dates.csv", header = T, stringsAsFactors = F)
  p_new$timepoint_date <- c(sapply(7:25, function(x) rep((timepoint_dates$Timepoint_Date[
      substr(colnames(d_new)[x],2,2) == timepoint_dates$Patient & 
      as.character(as.numeric(substr(colnames(d_new)[x],5,5)) + 1) == timepoint_dates$TimePoint
      ]), length(d_new[,x]))))
  p_new$timepoint_date <- as.Date(p_new$timepoint_date, format = "%m/%d/%y")
  p_new$date_of_transplant <- c(sapply(7:25, function(x) rep((timepoint_dates$Transplant_Date[
      substr(colnames(d_new)[x],2,2) == timepoint_dates$Patient & 
      as.character(as.numeric(substr(colnames(d_new)[x],5,5)) + 1) == timepoint_dates$TimePoint
      ]), length(d_new[,x]))))
  p_new$date_of_transplant <- as.Date(p_new$date_of_transplant, format = "%m/%d/%y")
  
  
    
  ######################################################################
  p1_T1vT2 <- read.csv(files[1]); p1_T1vT2 <- p1_T1vT2[1:609,]
  p1 <- data.frame(c(p1_T1vT2$s1_T1, p1_T1vT2$s1_T2), #spectral count
                   rep(p1_T1vT2$Person, 2), #patient
                   rep(p1_T1vT2$Description, 2), #protein
                   c(rep(1, dim(p1_T1vT2)[1]), rep(2, dim(p1_T1vT2)[1])), #discrete timepoint
                   c(as.Date((p1_T1vT2$T1), format = "%m/%d/%y"), as.Date((p1_T1vT2$T2), format = "%m/%d/%y")), #timepoint date
                   as.Date(rep(p1_T1vT2$DOT,2), format = "%m/%d/%Y"), #date of transplant
                   as.Date(rep(p1_T1vT2$DOR,2), format = "%m/%d/%Y") #date of rejection
  )
  colnames(p1) <- columnNames
  
  ######################################################################
  p3_T1vT2 <- read.csv(files[2])
  p3 <- data.frame(c(p3_T1vT2$s3_T1, p3_T1vT2$s3_T2), #spectral count
                   rep(p3_T1vT2$Person, 2), #patient
                   rep(p3_T1vT2$Description, 2), #protein
                   c(rep(1, dim(p3_T1vT2)[1]), rep(2, dim(p3_T1vT2)[1])), #discrete timepoint
                   c(as.Date((p3_T1vT2$T1), format = "%m/%d/%y"), as.Date((p3_T1vT2$T2), format = "%m/%d/%y")), #timepoint date
                   as.Date(rep(p3_T1vT2$DOT,2), format = "%m/%d/%Y"), #date of transplant
                   as.Date(rep(p3_T1vT2$DOR,2), format = "%m/%d/%Y") #date of rejection
  )
  colnames(p3) <- columnNames
  
  ######################################################################
  p3_T1vT3 <- read.csv(files[3])
  p3_t3 <- data.frame(p3_T1vT3$s3_T3, #spectral count
                      rep(p3_T1vT3$Person, 1), #patient
                      rep(p3_T1vT3$Description, 1), #protein
                      rep(3, dim(p3_T1vT3)[1]), #discrete timepoint
                      as.Date((p3_T1vT3$T3), format = "%m/%d/%y"), #timepoint date
                      as.Date(rep(p3_T1vT3$DOT,1), format = "%m/%d/%Y"), #date of transplant
                      as.Date(rep(p3_T1vT3$DOR,1), format = "%m/%d/%Y") #date of rejection
  )
  colnames(p3_t3) <- columnNames
  p3 <- rbind(p3, p3_t3)
  
  in13not12 <- setdiff(p3_T1vT3$Description, p3_T1vT2$Description)
  p3r_T1vT3 <- p3_T1vT3[sapply(1:length(in13not12), function(x) which(p3_T1vT3$Description == in13not12[x])),]
  p3_t1 <- data.frame(p3r_T1vT3$s3_T1, #spectral count
                      rep(p3r_T1vT3$Person, 1), #patient
                      rep(p3r_T1vT3$Description, 1), #protein
                      rep(1, dim(p3r_T1vT3)[1]), #discrete timepoint
                      as.Date((p3r_T1vT3$T3), format = "%m/%d/%y"), #timepoint date
                      as.Date(rep(p3r_T1vT3$DOT,1), format = "%m/%d/%Y"), #date of transplant
                      as.Date(rep(p3r_T1vT3$DOR,1), format = "%m/%d/%Y") #date of rejection
  )
  colnames(p3_t1) <- columnNames
  p3 <- rbind(p3, p3_t1)
  
  ######################################################################
  p5_T1vT2 <- read.csv(files[4]); p5_T1vT2 <- p5_T1vT2[1:260,]
  p5 <- data.frame(c(p5_T1vT2$s5_T1, p5_T1vT2$s5_T2), #spectral count
                   rep(p5_T1vT2$Person, 2), #patient
                   rep(p5_T1vT2$Description, 2), #protein
                   c(rep(1, dim(p5_T1vT2)[1]), rep(2, dim(p5_T1vT2)[1])), #discrete timepoint
                   c(as.Date((p5_T1vT2$T1), format = "%m/%d/%y"), as.Date((p5_T1vT2$T2), format = "%m/%d/%y")), #timepoint date
                   as.Date(rep(p5_T1vT2$DOT,2), format = "%m/%d/%Y"), #date of transplant
                   as.Date(rep(p5_T1vT2$DOR,2), format = "%m/%d/%Y") #date of rejection
  )
  colnames(p5) <- columnNames
  
  p5_T1vT3 <- read.csv(files[5]); p5_T1vT3 <- p5_T1vT3[1:312,]
  p5_t3 <- data.frame(p5_T1vT3$s5_T3, #spectral count
                      rep(p5_T1vT3$Person, 1), #patient
                      rep(p5_T1vT3$Description, 1), #protein
                      rep(3, dim(p5_T1vT3)[1]), #discrete timepoint
                      as.Date((p5_T1vT3$T3), format = "%m/%d/%y"), #timepoint date
                      as.Date(rep(p5_T1vT3$DOT,1), format = "%m/%d/%Y"), #date of transplant
                      as.Date(rep(p5_T1vT3$DOR,1), format = "%m/%d/%Y") #date of rejection
  )
  colnames(p5_t3) <- columnNames
  p5 <- rbind(p5, p5_t3)
  
  in13not12 <- setdiff(p5_T1vT3$Description, p5_T1vT2$Description)
  p5r_T1vT3 <- p5_T1vT3[sapply(1:length(in13not12), function(x) which(p5_T1vT3$Description == in13not12[x])),]
  p5_t1 <- data.frame(p5r_T1vT3$s5_T1, #spectral count
                      rep(p5r_T1vT3$Person, 1), #patient
                      rep(p5r_T1vT3$Description, 1), #protein
                      rep(1, dim(p5r_T1vT3)[1]), #discrete timepoint
                      as.Date((p5r_T1vT3$T3), format = "%m/%d/%y"), #timepoint date
                      as.Date(rep(p5r_T1vT3$DOT,1), format = "%m/%d/%Y"), #date of transplant
                      as.Date(rep(p5r_T1vT3$DOR,1), format = "%m/%d/%Y") #date of rejection
  )
  colnames(p5_t1) <- columnNames
  p5 <- rbind(p5, p5_t1)
  
  ######################################################################
  p6_T1vT2 <- read.csv(files[6]); p6_T1vT2 <- p6_T1vT2[1:185,]
  p6 <- data.frame(c(p6_T1vT2$s6_T1, p6_T1vT2$s6_T2), #spectral count
                   rep(p6_T1vT2$Person, 2), #patient
                   rep(p6_T1vT2$Description, 2), #protein
                   c(rep(1, dim(p6_T1vT2)[1]), rep(2, dim(p6_T1vT2)[1])), #discrete timepoint
                   c(as.Date((p6_T1vT2$T1), format = "%m/%d/%y"), as.Date((p6_T1vT2$T2), format = "%m/%d/%y")), #timepoint date
                   as.Date(rep(p6_T1vT2$DOT,2), format = "%m/%d/%Y"), #date of transplant
                   as.Date(rep(p6_T1vT2$DOR,2), format = "%m/%d/%Y") #date of rejection
  )
  colnames(p6) <- columnNames
  
  p6_T1vT3 <- read.csv(files[7]); p6_T1vT3 <- p6_T1vT3[1:220,]
  p6_t3 <- data.frame(p6_T1vT3$s6_T3, #spectral count
                      rep(p6_T1vT3$Person, 1), #patient
                      rep(p6_T1vT3$Description, 1), #protein
                      rep(3, dim(p6_T1vT3)[1]), #discrete timepoint
                      as.Date((p6_T1vT3$T3), format = "%m/%d/%y"), #timepoint date
                      as.Date(rep(p6_T1vT3$DOT,1), format = "%m/%d/%Y"), #date of transplant
                      as.Date(rep(p6_T1vT3$DOR,1), format = "%m/%d/%Y") #date of rejection
  )
  colnames(p6_t3) <- columnNames
  p6 <- rbind(p6, p6_t3)
  
  in13not12 <- setdiff(p6_T1vT3$Description, p6_T1vT2$Description)
  p6r_T1vT3 <- p6_T1vT3[sapply(1:length(in13not12), function(x) which(p6_T1vT3$Description == in13not12[x])),]
  p6_t1 <- data.frame(p6r_T1vT3$s6_T1, #spectral count
                      rep(p6r_T1vT3$Person, 1), #patient
                      rep(p6r_T1vT3$Description, 1), #protein
                      rep(1, dim(p6r_T1vT3)[1]), #discrete timepoint
                      as.Date((p6r_T1vT3$T3), format = "%m/%d/%y"), #timepoint date
                      as.Date(rep(p6r_T1vT3$DOT,1), format = "%m/%d/%Y"), #date of transplant
                      as.Date(rep(p6r_T1vT3$DOR,1), format = "%m/%d/%Y") #date of rejection
  )
  colnames(p6_t1) <- columnNames
  p6 <- rbind(p6, p6_t1)
  
  ######################################################################
  p7_T1vT2 <- read.csv(files[8]); p7_T1vT2 <- p7_T1vT2[1:438,]
  p7 <- data.frame(c(p7_T1vT2$s7_T1, p7_T1vT2$s7_T2), #spectral count
                   rep(p7_T1vT2$Person, 2), #patient
                   rep(p7_T1vT2$Description, 2), #protein
                   c(rep(1, dim(p7_T1vT2)[1]), rep(2, dim(p7_T1vT2)[1])), #discrete timepoint
                   c(as.Date((p7_T1vT2$T1), format = "%m/%d/%y"), as.Date((p7_T1vT2$T2), format = "%m/%d/%y")), #timepoint date
                   as.Date(rep(p7_T1vT2$DOT,2), format = "%m/%d/%Y"), #date of transplant
                   as.Date(rep(p7_T1vT2$DOR,2), format = "%m/%d/%Y") #date of rejection
  )
  colnames(p7) <- columnNames
  
  p7_T1vT3 <- read.csv(files[9]); p7_T1vT3 <- p7_T1vT3[1:505,]
  p7_t3 <- data.frame(p7_T1vT3$s7_T3, #spectral count
                      rep(p7_T1vT3$Person, 1), #patient
                      rep(p7_T1vT3$Description, 1), #protein
                      rep(3, dim(p7_T1vT3)[1]), #discrete timepoint
                      as.Date((p7_T1vT3$T3), format = "%m/%d/%y"), #timepoint date
                      as.Date(rep(p7_T1vT3$DOT,1), format = "%m/%d/%Y"), #date of transplant
                      as.Date(rep(p7_T1vT3$DOR,1), format = "%m/%d/%Y") #date of rejection
  )
  colnames(p7_t3) <- columnNames
  p7 <- rbind(p7, p7_t3)
  
  in13not12 <- setdiff(p7_T1vT3$Description, p7_T1vT2$Description)
  p7r_T1vT3 <- p7_T1vT3[sapply(1:length(in13not12), function(x) which(p7_T1vT3$Description == in13not12[x])),]
  p7_t1 <- data.frame(p7r_T1vT3$s7_T1, #spectral count
                      rep(p7r_T1vT3$Person, 1), #patient
                      rep(p7r_T1vT3$Description, 1), #protein
                      rep(1, dim(p7r_T1vT3)[1]), #discrete timepoint
                      as.Date((p7r_T1vT3$T3), format = "%m/%d/%y"), #timepoint date
                      as.Date(rep(p7r_T1vT3$DOT,1), format = "%m/%d/%Y"), #date of transplant
                      as.Date(rep(p7r_T1vT3$DOR,1), format = "%m/%d/%Y") #date of rejection
  )
  colnames(p7_t1) <- columnNames
  p7 <- rbind(p7, p7_t1)
  
  ######################################################################
  
  p9_T1vT2 <- read.csv(files[10]); p9_T1vT2 <- p9_T1vT2[1:1672,]
  p9 <- data.frame(c(p9_T1vT2$s9_T1, p9_T1vT2$s9_T2), #spectral count
                   rep(p9_T1vT2$Person, 2), #patient
                   rep(p9_T1vT2$Description, 2), #protein
                   c(rep(1, dim(p9_T1vT2)[1]), rep(2, dim(p9_T1vT2)[1])), #discrete timepoint
                   c(as.Date((p9_T1vT2$T1), format = "%m/%d/%y"), as.Date((p9_T1vT2$T2), format = "%m/%d/%y")), #timepoint date
                   as.Date(rep(p9_T1vT2$DOT,2), format = "%m/%d/%Y"), #date of transplant
                   as.Date(rep(p9_T1vT2$DOR,2), format = "%m/%d/%Y") #date of rejection
  )
  colnames(p9) <- columnNames
  
  p9_T1vT3 <- read.csv(files[11]); p9_T1vT3 <- p9_T1vT3[1:1792,]
  p9_t3 <- data.frame(p9_T1vT3$s9_T3, #spectral count
                      rep(p9_T1vT3$Person, 1), #patient
                      rep(p9_T1vT3$Description, 1), #protein
                      rep(3, dim(p9_T1vT3)[1]), #discrete timepoint
                      as.Date((p9_T1vT3$T3), format = "%m/%d/%y"), #timepoint date
                      as.Date(rep(p9_T1vT3$DOT,1), format = "%m/%d/%Y"), #date of transplant
                      as.Date(rep(p9_T1vT3$DOR,1), format = "%m/%d/%Y") #date of rejection
  )
  colnames(p9_t3) <- columnNames
  p9 <- rbind(p9, p9_t3)
  
  in13not12 <- setdiff(p9_T1vT3$Description, p9_T1vT2$Description)
  p9r_T1vT3 <- p9_T1vT3[sapply(1:length(in13not12), function(x) which(p9_T1vT3$Description == in13not12[x])),]
  p9_t1 <- data.frame(p9r_T1vT3$s9_T1, #spectral count
                      rep(p9r_T1vT3$Person, 1), #patient
                      rep(p9r_T1vT3$Description, 1), #protein
                      rep(1, dim(p9r_T1vT3)[1]), #discrete timepoint
                      as.Date((p9r_T1vT3$T3), format = "%m/%d/%y"), #timepoint date
                      as.Date(rep(p9r_T1vT3$DOT,1), format = "%m/%d/%Y"), #date of transplant
                      as.Date(rep(p9r_T1vT3$DOR,1), format = "%m/%d/%Y") #date of rejection
  )
  colnames(p9_t1) <- columnNames
  p9 <- rbind(p9, p9_t1)
  
  ######################################################################
  ######################### DATA IS READ IN NOW ########################
  ######################################################################
  ######################## AND ALSO MANIPULATED ########################
  ######################################################################
  
  if(!removeP9){
    d <- d.orig <- rbind(p1, p3, p5, p6, p7, p9)
  } else if(removeP9){
    d <- d.orig <- rbind(p1, p3, p5, p6, p7)
  }
  
  if(addNewPatients){
    d <- d.orig <- rbind(d, p_new)
  }
  
  #check overlap between old and new proteins
  old_proteins <- unique(rbind(p1, p3, p5, p6, p7, p9)$protein)
  new_proteins <- unique(p_new$protein)
  length(old_proteins) + length(new_proteins)
  length(unique(c(old_proteins, new_proteins)))
  
  if(dataSize == "full" || dataSize == "full_removeP9"){
    #populate the data with zero spectral counts where appropriate
    patients <- as.character(unique(d$patient))
    proteins <- as.character(unique(d$protein))
    
    for(i in 1:length(patients)){
      print(paste0("patient ",i))
      for(j in 1:(if(any(i==c(1,7,13))){2}else{3})){ #timepoints
        patient <- patients[i]
        timepoint <- j
        subData <- d[d$patient == patient & d$timepoint == timepoint,]
        proteinsHere <- subData$protein
        if(j == 1){print(paste0("num proteins in this patient = ", length(proteinsHere)))}
        proteinsNotHere <- setdiff(proteins, proteinsHere)
        numProteinsNotHere <- length(proteinsNotHere)
        print(paste0("num proteins added to timepoint ", timepoint, ": ", numProteinsNotHere))

        zeroesToAdd <- data.frame(rep(0, numProteinsNotHere), #spectral count
                                  rep(patient, numProteinsNotHere), #patient
                                  proteinsNotHere, #protein
                                  rep(timepoint, numProteinsNotHere), #discrete timepoint
                                  rep(subData$timepoint_date[1], numProteinsNotHere), #timepoint date
                                  rep(subData$date_of_transplant[1], numProteinsNotHere), #date of transplant
                                  rep(subData$date_of_rejection[1], numProteinsNotHere) #date of rejection
        )
        colnames(zeroesToAdd) <- columnNames
        d <- rbind(d, zeroesToAdd)
      }
    }
  }
  
  if(dataSize == "defaultfixed"){
    #populate the data with zero spectral counts where appropriate
    patients <- as.character(unique(d$patient))
    
    for(i in 1:length(patients)){
      print(paste0("patient ",i))
      for(j in 1:(if(i==1){2}else{3})){ #timepoints
        patient <- patients[i]
        timepoint <- j
        subData <- d[d$patient == patient & d$timepoint == timepoint,]
        proteins <- as.character(unique(d[d$patient == patient,]$protein))
        proteinsHere <- subData$protein
        proteinsNotHere <- setdiff(proteins, proteinsHere)
        numProteinsNotHere <- length(proteinsNotHere)
        print(numProteinsNotHere)
        
        zeroesToAdd <- data.frame(rep(0, numProteinsNotHere), #spectral count
                                  rep(patient, numProteinsNotHere), #patient
                                  proteinsNotHere, #protein
                                  rep(timepoint, numProteinsNotHere), #discrete timepoint
                                  rep(subData$timepoint_date[1], numProteinsNotHere), #timepoint date
                                  rep(subData$date_of_transplant[1], numProteinsNotHere), #date of transplant
                                  rep(subData$date_of_rejection[1], numProteinsNotHere) #date of rejection
        )
        colnames(zeroesToAdd) <- columnNames
        d <- rbind(d, zeroesToAdd)
      }
    }
  }
  
  #filter out keratin and "null" protein
  if(keratin == "removeCorrectly"){
    d <- d[!(grepl(pattern = "keratin", as.character(d$protein)) | grepl(pattern = "Keratin", as.character(d$protein)) | as.character(d$protein) == ""),]
  } else if (keratin == "removeIncorrectly"){
    d <- d[!(grepl(pattern = "keratin", as.character(d$protein)) | as.character(d$protein) == ""),]
  }
  d <- d[!(grepl(pattern = "CONTAMINANT", as.character(d$protein)) | grepl(pattern = "Reversed", as.character(d$protein)) | grepl(pattern = "Immunoglobulin", as.character(d$protein))),]
  
  
  if(dataSize == "sparse"){
    proteins <- as.character(unique(d$protein))
    numPatientsPerProtein <- sapply(1:length(proteins), function(x) length(unique(d[d$protein == proteins[x] & d$spectral_count != 0,]$patient)))
    proteinsShared <- proteins[numPatientsPerProtein >= 3]
    d <- d[d$protein %in% proteinsShared,]
  }
  
  #create other variables of interest
  d$rejection <- ifelse(is.na(d$date_of_rejection), 0, 1)
  d$timepoint1 <- as.numeric(d$timepoint == 1)
  d$timepoint2 <- as.numeric(d$timepoint == 2)
  d$timepoint3 <- as.numeric(d$timepoint == 3)
  d$daysRelativeToRejection <- as.numeric(d$timepoint_date - d$date_of_rejection)
  d$afterRejection <- as.numeric(d$daysRelativeToRejection > 0)
  d$daysRelativeToTransplant <- as.numeric(d$timepoint_date - d$date_of_transplant)
  d$adjDaysRelativeToTransplant <- d$daysRelativeToTransplant + mean(unique(as.numeric(d$date_of_transplant - d$date_of_rejection)))
  d$afterTransplant <- as.numeric(d$daysRelativeToTransplant > 0)
  d$daysBetwTrans_Rej <- as.numeric(d$date_of_rejection - d$date_of_transplant)
  d$transplantStdDaysRelativeToRejection <- d$daysRelativeToRejection / d$daysBetwTrans_Rej
  d$phase <- 1; d$phase <- d$phase + (d$daysRelativeToTransplant > 0) + (d$daysRelativeToRejection > 0); 
  table(d$phase)
  # d$pivot <- NA
  
  
  #coerce the protein id
  d$protein_id <- coerce_index(d$protein)
  d$patient_id <- coerce_index(d$patient)
  
  #EDA
  if(EDA == T){
    library(RColorBrewer)
    
    max(d$spectral_count)
    hist(d$daysRelativeToTransplant)
    plot(c(1,1), xlim = c(-800,3000), ylim = c(3,12), xlab = "days relative to transplant", ylab = "mean nonzero spectral count / timepoint", type = "n")
    par(xpd=TRUE); points(0, 12.55, pch = 25, bg = 2, col = 2); points(c(-18,18), c(12.62,12.62), pch = 16, col = 2, cex = 0.85); par(xpd=F)
    abline(v = 0, col = 2, lty = 2)
    cols <- brewer.pal(6, "Dark2")
    cols <- rep("red", 6)
    rej_i <- 1
    for(pat in 1:13){
        subd <- d[d$patient_id == unique(d$patient_id)[pat],]
        subd <- lapply(1:length(unique(subd$timepoint)), function(x) subd[subd$timepoint == x,]) 
        meanSC <- sapply(1:length(subd), function(x) mean(subd[[x]]$spectral_count[subd[[x]]$spectral_count != 0]))
    #    meanSC <- meanSC - meanSC[1] + 8 #toggle to fix mean at 8
        numprots <- sapply(1:length(subd), function(x) length(subd[[x]]$spectral_count[subd[[x]]$spectral_count != 0]))
        daysFromTransplant <- sapply(1:length(subd), function(x) mean(subd[[x]]$daysRelativeToTransplant))
        rejector <- subd[[1]]$rejection[1] 
        lines(x = daysFromTransplant, y = meanSC, col = ifelse(rejector, "red", "black"))
        # points(x = daysFromTransplant, y = meanSC, pch = 16, cex = ifelse(rejector, 1, 0.4), col = ifelse(rejector, cols[rej_i], "black"))
        points(x = daysFromTransplant, y = meanSC, pch = 16, cex = (numprots) / 600, col = ifelse(rejector, cols[rej_i], "black"))
        if(rejector){
          b_length <- 0.2
          h_pos <- subd[[1]]$daysBetwTrans_Rej[1]
          diffs <- daysFromTransplant - h_pos
          tps_rj <- c(max(which(diffs < 0)), min(which(diffs > 0)))
          cent_v <- meanSC[tps_rj[1]] + (meanSC[tps_rj[2]] - meanSC[tps_rj[1]] ) / (daysFromTransplant[tps_rj[2]] - daysFromTransplant[tps_rj[1]]) * -diffs[tps_rj[1]] 
          segments(x0 = h_pos, x1 = h_pos, y0 = cent_v - b_length / 2, y1 = cent_v + b_length / 2, lwd = 1.5, col = "red")
          if(subd[[1]]$patient[1] == "P7"){
            segments(x0 = h_pos, x1 = h_pos, y0 = meanSC[1] - b_length / 2, y1 = meanSC[1] + b_length / 2, lwd = 0.5, col = "red")
          }
            
          #abline(v = subd[[1]]$daysBetwTrans_Rej[1], col = cols[rej_i])
          # print(c(pat, cent_v))
          rej_i <- rej_i + 1;
        }
        # text(round(mean(numprots)), x = daysFromTransplant[length(daysFromTransplant)] + 90, y = meanSC[length(meanSC)])
        # text(subd[[1]]$patient[1], x = daysFromTransplant[1], y = meanSC[1] + 0.2)
        print(paste0(rejector, ": ", numprots))
    }
    legend("bottomright", legend = c("rejector", "non-rejector / \"control\""), col = c(2,1), lty = 1)
    legend("bottomleft", legend = c(100, 200, 400, 800), col = 2, pch = 16, pt.cex = c(100, 200, 400, 800) / 600, title = "# proteins / tp")
    text("l", col = 2, x = 1800, y = 12.1, cex = 1.3); text("represent rejection dates", col = 1, x = 2450, y = 12.1)
    title("Vertical Lines Represent Rejection Dates, Non-Red Colors Index Rejecting Patients,\n Numbers Represent Mean # of Unique Proteins Across All Available Timepoints")
    }
    
    numPatientsPerProtein <- sapply(1:length(proteins), function(x) length(unique(d[d$protein == proteins[x] & d$spectral_count != 0,]$patient)))
    proteinsShared <- proteins[numPatientsPerProtein > 5]
    cols <- brewer.pal(6, "Dark2")
    png(filename = paste0("allPatients.png"), width = 1000, height = 750)
    par(mfrow = c(2,3))
    for(i in 1:length(patients)){
      print(i)
      patient <- patients[i]
      # png(filename = paste0(patient, "SCtT.png"), width = 1500, height = 2000)
      plot(1,1,type = "l", col = "white", xlim = c(-1500, 1500), ylim = c(0,175), 
           xlab = "Days Relative to Acute Rejection Event", ylab = "Spectral Count")
      title(paste0("Patient ", patient))
      for(k in 1:length(proteins)){
        protein <- proteins[k]
        subData <- d[d$patient == patient & d$protein == protein,]
        subData <- subData[order(subData$timepoint),]
        times <- subData$daysRelativeToRejection
        counts <- subData$spectral_count
        lines(times, counts, col = cols[i])
        abline(v = times)
      }
      for(l in 1:length(proteinsShared)){
        protein <- proteinsShared[l]
        subData <- d[d$patient == patient & d$protein == protein,]
        subData <- subData[order(subData$timepoint),]
        times <- subData$daysRelativeToRejection
        counts <- subData$spectral_count
        lines(times, counts, col = 1)
        abline(v = times, lwd = 2)
      }
      abline(v = as.numeric(subData$date_of_transplant[1] - subData$date_of_rejection[1]), lwd = 3, col = 2)
      # dev.off()
    }
    dev.off()
    
    #look at protein specific graphs
    # for(i in 1:length(proteinsShared)){
    #   print(i)
    #   protein <- proteinsShared[i]
    #   png(filename = paste0("protein", i, ".png"), width = 1000, height = 1000)
    #   plot(1,1,type = "l", col = "white", xlim = c(-1500, 1500), ylim = c(0,200), 
    #        xlab = "Days Relative to Acute Rejection Event", ylab = "Spectral Count")
    #   title(paste0(protein))
    #   for(j in 1:length(patients)){
    #     patient <- patients[j]
    #     subData <- d[d$patient == patient & d$protein == protein,]
    #     subData <- subData[order(subData$timepoint),]
    #     times <- subData$daysRelativeToRejection
    #     counts <- subData$spectral_count
    #     lines(times, counts, col = cols[j])
    #     abline(v = times, lwd = 1, col = cols[j])
    #     abline(v = as.numeric(subData$date_of_transplant[1] - subData$date_of_rejection[1]), lwd = 3, col = cols[j])
    #   }
    #   dev.off()
    # }
  }
}
