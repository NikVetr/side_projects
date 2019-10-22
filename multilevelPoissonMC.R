sampleIterations <- 25000; warmupIterations <- 25000; nChains <- 1
for(i in 1:1){ #hacky way of loading everything
  
  library(rethinking)
  keratinIndex <- 1
  keratin <- c("removeCorrectly", "removeIncorrectly")[keratinIndex]
  dataSizeIndex <- 1
  dataSize <- c("full", "default", "sparse", "defaultfixed")[dataSizeIndex]
  EDA <- F
  
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
  
  d <- d.orig <- rbind(p1, p3, p5, p6, p7, p9)
  
  if(dataSize == "full"){
    #populate the data with zero spectral counts where appropriate
    patients <- as.character(unique(d$patient))
    proteins <- as.character(unique(d$protein))
    
    for(i in 1:length(patients)){
      print(paste0("patient ",i))
      for(j in 1:(if(i==1){2}else{3})){ #timepoints
        patient <- patients[i]
        timepoint <- j
        subData <- d[d$patient == patient & d$timepoint == timepoint,]
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
  
  if(dataSize == "sparse"){
    proteins <- as.character(unique(d$protein))
    numPatientsPerProtein <- sapply(1:length(proteins), function(x) length(unique(d[d$protein == proteins[x] & d$spectral_count != 0,]$patient)))
    proteinsShared <- proteins[numPatientsPerProtein >= 3]
    d <- d[d$protein %in% proteinsShared,]
  }
  
  #create other variables of interest
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
    numPatientsPerProtein <- sapply(1:length(proteins), function(x) length(unique(d[d$protein == proteins[x] & d$spectral_count != 0,]$patient)))
    proteinsShared <- proteins[numPatientsPerProtein > 5]
    library(RColorBrewer)
    cols <- brewer.pal(6, "Dark2")
    png(filename = paste0("allPatients.png"), width = 1500, height = 2000)
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
  
  #get appropriate object types
  d$spectral_count <- as.integer(d$spectral_count)
  d$daysRelativeToRejection <- as.integer(d$daysRelativeToRejection)
  d$afterRejection <- as.integer(d$afterRejection)
  d$daysRelativeToTransplant <- as.integer(d$daysRelativeToTransplant)
  d$afterTransplant <- as.integer(d$afterTransplant)
  d$daysBetwTrans_Rej <- as.integer(d$daysBetwTrans_Rej)
  
}

######################################################################
##################### NOW LET"S DO SOME MODELING #####################
######################################################################

#separate slopes for before and after the rejection event, seperate "intercepts" at the rejection event, intercept with zero slope before the transplant event,
#express slope in change per unit day between transplant and rejection event, express post-rejection intercept as displacement from pre-rejection intercept
#this is model 27 with tighter priors
nChains <- 1
m40 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0 + bR0*daysRelativeToTransplant*(1-afterRejection)*afterTransplant + bR1*daysRelativeToRejection*afterRejection, 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id] + a1*afterRejection,
    A0b ~ dnorm(0,1),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,0.25),
    
    a1 <- bR0*daysBetwTrans_Rej + A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,1),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,0.25),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,0.01),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,0.001),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,0.01),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,0.001)
    
  ) ,
  data= d, sample = F,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains, 
  start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, A0b=0, sigA0j=0.5, sigA0i=0.5, A1b=0, sigA1j=0.5, sigA1i=0.5)
)

write(m40$model, file = "~/cmdstan-2.18.0/transplant/m40.stan")
out_data <- list( spectral_count=d$spectral_count , 
                  daysRelativeToRejection=d$daysRelativeToRejection, 
                  afterRejection=d$afterRejection,
                  daysRelativeToTransplant=d$daysRelativeToTransplant,
                  afterTransplant=d$afterTransplant,
                  protein_id=d$protein_id,
                  patient_id=d$patient_id,
                  daysBetwTrans_Rej=d$daysBetwTrans_Rej
)
out_data$N <- nrow(d)
out_data$N_protein_id <- length(unique(d$protein_id))
out_data$N_patient_id <- length(unique(d$patient_id))
stan_rdump(ls(out_data), "~/cmdstan-2.18.0/transplant/transplant_input.R", envir = list2env(out_data))

nChains <- 1
d$daysRelativeToTransplantPos <- d$daysRelativeToTransplant 
d$daysRelativeToTransplantPos[d$daysRelativeToTransplantPos < 0] <- 1e4
d$daysRelativeToRejectionPos <- d$daysRelativeToRejection 
d$daysRelativeToRejectionPos[d$daysRelativeToRejectionPos < 0] <- 1e4
m41 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0 + bR0*sqrt(daysRelativeToTransplantPos)*(1-afterRejection)*afterTransplant + bR1*sqrt(daysRelativeToRejectionPos)*afterRejection, 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id] + a1*afterRejection,
    A0b ~ dnorm(0,1),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,0.25),
    
    a1 <- bR0*sqrt(daysBetwTrans_Rej) + A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,1),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,0.25),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,0.1),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,0.01),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,0.1),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,0.01)
    
  ) ,
  data= d, sample = F,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains, 
  start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, A0b=0, sigA0j=0.5, sigA0i=0.5, A1b=0, sigA1j=0.1, sigA1i=0.1)
)

write(m41$model, file = "~/cmdstan-2.18.0/transplantM41/m41.stan")
out_data <- list( spectral_count=d$spectral_count , 
                  daysRelativeToRejectionPos=d$daysRelativeToRejectionPos, 
                  afterRejection=d$afterRejection,
                  daysRelativeToTransplantPos=d$daysRelativeToTransplantPos,
                  afterTransplant=d$afterTransplant,
                  protein_id=d$protein_id,
                  patient_id=d$patient_id,
                  daysBetwTrans_Rej=d$daysBetwTrans_Rej
)
out_data$N <- nrow(d)
out_data$N_protein_id <- length(unique(d$protein_id))
out_data$N_patient_id <- length(unique(d$patient_id))
stan_rdump(ls(out_data), "~/cmdstan-2.18.0/transplantM41/transplant_input.R", envir = list2env(out_data))


#separate slopes for before and after the rejection event, seperate "intercepts" at the rejection event, intercept with zero slope before the transplant event,
#express slope in change per unit day between transplant and rejection event, express post-rejection intercept as displacement from pre-rejection intercept
#allow for exponential effect on relative change per day, learn exponent from data
nChains <- 1
d$daysRelativeToTransplantPos <- d$daysRelativeToTransplant 
d$daysRelativeToTransplantPos[d$daysRelativeToTransplantPos < 0] <- 1e4
d$daysRelativeToRejectionPos <- d$daysRelativeToRejection 
d$daysRelativeToRejectionPos[d$daysRelativeToRejectionPos < 0] <- 1e4
sampleIterations <- 100000; warmupIterations <- 50000; nChains <- 1
m35 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0 + bR0*(exp(1/p2tuner*log(daysRelativeToTransplantPos)))*(1-afterRejection)*afterTransplant + bR1*(exp(1/p3tuner*log(daysRelativeToRejectionPos)))*afterRejection,
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id] + a1*afterRejection,
    A0b ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,1),
    
    a1 <- bR0*(exp(1/p2tuner*log(daysBetwTrans_Rej))) + A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,1),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,1),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,1),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,1),
    
    p2tuner <- exp(tunerExp2),
    p3tuner <- exp(tunerExp3),
    c(tunerExp2, tunerExp3) ~ dnorm(0, 1.1513)
    
  ) ,
  data= d, sample = F,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains, 
  start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1, tunerExp2= 0, tunerExp3=0)
)

if(!dir.exists("~/cmdstan-2.18.0/m35")){dir.create("~/cmdstan-2.18.0/m35")}
write(m35$model, file = "~/cmdstan-2.18.0/m35/m35.stan")
out_data <- list( spectral_count=d$spectral_count , 
                  daysRelativeToRejectionPos=d$daysRelativeToRejectionPos, 
                  afterRejection=d$afterRejection,
                  daysRelativeToTransplantPos=d$daysRelativeToTransplantPos,
                  afterTransplant=d$afterTransplant,
                  protein_id=d$protein_id,
                  patient_id=d$patient_id,
                  daysBetwTrans_Rej=d$daysBetwTrans_Rej
)
out_data$N <- nrow(d)
out_data$N_protein_id <- length(unique(d$protein_id))
out_data$N_patient_id <- length(unique(d$patient_id))
stan_rdump(ls(out_data), "~/cmdstan-2.18.0/m35/transplant_input.R", envir = list2env(out_data))

