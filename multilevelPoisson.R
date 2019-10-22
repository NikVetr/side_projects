sampleIterations <- 25000; warmupIterations <- 25000; nChains <- 1
for(i in 1:1){ #hacky way of loading everything
  
  library(rethinking)
  keratinIndex <- 1
  keratin <- c("removeCorrectly", "removeIncorrectly")[keratinIndex]
  dataSizeIndex <- 1
  dataSize <- c("full", "default", "sparse", "defaultfixed")[dataSizeIndex]
  EDA <- F
  removeP9 <- T
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
  if(dataSize == "full" || dataSize == "full_removeP9"){
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
}

######################################################################
##################### NOW LET"S DO SOME MODELING #####################
######################################################################

nChains <- 1

#model has just intercept and one coefficient corresponding to an indicator variable signifying occurence of discrete rejection event
m1 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a + bR*afterRejection, 
    a ~ dnorm(0,100),
    bR ~ dnorm(0,1)
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = 4,
  start = list(a=0, bR=0)
)
save(file = paste0(dataSize, "/", "m1", dataSize), m1)
precis(m1)
#plot(m1)
#pairs(m1)

#model has a continuous parameter indicating the number of days relative to the transplant
m2 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a + bDR*daysRelativeToRejection, 
    a ~ dnorm(0,100),
    bDR ~ dnorm(0,1)
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(a=0, bDR=0)
)
save(file = paste0(dataSize, "/", "m2", dataSize), m2)
precis(m2)
#plot(m2)
#pairs(m2)

#let's try having them both in there
m3 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a + bDR*daysRelativeToRejection + bR*afterRejection, 
    a ~ dnorm(0,100),
    bDR ~ dnorm(0,1),
    bR ~ dnorm(0,1)
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(a=0, bR=0, bDR=0)
)
save(file = paste0(dataSize, "/", "m3", dataSize), m3)
precis(m3)
#plot(m3)
#pairs(m3)

compare(m1, m2, m3)

#ok, now let's have m1, but with patient-specific effects
m4 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a[patient_id] + bR[patient_id]*afterRejection, 
    a[patient_id] ~ dnorm(muA,sigA),
    muA ~ dnorm(0,5),
    sigA ~ dcauchy(0,3),
    bR[patient_id] ~ dnorm(muBR,sigbR),
    muBR ~ dnorm(0,1),
    sigbR ~ dcauchy(0,3)
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(muA=0, sigA=1, muBR=0, sigbR=1)
)
save(file = paste0(dataSize, "/", "m4", dataSize), m4)
precis(m4, depth = 2)
#plot(m4)
#pairs(m4)

#ok, now let's have m1, but with protein-specific effects
m5 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a[protein_id] + bR[protein_id]*afterRejection, 
    a[protein_id] ~ dnorm(muA,sigA),
    muA ~ dnorm(0,5),
    sigA ~ dcauchy(0,3),
    bR[protein_id] ~ dnorm(muBR,sigbR),
    muBR ~ dnorm(0,1),
    sigbR ~ dcauchy(0,3)
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(muA=0, sigA=1, muBR=0, sigbR=1)
)
save(file = paste0(dataSize, "/", "m5", dataSize), m5)
precis(m5, depth = 2)
#plot(m5)
#pairs(m5)

# #ok, now let's have m1, but with patient-specific and protein-specific effects; this model is nonidentifiable
# #because muAj and muAi only inflence the likelihood through their sum
# nChains <- 1
# m6 <- map2stan(
#   alist(
#     spectral_count ~ dpois( lambda ),
#     log(lambda) <- a + bR*afterRejection, 
#     a <- Ab + Aj[protein_id] + Ai[patient_id],
#     Ab ~ dnorm(0,20),
#     Aj[protein_id] ~ dnorm(muAj,sigAj),
#     Ai[patient_id] ~ dnorm(muAi,sigAi),
#     c(muAj, muAi) ~ dnorm(0,1),
#     c(sigAj, sigAi) ~ dcauchy(0,1),
#     bR <- bRb + bRj[protein_id] + bRi[patient_id],
#     bRb ~ dnorm(0,2),
#     bRj[protein_id] ~ dnorm(muBRj,sigbRj),
#     bRi[patient_id] ~ dnorm(muBRi,sigbRi),
#     c(muBRj, muBRi) ~ dnorm(0,1),
#     c(sigbRj, sigbRi) ~ dcauchy(0,1)
#   ) ,
#   data= d,
#   iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
#   start = list(muAj=0, muAi=0, sigAj=1, sigAi=1, muBRj=0, muBRi=0, sigbRj=1, sigbRi=1, bRb=0, Ab=0)
# )
# save(file = paste0(dataSize, "/", "m6", dataSize), m6)
# precis(m6, depth = 1)
# #plot(m6)
# #pairs(m6)



#ok, now let's have m1, but with patient-specific and protein-specific effects, but have them be deviations from some mean effect (i.e. an intercept)
nChains <- 1
m7 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a + bR*afterRejection, 
    a <- Ab + Aj[protein_id] + Ai[patient_id],
    Ab ~ dnorm(0,20),
    Aj[protein_id] ~ dnorm(0,sigAj),
    Ai[patient_id] ~ dnorm(0,sigAi),
    c(sigAj, sigAi) ~ dcauchy(0,1),
    bR <- bRb + bRj[protein_id] + bRi[patient_id],
    bRb ~ dnorm(0,2),
    bRj[protein_id] ~ dnorm(0,sigbRj),
    bRi[patient_id] ~ dnorm(0,sigbRi),
    c(sigbRj, sigbRi) ~ dcauchy(0,1)
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigAj=1, sigAi=1, sigbRj=1, sigbRi=1, bRb=0, Ab=0)
)
save(file = paste0(dataSize, "/", "m7", dataSize), m7)
precis(m7, depth = 2)
#plot(m7)
#pairs(m7)

#check overdispersion through posterior predictive simulation
m7p <- link(m7, n = 2000)
dens(d$spectral_count, xlim = c(0,50))
dens(m7p$lambda[1,], add = T)
hist(sapply(1:length(d$spectral_count), function(x) sum(rpois(2000,m7p$lambda[,x]) > d$spectral_count[x])/2000))

#let's try to plot some of these effects for m7, e.g. expected counts
m7s <- extract.samples(m7, clean.names = F)
proteinEffects <- m7s$bRj
d$protein_id <- coerce_index(d$protein)
proteinKey <- data.frame(as.character(d$protein), d$protein_id)
proteinKey <- proteinKey[order(proteinKey$d.protein_id),]
proteinKey <- unique(proteinKey)
dens(exp(m7s$bRb + proteinEffects[,1]) * 100, xlim = c(0,6e2), ylim = c(0,0.025), xlab = "Percent Change in Intercept with Rejection Event", col = 0)
num <- 0
for(j in 1:length(proteinEffects[1,])){
  percChange <- exp(m7s$bRb + proteinEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  if(probPos > 0.95){
    num <- num + 1
    cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
  }
}
num

#let's have a linear affect on proximity to the days relative to rejection variable -- i.e. reflection at the rejection event
# nChains <- 1
# m8 <- map2stan(
#   alist(
#     spectral_count ~ dpois( lambda ),
#     log(lambda) <- a + bRpiv*fabs(daysRelativeToRejection), 
#     
#     a <- Ab + Aj[protein_id] + Ai[patient_id],
#     Ab ~ dnorm(0,4),
#     Aj[protein_id] ~ dnorm(0,sigAj),
#     Ai[patient_id] ~ dnorm(0,sigAi),
#     c(sigAj, sigAi) ~ dcauchy(0,1),
#        
#     bRpiv <- bRpivb + bRpivj[protein_id] + bRpivi[patient_id],
#     bRpivj[protein_id] ~ dnorm(0,sigbRpivj),
#     bRpivi[patient_id] ~ dnorm(0,sigbRpivi),
#     
#     bRpivb ~ dnorm(0,0.03),
#     c(sigbRpivj, sigbRpivi) ~ dcauchy(0,0.1)
#   ) ,
#   data= d,
#   iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
#   start = list(sigAj=1, sigAi=1, sigbRpivj=0.01, sigbRpivi=0.01, bRpivb=0, Ab=0)
# )
# save(file = paste0(dataSize, "/", "m8", dataSize), m8)
# precis(m8, depth = 2)
# #plot(m8)
# #pairs(m8)

#linear effect of days through the rejection event, without the insane prior and starting state
nChains <- 1
m8.1 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a + bRpiv*daysRelativeToRejection, 
    
    a <- Ab + Aj[protein_id] + Ai[patient_id],
    Ab ~ dnorm(0,3),
    Aj[protein_id] ~ dnorm(0,sigAj),
    Ai[patient_id] ~ dnorm(0,sigAi),
    c(sigAj, sigAi) ~ dcauchy(0,1),
    
    bRpiv <- bRpivb + bRpivj[protein_id] + bRpivi[patient_id],
    bRpivj[protein_id] ~ dnorm(0,sigbRpivj),
    bRpivi[patient_id] ~ dnorm(0,sigbRpivi),
    
    bRpivb ~ dnorm(0,0.03),
    c(sigbRpivj, sigbRpivi) ~ dcauchy(0,0.1)
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigAj=1, sigAi=1, sigbRpivj=0.01, sigbRpivi=0.01, bRpivb=0, Ab=2)
)
save(file = paste0(dataSize, "/", "m8.1", dataSize), m8.1)
precis(m8.1, depth = 2)
#plot(m8.1)
#pairs(m8.1)
m8.1s <- extract.samples(m8.1)

#separate slopes for before and after the rejection event, but same intercept at the rejection event
nChains <- 1
m9 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a + bR0*daysRelativeToRejection*(1-afterRejection) + bR1*daysRelativeToRejection*afterRejection, 
    
    a <- Ab + Aj[protein_id] + Ai[patient_id],
    Ab ~ dnorm(0,4),
    Aj[protein_id] ~ dnorm(0,sigAj),
    Ai[patient_id] ~ dnorm(0,sigAi),
    c(sigAj, sigAi) ~ dcauchy(0,1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,0.03),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,0.1),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,0.03),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,0.1)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigAj=1, sigAi=1, sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, Ab=0)
)
save(file = paste0(dataSize, "/", "m9", dataSize), m9)
precis(m9, depth = 2)
#plot(m9)
#pairs(m9)

m9p <- link(m9, n = 2000)
hist(d$spectral_count, xlim = c(0,50), breaks = 200, col = rgb(0,0,1,alpha = 0.5))
hist(rpois(n = length(m9p$lambda[1,]), lambda = m9p$lambda[1,]), add = T, col = rgb(1,0,0,alpha = 0.5), breaks = 200)
hist(sapply(1:length(d$spectral_count), function(x) sum(rpois(2000,m9p$lambda[,x]) >= d$spectral_count[x])/2000))

#plot raw residuals from mean expectation
plot(sapply(1:length(m9p$lambda[1,]), function(x) mean(m9p$lambda[,x])), 
     sapply(1:length(m9p$lambda[1,]), function(x) mean(m9p$lambda[,x])) - d$spectral_count,
     col = rgb(0,0,0,0.3)) 

#plot pearson residuals?
plot(sapply(1:length(m9p$lambda[1,]), function(x) mean(m9p$lambda[,x])), 
     (sapply(1:length(m9p$lambda[1,]), function(x) mean(m9p$lambda[,x])) - d$spectral_count) / sqrt(sapply(1:length(m9p$lambda[1,]), function(x) mean(m9p$lambda[,x]))),
     xlab = "Predicted Rate", ylab = "Standardized Residual from Predicted Rate")
hist((sapply(1:length(m9p$lambda[1,]), function(x) mean(m9p$lambda[,x])) - d$spectral_count) / sqrt(sapply(1:length(m9p$lambda[1,]), function(x) mean(m9p$lambda[,x]))),
     breaks = 100, xlim = c(-5, 5), freq = F)
lines(seq(-5, 5, length.out = 200), dnorm(x = seq(-5, 5, length.out = 200), mean = 0, sd = 1))

#let's try to plot some of these effects for m9, e.g. expected counts per day before rejection
load("m9c")
m9s <- extract.samples(m9, clean.names = F)

proteinEffects <- m9s$bR0j
d$protein_id <- coerce_index(d$protein)
proteinKey <- data.frame(as.character(d$protein), d$protein_id)
proteinKey <- proteinKey[order(proteinKey$d.protein_id),]
proteinKey <- unique(proteinKey)
dens(exp(m9s$bR0b + proteinEffects[,1]) * 100, xlim = c(99.5,100.5), xlab = "Percent Change in Spectral Count per Day", col = 0, ylim = c(0,20))
num <- 0
for(j in 1:length(proteinEffects[1,])){
  percChange <- exp(m9s$bR0b + proteinEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  if(probPos > 0.95){
    num <- num + 1
    cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
  }
}
num

#just look at data
plot(sort(unique(d$daysRelativeToRejection)),
     sapply(1:length(sort(unique(d$daysRelativeToRejection))), function(x) mean(d$spectral_count[d$daysRelativeToRejection == sort(unique(d$daysRelativeToRejection))[x]])),
     xlab = "days relative to rejection", ylab = "mean spectral count")

#after rejection
proteinEffects <- m9s$bR1j
dens(exp(m9s$bR1b + proteinEffects[,1]) * 100, xlim = c(99.5,100.5), xlab = "Percent Change in Intercept with Day", col = 0)
num <- 0
for(j in 1:length(proteinEffects[1,])){
  percChange <- exp(m9s$bR1b + proteinEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  if(probPos < 0.95){
    num <- num + 1
    cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
  }
}
num

#individual effects
patientEffects <- m9s$bR0i
d$patient_id <- coerce_index(d$patient)
patientKey <- data.frame(as.character(d$patient), d$patient_id)
patientKey <- patientKey[order(patientKey$d.patient_id),]
patientKey <- unique(patientKey)
dens(exp(m9s$bR0b + patientEffects[,1]) * 100, xlim = c(99.8,100.2), xlab = "Percent Change in Spectral Count per Day Going Backward", col = 0, ylim = c(0,300))
for(j in 1:length(patientEffects[1,])){
  percChange <- exp(m9s$bR0b + patientEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  if(probPos > -1){
    cat(paste0(as.character(patientKey[patientKey$d.patient_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
    kde <- density(percChange)
    text(x = kde$x[which.max(kde$y)], y = max(kde$y)*1.05, labels = as.character(patientKey[patientKey$d.patient_id == j, 1]))
  }
}


patientEffects <- m9s$bR1i
dens(exp(m9s$bR1b + patientEffects[,1]) * 100, xlim = c(99.5,100.5), xlab = "Percent Change in Spectral Count per Day Going Forward", col = 0, ylim = c(0,200))
for(j in 1:length(patientEffects[1,])){
  percChange <- exp(m9s$bR1b + patientEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  if(probPos > -1){
    cat(paste0(as.character(patientKey[patientKey$d.patient_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
    kde <- density(percChange)
    text(x = kde$x[which.max(kde$y)], y = max(kde$y)*1.05, labels = as.character(patientKey[patientKey$d.patient_id == j, 1]))
  }
}

#plot lines for uncertainty in expected rate through time relative to rejection
str(m9s)
# range of days relative to rejection c(-1300, 1100)
i <- 1
plot(c(-1300,0), exp(m9s$Ab[i] + c(-1300,0) * m9s$bR0b[i]), type = "l", col = rgb(0,0,0,0), xlim = c(-1300, 1100), ylim = c(0,10), 
     xlab = "days relative to rejection event", ylab = "rate of spectral count for the average protein/patient")
abline(v = 0, col = 2)
for(i in 1:length(m9s$Ab)){
  lines(c(-1300,0), exp(m9s$Ab[i] + c(-1300,0) * m9s$bR0b[i]), type = "l", col = rgb(0,0,0,0.02), xlim = c(-1300, 1100), ylim = c(0,10))
  lines(c(0:1100), exp(m9s$Ab[i] + c(0:1100) * m9s$bR1b[i]), type = "l", col = rgb(0,0,0,0.02), xlim = c(-1300, 1100), ylim = c(0,10))
}
`%roundup%` <- function(e1, e2) {e1[e1 < e2] <- e2; e1}
points(d$daysRelativeToRejection + rnorm(length(d$daysRelativeToRejection), 0, 20), d$spectral_count + rnorm(length(d$spectral_count), 0, 1)  %roundup% 0, col = rgb(1,0,0,0.2))

#extracting stan code?
cat(m9@model)
m9_code <- 
  "
data{
int<lower=1> N;
int<lower=1> N_protein_id;
int<lower=1> N_patient;
int spectral_count[N];
real daysRelativeToRejection[N];
real afterRejection[N];
int patient[N];
int protein_id[N];
}
parameters{
real<lower=0> sigAj;
real<lower=0> sigAi;
real<lower=0> sigbR0j;
real<lower=0> sigbR0i;
real bR0b;
real<lower=0> sigbR1j;
real<lower=0> sigbR1i;
real bR1b;
real Ab;
vector[N_protein_id] Aj;
vector[N_patient] Ai;
vector[N_protein_id] bR0j;
vector[N_patient] bR0i;
vector[N_protein_id] bR1j;
vector[N_patient] bR1i;
}
model{
vector[N] bR0;
vector[N] a;
vector[N] lambda;
sigbR1i ~ cauchy( 0 , 0.1 );
sigbR1j ~ cauchy( 0 , 0.1 );
bR1b ~ normal( 0 , 0.03 );
bR1i ~ normal( 0 , sigbR1i );
bR1j ~ normal( 0 , sigbR1j );
for ( i in 1:N ) {
bR1[i] = bR1b + bR1j[protein_id[i]] + bR1i[patient[i]];
}
sigbR0i ~ cauchy( 0 , 0.1 );
sigbR0j ~ cauchy( 0 , 0.1 );
bR0b ~ normal( 0 , 0.03 );
bR0i ~ normal( 0 , sigbR0i );
bR0j ~ normal( 0 , sigbR0j );
for ( i in 1:N ) {
bR0[i] = bR0b + bR0j[protein_id[i]] + bR0i[patient[i]];
}
sigAi ~ cauchy( 0 , 1 );
sigAj ~ cauchy( 0 , 1 );
Ai ~ normal( 0 , sigAi );
Aj ~ normal( 0 , sigAj );
Ab ~ normal( 0 , 4 );
for ( i in 1:N ) {
a[i] = Ab + Aj[protein_id[i]] + Ai[patient[i]];
}
for ( i in 1:N ) {
lambda[i] = a[i] + bR0[i] * fabs(daysRelativeToRejection[i]) * (1 - afterRejection[i]) +      bR1[i] * fabs(daysRelativeToRejection[i]) * afterRejection[i];
}
spectral_count ~ poisson_log( lambda );
}
generated quantities{
vector[N] bR1;
vector[N] bR0;
vector[N] a;
vector[N] lambda;
real dev;
dev = 0;
for ( i in 1:N ) {
bR1[i] = bR1b + bR1j[protein_id[i]] + bR1i[patient[i]];
}
for ( i in 1:N ) {
bR0[i] = bR0b + bR0j[protein_id[i]] + bR0i[patient[i]];
}
for ( i in 1:N ) {
a[i] = Ab + Aj[protein_id[i]] + Ai[patient[i]];
}
for ( i in 1:N ) {
lambda[i] = a[i] + bR0[i] * fabs(daysRelativeToRejection[i]) * (1 - afterRejection[i]) +      bR1[i] * fabs(daysRelativeToRejection[i]) * afterRejection[i];
}
dev = dev + (-2)*poisson_log_lpmf( spectral_count | lambda );
}
"

m9mod_code <- 
  "data{
int<lower=1> N;
int<lower=1> N_protein_id;
int<lower=1> N_patient;
int spectral_count[N];
int piv[N];
real daysRelativeToRejection[N];
real afterRejection[N];
int patient[N];
int protein_id[N];
}
parameters{
real<lower=0> sigAj;
real<lower=0> sigAi;
real<lower=0> sigbR0j;
real<lower=0> sigbR0i;
real bR0b;
real<lower=0> sigbR1j;
real<lower=0> sigbR1i;
real bR1b;
real Ab;
vector[N_protein_id] Aj;
vector[N_patient] Ai;
vector[N_protein_id] bR0j;
vector[N_patient] bR0i;
vector[N_protein_id] bR1j;
vector[N_patient] bR1i;
}
model{
vector[N] bR1;
vector[N] bR0;
vector[N] a;
vector[N] lambda;
sigbR1i ~ cauchy( 0 , 0.1 );
sigbR1j ~ cauchy( 0 , 0.1 );
bR1b ~ normal( 0 , 0.03 );
bR1i ~ normal( 0 , sigbR1i );
bR1j ~ normal( 0 , sigbR1j );
for ( i in 1:N ) {
bR1[i] = bR1b + bR1j[protein_id[i]] + bR1i[patient[i]];
}
sigbR0i ~ cauchy( 0 , 0.1 );
sigbR0j ~ cauchy( 0 , 0.1 );
bR0b ~ normal( 0 , 0.03 );
bR0i ~ normal( 0 , sigbR0i );
bR0j ~ normal( 0 , sigbR0j );
for ( i in 1:N ) {
bR0[i] = bR0b + bR0j[protein_id[i]] + bR0i[patient[i]];
}
sigAi ~ cauchy( 0 , 1 );
sigAj ~ cauchy( 0 , 1 );
Ai ~ normal( 0 , sigAi );
Aj ~ normal( 0 , sigAj );
Ab ~ normal( 0 , 4 );
for ( i in 1:N ) {
a[i] = Ab + Aj[protein_id[i]] + Ai[patient[i]];
}
for ( i in 1:N ) {
lambda[i] = a[i] + bR0[i] * fabs(daysRelativeToRejection[i]) * (1 - afterRejection[i]) +      bR1[i] * fabs(daysRelativeToRejection[i]) * afterRejection[i];
}
spectral_count ~ poisson_log( lambda );
}
generated quantities{
vector[N] bR1;
vector[N] bR0;
vector[N] a;
vector[N] lambda;
real dev;
dev = 0;
for ( i in 1:N ) {
bR1[i] = bR1b + bR1j[protein_id[i]] + bR1i[patient[i]];
}
for ( i in 1:N ) {
bR0[i] = bR0b + bR0j[protein_id[i]] + bR0i[patient[i]];
}
for ( i in 1:N ) {
a[i] = Ab + Aj[protein_id[i]] + Ai[patient[i]];
}
for ( i in 1:N ) {
lambda[i] = a[i] + bR0[i] * fabs(daysRelativeToRejection[i]) * (1 - afterRejection[i]) +      bR1[i] * fabs(daysRelativeToRejection[i]) * afterRejection[i];
}
dev = dev + (-2)*poisson_log_lpmf( spectral_count | lambda );
}"

N <- 1000 # number of cases
N_miss <- 100 # number missing values
x_baserate <- 0.25 # prob x==1 in total sample
a <- 0 # intercept in y ~ N( a+b*x , 1 )
b <- 1 # slope in y ~ N( a+b*x , 1 )

# simulate data
x <- sample( 0:1 , size=N , replace=TRUE , prob=c(1-x_baserate,x_baserate) )
i_miss <- sample( 1:N , size=N_miss )
x_obs <- x
x_obs[i_miss] <- (-1)
x_NA <- x_obs
x_NA[i_miss] <- NA
x_miss <- ifelse( 1:N %in% i_miss , 1 , 0 )
y <- rnorm( N , a + b*x , 1 )

#separate slopes for before and after the rejection event, seperate intercepts at the rejection event
nChains <- 1
m13 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0*(1-afterRejection) + a1*afterRejection + bR0*daysRelativeToRejection*(1-afterRejection) + bR1*daysRelativeToRejection*afterRejection, 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id],
    A0b ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,1),
    
    a1 <- A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,0.03),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,0.1),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,0.03),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,0.1)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
save(file = paste0(dataSize, "/", "m13", dataSize), m13)
precis(m13, depth = 2)
#plot(m13)
#pairs(m13)

#let's try to plot some of these effects for m13, e.g. expected counts per day before rejection
load("default/m13default")
m13s <- extract.samples(m13, clean.names = F)
proteinEffects <- m13s$bR0j
proteinKey <- data.frame(as.character(d$protein), d$protein_id)
proteinKey <- proteinKey[order(proteinKey$d.protein_id),]
proteinKey <- unique(proteinKey)
dens(exp(m13s$bR0b + proteinEffects[,1]) * 100, xlim = c(99.5,100.5), xlab = "Percent Change in Spectral Count per Day", col = 0, ylim = c(0,20))
num <- 0
for(j in 1:length(proteinEffects[1,])){
  percChange <- exp(m13s$bR0b + proteinEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  if(probPos > 0.95 | probPos < 0.05){
    num <- num + 1
    cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
  }
}
num

#after rejection
proteinEffects <- m13s$bR1j
dens(exp(m13s$bR1b + proteinEffects[,1]) * 100, xlim = c(99.5,100.5), xlab = "Percent Change in Intercept with Day", col = 0, ylim = c(0,12))
num <- 0
for(j in 1:length(proteinEffects[1,])){
  percChange <- exp(m13s$bR1b + proteinEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  if(probPos > 0.95 | probPos < 0.05){
    num <- num + 1
    cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
  }
}
num

#individual effects
patientEffects <- m13s$bR0i
d$patient_id <- coerce_index(d$patient)
patientKey <- data.frame(as.character(d$patient), d$patient_id)
patientKey <- patientKey[order(patientKey$d.patient_id),]
patientKey <- unique(patientKey)
dens(exp(m13s$bR0b + patientEffects[,1]) * 100, xlim = c(99.8,100.2), xlab = "Percent Change in Spectral Count per Day Going Backward", col = 0, ylim = c(0,300))
for(j in 1:length(patientEffects[1,])){
  percChange <- exp(m13s$bR0b + patientEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  if(probPos > -1){
    cat(paste0(as.character(patientKey[patientKey$d.patient_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
    kde <- density(percChange)
    text(x = kde$x[which.max(kde$y)], y = max(kde$y)*1.05, labels = as.character(patientKey[patientKey$d.patient_id == j, 1]))
  }
}


patientEffects <- m13s$bR1i
dens(exp(m13s$bR1b + patientEffects[,1]) * 100, xlim = c(99.5,100.5), xlab = "Percent Change in Spectral Count per Day Going Forward", col = 0, ylim = c(0,200))
for(j in 1:length(patientEffects[1,])){
  percChange <- exp(m13s$bR1b + patientEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  if(probPos > -1){
    cat(paste0(as.character(patientKey[patientKey$d.patient_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
    kde <- density(percChange)
    text(x = kde$x[which.max(kde$y)], y = max(kde$y)*1.05, labels = as.character(patientKey[patientKey$d.patient_id == j, 1]))
  }
}

#plot lines for uncertainty in expected rate through time relative to rejection

#mean estimate
par(mar = c(4,4,2,2))
plot(c(-1300,0), exp(mean(m13s$A0b) + c(-1300,0) * mean(m13s$bR0b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), 
     xlab = "days relative to rejection event", ylab = "rate of spectral count for the average protein/patient", lwd = 2)
lines(c(0:1100), exp(mean(m13s$A1b) + c(0:1100) * mean(m13s$bR1b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
title("Expected Rates of Spectral Count Given Days Relative to Acute Rejection Event")
#89% hpdi
lines(c(-1300,0), exp(HPDI(prob = 0.5, m13s$A0b)[[1]] + c(-1300,0) * HPDI(prob = 0.5, m13s$bR0b)[[2]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(0:1100), exp(HPDI(prob = 0.5, m13s$A1b)[[1]] + c(0:1100) * HPDI(prob = 0.5, m13s$bR1b)[[1]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(-1300,0), exp(HPDI(prob = 0.5, m13s$A0b)[[2]] + c(-1300,0) * HPDI(prob = 0.5, m13s$bR0b)[[1]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(0:1100), exp(HPDI(prob = 0.5, m13s$A1b)[[2]] + c(0:1100) * HPDI(prob = 0.5, m13s$bR1b)[[2]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
#rejection date
abline(v = 0, col = 2)
#samples from posterior
for(i in seq(1, length(m13s$A0b), length.out = 1000)){
  lines(c(-1300,0), exp(m13s$A0b[i] + c(-1300,0) * m13s$bR0b[i]), type = "l", col = rgb(0,0,0,0.02), xlim = c(-1300, 1100), ylim = c(0,10))
  lines(c(0:1100), exp(m13s$A1b[i] + c(0:1100) * m13s$bR1b[i]), type = "l", col = rgb(0,0,0,0.02), xlim = c(-1300, 1100), ylim = c(0,10))
}
#boxplots of data
times <- sort(unique(d$daysRelativeToRejection))
for(i in 1:length(times)){
  boxplot(d$spectral_count[d$daysRelativeToRejection == times[i]], add = T, at = times[i], boxwex=20, outline = F, col.axis = rgb(0,0,0,0))
}


#same slopes for before and after the rejection event, seperate intercepts at the rejection event
nChains <- 1
m14 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0*(1-afterRejection) + a1*afterRejection + bR*daysRelativeToRejection, 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id],
    A0b ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,1),
    
    a1 <- A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,1),
    
    bR <- bRb + bRj[protein_id] + bRi[patient_id],
    bRj[protein_id] ~ dnorm(0,sigbRj),
    bRi[patient_id] ~ dnorm(0,sigbRi),
    bRb ~ dnorm(0,0.03),
    c(sigbRj, sigbRi) ~ dcauchy(0,0.1)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigbRj=0.01, sigbRi=0.01, bRb=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
save(file = paste0(dataSize, "/", "m14", dataSize), m14)
precis(m14, depth = 2)
#plot(m14)
#pairs(m14)

#similar slopes for before and after the rejection event, seperate intercepts at the rejection event
nChains <- 1
m15 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0*(1-afterRejection) + a1*afterRejection + bR*daysRelativeToRejection*(1-afterRejection) + bR*daysRelativeToRejection*afterRejection, 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id],
    A0b ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,1),
    
    a1 <- A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,1),
    
    bR <- bRb + bRj[protein_id] + bRi[patient_id] + bRAfter*afterRejection,
    bRj[protein_id] ~ dnorm(0,sigbRj),
    bRi[patient_id] ~ dnorm(0,sigbRi),
    bRb ~ dnorm(0,0.03),
    c(sigbRj, sigbRi) ~ dcauchy(0,0.1),
    bRAfter ~ dnorm(0, sigbRAfter),
    sigbRAfter ~ dcauchy(0,0.1)
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigbRj=0.01, sigbRi=0.01, sigbRAfter = 0.01, bRb=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
save(file = paste0(dataSize, "/", "m15", dataSize), m15)
precis(m15, depth = 2)
#plot(m15)
#pairs(m15)

#similar slopes for before and after the rejection event, but same intercept at the rejection event
nChains <- 1
m16 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a + bR*daysRelativeToRejection*(1-afterRejection) + bR*daysRelativeToRejection*afterRejection, 
    
    a <- Ab + Aj[protein_id] + Ai[patient_id],
    Ab ~ dnorm(0,4),
    Aj[protein_id] ~ dnorm(0,sigAj),
    Ai[patient_id] ~ dnorm(0,sigAi),
    c(sigAj, sigAi) ~ dcauchy(0,1),
    
    bR <- bRb + bRj[protein_id] + bRi[patient_id] + bRAfter*afterRejection,
    bRj[protein_id] ~ dnorm(0,sigbRj),
    bRi[patient_id] ~ dnorm(0,sigbRi),
    bRb ~ dnorm(0,0.03),
    c(sigbRj, sigbRi) ~ dcauchy(0,0.1),
    bRAfter ~ dnorm(0, sigbRAfter),
    sigbRAfter ~ dcauchy(0,0.1)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigAj=1, sigAi=1, sigbRj=0.01, sigbRi=0.01, sigbRAfter = 0.01, Ab=0)
)
save(file = paste0(dataSize, "/", "m16", dataSize), m16)
precis(m16, depth = 2)
#plot(m16)
#pairs(m16)

#separate slopes for before and after the rejection event, and the before takes a nonlinear (exponential) relationship with days
nChains <- 1
m10 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a*afterRejection + bR0*exp(daysRelativeToRejection*expSlopeTune)*(1-afterRejection) + bR1*daysRelativeToRejection*afterRejection, 
    expSlopeTune ~ dunif(0,2),
    
    a <- Ab + Aj[protein_id] + Ai[patient_id],
    Ab ~ dnorm(0,10),
    Aj[protein_id] ~ dnorm(0,sigAj),
    Ai[patient_id] ~ dnorm(0,sigAi),
    c(sigAj, sigAi) ~ dcauchy(0,1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,0.25),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,0.5),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,0.25),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,0.5)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigAj=1, sigAi=1, sigbR0j=1, sigbR0i=1, bR0b=0, sigbR1j=1, sigbR1i=1, bR1b=0, Ab=0)
)
save(file = paste0(dataSize, "/", "m10", dataSize), m10)
precis(m10, depth = 2)
#plot(m10)
#pairs(m10)

#separate slopes for before and after the rejection event, and the before takes a nonlinear (exponential) relationship with days, 
#and both relationships are constrained to have the same value at 0 (the rejection event)
nChains <- 1
m10.1 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- bR0*afterRejection + bR0*exp(daysRelativeToRejection*expSlopeTune)*(1-afterRejection) + bR1*daysRelativeToRejection*afterRejection, 
    expSlopeTune ~ dunif(0,1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,4),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,1),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,0.03),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,0.1)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigbR0j=.1, sigbR0i=.1, bR0b=2, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, expSlopeTune = 0.1)
)
save(file = paste0(dataSize, "/", "m10.1", dataSize), m10.1)
precis(m10.1, depth = 2)
#plot(m10.1)
#pairs(m10.1)

#separate slopes for before and after the rejection event, and the before takes a nonlinear (sqrt) relationship with days
nChains <- 1
m11 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a + bR0*sqrt(fabs(daysRelativeToRejection))*(1-afterRejection) + bR1*fabs(daysRelativeToRejection)*afterRejection, 
    
    a <- Ab + Aj[protein_id] + Ai[patient_id],
    Ab ~ dnorm(0,4),
    Aj[protein_id] ~ dnorm(0,sigAj),
    Ai[patient_id] ~ dnorm(0,sigAi),
    c(sigAj, sigAi) ~ dcauchy(0,1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,0.03),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,0.1),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,0.03),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,0.1)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigAj=1, sigAi=1, sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, Ab=0)
)
save(file = paste0(dataSize, "/", "m11", dataSize), m11)
precis(m11, depth = 2)
#plot(m11)
#pairs(m11)

#separate slopes for before and after the rejection event, and the before takes a nonlinear (squared) relationship with days
nChains <- 1
m12 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a + bR0*(fabs(daysRelativeToRejection)^2)*(1-afterRejection) + bR1*fabs(daysRelativeToRejection)*afterRejection, 
    
    a <- Ab + Aj[protein_id] + Ai[patient_id],
    Ab ~ dnorm(0,10),
    Aj[protein_id] ~ dnorm(0,sigAj),
    Ai[patient_id] ~ dnorm(0,sigAi),
    c(sigAj, sigAi) ~ dcauchy(0,1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,0.25),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,0.5),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,0.25),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,0.5)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigAj=1, sigAi=1, sigbR0j=1, sigbR0i=1, bR0b=0, sigbR1j=1, sigbR1i=1, bR1b=0, Ab=0)
)
save(file = paste0(dataSize, "/", "m12", dataSize), m12)
precis(m12, depth = 2)
#plot(m12)
#pairs(m12)

#separate slopes for before and after the rejection event, seperate "intercepts" at the rejection event, intercept with zero slope before the transplant event,
#express slope in change per unit day between transplant and rejection event
nChains <- 1
m23c2 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0*(1-afterRejection) + a1*afterRejection + bR0*daysRelativeToTransplant*(1-afterRejection)*afterTransplant + bR1*daysRelativeToRejection*afterRejection, 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id],
    A0b ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,1),
    
    a1 <- A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,0.03),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,0.01),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,0.03),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,0.01)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
save(file = paste0(dataSize, "/", "m23c2", dataSize), m23c2)
precis(m23c2, depth = 2)

#let's try to plot some of these effects for m23, e.g. expected counts per day before rejection
load("default/m23spdefault")

m23p <- link(m23sp, n = 2000)
hist(d$spectral_count, xlim = c(0,50), breaks = 200, col = rgb(0,0,1,alpha = 0.5))
hist(rpois(n = length(m23p$lambda[1,]), lambda = m23p$lambda[1,]), add = T, col = rgb(1,0,0,alpha = 0.5), breaks = 200)
hist(sapply(1:length(d$spectral_count), function(x) ppois(d$spectral_count[x], m23p$lambda[,x])))
hist(sapply(1:length(d$spectral_count), function(x) sum(rpois(2000,m23p$lambda[,x]) >= d$spectral_count[x])/2000))

#plot raw residuals from mean expectation
plot(sapply(1:length(m23p$lambda[1,]), function(x) mean(m23p$lambda[,x])), 
     sapply(1:length(m23p$lambda[1,]), function(x) mean(m23p$lambda[,x])) - d$spectral_count,
     col = rgb(0,0,0,0.3)) 

#plot pearson residuals?
plot(sapply(1:length(m23p$lambda[1,]), function(x) mean(m23p$lambda[,x])), 
     (sapply(1:length(m23p$lambda[1,]), function(x) mean(m23p$lambda[,x])) - d$spectral_count) / sqrt(sapply(1:length(m23p$lambda[1,]), function(x) mean(m23p$lambda[,x]))),
     xlab = "Predicted Rate", ylab = "Standardized Residual from Predicted Rate")
hist((sapply(1:length(m23p$lambda[1,]), function(x) mean(m23p$lambda[,x])) - d$spectral_count) / sqrt(sapply(1:length(m23p$lambda[1,]), function(x) mean(m23p$lambda[,x]))),
     breaks = 50, xlim = c(-5, 5), freq = F)
lines(seq(-5, 5, length.out = 200), dnorm(x = seq(-5, 5, length.out = 200), mean = 0, sd = 1))

m23s <- extract.samples(m23sp, clean.names = F)
# m23s <- extract.samples(m23, clean.names = F)
proteinEffects <- m23s$bR0j
proteinKey <- data.frame(as.character(d$protein), d$protein_id)
proteinKey <- proteinKey[order(proteinKey$d.protein_id),]
proteinKey <- unique(proteinKey)
hist(sapply(1:length(proteinEffects[1,]), function (x) mean(proteinEffects[,x])), breaks = 200)
dens(exp(m23s$bR0b + proteinEffects[,1]) * 100, xlim = c(98,102), xlab = "Percent Change in Spectral Count per Day", col = 0, ylim = c(0,2))
num <- 0
percChanges <- 0
for(j in 1:length(proteinEffects[1,])){
  percChange <- exp(m23s$bR0b + proteinEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  # cat(paste0(probPos, " "))
  if(probPos > 0.95 | probPos < 0.05){
    num <- num + 1
    cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
  }
}
num

#after rejection
proteinEffects <- m23s$bR1j
hist(sapply(1:length(proteinEffects[1,]), function (x) mean(proteinEffects[,x])), breaks = 200)
dens(exp(m23s$bR1b + proteinEffects[,1]) * 100, xlim = c(99.5,100.5), xlab = "Percent Change in Intercept with Day", col = 0, ylim = c(0,12))
num <- 0
for(j in 1:length(proteinEffects[1,])){
  percChange <- exp(m23s$bR1b + proteinEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  if(probPos > 0.95 | probPos < 0.05){
    num <- num + 1
    cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
  }
}
num

# probing specific proteins
targetProteins <- c("troponin", "prohibitin", "plakoglobin")
targetProteins <- proteinKey[unlist(sapply(1:length(targetProteins), function(x) which(grepl(targetProteins[x], proteinKey$as.character.d.protein)))),]
proteinEffects <- m23s$bR0j
for(j in 1:length(targetProteins[,1])){
  percChange <- exp(m23s$bR0b + proteinEffects[,targetProteins[j,2]]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == targetProteins[j,2], 1]), ", probability positive = ", probPos * 100, "%, mean change = ", round(mean(percChange), 3), "%\n"))
}
proteinEffects <- m23s$bR1j
for(j in 1:length(targetProteins[,1])){
  percChange <- exp(m23s$bR1b + proteinEffects[,targetProteins[j,2]]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == targetProteins[j,2], 1]), ", probability positive = ", probPos * 100, "%, mean change = ", round(mean(percChange), 3), "%\n"))
}

#individual effects
patientEffects <- m23s$bR0i
patientKey <- data.frame(as.character(d$patient), d$patient_id)
patientKey <- patientKey[order(patientKey$d.patient_id),]
patientKey <- unique(patientKey)
dens(exp(m23s$bR0b + patientEffects[,1]) * 100, xlim = c(99,101), xlab = "Percent Change in Spectral Count per Day Going Backward", col = 0, ylim = c(0,50))
for(j in 1:length(patientEffects[1,])){
  percChange <- exp(m23s$bR0b + patientEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  if(probPos > -1){
    cat(paste0(as.character(patientKey[patientKey$d.patient_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
    kde <- density(percChange)
    text(x = kde$x[which.max(kde$y)], y = max(kde$y)*1.05, labels = as.character(patientKey[patientKey$d.patient_id == j, 1]))
  }
}


patientEffects <- m23s$bR1i
dens(exp(m23s$bR1b + patientEffects[,1]) * 100, xlim = c(99.5,100.5), xlab = "Percent Change in Spectral Count per Day Going Forward", col = 0, ylim = c(0,200))
for(j in 1:length(patientEffects[1,])){
  percChange <- exp(m23s$bR1b + patientEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  if(probPos > -1){
    cat(paste0(as.character(patientKey[patientKey$d.patient_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
    kde <- density(percChange)
    text(x = kde$x[which.max(kde$y)], y = max(kde$y)*1.05, labels = as.character(patientKey[patientKey$d.patient_id == j, 1]))
  }
}

#plot lines for uncertainty in expected rate through time relative to rejection

#mean estimate
par(mar = c(4,4,2,2))

plot(c(-1300,-441), exp(mean(m23s$A0b) + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), 
     xlab = "days relative to rejection event", ylab = "rate of spectral count for the average protein/patient", lwd = 2)
lines(c(-441:0), exp(mean(m23s$A0b) + c(0:441) * mean(m23s$bR0b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
lines(c(0:1100), exp(mean(m23s$A1b) + c(0:1100) * mean(m23s$bR1b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
title("Expected Rates of Spectral Count Given Days Relative to Acute Rejection Event")
#89% hpdi
lines(c(-1300,-441), exp(HPDI(prob = 0.5, m23s$A0b)[[1]] + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(-441:0), exp(HPDI(prob = 0.5, m23s$A0b)[[1]] + c(0:441) * HPDI(prob = 0.5, m23s$bR0b)[[1]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(0:1100), exp(HPDI(prob = 0.5, m23s$A1b)[[1]] + c(0:1100) * HPDI(prob = 0.5, m23s$bR1b)[[1]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(-1300,-441), exp(HPDI(prob = 0.5, m23s$A0b)[[2]] + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(-441:0), exp(HPDI(prob = 0.5, m23s$A0b)[[2]] + c(0:441) * HPDI(prob = 0.5, m23s$bR0b)[[2]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(0:1100), exp(HPDI(prob = 0.5, m23s$A1b)[[2]] + c(0:1100) * HPDI(prob = 0.5, m23s$bR1b)[[2]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
#rejection date
abline(v = 0, col = 2, lwd = 2)
#transplant date
abline(v = -441, col = 2, lwd = 2)
#samples from posterior
for(i in seq(1, length(m23s$A0b), length.out = 2000)){
  lines(c(-1300,-441), exp(m23s$A0b[i] + c(0,0)), type = "l", col = rgb(0,0,0,0.03), xlim = c(-1300, 1100), ylim = c(0,10))
  lines(c(-441:0), exp(m23s$A0b[i] + c(0:441) * m23s$bR0b[i]), type = "l", col = rgb(0,0,0,0.03), xlim = c(-1300, 1100), ylim = c(0,10))
  lines(c(0:1100), exp(m23s$A1b[i] + c(0:1100) * m23s$bR1b[i]), type = "l", col = rgb(0,0,0,0.03), xlim = c(-1300, 1100), ylim = c(0,10))
}
#boxplots of data
times <- sort(unique(d$daysRelativeToRejection))
for(i in 1:length(times)){
  boxplot(d$spectral_count[d$daysRelativeToRejection == times[i]], add = T, at = times[i], boxwex=20, outline = F, col.axis = rgb(0,0,0,0))
}

#separate slopes for before and after the rejection event, seperate "intercepts" at the rejection event, intercept with zero slope before the transplant event,
#express slope in change per unit "percent" between transplant and rejection event
nChains <- 1
m24 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0*(1-afterRejection) + a1*afterRejection + bR0*daysRelativeToTransplant/daysBetwTrans_Rej*(1-afterRejection)*afterTransplant + bR1*daysRelativeToRejection*afterRejection, 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id],
    A0b ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,1),
    
    a1 <- A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,2),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,2),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,0.03),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,0.1)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=1, sigbR1i=1, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
save(file = paste0(dataSize, "/", "m24", dataSize), m24)
precis(m24, depth = 2)

#similar slopes for before and after the rejection event, seperate "intercepts" at the rejection event, intercept with zero slope before the transplant event,
#express slope in change per unit day between transplant and rejection event
nChains <- 1
m25 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0*(1-afterRejection) + a1*afterRejection + bR0*daysRelativeToTransplant*afterTransplant*(1-afterRejection) + bR0*daysRelativeToRejection*afterRejection, 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id],
    A0b ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,1),
    
    a1 <- A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id] + bR1*afterRejection,
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,0.03),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,0.1),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,0.03),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,0.1)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
save(file = paste0(dataSize, "/", "m25", dataSize), m25)
precis(m25, depth = 2)


#similar slopes for before and after the rejection event, seperate "intercepts" at the rejection event, intercept with zero slope before the transplant event,
#express slope in change per unit "percent" between transplant and rejection event
nChains <- 1
m26 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0*(1-afterRejection) + a1*afterRejection + bR0*daysRelativeToTransplant/daysBetwTrans_Rej*afterTransplant + bR1*daysRelativeToRejection*afterRejection, 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id],
    A0b ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,1),
    
    a1 <- A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,2),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,2),
    
    bR1 <- bR0/daysBetwTrans_Rej + bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,0.03),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,0.1)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=1, sigbR1i=1, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
save(file = paste0(dataSize, "/", "m26", dataSize), m26)
precis(m26, depth = 2)

#####################################################
#introducing the time relative to transplant variable -- rescaling approach
#####################################################

#separate slopes for before and after the rejection event, but same intercept at the rejection event; linear predictor on days relative to transplant
nChains <- 1
m9T <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a + bR0*((daysRelativeToRejection - daysRelativeToTransplant*(1-afterTransplant))/daysBetwTrans_Rej)*(1-afterRejection) + bR1*daysRelativeToRejection/daysBetwTrans_Rej*afterRejection, 
    
    a <- Ab + Aj[protein_id] + Ai[patient_id],
    Ab ~ dnorm(0,4),
    Aj[protein_id] ~ dnorm(0,sigAj),
    Ai[patient_id] ~ dnorm(0,sigAi),
    c(sigAj, sigAi) ~ dcauchy(0,2),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,2),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,2),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,2),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,2)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigAj=1, sigAi=1, sigbR0j=1, sigbR0i=1, bR0b=0, sigbR1j=1, sigbR1i=1, bR1b=0, Ab=0)
)
save(file = paste0(dataSize, "/", "m9T", dataSize), m9T)
precis(m9T, depth = 2)

#separate slopes for before and after the rejection event, seperate "intercepts" at the rejection event, intercept with zero slope before the transplant event,
#express slope in change per unit day between transplant and rejection event, express post-rejection intercept as displacement from pre-rejection intercept
nChains <- 1
m27 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0 + bR0*daysRelativeToTransplant*(1-afterRejection)*afterTransplant + bR1*daysRelativeToRejection*afterRejection, 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id] + a1*afterRejection,
    A0b ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,1),
    
    a1 <- bR0*daysBetwTrans_Rej + A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,0.03),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,0.01),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,0.03),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,0.01)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
save(file = paste0(dataSize, "/", "m27", dataSize), m27)
q(save = "no")
precis(m27, depth = 2)

# > compare(m23, m27, m28, m29, m30)
# WAIC  pWAIC  dWAIC weight     SE    dSE
# m27 46959.0 3397.0    0.0      1 322.14     NA
# m30 47675.4 3476.5  716.3      0 347.43  86.78
# m28 47692.1 3329.9  733.1      0 362.98 182.59
# m23 49535.2 4006.6 2576.1      0 401.56 259.33
# m29 50216.9 3959.0 3257.9      0 426.74 274.76


load("full/m27c2full")
m27 <- m27c2; rm(m27c2)
m27p <- link(m27, n = 2000)
hist(d$spectral_count, xlim = c(0,50), breaks = 200, col = rgb(0,0,1,alpha = 0.5), ylim = c(0,5000))
hist(rpois(n = length(m27p$lambda[1,]), lambda = m27p$lambda[1,]), add = T, col = rgb(1,0,0,alpha = 0.5), breaks = 200)
hist(sapply(1:length(d$spectral_count), function(x) ppois(d$spectral_count[x], m27p$lambda[,x])))
# hist(sapply(1:length(d$spectral_count), function(x) sum(rpois(2000,m27p$lambda[,x]) >= d$spectral_count[x])/2000))

#plot raw residuals from mean expectation
plot(sapply(1:length(m27p$lambda[1,]), function(x) mean(m27p$lambda[,x])), 
     sapply(1:length(m27p$lambda[1,]), function(x) mean(m27p$lambda[,x])) - d$spectral_count,
     col = rgb(0,0,0,0.3)) 

#plot pearson residuals -- mean lambda
plot(sapply(1:length(m27p$lambda[1,]), function(x) mean(m27p$lambda[,x])), 
     (sapply(1:length(m27p$lambda[1,]), function(x) mean(m27p$lambda[,x])) - d$spectral_count) / sqrt(sapply(1:length(m27p$lambda[1,]), function(x) mean(m27p$lambda[,x]))),
     xlab = "Predicted Rate", ylab = "Standardized Residual from Predicted Rate")
hist((sapply(1:length(m27p$lambda[1,]), function(x) mean(m27p$lambda[,x])) - d$spectral_count) / sqrt(sapply(1:length(m27p$lambda[1,]), function(x) mean(m27p$lambda[,x]))),
     breaks = 500, xlim = c(-5, 5), freq = F)
lines(seq(-5, 5, length.out = 200), dnorm(x = seq(-5, 5, length.out = 200), mean = 0, sd = 1))

#plot pearson residuals -- sampled lambda
numSamps <- 400
indices <- floor(seq(from = 1, to = length(m27p$lambda[,1]), length.out = numSamps))
resids <- rep(0, length(d$spectral_count)*length(indices))
for(x in 1:length(d$spectral_count)){
  if (x %% 100 == 0) {cat(paste(x, " "))}
  resids[((x-1)*numSamps+1):(x*numSamps)] <- (m27p$lambda[,x][indices] - d$spectral_count[x]) / sqrt(m27p$lambda[,x][indices])
}
hist(resids, breaks = 500, xlim = c(-5, 5), freq = F)
lines(seq(-5, 5, length.out = 200), dnorm(x = seq(-5, 5, length.out = 200), mean = 0, sd = 1))

m27s <- extract.samples(m27, clean.names = F)

#before rejection, after transplant
proteinEffects <- m27s$bR0j
proteinKey <- data.frame(as.character(d$protein), d$protein_id)
proteinKey <- proteinKey[order(proteinKey$d.protein_id),]
proteinKey <- unique(proteinKey)
hist(sapply(1:length(proteinEffects[1,]), function (x) mean(proteinEffects[,x])), breaks = 200)
dens(exp(m27s$bR0b + proteinEffects[,1]) * 100, xlim = c(98,102), xlab = "Percent Change in Spectral Count per Day", col = 0, ylim = c(0,2))
num <- 0
percChanges <- 0
proteinSpecificProbPositives <- rep(2, length(proteinEffects[1,]))
proteinSpecificMeanDailyChange <- rep(2, length(proteinEffects[1,]))
proteinList <- rep(2, length(proteinEffects[1,])) 
proteinsShared <- proteins[numPatientsPerProtein >= 3]
for(j in 1:length(proteinEffects[1,])){
  percChange <- exp(m27s$bR0b + proteinEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  proteinSpecificProbPositives[j] <- probPos
  proteinSpecificMeanDailyChange[j] <- mean(percChange)
  proteinList[j] <- as.character(proteinKey[proteinKey$d.protein_id == j, 1])
  # cat(paste0(probPos, " "))
  if(probPos > 0.95 | probPos < 0.05){
    num <- num + 1
    cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
  }
}
num
hist(proteinSpecificProbPositives, breaks = 50)
hist(proteinSpecificMeanDailyChange, breaks = 50)
numPatientsPerProtein <- sapply(1:length(proteinList), 
                                function(x) length(unique(d[d$protein == proteinList[x] & d$spectral_count != 0,]$patient)))
dataSummary <- cbind(proteinList, proteinSpecificProbPositives, proteinSpecificMeanDailyChange, numPatientsPerProtein) 
colnames(dataSummary) <- c("Protein Name", "Probability Change is > 0%", "Posterior Mean Effect", "Number of Patients w/ Non-Zero Count")
dataSummary <- dataSummary[order(dataSummary[,2], decreasing = T),]
write.csv(x = dataSummary, file = "model2_populated0s_beforeRejectionAfterTransplant_posteriorRelativeSlopeSummary.csv")

youWantToRestrictToJustPhase2 <- F
if(youWantToRestrictToJustPhase2 == T){
  numPatientsPerProtein <- sapply(1:length(proteinList), 
                                  function(x) length(unique(d[d$protein == proteinList[x] & d$spectral_count != 0 
                                                              & d$phase == 2,]$patient)))
  dataSummary <- cbind(proteinList, proteinSpecificProbPositives, proteinSpecificMeanDailyChange, numPatientsPerProtein) 
  colnames(dataSummary) <- c("Protein Name", "Probability Change is > 0%", "Posterior Mean Effect", "Number of Patients w/ Non-Zero Count")
  dataSummary <- dataSummary[order(dataSummary[,2], decreasing = T),]
  plot(dataSummary[,4], dataSummary[,2], ylab = "Probability Change is > 0%", xlab = "Number of Patients w/ Non-Zero Count in Phase 2, Specifically", xaxt = 'n')
  axis(side = 2, at = seq(0, 0.9, by = 0.1)); axis(side = 1, at = c(0,1,2,3))
}
plot(dataSummary[,2], dataSummary[,3], ylab = "Posterior Mean Effect", xlab = "Probability Change is > 0%")
plot(dataSummary[,4], dataSummary[,2], ylab = "Probability Change is > 0%", xlab = "Number of Patients w/ Non-Zero Count")
hist(exp(rnorm(n = 2e4, 0, 0.03))*100, breaks = 100, xlim = c(90,110), col = 2) #prior
hist(exp(m27s$bR0b)*100, breaks = 100, add = T, col = col2rgb("green", alpha = 0.5)) #posterior


#after rejection
proteinEffects <- m27s$bR1j
hist(sapply(1:length(proteinEffects[1,]), function (x) mean(proteinEffects[,x])), breaks = 200)
dens(exp(m27s$bR1b + proteinEffects[,1]) * 100, xlim = c(99.5,100.5), xlab = "Percent Change in Intercept with Day", col = 0, ylim = c(0,12))
num <- 0
proteinSpecificProbPositives <- rep(2, length(proteinEffects[1,]))
proteinSpecificMeanDailyChange <- rep(2, length(proteinEffects[1,]))
for(j in 1:length(proteinEffects[1,])){
  percChange <- exp(m27s$bR1b + proteinEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  proteinSpecificProbPositives[j] <- probPos
  proteinSpecificMeanDailyChange[j] <- mean(percChange)
  if(probPos > 0.95 | probPos < 0.05){
    num <- num + 1
    cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
  }
}
num
hist(proteinSpecificProbPositives, breaks = 50)
hist(proteinSpecificMeanDailyChange, breaks = 50)

dataSummary <- cbind(proteinList, proteinSpecificProbPositives, proteinSpecificMeanDailyChange, numPatientsPerProtein) 
colnames(dataSummary) <- c("Protein Name", "Probability Change is > 0%", "Posterior Mean Effect", "Number of Patients w/ Non-Zero Count")
dataSummary <- dataSummary[order(dataSummary[,2], decreasing = T),]
write.csv(x = dataSummary, file = "model2_populated0s_afterRejection_posteriorRelativeSlopeSummary.csv")

plot(dataSummary[,2], dataSummary[,3], ylab = "Posterior Mean Effect", xlab = "Probability Change is > 0%")
plot(dataSummary[,4], dataSummary[,2], ylab = "Probability Change is > 0%", xlab = "Number of Patients w/ Non-Zero Count")
hist(exp(rnorm(n = 2e4, 0, 0.03))*100, breaks = 100, xlim = c(90,110), col = 2)
hist(exp(m27s$bR1b)*100, breaks = 100, add = T, col = col2rgb("green", alpha = 0.5))


#discontinuity at rejection
proteinEffects <- m27s$A1j
hist(sapply(1:length(proteinEffects[1,]), function (x) mean(proteinEffects[,x])), breaks = 200)
dens(exp(m27s$A1b + proteinEffects[,1]) * 100, xlim = c(0,300), xlab = "Percent Change in Intercept at Rejection", col = 0, ylim = c(0,0.05))
num <- 0
proteinSpecificProbPositives <- rep(2, length(proteinEffects[1,]))
proteinSpecificMeanChange <- rep(2, length(proteinEffects[1,]))
for(j in 1:length(proteinEffects[1,])){
  percChange <- exp(m27s$A1b + proteinEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  proteinSpecificProbPositives[j] <- probPos
  proteinSpecificMeanChange[j] <- mean(percChange)
  # cat(paste(probPos, "\t"))
  if(probPos > 0.95 | probPos < 0.05){
    num <- num + 1
    cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
  }
}
num
hist(proteinSpecificProbPositives, breaks = 50)
hist(proteinSpecificMeanChange, breaks = 50)
dataSummary <- cbind(proteinList, proteinSpecificProbPositives, proteinSpecificMeanChange, numPatientsPerProtein) 
colnames(dataSummary) <- c("Protein Name", "Probability Change is > 0%", "Posterior Mean Effect", "Number of Patients w/ Non-Zero Count")
dataSummary <- dataSummary[order(as.numeric(dataSummary[,2]), decreasing = F),]
write.csv(x = dataSummary, file = "model2_populated0s_displacementAtRejection_posteriorRelativeChangeSummary.csv")

plot(dataSummary[,2], dataSummary[,3], ylab = "Posterior Mean Effect", xlab = "Probability Change is > 0%")
plot(dataSummary[,4], dataSummary[,2], ylab = "Probability Change is > 0%", xlab = "Number of Patients w/ Non-Zero Count")

# probing specific proteins
targetProteins <- c("troponin", "prohibitin", "plakoglobin")
targetProteins <- proteinKey[unlist(sapply(1:length(targetProteins), function(x) which(grepl(targetProteins[x], proteinKey$as.character.d.protein)))),]
proteinEffects <- m27s$bR0j
for(j in 1:length(targetProteins[,1])){
  percChange <- exp(m27s$bR0b + proteinEffects[,targetProteins[j,2]]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == targetProteins[j,2], 1]), ", probability positive = ", probPos * 100, "%, mean change = ", round(mean(percChange), 3), "%\n"))
}
proteinEffects <- m27s$bR1j
for(j in 1:length(targetProteins[,1])){
  percChange <- exp(m27s$bR1b + proteinEffects[,targetProteins[j,2]]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == targetProteins[j,2], 1]), ", probability positive = ", probPos * 100, "%, mean change = ", round(mean(percChange), 3), "%\n"))
}
proteinEffects <- m27s$A1j
for(j in 1:length(targetProteins[,1])){
  percChange <- exp(m27s$A1b + proteinEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == targetProteins[j,2], 1]), ", probability positive = ", probPos * 100, "%, mean change = ", round(mean(percChange), 3), "%\n"))
}

#individual effects
patientEffects <- m27s$bR0i
patientKey <- data.frame(as.character(d$patient), d$patient_id)
patientKey <- patientKey[order(patientKey$d.patient_id),]
patientKey <- unique(patientKey)
dens(exp(m27s$bR0b + patientEffects[,1]) * 100, xlim = c(99,101), xlab = "Percent Change in Spectral Count per Day Going Backward", col = 0, ylim = c(0,50))
for(j in 1:length(patientEffects[1,])){
  percChange <- exp(m27s$bR0b + patientEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  if(probPos > -1){
    cat(paste0(as.character(patientKey[patientKey$d.patient_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
    kde <- density(percChange)
    text(x = kde$x[which.max(kde$y)], y = max(kde$y)*1.05, labels = as.character(patientKey[patientKey$d.patient_id == j, 1]))
  }
}

patientEffects <- m27s$bR1i
dens(exp(m27s$bR1b + patientEffects[,1]) * 100, xlim = c(99.5,100.5), xlab = "Percent Change in Spectral Count per Day Going Forward", col = 0, ylim = c(0,200))
for(j in 1:length(patientEffects[1,])){
  percChange <- exp(m27s$bR1b + patientEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  if(probPos > -1){
    cat(paste0(as.character(patientKey[patientKey$d.patient_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
    kde <- density(percChange)
    text(x = kde$x[which.max(kde$y)], y = max(kde$y)*1.05, labels = as.character(patientKey[patientKey$d.patient_id == j, 1]))
  }
}

patientEffects <- m27s$A1i
hist(sapply(1:length(patientEffects[1,]), function (x) mean(patientEffects[,x])), breaks = 200)
dens(exp(m27s$A1b + patientEffects[,1]) * 100, xlim = c(0,300), xlab = "Percent Change in Intercept at Rejection", col = 0, ylim = c(0,0.08))
patientSpecificProbPositives <- rep(2, length(patientEffects[1,]))
patientSpecificMeanDailyChange <- rep(2, length(patientEffects[1,]))
for(j in 1:length(patientEffects[1,])){
  percChange <- exp(m27s$A1b + patientEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  patientSpecificProbPositives[j] <- probPos
  patientSpecificMeanDailyChange[j] <- mean(percChange)
  # cat(paste(probPos, "\t"))
  cat(paste0(as.character(patientKey[patientKey$d.patient_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
  dens(percChange, add = T)
  kde <- density(percChange)
  text(x = kde$x[which.max(kde$y)], y = max(kde$y)*1.05, labels = as.character(patientKey[patientKey$d.patient_id == j, 1]))
}


#plot lines for uncertainty in expected rate through time relative to rejection

#mean estimate
par(mar = c(4,4,2,2))
plot(c(-1300,-441), exp(mean(m27s$A0b) + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,3), 
     xlab = "days relative to rejection event", ylab = "rate of spectral count for the average protein/patient", lwd = 2)
lines(c(-441:0), exp(mean(m27s$A0b) + c(0:441) * mean(m27s$bR0b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
lines(c(0:1100), exp(mean(m27s$A1b) + mean(m27s$A0b) + 441*mean(m27s$bR0b) + c(0:1100) * mean(m27s$bR1b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
title("Expected Rates of Spectral Count Given Days Relative to Acute Rejection Event")
#50% hpdi
lines(c(-1300,-441), exp(HPDI(prob = 0.5, m27s$A0b)[[1]] + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(-441:0), exp(HPDI(prob = 0.5, m27s$A0b)[[1]] + c(0:441) * HPDI(prob = 0.5, m27s$bR0b)[[1]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(0:1100), exp(HPDI(prob = 0.5, m27s$A0b)[[1]] + 441*HPDI(prob = 0.5, m27s$bR0b)[[1]] + HPDI(prob = 0.5, m27s$A1b)[[1]] + c(0:1100) * HPDI(prob = 0.5, m27s$bR1b)[[1]]), 
      type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(-1300,-441), exp(HPDI(prob = 0.5, m27s$A0b)[[2]] + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(-441:0), exp(HPDI(prob = 0.5, m27s$A0b)[[2]] + c(0:441) * HPDI(prob = 0.5, m27s$bR0b)[[2]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(0:1100), exp(HPDI(prob = 0.5, m27s$A0b)[[2]] + 441*HPDI(prob = 0.5, m27s$bR0b)[[2]] + HPDI(prob = 0.5, m27s$A1b)[[2]] + c(0:1100) * HPDI(prob = 0.5, m27s$bR1b)[[2]]), 
      type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
#rejection date
abline(v = 0, col = 2, lwd = 2)
text(x = -25, y = 2, labels = "Rejection Date", srt = 90)
#transplant date
abline(v = -441, col = 2, lwd = 2)
text(x = -466, y = 2, labels = "Average Transplant Date", srt = 90)
#samples from posterior
for(i in seq(1, length(m27s$A0b), length.out = 2000)){
  lines(c(-1300,-441), exp(m27s$A0b[i] + c(0,0)), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
  lines(c(-441:0), exp(m27s$A0b[i] + c(0:441) * m27s$bR0b[i]), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
  lines(c(0:1100), exp(m27s$A0b[i] + 441*m27s$bR0b[i] + m27s$A1b[i] + c(0:1100) * m27s$bR1b[i]), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
}
#boxplots of data
times <- sort(unique(d$daysRelativeToRejection))
for(i in 1:length(times)){
  boxplot(d$spectral_count[d$daysRelativeToRejection == times[i]], add = T, at = times[i], boxwex=20, outline = F, col.axis = rgb(0,0,0,0))
}

#explore power law relationships for all slope coefficients, partial pooling on phases 2 and 3
nChains <- 1
d$daysRelativeToTransplantPos <- d$daysRelativeToTransplant 
d$daysRelativeToTransplantPos[d$daysRelativeToTransplantPos < 0] <- 1e4
d$daysRelativeToRejectionPos <- d$daysRelativeToRejection 
d$daysRelativeToRejectionPos[d$daysRelativeToRejectionPos < 0] <- 1e4
m28 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0 + bR0*log(daysRelativeToTransplantPos)*(1-afterRejection)*afterTransplant + bR1*log(daysRelativeToRejectionPos)*afterRejection, 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id] + a1*afterRejection,
    A0b ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,1),
    
    a1 <- bR0*log(daysBetwTrans_Rej) + A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,2),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,1),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,2),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,1)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigbR0j=1, sigbR0i=1, bR0b=0, sigbR1j=1, sigbR1i=1, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
save(file = paste0(dataSize, "/", "m28", dataSize), m28)
precis(m28, depth = 2)


load("full/m28full")
m28p <- link(m28, n = 2000)
hist(d$spectral_count, xlim = c(0,50), breaks = 200, col = rgb(0,0,1,alpha = 0.5), ylim = c(0,30000))
hist(rpois(n = length(m28p$lambda[1,]), lambda = m28p$lambda[1,]), add = T, col = rgb(1,0,0,alpha = 0.5), breaks = 200)

hist(d$spectral_count, xlim = c(0,18), breaks = 200, col = rgb(0,0,1,alpha = 0.5), ylim = c(0,2500))
hist(rpois(n = length(m28p$lambda[1,]), lambda = m28p$lambda[1,]), add = T, col = rgb(1,0,0,alpha = 0.5), breaks = 200)

hist(sapply(1:length(d$spectral_count), function(x) ppois(d$spectral_count[x], m28p$lambda[,x])))
# hist(sapply(1:length(d$spectral_count), function(x) sum(rpois(2000,m28p$lambda[,x]) >= d$spectral_count[x])/2000))

#plot raw residuals from mean expectation
plot(sapply(1:length(m28p$lambda[1,]), function(x) mean(m28p$lambda[,x])), 
     sapply(1:length(m28p$lambda[1,]), function(x) mean(m28p$lambda[,x])) - d$spectral_count,
     col = rgb(0,0,0,0.3)) 

#plot pearson residuals -- mean lambda
plot(sapply(1:length(m28p$lambda[1,]), function(x) mean(m28p$lambda[,x])), 
     (sapply(1:length(m28p$lambda[1,]), function(x) mean(m28p$lambda[,x])) - d$spectral_count) / sqrt(sapply(1:length(m28p$lambda[1,]), function(x) mean(m28p$lambda[,x]))),
     xlab = "Predicted Rate", ylab = "Standardized Residual from Predicted Rate")
hist((sapply(1:length(m28p$lambda[1,]), function(x) mean(m28p$lambda[,x])) - d$spectral_count) / sqrt(sapply(1:length(m28p$lambda[1,]), function(x) mean(m28p$lambda[,x]))),
     breaks = 50, xlim = c(-5, 5), freq = F)
lines(seq(-5, 5, length.out = 200), dnorm(x = seq(-5, 5, length.out = 200), mean = 0, sd = 1))

#plot pearson residuals -- sampled lambda
numSamps <- 400
indices <- floor(seq(from = 1, to = length(m28p$lambda[,1]), length.out = numSamps))
resids <- rep(0, length(d$spectral_count)*length(indices))
for(x in 1:length(d$spectral_count)){
  if (x %% 100 == 0) {cat(paste(x, " "))}
  resids[((x-1)*numSamps+1):(x*numSamps)] <- (m28p$lambda[,x][indices] - d$spectral_count[x]) / sqrt(m28p$lambda[,x][indices])
}
hist(resids, breaks = 500, xlim = c(-5, 5), freq = F)
lines(seq(-5, 5, length.out = 200), dnorm(x = seq(-5, 5, length.out = 200), mean = 0, sd = 1))

m28s <- extract.samples(m28, clean.names = F)

#before rejection, after transplant
proteinEffects <- m28s$bR0j
proteinKey <- data.frame(as.character(d$protein), d$protein_id)
proteinKey <- proteinKey[order(proteinKey$d.protein_id),]
proteinKey <- unique(proteinKey)
hist(sapply(1:length(proteinEffects[1,]), function (x) mean(proteinEffects[,x])), breaks = 200)
dens(exp(m28s$bR0b + proteinEffects[,1]) * 100, xlim = c(50,400), xlab = "Percent Change in Spectral Count per ln(Day)", col = 0, ylim = c(0,0.025))
num <- 0
percChanges <- 0
proteinSpecificProbPositives <- rep(2, length(proteinEffects[1,]))
proteinSpecificMeanLogDailyChange <- rep(2, length(proteinEffects[1,]))
proteinList <- rep(2, length(proteinEffects[1,]))
proteinsShared <- proteins[numPatientsPerProtein >= 3]
for(j in 1:length(proteinEffects[1,])){
  percChange <- exp(m28s$bR0b + proteinEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  proteinSpecificProbPositives[j] <- probPos
  proteinSpecificMeanLogDailyChange[j] <- mean(percChange)
  proteinList[j] <- as.character(proteinKey[proteinKey$d.protein_id == j, 1])
  # cat(paste0(probPos, " "))
  if(probPos > 0.99 | probPos < 0.01){
    num <- num + 1
    cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
  }
}
num
hist(proteinSpecificProbPositives, breaks = 50)
hist(proteinSpecificMeanLogDailyChange, breaks = 50)
numPatientsPerProtein <- sapply(1:length(proteinList), 
                                function(x) length(unique(d[d$protein == proteinList[x] & d$spectral_count != 0,]$patient)))
dataSummary <- cbind(proteinList, proteinSpecificProbPositives, proteinSpecificMeanLogDailyChange, numPatientsPerProtein) 
colnames(dataSummary) <- c("Protein Name", "Probability Change is > 0%", "Posterior Mean Effect", "Number of Patients w/ Non-Zero Count")
dataSummary <- dataSummary[order(dataSummary[,2], decreasing = T),]
write.csv(x = dataSummary, file = "model3_populated0s_beforeRejectionAfterTransplant_posteriorRelativeSlopeSummary.csv")

youWantToRestrictToJustPhase2 <- T
if(youWantToRestrictToJustPhase2 == T){
  numPatientsPerProtein <- sapply(1:length(proteinList), 
                                  function(x) length(unique(d[d$protein == proteinList[x] & d$spectral_count != 0 
                                                              & d$phase == 2,]$patient)))
  dataSummary <- cbind(proteinList, proteinSpecificProbPositives, proteinSpecificMeanLogDailyChange, numPatientsPerProtein) 
  colnames(dataSummary) <- c("Protein Name", "Probability Change is > 0%", "Posterior Mean Effect", "Number of Patients w/ Non-Zero Count")
  dataSummary <- dataSummary[order(dataSummary[,2], decreasing = T),]
  plot(dataSummary[,4], dataSummary[,2], ylab = "Probability Change is > 0%", xlab = "Number of Patients w/ Non-Zero Count in Phase 2, Specifically", xaxt = 'n')
  axis(side = 2, at = seq(0, 0.9, by = 0.1)); axis(side = 1, at = c(0,1,2,3))
}
plot(dataSummary[,2], dataSummary[,3], ylab = "Posterior Mean Effect", xlab = "Probability Change is > 0%")
plot(dataSummary[,4], dataSummary[,2], ylab = "Probability Change is > 0%", xlab = "Number of Patients w/ Non-Zero Count")
hist(exp(rnorm(n = 1e4, 0, 2))*100, breaks = 100000, xlim = c(0,400), col = 2)
hist(exp(m28s$bR0b)*100, breaks = 1000, add = T, col = col2rgb("green", alpha = 0.5))

#after rejection
proteinEffects <- m28s$bR1j
hist(sapply(1:length(proteinEffects[1,]), function (x) mean(proteinEffects[,x])), breaks = 200)
dens(exp(m28s$bR1b + proteinEffects[,1]) * 100, xlim = c(50,400), xlab = "Percent Change in Intercept with ln(Day)", col = 0, ylim = c(0,0.025))
num <- 0
proteinSpecificProbPositives <- rep(2, length(proteinEffects[1,]))
proteinSpecificMeanLogDailyChange <- rep(2, length(proteinEffects[1,]))
proteinList <- rep(2, length(proteinEffects[1,]))
for(j in 1:length(proteinEffects[1,])){
  percChange <- exp(m28s$bR1b + proteinEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  proteinSpecificProbPositives[j] <- probPos
  proteinSpecificMeanLogDailyChange[j] <- mean(percChange)
  proteinList[j] <- as.character(proteinKey[proteinKey$d.protein_id == j, 1])
  if(probPos > 0.95 | probPos < 0.05){
    num <- num + 1
    cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
  }
}
num
hist(proteinSpecificProbPositives, breaks = 50)
hist(proteinSpecificMeanLogDailyChange, breaks = 50)
numPatientsPerProtein <- sapply(1:length(proteinList), 
                                function(x) length(unique(d[d$protein == proteinList[x] & d$spectral_count != 0,]$patient)))
dataSummary <- cbind(proteinList, proteinSpecificProbPositives, proteinSpecificMeanLogDailyChange, numPatientsPerProtein) 
colnames(dataSummary) <- c("Protein Name", "Probability Change is > 0%", "Posterior Mean Effect", "Number of Patients w/ Non-Zero Count")
dataSummary <- dataSummary[order(dataSummary[,2], decreasing = T),]
write.csv(x = dataSummary, file = "model3_populated0s_afterRejection_posteriorRelativeSlopeSummary.csv")

plot(dataSummary[,2], dataSummary[,3], ylab = "Posterior Mean Effect", xlab = "Probability Change is > 0%")
plot(dataSummary[,4], dataSummary[,2], ylab = "Probability Change is > 0%", xlab = "Number of Patients w/ Non-Zero Count")
hist(exp(rnorm(n = 1e4, 0, 2))*100, breaks = 100000, xlim = c(0,400), col = 2)
hist(exp(m28s$bR1b)*100, breaks = 1000, add = T, col = col2rgb("green", alpha = 0.5))


#discontinuity at rejection
proteinEffects <- m28s$A1j
hist(sapply(1:length(proteinEffects[1,]), function (x) mean(proteinEffects[,x])), breaks = 200)
dens(exp(m28s$A1b + proteinEffects[,1]) * 100, xlim = c(0,200), xlab = "Percent Change in Intercept at Rejection", col = 0, ylim = c(0,0.4))
num <- 0
proteinSpecificProbPositives <- rep(2, length(proteinEffects[1,]))
proteinSpecificMeanChange <- rep(2, length(proteinEffects[1,]))
proteinList <- rep(2, length(proteinEffects[1,]))
for(j in 1:length(proteinEffects[1,])){
  percChange <- exp(m28s$A1b + proteinEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  proteinSpecificProbPositives[j] <- probPos
  proteinSpecificMeanChange[j] <- mean(percChange)
  proteinList[j] <- as.character(proteinKey[proteinKey$d.protein_id == j, 1])
  # cat(paste(probPos, "\t"))
  if(probPos > 0.99 | probPos < 0.01){
    num <- num + 1
    cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
  }
}
num
hist(proteinSpecificProbPositives, breaks = 50)
hist(proteinSpecificMeanChange, breaks = 50)
numPatientsPerProtein <- sapply(1:length(proteinList), 
                                function(x) length(unique(d[d$protein == proteinList[x] & d$spectral_count != 0,]$patient)))
dataSummary <- cbind(proteinList, proteinSpecificProbPositives, proteinSpecificMeanChange, numPatientsPerProtein) 
colnames(dataSummary) <- c("Protein Name", "Probability Change is > 0%", "Posterior Mean Effect", "Number of Patients w/ Non-Zero Count")
dataSummary <- dataSummary[order(dataSummary[,2], decreasing = T),]
write.csv(x = dataSummary, file = "model3_populated0s_displacementAtRejection_posteriorRelativeSlopeSummary.csv")

plot(dataSummary[,2], dataSummary[,3], ylab = "Posterior Mean Effect", xlab = "Probability Change is > 0%")
plot(dataSummary[,4], dataSummary[,2], ylab = "Probability Change is > 0%", xlab = "Number of Patients w/ Non-Zero Count")


# probing specific proteins
targetProteins <- c("troponin", "prohibitin", "plakoglobin")
targetProteins <- proteinKey[unlist(sapply(1:length(targetProteins), function(x) which(grepl(targetProteins[x], proteinKey$as.character.d.protein)))),]
proteinEffects <- m28s$bR0j
for(j in 1:length(targetProteins[,1])){
  percChange <- exp(m28s$bR0b + proteinEffects[,targetProteins[j,2]]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == targetProteins[j,2], 1]), ", probability positive = ", probPos * 100, "%, mean change = ", round(mean(percChange), 3), "%\n"))
}
proteinEffects <- m28s$bR1j
for(j in 1:length(targetProteins[,1])){
  percChange <- exp(m28s$bR1b + proteinEffects[,targetProteins[j,2]]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == targetProteins[j,2], 1]), ", probability positive = ", probPos * 100, "%, mean change = ", round(mean(percChange), 3), "%\n"))
}
proteinEffects <- m28s$A1j
for(j in 1:length(targetProteins[,1])){
  percChange <- exp(m28s$A1b + proteinEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == targetProteins[j,2], 1]), ", probability positive = ", probPos * 100, "%, mean change = ", round(mean(percChange), 3), "%\n"))
}

#individual effects
patientEffects <- m28s$bR0i
patientKey <- data.frame(as.character(d$patient), d$patient_id)
patientKey <- patientKey[order(patientKey$d.patient_id),]
patientKey <- unique(patientKey)
dens(exp(m28s$bR0b + patientEffects[,1]) * 100, xlim = c(50,400), xlab = "Percent Change in Spectral Count per ln(Day) Going Backward", col = 0, ylim = c(0,0.15))
for(j in 1:length(patientEffects[1,])){
  percChange <- exp(m28s$bR0b + patientEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  if(probPos > -1){
    cat(paste0(as.character(patientKey[patientKey$d.patient_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
    kde <- density(percChange)
    text(x = kde$x[which.max(kde$y)], y = max(kde$y)*1.05, labels = as.character(patientKey[patientKey$d.patient_id == j, 1]))
  }
}

patientEffects <- m28s$bR1i
dens(exp(m28s$bR1b + patientEffects[,1]) * 100, xlim = c(50,200), xlab = "Percent Change in Spectral Count per ln(Day) Going Forward", col = 0, ylim = c(0,0.5))
for(j in 1:length(patientEffects[1,])){
  percChange <- exp(m28s$bR1b + patientEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  if(probPos > -1){
    cat(paste0(as.character(patientKey[patientKey$d.patient_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
    kde <- density(percChange)
    text(x = kde$x[which.max(kde$y)], y = max(kde$y)*1.05, labels = as.character(patientKey[patientKey$d.patient_id == j, 1]))
  }
}

patientEffects <- m28s$A1i
hist(sapply(1:length(patientEffects[1,]), function (x) mean(patientEffects[,x])), breaks = 200)
dens(exp(m28s$A1b + patientEffects[,1]) * 100, xlim = c(0,300), xlab = "Percent Change in Intercept at Rejection", col = 0, ylim = c(0,0.6))
patientSpecificProbPositives <- rep(2, length(patientEffects[1,]))
patientSpecificMeanDailyChange <- rep(2, length(patientEffects[1,]))
for(j in 1:length(patientEffects[1,])){
  percChange <- exp(m28s$A1b + patientEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  patientSpecificProbPositives[j] <- probPos
  patientSpecificMeanDailyChange[j] <- mean(percChange)
  # cat(paste(probPos, "\t"))
  cat(paste0(as.character(patientKey[patientKey$d.patient_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
  dens(percChange, add = T)
  kde <- density(percChange)
  text(x = kde$x[which.max(kde$y)], y = max(kde$y)*1.05, labels = as.character(patientKey[patientKey$d.patient_id == j, 1]))
}


#plot lines for uncertainty in expected rate through time relative to rejection

#mean estimate
par(mar = c(4,4,2,2))
plot(c(-1300,-441), exp(mean(m28s$A0b) + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,3), 
     xlab = "days relative to rejection event", ylab = "rate of spectral count for the average protein/patient", lwd = 2)
lines(c(-440:0), exp(mean(m28s$A0b) + log(c(1:441)) * mean(m28s$bR0b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
lines(c(1:1100), exp(mean(m28s$A1b) + mean(m28s$A0b) + log(441)*mean(m28s$bR0b) + log(c(1:1100)) * mean(m28s$bR1b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
title("Expected Rates of Spectral Count Given Days Relative to Acute Rejection Event")
#50% hpdi
lines(c(-1300,-441), exp(HPDI(prob = 0.5, m28s$A0b)[[1]] + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(-441:0), exp(HPDI(prob = 0.5, m28s$A0b)[[1]] + log(c(0:441)) * HPDI(prob = 0.5, m28s$bR0b)[[1]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(0:1100), exp(HPDI(prob = 0.5, m28s$A0b)[[1]] + log(441)*HPDI(prob = 0.5, m28s$bR0b)[[1]] + HPDI(prob = 0.5, m28s$A1b)[[1]] + log(c(0:1100)) * HPDI(prob = 0.5, m28s$bR1b)[[1]]), 
      type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(-1300,-441), exp(HPDI(prob = 0.5, m28s$A0b)[[2]] + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(-440:0), exp(HPDI(prob = 0.5, m28s$A0b)[[2]] + log(c(1:441)) * HPDI(prob = 0.5, m28s$bR0b)[[2]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(1:1100), exp(HPDI(prob = 0.5, m28s$A0b)[[2]] + log(441)*HPDI(prob = 0.5, m28s$bR0b)[[2]] + HPDI(prob = 0.5, m28s$A1b)[[2]] + log(c(1:1100)) * HPDI(prob = 0.5, m28s$bR1b)[[2]]), 
      type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
#rejection date
abline(v = 0, col = 2, lwd = 2)
text(x = -25, y = 2, labels = "Rejection Date", srt = 90)
#transplant date
abline(v = -441, col = 2, lwd = 2)
text(x = -466, y = 2, labels = "Average Transplant Date", srt = 90)
#samples from posterior
for(i in seq(1, length(m28s$A0b), length.out = 2000)){
  lines(c(-1300,-441), exp(m28s$A0b[i] + c(0,0)), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
  lines(c(-440:0), exp(m28s$A0b[i] + log(c(1:441)) * m28s$bR0b[i]), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
  lines(c(1:1100), exp(m28s$A0b[i] + log(441)*m28s$bR0b[i] + m28s$A1b[i] + log(c(1:1100)) * m28s$bR1b[i]), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
}

#boxplots of data
times <- sort(unique(d$daysRelativeToRejection))
for(i in 1:length(times)){
  boxplot(d$spectral_count[d$daysRelativeToRejection == times[i]], add = T, at = times[i], boxwex=20, outline = F, col.axis = rgb(0,0,0,0))
}

#separate slopes for before and after the rejection event, seperate "intercepts" at the rejection event, intercept with zero slope before the transplant event,
#express slope in change per unit day between transplant and rejection event
nChains <- 1
d$daysRelativeToTransplantPos <- d$daysRelativeToTransplant 
d$daysRelativeToTransplantPos[d$daysRelativeToTransplantPos < 0] <- 1e4
d$daysRelativeToRejectionPos <- d$daysRelativeToRejection 
d$daysRelativeToRejectionPos[d$daysRelativeToRejectionPos < 0] <- 1e4
m29 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0*(1-afterRejection) + a1*afterRejection + bR0*log(daysRelativeToTransplantPos)*(1-afterRejection)*afterTransplant + bR1*log(daysRelativeToRejectionPos)*afterRejection, 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id],
    A0b ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,1),
    
    a1 <- A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,2),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,1),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,2),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,1)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigbR0j=1, sigbR0i=1, bR0b=0, sigbR1j=1, sigbR1i=1, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
save(file = paste0(dataSize, "/", "m29", dataSize), m29)
precis(m29, depth = 2)
plot(m29)


#explore power law relationships for post rejection slope coefficients, exponential for 3, partial pooling for intercepts on phases 2 and 3
nChains <- 1
d$daysRelativeToTransplantPos <- d$daysRelativeToTransplant 
d$daysRelativeToTransplantPos[d$daysRelativeToTransplantPos < 0] <- 1e4
d$daysRelativeToRejectionPos <- d$daysRelativeToRejection 
d$daysRelativeToRejectionPos[d$daysRelativeToRejectionPos < 0] <- 1e4
m30 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0 + bR0*daysRelativeToTransplantPos*(1-afterRejection)*afterTransplant + bR1*log(daysRelativeToRejectionPos)*afterRejection, 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id] + a1*afterRejection,
    A0b ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,1),
    
    a1 <- bR0*daysBetwTrans_Rej + A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,0.03),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,0.01),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,2),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,1)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigbR0j=.01, sigbR0i=.01, bR0b=0, sigbR1j=1, sigbR1i=1, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
save(file = paste0(dataSize, "/", "m30", dataSize), m30)
precis(m30, depth = 2)

#explore stairstep intercepts for all three phases, pooling across individuals, adjascent phases, and proteins
nChains <- 1
m31 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a, 
    
    a <- A0b + A0j[protein_id] + A0i[patient_id] + a1*afterTransplant + a2*afterRejection,
    A0b ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,1),
    
    a1 <- A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,1),
    
    a2 <- A2b + A2j[protein_id] + A2i[patient_id],
    A2b ~ dnorm(0,4),
    A2j[protein_id] ~ dnorm(0,sigA2j),
    A2i[patient_id] ~ dnorm(0,sigA2i),
    c(sigA2j, sigA2i) ~ dcauchy(0,1)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1, A2b=0, sigA2j=1, sigA2i=1)
)
save(file = paste0(dataSize, "/", "m31", dataSize), m31)
precis(m31, depth = 2)

#separate slopes for before and after the rejection event, seperate intercepts at the rejection event, 
#data before the transplant informs starting point at the transplant to the rejection; alternative parameterization that is not as principled
#but easier
nChains <- 1
m13T <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0*(1-afterRejection) + a1*afterRejection + 
      bR0*((daysRelativeToRejection - daysRelativeToTransplant*(1-afterTransplant))/daysBetwTrans_Rej)*(1-afterRejection) + 
      bR1*daysRelativeToRejection/daysBetwTrans_Rej*afterRejection, 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id],
    A0b ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,2),
    
    a1 <- A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,2),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,2),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,2),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,2),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,2)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigbR0j=1, sigbR0i=1, bR0b=0, sigbR1j=1, sigbR1i=1, bR1b=0, A1b=0, sigA1j=1, sigA1i=1)
)
save(file = paste0(dataSize, "/", "m13T", dataSize), m13T)
precis(m13T, depth = 2)


#separate slopes for before and after the rejection event, seperate "intercepts" at the rejection event, intercept with zero slope before the transplant event,
#express slope in change per unit day between transplant and rejection event, express post-rejection intercept as displacement from pre-rejection intercept
nChains <- 1
m32 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0 + bR0*daysRelativeToTransplant*(1-afterRejection)*afterTransplant + bR1*daysRelativeToRejection*afterRejection, 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id] + a1*afterRejection + a2*afterTransplant,
    A0b ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,1),
    
    a1 <- bR0*daysBetwTrans_Rej + A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,1),
    
    a2 <- A2b + A2j[protein_id] + A2i[patient_id],
    A2b ~ dnorm(0,4),
    A2j[protein_id] ~ dnorm(0,sigA2j),
    A2i[patient_id] ~ dnorm(0,sigA2i),
    c(sigA2j, sigA2i) ~ dcauchy(0,1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,0.03),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,0.01),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,0.03),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,0.01)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1, A2b=0, sigA2j=1, sigA2i=1)
)
save(file = paste0(dataSize, "/", "m32", dataSize), m32)
precis(m32, depth = 2)


#intercept before the rejection event, slope and interecept after the rejection event
nChains <- 1
m33 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a + bR*daysRelativeToRejection*afterRejection, 
    
    a <- Ab + Aj[protein_id] + Ai[patient_id],
    Ab ~ dnorm(0,4),
    Aj[protein_id] ~ dnorm(0,sigAj),
    Ai[patient_id] ~ dnorm(0,sigAi),
    c(sigAj, sigAi) ~ dcauchy(0,1),
    
    bR <- bRb + bRj[protein_id] + bRi[patient_id],
    bRj[protein_id] ~ dnorm(0,sigbRj),
    bRi[patient_id] ~ dnorm(0,sigbRi),
    bRb ~ dnorm(0,0.03),
    c(sigbRj, sigbRi) ~ dcauchy(0,0.01)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigbRj=0.01, sigbRi=0.01, bRb=0, Ab=0, sigAj=1, sigAi=1)
)
save(file = paste0(dataSize, "/", "m33", dataSize), m33)
precis(m33, depth = 2)

#model 27, but with tighter priors
nChains <- 1
m34 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0 + bR0*daysRelativeToTransplant*(1-afterRejection)*afterTransplant + bR1*daysRelativeToRejection*afterRejection, 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id] + a1*afterRejection,
    A0b ~ dnorm(0,1),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,1),
    
    a1 <- bR0*daysBetwTrans_Rej + A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,1),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,0.03),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,0.01),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,0.03),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,0.01)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
save(file = paste0(dataSize, "/", "m34", dataSize), m34)
precis(m34, depth = 2)


#separate slopes for before and after the rejection event, seperate "intercepts" at the rejection event, intercept with zero slope before the transplant event,
#express slope in change per unit day between transplant and rejection event, express post-rejection intercept as displacement from pre-rejection intercept
#allow for exponential effect on relative change per day, learn exponent from data
nChains <- 1
d$daysRelativeToTransplantPos <- d$daysRelativeToTransplant 
d$daysRelativeToTransplantPos[d$daysRelativeToTransplantPos < 0] <- 1e4
d$daysRelativeToRejectionPos <- d$daysRelativeToRejection 
d$daysRelativeToRejectionPos[d$daysRelativeToRejectionPos < 0] <- 1e4
sampleIterations <- 30000; warmupIterations <- 30000; nChains <- 1
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
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains, 
  start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1, tunerExp2= 0, tunerExp3=0)
)
save(file = paste0(dataSize, "/", "m35loosepriors", dataSize), m35)
q(save = "no")
precis(m35, depth = 2)

#explore power law relationships for all slope coefficients, complete independence between phases
nChains <- 1
d$daysRelativeToTransplantPos <- d$daysRelativeToTransplant 
d$daysRelativeToTransplantPos[d$daysRelativeToTransplantPos < 0] <- 1e4
d$daysRelativeToRejectionPos <- d$daysRelativeToRejection 
d$daysRelativeToRejectionPos[d$daysRelativeToRejectionPos < 0] <- 1e4
m36 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0*(1-afterTransplant) + bR0*log(daysRelativeToTransplantPos)*(1-afterRejection)*afterTransplant + bR1*log(daysRelativeToRejectionPos)*afterRejection, 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id],
    A0b ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,2),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,1),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,2),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,1)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigbR0j=1, sigbR0i=1, bR0b=0, sigbR1j=1, sigbR1i=1, bR1b=0, A0b=0, sigA0j=1, sigA0i=1)
)
save(file = paste0(dataSize, "/", "m36", dataSize), m36)
q(save = "no")
precis(m36, depth = 2)

#separate slopes for before and after the rejection event, seperate "intercepts" at the rejection event, intercept with zero slope before the transplant event,
#express slope in change per unit day between transplant and rejection event, express post-rejection intercept as displacement from pre-rejection intercept
#allow for exponential effect on relative change per day, learn exponent from data, but constraint it to a more reasonable, informative, regularizing range
nChains <- 1
d$daysRelativeToTransplantPos <- d$daysRelativeToTransplant 
d$daysRelativeToTransplantPos[d$daysRelativeToTransplantPos < 0] <- 1e4
d$daysRelativeToRejectionPos <- d$daysRelativeToRejection 
d$daysRelativeToRejectionPos[d$daysRelativeToRejectionPos < 0] <- 1e4
m37 <- map2stan(
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
    bR0b ~ dnorm(0,0.5),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,1),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,2),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,1),
    
    p2tuner ~ dgamma(10,20),
    p3tuner ~ dgamma(30,20)

  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains, 
  start = list(sigbR0j=1, sigbR0i=1, bR0b=0, sigbR1j=1, sigbR1i=1, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1, p2tuner= 0.5, p3tuner=1.5)
)
save(file = paste0(dataSize, "/", "m37", dataSize), m37)
q(save = "no")
precis(m37, depth = 2)

#TRYING TO EXAMINE THE PRIOR
#separate slopes for before and after the rejection event, seperate "intercepts" at the rejection event, intercept with zero slope before the transplant event,
#express slope in change per unit day between transplant and rejection event, express post-rejection intercept as displacement from pre-rejection intercept
#allow for exponential effect on relative change per day, learn exponent from data, but constraint it to a more reasonable, informative, regularizing range
nChains <- 1
d$daysRelativeToTransplantPos <- d$daysRelativeToTransplant 
d$daysRelativeToTransplantPos[d$daysRelativeToTransplantPos < 0] <- 1e4
d$daysRelativeToRejectionPos <- d$daysRelativeToRejection 
d$daysRelativeToRejectionPos[d$daysRelativeToRejectionPos < 0] <- 1e4
sampleIterations <- 200000; warmupIterations <- 100000; nChains <- 1
m38 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0 + bR0*(exp(1/p2tuner*log(daysRelativeToTransplantPos)))*(1-afterRejection)*afterTransplant + bR1*(exp(1/p3tuner*log(daysRelativeToRejectionPos)))*afterRejection,
    
    a0 <- A0b + A0j + A0i + a1*afterRejection,
    A0b ~ dnorm(0,4),
    A0j ~ dnorm(0,sigA0j),
    A0i ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,.1),
    
    a1 <- bR0*(exp(1/p2tuner*log(daysBetwTrans_Rej))) + A1b + A1j + A1i,
    A1b ~ dnorm(0,4),
    A1j ~ dnorm(0,sigA1j),
    A1i ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,.1),
    
    bR0 <- bR0b + bR0j + bR0i,
    bR0j ~ dnorm(0,sigbR0j),
    bR0i ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,0.5),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,.1),
    
    bR1 <- bR1b + bR1j + bR1i,
    bR1j ~ dnorm(0,sigbR1j),
    bR1i ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,2),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,.1),
    
    p2tuner ~ dgamma(10,20),
    p3tuner ~ dgamma(30,20)
    
  ) ,
  data= d[1:2,],
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains, verbose = T, 
  sample_file = paste0(dataSize, "/", "m38", dataSize), diagnostic_file = paste0(dataSize, "/", "m38diag_", dataSize),
  start = list(sigbR0j=1, sigbR0i=1, bR0b=0, sigbR1j=1, sigbR1i=1, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1, p2tuner= 0.5, p3tuner=1.5)
)
# save(file = paste0(dataSize, "/", "m38", dataSize), m38)
# q(save = "no")
precis(m38, depth = 2)

#using more informative priors??
#separate slopes for before and after the rejection event, seperate "intercepts" at the rejection event, intercept with zero slope before the transplant event,
#express slope in change per unit day between transplant and rejection event, express post-rejection intercept as displacement from pre-rejection intercept
#allow for exponential effect on relative change per day, learn exponent from data, but constraint it to a more reasonable, informative, regularizing range
nChains <- 1
d$daysRelativeToTransplantPos <- d$daysRelativeToTransplant 
d$daysRelativeToTransplantPos[d$daysRelativeToTransplantPos < 0] <- 1e4
d$daysRelativeToRejectionPos <- d$daysRelativeToRejection 
d$daysRelativeToRejectionPos[d$daysRelativeToRejectionPos < 0] <- 1e4
m39 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0 + bR0*(exp(1/p2tuner*log(daysRelativeToTransplantPos)))*(1-afterRejection)*afterTransplant + bR1*(exp(1/p3tuner*log(daysRelativeToRejectionPos)))*afterRejection,
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id] + a1*afterRejection,
    A0b ~ dnorm(0,2),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,0.1),
    
    a1 <- bR0*(exp(1/p2tuner*log(daysBetwTrans_Rej))) + A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,2),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,0.1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    bR0b ~ dnorm(0,0.5),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,0.1),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,2),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,0.1),
    
    p2tuner ~ dgamma(10,20),
    p3tuner ~ dgamma(30,20)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains, verbose = T, debug = T,
  start = list(sigbR0j=0.1, sigbR0i=0.1, bR0b=0, sigbR1j=0.1, sigbR1i=0.1, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1, p2tuner= 0.5, p3tuner=1.5)
)
save(file = paste0(dataSize, "/", "m39", dataSize), m39)
q(save = "no")
precis(m39, depth = 2)

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
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains, 
  start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, A0b=0, sigA0j=0.5, sigA0i=0.5, A1b=0, sigA1j=0.5, sigA1i=0.5)
)
save(file = paste0(dataSize, "/", "m40", dataSize), m40)
q(save = "no")
precis(m40, depth = 2)

#separate slopes for before and after the rejection event, seperate "intercepts" at the rejection event, intercept with zero slope before the transplant event,
#express slope in change per unit sqrt(day) between transplant and rejection event, express post-rejection intercept as displacement from pre-rejection intercept
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
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains, 
  start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, A0b=0, sigA0j=0.5, sigA0i=0.5, A1b=0, sigA1j=0.1, sigA1i=0.1)
)
save(file = paste0(dataSize, "/", "m41", dataSize), m41)
q(save = "no")
precis(m41, depth = 2)

# load("full_removeP9//m41full_removeP9")
load("full/m41full")
m41p <- link(m41, n = 2000)
hist(d$spectral_count, xlim = c(0,50), breaks = 200, col = rgb(0,0,1,alpha = 0.5), ylim = c(0,5000))
hist(rpois(n = length(m41p$lambda[1,]), lambda = m41p$lambda[1,]), add = T, col = rgb(1,0,0,alpha = 0.5), breaks = 200)
hist(sapply(1:length(d$spectral_count), function(x) ppois(d$spectral_count[x], m41p$lambda[,x])))
# hist(sapply(1:length(d$spectral_count), function(x) sum(rpois(2000,m41p$lambda[,x]) >= d$spectral_count[x])/2000))

#plot raw residuals from mean expectation
plot(sapply(1:length(m41p$lambda[1,]), function(x) mean(m41p$lambda[,x])), 
     sapply(1:length(m41p$lambda[1,]), function(x) mean(m41p$lambda[,x])) - d$spectral_count,
     col = rgb(0,0,0,0.3)) 

#plot pearson residuals -- mean lambda
plot(sapply(1:length(m41p$lambda[1,]), function(x) mean(m41p$lambda[,x])), 
     (sapply(1:length(m41p$lambda[1,]), function(x) mean(m41p$lambda[,x])) - d$spectral_count) / sqrt(sapply(1:length(m41p$lambda[1,]), function(x) mean(m41p$lambda[,x]))),
     xlab = "Predicted Rate", ylab = "Standardized Residual from Predicted Rate")
hist((sapply(1:length(m41p$lambda[1,]), function(x) mean(m41p$lambda[,x])) - d$spectral_count) / sqrt(sapply(1:length(m41p$lambda[1,]), function(x) mean(m41p$lambda[,x]))),
     breaks = 500, xlim = c(-5, 5), freq = F)
lines(seq(-5, 5, length.out = 200), dnorm(x = seq(-5, 5, length.out = 200), mean = 0, sd = 1))

#plot pearson residuals -- sampled lambda
numSamps <- 400
indices <- floor(seq(from = 1, to = length(m41p$lambda[,1]), length.out = numSamps))
resids <- rep(0, length(d$spectral_count)*length(indices))
for(x in 1:length(d$spectral_count)){
  if (x %% 100 == 0) {cat(paste(x, " "))}
  resids[((x-1)*numSamps+1):(x*numSamps)] <- (m41p$lambda[,x][indices] - d$spectral_count[x]) / sqrt(m41p$lambda[,x][indices])
}
hist(resids, breaks = 500, xlim = c(-5, 5), freq = F)
lines(seq(-5, 5, length.out = 200), dnorm(x = seq(-5, 5, length.out = 200), mean = 0, sd = 1))

#plot pearson residuals -- sampled lambda, nonzero counts only
numSamps <- 400
indices <- floor(seq(from = 1, to = length(m41p$lambda[,1]), length.out = numSamps))
resids <- rep(0, length(d$spectral_count)*length(indices))
for(x in 1:length(d$spectral_count)){
  if (x %% 100 == 0) {cat(paste(x, " "))}
  if(d$spectral_count[x] == 0) {resids[((x-1)*numSamps+1):(x*numSamps)] <- NA} else {
  resids[((x-1)*numSamps+1):(x*numSamps)] <- (m41p$lambda[,x][indices] - d$spectral_count[x]) / sqrt(m41p$lambda[,x][indices])}
}
hist(resids, breaks = 500, xlim = c(-5, 5), freq = F)
lines(seq(-5, 5, length.out = 200), dnorm(x = seq(-5, 5, length.out = 200), mean = 0, sd = 1))


m41s <- extract.samples(m41, clean.names = F)

#before rejection, after transplant
proteinEffects <- m41s$bR0j
proteinKey <- data.frame(as.character(d$protein), d$protein_id)
proteinKey <- proteinKey[order(proteinKey$d.protein_id),]
proteinKey <- unique(proteinKey)
hist(sapply(1:length(proteinEffects[1,]), function (x) mean(proteinEffects[,x])), breaks = 200)
dens(exp(m41s$bR0b + proteinEffects[,1]) * 100, xlim = c(90,120), xlab = "Percent Change in Spectral Count per sqrt(Day)", col = 0, ylim = c(0,0.5))
num <- 0
percChanges <- 0
proteinSpecificProbPositives <- rep(2, length(proteinEffects[1,]))
proteinSpecificMeanDailyChange <- rep(2, length(proteinEffects[1,]))
proteinList <- rep(2, length(proteinEffects[1,])) 
proteinsShared <- proteins[numPatientsPerProtein >= 3]
sdsEffect <- rep(2, length(proteinEffects[1,]))
for(j in 1:length(proteinEffects[1,])){
  percChange <- exp(m41s$bR0b + proteinEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  proteinSpecificProbPositives[j] <- probPos
  proteinSpecificMeanDailyChange[j] <- mean(percChange)
  proteinList[j] <- as.character(proteinKey[proteinKey$d.protein_id == j, 1])
  sdsEffect[j] <- sd(percChange)
  # cat(paste0(probPos, " "))
  if(probPos > 0.95 | probPos < 0.05){
    num <- num + 1
    cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
  }
}
num
hist(proteinSpecificProbPositives, breaks = 50)
hist(proteinSpecificMeanDailyChange, breaks = 50)
numPatientsPerProtein <- sapply(1:length(proteinList), 
                                function(x) length(unique(d[d$protein == proteinList[x] & d$spectral_count != 0,]$patient)))
dataSummary <- cbind(proteinList, proteinSpecificProbPositives, proteinSpecificMeanDailyChange, numPatientsPerProtein) 
colnames(dataSummary) <- c("Protein Name", "Probability Change is > 0%", "Posterior Mean Effect", "Number of Patients w/ Non-Zero Count")
dataSummary <- dataSummary[order(dataSummary[,2], decreasing = T),]
addPatientData <- T
if(addPatientData){
  dataSummary <- cbind(dataSummary, matrix(0, nrow = length(dataSummary[,1]), ncol = (3*5+1*2)))
  for(i in 1:length(dataSummary[,1])){
    protein <- dataSummary[i,1]
    proteinData <- d[d$protein == protein,]
    proteinData <- proteinData[order(proteinData[,2], proteinData[,4], decreasing = F),]
    spectralCounts <- proteinData$spectral_count
    #errorChecking
    if(i != 1 && names != paste0(proteinData[,2], "_timepoint_", proteinData[,4])){
      print("error error abort")
      break
    }
    names <- paste0(proteinData[,2], "_timepoint_", proteinData[,4])
    dataSummary[i, 5:21] <- spectralCounts
  }
  colnames(dataSummary) <- c("Protein Name", "Probability Change is > 0%", "Posterior Mean Effect", "Number of Patients w/ Non-Zero Count", names)
}
write.csv(x = dataSummary, file = "modelsqrt_populated0s_beforeRejectionAfterTransplant_posteriorRelativeSlopeSummary.csv")

youWantToRestrictToJustPhase2 <- F
if(youWantToRestrictToJustPhase2){
  numPatientsPerProtein <- sapply(1:length(proteinList), 
                                  function(x) length(unique(d[d$protein == proteinList[x] & d$spectral_count != 0 
                                                              & d$phase == 2,]$patient)))
  dataSummary <- cbind(proteinList, proteinSpecificProbPositives, proteinSpecificMeanDailyChange, numPatientsPerProtein) 
  colnames(dataSummary) <- c("Protein Name", "Probability Change is > 0%", "Posterior Mean Effect", "Number of Patients w/ Non-Zero Count")
  dataSummary <- dataSummary[order(dataSummary[,2], decreasing = T),]
  plot(dataSummary[,4], dataSummary[,2], ylab = "Probability Change is > 0%", xlab = "Number of Patients w/ Non-Zero Count in Phase 2, Specifically", xaxt = 'n')
  axis(side = 2, at = seq(0, 0.9, by = 0.1)); axis(side = 1, at = c(0,1,2,3))
} else {
  numPatientsPerProtein <- sapply(1:length(proteinList), 
                                  function(x) length(unique(d[d$protein == proteinList[x] & d$spectral_count != 0,]$patient)))
}
plot(dataSummary[,2], dataSummary[,3], ylab = "Posterior Mean Effect", xlab = "Probability Change is > 0%")
plot(dataSummary[,4], dataSummary[,2], ylab = "Probability Change is > 0%", xlab = "Number of Patients w/ Non-Zero Count")
plot(dataSummary[,4], sdsEffect, ylab = "Standard Deviation of Marginal Posterior", xlab = "Number of Patients w/ Non-Zero Count")
ggplot2::ggplot(data.frame(numPats = dataSummary[,4], sds = sdsEffect), aes(x = numPats, y = sds)) + theme_bw() + geom_violin() #compare standard deviations of marginal posterior
hist(exp(rnorm(n = 2e4, 0, 0.03))*100, breaks = 500, xlim = c(90,110), col = 2, ylim = c(0,500)) #prior
hist(exp(m41s$bR0b)*100, breaks = 1000, add = T, col = col2rgb("green", alpha = 0.5)) #posterior


#after rejection
proteinEffects <- m41s$bR1j
hist(sapply(1:length(proteinEffects[1,]), function (x) mean(proteinEffects[,x])), breaks = 200)
dens(exp(m41s$bR1b + proteinEffects[,1]) * 100, xlim = c(90,120), xlab = "Percent Change in Intercept with sqrt(Day)", col = 0, ylim = c(0,0.3))
num <- 0
proteinSpecificProbPositives <- rep(2, length(proteinEffects[1,]))
proteinSpecificMeanDailyChange <- rep(2, length(proteinEffects[1,]))
sdsEffect <- rep(2, length(proteinEffects[1,]))
for(j in 1:length(proteinEffects[1,])){
  percChange <- exp(m41s$bR1b + proteinEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  proteinSpecificProbPositives[j] <- probPos
  proteinSpecificMeanDailyChange[j] <- mean(percChange)
  sdsEffect[j] <- sd(percChange)
  if(probPos > 0.95 | probPos < 0.05){
    num <- num + 1
    cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
  }
}
num
hist(proteinSpecificProbPositives, breaks = 50)
hist(proteinSpecificMeanDailyChange, breaks = 50)
youWantToRestrictToJustPhase3 <- F
if(youWantToRestrictToJustPhase3){
  numPatientsPerProtein <- sapply(1:length(proteinList), 
                                  function(x) length(unique(d[d$protein == proteinList[x] & d$spectral_count != 0 
                                                              & d$phase == 3,]$patient)))
  dataSummary <- cbind(proteinList, proteinSpecificProbPositives, proteinSpecificMeanDailyChange, numPatientsPerProtein) 
  colnames(dataSummary) <- c("Protein Name", "Probability Change is > 0%", "Posterior Mean Effect", "Number of Patients w/ Non-Zero Count")
  dataSummary <- dataSummary[order(dataSummary[,2], decreasing = T),]
  plot(dataSummary[,4], dataSummary[,2], ylab = "Probability Change is > 0%", xlab = "Number of Patients w/ Non-Zero Count in Phase 2, Specifically", xaxt = 'n')
  axis(side = 2, at = seq(0, 0.9, by = 0.1)); axis(side = 1, at = c(0,1,2,3))
} else {
  numPatientsPerProtein <- sapply(1:length(proteinList), 
                                  function(x) length(unique(d[d$protein == proteinList[x] & d$spectral_count != 0,]$patient)))
}
dataSummary <- cbind(proteinList, proteinSpecificProbPositives, proteinSpecificMeanDailyChange, numPatientsPerProtein) 
colnames(dataSummary) <- c("Protein Name", "Probability Change is > 0%", "Posterior Mean Effect", "Number of Patients w/ Non-Zero Count")
dataSummary <- dataSummary[order(dataSummary[,2], decreasing = T),]
addPatientData <- T
if(addPatientData){
  dataSummary <- cbind(dataSummary, matrix(0, nrow = length(dataSummary[,1]), ncol = (3*5+1*2)))
  for(i in 1:length(dataSummary[,1])){
    protein <- dataSummary[i,1]
    proteinData <- d[d$protein == protein,]
    proteinData <- proteinData[order(proteinData[,2], proteinData[,4], decreasing = F),]
    spectralCounts <- proteinData$spectral_count
    #errorChecking
    if(i != 1 && names != paste0(proteinData[,2], "_timepoint_", proteinData[,4])){
      print("error error abort")
      break
    }
    names <- paste0(proteinData[,2], "_timepoint_", proteinData[,4])
    dataSummary[i, 5:21] <- spectralCounts
  }
  colnames(dataSummary) <- c("Protein Name", "Probability Change is > 0%", "Posterior Mean Effect", "Number of Patients w/ Non-Zero Count", names)
}
write.csv(x = dataSummary, file = "modelsqrt_populated0s_afterRejection_posteriorRelativeSlopeSummary.csv")


plot(dataSummary[,2], dataSummary[,3], ylab = "Posterior Mean Effect", xlab = "Probability Change is > 0%")
plot(dataSummary[,4], dataSummary[,2], ylab = "Probability Change is > 0%", xlab = "Number of Patients w/ Non-Zero Count")
ggplot2::ggplot(data.frame(numPats = dataSummary[,4], sds = sdsEffect), aes(x = numPats, y = sds)) + theme_bw() + geom_violin() #compare standard deviations of marginal posterior
hist(exp(rnorm(n = 2e4, 0, 0.03))*100, breaks = 100, xlim = c(90,110), col = 2)
hist(exp(m41s$bR1b)*100, breaks = 100, add = T, col = col2rgb("green", alpha = 0.5))


#discontinuity at rejection
proteinEffects <- m41s$A1j
hist(sapply(1:length(proteinEffects[1,]), function (x) mean(proteinEffects[,x])), breaks = 200)
dens(exp(m41s$A1b + proteinEffects[,1]) * 100, xlim = c(0,300), xlab = "Percent Change in Intercept at Rejection", col = 0, ylim = c(0,0.05))
num <- 0
proteinSpecificProbPositives <- rep(2, length(proteinEffects[1,]))
proteinSpecificMeanChange <- rep(2, length(proteinEffects[1,]))
for(j in 1:length(proteinEffects[1,])){
  percChange <- exp(m41s$A1b + proteinEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  proteinSpecificProbPositives[j] <- probPos
  proteinSpecificMeanChange[j] <- mean(percChange)
  # cat(paste(probPos, "\t"))
  if(probPos > 0.95 | probPos < 0.05){
    num <- num + 1
    cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
  }
}
num
hist(proteinSpecificProbPositives, breaks = 50)
hist(proteinSpecificMeanChange, breaks = 50)
dataSummary <- cbind(proteinList, proteinSpecificProbPositives, proteinSpecificMeanChange, numPatientsPerProtein) 
colnames(dataSummary) <- c("Protein Name", "Probability Change is > 0%", "Posterior Mean Effect", "Number of Patients w/ Non-Zero Count")
dataSummary <- dataSummary[order(as.numeric(dataSummary[,2]), decreasing = F),]
addPatientData <- T
if(addPatientData){
  dataSummary <- cbind(dataSummary, matrix(0, nrow = length(dataSummary[,1]), ncol = (3*5+1*2)))
  for(i in 1:length(dataSummary[,1])){
    protein <- dataSummary[i,1]
    proteinData <- d[d$protein == protein,]
    proteinData <- proteinData[order(proteinData[,2], proteinData[,4], decreasing = F),]
    spectralCounts <- proteinData$spectral_count
    #errorChecking
    if(i != 1 && names != paste0(proteinData[,2], "_timepoint_", proteinData[,4])){
      print("error error abort")
      break
    }
    names <- paste0(proteinData[,2], "_timepoint_", proteinData[,4])
    dataSummary[i, 5:21] <- spectralCounts
  }
  colnames(dataSummary) <- c("Protein Name", "Probability Change is > 0%", "Posterior Mean Effect", "Number of Patients w/ Non-Zero Count", names)
}
write.csv(x = dataSummary, file = "modelsqrt_populated0s_displacementAtRejection_posteriorRelativeChangeSummary.csv")

plot(dataSummary[,2], dataSummary[,3], ylab = "Posterior Mean Effect", xlab = "Probability Change is > 0%")
plot(dataSummary[,4], dataSummary[,2], ylab = "Probability Change is > 0%", xlab = "Number of Patients w/ Non-Zero Count")

# probing specific proteins
targetProteins <- c("troponin", "prohibitin", "plakoglobin")
targetProteins <- proteinKey[unlist(sapply(1:length(targetProteins), function(x) which(grepl(targetProteins[x], proteinKey$as.character.d.protein)))),]
proteinEffects <- m41s$bR0j
for(j in 1:length(targetProteins[,1])){
  percChange <- exp(m41s$bR0b + proteinEffects[,targetProteins[j,2]]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == targetProteins[j,2], 1]), ", probability positive = ", probPos * 100, "%, mean change = ", round(mean(percChange), 3), "%\n"))
}
proteinEffects <- m41s$bR1j
for(j in 1:length(targetProteins[,1])){
  percChange <- exp(m41s$bR1b + proteinEffects[,targetProteins[j,2]]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == targetProteins[j,2], 1]), ", probability positive = ", probPos * 100, "%, mean change = ", round(mean(percChange), 3), "%\n"))
}
proteinEffects <- m41s$A1j
for(j in 1:length(targetProteins[,1])){
  percChange <- exp(m41s$A1b + proteinEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  cat(paste0(as.character(proteinKey[proteinKey$d.protein_id == targetProteins[j,2], 1]), ", probability positive = ", probPos * 100, "%, mean change = ", round(mean(percChange), 3), "%\n"))
}

#looking for nonidentifiability in patient 9
cors <- rep(0, 1770)
for(i in 1:length(cors)){cors[i] <- cor(m41s$A0i[,6], m41s$A0j[,i])}
hist(cors)
for(i in 1:length(cors)){cors[i] <- cor(m41s$bR0i[,6], m41s$bR0j[,i])}
hist(cors)
for(i in 1:length(cors)){cors[i] <- cor(m41s$bR1i[,6], m41s$bR1j[,i])}
hist(cors)



#individual effects
patientEffects <- m41s$bR0i
patientKey <- data.frame(as.character(d$patient), d$patient_id)
patientKey <- patientKey[order(patientKey$d.patient_id),]
patientKey <- unique(patientKey)
dens(exp(m41s$bR0b + patientEffects[,1]) * 100, xlim = c(99,101), xlab = "Percent Change in Spectral Count per Day Going Backward", col = 0, ylim = c(0,50))
for(j in 1:length(patientEffects[1,])){
  percChange <- exp(m41s$bR0b + patientEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  if(probPos > -1){
    cat(paste0(as.character(patientKey[patientKey$d.patient_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
    kde <- density(percChange)
    text(x = kde$x[which.max(kde$y)], y = max(kde$y)*1.05, labels = as.character(patientKey[patientKey$d.patient_id == j, 1]))
  }
}

patientEffects <- m41s$bR1i
dens(exp(m41s$bR1b + patientEffects[,1]) * 100, xlim = c(99.5,100.5), xlab = "Percent Change in Spectral Count per Day Going Forward", col = 0, ylim = c(0,200))
for(j in 1:length(patientEffects[1,])){
  percChange <- exp(m41s$bR1b + patientEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  if(probPos > -1){
    cat(paste0(as.character(patientKey[patientKey$d.patient_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
    dens(percChange, add = T)
    kde <- density(percChange)
    text(x = kde$x[which.max(kde$y)], y = max(kde$y)*1.05, labels = as.character(patientKey[patientKey$d.patient_id == j, 1]))
  }
}

patientEffects <- m41s$A1i
hist(sapply(1:length(patientEffects[1,]), function (x) mean(patientEffects[,x])), breaks = 200)
dens(exp(m41s$A1b + patientEffects[,1]) * 100, xlim = c(0,300), xlab = "Percent Change in Intercept at Rejection", col = 0, ylim = c(0,0.08))
patientSpecificProbPositives <- rep(2, length(patientEffects[1,]))
patientSpecificMeanDailyChange <- rep(2, length(patientEffects[1,]))
for(j in 1:length(patientEffects[1,])){
  percChange <- exp(m41s$A1b + patientEffects[,j]) * 100
  probPos <- sum(percChange > 100)/length(percChange)
  patientSpecificProbPositives[j] <- probPos
  patientSpecificMeanDailyChange[j] <- mean(percChange)
  # cat(paste(probPos, "\t"))
  cat(paste0(as.character(patientKey[patientKey$d.patient_id == j, 1]), ", probability positive = ", probPos * 100, "%\n"))
  dens(percChange, add = T)
  kde <- density(percChange)
  text(x = kde$x[which.max(kde$y)], y = max(kde$y)*1.05, labels = as.character(patientKey[patientKey$d.patient_id == j, 1]))
}


#plot lines for uncertainty in expected rate through time relative to rejection

#mean estimate
par(mar = c(4,4,2,2))
plot(c(-1300,-441), exp(mean(m41s$A0b) + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,3), 
     xlab = "days relative to rejection event", ylab = "rate of spectral count for the average protein/patient", lwd = 2)
lines(c(-441:0), exp(mean(m41s$A0b) + sqrt(c(0:441)) * mean(m41s$bR0b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
lines(c(0:1100), exp(mean(m41s$A1b) + mean(m41s$A0b) + sqrt(441)*mean(m41s$bR0b) + sqrt(c(0:1100)) * mean(m41s$bR1b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
title("Expected Rates of Spectral Count Given Days Relative to Acute Rejection Event")
#50% hpdi
lines(c(-1300,-441), exp(HPDI(prob = 0.5, m41s$A0b)[[1]] + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(-441:0), exp(HPDI(prob = 0.5, m41s$A0b)[[1]] + sqrt(c(0:441)) * HPDI(prob = 0.5, m41s$bR0b)[[1]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(0:1100), exp(HPDI(prob = 0.5, m41s$A0b)[[1]] + sqrt(441)*HPDI(prob = 0.5, m41s$bR0b)[[1]] + HPDI(prob = 0.5, m41s$A1b)[[1]] + sqrt(c(0:1100)) * HPDI(prob = 0.5, m41s$bR1b)[[1]]), 
      type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(-1300,-441), exp(HPDI(prob = 0.5, m41s$A0b)[[2]] + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(-441:0), exp(HPDI(prob = 0.5, m41s$A0b)[[2]] + sqrt(c(0:441)) * HPDI(prob = 0.5, m41s$bR0b)[[2]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(0:1100), exp(HPDI(prob = 0.5, m41s$A0b)[[2]] + sqrt(441)*HPDI(prob = 0.5, m41s$bR0b)[[2]] + HPDI(prob = 0.5, m41s$A1b)[[2]] + sqrt(c(0:1100)) * HPDI(prob = 0.5, m41s$bR1b)[[2]]), 
      type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
#rejection date
abline(v = 0, col = 2, lwd = 2)
text(x = -25, y = 2, labels = "Rejection Date", srt = 90)
#transplant date
abline(v = -441, col = 2, lwd = 2)
text(x = -466, y = 2, labels = "Average Transplant Date", srt = 90)
#samples from posterior
for(i in seq(1, length(m41s$A0b), length.out = 2000)){
  lines(c(-1300,-441), exp(m41s$A0b[i] + c(0,0)), type = "l", col = rgb(0,0,0,0.025), xlim = c(-1300, 1100), ylim = c(0,10))
  lines(c(-441:0), exp(m41s$A0b[i] + sqrt(c(0:441)) * m41s$bR0b[i]), type = "l", col = rgb(0,0,0,0.025), xlim = c(-1300, 1100), ylim = c(0,10))
  lines(c(0:1100), exp(m41s$A0b[i] + sqrt(441)*m41s$bR0b[i] + m41s$A1b[i] + sqrt(c(0:1100)) * m41s$bR1b[i]), type = "l", col = rgb(0,0,0,0.025), xlim = c(-1300, 1100), ylim = c(0,10))
}
#boxplots of data
times <- sort(unique(d$daysRelativeToRejection))
for(i in 1:length(times)){
  boxplot(d$spectral_count[d$daysRelativeToRejection == times[i]], add = T, at = times[i], boxwex=20, outline = F, col.axis = rgb(0,0,0,0))
}


#same slopes for before and after the rejection event, seperate intercepts at the rejection event
nChains <- 1
m14T <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0*(1-afterRejection) + a1*afterRejection + bR*((daysRelativeToRejection - daysRelativeToTransplant*(1-afterTransplant))/daysBetwTrans_Rej), 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id],
    A0b ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,2),
    
    a1 <- A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,2),
    
    bR <- bRb + bRj[protein_id] + bRi[patient_id],
    bRj[protein_id] ~ dnorm(0,sigbRj),
    bRi[patient_id] ~ dnorm(0,sigbRi),
    bRb ~ dnorm(0,2),
    c(sigbRj, sigbRi) ~ dcauchy(0,2)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigbRj=1, sigbRi=1, bRb=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
save(file = paste0(dataSize, "/", "m14T", dataSize), m14T)
precis(m14T, depth = 2)
#plot(m14T)
#pairs(m14T)

#similar slopes for before and after the rejection event, seperate intercepts at the rejection event
nChains <- 1
m15T <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0*(1-afterRejection) + a1*afterRejection + bR*((daysRelativeToRejection - daysRelativeToTransplant*(1-afterTransplant))/daysBetwTrans_Rej), 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id],
    A0b ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,2),
    
    a1 <- A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,2),
    
    bR <- bRb + bRj[protein_id] + bRi[patient_id] + bRAfter*afterRejection,
    bRj[protein_id] ~ dnorm(0,sigbRj),
    bRi[patient_id] ~ dnorm(0,sigbRi),
    bRb ~ dnorm(0,2),
    c(sigbRj, sigbRi) ~ dcauchy(0,2),
    bRAfter ~ dnorm(0, sigbRAfter),
    sigbRAfter ~ dcauchy(0,2)
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigbRj=1, sigbRi=1, sigbRAfter = 1, bRb=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
save(file = paste0(dataSize, "/", "m15T", dataSize), m15T)
precis(m15T, depth = 2)
#plot(m15T)
#pairs(m15T)

#similar slopes for before and after the rejection event, but same intercept at the rejection event
nChains <- 1
m16T <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a + bR*((daysRelativeToRejection - daysRelativeToTransplant*(1-afterTransplant))/daysBetwTrans_Rej), 
    
    a <- Ab + Aj[protein_id] + Ai[patient_id],
    Ab ~ dnorm(0,4),
    Aj[protein_id] ~ dnorm(0,sigAj),
    Ai[patient_id] ~ dnorm(0,sigAi),
    c(sigAj, sigAi) ~ dcauchy(0,2),
    
    bR <- bRb + bRj[protein_id] + bRi[patient_id] + bRAfter*afterRejection,
    bRj[protein_id] ~ dnorm(0,sigbRj),
    bRi[patient_id] ~ dnorm(0,sigbRi),
    bRb ~ dnorm(0,2),
    c(sigbRj, sigbRi) ~ dcauchy(0,2),
    bRAfter ~ dnorm(0, sigbRAfter),
    sigbRAfter ~ dcauchy(0,2)
    
  ) ,
  data= d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
  start = list(sigAj=1, sigAi=1, sigbRj=1, sigbRi=1, sigbRAfter = 1, Ab=0)
)
save(file = paste0(dataSize, "/", "m16T", dataSize), m16T)
precis(m16T, depth = 2)
#plot(m16T)
#pairs(m16T)

#separate slopes for before and after the rejection event, seperate intercepts at the rejection event, 
#data before the transplant informs starting point at the transplant to the rejection
# nChains <- 1
# m17T <- map2stan(
#   alist(
#     spectral_count ~ dpois( lambda ),
#     log(lambda) <- aBT*(1-afterTransplant) + a0*(1-afterRejection)*afterTransplant + a1*afterRejection + 
#       bR0*daysRelativeToRejection/daysBetwTrans_Rej*(1-afterRejection)*afterTransplant + bR1*daysRelativeToRejection*afterRejection, 
#     
#     aBT <- aBTb + aBTj[protein_id] + aBTi[patient_id],
#     aBTb ~ dnorm(0,4),
#     aBTj[protein_id] ~ dnorm(0,sigaBTj),
#     aBTi[patient_id] ~ dnorm(0,sigaBTi),
#     c(sigaBTj, sigaBTi) ~ dcauchy(0,2),
#     
#     a0 <- aBT + bR0,
#     
#     a1 <- A1b + A1j[protein_id] + A1i[patient_id],
#     A1b ~ dnorm(0,4),
#     A1j[protein_id] ~ dnorm(0,sigA1j),
#     A1i[patient_id] ~ dnorm(0,sigA1i),
#     c(sigA1j, sigA1i) ~ dcauchy(0,2),
#     
#     bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id],
#     bR0j[protein_id] ~ dnorm(0,sigbR0j),
#     bR0i[patient_id] ~ dnorm(0,sigbR0i),
#     bR0b ~ dnorm(0,2),
#     c(sigbR0j, sigbR0i) ~ dcauchy(0,2),
#     
#     bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
#     bR1j[protein_id] ~ dnorm(0,sigbR1j),
#     bR1i[patient_id] ~ dnorm(0,sigbR1i),
#     bR1b ~ dnorm(0,03),
#     c(sigbR1j, sigbR1i) ~ dcauchy(0,0.1)
#     
#   ) ,
#   data= d,
#   iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains,
#   start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, aBTb=0, sigaBTj=1, sigaBTi=1, A1b=0, sigA1j=1, sigA1i=1)
# )
# save(file = paste0(dataSize, "/", "m17T", dataSize), m17T)
# precis(m17T, depth = 2)



#####################################################
#extroducing the time relative to transplant variable
#####################################################


############
#here is a model where you can have patient and protein-specific responses to the timepoint's distance in days from the rejection event
#or a linear effect, marginalizing over process uncertainty

N <- length(d$spectral_count)
N_miss <- N / 1.1
x <- rbinom( N , size=1 , prob=0.5 )
x[ sample(1:N,size=(N-2)) ] <- NA

y <- rnorm( N , 2*x , 1 )

f6 <- alist(
  y ~ dnorm( mu , sigma ),
  mu <- a + b*x,
  x ~ bernoulli( phi ),
  a ~ dnorm( 0 , 100 ),
  b ~ dnorm( 0  , 10 ),
  phi ~ beta( 1 , 1 ),
  sigma ~ dcauchy(0,2)
)
m6 <- map2stan( f6 , data=list(y=y,x=x) , constraints=list(phi="lower=0,upper=1"), do_discrete_imputation=TRUE )


sampleIterations <- 20
warmupIterations <- 10
nChains <- 1
pivot <- rep(NA, length(d$spectral_count)); pivot[1:2] <- c(0,1)
mX <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a + bRpiv*fabs(daysRelativeToRejection)*pivot + bRlin*daysRelativeToRejection*(1-pivot), 
    pivot ~ bernoulli(phi),
    phi ~ beta(1,1),
    a <- Ab,
    Ab ~ dnorm(0,5),
    bRlin <- bRlinb,
    bRpiv <- bRpivb,
    c(bRpivb, bRlinb) ~ dnorm(0,0.1)
  ) ,
  data = list(pivot=pivot, spectral_count=d$spectral_count, daysRelativeToRejection = d$daysRelativeToRejection, patient_id = d$patient_id, protein_id = d$protein_id),
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains, do_discrete_imputation=TRUE,
  start = list(sigAj=1, sigAi=1, sigbRj=0.01, sigbRi=0.01, bRb=0, Ab=0), constraints=list(phi="lower=0,upper=1")
)
precis(mX, depth = 2)
#plot(mX)
mX@model

"
data{
int<lower=1> N;
int<lower=1> N_pivot_whole;
int<lower=1> N_pivot;
int spectral_count[N];
int pivot[N_pivot];
int pivot_whole[N_pivot_whole];
real daysRelativeToRejection[N];
matrix[2,1] lambda_missmatrix;
}
parameters{
real sigAj;
real sigAi;
real sigbRj;
real sigbRi;
real bRb;
real Ab;
real<lower=0,upper=1> phi;
real bRpivb;
real bRlinb;
}
model{
vector[N] bRpiv;
vector[N] bRlin;
vector[N] a;
matrix[N,2] lambda_mxt;
matrix[N,2] lambda_mxtlm;
vector[N] lambda;
bRlinb ~ normal( 0 , 0.1 );
bRpivb ~ normal( 0 , 0.1 );
for ( i in 1:N ) {
bRpiv[i] = bRpivb;
}
for ( i in 1:N ) {
bRlin[i] = bRlinb;
}
Ab ~ normal( 0 , 5 );
for ( i in 1:N ) {
a[i] = Ab;
}
phi ~ beta( 1 , 1 );
pivot_whole ~ bernoulli( phi );
for ( i in 1:N ) {
if ( pivot[i]<0 ) {
for ( mmrow in 1:rows(lambda_missmatrix) ) {
vector[1] dimiproxy;
dimiproxy[1] = pivot[i];
if ( pivot[i]<0 ) dimiproxy[1] = lambda_missmatrix[mmrow,1];

lambda_mxtlm[i,mmrow] = a[i] + bRpiv[i] * fabs(daysRelativeToRejection[i]) * dimiproxy[1] + bRlin[i] * daysRelativeToRejection[i] *      (1 - dimiproxy[1]);
lambda_mxt[i,mmrow] = 0;
if ( lambda_missmatrix[mmrow,1]==1 )
lambda_mxt[i,mmrow] = lambda_mxt[i,mmrow] + log(phi);
else
lambda_mxt[i,mmrow] = lambda_mxt[i,mmrow] + log1m(phi);

}//mmrow
}//if
lambda[i] = a[i] + bRpiv[i] * fabs(daysRelativeToRejection[i]) * pivot[i] + bRlin[i] * daysRelativeToRejection[i] *      (1 - pivot[i]);
}
for ( i in 1:N ) {
if ( pivot[i]<0 ) {
for ( mmrow in 1:rows(lambda_missmatrix) )
lambda_mxt[i,mmrow] = lambda_mxt[i,mmrow] + poisson_log_lpmf( spectral_count[i] | lambda_mxtlm[i,mmrow] );
target += log_sum_exp(lambda_mxt[i]);
} else 
spectral_count[i] ~ poisson_log( lambda[i] );
}//i 
}
generated quantities{
vector[N] bRpiv;
vector[N] bRlin;
vector[N] a;
matrix[N,2] lambda_mxt;
matrix[N,2] lambda_mxtlm;
vector[N] lambda;
vector[N] pivot_impute;
for ( i in 1:N ) {
bRpiv[i] = bRpivb;
}
for ( i in 1:N ) {
bRlin[i] = bRlinb;
}
for ( i in 1:N ) {
a[i] = Ab;
}
for ( i in 1:N ) {
if ( pivot[i]<0 ) {
for ( mmrow in 1:rows(lambda_missmatrix) ) {
vector[1] dimiproxy;
dimiproxy[1] = pivot[i];
if ( pivot[i]<0 ) dimiproxy[1] = lambda_missmatrix[mmrow,1];

lambda_mxtlm[i,mmrow] = a[i] + bRpiv[i] * fabs(daysRelativeToRejection[i]) * dimiproxy[1] + bRlin[i] * daysRelativeToRejection[i] *      (1 - dimiproxy[1]);
lambda_mxt[i,mmrow] = 0;
if ( lambda_missmatrix[mmrow,1]==1 )
lambda_mxt[i,mmrow] = lambda_mxt[i,mmrow] + log(phi);
else
lambda_mxt[i,mmrow] = lambda_mxt[i,mmrow] + log1m(phi);

}//mmrow
}//if
lambda[i] = a[i] + bRpiv[i] * fabs(daysRelativeToRejection[i]) * pivot[i] + bRlin[i] * daysRelativeToRejection[i] *      (1 - pivot[i]);
}
for ( i in 1:N ) {
if ( pivot[i]<0 ) {
for ( mmrow in 1:rows(lambda_missmatrix) )
lambda_mxt[i,mmrow] = lambda_mxt[i,mmrow] + poisson_log_lpmf( spectral_count[i] | lambda_mxtlm[i,mmrow] );
pivot_impute[i] = pivot[i];
if ( pivot[i]<0 )pivot_impute[i] = dot_product( softmax(to_vector(lambda_mxt[i])) , lambda_missmatrix[:,1] );
} else {
pivot_impute[i] = pivot[i];
}
}//i 
}
"

#with pivot bernoulli
sampleIterations <- 20
warmupIterations <- 10
nChains <- 1
pivot <- rep(NA, length(d$spectral_count)); pivot[1:2] <- c(0,1)
mX <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a + bRpiv*fabs(daysRelativeToRejection)*pivot + bRlin*daysRelativeToRejection*(1-pivot), 
    pivot ~ bernoulli(phi),
    phi ~ beta(1,1),
    
    a <- Ab + Aj[protein_id] + Ai[patient_id],
    Ab ~ dnorm(0,5),
    Aj[protein_id] ~ dnorm(0,sigAj),
    Ai[patient_id] ~ dnorm(0,sigAi),
    c(sigAj, sigAi) ~ dcauchy(0,1),
    
    bRlin <- bRlinb + bRlinj[protein_id] + bRlini[patient_id],
    bRlinj[protein_id] ~ dnorm(0,sigbRlinj),
    bRlini[patient_id] ~ dnorm(0,sigbRlini),
    
    bRpiv <- bRpivb + bRpivj[protein_id] + bRpivi[patient_id],
    bRpivj[protein_id] ~ dnorm(0,sigbRpivj),
    bRpivi[patient_id] ~ dnorm(0,sigbRpivi),
    
    c(bRpivb, bRlinb) ~ dnorm(0,0.1),
    c(sigbRlinj, sigbRlini, sigbRpivj, sigbRpivi) ~ dcauchy(0,0.2)
  ) ,
  data = list(pivot=pivot, spectral_count=d$spectral_count, daysRelativeToRejection = d$daysRelativeToRejection, patient_id = d$patient_id, protein_id = d$protein_id),
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains, do_discrete_imputation=TRUE,
  start = list(sigAj=1, sigAi=1, sigbRj=0.01, sigbRi=0.01, bRb=0, Ab=0), constraints=list(phi="lower=0,upper=1")
)
precis(mX, depth = 2)
#plot(mX)
#pairs(mX)

#with pivot or not pivot weights
sampleIterations <- 2500
warmupIterations <- 1000
nChains <- 1
mX <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a + bRpiv*fabs(daysRelativeToRejection)*pivotWeightAct + bRlin*daysRelativeToRejection*(1-pivotWeightAct),
    pivotWeightAct ~ dbeta(shape1, 0.5-shape1),
    shape1 ~ dunif(0.1, 0.49),
    
    a <- Ab + Aj[protein_id] + Ai[patient_id],
    Ab ~ dnorm(0,5),
    Aj[protein_id] ~ dnorm(0,sigAj),
    Ai[patient_id] ~ dnorm(0,sigAi),
    c(sigAj, sigAi) ~ dcauchy(0,1),
    
    bRlin <- bRlinb + bRlinj[protein_id] + bRlini[patient_id],
    bRlinj[protein_id] ~ dnorm(0,sigbRlinj),
    bRlini[patient_id] ~ dnorm(0,sigbRlini),
    
    bRpiv <- bRpivb + bRpivj[protein_id] + bRpivi[patient_id],
    bRpivj[protein_id] ~ dnorm(0,sigbRpivj),
    bRpivi[patient_id] ~ dnorm(0,sigbRpivi),
    
    c(bRpivb, bRlinb) ~ dnorm(0,0.1),
    c(sigbRlinj, sigbRlini, sigbRpivj, sigbRpivi) ~ dcauchy(0,0.2)
  ) ,
  data = d,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains, 
  start = list(sigAj=1, sigAi=1, sigbRlinj=0.01, sigbRlini=0.01, sigbRpivj=0.01, sigbRpivi=0.01, bRlinb=0, bRpivb=0, Ab=2)
)
precis(mX, depth = 2)
#plot(mX)
#pairs(mX)

#model comparison
if(Sys.info()['sysname'] == "Darwin") {
  setwd("/Volumes/macOS/Users/nikolai/heart_transplant/")
} else {
  setwd("C:\\Users\\Nikolai\\Desktop\\Kate Dissertation\\CSVs")
}
dataSizeIndex <- 5
dataSize <- c("full", "default", "sparse", "defaultfixed", "full_removeP9")[dataSizeIndex]
output <- paste0(dataSize, "/", list.files(dataSize))
paste(1:length(output), output)
output <- output[c(13:22)]
output <- output[c(1,2,4)]
for(i in 1:length(output)){
  print(paste(i, "/", length(output))); load(output[i])
}
fitmodels <- ls()[startsWith(ls(), "m")]
eval(parse(text=paste0("compare(", paste0(paste0(fitmodels[-length(fitmodels)], ",", collapse = ""), fitmodels[length(fitmodels)]), ")")))

# default dataset analysis
# WAIC  pWAIC  dWAIC weight     SE    dSE
# m27   46959.0 3397.0    0.0   0.59 322.14     NA
# m27c2 46959.8 3396.4    0.8   0.40 322.23   3.87
# m32c2 46967.0 3509.6    7.9   0.01 319.54  24.62
# m32   46971.1 3511.8   12.0   0.00 319.47  25.27
# m30c2 47668.2 3469.0  709.2   0.00 348.27  88.13
# m30   47675.4 3476.5  716.3   0.00 347.43  86.78
# m28   47692.1 3329.9  733.1   0.00 362.98 182.59
# m29   50216.9 3959.0 3257.9   0.00 426.74 274.76
# m31   53306.2 3868.5 6347.1   0.00 488.16 320.81
# m31c2 53314.2 3871.2 6355.1   0.00 488.92 321.66

# full dataset analysis
# WAIC  pWAIC   dWAIC weight      SE    dSE
# m28c2 60332.4 4798.1     0.0      1  660.67     NA
# m27c2 60850.9 5272.5   518.5      0  649.18 354.81
# m23c2 75248.5 7011.1 14916.1      0  987.30 622.43
# m29c2 75675.7 6741.4 15343.4      0  979.85 584.76
# m31c2 81925.3 6226.8 21592.9      0 1116.88 725.62

samps <- extract.samples(m40)
sapply(1:14, function(x) effectiveSize(samps[[x]]))
sapply(1:14, function(x) geweke.diag(samps[[x]]))


#trace plots
par(mfrow = c(3,5))
for(i in 1:14){plot(1:25000, samps[i][[1]], type = "l", ylab = names(samps)[i])}
#bivariate scatterplots

par(mfrow = c(3,5))
for(j in 13:13){
par(mfrow = c(3,5))
for(i in 1:14){plot(samps[j][[1]][seq(1, 25000, by = 10)], samps[i][[1]][seq(1, 25000, by = 10)], ylab = names(samps)[i], xlab = names(samps)[j])}
}
