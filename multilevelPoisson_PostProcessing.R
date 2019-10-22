######################################################################
###################### Loading in the raw data #######################
######################################################################
if(!exists("d")){
for(i in 1:1){ #hacky way of loading everything

  library(rethinking)
  keratinIndex <- 1
  keratin <- c("removeCorrectly", "removeIncorrectly")[keratinIndex]
  dataSizeIndex <- 1
  dataSize <- c("full", "default", "sparse", "defaultfixed")[dataSizeIndex]
  EDA <- F
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
}
# ######################################################################
# ################## NOW LET"S DO SOME MODEL CHECKING ##################
# ######################################################################
# 
# #stairtep model, m31
# #stairstep intercepts for all three phases, pooling across individuals, adjascent phases, and proteins
# load("full/m31full")
# load("full/m31c2full")
# #postcheck(fit = m31, prob = 0.95, window =  1000)
# sflist <- c(m31@stanfit, m31c2@stanfit)
# m31comb <- sflist2stanfit(sflist); rm(sflist)
# meanDisagree <- ((summary(m31@stanfit)[[1]][,1] - summary(m31c2@stanfit)[[1]][,1])/summary(m31@stanfit)[[1]][,1])
# hist(as.vector(meanDisagree), breaks = 1000, xlim = c(-5, 5))
# max(meanDisagree); min(meanDisagree)
# 
# m31s <- extract.samples(m31)
# #mean estimate
# par(mar = c(4,4,2,2))
# plot(c(-1300,-441), exp(mean(m31s$A0b) + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,1.5),
#      xlab = "days relative to rejection event", ylab = "rate of spectral count for the average protein/patient", lwd = 2)
# lines(c(-441,0), exp(mean(m31s$A0b) + mean(m31s$A1b) + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
# lines(c(0,1100), exp(mean(m31s$A0b) + mean(m31s$A1b) + mean(m31s$A2b) + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
# title("Expected Rates of Spectral Count Given Days Relative to Acute Rejection Event")
# #50% hpdi
# lines(c(-1300,-441), exp(HPDI(prob = 0.5, m31s$A0b)[[1]] + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# lines(c(-441,0), exp(HPDI(prob = 0.5, m31s$A0b + m31s$A1b)[[1]]) + c(0,0), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# lines(c(0,1100), exp(HPDI(prob = 0.5, m31s$A0b + m31s$A1b + m31s$A2b)[[1]]) + c(0,0), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# lines(c(-1300,-441), exp(HPDI(prob = 0.5, m31s$A0b)[[2]]) + c(0,0), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# lines(c(-441,0), exp(HPDI(prob = 0.5, m31s$A0b + m31s$A1b)[[2]]) + c(0,0), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# lines(c(0,1100), exp(HPDI(prob = 0.5, m31s$A0b + m31s$A1b + m31s$A2b)[[2]]) + c(0,0), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# #rejection date
# abline(v = 0, col = 2, lwd = 2)
# text(x = -25, y = 1, labels = "Rejection Date", srt = 90)
# #transplant date
# abline(v = -441, col = 2, lwd = 2)
# text(x = -466, y = 1, labels = "Average Transplant Date", srt = 90)
# #samples from posterior
# for(i in seq(1, length(m31s$A0b), length.out = 2000)){
#   lines(c(-1300,-441), exp(m31s$A0b[i] + c(0,0)), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
#   lines(c(-441,0), exp(m31s$A0b[i] + m31s$A1b[i] + c(0,0)), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
#   lines(c(0,1100), exp(m31s$A0b[i] + m31s$A1b[i]+ m31s$A2b[i] + c(0,0)), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
# }
# 
# #linear on days model, m40
# #separate slopes for before and after the rejection event, seperate "intercepts" at the rejection event, intercept with zero slope before the transplant event,
# #express slope in change per unit day between transplant and rejection event, express post-rejection intercept as displacement from pre-rejection intercept
# #this is model 27 with tighter priors
# load("full/m40full")
# m40s <- extract.samples(m40)
# #mean estimate
# par(mar = c(4,4,2,2))
# plot(c(-1300,-441), exp(mean(m40s$A0b) + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,1.5),
#      xlab = "days relative to rejection event", ylab = "rate of spectral count for the average protein/patient", lwd = 2)
# lines(c(-441:0), exp(mean(m40s$A0b) + c(0:441) * mean(m40s$bR0b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
# lines(c(0:1100), exp(mean(m40s$A1b) + mean(m40s$A0b) + 441*mean(m40s$bR0b) + c(0:1100) * mean(m40s$bR1b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
# title("Expected Rates of Spectral Count Given Days Relative to Acute Rejection Event")
# #50% hpdi
# lines(c(-1300,-441), exp(HPDI(prob = 0.5, m40s$A0b)[[1]] + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# lines(c(-441:0), exp(HPDI(prob = 0.5, m40s$A0b)[[1]] + c(0:441) * HPDI(prob = 0.5, m40s$bR0b)[[1]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# lines(c(0:1100), exp(HPDI(prob = 0.5, m40s$A0b)[[1]] + 441*HPDI(prob = 0.5, m40s$bR0b)[[1]] + HPDI(prob = 0.5, m40s$A1b)[[1]] + c(0:1100) * HPDI(prob = 0.5, m40s$bR1b)[[1]]),
#       type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# lines(c(-1300,-441), exp(HPDI(prob = 0.5, m40s$A0b)[[2]] + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# lines(c(-441:0), exp(HPDI(prob = 0.5, m40s$A0b)[[2]] + c(0:441) * HPDI(prob = 0.5, m40s$bR0b)[[2]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# lines(c(0:1100), exp(HPDI(prob = 0.5, m40s$A0b)[[2]] + 441*HPDI(prob = 0.5, m40s$bR0b)[[2]] + HPDI(prob = 0.5, m40s$A1b)[[2]] + c(0:1100) * HPDI(prob = 0.5, m40s$bR1b)[[2]]),
#       type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# #rejection date
# abline(v = 0, col = 2, lwd = 2)
# text(x = -25, y = 1, labels = "Rejection Date", srt = 90)
# #transplant date
# abline(v = -441, col = 2, lwd = 2)
# text(x = -466, y = 1, labels = "Average Transplant Date", srt = 90)
# #samples from posterior
# for(i in seq(1, length(m40s$A0b), length.out = 2000)){
#   lines(c(-1300,-441), exp(m40s$A0b[i] + c(0,0)), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
#   lines(c(-441:0), exp(m40s$A0b[i] + c(0:441) * m40s$bR0b[i]), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
#   lines(c(0:1100), exp(m40s$A0b[i] + 441*m40s$bR0b[i] + m40s$A1b[i] + c(0:1100) * m40s$bR1b[i]), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
# }
# 
# 
# #sqrt model, m41
# #separate slopes for before and after the rejection event, seperate "intercepts" at the rejection event, intercept with zero slope before the transplant event,
# #express slope in change per unit sqrt(day) between transplant and rejection event, express post-rejection intercept as displacement from pre-rejection intercept
# load("full/m41full")
# library(RColorBrewer)
# patColors <- unique(as.character(d$patient))
# patColors <- cbind(patColors, brewer.pal(length(patColors), "Dark2"))
# patientColors <- sapply(1:length(d$patient), function(x) patColors[which(d$patient[x] == patColors[,1]),2])
# #postcheck(fit = m31, prob = 0.95, window =  300, col = patientColors)
# 
# m41s <- extract.samples(m41)
# #mean estimate
# par(mar = c(4,4,2,2))
# plot(c(-1300,-441), exp(mean(m41s$A0b) + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,1.5),
#      xlab = "days relative to rejection event", ylab = "rate of spectral count for the average protein/patient", lwd = 2)
# lines(c(-441:0), exp(mean(m41s$A0b) + sqrt(c(0:441)) * mean(m41s$bR0b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
# lines(c(0:1100), exp(mean(m41s$A1b) + mean(m41s$A0b) + sqrt(441)*mean(m41s$bR0b) + sqrt(c(0:1100)) * mean(m41s$bR1b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
# title("Expected Rates of Spectral Count Given Days Relative to Acute Rejection Event")
# #50% hpdi
# lines(c(-1300,-441), exp(HPDI(prob = 0.5, m41s$A0b)[[1]] + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# lines(c(-441:0), exp(HPDI(prob = 0.5, m41s$A0b)[[1]] + sqrt(c(0:441)) * HPDI(prob = 0.5, m41s$bR0b)[[1]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# lines(c(0:1100), exp(HPDI(prob = 0.5, m41s$A0b)[[1]] + sqrt(441)*HPDI(prob = 0.5, m41s$bR0b)[[1]] + HPDI(prob = 0.5, m41s$A1b)[[1]] + sqrt(c(0:1100)) * HPDI(prob = 0.5, m41s$bR1b)[[1]]),
#       type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# lines(c(-1300,-441), exp(HPDI(prob = 0.5, m41s$A0b)[[2]] + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# lines(c(-441:0), exp(HPDI(prob = 0.5, m41s$A0b)[[2]] + sqrt(c(0:441)) * HPDI(prob = 0.5, m41s$bR0b)[[2]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# lines(c(0:1100), exp(HPDI(prob = 0.5, m41s$A0b)[[2]] + sqrt(441)*HPDI(prob = 0.5, m41s$bR0b)[[2]] + HPDI(prob = 0.5, m41s$A1b)[[2]] + sqrt(c(0:1100)) * HPDI(prob = 0.5, m41s$bR1b)[[2]]),
#       type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# #rejection date
# abline(v = 0, col = 2, lwd = 2)
# text(x = -25, y = 1, labels = "Rejection Date", srt = 90)
# #transplant date
# abline(v = -441, col = 2, lwd = 2)
# text(x = -466, y = 1, labels = "Average Transplant Date", srt = 90)
# #samples from posterior
# for(i in seq(1, length(m41s$A0b), length.out = 2000)){
#   lines(c(-1300,-441), exp(m41s$A0b[i] + c(0,0)), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
#   lines(c(-441:0), exp(m41s$A0b[i] + sqrt(c(0:441)) * m41s$bR0b[i]), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
#   lines(c(0:1100), exp(m41s$A0b[i] + sqrt(441)*m41s$bR0b[i] + m41s$A1b[i] + sqrt(c(0:1100)) * m41s$bR1b[i]), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
# }
# 
# #m28, log relationships
# load("full/m28full")
# m28s <- extract.samples(m28)
# #mean estimate
# par(mar = c(4,4,2,2))
# plot(c(-1300,-441), exp(mean(m28s$A0b) + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,3),
#      xlab = "days relative to rejection event", ylab = "rate of spectral count for the average protein/patient", lwd = 2)
# lines(c(-440:0), exp(mean(m28s$A0b) + log(c(1:441)) * mean(m28s$bR0b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
# lines(c(1:1100), exp(mean(m28s$A1b) + mean(m28s$A0b) + log(441)*mean(m28s$bR0b) + log(c(1:1100)) * mean(m28s$bR1b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
# title("Expected Rates of Spectral Count Given Days Relative to Acute Rejection Event")
# #50% hpdi
# lines(c(-1300,-441), exp(HPDI(prob = 0.5, m28s$A0b)[[1]] + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# lines(c(-441:0), exp(HPDI(prob = 0.5, m28s$A0b)[[1]] + log(c(0:441)) * HPDI(prob = 0.5, m28s$bR0b)[[1]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# lines(c(0:1100), exp(HPDI(prob = 0.5, m28s$A0b)[[1]] + log(441)*HPDI(prob = 0.5, m28s$bR0b)[[1]] + HPDI(prob = 0.5, m28s$A1b)[[1]] + log(c(0:1100)) * HPDI(prob = 0.5, m28s$bR1b)[[1]]),
#       type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# lines(c(-1300,-441), exp(HPDI(prob = 0.5, m28s$A0b)[[2]] + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# lines(c(-440:0), exp(HPDI(prob = 0.5, m28s$A0b)[[2]] + log(c(1:441)) * HPDI(prob = 0.5, m28s$bR0b)[[2]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# lines(c(1:1100), exp(HPDI(prob = 0.5, m28s$A0b)[[2]] + log(441)*HPDI(prob = 0.5, m28s$bR0b)[[2]] + HPDI(prob = 0.5, m28s$A1b)[[2]] + log(c(1:1100)) * HPDI(prob = 0.5, m28s$bR1b)[[2]]),
#       type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
# #rejection date
# abline(v = 0, col = 2, lwd = 2)
# text(x = -25, y = 2, labels = "Rejection Date", srt = 90)
# #transplant date
# abline(v = -441, col = 2, lwd = 2)
# text(x = -466, y = 2, labels = "Average Transplant Date", srt = 90)
# #samples from posterior
# for(i in seq(1, length(m28s$A0b), length.out = 2000)){
#   lines(c(-1300,-441), exp(m28s$A0b[i] + c(0,0)), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
#   lines(c(-440:0), exp(m28s$A0b[i] + log(c(1:441)) * m28s$bR0b[i]), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
#   lines(c(1:1100), exp(m28s$A0b[i] + log(441)*m28s$bR0b[i] + m28s$A1b[i] + log(c(1:1100)) * m28s$bR1b[i]), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
# }
# 
# #comparing model fit
# fitmodels <- c("m28, m31", "m40", "m41")
# eval(parse(text=paste0("compare(", paste0(paste0(fitmodels[-length(fitmodels)], ",", collapse = ""), fitmodels[length(fitmodels)]), ")")))
# 
# ##############################################################################################################################
# ########################################### generating compound figure ##########################################
# ##############################################################################################################################
load("full/m28full")
load("full/m31full")
load("full/m40full")
load("full/m41full")

m28s <- extract.samples(m28)
m31s <- extract.samples(m31)
m40s <- extract.samples(m40)
m41s <- extract.samples(m41)

m28p <- link(m28, n = 2000)
m31p <- link(m31, n = 2000)
m40p <- link(m40, n = 2000)
m41p <- link(m41, n = 2000)

rm("m28")
rm("m31")
rm("m40")
rm("m41")

posteriorDraws <- F
modelText <- F
png(filename = "candidateModels.png", width = 1000, height = 1000)
cex.lab.plots = 1.5 #axis label size
redline.label.size = 1.5; redline.label.horizontal.offset = -20
title.size = 3
WAIC.size = 1.2
margins <- rep(4.3,4)
margins <- c(4.3,4.3,4.3,1)
#model 1
par(mar = margins, mfrow = c(2,2))
plot(c(-1300,-441), exp(mean(m31s$A0b) + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,1.5), cex.lab = cex.lab.plots,
     xlab = "days relative to rejection event", ylab = "rate of spectral count for the average protein/patient", lwd = 2)
lines(c(-441,0), exp(mean(m31s$A0b) + mean(m31s$A1b) + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
lines(c(0,1100), exp(mean(m31s$A0b) + mean(m31s$A1b) + mean(m31s$A2b) + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
title("Model 1", cex.main = title.size); text(paste0("WAIC = ", round(WAIC(m31)[1])), x = 900, y = 1.5, cex = WAIC.size)
#50% hpdi
lines(c(-1300,-441), exp(HPDI(prob = 0.5, m31s$A0b)[[1]] + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(-441,0), exp(HPDI(prob = 0.5, m31s$A0b + m31s$A1b)[[1]]) + c(0,0), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(0,1100), exp(HPDI(prob = 0.5, m31s$A0b + m31s$A1b + m31s$A2b)[[1]]) + c(0,0), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(-1300,-441), exp(HPDI(prob = 0.5, m31s$A0b)[[2]]) + c(0,0), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(-441,0), exp(HPDI(prob = 0.5, m31s$A0b + m31s$A1b)[[2]]) + c(0,0), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(0,1100), exp(HPDI(prob = 0.5, m31s$A0b + m31s$A1b + m31s$A2b)[[2]]) + c(0,0), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
#rejection date
abline(v = 0, col = 2, lwd = 2)
text(x = -25 + redline.label.horizontal.offset, y = 1, labels = "Rejection Date", srt = 90, cex = redline.label.size)
#transplant date
abline(v = -441, col = 2, lwd = 2)
text(x = -466 + redline.label.horizontal.offset, y = 1, labels = "Average Transplant Date", srt = 90, cex = redline.label.size)
#samples from posterior
if(posteriorDraws){
  for(i in seq(1, length(m31s$A0b), length.out = 2000)){
    lines(c(-1300,-441), exp(m31s$A0b[i] + c(0,0)), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
    lines(c(-441,0), exp(m31s$A0b[i] + m31s$A1b[i] + c(0,0)), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
    lines(c(0,1100), exp(m31s$A0b[i] + m31s$A1b[i]+ m31s$A2b[i] + c(0,0)), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
  }}

#model 2
par(mar = margins)
plot(c(-1300,-441), exp(mean(m40s$A0b) + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,1.5), cex.lab = cex.lab.plots,
     xlab = "days relative to rejection event", ylab = "rate of spectral count for the average protein/patient", lwd = 2)
lines(c(-441:0), exp(mean(m40s$A0b) + c(0:441) * mean(m40s$bR0b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
lines(c(0:1100), exp(mean(m40s$A1b) + mean(m40s$A0b) + 441*mean(m40s$bR0b) + c(0:1100) * mean(m40s$bR1b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
title("Model 2", cex.main = title.size); text(paste0("WAIC = ", round(WAIC(m40)[1])), x = 900, y = 1.5, cex = WAIC.size)
#50% hpdi
lines(c(-1300,-441), exp(HPDI(prob = 0.5, m40s$A0b)[[1]] + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(-441:0), exp(HPDI(prob = 0.5, m40s$A0b)[[1]] + c(0:441) * HPDI(prob = 0.5, m40s$bR0b)[[1]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(0:1100), exp(HPDI(prob = 0.5, m40s$A0b)[[1]] + 441*HPDI(prob = 0.5, m40s$bR0b)[[1]] + HPDI(prob = 0.5, m40s$A1b)[[1]] + c(0:1100) * HPDI(prob = 0.5, m40s$bR1b)[[1]]), 
      type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(-1300,-441), exp(HPDI(prob = 0.5, m40s$A0b)[[2]] + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(-441:0), exp(HPDI(prob = 0.5, m40s$A0b)[[2]] + c(0:441) * HPDI(prob = 0.5, m40s$bR0b)[[2]]), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
lines(c(0:1100), exp(HPDI(prob = 0.5, m40s$A0b)[[2]] + 441*HPDI(prob = 0.5, m40s$bR0b)[[2]] + HPDI(prob = 0.5, m40s$A1b)[[2]] + c(0:1100) * HPDI(prob = 0.5, m40s$bR1b)[[2]]), 
      type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2, lty = 2)
#rejection date
abline(v = 0, col = 2, lwd = 2)
text(x = -25 + redline.label.horizontal.offset, y = 1, labels = "Rejection Date", srt = 90, cex = redline.label.size)
#transplant date
abline(v = -441, col = 2, lwd = 2)
text(x = -466 + redline.label.horizontal.offset, y = 1, labels = "Average Transplant Date", srt = 90, cex = redline.label.size)
#samples from posterior
if(posteriorDraws){
  for(i in seq(1, length(m40s$A0b), length.out = 2000)){
    lines(c(-1300,-441), exp(m40s$A0b[i] + c(0,0)), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
    lines(c(-441:0), exp(m40s$A0b[i] + c(0:441) * m40s$bR0b[i]), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
    lines(c(0:1100), exp(m40s$A0b[i] + 441*m40s$bR0b[i] + m40s$A1b[i] + c(0:1100) * m40s$bR1b[i]), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
  }}

#model 3
par(mar = margins)
plot(c(-1300,-441), exp(mean(m41s$A0b) + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,1.5), cex.lab = cex.lab.plots,
     xlab = "days relative to rejection event", ylab = "rate of spectral count for the average protein/patient", lwd = 2)
lines(c(-441:0), exp(mean(m41s$A0b) + sqrt(c(0:441)) * mean(m41s$bR0b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
lines(c(0:1100), exp(mean(m41s$A1b) + mean(m41s$A0b) + sqrt(441)*mean(m41s$bR0b) + sqrt(c(0:1100)) * mean(m41s$bR1b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
title("Model 3", cex.main = title.size); text(paste0("WAIC = ", round(WAIC(m41)[1])), x = 900, y = 1.5, cex = WAIC.size)
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
text(x = -25 + redline.label.horizontal.offset, y = 1, labels = "Rejection Date", srt = 90, cex = redline.label.size)
#transplant date
abline(v = -441, col = 2, lwd = 2)
text(x = -466 + redline.label.horizontal.offset, y = 1, labels = "Average Transplant Date", srt = 90, cex = redline.label.size)
#samples from posterior
if(posteriorDraws){
  for(i in seq(1, length(m41s$A0b), length.out = 2000)){
    lines(c(-1300,-441), exp(m41s$A0b[i] + c(0,0)), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
    lines(c(-441:0), exp(m41s$A0b[i] + sqrt(c(0:441)) * m41s$bR0b[i]), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
    lines(c(0:1100), exp(m41s$A0b[i] + sqrt(441)*m41s$bR0b[i] + m41s$A1b[i] + sqrt(c(0:1100)) * m41s$bR1b[i]), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
  }}

#model 4
par(mar = margins)
plot(c(-1300,-441), exp(mean(m28s$A0b) + c(0,0)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,3), cex.lab = cex.lab.plots,
     xlab = "days relative to rejection event", ylab = "rate of spectral count for the average protein/patient", lwd = 2)
lines(c(-440:0), exp(mean(m28s$A0b) + log(c(1:441)) * mean(m28s$bR0b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
lines(c(1:1100), exp(mean(m28s$A1b) + mean(m28s$A0b) + log(441)*mean(m28s$bR0b) + log(c(1:1100)) * mean(m28s$bR1b)), type = "l", col = rgb(0,0,0,1), xlim = c(-1300, 1100), ylim = c(0,10), lwd = 2)
title("Model 4", cex.main = title.size); text(paste0("WAIC = ", round(WAIC(m28)[1])), x = 900, y = 3, cex = WAIC.size)
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
text(x = -25 + redline.label.horizontal.offset, y = 2, labels = "Rejection Date", srt = 90, cex = redline.label.size)
#transplant date
abline(v = -441, col = 2, lwd = 2)
text(x = -466 + redline.label.horizontal.offset, y = 2, labels = "Average Transplant Date", srt = 90, cex = redline.label.size)
#samples from posterior
if(posteriorDraws){
  for(i in seq(1, length(m28s$A0b), length.out = 2000)){
    lines(c(-1300,-441), exp(m28s$A0b[i] + c(0,0)), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
    lines(c(-440:0), exp(m28s$A0b[i] + log(c(1:441)) * m28s$bR0b[i]), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
    lines(c(1:1100), exp(m28s$A0b[i] + log(441)*m28s$bR0b[i] + m28s$A1b[i] + log(c(1:1100)) * m28s$bR1b[i]), type = "l", col = rgb(0,0,0,0.05), xlim = c(-1300, 1100), ylim = c(0,10))
  }}


if(modelText){
  library(latex2exp)
  #text for each model - m31, 40, 41, 21
  #m31
  m31t <-
    "spectral_count ~ poisson(lambda),
  log(lambda) <- a,
  
  a <- A0b + A0j[protein_index] + A0i[patient_index] + a1*afterTransplant + a2*afterRejection,
  A0b ~ normal(0,4),
  A0j[protein_index] ~ normal(0,sigA0j),
  A0i[patient_index] ~ normal(0,sigA0i),
  sigA0j, sigA0i ~ half-cauchy(0,1),
  
  a1 <- A1b + A1j[protein_index] + A1i[patient_index],
  A1b ~ normal(0,4),
  A1j[protein_index] ~ normal(0,sigA1j),
  A1i[patient_index] ~ normal(0,sigA1i),
  sigA1j, sigA1i ~ half-cauchy(0,1),
  
  a2 <- A2b + A2j[protein_index] + A2i[patient_index],
  A2b ~ normal(0,4),
  A2j[protein_index] ~ normal(0,sigA2j),
  A2i[patient_index] ~ normal(0,sigA2i),
  sigA2j, sigA2i ~ half-cauchy(0,1)"
  
  par(mar = c(0,0,0,0))
  plot.new()
  text(m31t, x = 0, y= 0.55, pos = 4, cex = 1.2)
  
  m40t <- 
    "spectral_count ~ poisson(lambda),
  log(lambda) <- a0 
  + bR0*daysRelativeToTransplant*(1-afterRejection)*afterTransplant 
  + bR1*daysRelativeToRejection*afterRejection, 
  
  a0 <- A0b + A0j[protein_index] + A0i[patient_index] + a1*afterRejection,
  A0b ~ normal(0,1),
  A0j[protein_index] ~ normal(0,sigA0j),
  A0i[patient_index] ~ normal(0,sigA0i),
  sigA0j, sigA0i ~ half-cauchy(0,0.25),
  
  a1 <- bR0*daysBetwTrans_Rej + A1b + A1j[protein_index] + A1i[patient_index],
  A1b ~ normal(0,1),
  A1j[protein_index] ~ normal(0,sigA1j),
  A1i[patient_index] ~ normal(0,sigA1i),
  sigA1j, sigA1i ~ half-cauchy(0,0.25),
  
  bR0 <- bR0b + bR0j[protein_index] + bR0i[patient_index],
  bR0j[protein_index] ~ normal(0,sigbR0j),
  bR0i[patient_index] ~ normal(0,sigbR0i),
  bR0b ~ normal(0,0.01),
  sigbR0j, sigbR0i ~ half-cauchy(0,0.001),
  
  bR1 <- bR1b + bR1j[protein_index] + bR1i[patient_index],
  bR1j[protein_index] ~ normal(0,sigbR1j),
  bR1i[patient_index] ~ normal(0,sigbR1i),
  bR1b ~ normal(0,0.01),
  sigbR1j, sigbR1i ~ half-cauchy(0,0.001)"
  
  plot.new()
  text(m40t, x = 0, y= 0.4, pos = 4, cex = 1.2)
  
  m41t <-
    "spectral_count ~ poisson(lambda),
  log(lambda) <- a0 
  + bR0*sqrt(daysRelativeToTransplant)*(1-afterRejection)*afterTransplant 
  + bR1*sqrt(daysRelativeToRejection)*afterRejection, 
  
  a0 <- A0b + A0j[protein_index] + A0i[patient_index] + a1*afterRejection,
  A0b ~ normal(0,1),
  A0j[protein_index] ~ normal(0,sigA0j),
  A0i[patient_index] ~ normal(0,sigA0i),
  sigA0j, sigA0i ~ half-cauchy(0,0.25),
  
  a1 <- bR0*sqrt(daysBetwTrans_Rej) + A1b + A1j[protein_index] + A1i[patient_index],
  A1b ~ normal(0,1),
  A1j[protein_index] ~ normal(0,sigA1j),
  A1i[patient_index] ~ normal(0,sigA1i),
  sigA1j, sigA1i ~ half-cauchy(0,0.25),
  
  bR0 <- bR0b + bR0j[protein_index] + bR0i[patient_index],
  bR0j[protein_index] ~ normal(0,sigbR0j),
  bR0i[patient_index] ~ normal(0,sigbR0i),
  bR0b ~ normal(0,0.1),
  sigbR0j, sigbR0i ~ half-cauchy(0,0.01),
  
  bR1 <- bR1b + bR1j[protein_index] + bR1i[patient_index],
  bR1j[protein_index] ~ normal(0,sigbR1j),
  bR1i[patient_index] ~ normal(0,sigbR1i),
  bR1b ~ normal(0,0.1),
  sigbR1j, sigbR1i ~ half-cauchy(0,0.01)"
  
  plot.new()
  text(m41t, x = 0, y= 0.4, pos = 4, cex = 1.2)
  
  m28t <- 
    "spectral_count ~ poisson(lambda),
  log(lambda) <- a0 
  + bR0*log(daysRelativeToTransplant)*(1-afterRejection)*afterTransplant 
  + bR1*log(daysRelativeToRejection)*afterRejection, 
  
  a0 <- A0b + A0j[protein_index] + A0i[patient_index] + a1*afterRejection,
  A0b ~ normal(0,4),
  A0j[protein_index] ~ normal(0,sigA0j),
  A0i[patient_index] ~ normal(0,sigA0i),
  sigA0j, sigA0i ~ half-cauchy(0,1),
  
  a1 <- bR0*log(daysBetwTrans_Rej) + A1b + A1j[protein_index] + A1i[patient_index],
  A1b ~ normal(0,4),
  A1j[protein_index] ~ normal(0,sigA1j),
  A1i[patient_index] ~ normal(0,sigA1i),
  sigA1j, sigA1i ~ half-cauchy(0,1),
  
  bR0 <- bR0b + bR0j[protein_index] + bR0i[patient_index],
  bR0j[protein_index] ~ normal(0,sigbR0j),
  bR0i[patient_index] ~ normal(0,sigbR0i),
  bR0b ~ normal(0,2),
  sigbR0j, sigbR0i ~ half-cauchy(0,1),
  
  bR1 <- bR1b + bR1j[protein_index] + bR1i[patient_index],
  bR1j[protein_index] ~ normal(0,sigbR1j),
  bR1i[patient_index] ~ normal(0,sigbR1i),
  bR1b ~ normal(0,2),
  sigbR1j, sigbR1i ~ half-cauchy(0,1)
  "
  plot.new()
  text(m28t, x = 0, y= 0.4, pos = 4, cex = 1.2)
  
  
  # #plot pearson residuals -- sampled lambda
  # numSamps <- 400
  # indices <- floor(seq(from = 1, to = length(m41p$lambda[,1]), length.out = numSamps))
  # resids <- rep(0, length(d$spectral_count)*length(indices))
  # for(x in 1:length(d$spectral_count)){
  #   if (x %% 100 == 0) {cat(paste(x, " "))}
  #   resids[((x-1)*numSamps+1):(x*numSamps)] <- (m41p$lambda[,x][indices] - d$spectral_count[x]) / sqrt(m41p$lambda[,x][indices])
  # }
  # hist(resids, breaks = 500, xlim = c(-5, 5), freq = F)
  # lines(seq(-5, 5, length.out = 200), dnorm(x = seq(-5, 5, length.out = 200), mean = 0, sd = 1), lwd = 2, col = 2)
  # 
  # #plot pearson residuals -- sampled lambda, nonzero counts only
  # numSamps <- 400
  # indices <- floor(seq(from = 1, to = length(m41p$lambda[,1]), length.out = numSamps))
  # resids <- rep(0, length(d$spectral_count)*length(indices))
  # for(x in 1:length(d$spectral_count)){
  #   if (x %% 100 == 0) {cat(paste(x, " "))}
  #   if(d$spectral_count[x] == 0) {resids[((x-1)*numSamps+1):(x*numSamps)] <- NA} else {
  #     resids[((x-1)*numSamps+1):(x*numSamps)] <- (m41p$lambda[,x][indices] - d$spectral_count[x]) / sqrt(m41p$lambda[,x][indices])}
  # }
  # hist(resids, breaks = 500, xlim = c(-5, 5), ylim = c(0,0.4), freq = F)
  # lines(seq(-5, 5, length.out = 200), dnorm(x = seq(-5, 5, length.out = 200), mean = 0, sd = 1), lwd = 2, col = 2)
}
dev.off()


#evaluating multiple chain diagnostics
load("full/m28full")
load("full/m28c2full")
sflist <- c(m28@stanfit, m28c2@stanfit)
m28comb <- sflist2stanfit(sflist); rm(sflist); rm(m28); rm(m28c2)
summary(m28comb)
min(summary(m28comb)$summary[,"n_eff"])
max(summary(m28comb)$summary[,"Rhat"])

load("full/m31full")
load("full/m31c2full")
sflist <- c(m31@stanfit, m31c2@stanfit)
m31comb <- sflist2stanfit(sflist); rm(sflist); rm(m31); rm(m31c2)
summary(m31comb)
min(summary(m31comb)$summary[,"n_eff"])
max(summary(m31comb)$summary[,"Rhat"])

load("full/m40c2full"); m40c2 <- m40; rm(m40)
load("full/m40full")
sflist <- c(m40@stanfit, m40c2@stanfit)
m40comb <- sflist2stanfit(sflist); rm(sflist); rm(m40); rm(m40c2)
summary(m40comb)
min(summary(m40comb)$summary[,"n_eff"])
max(summary(m40comb)$summary[,"Rhat"])

load("full/m41c2full"); m41c2 <- m41; rm(m41)
load("full/m41full")
sflist <- c(m41@stanfit, m41c2@stanfit)
m41comb <- sflist2stanfit(sflist); rm(sflist); rm(m41); rm(m41c2)
summary(m41comb)
min(summary(m41comb)$summary[,"n_eff"])
max(summary(m41comb)$summary[,"Rhat"])

m41s <- extract.samples(m41comb, clean.names = F)
