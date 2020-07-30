sampleIterations <- 20000; warmupIterations <- 10000; nChains <- 1
for(i in 1:1){ #hacky way of loading everything
  
  library(rethinking)
  thinDataDownForExploration <- T; min_n_pos <- 30
  keratinIndex <- 1
  keratin <- c("removeCorrectly", "removeIncorrectly")[keratinIndex]
  dataSizeIndex <- 1
  dataSize <- c("full", "default", "sparse", "defaultfixed")[dataSizeIndex]
  addNewPatients <- T
  EDA1 <- F
  EDA2 <- T
  exploreXFOLDchange <- F
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
  d_new_old <- read.csv("newData.csv", header = T)
  d_new <- read.csv("New-Data-Combined-KateFebruaryCorrections.csv", header = T)
  d_new <- read.csv("newData-H-2020Feb29.csv", header = T)
  n_timepoints <- 19
  p_new <- data.frame(matrix(NA, nrow = n_timepoints*nrow(d_new), ncol = length(columnNames))); colnames(p_new) <- columnNames
  p_new$spectral_count <- c(sapply(7:25, function(x) d_new[,x]))
  p_new$patient <- c(sapply(7:25, function(x) rep(substr(colnames(d_new)[x],1,2), length(d_new[,x]))))
  p_new$protein <- rep(trimws(d_new$Identified.Proteins), n_timepoints)
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
  
  #force overlap between old and new proteins
  library(stringr)
  d$protein <- as.character(d$protein)
  d$protein <- tolower(d$protein)
  d$protein <- str_replace(d$protein, pattern = " \\[homo sapiens\\]", replacement = "")
  d$protein <- str_replace(d$protein, pattern = " precursor", replacement = "")
  d$protein <- str_replace(d$protein, pattern = " preproprotein", replacement = "")
  d$protein <- sapply(1:length(d$protein), function(x) strsplit(d$protein[x], split = " os=")[[1]][1])
  # d$protein <- sapply(1:length(d$protein), function(x) strsplit(d$protein[x], split = " isoform")[[1]][1])
  d[as.numeric(strsplit("574 1183 1577 1954 2082 3151 3219 3815 3854 3942 4380 5448 7120 9043", split = " ")[[1]]),]
  d <- d[!is.na(d$protein),]
  d <- d[-grep(pattern = "contaminant", x = d$protein),]
  d$protein <- trimws(d$protein, which = "both")
  
  #filter out keratin and "null" protein
  if(keratin == "removeCorrectly"){
    d <- d[!(grepl(pattern = "keratin", as.character(d$protein)) | grepl(pattern = "Keratin", as.character(d$protein)) | as.character(d$protein) == ""),]
  } else if (keratin == "removeIncorrectly"){
    d <- d[!(grepl(pattern = "keratin", as.character(d$protein)) | as.character(d$protein) == ""),]
  }
  d <- d[!(grepl(pattern = "contaminant", as.character(d$protein)) | grepl(pattern = "reversed", as.character(d$protein)) | grepl(pattern = "immunoglobulin", as.character(d$protein))),]
  
  
  all_prots <- unique(d$protein)
  old_proteins <- unique(d$protein[!is.na(d$date_of_rejection)])
  new_proteins <- unique(d$protein[is.na(d$date_of_rejection)])
  uniq_old <- sort(setdiff(old_proteins, new_proteins))
  uniq_new <- sort(setdiff(new_proteins, old_proteins))
  overlap <- sort(intersect(new_proteins, old_proteins))
  
  write.csv("unique_proteins_newPatients.csv", x = uniq_new)
  write.csv("unique_proteins_oldPatients.csv", x = uniq_old)
  write.csv("overlapping_proteins.csv", x = overlap)
  
  

  
  #check overlap between old and new proteins
  
  # old_proteins <- as.character(unique(rbind(p1, p3, p5, p6, p7, p9)$protein))
  # new_proteins <- unique(p_new$protein)
  # new_proteins <- new_proteins[-grep(pattern = "CONTAMINANT", x = new_proteins)]
  # length(old_proteins) + length(new_proteins)
  # length(unique(c(old_proteins, new_proteins)))
  # new_proteins_new_names <- trimws(read.csv("newData_2.csv", header = T)$Identified.Proteins, which = "both")
  # 
  # length(old_proteins) + length(new_proteins_new_names)
  # length(unique(c(old_proteins, new_proteins_new_names)))
  # 
  # length(new_proteins) + length(new_proteins_new_names)
  # length(unique(c(new_proteins, new_proteins_new_names)))
  # 
  # head(sort(old_proteins))
  # head(sort(new_proteins_new_names))
  # head(sort(new_proteins))
  # 
  # old_proteins[grep(x = old_proteins, pattern = "60 kDa heat shock protein")]
  # new_proteins[grep(x = new_proteins, pattern = "60 kDa heat shock protein")]
  # 
  # old_proteins[grep(x = old_proteins, pattern = "prohibitin")]
  # new_proteins[grep(x = new_proteins, pattern = "proteasome")]
  # 
  # 
  # new_proteins_new_names[grep(x = new_proteins_new_names, pattern = "60 kDa heat shock protein")]
  # 
  
  if(dataSize == "full" || dataSize == "full_removeP9"){
    #populate the data with zero spectral counts where appropriate
    patients <- as.character(unique(d$patient))
    proteins <- as.character(unique(d$protein))
    
    unequalRedundancies <- d[1,]
    
    for(i in 1:length(patients)){
      print(paste0("patient ",i))
      n_tps <- length(unique(d[d$patient == patients[i],]$timepoint))
      for(j in 1:n_tps){ #timepoints
        patient <- patients[i]
        timepoint <- j
        subData <- d[d$patient == patient & d$timepoint == timepoint,]
        proteinsHere <- subData$protein
        
        protCounts <- (table(proteinsHere))
        repeatProteins <- attr(protCounts, which = "dimnames")$proteinsHere[which(as.integer(protCounts) > 1)]
        uneqRepProt <- repeatProteins[(sapply(1:length(repeatProteins), function(prot) (length(unique(subData[subData$protein == repeatProteins[prot],]$spectral_count)))) > 1)]
        unequalRedundancies <- rbind(unequalRedundancies, do.call(rbind, lapply(uneqRepProt, function(prot) subData[subData$protein == prot,])))
        
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
    write.csv(file = "unequalSpectralCounts_sameProteinPatientTimepoint.csv", x = unequalRedundancies[-1,], row.names = F)
  }
  
  if(dataSize == "defaultfixed"){
    #populate the data with zero spectral counts where appropriate
    patients <- as.character(unique(d$patient))
    
    for(i in 1:length(patients)){
      print(paste0("patient ",i))
      for(j in 1:(if(i==1 | i == 8){2}else{3})){ #timepoints
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
  d$afterRejection[is.na(d$afterRejection)] <- 0
  # d$pivot <- NA
  
  
  #coerce the protein id
  d$protein <- as.factor(d$protein)
  d$protein_id <- coerce_index(d$protein)
  d$patient <- as.factor(d$patient)
  d$patient_id <- coerce_index(d$patient)
  
  #EDA
  if(EDA2 == T){
    library(RColorBrewer)
    
    max(d$spectral_count)
    hist(d$daysRelativeToTransplant)
    plot(c(1,1), xlim = c(-800,3000), ylim = c(3,12), xlab = "days relative to transplant", ylab = "mean nonzero spectral count / timepoint", type = "n")
    par(xpd=TRUE); points(0, 12.55, pch = 25, bg = 2, col = 2); points(c(-14,14), c(12.62,12.62), pch = 16, col = 2, cex = 0.85); par(xpd=F)
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
    text("l", col = 2, x = 1985, y = 12.1, cex = 1.3); text("represent rejection dates", col = 1, x = 2450, y = 12.1)
    title("Vertical Lines Represent Rejection Dates, Non-Red Colors Index Rejecting Patients,\n Numbers Represent Mean # of Unique Proteins Across All Available Timepoints")
    }
    
    if(EDA1 == T){
    
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
    }
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
  #explore x-fold change for Kate
  if(exploreXFOLDchange){
    unique(d$patient)
    
    d$protein <- as.character(d$protein)
    d$patient <- as.character(d$patient)
    str(d)
    old <- d[sapply(1:length(d$patient), function(x) any(d$patient[x] == c("P1", "P3", "P5", "P6", "P7", "P9"))),]
    new <- d[sapply(1:length(d$patient), function(x) any(d$patient[x] == c("X1", "X2", "X3", "X4", "X5", "X6", "X7"))),]
    
    str(new)
    prots <- unique(new$protein)
    pats <- unique(new$patient)
    targ_prots <- prots[sapply(1:length(prots), function(prot) length(unique(new[new$protein == prots[prot] & new$spectral_count != 0,]$patient)) >= 3)]
    targ_prot <- 6
    xfold_2_1 <- t(sapply(1:length(targ_prots), function(targ_prot)
      sapply(1:length(pats), function(pat) new[new$protein == targ_prots[targ_prot] & new$patient == pats[pat] & new$timepoint2,]$spectral_count / 
               new[new$protein == targ_prots[targ_prot] & new$patient == pats[pat] & new$timepoint1,]$spectral_count)))
    dubs21 <- xfold_2_1 >= 2
    dubs21[is.na(dubs21)] <- FALSE
    targ_prots[apply(dubs21, 1, sum) >= 3]
    write.csv(cbind(gsub(x = targ_prots[apply(dubs21, 1, sum) >= 3], pattern = ",", replacement = ""), xfold_2_1[apply(dubs21, 1, sum) >= 3,]), 
              file = paste0("doubling_3orMorePatients_tp1-2_", Sys.Date(), ".csv"))
    halves21 <- xfold_2_1 <= 0.5
    halves21[is.na(halves21)] <- FALSE
    targ_prots[apply(halves21, 1, sum) >= 3]
    write.csv(cbind(gsub(x = targ_prots[apply(halves21, 1, sum) >= 3], pattern = ",", replacement = ""), xfold_2_1[apply(halves21, 1, sum) >= 3,]), 
              file = paste0("halving_3orMorePatients_tp1-2_", Sys.Date(), ".csv"))
    
    xfold_3_1 <- t(sapply(1:length(targ_prots), function(targ_prot)
      sapply(1:length(pats), function(pat) new[new$protein == targ_prots[targ_prot] & new$patient == pats[pat] & new$timepoint3,]$spectral_count / 
               new[new$protein == targ_prots[targ_prot] & new$patient == pats[pat] & new$timepoint1,]$spectral_count)))
    
    dubs31 <- xfold_3_1 >= 2
    dubs31[is.na(dubs31)] <- FALSE
    targ_prots[apply(dubs31, 1, sum) >= 3]
    write.csv(cbind(gsub(x = targ_prots[apply(dubs31, 1, sum) >= 3], pattern = ",", replacement = ""), xfold_3_1[apply(dubs31, 1, sum) >= 3,]), 
              file = paste0("doubling_3orMorePatients_tp1-3_", Sys.Date(), ".csv"))
    halves31 <- xfold_3_1 <= 0.5
    halves31[is.na(halves31)] <- FALSE
    targ_prots[apply(halves31, 1, sum) >= 3]
    write.csv(cbind(gsub(x = targ_prots[apply(halves31, 1, sum) >= 3], pattern = ",", replacement = ""), xfold_3_1[apply(halves31, 1, sum) >= 3,]), 
              file = paste0("halving_3orMorePatients_tp1-3_", Sys.Date(), ".csv"))
    
    plot(new$daysRelativeToTransplant, new$spectral_count, type = "l")
  }
  if(thinDataDownForExploration){
    patients <- as.character(unique(d$patient))
    proteins <- as.character(unique(d$protein))
    d_thin <- d[1,]
    for(i in 1:length(patients)){
      n_tps <- length(unique(d[d$patient == patients[i],]$timepoint))
      for(j in 1:n_tps){ #timepoints
        patient <- patients[i]
        timepoint <- j
        subData <- d[d$patient == patient & d$timepoint == timepoint,]
        sd_zeros <- subData[subData$spectral_count == 0,]
        sd_nonzeros <- subData[subData$spectral_count > 0,]
        d_thin <- rbind(d_thin, sd_nonzeros[sample(1:nrow(sd_nonzeros), min_n_pos, replace = F),])
      }
    }
    d_thin <- d_thin[-1,]
    thin_prots_incl <- unique(d_thin$protein)
    for(i in 1:length(patients)){
      n_tps <- length(unique(d[d$patient == patients[i],]$timepoint))
      for(j in 1:n_tps){ #timepoints
        patient <- patients[i]
        timepoint <- j
        subData <- d[d$patient == patient & d$timepoint == timepoint,]
        subData_thin <- d_thin[d_thin$patient == patient & d_thin$timepoint == timepoint,]
        notHereYet <- setdiff(thin_prots_incl, unique(subData_thin$protein))
        toAdd_inds <- sapply(1:length(notHereYet), function(prot) which(subData$protein == notHereYet[prot]))
        d_thin <- rbind(d_thin, subData[toAdd_inds,])
      }
    }
    d_thin$protein <- as.factor(as.character(d_thin$protein))
    d_thin$protein_id <- coerce_index(d_thin$protein)
    d_thin$patient <- as.factor(as.character(d_thin$patient))
    d_thin$patient_id <- coerce_index(d_thin$patient)
    
  }
}

#use persistent thinned dataset
# save(d_thin, file = "d_thin")
load("d_thin")

# model has just intercept and one coefficient corresponding to an indicator variable signifying occurence of discrete rejection event
# d_sub <- d_thin[c("spectral_count", "daysRelativeToTransplant", "afterRejection", "afterTransplant", 
#                   "daysRelativeToRejection", "protein_id", "patient_id", "daysBetwTrans_Rej")]
# ds <- d_thin[c("spectral_count", "patient_id", "afterRejection")]
# m1 <- map2stan(
#   alist(
#     spectral_count ~ dpois( lambda ),
#     log(lambda) <- a + bR*afterRejection,
#     a ~ dnorm(0,100),
#     bR ~ dnorm(0,1)
#   ) ,
#   data= ds,
#   iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = 4,
#   start = list(a=0, bR=0)
# )
# save(file = paste0(dataSize, "/", "m1", dataSize), m1)
# precis(m1)

nChains <- 1

d_thin$daysBetwTrans_Rej[d_thin$rejection == 0] <- 1
d_thin$daysRelativeToRejection[d_thin$rejection == 0] <- 1
d_sub <- d_thin[c("spectral_count", "daysRelativeToTransplant", "afterRejection", "afterTransplant", 
             "daysRelativeToRejection", "protein_id", "patient_id", "daysBetwTrans_Rej", "rejection")]
m2 <- map2stan(
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
  data= d_sub,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains
  # ,start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
save(file = paste0(dataSize, "/", "m2", dataSize), m2)

d_thin$daysBetwTrans_Rej[d_thin$rejection == 0] <- 1
d_thin$daysRelativeToRejection[d_thin$rejection == 0] <- 1
d_sub <- d_thin[c("spectral_count", "daysRelativeToTransplant", "afterRejection", "afterTransplant", 
                  "daysRelativeToRejection", "protein_id", "patient_id", "daysBetwTrans_Rej", "rejection")]
m3 <- map2stan(
  alist(
    spectral_count ~ dgampois( lambda , phi),
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
    c(sigbR1j, sigbR1i) ~ dcauchy(0,0.01),
    
    phi ~ dexp(0.5)
    
  ) ,
  data= d_sub,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains, init_r = 0.1
  # ,start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
precis(m3)

d_thin$daysBetwTrans_Rej[d_thin$rejection == 0] <- 1
d_thin$daysRelativeToRejection[d_thin$rejection == 0] <- 1
d_sub <- d_thin[c("spectral_count", "daysRelativeToTransplant", "afterRejection", "afterTransplant", 
                  "daysRelativeToRejection", "protein_id", "patient_id", "daysBetwTrans_Rej", "rejection")]

m4 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda),
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
  data= d_sub,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains, init_r = 0.1
  # ,start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
precis(m4)

#d_sub[(d_sub$rejection == 0),]

## m5 utterly failed to sample
d_thin$daysBetwTrans_Rej[d_thin$rejection == 0] <- 1
d_thin$daysRelativeToRejection[d_thin$rejection == 0] <- 1
d_sub <- d_thin[c("spectral_count", "daysRelativeToTransplant", "afterRejection", "afterTransplant", 
                  "daysRelativeToRejection", "protein_id", "patient_id", "daysBetwTrans_Rej", "rejection")]
m5 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda),
    log(lambda) <- a0 + bR0*daysRelativeToTransplant*(1-afterRejection)*afterTransplant + bR1*daysRelativeToRejection*afterRejection, 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id] + a1*afterRejection + a_rej*rejection,
    c(A0b,a_rej) ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,1),
    
    a1 <- bR0*daysBetwTrans_Rej + A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id] + b_rej0*rejection,
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    c(bR0b, b_rej0) ~ dnorm(0,0.03),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,0.01),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,0.03),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,0.01)
    
    
  ) ,
  data= d_sub,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains, init_r = 0.1
  # ,start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
precis(m5)

d_thin$daysBetwTrans_Rej[d_thin$rejection == 0] <- 1
d_thin$daysRelativeToRejection[d_thin$rejection == 0] <- 1
d_sub <- d_thin[c("spectral_count", "daysRelativeToTransplant", "afterRejection", "afterTransplant", 
                  "daysRelativeToRejection", "protein_id", "patient_id", "daysBetwTrans_Rej", "rejection")]
m6 <- map2stan(
  alist(
    spectral_count ~ dgampois( lambda , phi),
    log(lambda) <- a0 + bR0*daysRelativeToTransplant*(1-afterRejection)*afterTransplant + bR1*daysRelativeToRejection*afterRejection, 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id] + a1*afterRejection + a_rej,
    c(a_rej, A0b) ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dcauchy(0,1),
    
    a1 <- bR0*daysBetwTrans_Rej + A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dcauchy(0,1),
    
    bR0 <- bR0b + bR0j[protein_id] + bR0i[patient_id] + b_rej,
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    c(bR0b, b_rej) ~ dnorm(0,0.03),
    c(sigbR0j, sigbR0i) ~ dcauchy(0,0.01),
    
    bR1 <- bR1b + bR1j[protein_id] + bR1i[patient_id],
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,0.03),
    c(sigbR1j, sigbR1i) ~ dcauchy(0,0.01),
    
    phi ~ dexp(0.5)
    
  ) ,
  data= d_sub,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains, init_r = 0.1
  # ,start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
save(x = m6, file = "m6_centered")
precis(m6)

#let's try reparameterizing m6 to see if it improves sampling
#i.e. pull the mean and scale parameters out, so e.g.
# vector[K] beta
#beta ~ normal(mu_beta, sigma_beta)
#mu_beta ~ normal(0,4)
#sigma_beta ~ exponential(1)
## becomes ##
# vector[K] beta
#beta = mu_beta + sigma_beta * beta_raw;
#beta_raw ~ std_normal();
#mu_beta ~ normal(0,4)
#sigma_beta ~ exponential(1)

d_thin$daysBetwTrans_Rej[d_thin$rejection == 0] <- 1
d_thin$daysRelativeToRejection[d_thin$rejection == 0] <- 1
d_sub <- d_thin[c("spectral_count", "daysRelativeToTransplant", "afterRejection", "afterTransplant", 
                  "daysRelativeToRejection", "protein_id", "patient_id", "daysBetwTrans_Rej", "rejection")]

m7 <- map2stan(
  alist(
    spectral_count ~ dpois( lambda ),
    log(lambda) <- a0 + bR0*daysRelativeToTransplant*(1-afterRejection)*afterTransplant + bR1*daysRelativeToRejection*afterRejection, 
    
    a0 <- A0b + A0j[protein_id] + A0i[patient_id] + a1*afterRejection + a_rej,
    c(a_rej, A0b) ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,sigA0j),
    A0i[patient_id] ~ dnorm(0,sigA0i),
    c(sigA0j, sigA0i) ~ dexp(1),
    
    a1 <- bR0*daysBetwTrans_Rej + A1b + A1j[protein_id] + A1i[patient_id],
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,sigA1j),
    A1i[patient_id] ~ dnorm(0,sigA1i),
    c(sigA1j, sigA1i) ~ dexp(1),
    
    bR0 <- (bR0b + bR0j[protein_id] + bR0i[patient_id] + b_rej) * 0.03,
    bR0j[protein_id] ~ dnorm(0,sigbR0j),
    bR0i[patient_id] ~ dnorm(0,sigbR0i),
    c(bR0b, b_rej) ~ dnorm(0,1),
    c(sigbR0j, sigbR0i) ~ dexp(1),
    
    bR1 <- (bR1b + bR1j[protein_id] + bR1i[patient_id]) * 0.03,
    bR1j[protein_id] ~ dnorm(0,sigbR1j),
    bR1i[patient_id] ~ dnorm(0,sigbR1i),
    bR1b ~ dnorm(0,1),
    c(sigbR1j, sigbR1i) ~ dexp(1)
    
  ) ,
  data= d_sub,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains, init_r = 0.1
  # ,start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
precis(m7)

#further reparameterizatiton

d_thin$daysBetwTrans_Rej[d_thin$rejection == 0] <- 1
d_thin$daysRelativeToRejection[d_thin$rejection == 0] <- 1
d_sub <- d_thin[c("spectral_count", "daysRelativeToTransplant", "afterRejection", "afterTransplant", 
                  "daysRelativeToRejection", "protein_id", "patient_id", "daysBetwTrans_Rej", "rejection")]
m8 <- ulam(
  alist(
    spectral_count ~ dgampois( lambda , phi),
    log(lambda) <- a0 + bR0*daysRelativeToTransplant*(1-afterRejection)*afterTransplant + bR1*daysRelativeToRejection*afterRejection, 
    
    a0 <- A0b + A0j[protein_id]*sigA0j + A0i[patient_id]*sigA0i + a1*afterRejection + a_rej,
    c(a_rej, A0b) ~ dnorm(0,4),
    A0j[protein_id] ~ dnorm(0,1),
    A0i[patient_id] ~ dnorm(0,1),
    c(sigA0j, sigA0i) ~ dhalfnorm(0,1),
    
    a1 <- bR0*daysBetwTrans_Rej + A1b + A1j[protein_id]*sigA1j + A1i[patient_id]*sigA1i,
    A1b ~ dnorm(0,4),
    A1j[protein_id] ~ dnorm(0,1),
    A1i[patient_id] ~ dnorm(0,1),
    c(sigA1j, sigA1i) ~ dhalfnorm(0,1),
    
    bR0 <- (bR0b*2 + bR0j[protein_id]*sigbR0j + bR0i[patient_id]*sigbR0i + b_rej*2)*0.03,
    bR0j[protein_id] ~ dnorm(0,1),
    bR0i[patient_id] ~ dnorm(0,1),
    c(bR0b, b_rej) ~ dnorm(0,1),
    c(sigbR0j, sigbR0i) ~ dhalfnorm(0,1),
    
    bR1 <- (bR1b*2 + bR1j[protein_id]*sigbR1j + bR1i[patient_id]*sigbR1i)*0.03,
    bR1j[protein_id] ~ dnorm(0,1),
    bR1i[patient_id] ~ dnorm(0,1),
    bR1b ~ dnorm(0,1),
    c(sigbR1j, sigbR1i) ~ dhalfnorm(0,1),
    
    phi ~ dexp(0.5)
    
  ) ,
  data= d_sub,
  iter = sampleIterations, warmup = warmupIterations, chains = nChains, cores = nChains, init_r = 0.1, control = list(adapt_delta=0.98)
  # ,start = list(sigbR0j=0.01, sigbR0i=0.01, bR0b=0, sigbR1j=0.01, sigbR1i=0.01, bR1b=0, A0b=0, sigA0j=1, sigA0i=1, A1b=0, sigA1j=1, sigA1i=1)
)
save(m8, file = "m8_noncentered")
m8s <- extract.samples(m8)
coeftab(m8)
precis(m8)

