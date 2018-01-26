################################################################################
# Author: Jorge Gomez (j.gomez.robles@student.tue.nl)
################################################################################

################################################################################
################################  SET WORKSPACE ################################
################################################################################
currentDir <- getSrcDirectory(function(x) {x})
setwd(currentDir)

################################################################################
############################ LOAD REQUIRED LIBRARIES ###########################
################################################################################
library(stringr)
library(dplyr)
library(tidyr)

################################################################################
############################# SET GLOBAL VARIABLES #############################
################################################################################

# Assume we are working in the same directory where the script resides in
DIR_RELATIVE_PATH <- 'datasets/genderstats'

filesIdsS <- c("3", "4", "5", "6", "8a", "8b", "8c", "10", "20")

filesIdsNS <- c("24a", "24b", "24c", "25")


bySexSet <- NULL
noSexSet <- NULL

################################################################################
############################# DEFINE MAIN FUNCTION #############################
################################################################################
processFilesWithBySex <- function(){
  for(fileId in filesIdsS){
       path <- str_c(DIR_RELATIVE_PATH, "/", fileId, ".csv")  
       currentFile <- tbl_df( read.csv(path, stringsAsFactors = FALSE))
       names(currentFile)[1] <- "Indicator"
       names(currentFile)[3] <- "CountryCode"
       cleanFile <- select(currentFile, Indicator, Region, CountryCode, Country, Year, Sex, Value)
       #cleanFile <- cleanFile %>%
       #cleanFile <- group_by(cleanFile, Indicator, Country, Year)
       #cleanFile <- aggregate(state.x77, list(Region = state.region), mean)
       #  mutate(prop = Value[Sex=="Male"] )
        #mutate(prop = 1)
       #aggregate(cleanFile, list(as.factor(cleanFile$Indicator), as.factor(cleanFile$Country), as.factor(cleanFile$Indicator)), FUN= mean)
       #cleanFile <- aggregate(. ~ Indicator + Country, cleanFile, mean)
       # Combining all data
       # View(cleanFile, fileId)
       if( is.null(bySexSet) ){
         bySexSet <<- cleanFile
       }else{
         bySexSet <<- bind_rows(bySexSet, cleanFile)
       }
  }
  writeToCsv(bySexSet, "GenderStatisticsBySex")
}

################################################################################
############################# DEFINE MAIN FUNCTION #############################
################################################################################
processFilesWithoutBySex <- function(){
  for(fileId in filesIdsNS){
    path <- str_c(DIR_RELATIVE_PATH, "/", fileId, ".csv")  
    arrayId <- str_c("f", fileId)  
    currentFile <- tbl_df( read.csv(path, stringsAsFactors = FALSE))
    names(currentFile)[1] <- "Indicator"
    names(currentFile)[3] <- "CountryCode"
    cleanFile <- select(currentFile, Indicator, Region, CountryCode, Country, Year, Sex, Value)
    
    # Combining all data
    # View(cleanFile, fileId)
    if( is.null(noSexSet) ){
      noSexSet <<- cleanFile
    }else{
      noSexSet <<- bind_rows(noSexSet, cleanFile)
    }
  }
  writeToCsv(noSexSet, "GenderStatisticsFemaleOnly")
}


# Write the final dataset
writeToCsv <- function(dataset, fileName) {
  path <- str_c(DIR_RELATIVE_PATH, "/", fileName, ".csv")
  write.csv(dataset, path, row.names = FALSE)
}
