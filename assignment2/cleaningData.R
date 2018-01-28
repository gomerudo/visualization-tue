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

filesIdsS <- c("3", "5", "8a", "8b", "8c", "20")

filesIdsNS <- c("24a", "24b", "24c")


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
       
       if( is.null(bySexSet) ){
         bySexSet <<- cleanFile
       }else{
         bySexSet <<- bind_rows(bySexSet, cleanFile)
       }
  }
  getRatio()
  #writeToCsv(bySexSet, "GenderStatisticsBySex")
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
    if( is.null(noSexSet) ){
      noSexSet <<- cleanFile
    }else{
      noSexSet <<- bind_rows(noSexSet, cleanFile)
    }
  }
  getRatioFake()
  #writeToCsv(noSexSet, "GenderStatisticsParity")
}

getRatio <- function(){
  bySexSet <<- bySexSet %>% 
    group_by(Indicator, Country, Year)
  
  bySexSet <<- bySexSet[
    with(bySexSet, order(Indicator, Country, Year, Sex)),
    ]
  bySexSet <<- bySexSet %>% 
      mutate(SexRatio = customProportion(Value))
}

getRatioFake <- function(){
  noSexSet <<- noSexSet %>% 
    mutate(SexRatio = 1)
}


customProportion <- function(dataset){
  return(dataset[2] / dataset[3])
}



# Write the final dataset
writeToCsv <- function(dataset, fileName) {
  path <- str_c(DIR_RELATIVE_PATH, "/", fileName, ".csv")
  write.csv(dataset, path, row.names = FALSE)
}

main <- function(){
  processFilesWithBySex()
  processFilesWithoutBySex()
  jointSet <- bind_rows(bySexSet, noSexSet)
  writeToCsv(jointSet, "GenderStatisticsBySex")
}
