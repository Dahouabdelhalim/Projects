
library(RPostgreSQL)
library(graph)
library(igraph)
library(tidyverse)
library(ggplot2)
library(stringi)
library(stringr)

library(compiler)
# HelperFunctions.R
#SQUELCH = 10
## Loading required package: DBI
pg = dbDriver("PostgreSQL")
con = dbConnect(pg, user="smalec", password="mandarin",
                host="localhost", port=5432, dbname="smdb")


writeVec <- function(dat, fn) {
  write.table(
    unique(dat),
    file = fn,
    quote = FALSE,
    row.names = FALSE,
    eol = "\\n",
    col.names = FALSE,
    append = FALSE 
  )
}

cleanupResultsFile <- function(filename) {
  print("########### MONKEY CHILD ############")
  concepts <-
    as.data.frame(
      read.csv2(
        file = filename,
        sep = ":",
        stringsAsFactors = FALSE,
        header = FALSE
      )[2],
      stringAsFactors = FALSE
    )
  print("## CONCEPTS ###")
  print(concepts)
  write.table(
    file = filename,
    x = concepts,
    append = F,
    quote = F,
    row.names = F
  )
  return(concepts)
}
cleanupResultsFile_compiled <- cmpfun(cleanupResultsFile)

cleanFN <- function(x) {
  x <- 
    str_replace_all(string = x,
                    pattern = " ",
                    replacement = "_")
  x
}

parsex <- function(x) {
	xlist <- unlist(str_split(x, pattern = "/"))
	out = c()
	out$confounderSource <- xlist[1]
	out$varType <- xlist[2]
	out$rsID <- xlist[3]
	out$lev <- xlist[4]
	out$hoi <- xlist[5]
	out$hoiFN <- cleanFN(paste("utgz/", xlist[5], ".txt.gz", sep = ""))
	out$drug <- xlist[6]
	out$drugFN <- cleanFN(paste("utgz/", xlist[6], ".txt.gz", sep = ""))
	out
}
 # parsed <- parsex("confounders/omop/1/liver_failure/zidovudine")

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)),
        substring(s, 2),
        sep = "",
        collapse = " ")
}

writeVecAppend <- function(dat, fn) {
  write.table(
    unique(dat),
    file = fn,
    quote = FALSE,
    row.names = FALSE,
    eol = "\\n",
    col.names = FALSE,
    append = TRUE
  )
}

newResultsFolderRollout <- function(cl, drug, ade, adr, group, SQUELCH) {
  clDir <- paste("output/", cl, sep = "")
  if (!file.exists(clDir)) { system(paste("mkdir ", clDir, sep = "")) }
  groupDir <- paste(clDir, "/", group, sep = "")
  if (!file.exists(groupDir)) { system(paste("mkdir ", groupDir, sep = "")) }
  adrDir <- paste(groupDir, "/", adr, sep = "")
  if (!file.exists(adrDir)) { system(paste("mkdir ", adrDir, sep = "")) }
  adeDir <- paste(adrDir, "/", ade, sep = "")
  if (!file.exists(adeDir)) { system(paste("mkdir ", adeDir, sep = "")) }
  drugDir <- paste(adeDir, "/", drug, sep = "")
  if (!file.exists(drugDir)) { system(paste("mkdir ", drugDir, sep = "")) }
  finalDir <- paste(drugDir, "/", SQUELCH, sep = "")
  if (!file.exists(finalDir)) { system(paste("mkdir ", finalDir, sep = "")) }
  finalDir
}

outputResultsFolder <- function(cl, drug, ade, adr, group, SQUELCH) {
  clDir <- paste("confounderLists/", cl, SQUELCH, "/", cl, sep = "")
  #if (!file.exists(clDir)) { system(paste("mkdir ", clDir, sep = "")) }
  groupDir <- paste(clDir, "/", group, sep = "")
  #if (!file.exists(groupDir)) { system(paste("mkdir ", groupDir, sep = "")) }
  adrDir <- paste(groupDir, "/", adr, sep = "")
  #if (!file.exists(adrDir)) { system(paste("mkdir ", adrDir, sep = "")) }
  adeDir <- paste(adrDir, "/", ade, sep = "")
  #if (!file.exists(adeDir)) { system(paste("mkdir ", adeDir, sep = "")) }
  drugDir <- paste(adeDir, "/", drug, sep = "")
  #if (!file.exists(drugDir)) { system(paste("mkdir ", drugDir, sep = "")) }
  finalDir <- paste(drugDir, "/", SQUELCH, sep = "")
  #if (!file.exists(finalDir)) { system(paste("mkdir ", finalDir, sep = "")) }
  finalDir
}



constructDAG <- function(dat.names) {
  # function to build DAG modelstring from confounder list
  #dat.names <- c("fever", "pain")
  buildDAG <- function(dat.names) {
    bn.names.string <- ""
    bn.names.string.fin <- ""
    getNameString <- function(name.string) {
      n.str <- ""
      n.str <- paste("[", name.string, "]", sep = "")
    }
    bn.names.string <- sapply(dat.names, getNameString)
    bn.names.string.fin <- ""
    for (i in bn.names.string)
    {
      bn.names.string.fin <- paste(bn.names.string.fin, i, sep = "")
    }
    print(bn.names.string.fin)
  }
  buildDAG(dat.names)
  
  
  buildDAGade <- function(dat.names) {
    bn.names.string <- ""
    bn.names.string.fin <- ""
    getNameString <- function(name.string) {
      n.str <- ""
      n.str <- paste(name.string, ":", sep = "")
    }
    bn.names.string <- sapply(dat.names, getNameString)
    bn.names.string.fin <- "ade|"
    for (i in bn.names.string)
    {
      bn.names.string.fin <- paste(bn.names.string.fin, i, sep = "")
    }
    print(bn.names.string.fin)
  }
  buildDAGdrug <- function(dat.names) {
    bn.names.string <- ""
    bn.names.string.fin <- ""
    getNameString <- function(name.string) {
      n.str <- ""
      n.str <- paste(name.string, ":", sep = "")
    }
    bn.names.string <- sapply(dat.names, getNameString)
    bn.names.string.fin <- "medication|"
    for (i in bn.names.string)
    {
      bn.names.string.fin <- paste(bn.names.string.fin, i, sep = "")
    }
    print(bn.names.string.fin)
  }
  
  predSTR <- buildDAG(dat.names)
  adeSTR <-
    paste("[",
          gsubfn(
            buildDAGade(c("medication", dat.names)),
            pattern = ":$",
            replacement = ""
          ),
          "]",
          sep = "")
  drugSTR <-
    paste("[",
          gsubfn(
            buildDAGdrug(c(dat.names)),
            pattern = ":$",
            replacement = ""
          ),
          "]",
          sep = "")
  dagSTR <- paste(predSTR, drugSTR, adeSTR, sep = "")
  print("###################DAG #########################")
  
  print(dagSTR)
  dag <- bnlearn::model2network(dagSTR)
  plot(dag)
  return(dag)
}

isEmpty <- function(x) {
  return(length(x)==0)
}



