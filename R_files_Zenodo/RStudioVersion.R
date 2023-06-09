#written by IC, 08.04.2020

#Output current version of RStudio to a text file for reporting purposes.

vers <- as.character(RStudio.Version()$version)
writeLines(c(vers, "\\n"), "RStudioVersion.txt")
