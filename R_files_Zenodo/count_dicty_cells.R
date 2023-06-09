
## Script for analyzing fluorescent Dictyostelium cells via flow cytometry
## For use with R statistical computing environment (r-project.org)
## 
## jeff smith <jeffsmith@wustl.edu>, Dept of Biology, Washington Univ in St Louis
## 2011-2012


# 
# READ IN PACKAGES & DATA
# 

# If you don't have the packages: 
# 	source("http://bioconductor.org/biocLite.R")
# 	biocLite("flowCore")
# 	biocLite("flowViz")
# 	biocLite("flowStats")

print("Loading bioconductor.org flow cytometry packages ...", quote = FALSE)
require(flowCore)
require(flowViz)
require(flowStats)

print("Reading data from .fcs files in directory 'fcs_files' ...", quote = FALSE)
myData <- read.flowSet(path = "fcs_files", phenoData = list(sampleDescription="#SAMPLE", dataCell = "$SMNO", dateCollected = "$DATE", volume = "$VOL")) # script might not work properly if there's only one .fcs file
# pData(phenoData(myData)) # summary of metadata

myWorkFlow <- workFlow(myData, name = "myWorkFlow") # The Bioconductor flow cytometry packages use "workFlow" objects to manage transformations, gates, and naming of intermediate results 


# 
# PROCESSING AND CALCULATIONS
# 

print("Transforming data (asinh) ...", quote = FALSE)
myTransformation <- transformList(colnames(myData)[1:12], tfun = asinh, transformationId = "asinh")
add(myWorkFlow, myTransformation)

print("Screening out non-cell debris ...", quote = FALSE)
# So lymphGate() won't decide to choose debris instead of Dicty cells
debrisScreen <- rectangleGate("FSC-H" = c(11,18), "SSC-H" = c(10,20), filterId = "debrisScreen") 
add(myWorkFlow, debrisScreen, parent = "asinh")

print("Isolating Dicty cells ...", quote = FALSE)
# Fit an oval gate (bivariate normal) to the data in the region where we expect Dicty spores: 
sporeGate <- lymphGate(Data(myWorkFlow[["debrisScreen+"]]), channels = c("FSC-H","SSC-H"), preselection = list("FSC-H" = c(13.5, 16.5), "SSC-H" = c(10, 15)), eval = FALSE, filterId = "DictyCells")
add(myWorkFlow, sporeGate$n2gate, parent="debrisScreen+")

rfpBoundaryFilter <- boundaryFilter(filterId = "rfpBoundaryFilter", x = "FL3-H") # So that crap in the Dicty region doesn't throw off our numbers
add(myWorkFlow, rfpBoundaryFilter, parent = "DictyCells+")

print("Aligning rfp distribution peaks ...", quote = FALSE)
alignRFPpeaks <- normalization(normFun = function(x, parameters, ...) warpSet(x, parameters,...), parameters = c("FL3-H"), normalizationId = "alignRFPpeaks") # Helps rangeGate() create an appropriate breakpoint between rfp+ and rfp- cells
add(myWorkFlow, alignRFPpeaks, parent = "rfpBoundaryFilter+") 

print("Measuring rfp +/- ...", quote = FALSE)
rfpGate <- rangeGate(Data(myWorkFlow[["alignRFPpeaks"]]), "FL3-H", plot = FALSE, filterId = "rfp", sd = 2.5, refLine = 4) # Automagically create dividing line between rfp+ and rfp- cells
add(myWorkFlow, rfpGate, parent = "alignRFPpeaks")

print("Exporting results as tab-delimitted text file called 'results.txt' ...", quote = FALSE)
# summary(myWorkFlow[["rfp+"]]) 
results <- pData(phenoData(myData))
results <- cbind(results, "rfp+cells" = summary(myWorkFlow[["rfp+"]])$true, "rfp-cells" = summary(myWorkFlow[["rfp+"]])$false, "fractionRfp+" = summary(myWorkFlow[["rfp+"]])$p)
write.table(results, file = "results.txt", sep = "\\t", quote = FALSE, row.names = FALSE)


# 
# PLOT DATA AND GATES TO MAKE SURE EVERYTHING IS ACTING AS IT SHOULD
# 

print("Plotting data ...", quote = FALSE)
devAskNewPage(ask = TRUE)

# print(xyplot(`SSC-H` ~ `FSC-H`, myWorkFlow[["asinh"]], xlim = c(11,18), ylim = c(0,20), main = "Flow cytometry"))

print(xyplot(`SSC-H` ~ `FSC-H`, myWorkFlow[["DictyCells+"]], par.settings = list(gate = list(col = "red", fill = "red", alpha = 0.3)), xlim = c(11,18), ylim = c(8,20), main = "Gate for Dicty cells"))

print(densityplot(~`FL3-H`, Data(myWorkFlow[["rfpBoundaryFilter+"]]), main = "Fluorescence data", xlim = c(3,13) ), split = c(1,1,2,1), more = TRUE)

print(densityplot(~`FL3-H`, Data(myWorkFlow[["alignRFPpeaks"]]), refline = rfpGate@min, main = "Split between rfp+ and rfp- cells", xlim = c(3,13) ), split = c(2,1,2,1), more = FALSE)
