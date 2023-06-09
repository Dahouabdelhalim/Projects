# 104 individuals, 323 markers...
#
##
###
####
#####
######
############# ORIGINAL
# For PC
setwd("D:/PhD/Crossing/Analyses/Recomb.Map/Population/FINAL_DNA_PHENOTYPES")
getwd()

library(onemap)

# 104 individuals
LGs_1 <- read.outcross(,"Recombination_Mapping.txt")
LGs_2pt <- rf.2pts(LGs_1, LOD=8, max.rf=.4)
LGs_2pt_makeseq <- make.seq (LGs_2pt, "all")
marker.type(LGs_2pt_makeseq)
Grouped_LGs <- group(LGs_2pt_makeseq)
print(Grouped_LGs, detailed=FALSE)
print(Grouped_LGs)

set.map.fun(type="kosambi")

# First view and evaluate the first LG (LG_A) 
###### LG_A \\/ #########################
LG_A.frame1 <- make.seq(LGs_2pt, c(232,224,284))
LG_A.compare <-compare(LG_A.frame1)
LG_A.compare
LG_A.frame1 <- make.seq(LG_A.compare,1,1)
LG_A.frame1
### BACKBONE A.1 ### /\\ ####################

### IN BACKBONE \\/ ###
LG_A.frame_add224 <- try.seq(LG_A.frame2, 224, draw.try=TRUE)
LG_A.frame_add224 <- try.seq(LG_A.frame2, 224)
LG_A.frame_add224
LG_A.frame2 <- make.seq(LG_A.frame_add224,1,1)
LG_A.frame2
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add232 <- try.seq(LG_A.frame2, 232, draw.try=TRUE)
LG_A.frame_add232 <- try.seq(LG_A.frame2, 232)
LG_A.frame_add232
LG_A.frame2 <- make.seq(LG_A.frame_add232,1,1)
LG_A.frame2
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add284 <- try.seq(LG_A.frame1, 284, draw.try=TRUE)
LG_A.frame_add284 <- try.seq(LG_A.frame1, 284)
LG_A.frame_add284
LG_A.frame1 <- make.seq(LG_A.frame_add284,8,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])
### IN BACKBONE /\\ ###

##### Test other markers by adding one at a time to see if they fit the following criteria:
# Markers that did not show recombination frequencies monotonically increasing with distance 
# from the diagonal of the recombination frequency heat map were relocated using the try.seq 
# and make.seq functions, or removed. Once all markers within the LG displayed a monotonic 
# recombination frequency pattern, we forced each other marker initially grouped with those 
# markers onto the LG, one at a time, to determine if they fit soundly at any position along the 
# linage group. If forcing a marker onto the LG resulted in map expansion or violation of 
# monotony, we relocated or removed it. 

## Frame 1 \\/ ##############################
LG_A.frame_add78 <- try.seq(LG_A.frame1, 78, draw.try=TRUE)
LG_A.frame_add78 <- try.seq(LG_A.frame1, 78)
LG_A.frame_add78
LG_A.frame1 <- make.seq(LG_A.frame_add78,4,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add144 <- try.seq(LG_A.frame1, 144, draw.try=TRUE)
LG_A.frame_add144 <- try.seq(LG_A.frame1, 144)
LG_A.frame_add144
LG_A.frame1 <- make.seq(LG_A.frame_add144,4,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add251 <- try.seq(LG_A.frame1, 251, draw.try=TRUE)
LG_A.frame_add251 <- try.seq(LG_A.frame1, 251)
LG_A.frame_add251
LG_A.frame1 <- make.seq(LG_A.frame_add251,1,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add243 <- try.seq(LG_A.frame1, 243, draw.try=TRUE)
LG_A.frame_add243 <- try.seq(LG_A.frame1, 243)
LG_A.frame_add243
LG_A.frame1 <- make.seq(LG_A.frame_add243,7,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add73 <- try.seq(LG_A.frame1, 73, draw.try=TRUE)
LG_A.frame_add73 <- try.seq(LG_A.frame1, 73)
LG_A.frame_add73
LG_A.frame1 <- make.seq(LG_A.frame_add73,7,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

# \\/ LESS THAN 80% OF INDIVIDUALS SCORED \\/ #
LG_A.frame_add134 <- try.seq(LG_A.frame1, 134, draw.try=TRUE)
LG_A.frame_add134 <- try.seq(LG_A.frame1, 134)
LG_A.frame_add134
LG_A.frame1 <- make.seq(LG_A.frame_add134,6,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add262 <- try.seq(LG_A.frame1, 262, draw.try=TRUE)
LG_A.frame_add262 <- try.seq(LG_A.frame1, 262)
LG_A.frame_add262
LG_A.frame1 <- make.seq(LG_A.frame_add262,4,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add239 <- try.seq(LG_A.frame1, 239, draw.try=TRUE)
LG_A.frame_add239 <- try.seq(LG_A.frame1, 239)
LG_A.frame_add239
LG_A.frame1 <- make.seq(LG_A.frame_add239,5,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add257 <- try.seq(LG_A.frame1, 257, draw.try=TRUE)
LG_A.frame_add257 <- try.seq(LG_A.frame1, 257)
LG_A.frame_add257
LG_A.frame1 <- make.seq(LG_A.frame_add257,5,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add110 <- try.seq(LG_A.frame1, 110, draw.try=TRUE)
LG_A.frame_add110 <- try.seq(LG_A.frame1, 110)
LG_A.frame_add110
LG_A.frame1 <- make.seq(LG_A.frame_add110,11,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add309 <- try.seq(LG_A.frame1, 309, draw.try=TRUE)
LG_A.frame_add309 <- try.seq(LG_A.frame1, 309)
LG_A.frame_add309
LG_A.frame1 <- make.seq(LG_A.frame_add309,3,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add98 <- try.seq(LG_A.frame1, 98, draw.try=TRUE)
LG_A.frame_add98 <- try.seq(LG_A.frame1, 98)
LG_A.frame_add98
LG_A.frame1 <- make.seq(LG_A.frame_add98,15,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add303 <- try.seq(LG_A.frame1, 303, draw.try=TRUE)
LG_A.frame_add303 <- try.seq(LG_A.frame1, 303)
LG_A.frame_add303
LG_A.frame1 <- make.seq(LG_A.frame_add303,15,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add219 <- try.seq(LG_A.frame1, 219, draw.try=TRUE)
LG_A.frame_add219 <- try.seq(LG_A.frame1, 219)
LG_A.frame_add219
LG_A.frame1 <- make.seq(LG_A.frame_add219,15,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

rf.graph.table(LG_A.frame1, scale=2)
ripple.seq(LG_A.frame1, ws=6, LOD=1, tol=0.1)
maps <- list(LG_A.frame1)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=1.5)

## LG_A Frame 1 /\\ ################# TOTAL OF 98.27 CM ####
#
#
#
#
#
#
#
#
#
#
#

## LG_B 2 Frame 2 \\/ ##############################
LG_B.frame2 <- make.seq(LGs_2pt, c(164, 167,139,1))
LG_B.compare <-compare(LG_B.frame2)
LG_B.compare
LG_B.frame2 <- make.seq(LG_B.compare,1,1)
LG_B.frame2

### IN BACKBONE \\/ ###
LG_B.frame_add167 <- try.seq(LG_B.frame2, 167, draw.try=TRUE)
LG_B.frame_add167 <- try.seq(LG_B.frame2, 167)
LG_B.frame_add167
LG_B.frame2 <- make.seq(LG_B.frame_add167,8,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add164 <- try.seq(LG_B.frame2, 164, draw.try=TRUE)
LG_B.frame_add164 <- try.seq(LG_B.frame2, 164)
LG_B.frame_add164
LG_B.frame2 <- make.seq(LG_B.frame_add164,7,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add1 <- try.seq(LG_B.frame2, 1, draw.try=TRUE)
LG_B.frame_add1 <- try.seq(LG_B.frame2, 1)
LG_B.frame_add1
LG_B.frame2 <- make.seq(LG_B.frame_add1,5,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add139 <- try.seq(LG_B.frame2, 139, draw.try=TRUE)
LG_B.frame_add139 <- try.seq(LG_B.frame2, 139)
LG_B.frame_add139
LG_B.frame2 <- make.seq(LG_B.frame_add139,3,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])
### IN BACKBONE /\\ ###

## LG_B 2 \\/ ##############################
LG_B.frame_add206 <- try.seq(LG_B.frame2, 206, draw.try=TRUE)
LG_B.frame_add206 <- try.seq(LG_B.frame2, 206)
LG_B.frame_add206
LG_B.frame2 <- make.seq(LG_B.frame_add206,5,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add182 <- try.seq(LG_B.frame2, 182, draw.try=TRUE)
LG_B.frame_add182 <- try.seq(LG_B.frame2, 182)
LG_B.frame_add182
LG_B.frame2 <- make.seq(LG_B.frame_add182,1,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add178 <- try.seq(LG_B.frame2, 178, draw.try=TRUE)
LG_B.frame_add178 <- try.seq(LG_B.frame2, 178)
LG_B.frame_add178
LG_B.frame2 <- make.seq(LG_B.frame_add178,6,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add190 <- try.seq(LG_B.frame2, 190, draw.try=TRUE)
LG_B.frame_add190 <- try.seq(LG_B.frame2, 190)
LG_B.frame_add190
LG_B.frame2 <- make.seq(LG_B.frame_add190,2,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add203 <- try.seq(LG_B.frame2, 203, draw.try=TRUE)
LG_B.frame_add203 <- try.seq(LG_B.frame2, 203)
LG_B.frame_add203
LG_B.frame2 <- make.seq(LG_B.frame_add203,6,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add189 <- try.seq(LG_B.frame2, 189, draw.try=TRUE)
LG_B.frame_add189 <- try.seq(LG_B.frame2, 189)
LG_B.frame_add189
LG_B.frame2 <- make.seq(LG_B.frame_add189,9,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add102 <- try.seq(LG_B.frame2, 102, draw.try=TRUE)
LG_B.frame_add102 <- try.seq(LG_B.frame2, 102)
LG_B.frame_add102
LG_B.frame2 <- make.seq(LG_B.frame_add102,3,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add181 <- try.seq(LG_B.frame2, 181, draw.try=TRUE)
LG_B.frame_add181 <- try.seq(LG_B.frame2, 181)
LG_B.frame_add181
LG_B.frame2 <- make.seq(LG_B.frame_add181,8,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

rf.graph.table(LG_B.frame2, scale=2)
ripple.seq(LG_B.frame2, ws=6, LOD=1, tol=0.1)
maps <- list(LG_B.frame2)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=1.5)

## LG_B Frame 2 /\\ #### 83.67 CM ##########################
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
## LG_C Frame 3 \\/ ##############################
LG_C.frame3 <- make.seq(LGs_2pt, c(19,271,268))
LG_C.compare <-compare(LG_C.frame3)
LG_C.compare
LG_C.frame3 <- make.seq(LG_C.compare,2,1)
LG_C.frame3

### IN BACKBONE \\/ ####
LG_C.frame_add19 <- try.seq(LG_C.frame3, 19, draw.try=TRUE)
LG_C.frame_add19 <- try.seq(LG_C.frame3, 19)
LG_C.frame_add19
LG_C.frame3 <- make.seq(LG_C.frame_add19,4,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add271 <- try.seq(LG_C.frame3, 271, draw.try=TRUE)
LG_C.frame_add271 <- try.seq(LG_C.frame3, 271)
LG_C.frame_add271
LG_C.frame3 <- make.seq(LG_C.frame_add271,2,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add268 <- try.seq(LG_C.frame3, 268, draw.try=TRUE)
LG_C.frame_add268 <- try.seq(LG_C.frame3, 268)
LG_C.frame_add268
LG_C.frame3 <- make.seq(LG_C.frame_add268,5,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])
### IN BACKBONE /\\ ####
#
#
####### LG_C \\/ ##############################
LG_C.frame_add302 <- try.seq(LG_C.frame3, 302, draw.try=TRUE)
LG_C.frame_add302 <- try.seq(LG_C.frame3, 302)
LG_C.frame_add302
LG_C.frame3 <- make.seq(LG_C.frame_add302,4,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add149 <- try.seq(LG_C.frame3, 149, draw.try=TRUE)
LG_C.frame_add149 <- try.seq(LG_C.frame3, 149)
LG_C.frame_add149
LG_C.frame3 <- make.seq(LG_C.frame_add149,5,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add301 <- try.seq(LG_C.frame3, 301, draw.try=TRUE)
LG_C.frame_add301 <- try.seq(LG_C.frame3, 301)
LG_C.frame_add301
LG_C.frame3 <- make.seq(LG_C.frame_add301,6,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add27 <- try.seq(LG_C.frame3, 27, draw.try=TRUE)
LG_C.frame_add27 <- try.seq(LG_C.frame3, 27)
LG_C.frame_add27
LG_C.frame3 <- make.seq(LG_C.frame_add27,7,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add37 <- try.seq(LG_C.frame3, 37, draw.try=TRUE)
LG_C.frame_add37 <- try.seq(LG_C.frame3, 37)
LG_C.frame_add37
LG_C.frame3 <- make.seq(LG_C.frame_add37,1,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add273 <- try.seq(LG_C.frame3, 273, draw.try=TRUE)
LG_C.frame_add273 <- try.seq(LG_C.frame3, 273)
LG_C.frame_add273
LG_C.frame3 <- make.seq(LG_C.frame_add273,1,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add213 <- try.seq(LG_C.frame3, 213, draw.try=TRUE)
LG_C.frame_add213 <- try.seq(LG_C.frame3, 213)
LG_C.frame_add213
LG_C.frame3 <- make.seq(LG_C.frame_add213,1,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add306 <- try.seq(LG_C.frame3, 306, draw.try=TRUE)
LG_C.frame_add306 <- try.seq(LG_C.frame3, 306)
LG_C.frame_add306
LG_C.frame3 <- make.seq(LG_C.frame_add306,1,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add212 <- try.seq(LG_C.frame3, 212, draw.try=TRUE)
LG_C.frame_add212 <- try.seq(LG_C.frame3, 212)
LG_C.frame_add212
LG_C.frame3 <- make.seq(LG_C.frame_add212,1,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add107 <- try.seq(LG_C.frame3, 107, draw.try=TRUE)
LG_C.frame_add107 <- try.seq(LG_C.frame3, 107)
LG_C.frame_add107
LG_C.frame3 <- make.seq(LG_C.frame_add107,13,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add183 <- try.seq(LG_C.frame3, 183, draw.try=TRUE)
LG_C.frame_add183 <- try.seq(LG_C.frame3, 183)
LG_C.frame_add183
LG_C.frame3 <- make.seq(LG_C.frame_add183,14,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add123 <- try.seq(LG_C.frame3, 123, draw.try=TRUE)
LG_C.frame_add123 <- try.seq(LG_C.frame3, 123)
LG_C.frame_add123
LG_C.frame3 <- make.seq(LG_C.frame_add123,6,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

rf.graph.table(LG_C.frame3, scale=2)
ripple.seq(LG_C.frame3, ws=6, LOD=1, tol=0.1)
maps <- list(LG_C.frame3)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=.75)
####################################
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
## LG_D Frame 4 \\/ ##############################
LG_D.tryHARD2 <- make.seq(LGs_2pt, c(23,312,226,308,280,188,7,249,205,187,147,221,286,281,265,202,197,21,267,209,26,16,191,171)) 
LG_D_tryTD2 <- order.seq(LG_D.tryHARD2, subset.n.try=30, THRES=2, touchdown = TRUE)
LG_D_TD4 <- make.seq(LG_D_tryTD2)
print(LG_D_TD4)

LG_D.frame4 <- LG_D_TD4

rf.graph.table(LG_D_TD4, scale=2)
ripple.seq(LG_D_TD4, ws=6, LOD=1, tol=0.1)
maps <- list(LG_D_TD4)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=1)


### IN BACKBONE \\/ #### 316, 23, 7, 147, 21, 26, 16, 188, 191, 202, 209, 226, 265, 286
LG_D.frame_add316 <- try.seq(LG_D.frame4, 316, draw.try=TRUE)
LG_D.frame_add316 <- try.seq(LG_D.frame4, 316)
LG_D.frame_add316
LG_D.frame4 <- make.seq(LG_D.frame_add316,1,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add23 <- try.seq(LG_D.frame4, 23, draw.try=TRUE)
LG_D.frame_add23 <- try.seq(LG_D.frame4, 23)
LG_D.frame_add23
LG_D.frame4 <- make.seq(LG_D.frame_add23,1,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add7 <- try.seq(LG_D.frame4, 7, draw.try=TRUE)
LG_D.frame_add7 <- try.seq(LG_D.frame4, 7)
LG_D.frame_add7
LG_D.frame4 <- make.seq(LG_D.frame_add7,1,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add147 <- try.seq(LG_D.frame4, 147, draw.try=TRUE)
LG_D.frame_add147 <- try.seq(LG_D.frame4, 147)
LG_D.frame_add147
LG_D.frame4 <- make.seq(LG_D.frame_add147,10,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add21 <- try.seq(LG_D.frame4, 21, draw.try=TRUE)
LG_D.frame_add21 <- try.seq(LG_D.frame4, 21)
LG_D.frame_add21
LG_D.frame4 <- make.seq(LG_D.frame_add21,1,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add26 <- try.seq(LG_D.frame4, 26, draw.try=TRUE)
LG_D.frame_add26 <- try.seq(LG_D.frame4, 26)
LG_D.frame_add26
LG_D.frame4 <- make.seq(LG_D.frame_add26,9,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add16 <- try.seq(LG_D.frame4, 16, draw.try=TRUE)
LG_D.frame_add16 <- try.seq(LG_D.frame4, 16)
LG_D.frame_add16
LG_D.frame4 <- make.seq(LG_D.frame_add16,1,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add188 <- try.seq(LG_D.frame4, 188, draw.try=TRUE)
LG_D.frame_add188 <- try.seq(LG_D.frame4, 188)
LG_D.frame_add188
LG_D.frame4 <- make.seq(LG_D.frame_add188,5,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add191 <- try.seq(LG_D.frame4, 191, draw.try=TRUE)
LG_D.frame_add191 <- try.seq(LG_D.frame4, 191)
LG_D.frame_add191
LG_D.frame4 <- make.seq(LG_D.frame_add191,1,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add202 <- try.seq(LG_D.frame4, 202, draw.try=TRUE)
LG_D.frame_add202 <- try.seq(LG_D.frame4, 202)
LG_D.frame_add202
LG_D.frame4 <- make.seq(LG_D.frame_add202,15,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add209 <- try.seq(LG_D.frame4, 209, draw.try=TRUE)
LG_D.frame_add209 <- try.seq(LG_D.frame4, 209)
LG_D.frame_add209
LG_D.frame4 <- make.seq(LG_D.frame_add209,20,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add265 <- try.seq(LG_D.frame4, 265, draw.try=TRUE)
LG_D.frame_add265 <- try.seq(LG_D.frame4, 265)
LG_D.frame_add265
LG_D.frame4 <- make.seq(LG_D.frame_add265,14,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add286 <- try.seq(LG_D.frame4, 286, draw.try=TRUE)
LG_D.frame_add286 <- try.seq(LG_D.frame4, 286)
LG_D.frame_add286
LG_D.frame4 <- make.seq(LG_D.frame_add286,12,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add226 <- try.seq(LG_D.frame4, 226, draw.try=TRUE)
LG_D.frame_add226 <- try.seq(LG_D.frame4, 226)
LG_D.frame_add226
LG_D.frame4 <- make.seq(LG_D.frame_add226,4,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])


### IN BACKBONE /\\ #### 316, 23, 7, 147, 21, 26, 16, 188, 191, 202, 209, 226, 265, 286

LG_D.frame4 <- LG_D_TD4
LG_D.frame4

LG_D.frame_add221 <- try.seq(LG_D.frame4, 221, draw.try=TRUE)
LG_D.frame_add221 <- try.seq(LG_D.frame4, 221)
LG_D.frame_add221
LG_D.frame4 <- make.seq(LG_D.frame_add221,10,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add308 <- try.seq(LG_D.frame4, 308, draw.try=TRUE)
LG_D.frame_add308 <- try.seq(LG_D.frame4, 308)
LG_D.frame_add308
LG_D.frame4 <- make.seq(LG_D.frame_add308,7,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add281 <- try.seq(LG_D.frame4, 281, draw.try=TRUE)
LG_D.frame_add281 <- try.seq(LG_D.frame4, 281)
LG_D.frame_add281
LG_D.frame4 <- make.seq(LG_D.frame_add281,12,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add280 <- try.seq(LG_D.frame4, 280, draw.try=TRUE)
LG_D.frame_add280 <- try.seq(LG_D.frame4, 280)
LG_D.frame_add280
LG_D.frame4 <- make.seq(LG_D.frame_add280,7,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add267 <- try.seq(LG_D.frame4, 267, draw.try=TRUE)
LG_D.frame_add267 <- try.seq(LG_D.frame4, 267)
LG_D.frame_add267
LG_D.frame4 <- make.seq(LG_D.frame_add267,15,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add187 <- try.seq(LG_D.frame4, 187, draw.try=TRUE)
LG_D.frame_add187 <- try.seq(LG_D.frame4, 187)
LG_D.frame_add187
LG_D.frame4 <- make.seq(LG_D.frame_add187,14,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

rf.graph.table(LG_D.frame4, scale=2)
ripple.seq(LG_D.frame4, ws=6, LOD=1, tol=0.1)
maps <- list(LG_D.frame4)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=1)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
LG_E.tryHARD.2 <- make.seq(LGs_2pt, c(93,142,88,103))
LG_E_tryTD.2 <- order.seq(LG_E.tryHARD.2, subset.n.try=30, THRES=2, touchdown = TRUE)
LG_E_TD.2 <- make.seq(LG_E_tryTD.2)
print(LG_E_TD.2)

LG_E.frame5 <- LG_E_TD.2

rf.graph.table(LG_E.frame5, scale=2)
ripple.seq(LG_E.frame5, ws=6, LOD=1, tol=0.1)
maps <- list(LG_E.frame5)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=1)


### IN BACKBONE \\/ ####


LG_E.frame_add142 <- try.seq(LG_E.frame5, 142, draw.try=TRUE)
LG_E.frame_add142 <- try.seq(LG_E.frame5, 142)
LG_E.frame_add142
LG_E.frame5 <- make.seq(LG_E.frame_add142,29,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

LG_E.frame_add103 <- try.seq(LG_E.frame5, 103, draw.try=TRUE)
LG_E.frame_add103 <- try.seq(LG_E.frame5, 103)
LG_E.frame_add103
LG_E.frame5 <- make.seq(LG_E.frame_add103,31,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

LG_E.frame_add88 <- try.seq(LG_E.frame5, 88, draw.try=TRUE)
LG_E.frame_add88 <- try.seq(LG_E.frame5, 88)
LG_E.frame_add88
LG_E.frame5 <- make.seq(LG_E.frame_add88,31,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

LG_E.frame_add93 <- try.seq(LG_E.frame5, 93, draw.try=TRUE)
LG_E.frame_add93 <- try.seq(LG_E.frame5, 93)
LG_E.frame_add93
LG_E.frame5 <- make.seq(LG_E.frame_add93,26,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

### IN BACKBONE /\\ #### ########################
# 
#
#

LG_E.frame_add300 <- try.seq(LG_E.frame5, 300, draw.try=TRUE)
LG_E.frame_add300 <- try.seq(LG_E.frame5, 300)
LG_E.frame_add300
LG_E.frame5 <- make.seq(LG_E.frame_add300,2,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

LG_E.frame_add261 <- try.seq(LG_E.frame5, 261, draw.try=TRUE)
LG_E.frame_add261 <- try.seq(LG_E.frame5, 261)
LG_E.frame_add261
LG_E.frame5 <- make.seq(LG_E.frame_add261,6,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

LG_E.frame_add263 <- try.seq(LG_E.frame5, 263, draw.try=TRUE)
LG_E.frame_add263 <- try.seq(LG_E.frame5, 263)
LG_E.frame_add263
LG_E.frame5 <- make.seq(LG_E.frame_add263,4,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

LG_E.frame_add80 <- try.seq(LG_E.frame5, 80, draw.try=TRUE)
LG_E.frame_add80 <- try.seq(LG_E.frame5, 80)
LG_E.frame_add80
LG_E.frame5 <- make.seq(LG_E.frame_add80,7,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

LG_E.frame_add151 <- try.seq(LG_E.frame5, 151, draw.try=TRUE)
LG_E.frame_add151 <- try.seq(LG_E.frame5, 151)
LG_E.frame_add151
LG_E.frame5 <- make.seq(LG_E.frame_add151,6,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

LG_E.frame_add106 <- try.seq(LG_E.frame5, 106, draw.try=TRUE)
LG_E.frame_add106 <- try.seq(LG_E.frame5, 106)
LG_E.frame_add106
LG_E.frame5 <- make.seq(LG_E.frame_add106,8,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

LG_E.frame_add121 <- try.seq(LG_E.frame5, 121, draw.try=TRUE)
LG_E.frame_add121 <- try.seq(LG_E.frame5, 121)
LG_E.frame_add121
LG_E.frame5 <- make.seq(LG_E.frame_add121,7,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

rf.graph.table(LG_E.frame5, scale=2)
ripple.seq(LG_E.frame5, ws=6, LOD=1, tol=0.1)
maps <- list(LG_E.frame5)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=1)


################ END LG_E /\\ ######################
#
#
#
#
#
#
#
#
#
#
#
#
#
#
LG_F.tryHARD.2 <- make.seq(LGs_2pt, c(9,48,56,62,91,99,275))
LG_F_tryTD.2 <- order.seq(LG_F.tryHARD.2, subset.n.try=30, THRES=2, touchdown = TRUE)
LG_F_TD.2 <- make.seq(LG_F_tryTD.2)
print(LG_F_TD.2)

LG_F.frame6 <- LG_F_TD.2

rf.graph.table(LG_F.frame6, scale=2)
ripple.seq(LG_F.frame6, ws=6, LOD=1, tol=0.1)
maps <- list(LG_F.frame6)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=1)

## ## BACKBONE \\/ ###
LG_F.frame_add9 <- try.seq(LG_F.frame6, 9, draw.try=TRUE)
LG_F.frame_add9 <- try.seq(LG_F.frame6, 9)
LG_F.frame_add9
LG_F.frame6 <- make.seq(LG_F.frame_add9,1,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

LG_F.frame_add91 <- try.seq(LG_F.frame6, 91, draw.try=TRUE)
LG_F.frame_add91 <- try.seq(LG_F.frame6, 91)
LG_F.frame_add91
LG_F.frame6 <- make.seq(LG_F.frame_add91,1,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

LG_F.frame_add56 <- try.seq(LG_F.frame6, 56, draw.try=TRUE)
LG_F.frame_add56 <- try.seq(LG_F.frame6, 56)
LG_F.frame_add56
LG_F.frame6 <- make.seq(LG_F.frame_add56,1,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

LG_F.frame_add48 <- try.seq(LG_F.frame6, 48, draw.try=TRUE)
LG_F.frame_add48 <- try.seq(LG_F.frame6, 48)
LG_F.frame_add48
LG_F.frame6 <- make.seq(LG_F.frame_add48,1,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

LG_F.frame_add62 <- try.seq(LG_F.frame6, 62, draw.try=TRUE)
LG_F.frame_add62 <- try.seq(LG_F.frame6, 62)
LG_F.frame_add62
LG_F.frame6 <- make.seq(LG_F.frame_add62,1,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])


## ## BACKBONE /\\ ###



LG_F.frame_add246 <- try.seq(LG_F.frame6, 246, draw.try=TRUE)
LG_F.frame_add246 <- try.seq(LG_F.frame6, 246)
LG_F.frame_add246
LG_F.frame6 <- make.seq(LG_F.frame_add246,4,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

LG_F.frame_add253 <- try.seq(LG_F.frame6, 253, draw.try=TRUE)
LG_F.frame_add253 <- try.seq(LG_F.frame6, 253)
LG_F.frame_add253
LG_F.frame6 <- make.seq(LG_F.frame_add253,3,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

LG_F.frame_add274 <- try.seq(LG_F.frame6, 274, draw.try=TRUE)
LG_F.frame_add274 <- try.seq(LG_F.frame6, 274)
LG_F.frame_add274
LG_F.frame6 <- make.seq(LG_F.frame_add274,5,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

LG_F.frame_add290 <- try.seq(LG_F.frame6, 290, draw.try=TRUE)
LG_F.frame_add290 <- try.seq(LG_F.frame6, 290)
LG_F.frame_add290
LG_F.frame6 <- make.seq(LG_F.frame_add290,8,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

LG_F.frame_add248 <- try.seq(LG_F.frame6, 248, draw.try=TRUE)
LG_F.frame_add248 <- try.seq(LG_F.frame6, 248)
LG_F.frame_add248
LG_F.frame6 <- make.seq(LG_F.frame_add248,6,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

LG_F.frame_add242 <- try.seq(LG_F.frame6, 242, draw.try=TRUE)
LG_F.frame_add242 <- try.seq(LG_F.frame6, 242)
LG_F.frame_add242
LG_F.frame6 <- make.seq(LG_F.frame_add242,6,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

rf.graph.table(LG_F.frame6, scale=2)
ripple.seq(LG_F.frame6, ws=6, LOD=1, tol=0.1)
maps <- list(LG_F.frame6)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=1)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
LG_G.Backbone1 <- make.seq(LGs_2pt, c(20,117,18,15,59,193,140,217,304))
LG_G_tryTD.2 <- order.seq(LG_G.Backbone1, subset.n.try=30, touchdown = TRUE)
LG_G_TD.2 <- make.seq(LG_G_tryTD.2)
print(LG_G_TD.2)
LG_G.frame7 <- LG_G_TD.2

rf.graph.table(LG_G.frame7, scale=2)
ripple.seq(LG_G.frame7, ws=6, LOD=1, tol=0.1)
maps <- list(LG_G.frame7)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=.5)

## ## BACKBONE \\/ ### TRY 84, 76, 172, 163, 
LG_G.frame_add15 <- try.seq(LG_G.frame7, 15, draw.try=TRUE)
LG_G.frame_add15 <- try.seq(LG_G.frame7, 15)
LG_G.frame_add15
LG_G.frame7 <- make.seq(LG_G.frame_add15,1,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

LG_G.frame_add18 <- try.seq(LG_G.frame7, 18, draw.try=TRUE)
LG_G.frame_add18 <- try.seq(LG_G.frame7, 18)
LG_G.frame_add18
LG_G.frame7 <- make.seq(LG_G.frame_add18,1,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

LG_G.frame_add20 <- try.seq(LG_G.frame7, 20, draw.try=TRUE)
LG_G.frame_add20 <- try.seq(LG_G.frame7, 20)
LG_G.frame_add20
LG_G.frame7 <- make.seq(LG_G.frame_add20,1,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

LG_G.frame_add59 <- try.seq(LG_G.frame7, 59, draw.try=TRUE)
LG_G.frame_add59 <- try.seq(LG_G.frame7, 59)
LG_G.frame_add59
LG_G.frame7 <- make.seq(LG_G.frame_add59,1,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

LG_G.frame_add117 <- try.seq(LG_G.frame7, 117, draw.try=TRUE)
LG_G.frame_add117 <- try.seq(LG_G.frame7, 117)
LG_G.frame_add117
LG_G.frame7 <- make.seq(LG_G.frame_add117,1,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

LG_G.frame_add140 <- try.seq(LG_G.frame7, 140, draw.try=TRUE)
LG_G.frame_add140 <- try.seq(LG_G.frame7, 140)
LG_G.frame_add140
LG_G.frame7 <- make.seq(LG_G.frame_add140,1,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

LG_G.frame_add193 <- try.seq(LG_G.frame7, 193, draw.try=TRUE)
LG_G.frame_add193 <- try.seq(LG_G.frame7, 193)
LG_G.frame_add193
LG_G.frame7 <- make.seq(LG_G.frame_add193,1,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

LG_G.frame_add217 <- try.seq(LG_G.frame7, 217, draw.try=TRUE)
LG_G.frame_add217 <- try.seq(LG_G.frame7, 217)
LG_G.frame_add217
LG_G.frame7 <- make.seq(LG_G.frame_add217,1,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

LG_G.frame_add304 <- try.seq(LG_G.frame7, 304, draw.try=TRUE)
LG_G.frame_add304 <- try.seq(LG_G.frame7, 304)
LG_G.frame_add304
LG_G.frame7 <- make.seq(LG_G.frame_add304,1,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])
# ## ## BACKBONE /\\ ### ########## ##########

LG_G.frame_add76 <- try.seq(LG_G.frame7, 76, draw.try=TRUE)
LG_G.frame_add76 <- try.seq(LG_G.frame7, 76)
LG_G.frame_add76
LG_G.frame7 <- make.seq(LG_G.frame_add76,1,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

LG_G.frame_add84 <- try.seq(LG_G.frame7, 84, draw.try=TRUE)
LG_G.frame_add84 <- try.seq(LG_G.frame7, 84)
LG_G.frame_add84
LG_G.frame7 <- make.seq(LG_G.frame_add84,1,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

LG_G.frame_add172 <- try.seq(LG_G.frame7, 172, draw.try=TRUE)
LG_G.frame_add172 <- try.seq(LG_G.frame7, 172)
LG_G.frame_add172
LG_G.frame7 <- make.seq(LG_G.frame_add172,2,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

rf.graph.table(LG_G.frame7, scale=2)
ripple.seq(LG_G.frame7, ws=6, LOD=1, tol=0.1)
maps <- list(LG_G.frame7)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=.5)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# LG_H \\/ ###################
LG_H.tryHARD <- make.seq(LGs_2pt, c(66,95,116,133,154,222,228,247,252,259,266,278,287,289)) 
LG_H_tryTD <- order.seq(LG_H.tryHARD,  subset.n.try=30, touchdown = TRUE)
LG_H_TD <- make.seq(LG_H_tryTD)
print(LG_H_TD)

### BACKBONE \\/ ###########
LG_H.frame_add66 <- try.seq(LG_H.frame8, 66, draw.try=TRUE)
LG_H.frame_add66 <- try.seq(LG_H.frame8, 66)
LG_H.frame_add66
LG_H.frame8 <- make.seq(LG_H.frame_add66,1,1)
LG_H.frame8
dev.off(dev.list()["RStudioGD"])

LG_H.frame_add95 <- try.seq(LG_H.frame8, 95, draw.try=TRUE)
LG_H.frame_add95 <- try.seq(LG_H.frame8, 95)
LG_H.frame_add95
LG_H.frame8 <- make.seq(LG_H.frame_add95,1,1)
LG_H.frame8
dev.off(dev.list()["RStudioGD"])

LG_H.frame_add133 <- try.seq(LG_H.frame8, 133, draw.try=TRUE)
LG_H.frame_add133 <- try.seq(LG_H.frame8, 133)
LG_H.frame_add133
LG_H.frame8 <- make.seq(LG_H.frame_add133,1,1)
LG_H.frame8
dev.off(dev.list()["RStudioGD"])

LG_H.frame_add154 <- try.seq(LG_H.frame8, 154, draw.try=TRUE)
LG_H.frame_add154 <- try.seq(LG_H.frame8, 154)
LG_H.frame_add154
LG_H.frame8 <- make.seq(LG_H.frame_add154,1,1)
LG_H.frame8
dev.off(dev.list()["RStudioGD"])

LG_H.frame_add116 <- try.seq(LG_H.frame8, 116, draw.try=TRUE)
LG_H.frame_add116 <- try.seq(LG_H.frame8, 116)
LG_H.frame_add116
LG_H.frame8 <- make.seq(LG_H.frame_add116,1,1)
LG_H.frame8
dev.off(dev.list()["RStudioGD"])

LG_H.frame_add228 <- try.seq(LG_H.frame8, 228, draw.try=TRUE)
LG_H.frame_add228 <- try.seq(LG_H.frame8, 228)
LG_H.frame_add228
LG_H.frame8 <- make.seq(LG_H.frame_add228,1,1)
LG_H.frame8
dev.off(dev.list()["RStudioGD"])
#### BACKBONE /\\ ###########

LG_H.frame8 <- LG_H_TD

rf.graph.table(LG_H.frame8, scale=2)
ripple.seq(LG_H.frame8, ws=4, LOD=1, tol=0.1)
maps <- list(LG_H.frame8)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=.5)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
LG_I.Backbone1 <- make.seq(LGs_2pt, c(125,119,236,105))
LG_I_tryTD.2 <- order.seq(LG_I.Backbone1, THRES=2, subset.n.try=30, touchdown = TRUE)
LG_I_TD.2 <- make.seq(LG_I_tryTD.2)
print(LG_I_TD.2)
LG_I.frame8 <- LG_I_TD.2

LG_I.Backbone1 <- make.seq(LGs_2pt, c(125,119,236,105))

rf.graph.table(LG_I.frame8, scale=2)
ripple.seq(LG_I.frame8, ws=6, LOD=1, tol=0.1)
maps <- list(LG_I.frame8)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=.5)

####### BACKBONE \\/ ########### ########### add 244,310,44
LG_I.frame_add125 <- try.seq(LG_I.frame8, 125, draw.try=TRUE)
LG_I.frame_add125 <- try.seq(LG_I.frame8, 125)
LG_I.frame_add125
LG_I.frame8 <- make.seq(LG_I.frame_add125,1,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])

LG_I.frame_add119 <- try.seq(LG_I.frame8, 119, draw.try=TRUE)
LG_I.frame_add119 <- try.seq(LG_I.frame8, 119)
LG_I.frame_add119
LG_I.frame8 <- make.seq(LG_I.frame_add119,1,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])

LG_I.frame_add236 <- try.seq(LG_I.frame8, 236, draw.try=TRUE)
LG_I.frame_add236 <- try.seq(LG_I.frame8, 236)
LG_I.frame_add236
LG_I.frame8 <- make.seq(LG_I.frame_add236,1,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])

LG_I.frame_add105 <- try.seq(LG_I.frame8, 105, draw.try=TRUE)
LG_I.frame_add105 <- try.seq(LG_I.frame8, 105)
LG_I.frame_add105
LG_I.frame8 <- make.seq(LG_I.frame_add105,1,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])

####### BACKBONE /\\ ########### ########### add 244,310,44
LG_I.frame_add244 <- try.seq(LG_I.frame8, 244, draw.try=TRUE)
LG_I.frame_add244 <- try.seq(LG_I.frame8, 244)
LG_I.frame_add244
LG_I.frame8 <- make.seq(LG_I.frame_add244,1,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])

LG_I.frame_add44 <- try.seq(LG_I.frame8, 44, draw.try=TRUE)
LG_I.frame_add44 <- try.seq(LG_I.frame8, 44)
LG_I.frame_add44
LG_I.frame8 <- make.seq(LG_I.frame_add44,1,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])

LG_I.frame_add310 <- try.seq(LG_I.frame8, 310, draw.try=TRUE)
LG_I.frame_add310 <- try.seq(LG_I.frame8, 310)
LG_I.frame_add310
LG_I.frame8 <- make.seq(LG_I.frame_add310,2,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])

rf.graph.table(LG_I.frame8, scale=2)
ripple.seq(LG_I.frame8, ws=6, LOD=1, tol=0.1)
maps <- list(LG_I.frame8)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=.5)
############## LG_I /\\ ######################
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
############# LG_J ## # ###########
LG_J.Backbone1 <- make.seq(LGs_2pt, c(13,29,34,35,42,43,52,111,113,114,148,153,169,174,180,185,196,200,214,220,223,235,238,250,255,283,295,297,298,299,315))


LG_J.Backbone1 <- make.seq(LGs_2pt, c(13,196,223,200,29))

LG_J_tryTD.2 <- order.seq(LG_J.Backbone1, THRES=2, subset.n.try=30, touchdown = TRUE)
LG_J_TD.2 <- make.seq(LG_J_tryTD.2)
print(LG_J_TD.2)
LG_J.frame8 <- LG_J_TD.2


rf.graph.table(LG_J.frame8, scale=2)
ripple.seq(LG_J.frame8, ws=6, LOD=1, tol=0.1)
maps <- list(LG_J.frame8)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=.5)
######## BACKBONE \\/ ##########
LG_J.frame_add13 <- try.seq(LG_J.frame8, 13, draw.try=TRUE)
LG_J.frame_add13 <- try.seq(LG_J.frame8, 13)
LG_J.frame_add13
LG_J.frame8 <- make.seq(LG_J.frame_add13,1,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

LG_J.frame_add196 <- try.seq(LG_J.frame8, 196, draw.try=TRUE)
LG_J.frame_add196 <- try.seq(LG_J.frame8, 196)
LG_J.frame_add196
LG_J.frame8 <- make.seq(LG_J.frame_add196,1,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

LG_J.frame_add223 <- try.seq(LG_J.frame8, 223, draw.try=TRUE)
LG_J.frame_add223 <- try.seq(LG_J.frame8, 223)
LG_J.frame_add223
LG_J.frame8 <- make.seq(LG_J.frame_add223,1,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

LG_J.frame_add200 <- try.seq(LG_J.frame8, 200, draw.try=TRUE)
LG_J.frame_add200 <- try.seq(LG_J.frame8, 200)
LG_J.frame_add200
LG_J.frame8 <- make.seq(LG_J.frame_add200,1,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

LG_J.frame_add29 <- try.seq(LG_J.frame8, 29, draw.try=TRUE)
LG_J.frame_add29 <- try.seq(LG_J.frame8, 29)
LG_J.frame_add29
LG_J.frame8 <- make.seq(LG_J.frame_add29,1,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

######## BACKBONE /\\ ##########
LG_J.frame8 <- LG_J_TD.2

LG_J.frame_add174 <- try.seq(LG_J.frame8, 174, draw.try=TRUE)
LG_J.frame_add174 <- try.seq(LG_J.frame8, 174)
LG_J.frame_add174
LG_J.frame8 <- make.seq(LG_J.frame_add174,5,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

LG_J.frame_add185 <- try.seq(LG_J.frame8, 185, draw.try=TRUE)
LG_J.frame_add185 <- try.seq(LG_J.frame8, 185)
LG_J.frame_add185
LG_J.frame8 <- make.seq(LG_J.frame_add185,6,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

LG_J.frame_add220 <- try.seq(LG_J.frame8, 220, draw.try=TRUE)
LG_J.frame_add220 <- try.seq(LG_J.frame8, 220)
LG_J.frame_add220
LG_J.frame8 <- make.seq(LG_J.frame_add220,7,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

rf.graph.table(LG_J.frame8, scale=2)
ripple.seq(LG_J.frame8, ws=6, LOD=1, tol=0.1)
maps <- list(LG_J.frame8)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=.5)

############ LG_J /\\ #######################
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# 104 individuals, 323 markers...
#
##
###
####
#####
######
############# ALL TOGETHER
# For PC
setwd("D:/PhD/Crossing/Analyses/Recomb.Map/Population/FINAL_DNA_PHENOTYPES")
getwd()

library(onemap)

# 104 individuals
LGs_1 <- read.outcross(,"OneWay104.txt")
LGs_2pt <- rf.2pts(LGs_1, LOD=8, max.rf=.4)
LGs_2pt_makeseq <- make.seq (LGs_2pt, "all")
marker.type(LGs_2pt_makeseq)
Grouped_LGs <- group(LGs_2pt_makeseq)
print(Grouped_LGs, detailed=FALSE)
print(Grouped_LGs)

set.map.fun(type="kosambi")

###### LG_A \\/ #########################
LG_A.frame1 <- make.seq(LGs_2pt, c(232,224,284))
LG_A.compare <-compare(LG_A.frame1)
LG_A.compare
LG_A.frame1 <- make.seq(LG_A.compare,1,1)
LG_A.frame1
### BACKBONE A.1 ### /\\ ####################

### IN BACKBONE \\/ ###
LG_A.frame_add224 <- try.seq(LG_A.frame2, 224, draw.try=TRUE)
LG_A.frame_add224 <- try.seq(LG_A.frame2, 224)
LG_A.frame_add224
LG_A.frame2 <- make.seq(LG_A.frame_add224,1,1)
LG_A.frame2
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add232 <- try.seq(LG_A.frame2, 232, draw.try=TRUE)
LG_A.frame_add232 <- try.seq(LG_A.frame2, 232)
LG_A.frame_add232
LG_A.frame2 <- make.seq(LG_A.frame_add232,1,1)
LG_A.frame2
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add284 <- try.seq(LG_A.frame1, 284, draw.try=TRUE)
LG_A.frame_add284 <- try.seq(LG_A.frame1, 284)
LG_A.frame_add284
LG_A.frame1 <- make.seq(LG_A.frame_add284,8,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])
### IN BACKBONE /\\ ###

## Frame 1 \\/ ##############################
LG_A.frame_add78 <- try.seq(LG_A.frame1, 78, draw.try=TRUE)
LG_A.frame_add78 <- try.seq(LG_A.frame1, 78)
LG_A.frame_add78
LG_A.frame1 <- make.seq(LG_A.frame_add78,4,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add144 <- try.seq(LG_A.frame1, 144, draw.try=TRUE)
LG_A.frame_add144 <- try.seq(LG_A.frame1, 144)
LG_A.frame_add144
LG_A.frame1 <- make.seq(LG_A.frame_add144,4,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add251 <- try.seq(LG_A.frame1, 251, draw.try=TRUE)
LG_A.frame_add251 <- try.seq(LG_A.frame1, 251)
LG_A.frame_add251
LG_A.frame1 <- make.seq(LG_A.frame_add251,1,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add243 <- try.seq(LG_A.frame1, 243, draw.try=TRUE)
LG_A.frame_add243 <- try.seq(LG_A.frame1, 243)
LG_A.frame_add243
LG_A.frame1 <- make.seq(LG_A.frame_add243,7,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add73 <- try.seq(LG_A.frame1, 73, draw.try=TRUE)
LG_A.frame_add73 <- try.seq(LG_A.frame1, 73)
LG_A.frame_add73
LG_A.frame1 <- make.seq(LG_A.frame_add73,7,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

# \\/ LESS THAN 80% OF INDIVIDUALS SCORED \\/ #
LG_A.frame_add134 <- try.seq(LG_A.frame1, 134, draw.try=TRUE)
LG_A.frame_add134 <- try.seq(LG_A.frame1, 134)
LG_A.frame_add134
LG_A.frame1 <- make.seq(LG_A.frame_add134,6,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add262 <- try.seq(LG_A.frame1, 262, draw.try=TRUE)
LG_A.frame_add262 <- try.seq(LG_A.frame1, 262)
LG_A.frame_add262
LG_A.frame1 <- make.seq(LG_A.frame_add262,4,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add239 <- try.seq(LG_A.frame1, 239, draw.try=TRUE)
LG_A.frame_add239 <- try.seq(LG_A.frame1, 239)
LG_A.frame_add239
LG_A.frame1 <- make.seq(LG_A.frame_add239,5,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add257 <- try.seq(LG_A.frame1, 257, draw.try=TRUE)
LG_A.frame_add257 <- try.seq(LG_A.frame1, 257)
LG_A.frame_add257
LG_A.frame1 <- make.seq(LG_A.frame_add257,5,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add110 <- try.seq(LG_A.frame1, 110, draw.try=TRUE)
LG_A.frame_add110 <- try.seq(LG_A.frame1, 110)
LG_A.frame_add110
LG_A.frame1 <- make.seq(LG_A.frame_add110,11,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add309 <- try.seq(LG_A.frame1, 309, draw.try=TRUE)
LG_A.frame_add309 <- try.seq(LG_A.frame1, 309)
LG_A.frame_add309
LG_A.frame1 <- make.seq(LG_A.frame_add309,3,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add98 <- try.seq(LG_A.frame1, 98, draw.try=TRUE)
LG_A.frame_add98 <- try.seq(LG_A.frame1, 98)
LG_A.frame_add98
LG_A.frame1 <- make.seq(LG_A.frame_add98,15,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add303 <- try.seq(LG_A.frame1, 303, draw.try=TRUE)
LG_A.frame_add303 <- try.seq(LG_A.frame1, 303)
LG_A.frame_add303
LG_A.frame1 <- make.seq(LG_A.frame_add303,15,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add219 <- try.seq(LG_A.frame1, 219, draw.try=TRUE)
LG_A.frame_add219 <- try.seq(LG_A.frame1, 219)
LG_A.frame_add219
LG_A.frame1 <- make.seq(LG_A.frame_add219,15,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

rf.graph.table(LG_A.frame1, scale=2)
ripple.seq(LG_A.frame1, ws=6, LOD=1, tol=0.1)
maps <- list(LG_A.frame1)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=1.5)

Master <- LG_A.frame1
## LG_A Frame 1 /\\ ################# TOTAL OF 98.27 CM ####
#
#
# Add markers from LG_B
LG_B.frame_add182 <- try.seq(Master, 182, draw.try=TRUE)
LG_B.frame_add182 <- try.seq(Master, 182)
LG_B.frame_add182
LG_B.frame2 <- make.seq(LG_B.frame_add182,18,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add190 <- try.seq(LG_B.frame2, 190, draw.try=TRUE)
LG_B.frame_add190 <- try.seq(LG_B.frame2, 190)
LG_B.frame_add190
LG_B.frame2 <- make.seq(LG_B.frame_add190,19,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add102 <- try.seq(LG_B.frame2, 102, draw.try=TRUE)
LG_B.frame_add102 <- try.seq(LG_B.frame2, 102)
LG_B.frame_add102
LG_B.frame2 <- make.seq(LG_B.frame_add102,20,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add1 <- try.seq(LG_B.frame2, 1, draw.try=TRUE)
LG_B.frame_add1 <- try.seq(LG_B.frame2, 1)
LG_B.frame_add1
LG_B.frame2 <- make.seq(LG_B.frame_add1,21,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add203 <- try.seq(LG_B.frame2, 203, draw.try=TRUE)
LG_B.frame_add203 <- try.seq(LG_B.frame2, 203)
LG_B.frame_add203
LG_B.frame2 <- make.seq(LG_B.frame_add203,22,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add181 <- try.seq(LG_B.frame2, 181, draw.try=TRUE)
LG_B.frame_add181 <- try.seq(LG_B.frame2, 181)
LG_B.frame_add181
LG_B.frame2 <- make.seq(LG_B.frame_add181,23,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add139 <- try.seq(LG_B.frame2, 139, draw.try=TRUE)
LG_B.frame_add139 <- try.seq(LG_B.frame2, 139)
LG_B.frame_add139
LG_B.frame2 <- make.seq(LG_B.frame_add139,24,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add178 <- try.seq(LG_B.frame2, 178, draw.try=TRUE)
LG_B.frame_add178 <- try.seq(LG_B.frame2, 178)
LG_B.frame_add178
LG_B.frame2 <- make.seq(LG_B.frame_add178,25,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add189 <- try.seq(LG_B.frame2, 189, draw.try=TRUE)
LG_B.frame_add189 <- try.seq(LG_B.frame2, 189)
LG_B.frame_add189
LG_B.frame2 <- make.seq(LG_B.frame_add189,26,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add206 <- try.seq(LG_B.frame2, 206, draw.try=TRUE)
LG_B.frame_add206 <- try.seq(LG_B.frame2, 206)
LG_B.frame_add206
LG_B.frame2 <- make.seq(LG_B.frame_add206,27,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

Master <- LG_B.frame2
#
rf.graph.table(Master, scale=2)

ripple.seq(Master, ws=6, LOD=1, tol=0.1)
maps <- list(Master)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=1.5)

## Add LG_C \\/ ##############################



LG_C.frame_add183 <- try.seq(Master, 183, draw.try=TRUE)
LG_C.frame_add183 <- try.seq(Master, 183)
LG_C.frame_add183
LG_C.frame3 <- make.seq(LG_C.frame_add183,28,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add107 <- try.seq(LG_C.frame3, 107, draw.try=TRUE)
LG_C.frame_add107 <- try.seq(LG_C.frame3, 107)
LG_C.frame_add107
LG_C.frame3 <- make.seq(LG_C.frame_add107,29,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add27 <- try.seq(LG_C.frame3, 27, draw.try=TRUE)
LG_C.frame_add27 <- try.seq(LG_C.frame3, 27)
LG_C.frame_add27
LG_C.frame3 <- make.seq(LG_C.frame_add27,30,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add301 <- try.seq(LG_C.frame3, 301, draw.try=TRUE)
LG_C.frame_add301 <- try.seq(LG_C.frame3, 301)
LG_C.frame_add301
LG_C.frame3 <- make.seq(LG_C.frame_add301,31,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add149 <- try.seq(LG_C.frame3, 149, draw.try=TRUE)
LG_C.frame_add149 <- try.seq(LG_C.frame3, 149)
LG_C.frame_add149
LG_C.frame3 <- make.seq(LG_C.frame_add149,32,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add302 <- try.seq(LG_C.frame3, 302, draw.try=TRUE)
LG_C.frame_add302 <- try.seq(LG_C.frame3, 302)
LG_C.frame_add302
LG_C.frame3 <- make.seq(LG_C.frame_add302,33,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add268 <- try.seq(LG_C.frame3, 268, draw.try=TRUE)
LG_C.frame_add268 <- try.seq(LG_C.frame3, 268)
LG_C.frame_add268
LG_C.frame3 <- make.seq(LG_C.frame_add268,34,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add271 <- try.seq(LG_C.frame3, 271, draw.try=TRUE)
LG_C.frame_add271 <- try.seq(LG_C.frame3, 271)
LG_C.frame_add271
LG_C.frame3 <- make.seq(LG_C.frame_add271,35,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add19 <- try.seq(LG_C.frame3, 19, draw.try=TRUE)
LG_C.frame_add19 <- try.seq(LG_C.frame3, 19)
LG_C.frame_add19
LG_C.frame3 <- make.seq(LG_C.frame_add19,36,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add123 <- try.seq(LG_C.frame3, 123, draw.try=TRUE)
LG_C.frame_add123 <- try.seq(LG_C.frame3, 123)
LG_C.frame_add123
LG_C.frame3 <- make.seq(LG_C.frame_add123,37,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add37 <- try.seq(LG_C.frame3, 37, draw.try=TRUE)
LG_C.frame_add37 <- try.seq(LG_C.frame3, 37)
LG_C.frame_add37
LG_C.frame3 <- make.seq(LG_C.frame_add37,38,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add273 <- try.seq(LG_C.frame3, 273, draw.try=TRUE)
LG_C.frame_add273 <- try.seq(LG_C.frame3, 273)
LG_C.frame_add273
LG_C.frame3 <- make.seq(LG_C.frame_add273,39,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add213 <- try.seq(LG_C.frame3, 213, draw.try=TRUE)
LG_C.frame_add213 <- try.seq(LG_C.frame3, 213)
LG_C.frame_add213
LG_C.frame3 <- make.seq(LG_C.frame_add213,40,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add306 <- try.seq(LG_C.frame3, 306, draw.try=TRUE)
LG_C.frame_add306 <- try.seq(LG_C.frame3, 306)
LG_C.frame_add306
LG_C.frame3 <- make.seq(LG_C.frame_add306,41,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add212 <- try.seq(LG_C.frame3, 212, draw.try=TRUE)
LG_C.frame_add212 <- try.seq(LG_C.frame3, 212)
LG_C.frame_add212
LG_C.frame3 <- make.seq(LG_C.frame_add212,42,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])


Master <- LG_C.frame3
rf.graph.table(Master, scale=2)

####################################

## ADD LG_D ##############################

LG_D.frame_add16 <- try.seq(Master, 16, draw.try=TRUE)
LG_D.frame_add16 <- try.seq(Master, 16)
LG_D.frame_add16
LG_D.frame4 <- make.seq(LG_D.frame_add16,43,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add26 <- try.seq(LG_D.frame4, 26, draw.try=TRUE)
LG_D.frame_add26 <- try.seq(LG_D.frame4, 26)
LG_D.frame_add26
LG_D.frame4 <- make.seq(LG_D.frame_add26,44,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add209 <- try.seq(LG_D.frame4, 209, draw.try=TRUE)
LG_D.frame_add209 <- try.seq(LG_D.frame4, 209)
LG_D.frame_add209
LG_D.frame4 <- make.seq(LG_D.frame_add209,45,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add267 <- try.seq(LG_D.frame4, 267, draw.try=TRUE)
LG_D.frame_add267 <- try.seq(LG_D.frame4, 267)
LG_D.frame_add267
LG_D.frame4 <- make.seq(LG_D.frame_add267,46,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add191 <- try.seq(LG_D.frame4, 191, draw.try=TRUE)
LG_D.frame_add191 <- try.seq(LG_D.frame4, 191)
LG_D.frame_add191
LG_D.frame4 <- make.seq(LG_D.frame_add191,47,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add187 <- try.seq(LG_D.frame4, 187, draw.try=TRUE)
LG_D.frame_add187 <- try.seq(LG_D.frame4, 187)
LG_D.frame_add187
LG_D.frame4 <- make.seq(LG_D.frame_add187,48,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add281 <- try.seq(LG_D.frame4, 281, draw.try=TRUE)
LG_D.frame_add281 <- try.seq(LG_D.frame4, 281)
LG_D.frame_add281
LG_D.frame4 <- make.seq(LG_D.frame_add281,49,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add221 <- try.seq(LG_D.frame4, 221, draw.try=TRUE)
LG_D.frame_add221 <- try.seq(LG_D.frame4, 221)
LG_D.frame_add221
LG_D.frame4 <- make.seq(LG_D.frame_add221,50,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add286 <- try.seq(LG_D.frame4, 286, draw.try=TRUE)
LG_D.frame_add286 <- try.seq(LG_D.frame4, 286)
LG_D.frame_add286
LG_D.frame4 <- make.seq(LG_D.frame_add286,51,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add21 <- try.seq(LG_D.frame4, 21, draw.try=TRUE)
LG_D.frame_add21 <- try.seq(LG_D.frame4, 21)
LG_D.frame_add21
LG_D.frame4 <- make.seq(LG_D.frame_add21,52,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add226 <- try.seq(LG_D.frame4, 226, draw.try=TRUE)
LG_D.frame_add226 <- try.seq(LG_D.frame4, 226)
LG_D.frame_add226
LG_D.frame4 <- make.seq(LG_D.frame_add226,53,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add280 <- try.seq(LG_D.frame4, 280, draw.try=TRUE)
LG_D.frame_add280 <- try.seq(LG_D.frame4, 280)
LG_D.frame_add280
LG_D.frame4 <- make.seq(LG_D.frame_add280,54,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add147 <- try.seq(LG_D.frame4, 147, draw.try=TRUE)
LG_D.frame_add147 <- try.seq(LG_D.frame4, 147)
LG_D.frame_add147
LG_D.frame4 <- make.seq(LG_D.frame_add147,55,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add188 <- try.seq(LG_D.frame4, 188, draw.try=TRUE)
LG_D.frame_add188 <- try.seq(LG_D.frame4, 188)
LG_D.frame_add188
LG_D.frame4 <- make.seq(LG_D.frame_add188,56,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add23 <- try.seq(LG_D.frame4, 23, draw.try=TRUE)
LG_D.frame_add23 <- try.seq(LG_D.frame4, 23)
LG_D.frame_add23
LG_D.frame4 <- make.seq(LG_D.frame_add23,57,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add202 <- try.seq(LG_D.frame4, 202, draw.try=TRUE)
LG_D.frame_add202 <- try.seq(LG_D.frame4, 202)
LG_D.frame_add202
LG_D.frame4 <- make.seq(LG_D.frame_add202,58,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add265 <- try.seq(LG_D.frame4, 265, draw.try=TRUE)
LG_D.frame_add265 <- try.seq(LG_D.frame4, 265)
LG_D.frame_add265
LG_D.frame4 <- make.seq(LG_D.frame_add265,59,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

Master <- LG_D.frame4

rf.graph.table(Master, scale=2)
#
# ADD LG_E
#LG_E.frame_add93 <- try.seq(Master, 93, draw.try=TRUE)
LG_E.frame_add93 <- try.seq(Master, 93)
LG_E.frame_add93
LG_E.frame5 <- make.seq(LG_E.frame_add93,60,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

#LG_E.frame_add300 <- try.seq(LG_E.frame5, 300, draw.try=TRUE)
LG_E.frame_add300 <- try.seq(LG_E.frame5, 300)
LG_E.frame_add300
LG_E.frame5 <- make.seq(LG_E.frame_add300,61,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

#LG_E.frame_add88 <- try.seq(LG_E.frame5, 88, draw.try=TRUE)
LG_E.frame_add88 <- try.seq(LG_E.frame5, 88)
LG_E.frame_add88
LG_E.frame5 <- make.seq(LG_E.frame_add88,62,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

#LG_E.frame_add263 <- try.seq(LG_E.frame5, 263, draw.try=TRUE)
LG_E.frame_add263 <- try.seq(LG_E.frame5, 263)
LG_E.frame_add263
LG_E.frame5 <- make.seq(LG_E.frame_add263,63,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

#LG_E.frame_add103 <- try.seq(LG_E.frame5, 103, draw.try=TRUE)
LG_E.frame_add103 <- try.seq(LG_E.frame5, 103)
LG_E.frame_add103
LG_E.frame5 <- make.seq(LG_E.frame_add103,64,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

#LG_E.frame_add151 <- try.seq(LG_E.frame5, 151, draw.try=TRUE)
LG_E.frame_add151 <- try.seq(LG_E.frame5, 151)
LG_E.frame_add151
LG_E.frame5 <- make.seq(LG_E.frame_add151,65,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

#LG_E.frame_add142 <- try.seq(LG_E.frame5, 142, draw.try=TRUE)
LG_E.frame_add142 <- try.seq(LG_E.frame5, 142)
LG_E.frame_add142
LG_E.frame5 <- make.seq(LG_E.frame_add142,66,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

#LG_E.frame_add106 <- try.seq(LG_E.frame5, 106, draw.try=TRUE)
LG_E.frame_add106 <- try.seq(LG_E.frame5, 106)
LG_E.frame_add106
LG_E.frame5 <- make.seq(LG_E.frame_add106,67,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

#LG_E.frame_add80 <- try.seq(LG_E.frame5, 80, draw.try=TRUE)
LG_E.frame_add80 <- try.seq(LG_E.frame5, 80)
LG_E.frame_add80
LG_E.frame5 <- make.seq(LG_E.frame_add80,68,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

#LG_E.frame_add261 <- try.seq(LG_E.frame5, 261, draw.try=TRUE)
LG_E.frame_add261 <- try.seq(LG_E.frame5, 261)
LG_E.frame_add261
LG_E.frame5 <- make.seq(LG_E.frame_add261,69,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

Master <- LG_E.frame5
rf.graph.table(Master, scale=2)

#
# ADD LG_F
#LG_F.frame_add48 <- try.seq(Master, 48, draw.try=TRUE)
LG_F.frame_add48 <- try.seq(Master, 48)
LG_F.frame_add48
LG_F.frame6 <- make.seq(LG_F.frame_add48,70,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

#LG_F.frame_add290 <- try.seq(LG_F.frame6, 290, draw.try=TRUE)
LG_F.frame_add290 <- try.seq(LG_F.frame6, 290)
LG_F.frame_add290
LG_F.frame6 <- make.seq(LG_F.frame_add290,71,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

#LG_F.frame_add56 <- try.seq(LG_F.frame6, 56, draw.try=TRUE)
LG_F.frame_add56 <- try.seq(LG_F.frame6, 56)
LG_F.frame_add56
LG_F.frame6 <- make.seq(LG_F.frame_add56,72,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

#LG_F.frame_add246 <- try.seq(LG_F.frame6, 246, draw.try=TRUE)
LG_F.frame_add246 <- try.seq(LG_F.frame6, 246)
LG_F.frame_add246
LG_F.frame6 <- make.seq(LG_F.frame_add246,73,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

#LG_F.frame_add248 <- try.seq(LG_F.frame6, 248, draw.try=TRUE)
LG_F.frame_add248 <- try.seq(LG_F.frame6, 248)
LG_F.frame_add248
LG_F.frame6 <- make.seq(LG_F.frame_add248,74,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

#LG_F.frame_add242 <- try.seq(LG_F.frame6, 242, draw.try=TRUE)
LG_F.frame_add242 <- try.seq(LG_F.frame6, 242)
LG_F.frame_add242
LG_F.frame6 <- make.seq(LG_F.frame_add242,75,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

#LG_F.frame_add274 <- try.seq(LG_F.frame6, 274, draw.try=TRUE)
LG_F.frame_add274 <- try.seq(LG_F.frame6, 274)
LG_F.frame_add274
LG_F.frame6 <- make.seq(LG_F.frame_add274,76,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

#LG_F.frame_add62 <- try.seq(LG_F.frame6, 62, draw.try=TRUE)
LG_F.frame_add62 <- try.seq(LG_F.frame6, 62)
LG_F.frame_add62
LG_F.frame6 <- make.seq(LG_F.frame_add62,77,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

#LG_F.frame_add253 <- try.seq(LG_F.frame6, 253, draw.try=TRUE)
LG_F.frame_add253 <- try.seq(LG_F.frame6, 253)
LG_F.frame_add253
LG_F.frame6 <- make.seq(LG_F.frame_add253,78,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

#LG_F.frame_add91 <- try.seq(LG_F.frame6, 91, draw.try=TRUE)
LG_F.frame_add91 <- try.seq(LG_F.frame6, 91)
LG_F.frame_add91
LG_F.frame6 <- make.seq(LG_F.frame_add91,79,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

#LG_F.frame_add9 <- try.seq(LG_F.frame6, 9, draw.try=TRUE)
LG_F.frame_add9 <- try.seq(LG_F.frame6, 9)
LG_F.frame_add9
LG_F.frame6 <- make.seq(LG_F.frame_add9,80,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

Master <- LG_F.frame6
rf.graph.table(LG_F.frame6, scale=2)


#
# ADD LG_G
#LG_G.frame_add84 <- try.seq(Master, 84, draw.try=TRUE)
LG_G.frame_add84 <- try.seq(Master, 84)
LG_G.frame_add84
LG_G.frame7 <- make.seq(LG_G.frame_add84,81,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

#LG_G.frame_add76 <- try.seq(LG_G.frame7, 76, draw.try=TRUE)
LG_G.frame_add76 <- try.seq(LG_G.frame7, 76)
LG_G.frame_add76
LG_G.frame7 <- make.seq(LG_G.frame_add76,82,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

#LG_G.frame_add140 <- try.seq(LG_G.frame7, 140, draw.try=TRUE)
LG_G.frame_add140 <- try.seq(LG_G.frame7, 140)
LG_G.frame_add140
LG_G.frame7 <- make.seq(LG_G.frame_add140,83,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

#LG_G.frame_add20 <- try.seq(LG_G.frame7, 20, draw.try=TRUE)
LG_G.frame_add20 <- try.seq(LG_G.frame7, 20)
LG_G.frame_add20
LG_G.frame7 <- make.seq(LG_G.frame_add20,84,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

#LG_G.frame_add193 <- try.seq(LG_G.frame7, 193, draw.try=TRUE)
LG_G.frame_add193 <- try.seq(LG_G.frame7, 193)
LG_G.frame_add193
LG_G.frame7 <- make.seq(LG_G.frame_add193,85,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

#LG_G.frame_add117 <- try.seq(LG_G.frame7, 117, draw.try=TRUE)
LG_G.frame_add117 <- try.seq(LG_G.frame7, 117)
LG_G.frame_add117
LG_G.frame7 <- make.seq(LG_G.frame_add117,86,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

#LG_G.frame_add18 <- try.seq(LG_G.frame7, 18, draw.try=TRUE)
LG_G.frame_add18 <- try.seq(LG_G.frame7, 18)
LG_G.frame_add18
LG_G.frame7 <- make.seq(LG_G.frame_add18,87,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

#LG_G.frame_add217 <- try.seq(LG_G.frame7, 217, draw.try=TRUE)
LG_G.frame_add217 <- try.seq(LG_G.frame7, 217)
LG_G.frame_add217
LG_G.frame7 <- make.seq(LG_G.frame_add217,88,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

#LG_G.frame_add15 <- try.seq(LG_G.frame7, 15, draw.try=TRUE)
LG_G.frame_add15 <- try.seq(LG_G.frame7, 15)
LG_G.frame_add15
LG_G.frame7 <- make.seq(LG_G.frame_add15,89,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

#LG_G.frame_add304 <- try.seq(LG_G.frame7, 304, draw.try=TRUE)
LG_G.frame_add304 <- try.seq(LG_G.frame7, 304)
LG_G.frame_add304
LG_G.frame7 <- make.seq(LG_G.frame_add304,90,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

#LG_G.frame_add59 <- try.seq(LG_G.frame7, 59, draw.try=TRUE)
LG_G.frame_add59 <- try.seq(LG_G.frame7, 59)
LG_G.frame_add59
LG_G.frame7 <- make.seq(LG_G.frame_add59,91,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

Master <- LG_G.frame7

rf.graph.table(Master, scale=3)

#
############# Add LG_J ## # ###########

#LG_J.frame_add13 <- try.seq(Master, 13, draw.try=TRUE)
LG_J.frame_add13 <- try.seq(Master, 13)
LG_J.frame_add13
LG_J.frame8 <- make.seq(LG_J.frame_add13,92,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

#LG_J.frame_add196 <- try.seq(LG_J.frame8, 196, draw.try=TRUE)
LG_J.frame_add196 <- try.seq(LG_J.frame8, 196)
LG_J.frame_add196
LG_J.frame8 <- make.seq(LG_J.frame_add196,93,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

#LG_J.frame_add223 <- try.seq(LG_J.frame8, 223, draw.try=TRUE)
LG_J.frame_add223 <- try.seq(LG_J.frame8, 223)
LG_J.frame_add223
LG_J.frame8 <- make.seq(LG_J.frame_add223,94,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

#LG_J.frame_add200 <- try.seq(LG_J.frame8, 200, draw.try=TRUE)
LG_J.frame_add200 <- try.seq(LG_J.frame8, 200)
LG_J.frame_add200
LG_J.frame8 <- make.seq(LG_J.frame_add200,95,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

#LG_J.frame_add174 <- try.seq(LG_J.frame8, 174, draw.try=TRUE)
LG_J.frame_add174 <- try.seq(LG_J.frame8, 174)
LG_J.frame_add174
LG_J.frame8 <- make.seq(LG_J.frame_add174,96,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

#LG_J.frame_add185 <- try.seq(LG_J.frame8, 185, draw.try=TRUE)
LG_J.frame_add185 <- try.seq(LG_J.frame8, 185)
LG_J.frame_add185
LG_J.frame8 <- make.seq(LG_J.frame_add185,97,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

#LG_J.frame_add220 <- try.seq(LG_J.frame8, 220, draw.try=TRUE)
LG_J.frame_add220 <- try.seq(LG_J.frame8, 220)
LG_J.frame_add220
LG_J.frame8 <- make.seq(LG_J.frame_add220,98,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

#LG_J.frame_add29 <- try.seq(LG_J.frame8, 29, draw.try=TRUE)
LG_J.frame_add29 <- try.seq(LG_J.frame8, 29)
LG_J.frame_add29
LG_J.frame8 <- make.seq(LG_J.frame_add29,99,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

Master <- LG_J.frame8

rf.graph.table(Master, scale=3)
ripple.seq(LG_J.frame8, ws=6, LOD=1, tol=0.1)
maps <- list(LG_J.frame8)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=.5)

############ LG_J /\\ #######################
# Add LG_H \\/ ###################

#LG_H.frame_add116 <- try.seq(Master, 116, draw.try=TRUE)
LG_H.frame_add116 <- try.seq(Master, 116)
LG_H.frame_add116
LG_H.frame8 <- make.seq(LG_H.frame_add116,100,1)
LG_H.frame8
dev.off(dev.list()["RStudioGD"])

#LG_H.frame_add228 <- try.seq(LG_H.frame8, 228, draw.try=TRUE)
LG_H.frame_add228 <- try.seq(LG_H.frame8, 228)
LG_H.frame_add228
LG_H.frame8 <- make.seq(LG_H.frame_add228,101,1)
LG_H.frame8
dev.off(dev.list()["RStudioGD"])

#LG_H.frame_add133 <- try.seq(LG_H.frame8, 133, draw.try=TRUE)
LG_H.frame_add133 <- try.seq(LG_H.frame8, 133)
LG_H.frame_add133
LG_H.frame8 <- make.seq(LG_H.frame_add133,102,1)
LG_H.frame8
dev.off(dev.list()["RStudioGD"])

#LG_H.frame_add154 <- try.seq(LG_H.frame8, 154, draw.try=TRUE)
LG_H.frame_add154 <- try.seq(LG_H.frame8, 154)
LG_H.frame_add154
LG_H.frame8 <- make.seq(LG_H.frame_add154,103,1)
LG_H.frame8
dev.off(dev.list()["RStudioGD"])

#LG_H.frame_add95 <- try.seq(LG_H.frame8, 95, draw.try=TRUE)
LG_H.frame_add95 <- try.seq(LG_H.frame8, 95)
LG_H.frame_add95
LG_H.frame8 <- make.seq(LG_H.frame_add95,104,1)
LG_H.frame8
dev.off(dev.list()["RStudioGD"])

#LG_H.frame_add66 <- try.seq(LG_H.frame8, 66, draw.try=TRUE)
LG_H.frame_add66 <- try.seq(LG_H.frame8, 66)
LG_H.frame_add66
LG_H.frame8 <- make.seq(LG_H.frame_add66,105,1)
LG_H.frame8
dev.off(dev.list()["RStudioGD"])


Master <- LG_H.frame8 
rf.graph.table(Master, scale=3)
#

# ADD LG_I ######################
#LG_I.frame_add105 <- try.seq(Master, 105, draw.try=TRUE)
LG_I.frame_add105 <- try.seq(Master, 105)
LG_I.frame_add105
LG_I.frame8 <- make.seq(LG_I.frame_add105,106,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])

#LG_I.frame_add236 <- try.seq(LG_I.frame8, 236, draw.try=TRUE)
LG_I.frame_add236 <- try.seq(LG_I.frame8, 236)
LG_I.frame_add236
LG_I.frame8 <- make.seq(LG_I.frame_add236,107,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])

#LG_I.frame_add119 <- try.seq(LG_I.frame8, 119, draw.try=TRUE)
LG_I.frame_add119 <- try.seq(LG_I.frame8, 119)
LG_I.frame_add119
LG_I.frame8 <- make.seq(LG_I.frame_add119,108,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])

#LG_I.frame_add125 <- try.seq(LG_I.frame8, 125, draw.try=TRUE)
LG_I.frame_add125 <- try.seq(LG_I.frame8, 125)
LG_I.frame_add125
LG_I.frame8 <- make.seq(LG_I.frame_add125,109,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])

#LG_I.frame_add244 <- try.seq(LG_I.frame8, 244, draw.try=TRUE)
LG_I.frame_add244 <- try.seq(LG_I.frame8, 244)
LG_I.frame_add244
LG_I.frame8 <- make.seq(LG_I.frame_add244,110,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])

#LG_I.frame_add310 <- try.seq(LG_I.frame8, 310, draw.try=TRUE)
LG_I.frame_add310 <- try.seq(LG_I.frame8, 310)
LG_I.frame_add310
LG_I.frame8 <- make.seq(LG_I.frame_add310,111,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])

#LG_I.frame_add44 <- try.seq(LG_I.frame8, 44, draw.try=TRUE)
LG_I.frame_add44 <- try.seq(LG_I.frame8, 44)
LG_I.frame_add44
LG_I.frame8 <- make.seq(LG_I.frame_add44,112,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])

Master <- LG_I.frame8 
rf.graph.table(Master, scale=3)





#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# 104 individuals, 323 markers...
#
##
###
####
#####
######
############# ALL SEPARATE
# For PC
setwd("D:/PhD/Crossing/Analyses/Recomb.Map/Population/FINAL_DNA_PHENOTYPES")
getwd()

library(onemap)

# 104 individuals
LGs_1 <- read.outcross(,"OneWay104.txt")
LGs_2pt <- rf.2pts(LGs_1, LOD=8, max.rf=.4)
LGs_2pt_makeseq <- make.seq (LGs_2pt, "all")
marker.type(LGs_2pt_makeseq)
Grouped_LGs <- group(LGs_2pt_makeseq)
print(Grouped_LGs, detailed=FALSE)
print(Grouped_LGs)

set.map.fun(type="kosambi")

###### LG_A \\/ #########################
LG_A.frame1 <- make.seq(LGs_2pt, c(232,224,284))
LG_A.compare <-compare(LG_A.frame1)
LG_A.compare
LG_A.frame1 <- make.seq(LG_A.compare,1,1)
LG_A.frame1
### BACKBONE A.1 ### /\\ ####################

### IN BACKBONE \\/ ###
LG_A.frame_add224 <- try.seq(LG_A.frame2, 224, draw.try=TRUE)
LG_A.frame_add224 <- try.seq(LG_A.frame2, 224)
LG_A.frame_add224
LG_A.frame2 <- make.seq(LG_A.frame_add224,1,1)
LG_A.frame2
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add232 <- try.seq(LG_A.frame2, 232, draw.try=TRUE)
LG_A.frame_add232 <- try.seq(LG_A.frame2, 232)
LG_A.frame_add232
LG_A.frame2 <- make.seq(LG_A.frame_add232,1,1)
LG_A.frame2
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add284 <- try.seq(LG_A.frame1, 284, draw.try=TRUE)
LG_A.frame_add284 <- try.seq(LG_A.frame1, 284)
LG_A.frame_add284
LG_A.frame1 <- make.seq(LG_A.frame_add284,8,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])
### IN BACKBONE /\\ ###

## Frame 1 \\/ ##############################
LG_A.frame_add78 <- try.seq(LG_A.frame1, 78, draw.try=TRUE)
LG_A.frame_add78 <- try.seq(LG_A.frame1, 78)
LG_A.frame_add78
LG_A.frame1 <- make.seq(LG_A.frame_add78,4,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add144 <- try.seq(LG_A.frame1, 144, draw.try=TRUE)
LG_A.frame_add144 <- try.seq(LG_A.frame1, 144)
LG_A.frame_add144
LG_A.frame1 <- make.seq(LG_A.frame_add144,4,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add251 <- try.seq(LG_A.frame1, 251, draw.try=TRUE)
LG_A.frame_add251 <- try.seq(LG_A.frame1, 251)
LG_A.frame_add251
LG_A.frame1 <- make.seq(LG_A.frame_add251,1,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add243 <- try.seq(LG_A.frame1, 243, draw.try=TRUE)
LG_A.frame_add243 <- try.seq(LG_A.frame1, 243)
LG_A.frame_add243
LG_A.frame1 <- make.seq(LG_A.frame_add243,7,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add73 <- try.seq(LG_A.frame1, 73, draw.try=TRUE)
LG_A.frame_add73 <- try.seq(LG_A.frame1, 73)
LG_A.frame_add73
LG_A.frame1 <- make.seq(LG_A.frame_add73,7,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

# \\/ LESS THAN 80% OF INDIVIDUALS SCORED \\/ #
LG_A.frame_add134 <- try.seq(LG_A.frame1, 134, draw.try=TRUE)
LG_A.frame_add134 <- try.seq(LG_A.frame1, 134)
LG_A.frame_add134
LG_A.frame1 <- make.seq(LG_A.frame_add134,6,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add262 <- try.seq(LG_A.frame1, 262, draw.try=TRUE)
LG_A.frame_add262 <- try.seq(LG_A.frame1, 262)
LG_A.frame_add262
LG_A.frame1 <- make.seq(LG_A.frame_add262,4,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add239 <- try.seq(LG_A.frame1, 239, draw.try=TRUE)
LG_A.frame_add239 <- try.seq(LG_A.frame1, 239)
LG_A.frame_add239
LG_A.frame1 <- make.seq(LG_A.frame_add239,5,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add257 <- try.seq(LG_A.frame1, 257, draw.try=TRUE)
LG_A.frame_add257 <- try.seq(LG_A.frame1, 257)
LG_A.frame_add257
LG_A.frame1 <- make.seq(LG_A.frame_add257,5,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add110 <- try.seq(LG_A.frame1, 110, draw.try=TRUE)
LG_A.frame_add110 <- try.seq(LG_A.frame1, 110)
LG_A.frame_add110
LG_A.frame1 <- make.seq(LG_A.frame_add110,11,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add309 <- try.seq(LG_A.frame1, 309, draw.try=TRUE)
LG_A.frame_add309 <- try.seq(LG_A.frame1, 309)
LG_A.frame_add309
LG_A.frame1 <- make.seq(LG_A.frame_add309,3,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add98 <- try.seq(LG_A.frame1, 98, draw.try=TRUE)
LG_A.frame_add98 <- try.seq(LG_A.frame1, 98)
LG_A.frame_add98
LG_A.frame1 <- make.seq(LG_A.frame_add98,15,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add303 <- try.seq(LG_A.frame1, 303, draw.try=TRUE)
LG_A.frame_add303 <- try.seq(LG_A.frame1, 303)
LG_A.frame_add303
LG_A.frame1 <- make.seq(LG_A.frame_add303,15,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])

LG_A.frame_add219 <- try.seq(LG_A.frame1, 219, draw.try=TRUE)
LG_A.frame_add219 <- try.seq(LG_A.frame1, 219)
LG_A.frame_add219
LG_A.frame1 <- make.seq(LG_A.frame_add219,15,1)
LG_A.frame1
dev.off(dev.list()["RStudioGD"])


## LG_A Frame 1 /\\ ################# TOTAL OF 98.27 CM ####
#
#
# Add markers from LG_B
LG_B.frame2 <- make.seq(LGs_2pt, c(182,190,102))
LG_B.compare <-compare(LG_B.frame2)
LG_B.compare
LG_B.frame2 <- make.seq(LG_B.compare,2,1)
LG_B.frame2

## BACKBONE ########## ## # \\/ ## #################################
LG_B.frame_add182 <- try.seq(LG_B.frame2, 182, draw.try=TRUE)
LG_B.frame_add182 <- try.seq(LG_B.frame2, 182)
LG_B.frame_add182
LG_B.frame2 <- make.seq(LG_B.frame_add182,1,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add190 <- try.seq(LG_B.frame2, 190, draw.try=TRUE)
LG_B.frame_add190 <- try.seq(LG_B.frame2, 190)
LG_B.frame_add190
LG_B.frame2 <- make.seq(LG_B.frame_add190,19,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add102 <- try.seq(LG_B.frame2, 102, draw.try=TRUE)
LG_B.frame_add102 <- try.seq(LG_B.frame2, 102)
LG_B.frame_add102
LG_B.frame2 <- make.seq(LG_B.frame_add102,20,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])
## BACKBONE ########## ## # /\\ ## #############################

LG_B.frame_add1 <- try.seq(LG_B.frame2, 1, draw.try=TRUE)
LG_B.frame_add1 <- try.seq(LG_B.frame2, 1)
LG_B.frame_add1
LG_B.frame2 <- make.seq(LG_B.frame_add1,4,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add203 <- try.seq(LG_B.frame2, 203, draw.try=TRUE)
LG_B.frame_add203 <- try.seq(LG_B.frame2, 203)
LG_B.frame_add203
LG_B.frame2 <- make.seq(LG_B.frame_add203,5,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add181 <- try.seq(LG_B.frame2, 181, draw.try=TRUE)
LG_B.frame_add181 <- try.seq(LG_B.frame2, 181)
LG_B.frame_add181
LG_B.frame2 <- make.seq(LG_B.frame_add181,6,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add139 <- try.seq(LG_B.frame2, 139, draw.try=TRUE)
LG_B.frame_add139 <- try.seq(LG_B.frame2, 139)
LG_B.frame_add139
LG_B.frame2 <- make.seq(LG_B.frame_add139,7,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add178 <- try.seq(LG_B.frame2, 178, draw.try=TRUE)
LG_B.frame_add178 <- try.seq(LG_B.frame2, 178)
LG_B.frame_add178
LG_B.frame2 <- make.seq(LG_B.frame_add178,8,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add189 <- try.seq(LG_B.frame2, 189, draw.try=TRUE)
LG_B.frame_add189 <- try.seq(LG_B.frame2, 189)
LG_B.frame_add189
LG_B.frame2 <- make.seq(LG_B.frame_add189,9,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

LG_B.frame_add206 <- try.seq(LG_B.frame2, 206, draw.try=TRUE)
LG_B.frame_add206 <- try.seq(LG_B.frame2, 206)
LG_B.frame_add206
LG_B.frame2 <- make.seq(LG_B.frame_add206,10,1)
LG_B.frame2
dev.off(dev.list()["RStudioGD"])

maps <- list(LG_B.frame2)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=1)
###############################################

## Add LG_C \\/ ##############################
LG_C.frame3 <- make.seq(LGs_2pt, c(212,306,213))
LG_C.compare <-compare(LG_C.frame3)
LG_C.compare
LG_C.frame3 <- make.seq(LG_C.compare,1,1)
LG_C.frame3

## BACKBONE ########## ## # \\/ ##   ###################################
LG_C.frame_add213 <- try.seq(LG_C.frame3, 213, draw.try=TRUE)
LG_C.frame_add213 <- try.seq(LG_C.frame3, 213)
LG_C.frame_add213
LG_C.frame3 <- make.seq(LG_C.frame_add213,40,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add306 <- try.seq(LG_C.frame3, 306, draw.try=TRUE)
LG_C.frame_add306 <- try.seq(LG_C.frame3, 306)
LG_C.frame_add306
LG_C.frame3 <- make.seq(LG_C.frame_add306,41,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add212 <- try.seq(LG_C.frame3, 212, draw.try=TRUE)
LG_C.frame_add212 <- try.seq(LG_C.frame3, 212)
LG_C.frame_add212
LG_C.frame3 <- make.seq(LG_C.frame_add212,42,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])
## BACKBONE ########## ## # /\\ ########################################

LG_C.frame_add273 <- try.seq(LG_C.frame3, 273, draw.try=TRUE)
LG_C.frame_add273 <- try.seq(LG_C.frame3, 273)
LG_C.frame_add273
LG_C.frame3 <- make.seq(LG_C.frame_add273,4,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add37 <- try.seq(LG_C.frame3, 37, draw.try=TRUE)
LG_C.frame_add37 <- try.seq(LG_C.frame3, 37)
LG_C.frame_add37
LG_C.frame3 <- make.seq(LG_C.frame_add37,5,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add123 <- try.seq(LG_C.frame3, 123, draw.try=TRUE)
LG_C.frame_add123 <- try.seq(LG_C.frame3, 123)
LG_C.frame_add123
LG_C.frame3 <- make.seq(LG_C.frame_add123,6,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add19 <- try.seq(LG_C.frame3, 19, draw.try=TRUE)
LG_C.frame_add19 <- try.seq(LG_C.frame3, 19)
LG_C.frame_add19
LG_C.frame3 <- make.seq(LG_C.frame_add19,7,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add271 <- try.seq(LG_C.frame3, 271, draw.try=TRUE)
LG_C.frame_add271 <- try.seq(LG_C.frame3, 271)
LG_C.frame_add271
LG_C.frame3 <- make.seq(LG_C.frame_add271,8,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add268 <- try.seq(LG_C.frame3, 268, draw.try=TRUE)
LG_C.frame_add268 <- try.seq(LG_C.frame3, 268)
LG_C.frame_add268
LG_C.frame3 <- make.seq(LG_C.frame_add268,9,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add302 <- try.seq(LG_C.frame3, 302, draw.try=TRUE)
LG_C.frame_add302 <- try.seq(LG_C.frame3, 302)
LG_C.frame_add302
LG_C.frame3 <- make.seq(LG_C.frame_add302,10,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add149 <- try.seq(LG_C.frame3, 149, draw.try=TRUE)
LG_C.frame_add149 <- try.seq(LG_C.frame3, 149)
LG_C.frame_add149
LG_C.frame3 <- make.seq(LG_C.frame_add149,11,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add301 <- try.seq(LG_C.frame3, 301, draw.try=TRUE)
LG_C.frame_add301 <- try.seq(LG_C.frame3, 301)
LG_C.frame_add301
LG_C.frame3 <- make.seq(LG_C.frame_add301,12,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add27 <- try.seq(LG_C.frame3, 27, draw.try=TRUE)
LG_C.frame_add27 <- try.seq(LG_C.frame3, 27)
LG_C.frame_add27
LG_C.frame3 <- make.seq(LG_C.frame_add27,13,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add107 <- try.seq(LG_C.frame3, 107, draw.try=TRUE)
LG_C.frame_add107 <- try.seq(LG_C.frame3, 107)
LG_C.frame_add107
LG_C.frame3 <- make.seq(LG_C.frame_add107,14,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])

LG_C.frame_add183 <- try.seq(LG_C.frame3, 183, draw.try=TRUE)
LG_C.frame_add183 <- try.seq(LG_C.frame3, 183)
LG_C.frame_add183
LG_C.frame3 <- make.seq(LG_C.frame_add183,15,1)
LG_C.frame3
dev.off(dev.list()["RStudioGD"])


####################################

## ADD LG_D ##############################

LG_D.frame4 <- make.seq(LGs_2pt, c(265,202,23))
LG_D.compare <-compare(LG_D.frame4)
LG_D.compare
LG_D.frame4 <- make.seq(LG_D.compare,1,1)
LG_D.frame4

## BACKBONE ########## ## # \\/ ##   ###################################

LG_D.frame_add23 <- try.seq(LG_D.frame4, 23, draw.try=TRUE)
LG_D.frame_add23 <- try.seq(LG_D.frame4, 23)
LG_D.frame_add23
LG_D.frame4 <- make.seq(LG_D.frame_add23,58,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add202 <- try.seq(LG_D.frame4, 202, draw.try=TRUE)
LG_D.frame_add202 <- try.seq(LG_D.frame4, 202)
LG_D.frame_add202
LG_D.frame4 <- make.seq(LG_D.frame_add202,59,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add265 <- try.seq(LG_D.frame4, 265, draw.try=TRUE)
LG_D.frame_add265 <- try.seq(LG_D.frame4, 265)
LG_D.frame_add265
LG_D.frame4 <- make.seq(LG_D.frame_add265,60,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

## BACKBONE ########## ## # /\\ ########################################

LG_D.frame_add188 <- try.seq(LG_D.frame4, 188, draw.try=TRUE)
LG_D.frame_add188 <- try.seq(LG_D.frame4, 188)
LG_D.frame_add188
LG_D.frame4 <- make.seq(LG_D.frame_add188,4,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add147 <- try.seq(LG_D.frame4, 147, draw.try=TRUE)
LG_D.frame_add147 <- try.seq(LG_D.frame4, 147)
LG_D.frame_add147
LG_D.frame4 <- make.seq(LG_D.frame_add147,5,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add280 <- try.seq(LG_D.frame4, 280, draw.try=TRUE)
LG_D.frame_add280 <- try.seq(LG_D.frame4, 280)
LG_D.frame_add280
LG_D.frame4 <- make.seq(LG_D.frame_add280,6,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add226 <- try.seq(LG_D.frame4, 226, draw.try=TRUE)
LG_D.frame_add226 <- try.seq(LG_D.frame4, 226)
LG_D.frame_add226
LG_D.frame4 <- make.seq(LG_D.frame_add226,7,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add21 <- try.seq(LG_D.frame4, 21, draw.try=TRUE)
LG_D.frame_add21 <- try.seq(LG_D.frame4, 21)
LG_D.frame_add21
LG_D.frame4 <- make.seq(LG_D.frame_add21,8,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add286 <- try.seq(LG_D.frame4, 286, draw.try=TRUE)
LG_D.frame_add286 <- try.seq(LG_D.frame4, 286)
LG_D.frame_add286
LG_D.frame4 <- make.seq(LG_D.frame_add286,9,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add221 <- try.seq(LG_D.frame4, 221, draw.try=TRUE)
LG_D.frame_add221 <- try.seq(LG_D.frame4, 221)
LG_D.frame_add221
LG_D.frame4 <- make.seq(LG_D.frame_add221,10,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add281 <- try.seq(LG_D.frame4, 281, draw.try=TRUE)
LG_D.frame_add281 <- try.seq(LG_D.frame4, 281)
LG_D.frame_add281
LG_D.frame4 <- make.seq(LG_D.frame_add281,11,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add187 <- try.seq(LG_D.frame4, 187, draw.try=TRUE)
LG_D.frame_add187 <- try.seq(LG_D.frame4, 187)
LG_D.frame_add187
LG_D.frame4 <- make.seq(LG_D.frame_add187,12,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add191 <- try.seq(LG_D.frame4, 191, draw.try=TRUE)
LG_D.frame_add191 <- try.seq(LG_D.frame4, 191)
LG_D.frame_add191
LG_D.frame4 <- make.seq(LG_D.frame_add191,13,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add267 <- try.seq(LG_D.frame4, 267, draw.try=TRUE)
LG_D.frame_add267 <- try.seq(LG_D.frame4, 267)
LG_D.frame_add267
LG_D.frame4 <- make.seq(LG_D.frame_add267,14,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add209 <- try.seq(LG_D.frame4, 209, draw.try=TRUE)
LG_D.frame_add209 <- try.seq(LG_D.frame4, 209)
LG_D.frame_add209
LG_D.frame4 <- make.seq(LG_D.frame_add209,15,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add26 <- try.seq(LG_D.frame4, 26, draw.try=TRUE)
LG_D.frame_add26 <- try.seq(LG_D.frame4, 26)
LG_D.frame_add26
LG_D.frame4 <- make.seq(LG_D.frame_add26,16,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

LG_D.frame_add16 <- try.seq(LG_D.frame4, 16, draw.try=TRUE)
LG_D.frame_add16 <- try.seq(LG_D.frame4, 16)
LG_D.frame_add16
LG_D.frame4 <- make.seq(LG_D.frame_add16,17,1)
LG_D.frame4
dev.off(dev.list()["RStudioGD"])

maps <- list(LG_D.frame4)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=1)

#
#
#
#
#
#

# ADD LG_E


LG_E.frame5 <- make.seq(LGs_2pt, c(93,300,88))
LG_E.compare <-compare(LG_E.frame5)
LG_E.compare
LG_E.frame5 <- make.seq(LG_E.compare,1,1)
LG_E.frame5

## BACKBONE ########## ## # \\/ ##   ###################################

LG_E.frame_add93 <- try.seq(Master, 93, draw.try=TRUE)
LG_E.frame_add93 <- try.seq(Master, 93)
LG_E.frame_add93
LG_E.frame5 <- make.seq(LG_E.frame_add93,61,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

LG_E.frame_add300 <- try.seq(LG_E.frame5, 300, draw.try=TRUE)
LG_E.frame_add300 <- try.seq(LG_E.frame5, 300)
LG_E.frame_add300
LG_E.frame5 <- make.seq(LG_E.frame_add300,62,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

LG_E.frame_add88 <- try.seq(LG_E.frame5, 88, draw.try=TRUE)
LG_E.frame_add88 <- try.seq(LG_E.frame5, 88)
LG_E.frame_add88
LG_E.frame5 <- make.seq(LG_E.frame_add88,63,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

## BACKBONE ########## ## # /\\ ########################################


LG_E.frame_add263 <- try.seq(LG_E.frame5, 263, draw.try=TRUE)
LG_E.frame_add263 <- try.seq(LG_E.frame5, 263)
LG_E.frame_add263
LG_E.frame5 <- make.seq(LG_E.frame_add263,4,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

LG_E.frame_add103 <- try.seq(LG_E.frame5, 103, draw.try=TRUE)
LG_E.frame_add103 <- try.seq(LG_E.frame5, 103)
LG_E.frame_add103
LG_E.frame5 <- make.seq(LG_E.frame_add103,5,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

LG_E.frame_add151 <- try.seq(LG_E.frame5, 151, draw.try=TRUE)
LG_E.frame_add151 <- try.seq(LG_E.frame5, 151)
LG_E.frame_add151
LG_E.frame5 <- make.seq(LG_E.frame_add151,6,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

LG_E.frame_add142 <- try.seq(LG_E.frame5, 142, draw.try=TRUE)
LG_E.frame_add142 <- try.seq(LG_E.frame5, 142)
LG_E.frame_add142
LG_E.frame5 <- make.seq(LG_E.frame_add142,7,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

LG_E.frame_add106 <- try.seq(LG_E.frame5, 106, draw.try=TRUE)
LG_E.frame_add106 <- try.seq(LG_E.frame5, 106)
LG_E.frame_add106
LG_E.frame5 <- make.seq(LG_E.frame_add106,8,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

LG_E.frame_add80 <- try.seq(LG_E.frame5, 80, draw.try=TRUE)
LG_E.frame_add80 <- try.seq(LG_E.frame5, 80)
LG_E.frame_add80
LG_E.frame5 <- make.seq(LG_E.frame_add80,9,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

LG_E.frame_add261 <- try.seq(LG_E.frame5, 261, draw.try=TRUE)
LG_E.frame_add261 <- try.seq(LG_E.frame5, 261)
LG_E.frame_add261
LG_E.frame5 <- make.seq(LG_E.frame_add261,10,1)
LG_E.frame5
dev.off(dev.list()["RStudioGD"])

maps <- list(LG_E.frame5)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=1)

#####################################################################################
# ADD LG_F


LG_F.frame6 <- make.seq(LGs_2pt, c(9,91,253))
LG_F.compare <-compare(LG_F.frame6)
LG_F.compare
LG_F.frame6 <- make.seq(LG_F.compare,1,1)
LG_F.frame6

## BACKBONE ########## ## # \\/ ##   ###################################

LG_F.frame_add253 <- try.seq(LG_F.frame6, 253, draw.try=TRUE)
LG_F.frame_add253 <- try.seq(LG_F.frame6, 253)
LG_F.frame_add253
LG_F.frame6 <- make.seq(LG_F.frame_add253,20,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

LG_F.frame_add91 <- try.seq(LG_F.frame6, 91, draw.try=TRUE)
LG_F.frame_add91 <- try.seq(LG_F.frame6, 91)
LG_F.frame_add91
LG_F.frame6 <- make.seq(LG_F.frame_add91,81,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

LG_F.frame_add9 <- try.seq(LG_F.frame6, 9, draw.try=TRUE)
LG_F.frame_add9 <- try.seq(LG_F.frame6, 9)
LG_F.frame_add9
LG_F.frame6 <- make.seq(LG_F.frame_add9,82,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

## BACKBONE ########## ## # /\\ ########################################

LG_F.frame_add62 <- try.seq(LG_F.frame6, 62, draw.try=TRUE)
LG_F.frame_add62 <- try.seq(LG_F.frame6, 62)
LG_F.frame_add62
LG_F.frame6 <- make.seq(LG_F.frame_add62,4,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

LG_F.frame_add274 <- try.seq(LG_F.frame6, 274, draw.try=TRUE)
LG_F.frame_add274 <- try.seq(LG_F.frame6, 274)
LG_F.frame_add274
LG_F.frame6 <- make.seq(LG_F.frame_add274,5,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

LG_F.frame_add242 <- try.seq(LG_F.frame6, 242, draw.try=TRUE)
LG_F.frame_add242 <- try.seq(LG_F.frame6, 242)
LG_F.frame_add242
LG_F.frame6 <- make.seq(LG_F.frame_add242,6,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

LG_F.frame_add248 <- try.seq(LG_F.frame6, 248, draw.try=TRUE)
LG_F.frame_add248 <- try.seq(LG_F.frame6, 248)
LG_F.frame_add248
LG_F.frame6 <- make.seq(LG_F.frame_add248,7,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

LG_F.frame_add246 <- try.seq(LG_F.frame6, 246, draw.try=TRUE)
LG_F.frame_add246 <- try.seq(LG_F.frame6, 246)
LG_F.frame_add246
LG_F.frame6 <- make.seq(LG_F.frame_add246,8,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

LG_F.frame_add56 <- try.seq(LG_F.frame6, 56, draw.try=TRUE)
LG_F.frame_add56 <- try.seq(LG_F.frame6, 56)
LG_F.frame_add56
LG_F.frame6 <- make.seq(LG_F.frame_add56,9,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

LG_F.frame_add290 <- try.seq(LG_F.frame6, 290, draw.try=TRUE)
LG_F.frame_add290 <- try.seq(LG_F.frame6, 290)
LG_F.frame_add290
LG_F.frame6 <- make.seq(LG_F.frame_add290,10,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

LG_F.frame_add48 <- try.seq(LG_F.frame6, 48, draw.try=TRUE)
LG_F.frame_add48 <- try.seq(LG_F.frame6, 48)
LG_F.frame_add48
LG_F.frame6 <- make.seq(LG_F.frame_add48,11,1)
LG_F.frame6
dev.off(dev.list()["RStudioGD"])

maps <- list(LG_F.frame6)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=1)



#
# ADD LG_G

LG_G.frame7 <- make.seq(LGs_2pt, c(84,140,76))
LG_G.compare <-compare(LG_G.frame7)
LG_G.compare
LG_G.frame7 <- make.seq(LG_G.compare,1,1)
LG_G.frame7

## BACKBONE ########## ## # \\/ ##   ###################################

LG_G.frame_add84 <- try.seq(LG_G.frame7, 84, draw.try=TRUE)
LG_G.frame_add84 <- try.seq(LG_G.frame7, 84)
LG_G.frame_add84
LG_G.frame7 <- make.seq(LG_G.frame_add84,83,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

LG_G.frame_add76 <- try.seq(LG_G.frame7, 76, draw.try=TRUE)
LG_G.frame_add76 <- try.seq(LG_G.frame7, 76)
LG_G.frame_add76
LG_G.frame7 <- make.seq(LG_G.frame_add76,85,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

LG_G.frame_add140 <- try.seq(LG_G.frame7, 140, draw.try=TRUE)
LG_G.frame_add140 <- try.seq(LG_G.frame7, 140)
LG_G.frame_add140
LG_G.frame7 <- make.seq(LG_G.frame_add140,4,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

## BACKBONE ########## ## # /\\ ########################################



LG_G.frame_add20 <- try.seq(LG_G.frame7, 20, draw.try=TRUE)
LG_G.frame_add20 <- try.seq(LG_G.frame7, 20)
LG_G.frame_add20
LG_G.frame7 <- make.seq(LG_G.frame_add20,4,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

LG_G.frame_add193 <- try.seq(LG_G.frame7, 193, draw.try=TRUE)
LG_G.frame_add193 <- try.seq(LG_G.frame7, 193)
LG_G.frame_add193
LG_G.frame7 <- make.seq(LG_G.frame_add193,5,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

LG_G.frame_add117 <- try.seq(LG_G.frame7, 117, draw.try=TRUE)
LG_G.frame_add117 <- try.seq(LG_G.frame7, 117)
LG_G.frame_add117
LG_G.frame7 <- make.seq(LG_G.frame_add117,6,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

LG_G.frame_add18 <- try.seq(LG_G.frame7, 18, draw.try=TRUE)
LG_G.frame_add18 <- try.seq(LG_G.frame7, 18)
LG_G.frame_add18
LG_G.frame7 <- make.seq(LG_G.frame_add18,7,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

LG_G.frame_add217 <- try.seq(LG_G.frame7, 217, draw.try=TRUE)
LG_G.frame_add217 <- try.seq(LG_G.frame7, 217)
LG_G.frame_add217
LG_G.frame7 <- make.seq(LG_G.frame_add217,8,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

LG_G.frame_add15 <- try.seq(LG_G.frame7, 15, draw.try=TRUE)
LG_G.frame_add15 <- try.seq(LG_G.frame7, 15)
LG_G.frame_add15
LG_G.frame7 <- make.seq(LG_G.frame_add15,9,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

LG_G.frame_add304 <- try.seq(LG_G.frame7, 304, draw.try=TRUE)
LG_G.frame_add304 <- try.seq(LG_G.frame7, 304)
LG_G.frame_add304
LG_G.frame7 <- make.seq(LG_G.frame_add304,10,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])

LG_G.frame_add59 <- try.seq(LG_G.frame7, 59, draw.try=TRUE)
LG_G.frame_add59 <- try.seq(LG_G.frame7, 59)
LG_G.frame_add59
LG_G.frame7 <- make.seq(LG_G.frame_add59,11,1)
LG_G.frame7
dev.off(dev.list()["RStudioGD"])



#
############# Add LG_J ## # ###########

LG_J.frame8 <- make.seq(LGs_2pt, c(13,196,223))
LG_J.compare <-compare(LG_J.frame8)
LG_J.compare
LG_J.frame8 <- make.seq(LG_J.compare,1,1)
LG_J.frame8

## BACKBONE ########## ## # \\/ ##   ###################################

LG_J.frame_add13 <- try.seq(Master, 13, draw.try=TRUE)
LG_J.frame_add13 <- try.seq(Master, 13)
LG_J.frame_add13
LG_J.frame8 <- make.seq(LG_J.frame_add13,95,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

LG_J.frame_add196 <- try.seq(LG_J.frame8, 196, draw.try=TRUE)
LG_J.frame_add196 <- try.seq(LG_J.frame8, 196)
LG_J.frame_add196
LG_J.frame8 <- make.seq(LG_J.frame_add196,96,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

LG_J.frame_add223 <- try.seq(LG_J.frame8, 223, draw.try=TRUE)
LG_J.frame_add223 <- try.seq(LG_J.frame8, 223)
LG_J.frame_add223
LG_J.frame8 <- make.seq(LG_J.frame_add223,97,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

## BACKBONE ########## ## # /\\ ########################################



LG_J.frame_add200 <- try.seq(LG_J.frame8, 200, draw.try=TRUE)
LG_J.frame_add200 <- try.seq(LG_J.frame8, 200)
LG_J.frame_add200
LG_J.frame8 <- make.seq(LG_J.frame_add200,4,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

LG_J.frame_add174 <- try.seq(LG_J.frame8, 174, draw.try=TRUE)
LG_J.frame_add174 <- try.seq(LG_J.frame8, 174)
LG_J.frame_add174
LG_J.frame8 <- make.seq(LG_J.frame_add174,5,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

LG_J.frame_add185 <- try.seq(LG_J.frame8, 185, draw.try=TRUE)
LG_J.frame_add185 <- try.seq(LG_J.frame8, 185)
LG_J.frame_add185
LG_J.frame8 <- make.seq(LG_J.frame_add185,6,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

LG_J.frame_add220 <- try.seq(LG_J.frame8, 220, draw.try=TRUE)
LG_J.frame_add220 <- try.seq(LG_J.frame8, 220)
LG_J.frame_add220
LG_J.frame8 <- make.seq(LG_J.frame_add220,7,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

LG_J.frame_add29 <- try.seq(LG_J.frame8, 29, draw.try=TRUE)
LG_J.frame_add29 <- try.seq(LG_J.frame8, 29)
LG_J.frame_add29
LG_J.frame8 <- make.seq(LG_J.frame_add29,8,1)
LG_J.frame8
dev.off(dev.list()["RStudioGD"])

############ LG_J /\\ #######################

# Add LG_H \\/ ###################

LG_H.frame8 <- make.seq(LGs_2pt, c(66,95,154))
LG_H.compare <-compare(LG_H.frame8)
LG_H.compare
LG_H.frame8 <- make.seq(LG_H.compare,1,1)
LG_H.frame8

## BACKBONE ########## ## # \\/ ##   ###################################

LG_H.frame_add154 <- try.seq(LG_H.frame8, 154, draw.try=TRUE)
LG_H.frame_add154 <- try.seq(LG_H.frame8, 154)
LG_H.frame_add154
LG_H.frame8 <- make.seq(LG_H.frame_add154,106,1)
LG_H.frame8
dev.off(dev.list()["RStudioGD"])

LG_H.frame_add95 <- try.seq(LG_H.frame8, 95, draw.try=TRUE)
LG_H.frame_add95 <- try.seq(LG_H.frame8, 95)
LG_H.frame_add95
LG_H.frame8 <- make.seq(LG_H.frame_add95,107,1)
LG_H.frame8
dev.off(dev.list()["RStudioGD"])

LG_H.frame_add66 <- try.seq(LG_H.frame8, 66, draw.try=TRUE)
LG_H.frame_add66 <- try.seq(LG_H.frame8, 66)
LG_H.frame_add66
LG_H.frame8 <- make.seq(LG_H.frame_add66,108,1)
LG_H.frame8
dev.off(dev.list()["RStudioGD"])

## BACKBONE ########## ## # /\\ ########################################

LG_H.frame_add116 <- try.seq(LG_H.frame8, 116, draw.try=TRUE)
LG_H.frame_add116 <- try.seq(LG_H.frame8, 116)
LG_H.frame_add116
LG_H.frame8 <- make.seq(LG_H.frame_add116,4,1)
LG_H.frame8
dev.off(dev.list()["RStudioGD"])

LG_H.frame_add228 <- try.seq(LG_H.frame8, 228, draw.try=TRUE)
LG_H.frame_add228 <- try.seq(LG_H.frame8, 228)
LG_H.frame_add228
LG_H.frame8 <- make.seq(LG_H.frame_add228,5,1)
LG_H.frame8
dev.off(dev.list()["RStudioGD"])

LG_H.frame_add133 <- try.seq(LG_H.frame8, 133, draw.try=TRUE)
LG_H.frame_add133 <- try.seq(LG_H.frame8, 133)
LG_H.frame_add133
LG_H.frame8 <- make.seq(LG_H.frame_add133,6,1)
LG_H.frame8
dev.off(dev.list()["RStudioGD"])


#

# ADD LG_I ######################

LG_I.frame8 <- make.seq(LGs_2pt, c(44,310,244))
LG_I.compare <-compare(LG_I.frame8)
LG_I.compare
LG_I.frame8 <- make.seq(LG_I.compare,1,1)
LG_I.frame8

## BACKBONE ########## ## # \\/ ##   ###################################

LG_I.frame_add244 <- try.seq(LG_I.frame8, 244, draw.try=TRUE)
LG_I.frame_add244 <- try.seq(LG_I.frame8, 244)
LG_I.frame_add244
LG_I.frame8 <- make.seq(LG_I.frame_add244,113,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])

LG_I.frame_add310 <- try.seq(LG_I.frame8, 310, draw.try=TRUE)
LG_I.frame_add310 <- try.seq(LG_I.frame8, 310)
LG_I.frame_add310
LG_I.frame8 <- make.seq(LG_I.frame_add310,114,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])

LG_I.frame_add44 <- try.seq(LG_I.frame8, 44, draw.try=TRUE)
LG_I.frame_add44 <- try.seq(LG_I.frame8, 44)
LG_I.frame_add44
LG_I.frame8 <- make.seq(LG_I.frame_add44,115,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])

## BACKBONE ########## ## # /\\ ########################################

LG_I.frame_add125 <- try.seq(LG_I.frame8, 125, draw.try=TRUE)
LG_I.frame_add125 <- try.seq(LG_I.frame8, 125)
LG_I.frame_add125
LG_I.frame8 <- make.seq(LG_I.frame_add125,4,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])

LG_I.frame_add119 <- try.seq(LG_I.frame8, 119, draw.try=TRUE)
LG_I.frame_add119 <- try.seq(LG_I.frame8, 119)
LG_I.frame_add119
LG_I.frame8 <- make.seq(LG_I.frame_add119,5,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])

LG_I.frame_add236 <- try.seq(LG_I.frame8, 236, draw.try=TRUE)
LG_I.frame_add236 <- try.seq(LG_I.frame8, 236)
LG_I.frame_add236
LG_I.frame8 <- make.seq(LG_I.frame_add236,6,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])

LG_I.frame_add105 <- try.seq(LG_I.frame8, 105, draw.try=TRUE)
LG_I.frame_add105 <- try.seq(LG_I.frame8, 105)
LG_I.frame_add105
LG_I.frame8 <- make.seq(LG_I.frame_add105,7,1)
LG_I.frame8
dev.off(dev.list()["RStudioGD"])


maps <- list(LG_A.frame1,LG_B.frame2,LG_C.frame3,LG_D.frame4,LG_E.frame5,LG_F.frame6,LG_G.frame7,LG_H.frame8,LG_I.frame8,LG_J.frame8)
draw.map(maps, names=TRUE, grid=TRUE, cex.mrk=.75)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#