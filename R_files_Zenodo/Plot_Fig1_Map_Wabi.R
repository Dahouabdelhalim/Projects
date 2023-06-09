################################################################################
###                                                                          ###
###----- NAME: --------------------------------------------------------------###
###                                                                          ###
### Plot_Fig1_Map_Wabi.R                                                     ###
###                                                                          ###
###                                                                          ###
###----- AIM: ---------------------------------------------------------------###
###                                                                          ###
### Plot a world map with points for each city. Around the map add for       ###
### each city a small WABI indent with the traffic light categories. Add a   ###
### panel explaining the WABIs color scheme either to the right (Horizontal) ###
### or below (Vertical) the map.                                             ###
###                                                                          ###
###----- INPUT --------------------------------------------------------------###
###                                                                          ###
### 1.  "Map"      - excel sheet.                                            ###
###     A sheet from the excel file 'VelisEtAl2022_WABIs_Input.xlsx'.        ###
###     Contains the data required to plot the map (names, locations, WABI   ###
###     colors, and location of inset).                                      ###
###     . Contains the following columns:                                    ###
###   "Num"      - Running number.                                           ###
###   "casename" - A unique case name.                                       ###                    
###   "Country"  - Name of country.                                          ###                     
###   "ISO3CountryCode" - 3 letter ISO3 code for countries.                  ###               
###   "Y_WGS84" - Spatial location in WGS84.                                 ###                     
###   "X_WGS84" - Spatial location in WGS84.                                 ###                      
###   "IncomeCategory" - Income category of the country. Either "L", "L-M",  ###
###                      "U-M" or "H".                                       ###                
###   "City"           - Name of the city.                                   ###               
###   "WastePerCapita" - WABI I-0 - Waste generation rate [Kg*y^-1*p^-1].    ###               
###  "CollectionCoverage"     - WABI I-1.1 - Waste collection coverage.      ###           
###  "WasteQualityCollection" - WABI I-1C - Quality of waste collection.     ###
###                             service.                                     ###        
###  "ControlledDisposal" - WABI I-2 - Controlled treatment and disposal.    ###     
###  "QualityEnvironmentalProtection" - WABI I-2E - Environmental protection ###
###                                   in controlled treatment and disposal.  ###
###   "Generetation_col"       - Color for WABI I-0.                         ###              
###   "CollectionCoverage_col" _ Color for WABI I-1.1.                       ###       
###   "ColectionQuality_col"   _ Color for WABI I-1C.                        ###       
###   "DisposalCoverage_col"   _ Color for WABI I-2.                         ###      
###   "DisposalQuality_col"    _ Color for WABI I-2E.                        ###       
###   "Group"       - Location of inset:                                     ###
###                    1: top;  2: right;  3: Bottom;  4: Left.              ###                       
###   "RankInGroup" - Order of inset around the map (in each Group).         ###
###                                                                          ###
### 2.  "Fig1_Wabi_Horizontal.jpg"  ;  "Fig1_Wabi_Vertical.jpg"              ###
###    A JPEG image (300dpi) with the WABI side panel, that will appear on   ###
###    the right of (Horizontal), or below (Vertical) the map.               ###
###                                                                          ###
###----- PROCESS: -----------------------------------------------------------###
###                                                                          ###
### 1.  Packages.                                                            ###                                    
### 2.  Directories.                                                         ###                                       
### 3.  Controls.                                                            ###                                    
###   3.1.  General controls.                                                ###
###   3.2.  colors.                                                          ###
###   3.3.  Income category - data and colors.                               ###
###   3.4.  Controls - to console.                                           ###
###   3.5.  Controls - to file.                                              ###
### 4.  Read data.                                                           ###                                
### 5.  Prepare data for world map.                                          ###                                        
### 6.  Calculate inset and lines locations.                                 ###                                        
###   6.1.  Group 1.                                                         ###
###   6.2.  Group 2.                                                         ###
###   6.3.  Group 3.                                                         ###
###   6.4.  Group 4.                                                         ###
###   6.5.  Bind.                                                            ###
### 7.  Add colors and factor - Income category.                             ###                                       
### 8.  Plot the main map.                                                   ###                                     
###   8.1.  Start plot with world map.                                       ###
###   8.2.  Define coordinates limits.                                       ###
###   8.3.  Rect around map area (excluding inset area).                     ###
###   8.4.  Lines from points to inset + points.                             ###
###   8.5.  Set the theme for the legend.                                    ###
###   8.6.  Get the legend and insert it as inset.                           ###
###   8.7.  Add tag as text.                                                 ###
###   8.8.  Adjust scales, titles and theme.                                 ###
### 9.  Add the insets in a loop.                                            ###                                        
###   9.1.  Create the inset plot.                                           ###
###   9.2.  Change to grob and overlay on the correct location on the map.   ###
### 10.  Prepare the WABI panel and arrange panels in grob.                  ### 
### 11.  Write to console or file.                                           ###
###                                                                          ###
###----- OUTPUT: ------------------------------------------------------------###
###                                                                          ###
### 1.  A 2 panel figure with a map and a panel explaining the WABIs. The    ###
###     map includes inset for each city with the WABI colors.               ###
###                                                                          ###
###----- NOTES --------------------------------------------------------------###
###                                                                          ###
### 1.                                                                       ###
###                                                                          ###
###----- VERSIONS -----------------------------------------------------------###
###                                                                          ###
### 1.  16 Feb. 2022 - Original code from varioues earlier versions.         ###
###                                                                          ###
###----- DETAILS ------------------------------------------------------------###
###                                                                          ###
### 1. Author:      Gavish Yoni <gavishyoni@gmail.com>                       ###
###                                                                          ###
################################################################################


###--------------------------------------------------------------------------###
### 1.  Packages.                                         
###--------------------------------------------------------------------------###


library(readxl)
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(extrafont)           # for pdf output
library(ggmap)
library(ggplotify)
library(ggpubr)

library(sf)
library(sp)
library(raster)
library(jpeg)


###--------------------------------------------------------------------------###
### 2.  Directories.                                         
###--------------------------------------------------------------------------###


Dir1 <- "./Analysis/Input"
Dir2 <- "./Analysis/Figures"


###--------------------------------------------------------------------------###
### 3.  Controls.                                         
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
### 3.1.  General controls.

ToFile <- T
FormOut <- c(".png", ".pdf")     # formats of output files
#FormOut  <- ".png"                # formats of output files
Name     <- "VelisEtAL_Fig1_Map_WABI"
TagText <- "A:"

DoHorizontal <- F

###--------------------------------------------------------------------------###
### 3.2.  colors.

## inset
ColEdgeInset  <- "gray55" # border of inset shapes
FillInset     <- "white"  # fill of inset panels
BorderInset   <- "gray15" # border of inset panel

## map
ColMapRect     <- "black"  # color of box around map

###--------------------------------------------------------------------------###
### 3.3.  Income category - data and colors.

IncCat <-data.frame(IncomeCategory = c("L",   
                                       "L-M", 
                                       "U-M", 
                                       "H"),
                    Color = c("#F2E8F3",
                              "#BCA4B2",
                              "#AE8588",
                              "#71779B"),
                    Border = c("black", 
                               "black",
                               "black",
                               "black"),
                    Label = c("Low",   
                              "Lower-Middle", 
                              "Upper-Middle", 
                              "High"))
IncCat$IncomeCategory <- factor(IncCat$IncomeCategory, 
                                levels = IncCat$IncomeCategory)

# named vector for legend
IntCatCols        <- IncCat$Color
names(IntCatCols) <- IncCat$Label

###--------------------------------------------------------------------------###
### 3.4.  Controls - to console.

if(!ToFile){
  ###------------------------------------------------------------------------###
  ### General  - figure dimension/resolution
  Dim <- c(27, 15.5) # width/height cm
  PanelGrid_Size <- 0.2
  
  ###------------------------------------------------------------------------###
  ### map
  MapFill      <- "gray75"
  MapColor     <- "gray50"
  MapPoly_Size <- 0.3
  
  ###------------------------------------------------------------------------###
  ### inset plots layout Ranges (1 - top ; 2 - right ; 3 - bottom ; 4- left)
  
  RangeGroup1 <- c(-150, 150)
  RangeGroup2 <- c(-50, 70)
  RangeGroup3 <- c(-150, 150)
  RangeGroup4 <- c(-25, 40)
  
  ###------------------------------------------------------------------------###
  ### Points - size points on location and near inset on the margins 
  SizePoints       <- 3
  SizePointsMargin <- 5
  
  ###------------------------------------------------------------------------###
  ### inset - Title text size, size of inset on map, gap from margin
  SizeNameCase <- 2
  SizeInsetX   <- 10 # number of unit along the x axis from the inset center
  SizeInsetY   <- 10 # number of unit along the y axis from the inset center
  GapExpand    <- 1
  
  SizeLinesInset <- 0.5
  SizeBorderLine <- 0.5
  
  ###------------------------------------------------------------------------###
  ### inset lines to points - 
  ColSegments     <- "blue" 
  SizeLines    <- 0.5
  Seg_LineType <- "dotted"
  
  ###------------------------------------------------------------------------###
  ### legend - location, fill, color, title size, text size
  
  Legend_box      <- c(-160, -140, -67.5, -33)
  Leg_Fill        <- "white"
  Leg_Border      <- "black"
  Leg_Border_Size <- 0.4
  SizeLegTitle    <- 5
  SizeLegText     <- 4
  LegKeySize      <- 2 # legend key size in points
  
  ###------------------------------------------------------------------------###
  ### Axis
  SizeAxesText <- 5 
  
  ###------------------------------------------------------------------------###
  ### Tag
  TagXY <- c(-190, 100)
  SizeTag <- 8
}

###--------------------------------------------------------------------------###
### 3.5.  Controls - to file.

if(ToFile){
  ###------------------------------------------------------------------------###
  ### General  - figure dimension/resolution
  
  if(DoHorizontal){
    Dim <- c(35, 15.5) # width/height cm
  }
  
  if(!DoHorizontal){
    Dim <- c(27, 26.0) # width/height cm
  }
  
  PanelGrid_Size <- 0.2
  
  ###------------------------------------------------------------------------###
  ### map
  MapFill      <- "gray75"
  MapColor     <- "gray50"
  MapPoly_Size <- 0.3
  
  ###------------------------------------------------------------------------###
  ### inset plots layout Ranges (1 - top ; 2 - right ; 3 - bottom ; 4- left)
  
  RangeGroup1 <- c(-170, 168)
  RangeGroup2 <- c(-88, 82)
  RangeGroup3 <- c(-170, 168)
  RangeGroup4 <- c(-10, 70)
  
  ###------------------------------------------------------------------------###
  ### Points - size points on location and near inset on the margins 
  
  SizePoints       <- 2
  SizePointsMargin <- 3.5
  
  ###------------------------------------------------------------------------###
  ### inset - Title text size, size of inset on map, gap from margin
  SizeNameCase <- 1.7
  SizeInsetX   <- 15 # number of unit along the x axis from the inset center
  SizeInsetY   <- 15 # number of unit along the y axis from the inset center
  GapExpand    <- 0.0
  
  SizeLinesInset <- 0.3
  SizeBorderLine <- 0.3
  
  ###------------------------------------------------------------------------###
  ### inset lines to points - 
  
  ColSegments     <- "steelblue4"
  SizeLines    <- 0.25
  Seg_LineType <- "dashed"
  
  ###------------------------------------------------------------------------###
  ### legend - location, fill, color, title size, text size
  
  Legend_box      <- c(-170, -130, -60, -33)
  Leg_Fill        <- "white"
  Leg_Border      <- "black"
  Leg_Border_Size <- 0.3
  SizeLegTitle    <- 6
  SizeLegText     <- 5
  LegKeySize      <- 2 # legend key size in points
  
  ###------------------------------------------------------------------------###
  ### Axis
  SizeAxesText <- 5
  
  ###------------------------------------------------------------------------###
  ### Tag
  TagXY <- c(-200, 105)
  SizeTag <- 4.5
}


###--------------------------------------------------------------------------###
### 4.  Read data.                                         
###--------------------------------------------------------------------------###


# cities data
setwd(Dir1)

D0 <- read_excel(path = paste0(Dir1, "/VelisEtAl2022_WABIs_Input.xlsx"),
                   sheet = "Map")


# change to data.frame and factors for back compatibility 

D0 <- as.data.frame(D0, stringsAsFactors = T)

for(i in 1:ncol(D0)){
  if(class(D0[, i]) == "character"){
    D0[, i] <- factor(D0[, i])
  }
}


# WABI panel 
if(DoHorizontal){
  Img1 <- readJPEG("Fig1_Wabi_Horizontal.jpg")
}

if(!DoHorizontal){
  Img1 <- readJPEG("Fig1_Wabi_Vertical.jpg")
}


###--------------------------------------------------------------------------###
### 5.  Prepare data for world map.                                         
###--------------------------------------------------------------------------###


# world map
World <- map_data('world', wrap=c(-180,180))

# bounding box and final coord range
WorldBox <- c(-180, 180,-89.9, 83.6)
CordLim  <- c(WorldBox[1] - 2 * GapExpand - 2 * SizeInsetX, 
              WorldBox[2] + 2 * GapExpand + 2 * SizeInsetX,
              WorldBox[3] - 2 * GapExpand - 2 * SizeInsetY, 
              WorldBox[4] + 2 * GapExpand + 2 * SizeInsetY)


###--------------------------------------------------------------------------###
### 6.  Calculate inset and lines locations.                                         
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
### 6.1.  Group 1.

D1 <- D0[D0$Group == 1, ]
D1 <- D1[order(D1$RankInGroup),]
D1$InsetX   <- seq(RangeGroup1[1], RangeGroup1[2], length.out = nrow(D1)) 
D1$InsetY   <- WorldBox[4] + GapExpand + SizeInsetY
D1$SegmentX <- D1$InsetX
D1$SegmentY <- WorldBox[4]

###--------------------------------------------------------------------------###
### 6.2.  Group 2.

D2 <- D0[D0$Group == 2, ]
D2 <- D2[order(D2$RankInGroup),]
D2$InsetX   <- WorldBox[2] + GapExpand + SizeInsetX
D2$InsetY   <- seq(RangeGroup2[1], RangeGroup2[2], length.out = nrow(D2)) 
D2$SegmentY <- D2$InsetY
D2$SegmentX <- WorldBox[2]

###--------------------------------------------------------------------------###
### 6.3.  Group 3.

D3 <- D0[D0$Group == 3, ]
D3 <- D3[order(D3$RankInGroup),]
D3$InsetX   <- seq(RangeGroup3[1], RangeGroup3[2], length.out = nrow(D3)) 
D3$InsetY   <- WorldBox[3] - GapExpand - SizeInsetY
D3$SegmentX <- D3$InsetX
D3$SegmentY <- WorldBox[3]

###--------------------------------------------------------------------------###
### 6.4.  Group 4.

D4 <- D0[D0$Group == 4, ]
D4 <- D4[order(D4$RankInGroup),]
D4$InsetX   <- WorldBox[1] - GapExpand - SizeInsetX
D4$InsetY   <- seq(RangeGroup4[1], RangeGroup4[2], length.out = nrow(D4)) 
D4$SegmentY <- D4$InsetY
D4$SegmentX <- WorldBox[1]

###--------------------------------------------------------------------------###
### 6.5.  Bind.

Data <- rbind(D1, D2, D3, D4)
rm(D0, D1, D2, D3, D4)


###--------------------------------------------------------------------------###
### 7.  Add colors and factor - Income category.                                         
###--------------------------------------------------------------------------###


Data$IncCat_Fill <- IncCat[match(Data$IncomeCategory, IncCat$IncomeCategory), "Color"]

Data$IncomeCategory <- factor(Data$IncomeCategory, 
                              levels = IncCat$IncomeCategory, 
                              label  = IncCat$Label)


###--------------------------------------------------------------------------###
### 8.  Plot the main map.                                         
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
### 8.1.  Start plot with world map.

P1 <- ggplot(World, aes(long, lat))
P1 <- P1 + geom_map(map=World, aes(map_id=region), 
                    fill  = MapFill, 
                    color = MapColor, 
                    size  = MapPoly_Size)

###--------------------------------------------------------------------------###
### 8.2.  Define coordinates limits.
P1 <- P1 + coord_quickmap(xlim = c(CordLim[1], CordLim[2]),
                          ylim = c(CordLim[3], CordLim[4]), 
                          expand = F)

###--------------------------------------------------------------------------###
### 8.3.  Rect around map area (excluding inset area).

P1 <- P1 + geom_rect(aes(xmin = WorldBox[1], 
                         xmax = WorldBox[2], 
                         ymin = WorldBox[3], 
                         ymax = WorldBox[4]), 
                     fill  = NA, 
                     color = ColMapRect)

###--------------------------------------------------------------------------###
### 8.4.  Lines from points to inset + points.

P1 <- P1 + geom_segment(aes(x = SegmentX, 
                            y = SegmentY, 
                            xend = X_WGS84,
                            yend = Y_WGS84), 
                        color    = ColSegments,
                        linetype = Seg_LineType,
                        size     = SizeLines,
                        data     = Data)

P1 <- P1 + geom_point(aes(x = X_WGS84, 
                          y = Y_WGS84),
                      fill  = Data$IncCat_Fill,
                      color = ColSegments,
                      shape = 21,
                      size  = SizePoints,
                      data  = Data)

P1 <- P1 + geom_point(aes(x = SegmentX,
                          fill = IncomeCategory,
                          y = SegmentY),
                      color = ColSegments,
                      shape = 21,
                      size  = SizePointsMargin,
                      data  = Data)

###--------------------------------------------------------------------------###
### 8.5.  Set the theme for the legend. 

P1 <- P1 + scale_fill_manual(values = IntCatCols)
P1 <- P1 + guides(fill=guide_legend(title="Income category"))

P1 <- P1 + theme(
  legend.title = element_text(size = SizeLegTitle),
  legend.text  = element_text(size = SizeLegText),
  legend.background = element_rect(size  = Leg_Border_Size , 
                                   color = Leg_Border, 
                                   fill  = Leg_Fill),
  legend.key = element_rect(fill = "white"),
  legend.key.size = unit(LegKeySize, "points"))

###--------------------------------------------------------------------------###
### 8.6. Get the legend and insert it as inset.

P2 <- get_legend(P1)
P1 <- P1 + inset(P2, 
                 xmin = Legend_box[1], 
                 xmax = Legend_box[2], 
                 ymin = Legend_box[3], 
                 ymax = Legend_box[4])
rm(P2)

###--------------------------------------------------------------------------###
### 8.7. Add tag as text.

P1 <- P1 + geom_text(aes(x     = TagXY[1],
                         y     = TagXY[2],
                         label = TagText),
                     size  = SizeTag,
                     color = "black")

###--------------------------------------------------------------------------###
### 8.8.  Adjust scales, titles and theme.

P1 <- P1 + scale_x_continuous(breaks       = seq(-180, 180, 60),
                              minor_breaks = seq(-180, 180, 30))
P1 <- P1 + scale_y_continuous(breaks       = seq(-100, 100, 50),
                              minor_breaks = seq(-100, 100, 25))
P1 <- P1 + theme_bw()
P1 <- P1 + theme(
  legend.position = "None",
  axis.title = element_blank(), 
  panel.grid = element_line(size = PanelGrid_Size),
  axis.text  = element_text(size = SizeAxesText)
)


###--------------------------------------------------------------------------###
### 9.  Add the insets in a loop.                                         
###--------------------------------------------------------------------------###


for(i in 1:nrow(Data)){

    ###----------------------------------------------------------------------###
    ### 9.1.  Create the inset plot.
    P2 <- ggplot()
    
    P2 <- P2 + geom_rect(aes(xmin = 0, xmax = 1.45, ymin =0.1, ymax = 0.9), 
                         fill  = Data[i, "DisposalCoverage_col"], 
                         size  = SizeLinesInset,
                         color = ColEdgeInset)
    P2 <- P2 + geom_rect(aes(xmin = 1.55, xmax = 3.0, ymin =0.1, ymax = 0.9), 
                         fill  = Data[i, "DisposalQuality_col"],
                         size  = SizeLinesInset,
                         color = ColEdgeInset)
    P2 <- P2 + geom_rect(aes(xmin = 0, xmax = 1.45, ymin =1.1, ymax = 1.9), 
                         fill  = Data[i, "CollectionCoverage_col"],
                         size  = SizeLinesInset,
                         color = ColEdgeInset)
    P2 <- P2 + geom_rect(aes(xmin = 1.55, xmax = 3.0, ymin =1.1, ymax = 1.9), 
                         fill  = Data[i, "ColectionQuality_col"], 
                         size  = SizeLinesInset,
                         color = ColEdgeInset)
    P2 <- P2 + geom_rect(aes(xmin = 0, xmax = 3, ymin =2.1, ymax = 2.8),
                         fill  = Data[i, "Generetation_col"],
                         size  = SizeLinesInset,
                         color = ColEdgeInset)
    
    P2 <- P2 + geom_text(aes(x = 1.5, y = 3.2, label = Data[i, "casename"]), 
                         size = SizeNameCase)
   
    P2 <- P2 + ylim(0, 3.6) + xlim(-0.2, 3.2)
    
    P2 <- P2 + theme_minimal()
    P2 <- P2 + theme(
      axis.title       = element_blank(),
      axis.text        = element_blank(),
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = FillInset, color = BorderInset, size = SizeBorderLine)
    )
    
    ###----------------------------------------------------------------------###
    ### 9.2.  Change to grob and overlay on the correct location on the map.
    
    P2 <- as.grob(P2)
    
    CenterPointX <-Data[i, "InsetX"]
    CenterPointY <-Data[i, "InsetY"]
    
    P1 <- P1 + inset(P2, 
                     xmin = CenterPointX - SizeInsetX, 
                     xmax = CenterPointX + SizeInsetX, 
                     ymin = CenterPointY - SizeInsetY, 
                     ymax = CenterPointY + SizeInsetY)
    
    rm(P2, CenterPointX, CenterPointY)
    
}

rm(i)


###--------------------------------------------------------------------------###
### 10.  Prepare the WABI panel and arrange panels in grob. 
###--------------------------------------------------------------------------###


# create a grob from the raster
Img2 <- rasterGrob(Img1, interpolate=TRUE)


# arrange grob
if(DoHorizontal){
  G1 <- arrangeGrob(P1,Img2,
                    ncol    = 2,
                    nrow    = 1, 
                    widths  = c(15,5),
                    heights = 10) 
  Name = paste0(Name, "_Horizontal")
}

if(!DoHorizontal){
  G1 <- arrangeGrob(P1,Img2,
                    ncol    = 1,
                    nrow    = 2, 
                    widths  = 27,
                    heights = c(15.5, 11)) 
  
  Name = paste0(Name, "_Vertical")
}


###--------------------------------------------------------------------------###
### 11.  Write to console or file. 
###--------------------------------------------------------------------------###


if(!ToFile){
  plot.new()
  grid.draw(G1)
}

if(ToFile){
  setwd(Dir2)
  for(i in FormOut){
    if(i != ".pdf"){
      ggsave(file   = paste0(Name, i), 
             plot   = G1,
             width  = Dim[1], 
             height = Dim[2], 
             units  = c("cm"),
             dpi    = 300)
    }
    # use cairo for special characters. 
    if(i == ".pdf"){
      ggsave(file   = paste0(Name, i), 
             plot   = G1,
             width  = Dim[1], 
             height = Dim[2], 
             units  = c("cm"),
             dpi    = 300, 
             device = cairo_pdf)
    }
    
  }
  ###------------------------------------------------------------------------###
  ### clean
  rm(i)
  
}




################################################################################
###--------------------------------------------------------------------------###
###----------------------                              ----------------------### 
###-----------------                                        -----------------###
###-----------                                                    -----------###
###-----                            END CODE                            -----###
###-----------                                                    -----------###
###-----------------                                        -----------------###
###----------------------                              ----------------------### 
###--------------------------------------------------------------------------###
################################################################################






