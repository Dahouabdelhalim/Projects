library(data.table)
library(readxl)
library(arules)
library(ggrepel)
library(colorspace)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls(all=TRUE))

full.data <- read_xlsx("Data_supplement_Kastner_et_al_2021.xlsx","data_Fig2-Fig4", skip = 2, col_names = F)
full.data <- as.data.table(full.data)

setnames(full.data, c("Producer.country",	"Consumer.country", "Crop.type",
                      "kcal","area","species.loss","CS.forgone","def.area","def.em",
                      "area.no.trade","species.loss.no.trade","CS.forgeone.no.trade"))

plot.data <- full.data[!Producer.country == Consumer.country, 
                       .(kcal = sum(kcal), area  = sum(area),
                         species.loss = sum(species.loss), def.area = sum(def.area)), 
                       by = "Crop.type"]

temp.data <- full.data[, .(kcal.total = sum(kcal), area.total  = sum(area),
                         species.loss.total = sum(species.loss), def.area.total = sum(def.area)), 
                       by = "Crop.type"]

plot.data <- merge(plot.data,temp.data, by = "Crop.type")
rm(temp.data)
plot.data[, kcal.exp.share :=  kcal / kcal.total]
plot.data <- plot.data[!kcal == 0]
plot.data[,area:=area*100/1000/1000]


plot.data[, exp.share.classes := discretize(plot.data[, kcal.exp.share],
                                            method = "frequency", breaks = 5,
                                            labels = c("very low","low","medium","high","very high"))]

plot.data[Crop.type == "soyb", Crop.type := "Soybeans"]
plot.data[Crop.type == "whea", Crop.type := "Wheat"]
plot.data[Crop.type == "maiz", Crop.type := "Maize"]
plot.data[Crop.type == "oilp", Crop.type := "Oil palm"]
plot.data[Crop.type == "suga", Crop.type := "Sugar crops"]
plot.data[Crop.type == "barl", Crop.type := "Barley"]
plot.data[Crop.type == "rape", Crop.type := "Rapeseed"]
plot.data[Crop.type == "rice", Crop.type := "Rice"]
plot.data[Crop.type == "cnut", Crop.type := "Coconut"]
plot.data[Crop.type == "coff", Crop.type := "Coffee"]
plot.data[Crop.type == "coco", Crop.type := "Cocoa"]
plot.data[Crop.type == "rest", Crop.type := "Other crops"]


png(paste0("Fig4A.png"), width=8, height=6, units="in", res=600)

ggplot(data = plot.data, aes(x = kcal, y = species.loss, size = area, color = exp.share.classes)) +
  geom_point(alpha=0.5) +
  scale_size(range = c(1, 20), name="Area [Mha]") +
  geom_text_repel(aes(label=as.character(Crop.type)),
                  segment.color = 'black',
                  box.padding = 1,
                  size = 6,
                  color = "black")+
  scale_color_manual(values = rev(sequential_hcl(5, "Viridis")))+
  geom_point(alpha=1, shape = 21,color = "black")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        text = element_text(size=20)) +
  labs(color = "Export share")

dev.off()

png(paste0("Fig4B.png"), width=8, height=6, units="in", res=600)

ggplot(data = plot.data, aes(x = kcal, y = def.area, size = area, color = exp.share.classes)) +
  geom_point(alpha=0.5) +
  scale_size(range = c(1, 20), name="Area [Mha]") +
  geom_text_repel(aes(label=as.character(Crop.type)),
                  segment.color = 'black',
                  box.padding = 1,
                  size = 6,
                  color = "black")+
  scale_color_manual(values = rev(sequential_hcl(5, "Viridis")))+
  geom_point(alpha=1, shape = 21,color = "black")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        text = element_text(size=20)) +
  labs(color = "Export share")

dev.off()