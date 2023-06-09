library(ggplot2)
library(forcats)
library(RColorBrewer)
library(viridis)
library(scales)

#setwd
cp.data<-read.csv("Canonical Pathways_Anterior-Posterior-WingsvLegsNA.csv", header=T)
cp.data$X.log..p.value.
pa.cpdata<-read.csv("Posterior vs Anterior Canonical Pathways Significantly Enriched.csv", header=TRUE)
?reorder


plot<-ggplot(data=cp.data, aes(x=fct_rev(Tissue), y=X.log..p.value., fill = z.score)) +geom_bar(position = position_dodge(),stat="identity", alpha=0.7) + coord_flip() + scale_fill_viridis_c(option="viridis")
plot + facet_grid(Canonical.Pathway~.) + ylim(0, 40) 

pl<-ggplot(data=cp.data, aes(x=Tissue, y=X.log..p.value., fill = z.score))  + facet_grid(.~Canonical.Pathway) + geom_col(position = "dodge") 

?scale_fill_viridis_c



cp.data$Canonical.Pathway<-reorder(cp.data$Canonical.Pathway, cp.data$X.log..p.value.)
cp.data$Canonical.Pathway<-fct_rev(cp.data$Canonical.Pathway)

tiff("Canonical Pathways Ant, Post, Wing vs Legs Full Names.tiff", units="in", width=7, height=10.5, res=200)
plot2<-ggplot(data=cp.data, aes(x=fct_rev(Tissue), y=X.log..p.value., fill = z.score)) +geom_bar(position = position_dodge(),stat="identity", alpha=0.7) + coord_flip() + scale_fill_gradient2(
  low = muted("red"),
  mid = "white",
  high = muted("blue"),
  midpoint = 0,
  space = "Lab",
  na.value = "grey50",
  guide = "colourbar",
  aesthetics = "fill"
) + 
  labs(y="-log(p-value)", x="Tissue Compartment")
plot2 + facet_grid(Canonical.Pathway~., labeller = label_wrap_gen(width = 10)) + ylim(0, 40) + theme_bw() + geom_hline(yintercept=1.3, linetype="dashed", color="black")

dev.off()


pa.cpdata$Canonical.Pathway<-reorder(pa.cpdata$Canonical.Pathway, pa.cpdata$X.log.p.value.)
#pa.cpdata$Canonical.Pathway<-fct_rev(pa.cpdata$Canonical.Pathway)

tiff("Canonical Pathways Ant vs Post.tiff", units="in", width=8, height=5, res=200)
plot3<-ggplot(data=pa.cpdata, aes(x=Canonical.Pathway, y=X.log.p.value., fill = z.score)) +geom_bar(position = position_dodge(),stat="identity", alpha=0.7) + coord_flip() + scale_fill_gradient2(
  low = muted("red"),
  mid = "white",
  high = muted("blue"),
  midpoint = 0,
  space = "Lab",
  na.value = "grey50",
  guide = "colourbar",
  aesthetics = "fill"
) + 
  labs(y="-log(p-value)", x="Canonical Pathway")
plot3 + ylim(0, 12) + theme_bw() + geom_hline(yintercept=1.3, linetype="dashed", color="black")
dev.off()

