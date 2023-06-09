data<-read.table(file="HxSPhenotypingDataDryad.txt",sep = '\\t',header=TRUE)
library(ggplot2)

data$Plant<-as.character(data$Plant)

# Generate BoxPlots of all phenotypes measured
pdf(file="final/results/HxSPhenotypes_color.pdf",height=4, width=6)
ggplot(data, aes(x=combo, y=Branch_mm, fill=combo)) + geom_boxplot()+ theme_bw()+ scale_fill_manual(values=c("#7B1445","#E7298A", "#33531e","#66a60F", "#0F3C5A","#1F78B4","#735601","#E6AB02"))
ggplot(data, aes(x=combo, y=Nod_Tot, fill=combo)) + geom_boxplot()+ theme_bw()+ scale_fill_manual(values=c("#7B1445","#E7298A", "#33531e","#66a60F", "#0F3C5A","#1F78B4","#735601","#E6AB02"))
ggplot(data, aes(x=combo, y=Biomass_veg, fill=combo)) + geom_boxplot()+ theme_bw()+ scale_fill_manual(values=c("#7B1445","#E7298A", "#33531e","#66a60F", "#0F3C5A","#1F78B4","#735601","#E6AB02"))
ggplot(data, aes(x=combo, y=mg_per_Nod, fill=combo)) + geom_boxplot()+ theme_bw()+ scale_fill_manual(values=c("#7B1445","#E7298A", "#33531e","#66a60F", "#0F3C5A","#1F78B4","#735601","#E6AB02"))
ggplot(data, aes(x=combo, y=Nod_Tot*mg_per_Nod, fill=combo)) + geom_boxplot()+ theme_bw()+ scale_fill_manual(values=c("#7B1445","#E7298A", "#33531e","#66a60F", "#0F3C5A","#1F78B4","#735601","#E6AB02"))
ggplot(data, aes(x=combo, y=Biomass_veg/Nod_Tot, fill=combo)) + geom_boxplot()+ theme_bw()+ scale_fill_manual(values=c("#7B1445","#E7298A", "#33531e","#66a60F", "#0F3C5A","#1F78B4","#735601","#E6AB02"))
dev.off()

#### Run Full Anova models of each trait with Host (Plant), Symbiont (Strain), Strain*Plant as predictors and Tray and Edge position as covariates ####

anova(lm(Nod_Tot~ Plant + Strain + Strain*Plant + Tray + Edge, data=data))
anova(lm(Biomass_veg/Nod_Tot~ Plant + Strain + Strain*Plant + Tray + Edge, data=data))
anova(lm(mg_per_Nod~ Plant + Strain + Strain*Plant + Tray + Edge, data=data))
anova(lm(Biomass_veg~ Plant + Strain + Strain*Plant + Tray + Edge, data=data))

# Subset into strain comparisons within each Host for the four focal traits in the paper

anova(lm(Biomass_veg~ Strain + Tray + Edge, data=data[data$Plant=="56",]))
anova(lm(Biomass_veg~ Strain + Tray + Edge, data=data[data$Plant=="101",]))
anova(lm(Biomass_veg~ Strain + Tray + Edge, data=data[data$Plant=="34",]))
anova(lm(Biomass_veg~ Strain + Tray + Edge, data=data[data$Plant=="340",]))

anova(lm(Nod_Tot~ Strain + Tray + Edge, data=data[data$Plant=="56",]))
anova(lm(Nod_Tot~ Strain + Tray + Edge, data=data[data$Plant=="101",]))
anova(lm(Nod_Tot~ Strain + Tray + Edge, data=data[data$Plant=="34",]))
anova(lm(Nod_Tot~ Strain + Tray + Edge, data=data[data$Plant=="340",]))

anova(lm(mg_per_Nod~ Strain + Tray + Edge, data=data[data$Plant=="56",]))
anova(lm(mg_per_Nod~ Strain + Tray + Edge, data=data[data$Plant=="101",]))
anova(lm(mg_per_Nod~ Strain + Tray + Edge, data=data[data$Plant=="34",]))
anova(lm(mg_per_Nod~ Strain + Tray + Edge, data=data[data$Plant=="340",]))

anova(lm((Biomass_veg/Nod_Tot)~ Strain + Tray + Edge, data=data[data$Plant=="56",]))
anova(lm((Biomass_veg/Nod_Tot)~Strain + Tray + Edge, data=data[data$Plant=="101",]))
anova(lm((Biomass_veg/Nod_Tot)~ Strain + Tray + Edge, data=data[data$Plant=="34",]))
anova(lm((Biomass_veg/Nod_Tot)~ Strain + Tray + Edge, data=data[data$Plant=="340",]))
