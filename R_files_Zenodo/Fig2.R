library(openxlsx)
library(ggplot2)
library(data.table)

cits = as.data.table(read.csv("Fig2-Hardin_citations.csv"))
cits = cits[Publication.Year < 2021,]
#to count citations later
cits[, row:= 1]

#when is coexist/coexisting/coexistence mentioned in these papers? (In either title or abstract)
cits[,titleabstract := paste(Article.Title," ", Abstract)]
table(grepl("coexist", cits$titleabstract, ignore.case = T))

#these are the papers citing Hardin since 2000 that mention "coexist"
cits1 = cits[!is.na(Abstract)&grepl("coexist", titleabstract, ignore.case = T)]
cits1 = cits1[order(Publication.Year), .(sum_coexist = sum(row)), by = Publication.Year]
cits1[, cumsum_coexist := cumsum(sum_coexist)]

#this is the cumulative sum of all papers since 2000 citing Hardin
cits2 = cits[order(Publication.Year), .(sum_all = sum(row)), by = Publication.Year]
cits2[, cumsum_all := cumsum(sum_all)]

###merge 2 groups and melt to create factors with legend in ggplot
cits3 = merge(cits1, cits2, by = 'Publication.Year')
cits3 = cits3[,list(Publication.Year, cumsum_all, cumsum_coexist)]
colnames(cits3) = c("Publication.Year", "All", "About coexistence")
cits3 = melt(cits3, id = "Publication.Year")
#Add cumulative sum labels to last value in each group
cits3[, label := NA_integer_]
cits3[Publication.Year == 2020 & variable == "All", label := 1359]
cits3[Publication.Year == 2020 & variable == "About coexistence", label := 438]

#cumulative sum figure
ggsave(paste0(Sys.Date(), "_Hardin-citations-cumulative.pdf"), width = 3.15, height = 3.15, units = "in", dpi = 800)
ggplot(data = cits3, aes(x = Publication.Year, y = value, color = factor(variable)))+
  scale_color_manual(values = c("black", "gray50"))+ #, guide = guide_legend(reverse = TRUE))+
  geom_step(size = 1)+
  theme_bw()+
  theme(text = element_text(size = 11), axis.text = element_text(size = 11))+
  xlab("Publication Year")+ylab("Cumulative Citations")+
  theme(legend.position = "bottom", legend.text = element_text(size = 11))+
  labs(color = "")+ 
  geom_text(data = cits3, aes(label = label), nudge_x = -0.7, nudge_y = 70, 
             na.rm = T, show.legend = F)
dev.off()
