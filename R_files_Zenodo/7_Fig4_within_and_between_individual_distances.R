# Perofsky et al. "Social groups constrain the spatiotemporal dynamics of wild sifaka gut microbiomes"
# Contact: Amanda Perofsky, amanda.perofsky@nih.gov

rm(list=ls())
graphics.off()

## set your working directory to where files are located
## input files are in "Rdata" folder
## figures are saved in "figures" folder
dir <- "~/OneDrive - National Institutes of Health/NIH_Laptop_Updates_Post_Damage/Documents/Sifaka_KMNP_2016/Final_Code/Cleaned_Scripts/" ## set your working directory to where files are located
## set your working directory to where files are located
setwd(dir)

library(ggplot2)
library(cowplot)
library(dplyr)
library(rcompanion)
library(ggpubr)
library(ggpubfigs)

red_df = readRDS("Rdata/within_and_bw_host_df_lim.rds")
head(red_df)
names(red_df)
#######################################################################
## Figure 4: within and between individual pairwise distances, within and across years
#######################################################################
my_comparisons = compare_means(bray_curtis ~ Value,  data = red_df,p.adjust.method="BH")

PT = pairwisePermutationTest(bray_curtis ~ Value,
                             data   = red_df,
                             method = "BH")

pwc <- red_df %>% 
  rstatix::dunn_test(bray_curtis ~ Value, p.adjust.method = "BH") 
pwc
pwc$group_name = paste(pwc$group1,pwc$group2,sep="-")

cld_fun <- function(sub_df){
  ret_df = cldList(p.adj ~ group_name,data = sub_df,threshold = 0.05)
  return(ret_df)
}
dunn_df = cld_fun(pwc)
names(dunn_df)[1]<-"Value"

max_values = red_df  %>%
  summarize(max_y = max(bray_curtis)+0.1)

dunn_df$max_y = max_values$max_y 

pal = friendly_pal("muted_nine",9)


p <- ggplot(red_df,aes(x=Value,y=bray_curtis,fill=Value))+
  geom_violin(trim=F,alpha=0.7)+
  geom_boxplot(width=0.1) + 
  ylab("Microbial Dissimilarity\\n(Bray-Curtis)")+
  geom_text(data = dunn_df, aes(y = max_y, label = Letter), 
            position = position_dodge(width = .75),size=5, color="black")+
  xlab("")+
  theme_cowplot(font_size = 14)+
  background_grid("xy")+
  scale_x_discrete(breaks=c("Same_Ind_Same_Year_Same_Group",
                            "Diff_Ind_Same_Year_Same_Group",
                            "Diff_Ind_Same_Year_Diff_Group",
                            "Same_Ind_Diff_Year_Same_Group",
                            "Diff_Ind_Diff_Year_Same_Group",
                            "Diff_Ind_Diff_Year_Diff_Group"),
                   labels=c("Same individual \\n Within season", 
                            "Group members \\n Within season", 
                            "Non-group members \\n Within season",
                            "Same individual \\n Different Years",
                            "Group members \\n Different Years",
                            "Non-group members \\n Different Years")) + # ggtitle("Pairwise sample comparisons, all years")
  # scale_fill_manual(values=c("#9A8822" ,"#9A8822","#74A089","#74A089", "#F8AFA8" ,"#F8AFA8"))+
  scale_fill_manual(values=c("#CC6677","#CC6677", "#88CCEE" , "#88CCEE" ,"#999933","#999933"))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))
p
save_plot(p, filename = "figures/Fig4_individual_vs_group_and_non_group_members_pairwise_comparisons_all_years.pdf",base_width = 12,base_height = 5)
save_plot(p, filename = "figures/Fig4_individual_vs_group_and_non_group_members_pairwise_comparisons_all_years.png",dpi=600,base_width = 12,base_height = 5)

#### weighted unifrac

pwc <- red_df %>% 
  rstatix::dunn_test(w_unifrac ~ Value, p.adjust.method = "BH") 
pwc
pwc$group_name = paste(pwc$group1,pwc$group2,sep="-")

cld_fun <- function(sub_df){
  ret_df = cldList(p.adj ~ group_name,data = sub_df,threshold = 0.05)
  return(ret_df)
}
dunn_df = cld_fun(pwc)
names(dunn_df)[1]<-"Value"

max_values = red_df  %>%
  summarize(max_y = max(w_unifrac)+0.1)

dunn_df$max_y = max_values$max_y 

p <- ggplot(red_df,aes(x=Value,y=w_unifrac,fill=Value))+
  geom_violin(trim=F,alpha=0.7)+
  geom_boxplot(width=0.1) + 
  ylab("Microbial Dissimilarity\\n(Weighted Unifrac)")+
  geom_text(data = dunn_df, aes(y = max_y, label = Letter), 
            position = position_dodge(width = .75),size=5, color="black")+
  xlab("")+
  theme_cowplot(font_size = 14)+
  background_grid("xy")+
  scale_x_discrete(breaks=c("Same_Ind_Same_Year_Same_Group",
                            "Diff_Ind_Same_Year_Same_Group",
                            "Diff_Ind_Same_Year_Diff_Group",
                            "Same_Ind_Diff_Year_Same_Group",
                            "Diff_Ind_Diff_Year_Same_Group",
                            "Diff_Ind_Diff_Year_Diff_Group"),
                   labels=c("Same individual \\n Within season", 
                            "Group members \\n Within season", 
                            "Non-group members \\n Within season",
                            "Same individual \\n Different Years",
                            "Group members \\n Different Years",
                            "Non-group members \\n Different Years")) + # ggtitle("Pairwise sample comparisons, all years")
  scale_fill_manual(values=c("#9A8822" ,"#9A8822","#74A089","#74A089", "#F8AFA8" ,"#F8AFA8"))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))
p

### differences between samples collected in 2015 and 2016 (one year apart)
red_df_one_year = red_df %>% filter(Year1 %in% c("2015","2016") & Year2 %in%c("2015","2016")) %>% droplevels()
df_comp = compare_means(bray_curtis ~ Value,  data = red_df_one_year,method="wilcox.test", p.adjust.method = "BH")
df_comp %>% filter(p.adj>0.05)

PT = pairwisePermutationTest(bray_curtis ~ Value,
                             data   = red_df_one_year,
                             method = "BH")
PT

pwc <- red_df_one_year %>% 
  rstatix::dunn_test(bray_curtis ~ Value, p.adjust.method = "BH") 
pwc
pwc$group_name = paste(pwc$group1,pwc$group2,sep="-")


cld_fun <- function(sub_df){
  ret_df = cldList(p.adj ~ group_name,data = sub_df,threshold = 0.05)
  return(ret_df)
}
dunn_df = cld_fun(pwc)
names(dunn_df)[1]<-"Value"

max_values = red_df_one_year  %>%
  summarize(max_y = max(bray_curtis)+0.1)

dunn_df$max_y = max_values$max_y 

p <- ggplot(red_df_one_year,aes(x=Value,y=bray_curtis,fill=Value))+
  geom_violin(trim=F,alpha=0.7)+
  geom_boxplot(width=0.1) + 
  geom_text(data = dunn_df, aes(y = max_y, label = Letter), 
            position = position_dodge(width = .75),size=5, color="black")+
  ylab("Microbial Dissimilarity\\n(Bray-Curtis)")+
  xlab("")+
  theme_cowplot()+
  background_grid("xy")+
  scale_x_discrete(breaks=c("Same_Ind_Same_Year_Same_Group",
                            "Diff_Ind_Same_Year_Same_Group",
                            "Diff_Ind_Same_Year_Diff_Group",
                            "Same_Ind_Diff_Year_Same_Group",
                            "Diff_Ind_Diff_Year_Same_Group",
                            "Diff_Ind_Diff_Year_Diff_Group"),
                   labels=c("Same individual \\n Within season", 
                            "Group members \\n Within season", 
                            "Non-group members \\n Within season",
                            "Same individual \\n Different Years",
                            "Group members \\n Different Years",
                            "Non-group members \\n Different Years")) + # ggtitle("Pairwise sample comparisons, all years")
  scale_fill_manual(values=c("#9A8822" ,"#9A8822","#74A089","#74A089", "#F8AFA8" ,"#F8AFA8"))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))
p

pwc <- red_df_one_year %>% 
  rstatix::dunn_test(w_unifrac ~ Value, p.adjust.method = "BH") 
pwc
pwc$group_name = paste(pwc$group1,pwc$group2,sep="-")

cld_fun <- function(sub_df){
  ret_df = cldList(p.adj ~ group_name,data = sub_df,threshold = 0.05)
  return(ret_df)
}
dunn_df = cld_fun(pwc)
names(dunn_df)[1]<-"Value"

max_values = red_df_one_year  %>%
  summarize(max_y = max(w_unifrac)+0.1)

dunn_df$max_y = max_values$max_y 

q <- ggplot(red_df_one_year,aes(x=Value,y=w_unifrac,fill=Value))+
  geom_violin(trim=F,alpha=0.7)+
  geom_boxplot(width=0.1) + 
  geom_text(data = dunn_df, aes(y = max_y, label = Letter), 
            position = position_dodge(width = .75),size=5, color="black")+
  ylab("Microbial Dissimilarity\\n(Weighted Unifrac)")+
  xlab("")+
  theme_cowplot()+
  background_grid("xy")+
  scale_x_discrete(breaks=c("Same_Ind_Same_Year_Same_Group",
                            "Diff_Ind_Same_Year_Same_Group",
                            "Diff_Ind_Same_Year_Diff_Group",
                            "Same_Ind_Diff_Year_Same_Group",
                            "Diff_Ind_Diff_Year_Same_Group",
                            "Diff_Ind_Diff_Year_Diff_Group"),
                   labels=c("Same individual \\n Within season", 
                            "Group members \\n Within season", 
                            "Non-group members \\n Within season",
                            "Same individual \\n Different Years",
                            "Group members \\n Different Years",
                            "Non-group members \\n Different Years")) + # ggtitle("Pairwise sample comparisons, all years")
  scale_fill_manual(values=c("#9A8822" ,"#9A8822","#74A089","#74A089", "#F8AFA8" ,"#F8AFA8"))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))
q

one_year_grid = plot_grid(p,q,labels = c("A","B"),nrow=2)
one_year_grid
