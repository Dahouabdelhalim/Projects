##########################################
#Pie chart top guidelines evaluated per db



###########################
# import
Overlap_SCOPUS_TOP=readRDS("C:/Users/ipatarc/Desktop/PR_TOP/Overlap_SCOPUS_TOP.rds")
library("cowplot")
library("gridExtra")
library(ggridges)
library(pheatmap)
library(ggplot2)
library(stringr)



########
#1 creating Multidisciplinary journal category
#reviewing top category of sciences in the top guidelines
    Sci_categories=c("Top.level...Life.Sciences","Top.level...Social.Sciences","Top.level...Physical.Sciences",
                     "Top.level...Health.Sciences")
    table(apply(Overlap_SCOPUS_TOP[,Sci_categories],1,paste0,collapse = ""))
    # RESULTS: health science 105, life science 67, social science 423, phys science 54, NA 10
    
    #hard coding categories of science
    Overlap_SCOPUS_TOP$SCIENCE <- str_replace_all(apply(Overlap_SCOPUS_TOP[,Sci_categories],1,paste0,collapse = ""),"Sciences","")
    
    Overlap_SCOPUS_TOP$SCIENCE[Overlap_SCOPUS_TOP$SCIENCE%in%c("Social Health ",
                                          "Life Social ",
                                          "Life Social Health ",
                                          "Social Physical ",
                                          "Social Physical Health ",
                                          "Life Social Physical ",
                                          "Life Social Physical Health ",
                                          "Life Health ",
                                          "Life Physical ",
                                          "Physical Health ",
                                          "Life Physical Health ")]="Multidisc."
    
    # Overlap_SCOPUS_TOP$SCIENCE[Overlap_SCOPUS_TOP$SCIENCE%in%c("Life Health ","Life Physical ","Physical Health ","Life Physical Health ")]="Multidisc.(-Social)"

#results statistics
    #creating 6 categories
    table(Overlap_SCOPUS_TOP$SCIENCE)
    
    # Health       Life  Multidisc.  Physical  
    # 10        297        162        642        118 
    # Social  
    # 595
   
    round(100*table(Overlap_SCOPUS_TOP$SCIENCE)/1824)
    # Health       Life  Multidisc.  Physical  
    # 1         16          9         35          6 
    # Social  
    # 33 
    
     
#removing 10 journals that do not have affiliation with field of science
    Overlap_SCOPUS_TOP=Overlap_SCOPUS_TOP[Overlap_SCOPUS_TOP$SCIENCE!="",]

    Sci_field_overview=as.data.frame(table(Overlap_SCOPUS_TOP$SCIENCE))
    colnames(Sci_field_overview)=c("Science","Freq")


#2. plotting results - PIE chart
      #saveRDS(Overlap_SCOPUS_TOP,"C:/Users/ipatarc/Desktop/PR_TOP/Overlap_SCOPUS_TOP_SCIENCE_added.rds")
      # Barplot
      bp<- ggplot(Sci_field_overview, aes(x="", y=Freq, fill=Science))+
        geom_bar(width = 1, stat = "identity")+
        ggtitle("Journals in the osf.io across science disciplines")
      bp
      
      blank_theme <- theme_minimal()+
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank(),
          axis.ticks = element_blank(),
          plot.title=element_text(size=14, face="bold")
        )
      
      pie <- bp + coord_polar("y")+
        blank_theme +
        theme(axis.text.x=element_blank())
      pie
      
      png("C:/Users/ipatarc/Desktop/PR_TOP/TOP_guidelines_Across_Journalspie.png"
          ,width = 477,height = 344)
      pie
      dev.off()
      #sizes#sizes of the saved image 477/344
      
      
      

######################################
# IDEA comparing results with the Scopus database - checking if sample equals (TOP Factor) to population (Scopus database)
      
#import SCOPUS db
      SCOPUS_content=read.csv("C:/Users/ipatarc/Desktop/PR_TOP/SCOPUS_content_July2022.csv",
                              header = T) # number of rows seems to be ok!
      
      SCOPUS_content_active=SCOPUS_content[SCOPUS_content$Active.or.Inactive=="Active",]
      dim(SCOPUS_content) #43016    53 number of all journals
      dim(SCOPUS_content_active) #27253    53 number of active journals
      
      Sci_categories=c("Top.level...Life.Sciences","Top.level...Social.Sciences",
      "Top.level...Physical.Sciences","Top.level...Health.Sciences" )
      
      SCOPUS_category <- str_replace_all(apply(SCOPUS_content_active[,Sci_categories],1,paste0,collapse = ""),"Sciences","")
      
      SCOPUS_category[SCOPUS_category%in%c("Social Health ","Life Social ",
                                                                 "Life Social Health ",
                                                                 "Social Physical ",
                                                                 "Social Physical Health ",
                                                                 "Life Social Physical ",
                                                                 "Life Social Physical Health ",
                                           "Life Health ","Life Physical ","Physical Health ","Life Physical Health ")]="Multidisc."
      
      # SCOPUS_category[SCOPUS_category%in%c("Life Health ","Life Physical ","Physical Health ","Life Physical Health ")]="Multidisc.(-Social)"

#creating 6 categories
    table(SCOPUS_category)
    
    # SCOPUS_category
    # Health       Life  Multidisc.  Physical     Social  
    # 93       4603       1937       5973       5848       8799
    # SCOPUS_category

        # Health       Life  Multidisc.  Physical     Social  
    # 93       4603       1937       5973       5848       8799 
round(100*table(SCOPUS_category)/length(SCOPUS_category))

# SCOPUS_category
# Health       Life  Multidisc.  Physical     Social  
# 93       4603       1937       5973       5848       8799 

# 2. calculating chi-square test to see if distiributions of sample and data are the same

    x=table(Overlap_SCOPUS_TOP$SCIENCE)
    y=table(SCOPUS_category)
    y=y[!names(y)==""] #removing category without name
    df=as.data.frame(cbind(x,y))
    chisq.test(df)
    
    # Pearson's Chi-squared test
    # 
    # data:  df
    # X-squared = 353.09, df = 5, p-value < 2.2e-16
    
    df=as.data.frame(cbind(100*x/sum(x),100*y/sum(y)))
    colnames(df)=c("TOP_guidelines","Scopus")
    df$SCIENCE=rownames(df)
    df

##################################
    
#3. plot
    
    colnames(Overlap_SCOPUS_TOP)
    
    
    policies=which(str_detect(colnames(Overlap_SCOPUS_TOP),"score"))
    
    Discipline=unique(Overlap_SCOPUS_TOP$SCIENCE)
   
    # creating data for statistics
        for_list=list()
        
         for (i in (Discipline)){
         
          filtered.j=Overlap_SCOPUS_TOP[Overlap_SCOPUS_TOP$SCIENCE==i,policies]
          
               per_score_policy=rbind(apply(filtered.j,2,function(x){length(which(x==1))}),
                                 apply(filtered.j,2,function(x){length(which(x==2))}),
                                 apply(filtered.j,2,function(x){length(which(x==3))}))
               
               per_score_policy=as.data.frame(colSums(per_score_policy))
          
               per_score_policy$Discipline=i
               per_score_policy$TOP=rownames(per_score_policy)
              colnames(per_score_policy)[1]="value"
               for_list[[i]]=per_score_policy
        }
        
          Per_discipline_per_guideline_stat=do.call("rbind", for_list)
          
          
          
          # value Discipline                                          TOP
          # Multidisc..Data.citation.score                            422 Multidisc.                          Data.citation.score
          # Multidisc..Data.transparency.score                        277 Multidisc.                      Data.transparency.score
          # Multidisc..Analysis.code.transparency.score               181 Multidisc.             Analysis.code.transparency.score
          # Multidisc..Materials.transparency.score                    81 Multidisc.                 Materials.transparency.score
          # Multidisc..Design...analysis.reporting.guidelines.score   252 Multidisc. Design...analysis.reporting.guidelines.score
          # Multidisc..Study.preregistration.score                     56 Multidisc.                  Study.preregistration.score
          # Multidisc..Analysis.plan.preregistration.score             55 Multidisc.          Analysis.plan.preregistration.score
          # Multidisc..Replication.score                               94 Multidisc.                            Replication.score
          # Multidisc..Registered.reports...publication.bias.score     56 Multidisc.  Registered.reports...publication.bias.score
          # Multidisc..Open.science.badges.score                       31 Multidisc.                    Open.science.badges.score
          # Social .Data.citation.score                               259    Social                           Data.citation.score
          # Social .Data.transparency.score                           220    Social                       Data.transparency.score
          # Social .Analysis.code.transparency.score                  107    Social              Analysis.code.transparency.score
          # Social .Materials.transparency.score                       81    Social                  Materials.transparency.score
          # Social .Design...analysis.reporting.guidelines.score      101    Social  Design...analysis.reporting.guidelines.score
          # Social .Study.preregistration.score                        52    Social                   Study.preregistration.score
          # Social .Analysis.plan.preregistration.score                42    Social           Analysis.plan.preregistration.score
          # Social .Replication.score                                  87    Social                             Replication.score
          # Social .Registered.reports...publication.bias.score        73    Social   Registered.reports...publication.bias.score
          # Social .Open.science.badges.score                          55    Social                     Open.science.badges.score
          # Physical .Data.citation.score                              85  Physical                           Data.citation.score
          # Physical .Data.transparency.score                          61  Physical                       Data.transparency.score
          # Physical .Analysis.code.transparency.score                 36  Physical              Analysis.code.transparency.score
          # Physical .Materials.transparency.score                     25  Physical                  Materials.transparency.score
          # Physical .Design...analysis.reporting.guidelines.score     25  Physical  Design...analysis.reporting.guidelines.score
          # Physical .Study.preregistration.score                      15  Physical                   Study.preregistration.score
          # Physical .Analysis.plan.preregistration.score              15  Physical           Analysis.plan.preregistration.score
          # Physical .Replication.score                                16  Physical                             Replication.score
          # Physical .Registered.reports...publication.bias.score       2  Physical   Registered.reports...publication.bias.score
          # Physical .Open.science.badges.score                         7  Physical                     Open.science.badges.score
          # Life .Data.citation.score                                 131      Life                           Data.citation.score
          # Life .Data.transparency.score                              94      Life                       Data.transparency.score
          # Life .Analysis.code.transparency.score                     70      Life              Analysis.code.transparency.score
          # Life .Materials.transparency.score                         37      Life                  Materials.transparency.score
          # Life .Design...analysis.reporting.guidelines.score         84      Life  Design...analysis.reporting.guidelines.score
          # Life .Study.preregistration.score                          30      Life                   Study.preregistration.score
          # Life .Analysis.plan.preregistration.score                  28      Life           Analysis.plan.preregistration.score
          # Life .Replication.score                                    34      Life                             Replication.score
          # Life .Registered.reports...publication.bias.score          10      Life   Registered.reports...publication.bias.score
          # Life .Open.science.badges.score                             6      Life                     Open.science.badges.score
          # Health .Data.citation.score                               222    Health                           Data.citation.score
          # Health .Data.transparency.score                           178    Health                       Data.transparency.score
          # Health .Analysis.code.transparency.score                  122    Health              Analysis.code.transparency.score
          # Health .Materials.transparency.score                       24    Health                  Materials.transparency.score
          # Health .Design...analysis.reporting.guidelines.score      184    Health  Design...analysis.reporting.guidelines.score
          # Health .Study.preregistration.score                        21    Health                   Study.preregistration.score
          # Health .Analysis.plan.preregistration.score                21    Health           Analysis.plan.preregistration.score
          # Health .Replication.score                                  30    Health                             Replication.score
          # Health .Registered.reports...publication.bias.score        16    Health   Registered.reports...publication.bias.score
          # Health .Open.science.badges.score                           5    Health                     Open.science.badges.score
          #
          
          
          #############
          # statistics for table - per discipline per standard statistics
          Per_discipline_per_guideline_stat.c=do.call("cbind", for_list)
          t=Per_discipline_per_guideline_stat.c[,c(1,4,7,10,13)]
          
          
          
          library(pheatmap)
          
          (apply(t,2,function(x)round(100*x/sum(x))))
          
          # Multidisc..value Social .value
          # Data.citation.score                                        28            24
          # Data.transparency.score                                    18            20
          # Analysis.code.transparency.score                           12            10
          # Materials.transparency.score                                5             8
          # Design...analysis.reporting.guidelines.score               17             9
          # Study.preregistration.score                                 4             5
          # Analysis.plan.preregistration.score                         4             4
          # Replication.score                                           6             8
          # Registered.reports...publication.bias.score                 4             7
          # Open.science.badges.score                                   2             5
          # Physical .value Life .value
          # Data.citation.score                                       30          25
          # Data.transparency.score                                   21          18
          # Analysis.code.transparency.score                          13          13
          # Materials.transparency.score                               9           7
          # Design...analysis.reporting.guidelines.score               9          16
          # Study.preregistration.score                                5           6
          # Analysis.plan.preregistration.score                        5           5
          # Replication.score                                          6           6
          # Registered.reports...publication.bias.score                1           2
          # Open.science.badges.score                                  2           1
          # Health .value
          # Data.citation.score                                     27
          # Data.transparency.score                                 22
          # Analysis.code.transparency.score                        15
          # Materials.transparency.score                             3
          # Design...analysis.reporting.guidelines.score            22
          # Study.preregistration.score                              3
          # Analysis.plan.preregistration.score                      3
          # Replication.score                                        4
          # Registered.reports...publication.bias.score              2
          # Open.science.badges.score                                1
          
          for.heatmap=(apply(t,2,function(x)round(100*x/sum(x))))
          colnames(for.heatmap)=str_replace(colnames(for.heatmap)," \\\\.value","")
          rownames(for.heatmap)=str_replace(rownames(for.heatmap)," \\\\.score","")
          rownames(for.heatmap)=str_replace(rownames(for.heatmap),"\\\\.score","")
          rownames(for.heatmap)=str_replace(rownames(for.heatmap),"\\\\."," ")
          rownames(for.heatmap)[rownames(for.heatmap)=="Registered reports...publication.bias"]="Registered reports"
          rownames(for.heatmap)[rownames(for.heatmap)=="Design ..analysis.reporting.guidelines"]="Design analysis reporting.guidelines"
          colnames(for.heatmap)[1]="Multidisc."
           
          png("C:/Users/ipatarc/Desktop/PR_TOP/Supp_Figure_Heatmap.png", 
              width =700,height = 900)
          
           pheatmap( for.heatmap,display_numbers = T,cluster_rows = F, fontsize_row=20,fontsize_col = 20,fontsize_number = 17)
          
          dev.off()
          
          # DOING cosmetics
              colnames(Per_discipline_per_guideline_stat)=c("Journals","Discipline","TOP_policy")
              
              Per_discipline_per_guideline_stat$TOP_policy=str_replace(Per_discipline_per_guideline_stat$TOP_policy,
                                                                       "...publication.bias.score","")
              
              Per_discipline_per_guideline_stat$TOP_policy=str_replace(Per_discipline_per_guideline_stat$TOP_policy,
                                                                       ".score","")
              
              Per_discipline_per_guideline_stat$TOP_policy=str_replace(Per_discipline_per_guideline_stat$TOP_policy,
                                                                       "\\\\."," ")
              Per_discipline_per_guideline_stat$TOP_policy=str_replace(Per_discipline_per_guideline_stat$TOP_policy,
                                                                       "\\\\.\\\\.","\\n")
              
              Per_discipline_per_guideline_stat$TOP_policy=str_replace(Per_discipline_per_guideline_stat$TOP_policy,
                                                                       "\\\\.","\\n")
              
              Per_discipline_per_guideline_stat$TOP_policy=paste0(rep(LETTERS[10:1],each=1),
                                                                  ".",Per_discipline_per_guideline_stat$TOP_policy)
              
              # Per_discipline_per_guideline_stat$Level=1:3
              
              
 png("C:/Users/ipatarc/Desktop/PR_TOP/Fig3_TOP_policies_across_TOP_disciplines_barplot.png", 
     width =800,height = 1000)
    
    ggplot(data = Per_discipline_per_guideline_stat)+
      geom_bar(aes(x = Journals,
                   # value=as.factor(Discipline),
                   y=as.factor(TOP_policy),
                   fill = as.factor(Discipline)),  
               position = "dodge",
               stat="identity")  +
      #scale_fill_discrete(breaks=c("1","2","3"))+  
      theme_bw()+
      labs(y="Standards of the TOP Guidelines",
           x="Number of journals")+#,
          #title="Number of journals within each requirement level of the TOP policy") + 
      theme(axis.text=element_text(size=16,face="bold"),
            axis.title=element_text(size=20,face="bold"),
            plot.title=element_text(size=16, face="bold"),
            legend.title = element_text(size=18, face="bold"), 
            legend.text = element_text(size=18, face="bold"))+
      guides(fill=guide_legend(title="Disciplines"))
      
   dev.off()
    
    
    
    
   
   
   
   
   
    
###############################################################
#############################################################
# OPTIONAL ANALYSIS: plotting density plots for Scores for TOP Guidelines


    png("C:/Users/ipatarc/Desktop/top_guidelines_distributions.png"
        ,width = 900,height = 900)
    plot3 <- ggplot(Overlap_SCOPUS_TOP, aes(x=Total, fill = SCIENCE, y = SCIENCE)) + 
      geom_density_ridges() + labs(x="\\n Max score (TOP guidelines)") +
      scale_x_continuous()+#trans = "log10")+
      theme(legend.position="none",
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 90),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.title = element_text(size=40,face="bold"),
            axis.text =element_text(size=25,face="bold"),
            text = element_text(size=25,face="bold", vjust = 2))
    plot3
    dev.off()
    
    
    stat_per_science=t(table(Overlap_SCOPUS_TOP$Total,Overlap_SCOPUS_TOP$SCIENCE))
    stat_per_science_per=round(100*stat_per_science/rowSums(stat_per_science))
    pheatmap(stat_per_science_per,cluster_rows = F, cluster_cols = F,display_numbers = T)



#############################################################
# OPTIONAL ANALYSIS: plotting density plots for Scores for TOP Guidelines



    plot3 <- ggplot(Overlap_SCOPUS_TOP, aes(x=X2021.CiteScore, fill = SCIENCE, y = SCIENCE)) + 
      geom_density_ridges() + labs(x="\\n Max score (TOP guidelines)") +
      scale_x_continuous()+#trans = "log10")+
      theme(legend.position="none",
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 90),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.title = element_text(size=40,face="bold"),
            axis.text =element_text(size=25,face="bold"),
            text = element_text(size=25,face="bold", vjust = 2))
    plot3
    
    
    
#############################################################
# OPTIONAL ANALYSIS: how does citation factor regresses with TOP score

ggplot(Overlap_SCOPUS_TOP, aes(Total,X2021.CiteScore )) + 
  geom_point() + facet_grid(. ~ SCIENCE) + stat_smooth(method = "lm") #+
 # background_grid(major = &#39;y&#39;, minor = "none") + # add thin horizontal lines 
                    panel_border() # and a border around each panel
                  # plot.mpt and plot.diamonds were defined earlier
                  ggdraw() +
                    draw_plot(plot.iris, 0, .5, 1, .5) +
                    draw_plot(sp, 0, 0, .5, .5) +
                    draw_plot(bp, .5, 0, .5, .5) +
                    draw_plot_label(c("A", "B", "C"), c(0, 0, 0.5), c(1, 0.5, 0.5), size = 15)
                  
                  
                  
                  
                  
                  sort(table(Overlap_SCOPUS_TOP$Publisher))