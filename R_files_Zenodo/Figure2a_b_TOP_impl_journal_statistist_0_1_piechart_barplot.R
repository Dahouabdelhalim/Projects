#import top guidelines
top=read.csv("C:/Users/ipatarc/Desktop/PR_TOP/top-factor -v33.csv")


#program information why journal names have higher overlap than ISSN


library(dplyr)
library(stringr)
library(ggplot2)

###################
#here is the idea to get a statistics on how many journals in the database are actually having a sum of top guideline scores =0 or >0


    stat_top_g1=(table(top$Total!=0))
    #360 journals have a total of zero top guidelines implementation, whereas 854 have a least one guidelines implemented
    
    names(stat_top_g1)=c("No","Yes") # adjusting the names
    stat_top_g1=as.data.frame(stat_top_g1)
    colnames(stat_top_g1)=c('Implementation',"Freq")
    
    
    # Barplot
    bp<- ggplot(stat_top_g1, aes(x="", y=Freq, fill=Implementation))+
      geom_bar(width = 1, stat = "identity")+
      ggtitle("Adoption of at least one TOP standard")
    bp
    
    blank_theme <- theme_minimal()+
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=24, face="bold"),
        legend.title = element_text(size=24), #change legend title font size
        legend.text = element_text(size=24)
      )
    
    pie <- bp + coord_polar("y")+
      blank_theme +
      theme(axis.text.x=element_blank())
 
    pie
    
    png("C:/Users/ipatarc/Desktop/PR_TOP/top_guidelines_implem_y_n.png"
        ,width = 477,height = 344)
    pie
    dev.off()
    #sizes of the saved image 477/344



########################
# extracting how many policies each journal takes into account


    x=top[1,colnames(top)[str_detect(colnames(top),"score")]]
    
    
    # writing function that will report me number of TOP policies that were adopted
    # by the journal. regardless of the stringency level
    
    n_adopt_policies= function(x){
      x=x[!is.na(x)]
      return(sum(x!=0))
      
    }
    
    Policies_per_journal=apply(top[,colnames(top)[str_detect(colnames(top),"score")]],1,n_adopt_policies)
    
    
    data=as.data.frame(table(Policies_per_journal)[-1])
    
    
    plot=ggplot(data, aes(x=Policies_per_journal, 
                          y=Freq#, 
                          # color="00BBDB",
                          # fill="00BBDB")) + 
    ))+
      geom_bar(stat = "identity")
    
    
    
    png("C:/Users/ipatarc/Desktop/PR_TOP/Adoption_TOP_per_journal_barplot.png"
        ,width = 350,height = 250)
    
    plot + theme_classic()+
      labs(#title = "",
            # subtitle = "",
             #caption = "D",
             x = "Number of adopted TOP standards", y = "Number of journals"#,
             #tag = "B"
             )+ theme(
              # plot.title = element_text(color="red", size=14, face="bold.italic"),
               axis.text=element_text(size=12),
               axis.title.x = element_text(color="black", size=14, face="bold"),
               axis.title.y = element_text(color="black", size=14, face="bold")
             )+theme(legend.position="none")
    
    dev.off()


    
    
    
    
    
###############################################################    
#######
# statistics for the paper
    summary(Policies_per_journal)
    table(Policies_per_journal)
    
    # 0   1   2   3   4   5   6   7   8   9  10 
    # 455 561 263 126 337  59  19  60  78  25  17 
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    #0.00    1.00    1.00    2.33    4.00   10.00 




#######################
# Analysis: what happens to journals that have 4 adopted policies
    journals_with_4_policy=top[which(Policies_per_journal==4),
        colnames(top)[str_detect(colnames(top),"score")]]
    
    
    combinations=t(apply(journals_with_4_policy,1,function(x){which(x!=0)}))
    
    
    str_replace(colnames(journals_with_4_policy),"score","")
    
    combinations[combinations=="10"]="Open science badges"
    combinations[combinations=="9"]="Registered reports/publication bias"
    combinations[combinations=="8"]="Replication"
    combinations[combinations=="7"]="Analysis plan preregistration."
    combinations[combinations=="6"]="Study preregistration"
    combinations[combinations=="5"]="Design analysis reporting guidelines"
    combinations[combinations=="4"]="Materials transparency"
    combinations[combinations=="3"]="Analysis code transparency"
    combinations[combinations=="2"]="Data transparency"
    combinations[combinations=="1"]="Data citation"
    

    #there are 20 combinations
    length(sort(table(apply(combinations,1,paste0,collapse=""))))
    
    # 24500/337= 72.3% combo 1235 which corresponds to "data citation",
    # "data transparency", "analysis code transpareny", "design analysis reporting"
    
    #there are 20 combinations
    as.data.frame(sort(table(apply(combinations,1,paste0,collapse="+"))))
    # this was manually entered in the data table in the article


#######################
# Analysis: what happens to journals that have 1 adopted policy
      journals_with_1_policy=top[which(Policies_per_journal==1),
                                 colnames(top)[str_detect(colnames(top),"score")]]
      
      
      combinations1=t(apply(journals_with_1_policy,1,function(x){which(x!=0)}))
      combinations1[combinations1=="10"]="Open science badges"
      combinations1[combinations1=="9"]="Registered reports/publication bias"
      combinations1[combinations1=="8"]="Replication"
      combinations1[combinations1=="7"]="Analysis plan preregistration."
      combinations1[combinations1=="6"]="Study preregistration"
      combinations1[combinations1=="5"]="Design analysis reporting guidelines"
      combinations1[combinations1=="4"]="Materials transparency"
      combinations1[combinations1=="3"]="Analysis code transparency"
      combinations1[combinations1=="2"]="Data transparency"
      combinations1[combinations1=="1"]="Data citation"
      
      as.data.frame(sort(table(combinations1)))
      
      #as.data.frame(sort(table(combinations1)))
      # combinations1 Freq
      # 1           Analysis code transparency    3
      # 2                          Replication   10
      # 3  Registered reports/publication bias   20
      # 4                  Open science badges   21
      # 5                    Data transparency   51
      # 6 Design analysis reporting guidelines   62
      # 7                        Data citation  394
      
      
      #######################
      # Analysis: what happens to journals that have 2 adopted policy
      journals_with_4_policy=top[which(Policies_per_journal==2),
                                 colnames(top)[str_detect(colnames(top),"score")]]
      
     
      combinations=t(apply(journals_with_4_policy,1,function(x){which(x!=0)}))
      
      
      str_replace(colnames(journals_with_4_policy),"score","")
      
      combinations[combinations=="10"]="Open science badges"
      combinations[combinations=="9"]="Registered reports/publication bias"
      combinations[combinations=="8"]="Replication"
      combinations[combinations=="7"]="Analysis plan preregistration."
      combinations[combinations=="6"]="Study preregistration"
      combinations[combinations=="5"]="Design analysis reporting guidelines"
      combinations[combinations=="4"]="Materials transparency"
      combinations[combinations=="3"]="Analysis code transparency"
      combinations[combinations=="2"]="Data transparency"
      combinations[combinations=="1"]="Data citation"
      
      
      #there are 20 combinations
      length(sort(table(apply(combinations,1,paste0,collapse=""))))
      
      # 24500/337= 72.3% combo 1235 which corresponds to "data citation",
      # "data transparency", "analysis code transpareny", "design analysis reporting"
      
      #there are 20 combinations
      as.data.frame(sort(table(apply(combinations,1,paste0,collapse="+"))))  
      
      
      
      # Var1 Freq
      # 1        Analysis plan preregistration.+Registered reports/publication bias    1
      # 2                                  Data transparency+Materials transparency    1
      # 3                     Data transparency+Registered reports/publication bias    1
      # 4                                             Data transparency+Replication    1
      # 5       Design analysis reporting guidelines+Analysis plan preregistration.    1
      # 6  Design analysis reporting guidelines+Registered reports/publication bias    1
      # 7               Materials transparency+Design analysis reporting guidelines    1
      # 8                                           Replication+Open science badges    1
      # 9                  Design analysis reporting guidelines+Open science badges    3
      # 10                         Design analysis reporting guidelines+Replication    3
      # 11                                    Data transparency+Open science badges    4
      # 12                  Registered reports/publication bias+Open science badges    4
      # 13                                                Data citation+Replication    5
      # 14                                        Data citation+Open science badges    6
      # 15                        Data citation+Registered reports/publication bias    8
      # 16                          Replication+Registered reports/publication bias   13
      # 17                             Data transparency+Analysis code transparency   16
      # 18                   Data transparency+Design analysis reporting guidelines   21
      # 19                       Data citation+Design analysis reporting guidelines   49
      # 20                                          Data citation+Data transparency  123
      
      # nrow(journals_with_4_policy)-123-49-21-16-13
      
      #######################
      # Analysis: what happens to journals that have 3 adopted policy
      journals_with_4_policy=top[which(Policies_per_journal==3),
                                 colnames(top)[str_detect(colnames(top),"score")]]
      
      
      combinations=t(apply(journals_with_4_policy,1,function(x){which(x!=0)}))
      
      
      str_replace(colnames(journals_with_4_policy),"score","")
      
      combinations[combinations=="10"]="Open science badges"
      combinations[combinations=="9"]="Registered reports/publication bias"
      combinations[combinations=="8"]="Replication"
      combinations[combinations=="7"]="Analysis plan preregistration."
      combinations[combinations=="6"]="Study preregistration"
      combinations[combinations=="5"]="Design analysis reporting guidelines"
      combinations[combinations=="4"]="Materials transparency"
      combinations[combinations=="3"]="Analysis code transparency"
      combinations[combinations=="2"]="Data transparency"
      combinations[combinations=="1"]="Data citation"
      
      
      #there are 20 combinations
      length(sort(table(apply(combinations,1,paste0,collapse=""))))
      
      # 24500/337= 72.3% combo 1235 which corresponds to "data citation",
      # "data transparency", "analysis code transpareny", "design analysis reporting"
      
      #there are 20 combinations
      as.data.frame(sort(table(apply(combinations,1,paste0,collapse="+"))))  
      
      
      
      # Var1 Freq
      # 1                           Data citation+Data transparency+Registered reports/publication bias    1
      # 2                        Data citation+Design analysis reporting guidelines+Open science badges    1
      # 3              Data transparency+Analysis code transparency+Registered reports/publication bias    1
      # 4                    Data transparency+Design analysis reporting guidelines+Open science badges    1
      # 5                  Data transparency+Design analysis reporting guidelines+Study preregistration    1
      # 6                  Data transparency+Materials transparency+Registered reports/publication bias    1
      # 7                     Data transparency+Registered reports/publication bias+Open science badges    1
      # 8                                             Data transparency+Replication+Open science badges    1
      # 9     Design analysis reporting guidelines+Study preregistration+Analysis plan preregistration.    1
      # 10                       Design analysis reporting guidelines+Study preregistration+Replication    1
      # 11     Study preregistration+Analysis plan preregistration.+Registered reports/publication bias    1
      # 12                               Data citation+Design analysis reporting guidelines+Replication    2
      # 13 Design analysis reporting guidelines+Registered reports/publication bias+Open science badges    2
      # 14                         Design analysis reporting guidelines+Replication+Open science badges    2
      # 15                          Replication+Registered reports/publication bias+Open science badges    2
      # 16   Data transparency+Design analysis reporting guidelines+Registered reports/publication bias    3
      # 17                Data transparency+Materials transparency+Design analysis reporting guidelines    3
      # 18         Design analysis reporting guidelines+Replication+Registered reports/publication bias    3
      # 19                                   Data citation+Data transparency+Analysis code transparency    4
      # 20       Data citation+Design analysis reporting guidelines+Registered reports/publication bias    5
      # 21            Data transparency+Analysis code transparency+Design analysis reporting guidelines    5
      # 22                           Data transparency+Analysis code transparency+Study preregistration    5
      # 23                                Data citation+Replication+Registered reports/publication bias    6
      # 24                            Data transparency+Replication+Registered reports/publication bias    6
      # 25                                          Data citation+Data transparency+Open science badges    9
      # 26                          Data transparency+Analysis code transparency+Materials transparency   27
      # 27                         Data citation+Data transparency+Design analysis reporting guidelines   31
      
      
