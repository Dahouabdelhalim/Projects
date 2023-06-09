#################
# This is a script that overlaps TOP guidelines df and SCOPUS database
#it does it in two steps: 1) journal name and 2) print and e.issn
# then it rinds results and selects only some columns of interest and saves the document


#import SCOPUS db
SCOPUS_content=read.csv("C:/Users/ipatarc/Desktop/PR_TOP/SCOPUS_content_July2022.csv",
                        header = T) # number of rows seems to be ok!

#import top guidelines
top=read.csv("C:/Users/ipatarc/Desktop/PR_TOP/top-factor -v33.csv")


#program information why journal names have higher overlap than ISSN


library(dplyr)
library(stringr)


#since top guidelines df has ISSN number issued with xxxx-xxxx but SCOPUS has no - xxxxxxxx, than I will remove - sign and do overlap as such

top$Issn_A = str_replace(top$Issn,"-","")

#overlap based on Journal's names (other option possible Issn)
Overlap_SCOPUS_TOP_names <- merge(top,
                                  SCOPUS_content,
                                  by.x = "ï..Journal",
                                  by.y = "Source.Title..Medline.sourced.journals.are.indicated.in.Green.",
                                  all=F)


Overlap_SCOPUS_TOP_PISSN <- merge(top,
                                  SCOPUS_content[SCOPUS_content$Print.ISSN!="",],
                                  by.x = "Issn_A",
                                  by.y = "Print.ISSN")#"E.ISSN" or "Print.ISSN"



SCOPUS_content_EISS=SCOPUS_content[SCOPUS_content$E.ISSN!="",]
head(SCOPUS_content_EISS)
Overlap_SCOPUS_TOP_EISSN <- merge(top,
                                  SCOPUS_content_EISS,
                                  by.x = "Issn_A",
                                  by.y = "E.ISSN")#"E.ISSN" or "Print.ISSN"

#sum(unique(Overlap_SCOPUS_TOP_names$Journal)%in%unique(c(Overlap_SCOPUS_TOP_EISSN$Journal,Overlap_SCOPUS_TOP_PISSN$Journal))
    #overlap 318
 #   sum(unique(Overlap_SCOPUS_TOP_names$Journal)%in%unique(Overlap_SCOPUS_TOP_PISSN$Journal))
    #overlap 273
    
    #which journals do i get on top of overlap with the names
  #  sum(unique(Overlap_SCOPUS_TOP_names$Journal)%in%unique(c(Overlap_SCOPUS_TOP_EISSN$Journal,Overlap_SCOPUS_TOP_PISSN$Journal)))
   # sum(!unique(c(Overlap_SCOPUS_TOP_EISSN$Journal,Overlap_SCOPUS_TOP_PISSN$Journal))%in%unique(Overlap_SCOPUS_TOP_names$Journal))
    
    
    #in nutshell I have 1219 journals added into the top guidelines
    #879 I can overlap with SCOPUS Journals based on their name
    #additional 185 I can overlap with scopus Journals based on PISSN or EISSN
    #this is 1064 journals or 87.3%
    
    
    ############################
    # OVERLAPING scopus and TOP guidelines
    # since cols are different, I need to remove first column for top guidelines
    # idea: merge data frames and then select unique entries
    #Overlap_SCOPUS_TOP_names=Overlap_SCOPUS_TOP_names[,-1]
    
    colnames(Overlap_SCOPUS_TOP_names)[order(colnames(Overlap_SCOPUS_TOP_names))]==colnames(Overlap_SCOPUS_TOP_EISSN)[order(colnames(Overlap_SCOPUS_TOP_EISSN))]
    
    
     Columns_of_Interest=c("ï..Journal","Publisher","Societies","Author.guideline.url",
                          "Data.citation.score","Data.citation.justification","Data.transparency.score",
                          "Data.transparency.justification","Analysis.code.transparency.score",
                          "Analysis.code.transparency.justification","Materials.transparency.score",
                          "Materials.transparency.justification","Design...analysis.reporting.guidelines.score",
                          "Design...analysis.reporting.guidelines.justification","Study.preregistration.score",
                          "Study.preregistration.justification","Analysis.plan.preregistration.score",
                          "Analysis.plan.preregistration.justification","Replication.score","Replication.justification",
                          "Registered.reports...publication.bias.score","Registered.reports...publication.bias.justification",
                          "Open.science.badges.score","Open.science.badges.justification","Total", "Active.or.Inactive",
                          "Coverage","Titles.discontinued.by.Scopus.due.to.quality.issues",
                          "Article.language.in.source..three.letter.ISO.language.codes.","X2021.CiteScore",
                          "Medline.sourced.Title...see.more.info.under.separate.tab.","Open.Access.status",
                          "Articles.in.Press.included.","Added.to.list.May.2022","Source.Type","Title.history.indication",
                          "Related.title.to.title.history.indication","Publisher.s.Name",
                          "Publisher.imprints.grouped.to.main.Publisher","All.Science.Journal.Classification.Codes..ASJC.", 
                          "Top.level...Life.Sciences","Top.level...Social.Sciences","Top.level...Physical.Sciences",
                          "Top.level...Health.Sciences")

    
    
    Overlap_SCOPUS_TOP=unique(rbind(Overlap_SCOPUS_TOP_names[,Columns_of_Interest],
                             Overlap_SCOPUS_TOP_EISSN[,Columns_of_Interest],
                             Overlap_SCOPUS_TOP_PISSN[,Columns_of_Interest]))
    
 
    saveRDS(Overlap_SCOPUS_TOP,"C:/Users/ipatarc/Desktop/PR_TOP/Overlap_SCOPUS_TOP.rds")
    
    
    #1824/2000
    