# BEDE Network Survey, 2019-2020 ####
# Nate Emery & Kait Farrell

# For survey questions, see SurveyQuestions.md

# Load packages and data, set graphing theme ####
#install.packages('pacman')
pacman::p_load(tidyverse, corrplot, Hmisc, viridis, 
               rnaturalearth, cowplot, svglite)

# Load clean data file
BEDE_data <- read_csv('./data/survey_data_clean.csv') %>% 
  mutate_if(is.character, as.factor)

# ggplot theme
mytheme <- theme(panel.grid.major = element_line(colour = 'gray80'), 
                 panel.grid.minor = element_blank(),  
                 axis.line.x = element_line(colour = "black"), 
                 axis.line.y = element_line(colour = "black"), 
                 axis.text.x=element_text(size=12, colour='black'), 
                 axis.text.y=element_text(size=12, colour='black'), 
                 axis.title.x=element_text(size=14), 
                 axis.title.y=element_text(size=14),
                 strip.text.x = element_text(size=10), 
                 strip.text.y = element_text(size=10),
                 panel.background = element_rect(fill = NA, color = "black"),
                 strip.background = element_rect(fill = NA, color = 'black'),
                 legend.title=element_text(size=12), 
                 legend.text=element_text(size=10),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.major.y = element_blank())

# Correlation among numeric responses ####
BEDE_Corr <- BEDE_data %>% 
    select(Student_Learn_CourseNonDept : Student_Learn_BeyondInst,
           InterestTraining_DataManage : InterestTraining_Reproducibility,
           BigBarrier_LackOfDeptSupport : BigBarrier_StudentInterest, YearOfDegree)

CorPlot1 <- rcorr(as.matrix(BEDE_Corr))

corrplot(CorPlot1$r, method="ellipse", type="upper", p.mat=CorPlot1$P, insig = "blank",
         number.cex=0.8, tl.cex=0.8, tl.col = 'black')

# Summarize count and percent of demographic category responses ####
demographics <- BEDE_data %>% 
  mutate(Decade = ifelse(YearOfDegree < 1970, '1960s', 
                         ifelse(YearOfDegree < 1980, '1970s',
                                ifelse(YearOfDegree < 1990, '1980s', 
                                       ifelse(YearOfDegree < 2000, '1990s',
                                              ifelse(YearOfDegree < 2010, '2000s',
                                                     ifelse(YearOfDegree < 2020, '2010s', '2020s')))))),
         Decade = factor(Decade)) %>% 
  select(-YearOfDegree) %>% 
  pivot_longer(CarnegieClass : Decade, names_to = "Factor", values_to = "Response") %>% 
  select(Factor, Response) %>% na.omit()

Demographic_Summaries <- demographics %>% 
  group_by(Factor, Response) %>% 
  filter(Factor != 'YearOfDegree') %>% 
  dplyr::summarize(Count = n()) %>% 
  mutate(Percent = round((Count/sum(Count)*100),2)) %>% 
  write_csv('./output/Demographic_Summaries.csv')

# Year of degree
DegreeYear <- BEDE_data %>% select(YearOfDegree)

# Print mean, median degree year
DegreeYear_Summary <- DegreeYear %>% 
  dplyr::summarize(Mean = mean(YearOfDegree, na.rm=T), 
                   Median = median(YearOfDegree, na.rm=T))

ggplot(DegreeYear, aes(YearOfDegree)) + mytheme +
  geom_histogram(binwidth=2) +
  scale_x_continuous("Year terminal degree completed", limits=c(1960,2020),
                     breaks=seq(1960,2020,10))+
  scale_y_continuous("Respondents", limits=c(0,12), breaks=seq(0,12,3)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = DegreeYear_Summary$Mean, lty= 2)+
  annotate("text", x = 1999, y = 12, label = "Mean = 2004")
# 2 NA values dropped for respondents who did not indicate degree earned

# Figure 1: Panel figure of demographics ####
panelTheme <- theme(panel.grid.major = element_line(colour = 'gray80'), 
                 panel.grid.minor = element_blank(),  
                 axis.line.x = element_line(colour = "black"), 
                 axis.line.y = element_line(colour = "black"), 
                 axis.text.x=element_text(size=8, colour='black'), 
                 axis.text.y=element_text(size=8, colour='black'), 
                 axis.title.x=element_text(size=10), 
                 axis.title.y=element_text(size=10),
                 panel.background = element_rect(fill = NA, color = "black"),
                 legend.title=element_text(size=8), 
                 legend.text=element_text(size=8),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.major.y = element_blank())

degree <- ggplot(BEDE_data, aes(YearOfDegree)) + panelTheme +
  geom_histogram(binwidth=1, fill = "#238A8DFF") + 
  xlab("Year degree earned")+
  scale_y_continuous("Respondents", limits=c(0,10), breaks=seq(0,10,2))

carnPlot <- BEDE_data %>% 
  mutate(CarnegieClass = factor(CarnegieClass, levels=c('Baccalaureate College',
                                                        "Master’s College or University",
                                                        'Doctoral University', 'I do not know'),
                                labels = c('Baccalaureate \\nCollege',
                                           "Master’s College \\nor University",
                                           'Doctoral \\nUniversity','I do not know')),
         NumStudents= factor(NumStudents, levels=c('< 5,000 students','5,000 - 15,000 students',
                                                   '> 15,000 students', 'I do not know'),
                             labels=c('< 5,000', '5,000-15,000','> 15,000', 'I do not know')))

carnegie <- ggplot(carnPlot, aes(CarnegieClass)) + panelTheme +
  geom_bar(aes(fill=NumStudents)) + 
  scale_y_continuous("Respondents", limits=c(0,50), breaks=seq(0,50,10))+
  xlab("") +
  scale_fill_viridis("Students",discrete=TRUE) +
  theme(panel.grid.major.y=element_line(),
        panel.grid.major.x=element_blank(),
        #legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = c(.01,.99), legend.justification = c(0,.99))

apptPlot <- BEDE_data %>% 
  mutate(ApptType = factor(ApptType,levels = c('Part-time staff', 'Full-time staff', 
                                                 'Visiting/temporary/adjunct faculty', 
                                               'Tenure-track faculty','Tenured faculty', 
                                               'Other (please specify)'),
                             labels = c('Part-time staff', 'Full-time staff', 'Temporary/\\nadjunct faculty', 
                                        'Tenure-track faculty','Tenured faculty', 'Other')),
         DegreeEarned = factor(DegreeEarned, levels = c("B.S. (or equivalent)","Ph.D. (or equivalent)",
                                                        "Professional degree (e.g. M.D.)", "Other (please specify)"),
                               labels= c("B.S. (or equivalent)", "Ph.D. (or equivalent)",
                               "Prof. degree \\n(e.g. M.D.)", "Other")))

appt <- ggplot(apptPlot, aes(ApptType)) + panelTheme +
  geom_bar(aes(fill=DegreeEarned)) + 
  ylab("Respondents") +
  xlab("") +
  scale_fill_viridis("Degree",discrete=TRUE, direction=1) +
  theme(panel.grid.major.y=element_line(),
        panel.grid.major.x=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = c(.01,.99), legend.justification = c(0,.99))


# Visualize participant locations
mapPoints <- read_csv('./data/Locations.csv') %>% na.omit()
world <- ne_countries(scale = "medium", returnclass = "sf")

locs <- ggplot(data = world) + panelTheme +
  geom_sf(col='gray70') +
  labs(x = "Longitude", y = "Latitude") +
  geom_point(data = mapPoints, aes(x = long, y = lat), size = 2, 
             shape=21, fill = "#238A8DFF", alpha=0.7) +
  coord_sf(xlim = c(-130.00, 125.00), ylim = c(10.00, 80.00), 
           expand = F)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

demoPlot <- plot_grid(degree, locs, appt, carnegie, labels="AUTO", align='v',
                      rel_heights = c(1,1.5))

demoPlot
ggsave('./output/vector/Figure1.svg', width = 8, height = 7, units = "in")

#jpeg(filename = "./output/demographics.jpeg", width = 8, height = 7, units = "in", res = 500)
#demoPlot
#dev.off()

# Figures 2, 3: Frequency of Using Data Skills ####
Use <- BEDE_data %>% 
  select(YearOfDegree, CarnegieClass, UseDS_DataManage : UseDS_Reproducibility) %>% 
  mutate(YearOfDegree = as.numeric(YearOfDegree)) %>% 
  pivot_longer(UseDS_DataManage : UseDS_Reproducibility, 
               names_to = 'Factor', values_to = 'Response') %>% 
  filter(Response != "NA") %>% 
  mutate_if(is.character, as.factor) %>% 
  mutate(Factor = factor(Factor, levels=c('UseDS_Reproducibility','UseDS_Modeling','UseDS_Code',
                                          'UseDS_DataManage','UseDS_DataVis','UseDS_DataAnalysis'), 
                         labels=c('Reproducibility','Modeling','Code',
                                  'Data Management','Data Visualization','Data Analysis')),
         SimpleResponse = ifelse(Response %in% c('Daily','Once to twice per week'), "Frequently", 
                           ifelse(Response == "Once to twice per month", "Often", "Rarely")),
         SimpleResponse = factor(SimpleResponse, levels=c('Rarely','Often','Frequently')),
         CarnegieClass = factor(CarnegieClass, levels=c('Baccalaureate College',
                                                        "Master’s College or University",
                                                        'Doctoral University'),
                                labels = c('Baccalaureate \\nCollege',
                                           "Master’s College \\nor University",
                                           'Doctoral \\nUniversity')))

# Calculate median degree year by frequency of use
use_groups_year <- Use %>% 
  group_by(Factor, SimpleResponse) %>% 
  dplyr::summarize(Count = n(), medianYear = median(YearOfDegree, na.rm=T))

Use_percents <- Use %>% na.omit() %>% 
  group_by(CarnegieClass, Factor, SimpleResponse) %>% 
  dplyr::summarize(Count = n()) %>% 
  mutate(Percent = round((Count/sum(Count)*100),2)) %>% ungroup() 

# Stacked bar plot of frequency of use
use_plot <- ggplot(Use_percents, aes(x = Factor, y = Percent, fill = SimpleResponse)) + mytheme +
  geom_bar(stat="identity", position="stack", col='black') +
  guides(fill = guide_legend(reverse=TRUE)) +
  coord_flip() +
  facet_grid(CarnegieClass ~ .)+
  scale_y_continuous(limits=c(0,100.01), breaks=seq(0,100,25)) +
  labs(y= "Percent of respondents", x = "Data science skill", fill = 'Frequency') +
  scale_fill_viridis(discrete = T)

use_plot
ggsave('./output/vector/Figure2.svg', width = 6.5, height = 7, units = "in")

#jpeg(filename = "./output/Use_frequency.jpeg", width = 6.5, height = 7, units = "in", res = 500)
#use_plot
#dev.off()

# Boxplot of frequency of use by degree year
Use_box <- Use %>% 
  mutate(Factor = factor(Factor, levels=c('Data Analysis','Data Visualization','Data Management',
                                          'Code','Modeling','Reproducibility')))

use_boxplot <- ggplot(Use_box, aes(x = SimpleResponse, y = YearOfDegree, fill = SimpleResponse)) + mytheme +
  geom_jitter(width=.1, cex=1, col='gray60') +
  geom_boxplot(size=.75, width=.3, alpha=.5, outlier.shape=NA) + 
  coord_flip() +
  facet_wrap(Factor ~ .) +
  scale_fill_viridis(discrete = T) +
  scale_y_continuous(limits=c(1965,2020), breaks=seq(1970,2020,10)) +
  labs(y= "Year of degree", x = "Frequency of use")+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1), 
        legend.position='none', strip.text.x = element_text(size=12))

use_boxplot
ggsave('./output/vector/Figure3.svg', width = 6, height = 4, units = "in")

jpeg(filename = "./output/Use_boxplot.jpeg", width = 6, 
     height = 4, units = "in", res = 500)
use_boxplot
dev.off()

# Figure 4: Importance of Data Science Skills for Undergrads ####
Importance <- BEDE_data %>% 
  select(ResponseId : Importance_Learn_Reproducibility) %>% 
  pivot_longer(Importance_Learn_DataManage : Importance_Learn_Reproducibility, 
               names_to = 'Factor', values_to = 'Response') %>% 
  mutate_if(is.character, as.factor) %>% 
  group_by(Factor, Response) %>% 
  dplyr::summarize(Count = n()) %>% 
  mutate(Percent = round((Count/sum(Count)*100),2)) %>% ungroup() %>% 
  mutate(Factor = factor(Factor, levels=c('Importance_Learn_Reproducibility',
                                          'Importance_Learn_Modeling','Importance_Learn_Code',
                                          'Importance_Learn_DataManage',
                                         'Importance_Learn_DataVis', 'Importance_Learn_DataAnalysis'), 
                         labels=c('Reproducibility','Modeling','Code','Data Management','Data Visualization','Data Analysis')),
         Response = factor(Response, levels=c("Not at all important", "Slightly important", 
                                              "Moderately important", "Very important",
                                              "Extremely important")))

# Stacked bar plot of perceived importance of different skills
imp_plot <- ggplot(Importance, aes(x = Factor, y = Percent, fill = Response)) + mytheme +
  geom_bar(stat="identity", position="stack", col='black') +
  guides(fill = guide_legend(reverse=TRUE)) +
  coord_flip() +
  scale_y_continuous(limits=c(0,100), breaks=seq(0,100,25)) +
  labs(y= "Percent of respondents", x = "Data science skill", fill = 'Importance') +
  scale_fill_viridis(discrete=T)

imp_plot
ggsave('./output/vector/Figure4.svg', width = 6.5, height = 4, units = "in")

#jpeg(filename = "./output/Importance.jpeg", width = 6.5, height = 4, units = "in", res = 500)
#imp_plot
#dev.off()

# Figure 5: Where are students receiving training ####
Learning <- BEDE_data %>% 
  select(CarnegieClass, Student_Learn_CourseNonDept : Student_Learn_BeyondInst) %>% 
  pivot_longer(Student_Learn_CourseNonDept : Student_Learn_BeyondInst, 
               names_to = 'Factor', values_to = 'Response') %>% na.omit %>%
  mutate_if(is.character, as.factor) %>% 
  group_by(CarnegieClass, Factor, Response) %>% 
  dplyr::summarize(Count = n()) %>% 
  mutate(Percent = round((Count/sum(Count)*100),2)) %>% ungroup() %>% 
  mutate(CarnegieClass = factor(CarnegieClass, levels=c('Baccalaureate College',
                                                        "Master’s College or University",
                                                        'Doctoral University'),
                                labels = c('Baccalaureate \\nCollege',
                                           "Master’s College \\nor University",
                                           'Doctoral \\nUniversity')),
         Factor = factor(Factor, levels=c("Student_Learn_BeyondInst","Student_Learn_NonCourse",
                                          "Student_Learn_CourseNonDept", "Student_Learn_ElecCourse",
                                          "Student_Learn_Req_Course"),
                           labels=c("Outside of institution", "Outside of coursework",
                                    "Course in another department", "Elective course",
                                    "Required course")),
         Response = factor(Response, levels=c('5','4','3','2','1'),
                           labels = c('5- Least likely','4','3','2','1- Most likely'))) %>% 
  filter(CarnegieClass != "I do not know")

# Stacked bar plot of where students learning skills by carnegie classification
learning_plot <- ggplot(Learning, aes(x = Factor, y = Percent, fill = Response)) + mytheme +
  geom_bar(stat="identity", position="stack", col='black') +
  guides(fill = guide_legend(reverse=TRUE)) +
  facet_grid(CarnegieClass ~ .)+
  coord_flip() +
  scale_y_continuous(limits=c(0,100.01), breaks=seq(0,100,25)) +
  labs(y= "Percent of respondents", 
       x = "Where are student majors learning data science?", fill = 'Ranking') +
  scale_fill_viridis(discrete=T, direction=1)

learning_plot
ggsave('./output/vector/Figure5.svg', width = 7.5, height = 6.5, units = "in")

#jpeg(filename = "./output/Learning.jpeg", width = 7.5, height = 6.5, units = "in", res = 500)
#learning_plot
#dev.off()

# Figure 6: Teaching Data Science Skills ####
Teaching <- BEDE_data %>% 
  select(ResponseId, DoYouTeach_DataManage : DoYouTeach_Reproducibility) %>% 
  pivot_longer(DoYouTeach_DataManage : DoYouTeach_Reproducibility, 
               names_to = 'Factor', values_to = 'Response') %>% na.omit %>%
  mutate_if(is.character, as.factor) %>% 
  group_by(Factor, Response) %>% 
  dplyr::summarize(Count = n()) %>% 
  mutate(Percent = round((Count/sum(Count)*100),2)) %>% ungroup() %>% 
  mutate(Factor = factor(Factor, levels=c('DoYouTeach_Reproducibility','DoYouTeach_Modeling',
                                          'DoYouTeach_Code','DoYouTeach_DataManage',
                                          'DoYouTeach_DataVis', 'DoYouTeach_DataAnalysis'), 
                         labels=c('Reproducibility','Modeling','Code','Data Management','Data Visualization','Data Analysis')),
         Response = factor(Response, levels=c("I don't teach or intend to teach this", 
                                              "I want to teach this but don't know how",
                                             "I intend to teach this", "I teach this"),
                           labels=c("Don't teach or intend to","I want to but don't know how", 
                                              "I intend to teach this", "I teach this")))

# Stacked bar plot of teaching of different skills
teach_plot <- ggplot(Teaching, aes(x = Factor, y = Percent, fill = Response)) + mytheme +
  geom_bar(stat="identity", position="stack", col='black') +
  guides(fill = guide_legend(reverse=TRUE)) +
  coord_flip() +
  scale_y_continuous(limits=c(0,100), breaks=seq(0,100,25)) +
  labs(y= "Percent of respondents", x = "Data science skill", fill = 'Teaching') +
  scale_fill_viridis(discrete=T)

teach_plot
ggsave('./output/vector/Figure6.svg', width = 6.5, height = 3, units = "in")

#jpeg(filename = "./output/Teaching.jpeg", width = 6.5, height = 3, units = "in", res = 500)
#teach_plot
#dev.off()

# Figure 7: Teaching status vs. Year of Degree ####
TeachYear <- BEDE_data %>% 
  select(YearOfDegree, DoYouTeach_DataManage : DoYouTeach_Reproducibility) %>% 
  mutate(YearOfDegree = as.numeric(YearOfDegree)) %>% 
  pivot_longer(DoYouTeach_DataManage : DoYouTeach_Reproducibility, 
               names_to = 'Factor', values_to = 'Response') %>% na.omit %>%
  mutate_if(is.character, as.factor) %>% 
  mutate(Factor = factor(Factor, levels=c('DoYouTeach_DataAnalysis','DoYouTeach_DataVis',
                                         'DoYouTeach_DataManage', 'DoYouTeach_Reproducibility',
                                         'DoYouTeach_Modeling','DoYouTeach_Code'), 
                         labels=c('Data Analysis','Data Visualization','Data Management',
                                  'Reproducibility','Modeling','Code')),
         Response = factor(Response, levels=c("I don't teach or intend to teach this", 
                                              "I want to teach this but don't know how",
                                              "I intend to teach this", "I teach this"),
                           labels=c("Don't teach or \\nintend to","I want to but \\ndon't know how", 
                                    "I intend to \\nteach this", "I teach this")))

# Stacked bar plot of teaching intention vs. year of degree
teachyear_plot <- ggplot(TeachYear, aes(x = Response, y = YearOfDegree, fill = Response)) + mytheme +
  geom_jitter(width=.1, cex=1, col='gray60') +
  geom_boxplot(size=.75, width=.3, alpha=.5, outlier.shape=NA) +
  coord_flip() +
  facet_wrap(Factor ~ .) +
  scale_fill_viridis(discrete=T) +
  scale_y_continuous(limits=c(1965,2020), breaks=seq(1970,2020,10)) +
  labs(y= "Year of degree", x = "Teaching status")+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1), 
        legend.position='none', strip.text.x = element_text(size=12),
        axis.text.y = element_text(size=11))

teachyear_plot
ggsave('./output/vector/Figure7.svg', width = 6.5, height = 5, units = "in")

#jpeg(filename = "./output/Teach_by_YOD.jpeg", width = 6.5, height = 5, units = "in", res = 500)
#teachyear_plot
#dev.off()

# Figure 8: Source of teaching materials ####
Source <- BEDE_data %>% 
  select(CarnegieClass, SourceTeach_DataManage : SourceTeach_Reproducibility) %>% 
  pivot_longer(SourceTeach_DataManage : SourceTeach_Reproducibility, 
               names_to = 'Factor', values_to = 'Response') %>% na.omit %>%
  mutate_if(is.character, as.factor) %>% 
  group_by(CarnegieClass, Response) %>% 
  dplyr::summarize(Count = n()) %>% 
  mutate(Percent = round((Count/sum(Count)*100),2)) %>% ungroup() %>% 
  mutate(CarnegieClass = factor(CarnegieClass, levels=c('Doctoral University',
                                                        "Master’s College or University",
                                                        'Baccalaureate College'),
                                labels = c('Doctoral \\nUniversity',
                                                  "Master’s College \\nor University",
                                           'Baccalaureate \\nCollege')),
         Response = factor(Response, levels=c("Materials developed at my institution", 
                                              "Proprietary materials (require paid service/company)",
                                              "Open source online materials", "My own materials"),
                           labels=c("Materials developed at \\nmy institution", "Proprietary materials",
                                    "Open source \\nonline materials", "My own materials"))) %>% 
  filter(CarnegieClass != "I do not know")

# Stacked bar plot of source of materials by carnegie classification
source_plot <- ggplot(Source, aes(x = CarnegieClass, y = Percent, fill = Response)) + mytheme +
  geom_bar(stat="identity", position="stack", col='black') +
  guides(fill = guide_legend(reverse=TRUE)) +
  coord_flip() +
  scale_y_continuous(limits=c(0,100), breaks=seq(0,100,25)) +
  labs(y= "Percent of respondents", x = "Carnegie classification", fill = 'Source') +
  scale_fill_viridis(discrete=T)

source_plot
ggsave('./output/vector/Figure8.svg', width = 6.5, height = 2.5, units = "in")

#jpeg(filename = "./output/Source.jpeg", width = 6.5, height = 2.5, units = "in", res = 500)
#source_plot
#dev.off()

# Figure 9: Instructor interest in receiving training ####
Interest <- BEDE_data %>% 
  select(ResponseId, InterestTraining_DataManage : InterestTraining_Reproducibility) %>% 
  pivot_longer(InterestTraining_DataManage : InterestTraining_Reproducibility, 
               names_to = 'Factor', values_to = 'Response') %>% 
  mutate_if(is.character, as.factor) %>% 
  group_by(Factor, Response) %>% 
  dplyr::summarize(Count = n()) %>% 
  mutate(Percent = round((Count/sum(Count)*100),2)) %>% ungroup() %>% na.omit() %>% 
  mutate(Factor = factor(Factor, levels=c('InterestTraining_Reproducibility','InterestTraining_Modeling',
                                          'InterestTraining_Code','InterestTraining_DataManage',
                                          'InterestTraining_DataVis', 'InterestTraining_DataAnalysis'), 
                         labels=c('Reproducibility','Modeling','Code','Data Management','Data Visualization','Data Analysis')),
         Response = factor(Response, levels=c('6','5','4','3','2','1'),
                           labels=c('6- Lowest interest','5','4','3','2','1- Highest interest')))

# Stacked bar plot of perceived importance of different skills
interest_plot <- ggplot(Interest, aes(x = Factor, y = Percent, fill = Response)) + mytheme +
  geom_bar(stat="identity", position="stack", col='black') +
  coord_flip() +
  scale_y_continuous(limits=c(0,100), breaks=seq(0,100,25)) +
  labs(y= "Percent of respondents", x = "Data science skill", fill = 'Ranked interest \\nin training') +
  scale_fill_viridis(discrete=T)+
  guides(fill = guide_legend(reverse=T))

interest_plot
ggsave('./output/vector/Figure9.svg', width = 7.5, height = 4, units = "in")

#jpeg(filename = "./output/Interest_Training.jpeg", width = 7.5, height = 4, units = "in", res = 500)
#interest_plot
#dev.off()

# Figure 10: Gap between importance and teaching ####
high_importance <- Importance %>% 
  group_by(Factor) %>% 
  mutate(Skill_Sum = sum(Count)) %>% 
  filter(Response %in% c("Very important","Extremely important")) %>% 
  group_by(Factor, Skill_Sum) %>% 
  dplyr::summarize(Count = sum(Count)) %>% 
  mutate(Percent = round((Count/Skill_Sum*100),2),
         High_Imp = Percent) %>% 
  select(Factor, High_Imp)

actual_teaching <- Teaching %>% 
  group_by(Factor) %>% 
  mutate(Skill_Sum = sum(Count)) %>% 
  filter(Response == "I teach this")%>% 
  group_by(Factor, Skill_Sum) %>% 
  dplyr::summarize(Count = sum(Count)) %>% 
  mutate(Percent = round((Count/Skill_Sum*100),2),
         Teach_This = Percent)  %>% 
  select(Factor, Teach_This)

teaching_gap <- left_join(high_importance, actual_teaching)

gap_plot <- ggplot(teaching_gap, aes(x=Factor, y=High_Imp)) + 
  geom_segment(aes(x=Factor, xend=Factor, y=Teach_This, yend=High_Imp), color="black") +
  geom_point(size=5, color="#20A387FF") +
  geom_point(aes(y=Teach_This), size=5, color="#404788FF") +
  coord_flip() + mytheme +
  scale_y_continuous(limits=c(0,100), breaks=seq(0,100,25)) +
  labs(y= "Percent of respondents", x = "Data science skill")+
  annotate("text", x = 1.1, y = 8, label = "Currently \\nteach")+
  annotate("text", x = 1, y = 54, label= "Consider 'very' or \\n'extremely' important")

gap_plot
ggsave('./output/vector/Figure10.svg', width = 6.5, height = 4, units = "in")

#jpeg(filename = "./output/Teaching_Gap.jpeg", width = 6.5, height = 4, units = "in", res = 500)
#gap_plot
#dev.off()