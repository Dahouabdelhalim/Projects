#--------------------------------------------------------------
# Chapter 2: Data science in biomedicine
#--------------------------------------------------------------
# Install and load the packages

# install.packages("gtrendsR")
# install.packages("ggplot2")
# install.packages("gridExtra")
library(gtrendsR)
library(ggplot2)
library(gridExtra)

#------------------------------------------------------------- 
# Figure 2.2: Google trends for the terms "Data Science" (red), 
# "Big Data" (green), and "Cloud Computing" (blue) for global 
# queries.
#------------------------------------------------------------- 

res <- gtrends(c("Cloud Computing", "Big Data", 
                 "Data Science"),
               time = "2004-01-01 2019-12-31")
plot(res, main = "")

windows()
world <- res$interest_over_time
world$hits[world$hits == "<1"] <- 0
world$hits<- as.numeric(as.character(world$hits))
world$date <- as.Date(world$date)
world$keyword <- factor(world$keyword, 
                        levels = c("Data Science", 
                                   "Big Data", 
                                   "Cloud Computing"))

#---

gg_world <- ggplot(world, aes(x = date, y = hits, 
                              color = keyword)) + 
  geom_line(size = 1) + geom_point() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Date", y = "Interest", color = "Term") +
  theme(legend.position = "none")
plot1 <- print(gg_world + ggtitle("World")) 

#------------------------------------------------------------- 
# Figure 2.3: Google trends for the terms "Data Science" (red), 
# "Big Data" (green), and "Cloud Computing" (blue) for some 
# countries of Europe
#------------------------------------------------------------- 

res_spain <- gtrends(c("Cloud Computing", "Big Data", 
                       "Data Science"), 
                     time = "2004-01-01 2019-12-31", 
                     geo = "ES")
res_germany <- gtrends(c("Cloud Computing", "Big Data", 
                         "Data Science"), 
                       time = "2004-01-01 2019-12-31",
                       geo = "DE")
res_unitedkingdom <- gtrends(c("Cloud Computing", "Big Data", 
                               "Data Science"), 
                             time = "2004-01-01 2019-12-31",
                             geo = "GB")
res_italy <- gtrends(c("Cloud Computing", "Big Data", 
                       "Data Science"), 
                     time = "2004-01-01 2019-12-31", 
                     geo = "IT")
sp <- res_spain$interest_over_time
sp$hits[sp$hits == "<1"] <- 0
sp$hits<- as.numeric(as.character(sp$hits))
de <- res_germany$interest_over_time
uk <- res_unitedkingdom$interest_over_time
it <- res_italy$interest_over_time
sp$date <- as.Date(sp$date)
de$date <- as.Date(de$date)
uk$date <- as.Date(uk$date)
it$date <- as.Date(it$date)

sp$keyword <- factor(sp$keyword, 
                     levels = c("Data Science", 
                                "Big Data", 
                                "Cloud Computing"))
de$keyword <- factor(de$keyword, 
                     levels = c("Data Science", 
                                "Big Data", 
                                "Cloud Computing"))
uk$keyword <- factor(uk$keyword, 
                     levels = c("Data Science", 
                                "Big Data", 
                                "Cloud Computing"))
it$keyword <- factor(it$keyword, 
                     levels = c("Data Science", 
                                "Big Data", 
                                "Cloud Computing"))
#---

gg_sp <- ggplot(sp, aes(x = date, y = hits, 
                        color = keyword)) + 
  geom_line(size = 1) + geom_point() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none") + 
  labs(x = "Date", y = "Interest", color = "Term")

gg_de <- ggplot(de, aes(x = date, y = hits, 
                        color = keyword)) + 
  geom_line(size = 1) + geom_point() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none") +
  labs(x = "Date", y = "Interest", color = "Term")

gg_uk <- ggplot(uk, aes(x = date, y = hits, 
                        color = keyword)) + 
  geom_line(size = 1) + geom_point() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none") +
  labs(x = "Date", y = "Interest", color = "Term")

gg_it <- ggplot(it, aes(x = date, y = hits, 
                        color = keyword)) + 
  geom_line(size = 1) + geom_point() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none") +
  labs(x = "Date", y = "Interest", color = "Term")


plot1 <- gg_sp + ggtitle("Spain") 
plot2 <- gg_de + ggtitle("Germany")
plot3 <- gg_uk + ggtitle("United Kingdom")
plot4 <- gg_it + ggtitle("Italy")
windows()
grid.arrange(plot1, plot2, plot3, plot4, ncol = 2) 

plots <- list(x = plot1, y = plot2, z = plot3, 
              t = plot4)


#------------------------------------------------------------- 
# Figure 2.4: Google trends for the terms "Data Science" (red), 
# "Big Data" (green), and "Cloud Computing" (blue) for United 
# States and some of its states
#------------------------------------------------------------- 

res_usa <- gtrends(c("Cloud Computing", "Big Data", 
                     "Data Science"), 
                   time = "2004-01-01 2019-12-31",  
                   geo = "US")
res_usa_ma <- gtrends(c("Cloud Computing", "Big Data", 
                        "Data Science"), 
                      time = "2004-01-01 2019-12-31", 
                      geo = "US-MA")
res_usa_ca <- gtrends(c("Cloud Computing", "Big Data", 
                        "Data Science"), 
                      time = "2004-01-01 2019-12-31", 
                      geo = "US-CA")
res_usa_wa <- gtrends(c("Cloud Computing", "Big Data", 
                        "Data Science"), 
                      time = "2004-01-01 2019-12-31", 
                      geo = "US-WA")

usa <- res_usa$interest_over_time
usa$hits[usa$hits == "<1"] <- 0
usa$hits<- as.numeric(as.character(usa$hits))
usa_ma <- res_usa_ma$interest_over_time
usa_ca <- res_usa_ca$interest_over_time
usa_wa <- res_usa_wa$interest_over_time
usa$date <- as.Date(usa$date)
usa_ma$date <- as.Date(usa_ma$date)
usa_ca$date <- as.Date(usa_ca$date)
usa_wa$date <- as.Date(usa_wa$date)

usa$keyword <- factor(usa$keyword, 
                      levels = c("Data Science",
                                 "Big Data", 
                                 "Cloud Computing"))
usa_ma$keyword <- factor(usa_ma$keyword, 
                         levels = c("Data Science",
                                    "Big Data", 
                                    "Cloud Computing"))
usa_ca$keyword <- factor(usa_ca$keyword, 
                         levels = c("Data Science", 
                                    "Big Data", 
                                    "Cloud Computing"))
usa_wa$keyword <- factor(usa_wa$keyword, 
                         levels = c("Data Science",
                                    "Big Data", 
                                    "Cloud Computing"))


gg_usa <- ggplot(usa, aes(x = date, y = hits, 
                          color = keyword)) + 
  geom_line(size = 1) + geom_point() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none") +
  labs(x = "Date", y = "Interest", color = "Term")

gg_usa_ma <- ggplot(usa_ma, aes(x = date, y = hits, 
                                color = keyword)) + 
  geom_line(size = 1) + geom_point() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none") + 
  labs(x = "Date", y = "Interest", color = "Term")

gg_usa_ca <- ggplot(usa_ca, aes(x = date, y = hits, 
                                color = keyword)) + 
  geom_line(size = 1) + geom_point() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none") +
  labs(x = "Date", y = "Interest", color = "Term")

gg_usa_wa <- ggplot(usa_wa, aes(x = date, y = hits, 
                                color = keyword)) + 
  geom_line(size = 1) + geom_point() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none") +
  labs(x = "Date", y = "Interest", color = "Term")


plot5 <- gg_usa + ggtitle("United States")
plot6 <- gg_usa_ma + ggtitle("US - Massachusetts")
plot7 <- gg_usa_ca + ggtitle("US - California")
plot8 <- gg_usa_wa + ggtitle("US - Washington")
windows()
grid.arrange(plot5, plot6, plot7, plot8, ncol = 2)

#-------------------------------------------------------------
# Figure 2.5: Google trends for the terms "Data Science" (red), 
# "Big Data" (green), and "Cloud Computing" (blue) in some 
# countries of Asia and in Australia
#------------------------------------------------------------- 

res_china <- gtrends(c("data science", "cloud computing", 
                       "big data"), 
                     time = "2004-01-01 2019-12-31", 
                     geo = "CN")
res_japan <- gtrends(c("data science", "cloud computing", 
                       "big data"), 
                     time = "2004-01-01 2019-12-31",
                     geo = "JP")
res_india <- gtrends(c("data science", "cloud computing", 
                       "big data"), 
                     time = "2004-01-01 2019-12-31",
                     geo = "IN")
res_australia <- gtrends(c("data science", "cloud computing", 
                           "big data"), 
                         time = "2004-01-01 2019-12-31",
                         geo = "AU")

cn <- res_china$interest_over_time
cn$hits[cn$hits == "<1"] <- 0
cn$hits<- as.numeric(as.character(cn$hits))
ja <- res_japan$interest_over_time
india <- res_india$interest_over_time
au <- res_australia$interest_over_time
cn$date <- as.Date(cn$date)
ja$date <- as.Date(ja$date)
india$date <- as.Date(india$date)
au$date <- as.Date(au$date)

cn$keyword <- as.factor(cn$keyword)
ja$keyword <- as.factor(ja$keyword)
india$keyword <- as.factor(india$keyword)
au$keyword <- as.factor(au$keyword)
levels(cn$keyword) <- c("Data Science", "Big Data", 
                        "Cloud Computing")
levels(ja$keyword) <- c("Data Science", "Big Data", 
                        "Cloud Computing")
levels(india$keyword) <- c("Data Science", "Big Data", 
                           "Cloud Computing")
levels(au$keyword) <- c("Data Science", "Big Data", 
                        "Cloud Computing")

gg_cn <- ggplot(cn, aes(x = date, y = hits, 
                        color = keyword)) + 
  geom_line(size = 1) + geom_point() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none") +
  labs(x = "Date", y = "Interest", color = "Term")

gg_ja <- ggplot(ja, aes(x = date, y = hits, 
                        color = keyword)) + 
  geom_line(size = 1) + geom_point() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none") +
  labs(x = "Date", y = "Interest", color = "Term")
  labs(x = "Date", y = "Interest", color = "Term", 
       legend.position = "none")

gg_india <- ggplot(india, aes(x = date, y = hits, 
                              color = keyword)) + 
  geom_line(size = 1) + geom_point() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none") +
  labs(x = "Date", y = "Interest", color = "Term")
  labs(x = "Date", y = "Interest", color = "Term", 
       legend.position = "none")

gg_au <- ggplot(au, aes(x = date, y = hits, 
                        color = keyword)) + 
  geom_line(size = 1) + geom_point() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none") +
  labs(x = "Date", y = "Interest", color = "Term")
  labs(x = "Date", y = "Interest", color = "Term", 
       legend.position = "none")

plot9 <- gg_cn + ggtitle("China")
plot10 <- gg_ja + ggtitle("Japan")
plot11 <- gg_india + ggtitle("India")
plot12 <- gg_au + ggtitle("Australia")
windows()
grid.arrange(plot9, plot10, plot11, plot12, ncol = 2)

