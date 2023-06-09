library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Load data
	dat <- read.csv("MPB_full_clean.csv")

# Abbreviate authors and titles to fit in plot
	df <- data.frame(	no = dat$Number,
						title = substr(x = dat$Title_abbr, start = 1, stop = 10),
						last = substr(x = dat$Last, start = 1, stop = 10),
						year =  substr(x = dat$Year, start = 3, stop = 4),
						comp = dat$Competition,
						antag = dat$Antagonism,
						mut = dat$Mutualism)

# Cumulative sums
	df$comp_cum <- cumsum(df$comp)
	df$antag_cum <- cumsum(df$antag)
	df$mut_cum <- cumsum(df$mut)

# Cumulative w/o Bascompte and Jordano
	bascompte_id <- which(dat$Last == "Bascompte")
	df$mut_cum_noBas <- df$mut_cum
	df$mut_cum_noBas[bascompte_id:nrow(df)] <- df$mut_cum[bascompte_id:nrow(df)] - df$mut[bascompte_id]

# Prep df 
	df_gg <- df %>% pivot_longer(comp_cum:mut_cum_noBas, names_to = "int_type", values_to = "cumsum") %>%
	  group_by(int_type) %>%
	  arrange(no)

# Add max cum sum as label
	df_gg<-as.data.frame(df_gg)
	df_gg$label<-rep("", nrow(df_gg))
	df_gg[df_gg$no==max(df_gg$no), "label"]<-df_gg[df_gg$no==max(df_gg$no), "cumsum"]

# Plot
	colors <- c("antag_cum" = "#b38811", "comp_cum" = "#b35806", "mut_cum" = "#8073ac", "mut_cum_noBas" = "#b2abd2")
	ggplot(df_gg, aes(x=no, y=cumsum, color=as.factor(int_type)))+
	  geom_step()+
	  theme_bw()+
	  labs(x= "Ordered Princeton Population Monographs in Population Biology (issue no.)", 
	       y = "Approx. no. of pages", color = "Interaction type")+
	  scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 50, 60)) +
	  scale_color_manual(values = colors, name = , 
	                     labels = c("antagonism", "competition", "mutualism", "mutualism w/o Bascompte and Jordano"))+
	  theme(axis.text=element_text(size = 11),
	        axis.title=element_text(size = 11),
	        legend.position=c(.25, .65),
	        legend.background=element_rect(fill = "white", colour = "black", size = 0.25))+
	  geom_text_repel(data = df_gg, aes(label = label), nudge_y = 100, nudge_x = 2, na.rm = T, show.legend = F)+
	  ylim(0,7200)
	
	ggsave("fig3.png", units = "cm", height = 11.3, width=17, dpi=600)	
	