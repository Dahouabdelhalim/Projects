packages <- c("dplyr","lubridate","xtable","tidyr","qdapTools","gghighlight",
              "hexbin","ggpubr","censReg","VGAM","stargazer","gridExtra","reshape2")
lapply(packages, require, character.only = TRUE)


library(foreign)

dataCore=read.table("MASZK_dataCore.csv", sep=",", header = T)
data <- read.spss("MASZK_2021_05.sav", to.data.frame=TRUE)

d=data[,c("K24", "K25_1_1", "K25_3", "K25_7",
          "K25_8C1", "K25_8C2", "K25_8C3", "K25_8C4", "K25_8C5", "K25_8C6",
          "K26_2C1", "K26_2C2", "K26_2C3", "K26_2C4", "K26_2C5", "K26_2C6",
          "K26_3",
          "K88_2C1", "K88_2C2", "K88_2C3", "K88_2C4", "K88_2C5", "K88_2C6", "K88_2C7", "K88_2C8",
          "K88_3C1", "K88_3C2", "K88_3C3", "K88_3C4", "K88_3C5", "K88_3C6", "K88_3C7")]

names(d)=c("vaccinated", "month", "vaccine_type", "reject",
           "rej_Pfi", "rej_Mod", "rej_Ast", "rej_Spu", "rej_Sin", "rej_Joh",
           "pref_Pfi", "pref_Mod", "pref_Ast", "pref_Spu", "pref_Sin", "pref_Joh",
           "wait",
           "adv_doc", "adv_sci", "adv_fake", "adv_pol", "adv_fam", "adv_fri", "adv_jou", "adv_cel",
           "source_int", "source_soc", "source_print", "source_radio", "source_tv", "source_family", "source_friends")

h=data.frame(data$RecordNo, d, 
             data$K23_5_1, data$K23_5_2, data$K23_5_3, data$K23_5_4, data$K23_5_5, data$K23_5_6,
             data$K23_6_1, data$K23_6_1O, data$K23_6_2, data$K23_6_2O, data$K23_6_3,
             data$K23_6_3O, data$K23_6_4, data$K23_6_4O, data$K23_6_5, data$K23_6_5O, data$K23_6_6, data$K23_6_6O)

names(h)=c("Id",  "vaccinated", "month", "vaccine_type", "reject",
           "rej_Pfi", "rej_Mod", "rej_Ast", "rej_Spu", "rej_Sin", "rej_Joh",
           "pref_Pfi", "pref_Mod", "pref_Ast", "pref_Spu", "pref_Sin", "pref_Joh",
           "wait",
           "adv_doc", "adv_sci", "adv_fake", "adv_pol", "adv_fam", "adv_fri", "adv_jou", "adv_cel",
           "source_int", "source_soc", "source_print", "source_radio", "source_tv", "source_family", "source_friends",
           "ass_pf", "ass_mo", "ass_as", "ass_sp", "ass_si", "ass_ja",
           "ass_ch_pf", "ass_ch_pf_v","ass_ch_mo", "ass_ch_mo_v","ass_ch_as", "ass_ch_as_v",
           "ass_ch_sp", "ass_ch_sp_v","ass_ch_si", "ass_ch_si_v","ass_ch_ja", "ass_ch_ja_v")


fig1Data = dataCore %>% select(c(Id,VaccineType,Vaccine1Date,Age,Group)) %>% drop_na()

fig1Data$reject=0
fig1Data$reject[fig1Data$Group=="Chr=0&Rej=1" | fig1Data$Group=="Chr=1&Rej=1"]=1

table(fig1Data$reject)

h$pf_deny=0
h$pf_deny[h$ass_pf=="Elfogadhatatlan"]=1
h$mo_deny=0
h$mo_deny[h$ass_mo=="Elfogadhatatlan"]=1
h$as_deny=0
h$as_deny[h$ass_as=="Elfogadhatatlan"]=1
h$sp_deny=0
h$sp_deny[h$ass_sp=="Elfogadhatatlan"]=1
h$si_deny=0
h$si_deny[h$ass_si=="Elfogadhatatlan"]=1


h$sex=data$nem
h$birthyear=data$K1
h$age_d="80<"
h$age_d[h$birthyear>1941 & h$birthyear<1960]="60<80"
h$age_d[h$birthyear>1959 & h$birthyear<1980]="40<60"
h$age_d[h$birthyear>1979 & h$birthyear<2000]="20<40"
h$age_d[h$birthyear>1999]="<20"

h$age_d=as.factor(h$age_d)

h$chronic=data$K8 # chronic illness
h$citytype=data$teltip # city type
h$covid=data$K18_4 # COVID positive test
h$covidinfamily=data$K21 # COVID in family
h$ser_covid=data$K23_2 # in hospital with COVID
h$ser_covid_family=data$K23_3 # family member in hospital with COVID
h$schooling=data$K10 # schooling
h$schooling=as.character(h$schooling)

h$chronic_family="Nem" # chronic illness in family
h$chronic_family[data$K17_2_1=="Igen" | data$K17_2_2=="Igen" | data$K17_2_3=="Igen" | data$K17_2_4=="Igen" | data$K17_2_5=="Igen" | data$K17_2_6=="Igen" | data$K17_2_7=="Igen"]="Igen"

# schooling
h$schooling[h$schooling=="Ã‰rettsÃ©gire Ã©pÃ¼lÅ‘ szakkÃ©pesÃ­tÅ‘ bizonyÃ­tvÃ¡ny"]="Közép fokú"
h$schooling[h$schooling=="Ã‰rettsÃ©gi bizonyÃ­tvÃ¡ny szakkÃ©pesÃ­tÃ©s nÃ©lkÃ¼l"] ="Közép fokú"
h$schooling[h$schooling=="Egyetemi (vagy azzal egyenÃ©rtÃ©kÅ±, pl. MA/MSc) oklevÃ©l"]="Felsõ fokú"
h$schooling[h$schooling=="FÅ‘iskolai (vagy azzal egyenÃ©rtÃ©kÅ±, pl. BA/BSc) oklevÃ©l"]="Felsõ fokú"
h$schooling[h$schooling=="FelsÅ‘fokÃº (akkreditÃ¡lt is) szakkÃ©pesÃ­tÅ‘ bizonyÃ­tvÃ¡ny"]="Felsõ fokú"
h$schooling[h$schooling=="Szakiskolai oklevÃ©l, bizonyÃ­tvÃ¡ny"]="Közép fokú"
h$schooling[h$schooling=="SzakmunkÃ¡skÃ©pzÅ‘ iskolai bizonyÃ­tvÃ¡ny"]="Közép fokú"
h$schooling[h$schooling=="Ã‰rettsÃ©gi bizonyÃ­tvÃ¡ny szakkÃ©pesÃ­tÃ©ssel, kÃ©pesÃ­tÅ‘ bizonyÃ­tvÃ¡ny (az Ã©rettsÃ©givel egyÃ¼tt szerzett szakma)"]="Közép fokú"
h$schooling[h$schooling=="Ã\\u0081ltalÃ¡nos (elemi, polgÃ¡ri) iskola 8. osztÃ¡ly, Ã©vfolyam"]="Általános"
h$schooling[h$schooling=="Ã\\u0081ltalÃ¡nos (elemi, polgÃ¡ri) iskola 8. osztÃ¡lynÃ¡l, Ã©vfolyamnÃ¡l alacsonyabb"]="Általános"

# handle NA
h$covid=as.character(h$covid)
h$covid[h$covid==0]="Nem"
h$covid[is.na(h$covid)]="Nem"


h$ser_covid=as.character(h$ser_covid)
h$ser_covid[h$ser_covid=="0"]="Nem"
h$ser_covid[is.na(h$ser_covid)]="Nem"

h$Id=as.numeric(h$Id)

h$Id=as.integer(h$Id)

m=h[,c(1,5,60)]

fig1Data=merge(fig1Data, m, by="Id", all.x=T, all.y=F)

fig1Data$col="blue"
fig1Data$col[fig1Data$chronic=="Van"]="red"

fig1Data$shape=1
fig1Data$shape[fig1Data$reject.x==1]=19

fig1Data$Vaccine1Date_t=as.Date(fig1Data$Vaccine1Date)

fig1Data$lwd=1
fig1Data$lwd[fig1Data$chronic=="Van"]=2

fig1Data$index=1:nrow(fig1Data)

png("fig1_new.png", width=1200, height=600) 
plot(fig1Data$Vaccine1Date_t, fig1Data$Age, col = fig1Data$col, pch = fig1Data$shape, cex=2,
     xlab="Months", ylab="Age", lwd=2)
lines(lowess(fig1Data$Vaccine1Date_t[fig1Data$chronic=="Nincs"], fig1Data$Age[fig1Data$chronic=="Nincs"]), 
      col = "blue", lwd=3, lty=2)
lines(lowess(fig1Data$Vaccine1Date_t[fig1Data$chronic=="Van"], fig1Data$Age[fig1Data$chronic=="Van"]), 
      col = "red", lwd=3, lty=2)
legend("topleft", 
       legend=c("Chronic NO Reject",  "Not Chronic NO Reject", "Chronic Reject", "Not Chronic Reject"),
       text.col="black",
       col=c("red", "blue", "red", "blue"),
       cex=1.2, bty="n",
       fill="white",
       border='white',
       #      title= "",
       #lwd=c(2,1,2,1),
       pch=c(1,1,19,19),
       pt.bg='white')
dev.off()