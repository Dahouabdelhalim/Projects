# This script reads a CSV file in GNU R.
# While reading this file, comments will be created for all variables.
# The comments for values will be stored as attributes (attr) as well.

#ds_file = file.choose()
 setwd("./")
 ds_file = "rdata_anon_survey2020.csv"

ds = read.table(
  file=ds_file, encoding="UTF-8",
  header = FALSE, sep = "\\t", quote = "\\"",
  dec = ".", row.names = "CASE",
  col.names = c(
    "CASE","SERIAL","REF","QUESTNNR","MODE","STARTED","DE01","DE04","DE05","DE07",
    "DE08","DE09","DE10","DE11","OP01_01","OP01_02","OP01_03","OP01_04","OP01_05",
    "OP01_06","OP01_07","OP01_08","OP01_09","OP01_10","OP01_11","OP01_12","OP01_13",
    "OP01_14","OP01_15","OP04_01","TIME001","TIME002","TIME003","TIME_SUM",
    "MAILSENT","LASTDATA","FINISHED","Q_VIEWER","LASTPAGE","MAXPAGE","MISSING",
    "MISSREL","TIME_RSI","DEG_TIME"
  ),
  as.is = TRUE,
  colClasses = c(
    CASE="numeric", SERIAL="character", REF="character", QUESTNNR="character",
    MODE="character", STARTED="POSIXct", DE01="numeric", DE04="numeric",
    DE05="numeric", DE07="numeric", DE08="numeric", DE09="numeric",
    DE10="numeric", DE11="numeric", OP01_01="numeric", OP01_02="numeric",
    OP01_03="numeric", OP01_04="numeric", OP01_05="numeric", OP01_06="numeric",
    OP01_07="numeric", OP01_08="numeric", OP01_09="numeric", OP01_10="numeric",
    OP01_11="numeric", OP01_12="numeric", OP01_13="numeric", OP01_14="numeric",
    OP01_15="numeric", OP04_01="character", TIME001="integer",
    TIME002="integer", TIME003="integer", TIME_SUM="integer",
    MAILSENT="POSIXct", LASTDATA="POSIXct", FINISHED="logical",
    Q_VIEWER="logical", LASTPAGE="numeric", MAXPAGE="numeric",
    MISSING="numeric", MISSREL="numeric", TIME_RSI="numeric", DEG_TIME="numeric"
  ),
  skip = 1,
  check.names = TRUE, fill = TRUE,
  strip.white = FALSE, blank.lines.skip = TRUE,
  comment.char = "",
  na.strings = ""
)

rm(ds_file)

attr(ds, "project") = "prereqs2020"
attr(ds, "description") = "WEBDB2020"
attr(ds, "date") = "2020-04-17 10:28:15"
attr(ds, "server") = "https://www.soscisurvey.de"

# Variable und Value Labels
ds$DE01 = factor(ds$DE01, levels=c("1","2","3","4","-9"), labels=c("Female","Male","Other","Do not disclose","[NA] Not answered"), ordered=FALSE)
ds$DE04 = factor(ds$DE04, levels=c("1","2","3","4","-1","-9"), labels=c("Less than 25%","Between 25% and 50%","Between 50% and 75%","Between 75% and 100%","[NA] I am not taking this course","[NA] Not answered"), ordered=FALSE)
ds$DE05 = factor(ds$DE05, levels=c("1","2","3","4","-1","-9"), labels=c("Less than 25%","Between 25% and 50%","Between 50% and 75%","Between 75% and 100%","[NA] I am not taking this course","[NA] Not answered"), ordered=FALSE)
ds$DE07 = factor(ds$DE07, levels=c("1","2","3","4","5","-1","-9"), labels=c("Very good","Good","Average","Bad","Very bad","[NA] Do not disclose","[NA] Not answered"), ordered=FALSE)
ds$DE08 = factor(ds$DE08, levels=c("1","2","3","4","5","-1","-9"), labels=c("Very little","Little","Average","Much","Very much","[NA] Do not disclose","[NA] Not answered"), ordered=FALSE)
ds$DE09 = factor(ds$DE09, levels=c("1","2","3","4","-1","-9"), labels=c("20 or younger","21-25","26-30","older than 30","[NA] Do not disclose","[NA] Not answered"), ordered=FALSE)
ds$DE10 = factor(ds$DE10, levels=c("1","2","-1","-9"), labels=c("Yes","No","[NA] Do not disclose","[NA] Not answered"), ordered=FALSE)
ds$DE11 = factor(ds$DE11, levels=c("1","2","-1","-9"), labels=c("Yes","No","[NA] Do not disclose","[NA] Not answered"), ordered=FALSE)
attr(ds$OP01_01,"1") = "strongly disagree [0]"
attr(ds$OP01_01,"2") = "[1]"
attr(ds$OP01_01,"3") = "[2]"
attr(ds$OP01_01,"4") = "[3]"
attr(ds$OP01_01,"5") = "strongly agree [4]"
attr(ds$OP01_01,"-1") = "don’t know"
attr(ds$OP01_02,"1") = "strongly disagree [0]"
attr(ds$OP01_02,"2") = "[1]"
attr(ds$OP01_02,"3") = "[2]"
attr(ds$OP01_02,"4") = "[3]"
attr(ds$OP01_02,"5") = "strongly agree [4]"
attr(ds$OP01_02,"-1") = "don’t know"
attr(ds$OP01_03,"1") = "strongly disagree [0]"
attr(ds$OP01_03,"2") = "[1]"
attr(ds$OP01_03,"3") = "[2]"
attr(ds$OP01_03,"4") = "[3]"
attr(ds$OP01_03,"5") = "strongly agree [4]"
attr(ds$OP01_03,"-1") = "don’t know"
attr(ds$OP01_04,"1") = "strongly disagree [0]"
attr(ds$OP01_04,"2") = "[1]"
attr(ds$OP01_04,"3") = "[2]"
attr(ds$OP01_04,"4") = "[3]"
attr(ds$OP01_04,"5") = "strongly agree [4]"
attr(ds$OP01_04,"-1") = "don’t know"
attr(ds$OP01_05,"1") = "strongly disagree [0]"
attr(ds$OP01_05,"2") = "[1]"
attr(ds$OP01_05,"3") = "[2]"
attr(ds$OP01_05,"4") = "[3]"
attr(ds$OP01_05,"5") = "strongly agree [4]"
attr(ds$OP01_05,"-1") = "don’t know"
attr(ds$OP01_06,"1") = "strongly disagree [0]"
attr(ds$OP01_06,"2") = "[1]"
attr(ds$OP01_06,"3") = "[2]"
attr(ds$OP01_06,"4") = "[3]"
attr(ds$OP01_06,"5") = "strongly agree [4]"
attr(ds$OP01_06,"-1") = "don’t know"
attr(ds$OP01_07,"1") = "strongly disagree [0]"
attr(ds$OP01_07,"2") = "[1]"
attr(ds$OP01_07,"3") = "[2]"
attr(ds$OP01_07,"4") = "[3]"
attr(ds$OP01_07,"5") = "strongly agree [4]"
attr(ds$OP01_07,"-1") = "don’t know"
attr(ds$OP01_08,"1") = "strongly disagree [0]"
attr(ds$OP01_08,"2") = "[1]"
attr(ds$OP01_08,"3") = "[2]"
attr(ds$OP01_08,"4") = "[3]"
attr(ds$OP01_08,"5") = "strongly agree [4]"
attr(ds$OP01_08,"-1") = "don’t know"
attr(ds$OP01_09,"1") = "strongly disagree [0]"
attr(ds$OP01_09,"2") = "[1]"
attr(ds$OP01_09,"3") = "[2]"
attr(ds$OP01_09,"4") = "[3]"
attr(ds$OP01_09,"5") = "strongly agree [4]"
attr(ds$OP01_09,"-1") = "don’t know"
attr(ds$OP01_10,"1") = "strongly disagree [0]"
attr(ds$OP01_10,"2") = "[1]"
attr(ds$OP01_10,"3") = "[2]"
attr(ds$OP01_10,"4") = "[3]"
attr(ds$OP01_10,"5") = "strongly agree [4]"
attr(ds$OP01_10,"-1") = "don’t know"
attr(ds$OP01_11,"1") = "strongly disagree [0]"
attr(ds$OP01_11,"2") = "[1]"
attr(ds$OP01_11,"3") = "[2]"
attr(ds$OP01_11,"4") = "[3]"
attr(ds$OP01_11,"5") = "strongly agree [4]"
attr(ds$OP01_11,"-1") = "don’t know"
attr(ds$OP01_12,"1") = "strongly disagree [0]"
attr(ds$OP01_12,"2") = "[1]"
attr(ds$OP01_12,"3") = "[2]"
attr(ds$OP01_12,"4") = "[3]"
attr(ds$OP01_12,"5") = "strongly agree [4]"
attr(ds$OP01_12,"-1") = "don’t know"
attr(ds$OP01_13,"1") = "strongly disagree [0]"
attr(ds$OP01_13,"2") = "[1]"
attr(ds$OP01_13,"3") = "[2]"
attr(ds$OP01_13,"4") = "[3]"
attr(ds$OP01_13,"5") = "strongly agree [4]"
attr(ds$OP01_13,"-1") = "don’t know"
attr(ds$OP01_14,"1") = "strongly disagree [0]"
attr(ds$OP01_14,"2") = "[1]"
attr(ds$OP01_14,"3") = "[2]"
attr(ds$OP01_14,"4") = "[3]"
attr(ds$OP01_14,"5") = "strongly agree [4]"
attr(ds$OP01_14,"-1") = "don’t know"
attr(ds$OP01_15,"1") = "strongly disagree [0]"
attr(ds$OP01_15,"2") = "[1]"
attr(ds$OP01_15,"3") = "[2]"
attr(ds$OP01_15,"4") = "[3]"
attr(ds$OP01_15,"5") = "strongly agree [4]"
attr(ds$OP01_15,"-1") = "don’t know"
attr(ds$FINISHED,"F") = "Canceled"
attr(ds$FINISHED,"T") = "Finished"
attr(ds$Q_VIEWER,"F") = "Respondent"
attr(ds$Q_VIEWER,"T") = "Spectator"
comment(ds$SERIAL) = "Serial number (if provided)"
comment(ds$REF) = "Reference (if provided in link)"
comment(ds$QUESTNNR) = "Questionnaire that has been used in the interview"
comment(ds$MODE) = "Interview mode"
comment(ds$STARTED) = "Time the interview has started (Europe/Berlin)"
comment(ds$DE01) = "gender"
comment(ds$DE04) = "knowVeff"
comment(ds$DE05) = "knowGag"
comment(ds$DE07) = "gradesHighSchool"
comment(ds$DE08) = "activityHighSchool"
comment(ds$DE09) = "ageGroup"
comment(ds$DE10) = "previousDegrees"
comment(ds$DE11) = "taekniskoli"
comment(ds$OP01_01) = "feelingsAgree: I regularly feel like I don’t belong in my current studies."
comment(ds$OP01_02) = "feelingsAgree: I feel that I need to be a good role model for students with similar background/characteristics."
comment(ds$OP01_03) = "feelingsAgree: I feel that my grades are according to my own expectations."
comment(ds$OP01_04) = "feelingsAgree: I regularly feel stupid compared to other students in my courses."
comment(ds$OP01_05) = "feelingsAgree: I regularly feel intimidated by other, more knowledgeable students."
comment(ds$OP01_06) = "feelingsAgree: I often do not dare to ask open questions (in class/on Piazza)."
comment(ds$OP01_07) = "feelingsAgree: I regularly feel that my grades don’t represent my actual performance."
comment(ds$OP01_08) = "feelingsAgree: Working on assignments together with knowledgeable students intimidates me."
comment(ds$OP01_09) = "feelingsAgree: Working on assignments together with knowledgeable students benefits my understanding."
comment(ds$OP01_10) = "feelingsAgree: I am working hard for my university courses."
comment(ds$OP01_11) = "feelingsAgree: I often feel that teachers underestimate my skills/knowledge."
comment(ds$OP01_12) = "feelingsAgree: I regularly expect my (written) exam grades to be worse than they are."
comment(ds$OP01_13) = "feelingsAgree: Working on assignments with vocal students gives me comfort."
comment(ds$OP01_14) = "feelingsAgree: I see myself as a vocal student because I regularly ask questions."
comment(ds$OP01_15) = "feelingsAgree: I feel that my assignments are often graded unfairly."
comment(ds$OP04_01) = "clarText: [01]"
comment(ds$TIME001) = "Time spent on page 1"
comment(ds$TIME002) = "Time spent on page 2"
comment(ds$TIME003) = "Time spent on page 3"
comment(ds$TIME_SUM) = "Time spent overall (except outliers)"
comment(ds$MAILSENT) = "Time when the invitation mailing was sent (personally identifiable recipients, only)"
comment(ds$LASTDATA) = "Time when the data was most recently updated"
comment(ds$FINISHED) = "Has the interview been finished (reached last page)?"
comment(ds$Q_VIEWER) = "Did the respondent only view the questionnaire, omitting mandatory questions?"
comment(ds$LASTPAGE) = "Last page that the participant has handled in the questionnaire"
comment(ds$MAXPAGE) = "Hindmost page handled by the participant"
comment(ds$MISSING) = "Missing answers in percent"
comment(ds$MISSREL) = "Missing answers (weighted by relevance)"
comment(ds$TIME_RSI) = "Degradation points for being very fast"
comment(ds$DEG_TIME) = "Degradation points for being very fast"



# Assure that the comments are retained in subsets
as.data.frame.avector = as.data.frame.vector
`[.avector` <- function(x,i,...) {
  r <- NextMethod("[")
  mostattributes(r) <- attributes(x)
  r
}
ds_tmp = data.frame(
  lapply(ds, function(x) {
    structure( x, class = c("avector", class(x) ) )
  } )
)
mostattributes(ds_tmp) = attributes(ds)
ds = ds_tmp
rm(ds_tmp)

