############################################################ #
#
# Main function of the nCov school-testing simulator
#
############################################################ #
nCov.simulator.Testing<-function(n.teachers, 
                                 lambda.b, lambda.w, 
                                 rho.ch, rho.ad, 
                                 alpha.as, 
                                 R.s, 
                                 classes.comp, 
                                 seeds.ch,seeds.ad, 
                                 test.sens.as, test.sens.s, 
                                 test.delay.s, test.delay.as, 
                                 suscep.ch, 
                                 pdiagn.ch, pdiagn.ad, 
                                 t.detection, t.isolation, 
                                 timewindow.closure, threshold.school, 
                                 prop.immune.ch,
                                 prop.immune.ad,
                                 reseeding.ch,reseeding.ad, 
                                 t.stop,
                                 screening.delay,  
                                 threshold.class,
                                 n.test.week,
                                 strategy,
                                 reg.screening.type,
                                 bool.log.contacts,
                                 R.ctcs,
                                 compliance,
                                 nrs,
                                 seeding.time,
                                 variant){
  
  
  ## SETUP POPULATION ----
  
  n.class<-length(classes.comp)
  n<-n.teachers+sum(classes.comp)
  
  # matrix containing information about the state of the individuals
  status.matrix <- data.frame(infected          = rep(0,n), # 1
                              time.of.infection = NA,       # 2
                              infector          = NA,       # 3
                              severity          = -1,       # 4
                              TimeSymptomOnset  = Inf,      # 5
                              TimeToDetection   = NA,       # 6
                              Detected          = 0,        # 7
                              Category          = 0,        # 8
                              NextTest          = NA,       # 9
                              TeacherAssign     = 0,        # 10
                              Compliance        = 1,        # 11
                              IndexCase         = 0,        # 11
                              TimePosTest       = Inf,
                              Nscreening        = 0)      # 12

  group.label<-1:n.class

class_id<-NULL
for (i_id in 1:length(classes.comp)) {
  class_id<-c(class_id,rep(i_id,classes.comp[i_id]))
}

status.matrix$Category[1:sum(classes.comp)]     <-class_id # we enumerate individuals, the one with the same number will belong to the same class
status.matrix$Category[(sum(classes.comp)+1):n] <-0 #teachers have number 0
status.matrix<-assignTeacherToClasses(status.matrix = status.matrix,classes.comp = classes.comp,n.teachers = n.teachers)

#Proportion of immune
status.matrix$infected[sample(1:sum(classes.comp),floor(prop.immune.ch*sum(classes.comp)))] <- -2 #children that are immune
status.matrix$infected[sample((sum(classes.comp)+1):n,floor(prop.immune.ad*n.teachers))] <- -2 #children that are immune

# other population features ----
status.matrix$recovery.vector       <- rep(NA,n) #vector of recovery times
status.matrix$detection.day         <- rep(Inf,n)
status.matrix$infectives            <- rep(0,n) # vector that indicates who is infectious at the current time: 1 infectious 0 non infectious
status.matrix$isolation             <- rep(0,n) # vector that says who is at in quarantine at the current time. 0 means no isolation, 1 isolation
status.matrix$index.contact.within  <- rep(0,n) # vector that selects the individuals that have to propose a new social contact(global) - 1 yes 0 no
status.matrix$index.contact.between <- rep(0,n) # vector that selects the individuals that have to propose a new social contact(global) - 1 yes 0 no
status.matrix$positive.per.class    <- rep(0,n) # this is actually a vector of size n.class (warning!)
status.matrix$testing.day           <- rep(Inf,n) 
status.matrix$isolation.day         <- rep(Inf,n)
status.matrix$back.from.isolation   <- rep(Inf,n)
status.matrix$det.timewindow        <- rep(0,n) #individuals who test positive in a specific time window (to compute threshold)
status.matrix$school.days.losts     <-rep(0,n)

notcomply<-n-floor(n*compliance)
status.matrix$Compliance[sample(1:n,notcomply)]<-0


# log contacts
status.matrix$log.cnt.within.ch     <- ''
status.matrix$log.cnt.within.ad     <- ''
status.matrix$log.cnt.between.ch    <- ''


# transmission ----
# transmission parameters, per individual
status.matrix$id                      <- 1:n
status.matrix$q                       <- NA
status.matrix$total_infectionPeriod   <- NA
status.matrix$contact_rate_within     <- lambda.w
status.matrix$contact_rate_between    <- lambda.b
status.matrix$susceptibility          <- 1

# set age-specific attributes
status.matrix$susceptibility[status.matrix$Category!=0] <- suscep.ch
status.matrix$rho                               <- rho.ad
status.matrix$rho[status.matrix$Category!=0]    <- rho.ch
status.matrix$pdiagn                            <- pdiagn.ad
status.matrix$pdiagn[status.matrix$Category!=0] <- pdiagn.ch

# adjust contact rate for adults (which do not have between class contacts)
status.matrix$contact_rate_between[status.matrix$Category==0]    <- 0
# print(paste('cnt_rate_between adult',unique(status.matrix$contact_rate_between[status.matrix$Category==0] )))



# contacts ----
# info on the proposed time of the next contact (first colum) and the contact individual (second column)
status.matrix$contact.time.within.pr.ctc       <- NA
status.matrix$contact.time.within.pr.infectee  <- NA
status.matrix$contact.time.between.pr.ctc      <- NA
status.matrix$contact.time.between.pr.infectee <- NA

## SETUP EVENTS ----
# setup the structure to keep track of the events
current.time    <- 0
time.events     <<- matrix(NA,1,3)

events<-data.frame(NextCtc        = Inf,
                   Detection      = Inf,
                   Isolation      = Inf, 
                   BackfromIso    = Inf,
                   Recovery       = Inf, 
                   ThresholdCheck = Inf, 
                   NewSeeding     = Inf,
                   TestingDay     = Inf,
                   ScheduleTest   = Inf)

# initialize events ----
events$NewSeeding     <- 0                    # seed infected to start the simulation
events$ThresholdCheck <- timewindow.closure   # first threshold check 

# repetitive testing options
test.window.schedule <- ifelse(n.test.week>0,7/n.test.week,Inf)  

# RepS.A: test all children of all classes at the start of test window
if(!is.na(reg.screening.type) && reg.screening.type == 'A'){
  events$ScheduleTest  <- test.window.schedule
  test.prop.classes.per.school  <- 1
  test.prop.childeren.per.class <- 1
}


# RepS.B: half of the classes at t and the other halve at t+W/2
if(!is.na(reg.screening.type) && reg.screening.type == 'B'){
  events$ScheduleTest  <- test.window.schedule
  test.prop.classes.per.school <- 0.5
  test.prop.childeren.per.class <- 1
}

# RepS.C: half the children of each classes at t and the other halve at t+W/2
if(!is.na(reg.screening.type) && reg.screening.type == 'C'){
  events$ScheduleTest  <- test.window.schedule
  test.prop.classes.per.school  <- 1
  test.prop.childeren.per.class <- 0.5
}



# between-classes contact
between.classes.contacts<-list()
for (i in 1:n){
  between.classes.contacts[[i]]<--1
}
is.between<-FALSE


## RUN SIMULATION ----

# log variables
err                 <- 0
amountOfTest        <- 0

n.detectedpositive <-0
n.schoolclosure <-0
n.classclosure<-0

#  while(sum(status.matrix$infectives, na.rm = TRUE) > 0){ #while there are still infectives
while(current.time < t.stop){ # while current.time is less than t.stop 
  
  # set events ----
  # Phase 1: individuals who can, propose a new social contact(s)
  for (i in which(status.matrix$index.contact.within==1) ){ # for all the individuals that has to propose a within class contact
    status.matrix$contact.time.within.pr.ctc[i]<-rexp(1,status.matrix$contact_rate_within[i])+current.time# I generate the next interarrival time for individual i
    status.matrix$index.contact.within[i]<-0
  }
  
  for (i in which(status.matrix$index.contact.between==1) ){ # for all the individuals that has to propose a between class contact
    status.matrix$contact.time.between.pr.ctc[i]<-rexp(1,status.matrix$contact_rate_between[i])+current.time# I generate the next interarrival time for individual i
    status.matrix$index.contact.between[i]<-0
  }
  
  # aggregate contact times
  contact.time.overall<-c(status.matrix$contact.time.within.pr.ctc, status.matrix$contact.time.between.pr.ctc) #overall contact times
  
  # note: suppressWarnings() for vectors that contain all NA at the start
  events$NextCtc     <- suppressWarnings(min(contact.time.overall, na.rm = T))  #LW: if all 'NA'   ==> the minimum is 'Inf'
  events$Detection   <- min(status.matrix$detection.day)              #LW: if all 'Inf'  ==> the minimum is 'Inf'
  events$Isolation   <- min(status.matrix$isolation.day)              #LW: if all 'Inf'  ==> the minimum is 'Inf'
  events$BackfromIso <- min(status.matrix$back.from.isolation)        #LW: if all 'Inf'  ==> the minimum is 'Inf'
  events$Recovery    <- suppressWarnings(min(status.matrix$recovery.vector, na.rm = T)) # minimum among the recovery times
  events$TestingDay  <- min(status.matrix$testing.day)
  
  # select the next event 
  next.evts<-colnames(events)[min(events)==events]
  if (length(next.evts)>1){
    next.evts<-sample(next.evts,1)
  }
  
  # event: newSeeding ----
  if (next.evts=="NewSeeding"){ #print("NewSeeding")
    current.time<-events$NewSeeding
    events$NewSeeding<-current.time+seeding.time
    
    nSeeds.child <-ifelse(current.time == 0, seeds.ch, reseeding.ch)
    nSeeds.adult <-ifelse(current.time == 0, seeds.ad, reseeding.ad)
    
    temp.suscept.ch<-which(status.matrix$Category!=0 & status.matrix$infected == 0)
    temp.seeds.ch <- sample_cases(id.potential = temp.suscept.ch, num.infections = nSeeds.child)
    
    temp.suscept.ad<-which(status.matrix$Category==0 & status.matrix$infected == 0)
    temp.seeds.ad <- sample_cases(id.potential = temp.suscept.ad, num.infections = nSeeds.adult)
    
    status.matrix  <- register_infected_cases(id.vector=c(temp.seeds.ch,temp.seeds.ad),
                                              current.time,
                                              status.matrix,
                                              id.infector = NA)
    status.matrix$IndexCase[temp.seeds.ch]<-1
    status.matrix$IndexCase[temp.seeds.ad]<-1
    
  } 
  
  # event: ThresholdCheck ---- 
  if (next.evts=="ThresholdCheck"){ #print("ThresholdCheck")
    current.time           <- events$ThresholdCheck
    events$ThresholdCheck  <- events$ThresholdCheck+timewindow.closure
    for (x in 1:n.class){
      temp.students<-which(status.matrix$Category==x)
      temp.teachers<-which(status.matrix$TeacherAssign==x)
      if (sum(status.matrix$det.timewindow[c(temp.students,temp.teachers)])>=threshold.class){
        #if (status.matrix$positive.per.class[x]>threshold.class){
        status.matrix$isolation.day[temp.students]<-current.time
        time.events <<- rbind(time.events,c(current.time,'ClassIsolation',NA))
        n.classclosure<-n.classclosure+1
      }
    }
    if (sum(status.matrix$det.timewindow)>threshold.school){
      status.matrix$isolation.day <- current.time            # set school in isolation
      time.events <<- rbind(time.events,c(current.time,'SchoolIsolation',NA))
      n.schoolclosure<-n.schoolclosure+1
    }
    #status.matrix$positive.per.class<-rep(0,n) ## reset counter to 0
    status.matrix$det.timewindow<-rep(0,n) ## reset counter to 0
    time.events <<- rbind(time.events,c(current.time,'ThresholdCheck',NA))
  }
  
  # event: TestingDay ---- 
  if (next.evts=="TestingDay"){
    current.time <- events$TestingDay
    person.test  <- which(events$TestingDay==status.matrix$testing.day) # can be that more than one are  quarantined at the same time
    status.matrix$testing.day[person.test] <- Inf
    time.events <<- rbind(time.events,c(current.time,'TestingDay',NA))
    for (s in person.test){
      if (status.matrix$Nscreening[s]>0){
        status.matrix$Nscreening[s]<-status.matrix$Nscreening[s]-1
        status.matrix$testing.day[s]<-current.time+t.detection
      }
      if (status.matrix$Compliance[s]==1){
        if (is.positive(individual     = s, 
                        sensitivity.as = test.sens.as, 
                        sensitivity.s  = test.sens.s, 
                        t.detection    = t.detection, 
                        current.time   = current.time, 
                        status.matrix  = status.matrix )){
          status.matrix$detection.day[s]<-Inf
          status.matrix$TimePosTest[s]<-current.time
          if (status.matrix$TimeSymptomOnset[s] <= current.time){
            status.matrix$isolation.day[s]<-current.time+test.delay.s #now set for PCR          
          }else{
            status.matrix$isolation.day[s]<-current.time+test.delay.as #now set for PCR          
          } 
          #  if (status.matrix$Category[s]>0){
          #status.matrix$positive.per.class[status.matrix$Category[s]]<-status.matrix$positive.per.class[status.matrix$Category[s]]+1
          status.matrix$det.timewindow[s]<-1
          #n.detectedpositive<-n.detectedpositive+1
          #  }
          #else{
          #  temp.cl.teach<-unique(status.matrix[which(status.matrix[,10]==s),8])
          #  positive.per.class[temp.cl.teach]<-positive.per.class[temp.cl.teach]+rep(1,length(temp.cl.teach))
          #}
        }
        amountOfTest<-amountOfTest+1
      }
    }
  }
  
  # event: ScheduleTest ----
  if (next.evts=="ScheduleTest"){
    current.time        <- events$ScheduleTest
    events$ScheduleTest <- events$ScheduleTest+test.window.schedule
    time.events <- rbind(time.events,c(current.time,'ScheduleTest',NA))
    ind.test.beg        <- NULL
    half.classes        <- sample(group.label,n.class*test.prop.classes.per.school)
    for (id.class in half.classes) {
      temp.classmates<-c(which(status.matrix$Category==id.class),which(status.matrix$TeacherAssign==id.class))
      ind.test.beg <-c(ind.test.beg, sample(temp.classmates,floor(length(temp.classmates)*test.prop.childeren.per.class)))
      #temp.teach<-which(status.matrix[,8]==0)
      #ind.test.beg<-c(ind.test.beg,temp.teach)
    }
    #ind.test.beg              <- c(ind.test.beg,which(status.matrix$Category==0)) # add teachers
    status.matrix$testing.day[ind.test.beg] <- current.time+0.5
    
    other.half                <- setdiff(1:n,ind.test.beg)
    status.matrix$testing.day[other.half]   <- current.time + (test.window.schedule*0.5)
  }
  
  # event: Detection ----
  if(next.evts=="Detection"){ #print("Detection")
    current.time         <- events$Detection
    person.detected.all  <- which(events$Detection==status.matrix$detection.day) # can be that multiple individuals are quarantined at the same time
    status.matrix$detection.day[person.detected.all]<-Inf #rest
    
    for(person.detected in person.detected.all){
      if (status.matrix$isolation[person.detected]==0 & status.matrix$Detected[person.detected]==0){
        status.matrix$Detected[person.detected] <- 1
        amountOfTest                            <- amountOfTest+1
        if (is.positive(individual     = person.detected, 
                        sensitivity.as = test.sens.as, 
                        sensitivity.s  = test.sens.s, 
                        t.detection    = t.detection, 
                        current.time   = current.time, 
                        status.matrix  = status.matrix )){ ##LW: symptomatic cases are less likely to be detected?
          time.events <<- rbind(time.events,c(current.time,'Detection',person.detected))
          n.detectedpositive<-n.detectedpositive+1
          status.matrix$TimePosTest[person.detected]<-current.time
          # if symptomatic
          #           if (status.matrix$TimeSymptomOnset[person.detected] >= current.time){ 
          if (status.matrix$TimeSymptomOnset[person.detected] <= current.time){ #AT: if symptomatic when taking the test we do a symptomatic test
            status.matrix$isolation.day[person.detected] <- current.time+test.delay.s 
            test.delay.classmates          <- test.delay.s
            #  }else{ # if not (yet) symptomatic
          }else{ # if asymptomatic when taking the test
            status.matrix$isolation.day[person.detected] <- current.time+test.delay.as 
            test.delay.classmates          <- test.delay.as
          }
          
          if (status.matrix$Category[person.detected]>0){ # if child
            class.mates <- c(which(status.matrix$Category[person.detected]==status.matrix$Category))
            #status.matrix$Detected[class.mates] <- 2
            if(strategy == 'PI'){
              status.matrix$isolation.day[class.mates] <- current.time+test.delay.classmates
            }
            if(strategy == 'RS'){
              #status.matrix$testing.day[class.mates] <- current.time+screening.delay
              status.matrix$testing.day[class.mates] <- current.time+test.delay.classmates+screening.delay #AT: day of testing schoolmates: current.time+test.delay.s+screening.delay
             # ifelse(status.matrix$testing.day[class.mates[1]]<Inf,status.matrix$testing.day[class.mates]<-status.matrix$testing.day[class.mates],status.matrix$testing.day[class.mates]<-current.time+test.delay.classmates+screening.delay)
              status.matrix$Nscreening[class.mates]<-nrs
            }
            #LW: class-id in population status.matrix ==>> dangerous (to be continued)
            #status.matrix$positive.per.class[status.matrix$Category[person.detected]]<-status.matrix$positive.per.class[status.matrix$Category[person.detected]]+1
            status.matrix$det.timewindow[person.detected]<-1
          }else{
            temp.cls<-status.matrix$TeacherAssign[person.detected]
            class.mates<-c(which(status.matrix$Category==temp.cls),which(status.matrix$TeacherAssign==temp.cls))
            class.mates<-setdiff(class.mates,person.detected)
            #status.matrix$Detected[class.mates] <- 2
            if(strategy == 'PI'){
              status.matrix$isolation.day[class.mates] <- current.time+test.delay.classmates
            }
            if(strategy == 'RS'){
              #status.matrix$testing.day[class.mates] <- current.time+screening.delay
              status.matrix$testing.day[class.mates] <- current.time+test.delay.classmates+screening.delay #AT: day of testing schoolmates: current.time+test.delay.s+screening.delay
              #ifelse(status.matrix$testing.day[class.mates[1]]<Inf,status.matrix$testing.day[class.mates]<-status.matrix$testing.day[class.mates],status.matrix$testing.day[class.mates]<-current.time+test.delay.classmates+screening.delay)
              status.matrix$Nscreening[class.mates]<-nrs
            }
            status.matrix$det.timewindow[person.detected]<-1
          }
          
        }
      }
    }
  }
  
  # event: NextCtc ----
  if (next.evts=="NextCtc"){ #print("NextCtc")
    current.time<-events$NextCtc
    if (length(which(events$NextCtc==status.matrix$contact.time.within.pr.ctc))>0){ #if it is a within class contact
      
      infector <- which(status.matrix$contact.time.within.pr.ctc ==events$NextCtc)
      status.matrix$index.contact.within[infector]        <- 1  # re-enable to schedule a contact event
      status.matrix$contact.time.within.pr.ctc[infector]  <- NA # re-enable to schedule a contact event
      
      if (status.matrix$Category[infector]==0){ # if adult
        # get all children
        temp.cls<-status.matrix$TeacherAssign[infector]
        infectee.pool <- c(which(status.matrix$Category==temp.cls,which(status.matrix$TeacherAssign==temp.cls))) 
        if(bool.log.contacts) status.matrix$log.cnt.within.ad[infector] <- paste(status.matrix$log.cnt.within.ad[infector], current.time)
      }else{ # else, if child
        # children in the same class + all the teachers
        #infectee.pool <- which(status.matrix$Category==status.matrix$Category[infector] | status.matrix$Category==0) # if too many teachers are assigned, contacts can happen rarely between children
        infectee.pool <- c(which(status.matrix$Category==status.matrix$Category[infector]),which(status.matrix$TeacherAssign==status.matrix$Category[infector]))
        if(bool.log.contacts) status.matrix$log.cnt.within.ch[infector] <- paste(status.matrix$log.cnt.within.ch[infector], current.time)
      }
      
    }else{ # between class
      infector <- which(status.matrix$contact.time.between.pr.ctc == events$NextCtc) 
      status.matrix$index.contact.between[infector]       <- 1   # re-enable to schedule a contact event
      status.matrix$contact.time.between.pr.ctc[infector] <- NA  # re-enable to schedule a contact event
      # children in other classes
      is.between<-TRUE
      infectee.pool  <- which(status.matrix$Category!=status.matrix$Category[infector] & status.matrix$Category!=0)
      #infectee.pool  <- which(status.matrix$Category!=status.matrix$Category[infector]) #including teachers in the between contacts
      if(bool.log.contacts) status.matrix$log.cnt.between.ch[infector] <- paste(status.matrix$log.cnt.between.ch[infector], current.time)
    }
    
    #pick a random individual (make sure it is not the infector)
    infectee <- sample(setdiff(infectee.pool,infector),1)
    if (is.between==TRUE){
      between.classes.contacts[[infector]]<-c(between.classes.contacts[[infector]],infectee)
      is.between<-FALSE
    }
    acc.rate <- nCov.InfMeasure(t=(current.time-status.matrix$time.of.infection[infector]), 
                                lengthI = status.matrix$total_infectionPeriod[infector])*
      status.matrix$q[infector] *
      status.matrix$susceptibility[infectee]  
    
    if (status.matrix$isolation[infectee]==1){acc.rate<-0}
    if (acc.rate>1){err<-err+1}
    if (status.matrix$infected[infectee]==0 & runif(1)<acc.rate){
      
      status.matrix  <- register_infected_cases(id.vector=infectee,
                                                current.time,
                                                status.matrix,
                                                id.infector = infector)
      
      time.events <<- rbind(time.events,c(current.time,'Infector',infector))
    }
  }
  
  # event: Isolation ----
  if (next.evts=="Isolation"){ #print("Isolation")
    current.time         <- events$Isolation
    isolated.individuals <- which(status.matrix$isolation.day==events$Isolation)
    status.matrix$isolation.day[isolated.individuals] <- Inf
    
    for (k in isolated.individuals){
      if (status.matrix$isolation[k]==0){
        status.matrix$isolation[k]                   <- 1
        status.matrix$back.from.isolation[k]         <- current.time+t.isolation
        status.matrix$contact.time.between.pr.ctc[k] <- NA
        status.matrix$contact.time.within.pr.ctc[k]  <- NA
        status.matrix$index.contact.between[k]       <- 0
        status.matrix$index.contact.within[k]        <- 0
        time.events <<- rbind(time.events,c(current.time,'Isolation',k))
        status.matrix$school.days.losts[k]<-status.matrix$school.days.losts[k]+t.isolation
      }
    }
  }
  
  # event: BackfromIso ----
  if (next.evts=="BackfromIso"){ #print("BackfromIso")
    current.time              <- events$BackfromIso
    individuals.back.from.iso <- which(status.matrix$back.from.isolation==events$BackfromIso)
    status.matrix$back.from.isolation[individuals.back.from.iso]<-Inf
    
    for (w in individuals.back.from.iso) {
      
      if (status.matrix$infected[w]==1){
        # if infected and symptomatic => get out of isolation and re-start making contacts
        #if (status.matrix$TimeSymptomOnset[w]<current.time){
        if ((status.matrix$TimeSymptomOnset[w]<(current.time-t.isolation)) | (status.matrix$TimeSymptomOnset[w]>current.time)){ #AT: individual back from Iso showing symptoms before starting isolation, or asymptomatic are sent back to school
          status.matrix$isolation[w]            <- 0
          status.matrix$Detected[w]             <- 0
          status.matrix$index.contact.within[w] <- 1
          if (status.matrix$Category[w]!=0){ status.matrix$index.contact.between[w]<-1 }
          time.events <- rbind(time.events,c(current.time,'BackfromIso',w))
          #  }else{ # if infected but not (yet) symptomatic => remain in isolation
        }else{ # AT: individuals showing symptoms during isolation are kept in isolation for another iso.period
          status.matrix$back.from.isolation[w]<-current.time+t.isolation
          status.matrix$school.days.losts[w]<-status.matrix$school.days.losts[w]+t.isolation
          time.events <- rbind(time.events,c(current.time,'StillInIso',w))
        }
      }else{
        status.matrix$isolation[w]<-0
        status.matrix$Detected[w]<-0
      }
    }
  }
  
  # event: Recovery ----
  if (next.evts=="Recovery"){ #print("Recovery")
    current.time               <- events$Recovery
    recovered                  <- which(status.matrix$recovery.vector==events$Recovery)
    status.matrix$recovery.vector[recovered] <- NA
    
    # set as recovered (= not detectible anymore)
    status.matrix$infected[recovered] <- -1
    status.matrix$Detected[recovered] <- 0
    time.events <<- rbind(time.events,c(current.time,-1,recovered))
    
    # reset all infection related event features
    status.matrix$contact.time.between.pr.ctc[recovered] <- NA
    status.matrix$contact.time.within.pr.ctc[recovered]  <- NA
    status.matrix$index.contact.between[recovered]       <- 0
    status.matrix$index.contact.within[recovered]        <- 0
    status.matrix$infectives[recovered]                  <- 0
    status.matrix$detection.day[recovered]               <- Inf
    status.matrix$isolation.day[recovered]               <- Inf
    status.matrix$isolation[recovered]                   <- 0
  }
}

## GET MODEL OUTPUT ----

#When also the other pathogen is present.
timev.name<-c("time","event","who")
dimnames(time.events)<-list(NULL,timev.name)

return(list(time.events=time.events, status.matrix=status.matrix,err=err, n.tests=amountOfTest, school.days.losts=sum(status.matrix$school.days.losts), n.schoolclosure=n.schoolclosure, n.classclosure=n.classclosure, n.detectedpositive=n.detectedpositive, between.classes.contacts=between.classes.contacts))
}

# function to sample from the given list, with defensive programming checks 
sample_cases <- function(id.potential,num.infections){
  
  if(num.infections==0){
    return(NULL)
  }
  
  if (length(id.potential)==1){
    id.vector <- id.potential
  } else {
    id.vector<-sample(id.potential,min(num.infections,length(id.potential)), replace = FALSE)
  }
  
  return(id.vector)
}

# function to start infection for the individuals of the given id's at current.time
register_infected_cases <- function(id.vector,current.time,status.matrix,id.infector){
  
  for (i_id in id.vector){
    status.matrix$infected[i_id]             <- 1 
    status.matrix$time.of.infection[i_id]    <- current.time
    status.matrix$recovery.vector[i_id]      <- current.time+nCov.SymptmPeriod(mu.IP = mu.IP,variant=variant) # the total length since infection (Exposed+IP) 
    status.matrix$TimeSymptomOnset[i_id]     <- current.time+nCov.onset.sympt(lengthIP = (status.matrix$recovery.vector[i_id]-current.time))
    status.matrix$total_infectionPeriod[i_id]<- status.matrix$recovery.vector[i_id]-current.time
    status.matrix$infectives[i_id]<-1
    
    status.matrix$infector[i_id] <- id.infector
    
    lambda.w <- status.matrix$contact_rate_within[i_id]
    lambda.b <- status.matrix$contact_rate_between[i_id]
    
    # default severity
    status.matrix$severity[i_id]     <- 1
    
    if (runif(1)<status.matrix$rho[i_id]){               # if the indiviual is going to be symptomatic
      status.matrix$q[i_id]  <- R.s/R.ctcs  # A single q parameter for everyone
      if (runif(1)<status.matrix$pdiagn[i_id]){          # if the individual is going to be diagnosed
        status.matrix$detection.day[i_id]          <- status.matrix$TimeSymptomOnset[i_id]  ##LW: timeToDetection?
        status.matrix$severity[i_id] <- 2
      }
    }else{ # not symptomatic
      status.matrix$TimeSymptomOnset[i_id]<-Inf
      status.matrix$q[i_id]<-R.s/R.ctcs*alpha.as # relative infectiousness of asymptomatic
    }
    time.events<<-rbind(time.events,c(current.time,status.matrix$severity[i_id],i_id))
    
    # register contacts 
    status.matrix$index.contact.within[i_id] <- 1  
    if (status.matrix$Category[i_id]!=0){ 
      status.matrix$index.contact.between[i_id]<-1  
    }
    
  } # end for-loop to seed infected
  
  return(status.matrix)
}

#function to assign teachers to classes
assignTeacherToClasses<-function(status.matrix,classes.comp,n.teachers){
  teachers.id<-which(status.matrix$Category==0)
  count<-1
  for (j in teachers.id){
    if (count==(length(classes.comp)+1)){count<-1}
    status.matrix$TeacherAssign[j]<- count
    count<-count+1
  }
  return(status.matrix)
}



