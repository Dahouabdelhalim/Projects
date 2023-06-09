library(shiny)
library(ggplot2)
library(smatr)
library(shinyjs)

function(input, output){

  observeEvent(input$resetAll, {
    reset("G")
  })
  
  observeEvent(input$resetAll, {
    reset("sG")
  })
  
  observeEvent(input$resetAll, {
    reset("i1")
  })
  
  observeEvent(input$resetAll, {
    reset("si1")
  })
  
  observeEvent(input$resetAll, {
    reset("i2")
  })
  
  observeEvent(input$resetAll, {
    reset("si2")
  })
  
  observeEvent(input$resetAll, {
    reset("k1")
  })
  
  observeEvent(input$resetAll, {
    reset("sk1")
  })
  
  observeEvent(input$resetAll, {
    reset("k2")
  })
  
  observeEvent(input$resetAll, {
    reset("sk2")
  })
  
  
  pop.scaling.ggplot<-function(G, sG, i1, si1, i2, si2, k1, sk1, k2, sk2){
    #generate points
    t1all<-NA
    t2all<-NA
    for (i in 1:500){
      Gsample<-rnorm(1,G,sG)
      i1sample<-rnorm(1,i1,si1)
      k1sample<-rnorm(1,k1,sk1)
      i2sample<-rnorm(1,i2,si2)
      k2sample<-rnorm(1,k2,sk2)
      
      t1<-(Gsample*k1sample+i1sample)
      t2<-(Gsample*k2sample+i2sample)
      
      t1all<-rbind(t1all, t1)
      t2all<-rbind(t2all, t2)
    }
    t1t2<-cbind(t1all,t2all)
    t1t2<-na.omit(t1t2)
    t1t2<-as.data.frame(t1t2)
    
    #Line fitting
    st1<-sqrt(k1^2*sG^2+G^2*sk1^2+si1^2+sG^2*sk1^2)
    st2<-sqrt(k2^2*sG^2+G^2*sk2^2+si2^2+sG^2*sk2^2)
    cov<-sG^2*k1*k2
    r<-cov/(st1*st2)
    
    OLS<-cov/(st1^2)
    SMA<-st2/st1
    d<-(OLS^2-r^2)/(r^2*OLS)
    MA<-(d+sqrt(d^2+4))/2
    mut2<-G*k2+i2
    mut1<-G*k1+i1
    TRU<-k2/k1
    
    intOLS<-mut2-OLS*mut1
    intMA<-mut2-MA*mut1
    intSMA<-mut2-SMA*mut1
    intTRU<-mut2-(k2/k1)*mut1
    
    OLSpop<-with(t1t2, line.cis(t1t2[,2],t1t2[,1], method = "OLS", intercept = TRUE))
    SMApop<-with(t1t2, line.cis(t1t2[,2],t1t2[,1], method = "SMA", intercept = TRUE))
    MApop<-with(t1t2, line.cis(t1t2[,2],t1t2[,1], method = "MA", intercept = TRUE))
    
    ggplot(t1t2, aes_string(colnames(t1t2)[1], colnames(t1t2)[2])
    ) + geom_point(
    ) + geom_abline(slope=OLS, intercept=intOLS, colour="red"
    ) + geom_abline(slope=OLSpop[2,1], intercept=OLSpop[1,1], colour="red",linetype="3313" 
    ) + geom_abline(slope=MA, intercept=intMA, colour="blue"
    ) + geom_abline(slope=MApop[2,1], intercept=MApop[1,1], colour="blue",linetype="3313"                  
    ) + geom_abline(slope=SMA, intercept=intSMA, colour="orange"
    ) + geom_abline(slope=SMApop[2,1], intercept=SMApop[1,1], colour="orange",linetype="3313" 
    ) + geom_abline(slope=TRU, intercept=intTRU, colour="green")+
      xlab(expression(T[x]))+
      ylab(expression(T[y]))
  }
  
  
  
  OLSfun<-function(G, sG, i1, si1, i2, si2, k1, sk1, k2, sk2){
    st1<-sqrt(k1^2*sG^2+G^2*sk1^2+si1^2+sG^2*sk1^2)
    st2<-sqrt(k2^2*sG^2+G^2*sk2^2+si2^2+sG^2*sk2^2)
    cov<-sG^2*k1*k2
    OLS<-cov/(st1^2)
    return(OLS)
  }
  
  
  SMAfun<-function(G, sG, i1, si1, i2, si2, k1, sk1, k2, sk2){
    st1<-sqrt(k1^2*sG^2+G^2*sk1^2+si1^2+sG^2*sk1^2)
    st2<-sqrt(k2^2*sG^2+G^2*sk2^2+si2^2+sG^2*sk2^2)
    SMA<-st2/st1
    return(SMA)
  }
  
  MAfun<-function(G, sG, i1, si1, i2, si2, k1, sk1, k2, sk2){
    st1<-sqrt(k1^2*sG^2+G^2*sk1^2+si1^2+sG^2*sk1^2)
    st2<-sqrt(k2^2*sG^2+G^2*sk2^2+si2^2+sG^2*sk2^2)
    cov<-sG^2*k1*k2
    OLS<-cov/(st1^2)
    r<-cov/(st1*st2)
    d<-(OLS^2-r^2)/(r^2*OLS)
    MA<-(d+sqrt(d^2+4))/2
    return(MA)
  }  

  trufun<-function(G, sG, i1, si1, i2, si2, k1, sk1, k2, sk2){
    tru<-k2/k1
    return(tru)
  }  
  
  
  slope.against.k1<-function(G, sG, i1, si1, i2, si2, k1, sk1, k2, sk2){
    fun.1<-function(x) MAfun(G, sG, i1, si1, i2, si2, x, sk1, k2, sk2)
    fun.2<-function(x) SMAfun(G, sG, i1, si1, i2, si2, x, sk1, k2, sk2)
    fun.3<-function(x) OLSfun(G, sG, i1, si1, i2, si2, x, sk1, k2, sk2)
    fun.4<-function(x) trufun(G, sG, i1, si1, i2, si2, x, sk1, k2, sk2)
    p<-ggplot(data = data.frame(x = c(0.01, 2)), aes(x))
    p+stat_function(fun=fun.1, color="blue")+
      stat_function(fun=fun.2, color="orange")+
      stat_function(fun=fun.3, color="red")+
      stat_function(fun=fun.4, color="green") +
      geom_vline(xintercept=k2, color="black",linetype="3313") +
      geom_vline(xintercept= k1, color="black")+
      xlab(expression(mu[k[x]]))+
      ylab("slope")
  }  
  
  slope.against.k2<-function(G, sG, i1, si1, i2, si2, k1, sk1, k2, sk2){
    fun.1<-function(x) MAfun(G, sG, i1, si1, i2, si2, k1, sk1, x, sk2)
    fun.2<-function(x) SMAfun(G, sG, i1, si1, i2, si2, k1, sk1, x, sk2)
    fun.3<-function(x) OLSfun(G, sG, i1, si1, i2, si2, k1, sk1, x, sk2)
    fun.4<-function(x) trufun(G, sG, i1, si1, i2, si2, k1, sk1, x, sk2)
    p<-ggplot(data = data.frame(x = c(0.01, 2)), aes(x))
    p+stat_function(fun=fun.1, color="blue")+
    stat_function(fun=fun.2, color="orange")+
    stat_function(fun=fun.3, color="red")+
    stat_function(fun=fun.4, color="green") +
    geom_vline(xintercept=k1, color="white",linetype="3313") +
    geom_vline(xintercept= k2, color="white")+
      xlab(expression(mu[k[y]]))+
    ylab("slope")
  } 
  
  
  slope.against.si2<-function(G, sG, i1, si1, i2, si2, k1, sk1, k2, sk2){
    fun.1<-function(x) MAfun(G, sG, i1, si1, i2, x, k1, sk1, k2, sk2)
    fun.2<-function(x) SMAfun(G, sG, i1, si1, i2, x, k1, sk1, k2, sk2)
    fun.3<-function(x) OLSfun(G, sG, i1, si1, i2, x, k1, sk1, k2, sk2)
    fun.4<-function(x) trufun(G, sG, i1, si1, i2, x, k1, sk1, k2, sk2)
    p<-ggplot(data = data.frame(x = c(0.01, 1)), aes(x))
    p+stat_function(fun=fun.1, color="blue")+
      stat_function(fun=fun.2, color="orange")+
      stat_function(fun=fun.3, color="red")+
      stat_function(fun=fun.4, color="green") +
      geom_vline(xintercept=si1, color="black",linetype="3313") +
      geom_vline(xintercept= si2, color="black")+
      xlab(expression(sigma[i[y]]))+
      ylab("slope")
  }
  
  slope.against.si1<-function(G, sG, i1, si1, i2, si2, k1, sk1, k2, sk2){
    fun.1<-function(x) MAfun(G, sG, i1, x, i2, si2, k1, sk1, k2, sk2)
    fun.2<-function(x) SMAfun(G, sG, i1, x, i2, si2, k1, sk1, k2, sk2)
    fun.3<-function(x) OLSfun(G, sG, i1, x, i2, si2, k1, sk1, k2, sk2)
    fun.4<-function(x) trufun(G, sG, i1, x, i2, si2, k1, sk1, k2, sk2)
    p<-ggplot(data = data.frame(x = c(0.01, 1)), aes(x))
    p+stat_function(fun=fun.1, color="blue")+
      stat_function(fun=fun.2, color="orange")+
      stat_function(fun=fun.3, color="red")+
      stat_function(fun=fun.4, color="green") +
      geom_vline(xintercept=si2, color="black",linetype="3313") +
      geom_vline(xintercept= si1, color="black")+
      xlab(expression(sigma[i[x]]))+
      ylab("slope")
  }
  
  slope.against.sk1<-function(G, sG, i1, si1, i2, si2, k1, sk1, k2, sk2){
    fun.1<-function(x) MAfun(G, sG, i1, si1, i2, si2, k1, x, k2, sk2)
    fun.2<-function(x) SMAfun(G, sG, i1, si1, i2, si2, k1, x, k2, sk2)
    fun.3<-function(x) OLSfun(G, sG, i1, si1, i2, si2, k1, x, k2, sk2)
    fun.4<-function(x) trufun(G, sG, i1, si1, i2, si2, k1, x, k2, sk2)
    p<-ggplot(data = data.frame(x = c(0.01, 1)), aes(x))
    p+stat_function(fun=fun.1, color="blue")+
      stat_function(fun=fun.2, color="orange")+
      stat_function(fun=fun.3, color="red")+
      stat_function(fun=fun.4, color="green") +
      geom_vline(xintercept=sk2, color="black",linetype="3313") +
      geom_vline(xintercept= sk1, color="black")+
      xlab(expression(sigma[k[x]]))+
      ylab("slope")
  }
  
  slope.against.sk2<-function(G, sG, i1, si1, i2, si2, k1, sk1, k2, sk2){
    fun.1<-function(x) MAfun(G, sG, i1, si1, i2, si2, k1, sk1, k2, x)
    fun.2<-function(x) SMAfun(G, sG, i1, si1, i2, si2, k1, sk1, k2, x)
    fun.3<-function(x) OLSfun(G, sG, i1, si1, i2, si2, k1, sk1, k2, x)
    fun.4<-function(x) trufun(G, sG, i1, si1, i2, si2, k1, sk1, k2, x)
    p<-ggplot(data = data.frame(x = c(0.01, 1)), aes(x))
    p+stat_function(fun=fun.1, color="blue")+
      stat_function(fun=fun.2, color="orange")+
      stat_function(fun=fun.3, color="red")+
      stat_function(fun=fun.4, color="green") +
      geom_vline(xintercept=sk1, color="black",linetype="3313") +
      geom_vline(xintercept= sk2, color="black")+
      xlab(expression(sigma[k[y]]))+
      ylab("slope")
  }
  
  slope.against.sG<-function(G, sG, i1, si1, i2, si2, k1, sk1, k2, sk2){
    fun.1<-function(x) MAfun(G, x, i1, si1, i2, si2, k1, sk1, k2, sk2)
    fun.2<-function(x) SMAfun(G, x, i1, si1, i2, si2, k1, sk1, k2, sk2)
    fun.3<-function(x) OLSfun(G, x, i1, si1, i2, si2, k1, sk1, k2, sk2)
    fun.4<-function(x) trufun(G, x, i1, si1, i2, si2, k1, sk1, k2, sk2)
    p<-ggplot(data = data.frame(x = c(0.01, 1)), aes(x))
    p+stat_function(fun=fun.1, color="blue")+
      stat_function(fun=fun.2, color="orange")+
      stat_function(fun=fun.3, color="red")+
      stat_function(fun=fun.4, color="green") +
      xlab(expression(sigma[S]))+
      ylab("slope")
  }  
  
  slope.against.G<-function(G, sG, i1, si1, i2, si2, k1, sk1, k2, sk2){
    fun.1<-function(x) MAfun(x, sG, i1, si1, i2, si2, k1, sk1, k2, sk2)
    fun.2<-function(x) SMAfun(x, sG, i1, si1, i2, si2, k1, sk1, k2, sk2)
    fun.3<-function(x) OLSfun(x, sG, i1, si1, i2, si2, k1, sk1, k2, sk2)
    fun.4<-function(x) trufun(x, sG, i1, si1, i2, si2, k1, sk1, k2, sk2)
    p<-ggplot(data = data.frame(x = c(-1, 1)), aes(x))
    p+stat_function(fun=fun.1, color="blue")+
      stat_function(fun=fun.2, color="orange")+
      stat_function(fun=fun.3, color="red")+
      stat_function(fun=fun.4, color="green") +
      xlab(expression(mu[S]))+
      ylab("slope")
  } 
  
  
  
  #sk1.v.sk2.plot<-function(G, sG, i1, si1, i2, si2, k1, sk1, k2, sk2,min,max){
    #table<-data.frame("sk1" = seq(min,max, by = (max-min)/100), "sk2" = rep(seq(min,max,by=(max-min)/100), each=110))
  sk1.v.sk2.plot<-function(G, sG, i1, si1, i2, si2, k1, sk1, k2, sk2){  
    table<-data.frame("sk1" = seq(0,1, by = 1/100),"sk2" = rep(seq(0,1, by = 1/100), each =110))
    table$OLS<-OLSfun(G, sG, i1, si1, i2, si2, k1,table$sk1,k2,table$sk2)
    table$SMA<-SMAfun(G, sG, i1, si1, i2, si2, k1,table$sk1,k2,table$sk2)
    table$MA<-MAfun(G, sG, i1, si1, i2, si2, k1,table$sk1,k2,table$sk2)
    table$OLSdiff<-abs(table$OLS-(k2/k1))
    table$SMAdiff<-abs(table$SMA-(k2/k1))
    table$MAdiff<-abs(table$MA-(k2/k1))
    SMAhigh<-(table$SMAdiff<table$MAdiff)&(table$SMAdiff<table$OLSdiff)
    table[SMAhigh,"BestFit"] <-"SMA"                                                          
    OLShigh<-(table$OLSdiff<table$MAdiff)&(table$OLSdiff<table$SMAdiff)
    table[OLShigh,"BestFit"] <-"OLS"
    MAhigh<-(table$MAdiff<table$SMAdiff)&(table$MAdiff<table$OLSdiff)
    table[MAhigh,"BestFit"] <-"MA" 
    ggplot(table,aes(sk1,sk2))+
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      theme_bw()+
      theme(axis.line = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())+
      scale_fill_manual(values = c("MA"="blue", "OLS"="red", "SMA"="orange"))+
      geom_raster(aes(fill = BestFit))+
      geom_hline(yintercept=sk2, color="white") +
      geom_vline(xintercept=sk1, color="white")+
      xlab(expression(sigma[k[x]]))+
      ylab(expression(sigma[k[y]]))
  }
  
  #si1.v.si2.plot<-function(G, sG, i1, si1, i2, si2, k1, sk1, k2, sk2,min,max){
    #table<-data.frame("si1" = seq(min,max, by = (max-min)/100), "si2" = rep(seq(min,max,by=(max-min)/100), each=110))
 si1.v.si2.plot<-function(G, sG, i1, si1, i2, si2, k1, sk1, k2, sk2){
    table<-data.frame("si1" = seq(0,1, by = 1/100),"si2" = rep(seq(0,1, by = 1/100), each =110))
    table$OLS<-OLSfun(G, sG, i1, table$si1, i2, table$si2, k1,sk1,k2,sk2)
    table$SMA<-SMAfun(G, sG, i1, table$si1, i2, table$si2, k1,sk1,k2,sk2)
    table$MA<-MAfun(G, sG, i1, table$si1, i2, table$si2, k1,sk1,k2,sk2)
    table$OLSdiff<-abs(table$OLS-(k2/k1))
    table$SMAdiff<-abs(table$SMA-(k2/k1))
    table$MAdiff<-abs(table$MA-(k2/k1))
    SMAhigh<-(table$SMAdiff<table$MAdiff)&(table$SMAdiff<table$OLSdiff)
    table[SMAhigh,"BestFit"] <-"SMA"                                                          
    OLShigh<-(table$OLSdiff<table$MAdiff)&(table$OLSdiff<table$SMAdiff)
    table[OLShigh,"BestFit"] <-"OLS"
    MAhigh<-(table$MAdiff<table$SMAdiff)&(table$MAdiff<table$OLSdiff)
    table[MAhigh,"BestFit"] <-"MA" 
    ggplot(table,aes(si1,si2))+
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      theme_bw()+
      theme(axis.line = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())+
      scale_fill_manual(name="Line\\nFitting\\nMethod", values = c("MA"="blue", "OLS"="red", "SMA"="orange"))+
      geom_raster(aes(fill = BestFit))+
      geom_hline(yintercept=si2, color="white") +
      geom_vline(xintercept=si1, color="white")+
    xlab(expression(sigma[i[x]]))+
      ylab(expression(sigma[i[y]]))
  }
  
  
  
  output$pop.scaling<-renderPlot({pop.scaling.ggplot(input$G,input$sG,input$i1,input$si1,input$i2,input$si2,input$k1,input$sk1,input$k2,input$sk2)})
  output$slope.against.k1<-renderPlot({slope.against.k1(input$G,input$sG,input$i1,input$si1,input$i2,input$si2,input$k1,input$sk1,input$k2,input$sk2)})
  output$slope.against.k2<-renderPlot({slope.against.k2(input$G,input$sG,input$i1,input$si1,input$i2,input$si2,input$k1,input$sk1,input$k2,input$sk2)})
  output$slope.against.si1<-renderPlot({slope.against.si1(input$G,input$sG,input$i1,input$si1,input$i2,input$si2,input$k1,input$sk1,input$k2,input$sk2)})
  output$slope.against.si2<-renderPlot({slope.against.si2(input$G,input$sG,input$i1,input$si1,input$i2,input$si2,input$k1,input$sk1,input$k2,input$sk2)})
  output$slope.against.sk1<-renderPlot({slope.against.sk1(input$G,input$sG,input$i1,input$si1,input$i2,input$si2,input$k1,input$sk1,input$k2,input$sk2)})
  output$slope.against.sk2<-renderPlot({slope.against.sk2(input$G,input$sG,input$i1,input$si1,input$i2,input$si2,input$k1,input$sk1,input$k2,input$sk2)})
  output$slope.against.sG<-renderPlot({slope.against.sG(input$G,input$sG,input$i1,input$si1,input$i2,input$si2,input$k1,input$sk1,input$k2,input$sk2)})
  output$slope.against.G<-renderPlot({slope.against.G(input$G,input$sG,input$i1,input$si1,input$i2,input$si2,input$k1,input$sk1,input$k2,input$sk2)})
  #output$si1.v.si2.plot<-renderPlot({si1.v.si2.plot(input$G,input$sG,input$i1,input$si1,input$i2,input$si2,input$k1,input$sk1,input$k2,input$sk2,input$Range[1],input$Range[2])})
  #output$sk1.v.sk2.plot<-renderPlot({sk1.v.sk2.plot(input$G,input$sG,input$i1,input$si1,input$i2,input$si2,input$k1,input$sk1,input$k2,input$sk2,input$Range[1],input$Range[2])})
  output$si1.v.si2.plot<-renderPlot({si1.v.si2.plot(input$G,input$sG,input$i1,input$si1,input$i2,input$si2,input$k1,input$sk1,input$k2,input$sk2)})
  output$sk1.v.sk2.plot<-renderPlot({sk1.v.sk2.plot(input$G,input$sG,input$i1,input$si1,input$i2,input$si2,input$k1,input$sk1,input$k2,input$sk2)})
}