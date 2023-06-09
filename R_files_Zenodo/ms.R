## ----packages,include=FALSE,cache=FALSE---------------------------------------
library(tidyverse)
library(ggtree)
library(pomp)
library(phylopomp)
library(cowplot)
stopifnot(getRversion() >= "4.1")
stopifnot(packageVersion("pomp")>="3.2")
stopifnot(packageVersion("phylopomp")>="0.0.23.1")
theme_set(theme_bw(base_family="serif"))
set.seed(1159254136)


## ----tree1,fig.dim=c(5,5),out.width="70%",results="hide"----------------------
freeze(
  seed=52212439,
  playSIRwS(
    Beta=1,gamma=0.1,psi=0.07,
    S0=20,I0=2,R0=0,
    times=4,
    tree=TRUE,
    ill=TRUE
  )
) -> x
x %>% plot(points=TRUE,diagram=TRUE)
## x %>% plot(points=TRUE,diagram=TRUE) -> pl
## pl[[1]]+
##   geom_tiplab(
##     aes(alpha=nodecol=="black"),
##     hjust=1,offset=0.25
##   )+
##   geom_nodelab(
##     aes(alpha=nodecol=="green")
##   )


## ----bds0,include=FALSE-------------------------------------------------------
freeze(
  seed=522390503,
  playSIRwS(
    Beta=1,gamma=0.1,psi=0.05,S0=8,I0=1,R0=0,
    times=c(
      2.5,2.76,
      4.25,4.44,
      5.3,5.345,
      15,15.1
    ),
    tree=TRUE,
    ill=TRUE
  )
)-> x

## ----bds1,fig.dim=c(4,4.6),out.width="70%",results="hide"---------------------
plot_grid(
  plotlist=c(
    x[c(1,2),] %>% plot(points=TRUE,diagram=TRUE),
    x[c(3,4),] %>% plot(points=TRUE,diagram=TRUE)
  ),
  labels="AUTO",
  ncol=2
)


## ----bds2,fig.dim=c(4,4.6),out.width="70%",results="hide"---------------------
plot_grid(
  plotlist=c(
    x[c(5,6),] %>% plot(points=TRUE,diagram=TRUE),
    x[c(7,8),] %>% plot(points=TRUE,diagram=TRUE)
  ),
  labels="AUTO",
  ncol=2
)


## ----prune1,fig.dim=c(8,3),out.width="100%",results="hide"--------------------
freeze(
  seed=522390503,
  playSIRwS(
    Beta=1,gamma=0.1,psi=0.05,S0=8,I0=1,R0=0,
    times=10,
    tree=TRUE,
    ill=TRUE
  )
) -> x

x %>%
  getInfo() %>%
  getElement("lineages") %>%
  covariate_table(times="time",order="constant") %>%
  lookup(t=seq(0,10,by=0.01)) %>%
  arrange(t) %>%
  ggplot(aes(x=t,y=lineages))+
  geom_step()+
  expand_limits(x=c(0,10))+
  labs(y=quote("\\u2113(t)"))+
  theme_classic() -> inset 

plot_grid(
  plotlist=c(
    A=x %>% getInfo(prune=FALSE) %>%
      plot(points=TRUE,diagram=TRUE),
    B=x %>% getInfo(prune=TRUE,compact=FALSE) %>%
      plot(points=TRUE,diagram=TRUE) %>%
      {
        .[[1]] <- .[[1]]+expand_limits(x=10)
        .
      },
    C=x %>% getInfo(prune=TRUE,compact=TRUE) %>%
      plot(points=TRUE,diagram=FALSE) %>%
      {
        .[[1]] <- .[[1]]+
          expand_limits(x=10)+
          annotation_custom(
            grob=ggplotGrob(inset),
            xmin=-2,xmax=10,
            ymin=-1,ymax=2
          )+
          expand_limits(y=c(-1,3))
        .
      }
  ),
  nrow=1,
  labels="AUTO"
)


## ----vischain,fig.dim=c(8,2.5),out.width="100%",results="hide"----------------
freeze(
  seed=292791323,
  {
    playSIRwS(
      Beta=1,gamma=0.1,psi=0.2,S0=20,I0=1,R0=0,
      times=5
    )-> x1
    playSIRwS(x1,times=5.5) -> x2
    playSIRwS(x2,times=6) -> x3
  }
)

x1 %>% getInfo(prune=TRUE,compact=TRUE) -> x1
x2 %>% getInfo(prune=TRUE,compact=TRUE) -> x2
x3 %>% getInfo(prune=TRUE,compact=TRUE) -> x3

plot_grid(
  plotlist=list(
  {
    x1 %>% getElement("lineages") -> ell
    x1 %>% plot(points=TRUE,ladderize=FALSE) -> pl
    pl[[1]]+
      expand_limits(x=6)+
      lims(y=c(2,NA))+
      geom_vline(xintercept=head(ell$time,-1),alpha=0.2)+
      geom_vline(xintercept=ell$time[c(6,9)],alpha=0.5,linetype=2)
  },
  {
    x2 %>% getElement("lineages") -> ell
    x2 %>% plot(points=TRUE,ladderize=FALSE) -> pl
    pl[[1]]+
      expand_limits(x=6)+
      lims(y=c(2,NA))+
      geom_vline(xintercept=head(ell$time,-1),alpha=0.2)+
      geom_vline(xintercept=ell$time[c(4,11)],alpha=0.5,linetype=2)
  },
  {
    x3 %>% getElement("lineages") -> ell
    x3 %>% plot(points=TRUE,ladderize=FALSE) -> pl
    pl[[1]]+
      expand_limits(x=6)+
      lims(y=c(2,NA))+
      geom_vline(xintercept=head(ell$time,-1),alpha=0.2)+
      geom_vline(xintercept=ell$time[c(6,12)],alpha=0.5,linetype=2)
  }
  ),
  labels="AUTO",
  nrow=1
)


## ----lbdp1--------------------------------------------------------------------
bake(
  file="lbdp1.rds",
  seed=728604304L,
  playLBDP(lambda=1.5,mu=0.8,psi=1,times=5,n0=1)
) -> x

x |>
  getInfo() |>
  getElement("tree") |>
  newick2df(time=5) -> lbdp_dat

bake(
  file="lbdp2.rds",{
    library(doParallel)
    library(doRNG)
    registerDoParallel()
    registerDoRNG(728604304)
    expand_grid(
      rep=1:16,
      lambda=seq(0.7,2.5,by=0.1),
      mu=0.8,
      psi=1,
      times=5,
      n0=1,
      Np=c(1000,10000)
    ) -> params
    
    foreach (
      p=iter(params,"row"),
      .inorder=FALSE,
      .packages=c("magrittr","phylopomp","pomp"),
      .combine=bind_rows
    ) %dopar% {
      p %$% {
        lbdp_dat |>
          lbdp_exact(lambda=lambda,mu=mu,psi=psi)
      } -> ll1
      p %$% {
        lbdp_dat |>
          lbdp_pomp(lambda=lambda,mu=mu,psi=psi,n0=n0) |>
          pfilter(Np=Np) |>
          logLik()
      } -> ll2
      bind_cols(p,exact=ll1,pf=ll2)
    }
  }
)-> params


## ----lbdp1_plot,fig.dim=c(8,2.8),out.width="100%"-----------------------------
x |>
  getInfo(compact=TRUE) |>
  plot(points=TRUE) |>
  getElement(1L) -> pl1

params |>
  pivot_longer(c(exact,pf)) |>
  unite(name,name,Np) |>
  mutate(
    name=if_else(grepl("exact",name),"exact",name),
    name=gsub("pf_","",name)
  ) |>
  group_by(lambda,mu,psi,times,n0,name) |>
  summarize(
    type=c("logLik","logLik_se"),
    value=logmeanexp(value,se=TRUE)
  ) |>
  ungroup() |>
  pivot_wider(names_from=type) |>
  mutate(
    y=logLik,
    ymax=logLik+2*logLik_se,
    ymin=logLik-2*logLik_se
  ) |>
  filter(logLik>max(logLik)-16) |>
  ggplot(aes(x=lambda,group=name,color=name,
             y=y,ymin=ymin,ymax=ymax))+
  geom_errorbar(
    position="dodge"
  )+
  scale_color_manual(
    labels=c(
      exact="exact",
      `1000`="1000 particles",
      `10000`="10000 particles"
    ),
    values=c(
      exact="black",
      `1000`="blue",
      `10000`="red"
    )
  )+
  labs(
    color=character(0),
    y="log likelihood",
    x=expression(lambda)
  )+
  theme(
    legend.position=c(0.5,0.27)
  ) -> pl2

plot_grid(
  A=pl1,
  B=pl2,
  labels="AUTO",
  nrow=1,
  rel_widths=c(1,1)
)


## ----sir1---------------------------------------------------------------------
bake(
  file="sir1.rds",
  seed=328168304L,
  playSIRwS(
    Beta=4,gamma=1,psi=1,
    S0=97,I0=3,R0=0,t0=0,times=40,
    tree=TRUE
  )
) -> x

x |>
  getInfo() |>
  getElement("tree") |>
  newick2df(time=8) -> sir_dat

bake(
  file="sir2.rds",
  dependson="sir_dat",
  {
    library(doParallel)
    library(doRNG)
    registerDoParallel()
    registerDoRNG(728604304)
    expand_grid(
      Beta=seq(2,9,by=0.2),
      gamma=1,
      psi=1,
      S0=97,I0=3,R0=0,
      t0=0,
      rep=1:16,
      Np=10000,
      ) -> params
    
    foreach (
      p=iter(params,"row"),
      .inorder=FALSE,
      .packages=c("magrittr","phylopomp","pomp"),
      .combine=bind_rows
    ) %dopar% {
      p %$% {
        sir_dat |>
          sir_pomp(
            Beta=Beta,gamma=gamma,psi=psi,
            S0=S0,I0=I0,R0=R0,t0=0
          ) |>
          pfilter(Np=Np) |>
          logLik()
      } -> ll
      bind_cols(p,logLik=ll)
    } -> params
  }
)-> params


## ----sir1_plot,fig.dim=c(8,2.8),out.width="100%"------------------------------
x |>
  getInfo(compact=TRUE) |>
  plot(points=TRUE) |>
  getElement(1L) -> pl1

params |>
  group_by(Beta,gamma,psi,S0,I0,R0,Np) |>
  summarize(
    name=c("logLik","logLik_se"),
    value=logmeanexp(logLik,se=TRUE)
  ) |>
  ungroup() |>
  pivot_wider() |>
  mutate(
    y=logLik,
    ymax=logLik+2*logLik_se,
    ymin=logLik-2*logLik_se
  ) |>
  filter(logLik>max(logLik)-16) |>
  ggplot(aes(x=Beta/(S0+I0+R0),group=Np,
    y=y,ymin=ymin,ymax=ymax))+
  geom_errorbar(
    position="dodge"
  )+
  labs(
    color=character(0),
    y="log likelihood",
    x=expression(b)
  )+
  theme(
    legend.position=c(0.5,0.27)
  ) -> pl2

plot_grid(
  A=pl1,
  B=pl2,
  labels="AUTO",
  nrow=1,
  rel_widths=c(1,1)
)


## ----sirs1--------------------------------------------------------------------
bake(
  file="sirs1.rds",
  seed=328168304L,
  playSIRS(
    Beta=4,gamma=2,psi=1,Delta=1,
    S0=97,I0=3,R0=0,t0=0,times=40,
    tree=TRUE
  ) -> x
) -> x

x |>
  getInfo() |>
  getElement("tree") |>
  newick2df(time=40) -> sirs_dat

bake(
  file="sirs2.rds",
  dependson="sirs_dat",
  {
    library(doParallel)
    library(doRNG)
    registerDoParallel()
    registerDoRNG(728604304)
    expand_grid(
      Beta=4,
      gamma=2,
      psi=1,
      Delta=seq(0.5,2,by=0.05),
      S0=97,I0=3,R0=0,
      t0=0,
      rep=1:16,
      Np=5000,
      ) -> params
    
    foreach (
      p=iter(params,"row"),
      .inorder=FALSE,
      .packages=c("magrittr","phylopomp","pomp"),
      .combine=bind_rows
    ) %dopar% {
      p %$% {
        sirs_dat |>
          sirs_pomp(
            Beta=Beta,gamma=gamma,psi=psi,Delta=Delta,
            S0=S0,I0=I0,R0=R0,t0=0
          ) |>
          pfilter(Np=Np) |>
          logLik()
      } -> ll
      bind_cols(p,logLik=ll)
    } -> params
  }
)-> params


## ----sirs1_plot,fig.dim=c(8,2.8),out.width="100%"-----------------------------
x |>
  getInfo(compact=TRUE) |>
  plot(points=FALSE) |>
  getElement(1L) -> pl1

params |>
  group_by(Beta,gamma,psi,Delta,S0,I0,R0,Np) |>
  summarize(
    name=c("logLik","logLik_se"),
    value=logmeanexp(logLik,se=TRUE)
  ) |>
  ungroup() |>
  pivot_wider() |>
  mutate(
    y=logLik,
    ymax=logLik+2*logLik_se,
    ymin=logLik-2*logLik_se
  ) |>
  filter(round(Beta,2)==4) |>
  filter(logLik>max(logLik)-16) |>
  ggplot(aes(x=Delta,group=Np,
             y=y,ymin=ymin,ymax=ymax))+
  geom_errorbar(
    position="dodge"
  )+
  labs(
    color=character(0),
    y="log likelihood",
    x=expression(delta)
  )+
  theme(
    legend.position=c(0.5,0.27)
  ) -> pl2

plot_grid(
  A=pl1,
  B=pl2,
  labels="AUTO",
  nrow=1,
  rel_widths=c(1,1)
)


## ----sessioninfo,include=FALSE------------------------------------------------
sessionInfo()

