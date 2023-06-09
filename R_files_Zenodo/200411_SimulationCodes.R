## Simulation codes

## packages
{
  library("animation")
  library("jpeg")
  library("tcltk")
}

## Parameters
{
  #### individual parameters (species-specific)
  ## c(excavation prob, waiting prob, shape1, shape2)
  ## shape parameters are for beta distributions (back distance)
  Hetero_para <- c(0.5, 0.4, 0.483, 0.232)
  Reticuli_para <- c(0.26, 0.56, 0.567, 0.155)
  
  #### individual parameters (species-specific) fir P. sim.
  ## c(transport prob, excavation prob, waiting prob, shape1, shape2, )
  ## shape parameters are for beta distributions (back distance)
  Paraneo_para <- c(0.59, 0.06, 0.30, 0.541, 0.765)
  
  ## Arena setting
  N <- 20 # number of individuals
  length_lim <- 10 + 2 # maximum length of tunnels (cells + 2 for calcuration)
  num_lim <- 40 # maximum number of branches
}

## For R. tibialis and H. aureus
{
  ## Main fuctions
  {
    #### plot tunnel branching patterns (fixed branching angle = 45 degree)
    plot_tunnel <- function(Arena, ind_x, ind_y, tunnel_connect, dual){
      par(mfrow=c(1,1), pin=c(4,4))
      if(dual==1){
        par(mfrow=c(1,2), pin=c(5,5))
        image(1:length_lim, 1:num_lim, Arena, col=c("white","grey","red","green","blue","orange","black","purple"),
              xlab="Tunnel length (mm)", ylab="Tunnel number")
      }
      tunnel_length <- apply(Arena[2:length_lim,]!=1, 2, sum)
      minx <- NULL
      for(i in 1:length((1:num_lim)[apply(Arena[2:length_lim,]!=1, 2, sum)>0])){
        minx <- c(minx, min( (0:length_lim)[c(F,Arena[2:length_lim,i]!=1)] ))
      }
      minx[1] <- 0
      
      beginx = beginy = angle = rep(0,num_lim)
      shita <- pi/4
      note <- shita
      for(i in 1:num_lim){
        if(tunnel_length[i]>0){
          if(tunnel_connect[i] == 0){
            angle[i] = 0;
          } else if(tunnel_connect[i]==1){
            angle[i] = note;
            note = note*(-1);
          } else {
            angle[i] = angle[tunnel_connect[i]] + angle[tunnel_connect[i]]/abs(angle[tunnel_connect[i]])*shita
          }
        }
      }
      for(i in 1:num_lim){
        if(tunnel_length[i]>0){
          if(tunnel_connect[i]==0){
            beginx[i] = 0; beginy[i] = 0;
          } else {
            beginx[i] = beginx[tunnel_connect[i]] + (minx[i]-minx[tunnel_connect[i]])*cos(angle[tunnel_connect[i]]);
            beginy[i] = beginy[tunnel_connect[i]] + (minx[i]-minx[tunnel_connect[i]])*sin(angle[tunnel_connect[i]]);
          }
        }
      }
      beginx[beginx>0] <- beginx[beginx>0]-1
      plot(0, type="n", xlim=c(-5,15), ylim=c(-10,10),
           xlab="x (mm)", ylab="y (mm)", las=1)
      
      for(i in 1:num_lim){
        if(tunnel_length[i]>0){
          endx = beginx[i] + tunnel_length[i]*cos(angle[i])
          endy = beginy[i] + tunnel_length[i]*sin(angle[i])
          arrows(beginx[i],beginy[i],endx,endy, length=0)
        }
      }
      
      for(i in 1:N){
        if(ind_y[i]==1){
          points(ind_x[i]-1,0, col=c("grey","red","green","blue","orange","black","purple")[Arena[ind_x[i],ind_y[i]]], pch=19)
        }else{
          posx <- beginx[ind_y[i]] + (ind_x[i]-minx[ind_y[i]])*cos(angle[ind_y[i]])
          posy <- beginy[ind_y[i]] + (ind_x[i]-minx[ind_y[i]])*sin(angle[ind_y[i]])
          points(posx,posy, col=c("grey","red","green","blue","orange","black","purple")[Arena[ind_x[i],ind_y[i]]], pch=19)
        }
      }
      legend(x=-5, y=-3, legend=c("Move forward", "Excavation", "Backing", "Waiting", "loadong"), col=c("red","green","blue","orange", "purple"), pch=19)
    }
  
    #### Run simulations
    RunSimulation <- function(N, exrate, wrate, shape1, shape2, plot_image, dual){
      ## plot_image = 1: with Sys.sleep(0.1)
      ## plot_image = 2: without Sys.sleep() (use for Gif creation)
      
      wrate <- wrate + exrate
      ## Setting
      ind_status <- rep(0,N)    # go ahead (0), back (1), excavating (2), waiting (3)
      ind_x <- rep(1,N)         # individual position x
      ind_y <- rep(1,N)         # individual position tunnel (=y)
      ind_back <- rep(0,N)      # back distance (sample from random variable)
      ind_side <- rep(0,N)      # side digging or not
      ind_skip <- rep(FALSE, N) # skip or not (used for wait and swap calcuration)
      ind_load <- rep(0,N)      # loading (1) or not (0)
      
      # Arena
      # state: 0:empty, 1:sand, 2:going ahead individual, 3: excavating individual, 
      # 4: backing individuals, 5:waiting individuals, 6:not accessible, 7:loading (back)
      Arena <- matrix(1, length_lim, num_lim)
      Arena[1,1] <- 2
      Arena[1,2:num_lim] <- 6
      Arena[1,num_lim:(num_lim-7)] <- 0:7
      Arena[2,1] <- 0
      tunnel_connect <- rep(0,num_lim)
      
      new_tunnel <- 2 # id of next tunnel
      for(j in 1:1000){
        activeN <- sample(1:N, N, replace=F)
        for(i in activeN){
          if(ind_skip[i]){
            ind_skip[i] <- FALSE;
            next;
          }
          
          # in case out of tunnels
          if(ind_x[i] == 1){
            if(Arena[2,1]==0 && (Arena[3,1]!=4 && Arena[3,1]!=7) ){
              ind_x[i] = 2;
              Arena[ind_x[i], ind_y[i]] = 2;
            }
            
            # in case inside tunnel
          } else {
            
            # ind i is going ahead
            if(ind_status[i] == 0){
              # cell infront of the individual is
              # 0: empty (or side is empty)
              if( (Arena[ind_x[i]+1, ind_y[i]] == 0) ||
                  sum(Arena[ind_x[i]+1,]==0 & tunnel_connect == ind_y[i] & Arena[ind_x[i],] == 1)>0){
                Arena[ind_x[i], ind_y[i]] = 0;
                
                # only front is empty
                if(sum(Arena[ind_x[i]+1,]==0 & tunnel_connect == ind_y[i] & Arena[ind_x[i],] == 1) == 0){
                  ind_x[i] = ind_x[i]+1
                # both empty
                } else if (Arena[ind_x[i]+1,ind_y[i]]==0){
                  r <- runif(1,0,1)
                  if(r < 0.5){
                    ind_x[i] = ind_x[i]+1
                  } else {
                    ind_y[i] = (1:num_lim)[Arena[ind_x[i]+1,]==0 & tunnel_connect == ind_y[i] & Arena[ind_x[i],] == 1]
                    ind_x[i] = ind_x[i]+1
                  }
                # only side is empty
                } else {
                  ind_y[i] = (1:num_lim)[Arena[ind_x[i]+1,]==0 & tunnel_connect == ind_y[i] & Arena[ind_x[i],] == 1]
                  ind_x[i] = ind_x[i]+1
                }
                Arena[ind_x[i], ind_y[i]] = 2;
                
              # 1: sand
              } else if(Arena[ind_x[i]+1, ind_y[i]] == 1){
                if(ind_status[i] == 0){
                  ind_status[i] = 2;
                  Arena[ind_x[i], ind_y[i]] = 3;
                }
              
              # 2: going ahead individual
              } else if(Arena[ind_x[i]+1, ind_y[i]] == 2){
                
              # 3: excavating individual
              } else if(Arena[ind_x[i]+1, ind_y[i]] == 3){
                # in case no side tunnel
                if(sum(Arena[ind_x[i]+1,]!=1 & tunnel_connect == ind_y[i] & Arena[ind_x[i],] == 1) == 0){
                  r <- runif(1,0,1)
                  if(r < exrate){
                    ind_status[i] = 2;
                    ind_side[i] = new_tunnel;
                    tunnel_connect[new_tunnel] = ind_y[i];
                    new_tunnel <- new_tunnel + 1
                    Arena[ind_x[i], ind_y[i]] = 3;
                  } else if(r < wrate){
                    ind_status[i] = 3;
                    Arena[ind_x[i], ind_y[i]] = 5;
                  }
                }
                
              # 4: backing individual  
              } else if(Arena[ind_x[i]+1, ind_y[i]] == 4 || Arena[ind_x[i]+1, ind_y[i]] == 7){
                Arena[ind_x[i], ind_y[i]] = 4;
                ind_status[i] = 1;
                ind_back[i] = ceiling(rbeta(1, shape1, shape2)*ind_x[i])
              } 
              
            # ind i is backing
            } else if (ind_status[i] == 1){
              if(ind_back[i] > 0){
                if (Arena[ind_x[i]-1,ind_y[i]] == 1){
                  # back to the end of tunnel
                  ind_side[i] = 1;
                }
                
                # back straight
                if(ind_side[i] == 0){
                  
                  # back is empty
                  if(Arena[ind_x[i]-1, ind_y[i]] == 0 || (ind_x[i]-1) == 1){
                    Arena[ind_x[i], ind_y[i]] = 0;
                    ind_x[i] = ind_x[i]-1
                    if(ind_load[i] == 0){
                      Arena[ind_x[i], ind_y[i]] = 4;
                    } else {
                      Arena[ind_x[i], ind_y[i]] = 7;
                    }
                    ind_back[i] = ind_back[i]-1
                    
                    if(ind_x[i] == 1){
                      # come out of tunnel
                      Arena[ind_x[i], ind_y[i]] = 2;
                      ind_status[i] = 0;
                      ind_back[i] = 0;
                      ind_load[i] = 0;
                    } 
                    
                    # back is waiting ind
                  } else if (Arena[ind_x[i]-1, ind_y[i]] == 5){
                    swap <- (1:N)[(ind_x == ind_x[i]-1) & (ind_y == ind_y[i])]
                    ind_x[i] = ind_x[i]-1
                    ind_x[swap] = ind_x[swap]+1
                    if(ind_load[i] == 0){
                      Arena[ind_x[i], ind_y[i]] = 4;
                    } else {
                      Arena[ind_x[i], ind_y[i]] = 7;
                    } 
                    Arena[ind_x[swap], ind_y[swap]] = 2;
                    ind_back[i] = ind_back[i]-1
                    ind_skip[swap] <- TRUE
                    ind_status[swap] <- 0
                  }
                  
                  # back side
                } else {
                  # back is empty
                  if(Arena[ind_x[i]-1,tunnel_connect[ind_y[i]]]==0){
                    Arena[ind_x[i], ind_y[i]] = 0;
                    ind_y[i] = tunnel_connect[ind_y[i]];
                    ind_x[i] = ind_x[i]-1
                    if(ind_load[i] == 0){
                      Arena[ind_x[i], ind_y[i]] = 4;
                    } else {
                      Arena[ind_x[i], ind_y[i]] = 7;
                    }
                    ind_back[i] = ind_back[i]-1
                    ind_side[i] = 0;
                    
                    # back is waiting ind
                  } else if (Arena[ind_x[i]-1,tunnel_connect[ind_y[i]]]==5){
                    swap <- (1:N)[(ind_x == ind_x[i]-1) & (ind_y == tunnel_connect[ind_y[i]])]
                    ind_y[swap] = ind_y[i]
                    ind_y[i] = tunnel_connect[ind_y[i]];
                    ind_x[i] = ind_x[i]-1
                    ind_x[swap] = ind_x[swap]+1
                    if(ind_load[i] == 0){
                      Arena[ind_x[i], ind_y[i]] = 4;
                    } else {
                      Arena[ind_x[i], ind_y[i]] = 7;
                    }
                    Arena[ind_x[swap], ind_y[swap]] = 2;
                    ind_back[i] = ind_back[i]-1
                    ind_skip[swap] <- TRUE
                    ind_status[swap] <- 0
                    ind_side[i] = 0;
                  }
                }
              } else {
                Arena[ind_x[i], ind_y[i]] = 2;
                ind_status[i] = 0;
                ind_load[i] = 0;
              }
              
              # ind i is excavating
            } else if (ind_status[i] == 2){
              if(ind_side[i] == 0){
                # excavate forward
                if(Arena[ind_x[i]+1, ind_y[i]] == 1){
                  Arena[ind_x[i]+1, ind_y[i]] = 0;
                  Arena[ind_x[i], ind_y[i]] = 7;
                  ind_back[i] = ceiling(rbeta(1, shape1, shape2)*ind_x[i])
                  ind_status[i] = 1;
                  ind_load[i] = 1;
                }
              } else {
                # excavate side
                if(Arena[ind_x[i]+1, ind_side[i]] == 1){
                  Arena[ind_x[i]+1, ind_side[i]] = 0;
                  Arena[ind_x[i], ind_y[i]] = 7;
                  ind_back[i] = ceiling(rbeta(1, shape1, shape2)*ind_x[i])
                  ind_status[i] = 1;
                  ind_side[i] = 0; 
                  ind_load[i] = 1;
                }
              }
              
              # ind i is waiting
            } else if (ind_status[i] == 3){ 
              
            }
          }
        }
        
        if(plot_image == 1){
          par(mfrow=c(1,1), pin=c(5,5))
          plot_tunnel(Arena, ind_x, ind_y, tunnel_connect, dual); Sys.sleep(0.1); 
        } else if (plot_image == 2){
          par(mfrow=c(1,1), pin=c(5,5))
          plot_tunnel(Arena, ind_x, ind_y, tunnel_connect, dual);
        }
        if(sum(Arena[length_lim,]==0) > 0){
          break;
        }
        options(warn=2)
      }
      #plot_tunnel(Arena, ind_x, ind_y, tunnel_connect, dual);
      
      return(list(Arena, ind_x, ind_y, tunnel_connect))
    }
  }
  
  ## Calcuration
  ## for number of branches
  { 
    Res <- NULL
    for(iii in 1:1000){
      iter <- 15
      tipnum1 = tipnum2 = sumlength1 = sumlength2 <- NULL
      #par(mfrow=c(5,10), pin=c(1,1))
      Parameter <- Reticuli_para
      for(replicates in 1:iter){
        res <- RunSimulation(N, Parameter[1], Parameter[2], Parameter[3], Parameter[4],
                             plot_image = 0, dual = 1)
        
        tunnel_length <- apply(res[[1]][2:length_lim,]!=1, 2, sum)
        tipnum1 <- c(tipnum1, length(tunnel_length[tunnel_length>0]))
        tipnum2 <- c(tipnum2, length(tunnel_length[tunnel_length>1]))
        sumlength1 <- c(sumlength1, sum(tunnel_length[tunnel_length>0]))
        sumlength2 <- c(sumlength2, sum(tunnel_length[tunnel_length>1]))
        
        #plot_tunnel(res[[1]], res[[2]], res[[3]], res[[4]], dual=1)
      }
      Res<-c(Res,mean(tipnum1))
    }
  
    hist(Res)
    
    #Res_tip <- Res
    #save(Res_tip, file="Res_tip.Rdata")
    #write.table(R_Res, "E:Reticuli-tip.txt")
    
    #Het_tip <- Res
    #save(Het_tip, file="Het_tip.Rdata")
    #write.table(H_Res, "E:Hetero-tip.txt")
  }
  
  ## for tunnel length by branching segments
  {
    iter <- 1000
    branch1 = branch2 = branch3 <- NULL
    Parameter <- Hetero_para
    for(replicates in 1:iter){
      res <- RunSimulation(N, Parameter[1], Parameter[2], Parameter[3], Parameter[4],
                           plot_image = 0, dual = 1)
      res2 <- res[[1]][3:12,]
      res2[res2>1] <- 0
      res2
      for(i in (1:num_lim)[res[[4]]>0]){
        if(sum(res2[,i]==0)==0){
          res[[4]][i]<-0
        }
      }
      if(sum(res[[4]])==0){
        branch3 <- c(branch3, 10)
      } else if (sum(res[[4]])==1) {
        b1st <- min((1:10)[apply(res2[,2:num_lim]==0, 1, sum)>0]-1)
        branch1 <- c(branch1, b1st)
        b2nd <- apply(res2[(b1st+1):10,]==0, 2, sum)
        branch3 <- c(branch3, b2nd[b2nd>0])
      } else {
        for(i in 1:20){
          if(sum(res2[,i]==0)>0){
            if(sum(res[[4]]==i)>0){
              karimat <- matrix(1:10, 10, sum(res[[4]]==i)) * (res2[,res[[4]]==i] == 0)
              karimat[karimat==0] <- 100
              cut <- apply(karimat, 2, min)-1
              cut <- sort(cut)
              if(length(cut)==1){
                if(i==1){
                  branch1 <- c(branch1, cut)
                } else {
                  branch2 <- c(branch2, cut - (min((1:10)[res2[,i]==0])-1) )
                }
                branch3 <- c(branch3, max((1:10)[res2[,i]==0])-cut)
              } else {
                for(j in 1:length(cut)){
                  if(i == 1 && j == 1){
                    branch1 <- c(branch1, cut[1])
                  } else if (j == 1){
                    branch2 <- c(branch2, cut[1]-(min((1:10)[res2[,i]==0])-1) )
                  } else {
                    branch2 <- c(branch2, cut[j]-max(cut[1:(j-1)]))
                  }
                }
                branch3 <- c(branch3, max((1:10)[res2[,i]==0])-max(cut))
              }
            } else {
              branch3 <- c(branch3, sum(res2[,i]==0))
            }
          }
        }
      }
      #branch1
      #branch2
      #branch3
      #plot_tunnel(res[[1]], res[[2]], res[[3]], res[[4]], dual=1)
    }
    
    branch = c( rep(1,length(branch1)), rep(2,length(branch2)), rep(3,length(branch3)))
    length = c(branch1, branch2, branch3)
    d <- data.frame(branch, length)
    posn.d <- position_dodge(width = 0.5)
    ggplot(d, aes(branch, length)) + 
      stat_summary(geom = 'errorbar', position = posn.d, fun.data = mean_se, width=0.25) +
      stat_summary(geom = 'point', position = posn.d, fun.y = mean, size = 3) +
      theme_bw() + theme(aspect.ratio = 1)
    
    Res_branch <- d
    #save(Res_branch, file="Res_branch.Rdata")
    
    Het_branch <- d
    #save(Het_branch, file="Het_branch.Rdata")
    
    #write.table(d, "E:\\\\Reticlu_branch.txt")
    #write.table(d, "E:\\\\Hetero_branch.txt")
  }
  
  ## to make gif animation
  if(F){
    Parameter <- Reticuli_para
    saveGIF(expr={ RunSimulation(N, Parameter[1], Parameter[2], Parameter[3], Parameter[4],
                   plot_image = 2, dual = 0)
      }, interval = 0.05, movie.name = "E:Reticuli.gif", loop=1)
  }
}

## For P. simplicicornis
{
  ## Main fuction
  {
    plot_tunnel <- function(Arena, ind_x, ind_y, tunnel_connect, dual){
    if(dual==1){
      par(mfrow=c(1,2), pin=c(5,5))
      image(1:length_lim, 1:num_lim, Arena, col=c("white","grey","red","green","blue","orange","black","purple","purple2","brown"),
            xlab="Tunnel length (mm)", ylab="Tunnel number")
    }
    tunnel_length <- apply(Arena[2:length_lim,]!=1, 2, sum)
    minx <- NULL
    for(i in 1:length((1:num_lim)[apply(Arena[2:length_lim,]!=1, 2, sum)>0])){
      minx <- c(minx, min( (0:length_lim)[c(F,Arena[2:length_lim,i]!=1)] ))
    }
    minx[1] <- 0
    
    beginx = beginy = angle = rep(0,num_lim)
    shita <- pi/4
    note <- shita
    for(i in 1:num_lim){
      if(tunnel_length[i]>0){
        if(tunnel_connect[i] == 0){
          angle[i] = 0;
        } else if(tunnel_connect[i]==1){
          angle[i] = note;
          note = note*(-1);
        } else {
          angle[i] = angle[tunnel_connect[i]] + angle[tunnel_connect[i]]/abs(angle[tunnel_connect[i]])*shita
        }
      }
    }
    for(i in 1:num_lim){
      if(tunnel_length[i]>0){
        if(tunnel_connect[i]==0){
          beginx[i] = 0; beginy[i] = 0;
        } else {
          beginx[i] = beginx[tunnel_connect[i]] + (minx[i]-minx[tunnel_connect[i]])*cos(angle[tunnel_connect[i]]);
          beginy[i] = beginy[tunnel_connect[i]] + (minx[i]-minx[tunnel_connect[i]])*sin(angle[tunnel_connect[i]]);
        }
      }
    }
    beginx[beginx>0] <- beginx[beginx>0]-1
    plot(0, type="n", xlim=c(-15,25), ylim=c(-20,20),
         xlab="x (mm)", ylab="y (mm)", las=1)
    
    for(i in 1:num_lim){
      if(tunnel_length[i]>0){
        endx = beginx[i] + tunnel_length[i]*cos(angle[i])
        endy = beginy[i] + tunnel_length[i]*sin(angle[i])
        arrows(beginx[i],beginy[i],endx,endy, length=0)
      }
    }
    
    for(i in 1:N){
      if(ind_y[i]==1){
        points(ind_x[i]-1,0, col=c("grey","red","green","blue","orange","black","purple")[Arena[ind_x[i],ind_y[i]]], pch=19)
      }else{
        posx <- beginx[ind_y[i]] + (ind_x[i]-minx[ind_y[i]])*cos(angle[ind_y[i]])
        posy <- beginy[ind_y[i]] + (ind_x[i]-minx[ind_y[i]])*sin(angle[ind_y[i]])
        points(posx,posy, col=c("grey","red","green","blue","orange","black","purple")[Arena[ind_x[i],ind_y[i]]], pch=19)
      }
    }
    legend(x=-10, y=-7.5, legend=c("Move forward", "Excavation", "Backing", "Waiting", "loadong"), col=c("red","green","blue","orange", "purple"), pch=19)
  }
    RunSimulation <- function(N, trate, exrate, wrate, shape1, shape2, plot_image, dual){
      exrate <- trate + exrate
      wrate <- wrate + exrate
      ## Setting
      ind_status <- rep(0,N)    # go ahead (0), back (1), excavating (2), waiting (3)
      ind_x <- rep(1,N)         # individual position x
      ind_y <- rep(1,N)         # individual position tunnel (=y)
      ind_back <- rep(0,N)      # back distance (sample from random variable)
      ind_side <- rep(0,N)      # side digging or not
      ind_skip <- rep(FALSE, N) # skip or not (used for wait and swap calcuration)
      ind_load <- rep(0,N)      # loading (1) or not (0)
      
      # Arena
      # state: 0:empty, 1:sand, 2:going ahead individual, 3: excavating individual, 
      # 4: backing individuals, 5:waiting individuals, 6:not accessible, 7:loading (back)
      Arena <- matrix(1, length_lim, num_lim)
      Arena[1,1] <- 2
      Arena[1,2:num_lim] <- 6
      Arena[1,num_lim:(num_lim-9)] <- 0:9
      Arena[2,1] <- 0
      tunnel_connect <- rep(0,num_lim)
      
      new_tunnel <- 2 # id of next tunnel
      for(j in 1:1000){
        activeN <- sample(1:N, N, replace=F)
        for(i in activeN){
          if(ind_skip[i]){
            ind_skip[i] <- FALSE;
            next;
          }
          
          # in case out of tunnels
          if(ind_x[i] == 1){
            if(Arena[2,1]==0 && (Arena[3,1]!=4 && Arena[3,1]!=7) ){
              ind_x[i] = 2;
              Arena[ind_x[i], ind_y[i]] = 2;
            } else if (Arena[2,1]==9){
              Arena[2,1]=0;
            }
            
            # in case inside tunnel
          } else {
            
            # ind i is going ahead
            if(ind_status[i] == 0){
              # cell infront of the individual is
              # 0: empty (or side is empty)
              if( (Arena[ind_x[i]+1, ind_y[i]] == 0) ||
                  sum(Arena[ind_x[i]+1,]==0 & tunnel_connect == ind_y[i] & Arena[ind_x[i],] == 1)>0){
                if(Arena[ind_x[i], ind_y[i]] == 9){
                  Arena[ind_x[i], ind_y[i]] = 9;
                } else {
                  Arena[ind_x[i], ind_y[i]] = 0;
                }
                
                # only front is empty
                if(sum(Arena[ind_x[i]+1,]==0 & tunnel_connect == ind_y[i] & Arena[ind_x[i],] == 1) == 0){
                  ind_x[i] = ind_x[i]+1
                  # both empty
                } else if (Arena[ind_x[i]+1,ind_y[i]]==0){
                  r <- runif(1,0,1)
                  if(r < 0.5){
                    ind_x[i] = ind_x[i]+1
                  } else {
                    ind_y[i] = (1:num_lim)[Arena[ind_x[i]+1,]==0 & tunnel_connect == ind_y[i] & Arena[ind_x[i],] == 1]
                    ind_x[i] = ind_x[i]+1
                  }
                  # only side is empty
                } else {
                  ind_y[i] = (1:num_lim)[Arena[ind_x[i]+1,]==0 & tunnel_connect == ind_y[i] & Arena[ind_x[i],] == 1]
                  ind_x[i] = ind_x[i]+1
                }
                Arena[ind_x[i], ind_y[i]] = 2;
                
                # 9: loaded sand
              }else if(Arena[ind_x[i]+1, ind_y[i]] == 9 ||
                       sum(Arena[ind_x[i]+1,]==9 & tunnel_connect == ind_y[i] & Arena[ind_x[i],] == 1)>0){
                if(Arena[ind_x[i]+1, ind_y[i]] == 9){
                  ind_load[i] = 1;
                  ind_status[i] = 1;
                  ind_back[i] = ceiling(rbeta(1, shape1, shape2)*ind_x[i]);
                  Arena[ind_x[i], ind_y[i]] = 7;
                  Arena[ind_x[i]+1, ind_y[i]] = 0;
                } else {
                  ind_load[i] = 1;
                  ind_status[i] = 1;
                  ind_back[i] = ceiling(rbeta(1, shape1, shape2)*ind_x[i]);
                  Arena[ind_x[i], ind_y[i]] = 7;
                  Arena[ind_x[i]+1, (1:num_lim)[Arena[ind_x[i]+1,]==9 & tunnel_connect == ind_y[i] & Arena[ind_x[i],] == 1]] = 0;
                }
                # 1: sand
              } else if(Arena[ind_x[i]+1, ind_y[i]] == 1){
                if(ind_status[i] == 0){
                  ind_status[i] = 2;
                  Arena[ind_x[i], ind_y[i]] = 3;
                }
                
                # 2: going ahead individual
              } else if(Arena[ind_x[i]+1, ind_y[i]] == 2){
                
                # 3: excavating individual
              } else if(Arena[ind_x[i]+1, ind_y[i]] == 3){
                # in case no side tunnel
                if(sum(Arena[ind_x[i]+1,]!=1 & tunnel_connect == ind_y[i] & Arena[ind_x[i],] == 1) == 0){
                  r <- runif(1,0,1)
                  if(r < trate){
                  } else if(r < exrate){
                    ind_status[i] = 2;
                    ind_side[i] = new_tunnel;
                    tunnel_connect[new_tunnel] = ind_y[i];
                    new_tunnel <- new_tunnel + 1
                    Arena[ind_x[i], ind_y[i]] = 3;
                  } else if(r < wrate){
                    ind_status[i] = 3;
                    Arena[ind_x[i], ind_y[i]] = 5;
                  }
                }
                
                # 4: backing individual  
              } else if(Arena[ind_x[i]+1, ind_y[i]] == 4){
                Arena[ind_x[i], ind_y[i]] = 4;
                ind_status[i] = 1;
                ind_back[i] = ceiling(rbeta(1, shape1, shape2)*ind_x[i])
                
                # 7: loading individual
              } else if(Arena[ind_x[i]+1, ind_y[i]] == 7){
                swap <- (1:N)[(ind_x == ind_x[i]+1) & (ind_y == ind_y[i])]
                ind_load[swap] = 0;
                ind_status[swap] = 0;
                ind_back[swap] = 0;
                Arena[ind_x[swap], ind_y[swap]] = 2;
                
                ind_load[i] = 1;
                ind_status[i] = 1;
                ind_back[i] = ceiling(rbeta(1, shape1, shape2)*ind_x[i]);
                Arena[ind_x[i], ind_y[i]] = 7;
                ind_skip[i] <- TRUE
              }
              
              # ind i is backing
            } else if (ind_status[i] == 1){
              if(ind_back[i] > 0){
                if (Arena[ind_x[i]-1,ind_y[i]] == 1){
                  # back to the end of tunnel
                  ind_side[i] = 1;
                }
                
                # back straight
                if(ind_side[i] == 0){
                  
                  # back is empty
                  if(Arena[ind_x[i]-1, ind_y[i]] == 0 || (ind_x[i]-1) == 1){
                    Arena[ind_x[i], ind_y[i]] = 0;
                    ind_x[i] = ind_x[i]-1
                    if(ind_load[i] == 0){
                      Arena[ind_x[i], ind_y[i]] = 4;
                    } else {
                      Arena[ind_x[i], ind_y[i]] = 7;
                    }
                    ind_back[i] = ind_back[i]-1
                    
                    if(ind_x[i] == 1){
                      # come out of tunnel
                      Arena[ind_x[i], ind_y[i]] = 2;
                      ind_status[i] = 0;
                      ind_back[i] = 0;
                      ind_load[i] = 0;
                    } 
                    
                    # back is waiting ind
                  } else if (Arena[ind_x[i]-1, ind_y[i]] == 5){
                    swap <- (1:N)[(ind_x == ind_x[i]-1) & (ind_y == ind_y[i])]
                    ind_x[i] = ind_x[i]-1
                    ind_x[swap] = ind_x[swap]+1
                    if(ind_load[i] == 0){
                      Arena[ind_x[i], ind_y[i]] = 4;
                    } else {
                      Arena[ind_x[i], ind_y[i]] = 7;
                    }
                    Arena[ind_x[swap], ind_y[swap]] = 2;
                    ind_back[i] = ind_back[i]-1
                    ind_skip[swap] <- TRUE
                    ind_status[swap] <- 0
                  }
                  
                  # back side
                } else {
                  # back is empty
                  if(Arena[ind_x[i]-1,tunnel_connect[ind_y[i]]]==0){
                    Arena[ind_x[i], ind_y[i]] = 0;
                    ind_y[i] = tunnel_connect[ind_y[i]];
                    ind_x[i] = ind_x[i]-1
                    if(ind_load[i] == 0){
                      Arena[ind_x[i], ind_y[i]] = 4;
                    } else {
                      Arena[ind_x[i], ind_y[i]] = 7;
                    }
                    ind_back[i] = ind_back[i]-1
                    ind_side[i] = 0;
                    
                    # back is waiting ind
                  } else if (Arena[ind_x[i]-1,tunnel_connect[ind_y[i]]]==5){
                    swap <- (1:N)[(ind_x == ind_x[i]-1) & (ind_y == tunnel_connect[ind_y[i]])]
                    ind_y[swap] = ind_y[i]
                    ind_y[i] = tunnel_connect[ind_y[i]];
                    ind_x[i] = ind_x[i]-1
                    ind_x[swap] = ind_x[swap]+1
                    if(ind_load[i] == 0){
                      Arena[ind_x[i], ind_y[i]] = 4;
                    } else {
                      Arena[ind_x[i], ind_y[i]] = 7;
                    }
                    Arena[ind_x[swap], ind_y[swap]] = 2;
                    ind_back[i] = ind_back[i]-1
                    ind_skip[swap] <- TRUE
                    ind_status[swap] <- 0
                    ind_side[i] = 0;
                  }
                }
              } else {
                Arena[ind_x[i], ind_y[i]] = 9;
                ind_status[i] = 0;
                ind_load[i] = 0;
              }
              
              # ind i is excavating
            } else if (ind_status[i] == 2){
              if(ind_side[i] == 0){
                # excavate forward
                if(Arena[ind_x[i]+1, ind_y[i]] == 1){
                  Arena[ind_x[i]+1, ind_y[i]] = 0;
                  ind_load[i] = 1;
                  Arena[ind_x[i], ind_y[i]] = 7;
                  ind_back[i] = ceiling(rbeta(1, shape1, shape2)*ind_x[i])
                  ind_status[i] = 1;
                }
              } else {
                # excavate side
                if(Arena[ind_x[i]+1, ind_side[i]] == 1){
                  Arena[ind_x[i]+1, ind_side[i]] = 0;
                  ind_load[i] = 1;
                  Arena[ind_x[i], ind_y[i]] = 7;
                  ind_back[i] = ceiling(rbeta(1, shape1, shape2)*ind_x[i])
                  ind_status[i] = 1;
                  ind_side[i] = 0; 
                }
              }
              
              # ind i is waiting
            } else if (ind_status[i] == 3){ 
              
            }
          }
        }
        
        if(plot_image == 1){
          par(mfrow=c(1,1), pin=c(5,5))
          plot_tunnel(Arena, ind_x, ind_y, tunnel_connect, dual); Sys.sleep(0.1); 
        } else if (plot_image == 2){
          par(mfrow=c(1,1), pin=c(5,5))
          plot_tunnel(Arena, ind_x, ind_y, tunnel_connect, dual);
        }
        if(sum(Arena[length_lim,]==0) > 0){
          break;
        }
        options(warn=2)
      }
      #plot_tunnel(Arena, ind_x, ind_y, tunnel_connect, dual);
      
      return(list(Arena, ind_x, ind_y, tunnel_connect))
    }
  }
  
  ## Calcuration
  ## for number of branches
  {
    Res <- NULL
    for(iii in 1:1000){
      iter <- 16
      tipnum1 = tipnum2 = sumlength1 = sumlength2 <- NULL
      branch1 = branch2 = branch3 <- NULL
      #par(mfrow=c(5,10), pin=c(1,1))
      
      for(replicates in 1:iter){
        Parameter <- Paraneo_para
        res <- RunSimulation(N, Parameter[1], Parameter[2], Parameter[3], Parameter[4], Parameter[5],
                             plot_image = 0, dual = 1)
        
        tunnel_length <- apply(res[[1]][2:length_lim,]!=1, 2, sum)
        tipnum1 <- c(tipnum1, length(tunnel_length[tunnel_length>0]))
        tipnum2 <- c(tipnum2, length(tunnel_length[tunnel_length>1]))
        sumlength1 <- c(sumlength1, sum(tunnel_length[tunnel_length>0]))
        sumlength2 <- c(sumlength2, sum(tunnel_length[tunnel_length>1]))
        
        #plot_tunnel(res[[1]], res[[2]], res[[3]], res[[4]], dual=1)
      }
      Res<-c(Res,mean(tipnum1))
    }
    truehist(Res,breaks=seq(1,2,0.1))
    Para_tip <- Res
    #save(Para_tip, file="Para_tip.Rdata")
    #write.table(Res, "E:\\\\Para-tipnum.txt")
  }
  
  ## for tunnel length by branching segments
  {
    iter <- 1000
    branch1 = branch2 = branch3 <- NULL
    for(replicates in 1:iter){
      Parameter <- Paraneo_para
      res <- RunSimulation(N, Parameter[1], Parameter[2], Parameter[3], Parameter[4], Parameter[5],
                           plot_image = 0, dual = 1)
      res2 <- res[[1]][3:12,]
      res2[res2>1] <- 0
      res2
      for(i in (1:num_lim)[res[[4]]>0]){
        if(sum(res2[,i]==0)==0){
          res[[4]][i]<-0
        }
      }
      if(sum(res[[4]])==0){
        branch3 <- c(branch3, 10)
      } else if (sum(res[[4]])==1) {
        b1st <- min((1:10)[apply(res2[,2:num_lim]==0, 1, sum)>0]-1)
        branch1 <- c(branch1, b1st)
        b2nd <- apply(res2[(b1st+1):10,]==0, 2, sum)
        branch3 <- c(branch3, b2nd[b2nd>0])
      } else {
        for(i in 1:20){
          if(sum(res2[,i]==0)>0){
            if(sum(res[[4]]==i)>0){
              karimat <- matrix(1:10, 10, sum(res[[4]]==i)) * (res2[,res[[4]]==i] == 0)
              karimat[karimat==0] <- 100
              cut <- apply(karimat, 2, min)-1
              cut <- sort(cut)
              if(length(cut)==1){
                if(i==1){
                  branch1 <- c(branch1, cut)
                } else {
                  branch2 <- c(branch2, cut - (min((1:10)[res2[,i]==0])-1) )
                }
                branch3 <- c(branch3, max((1:10)[res2[,i]==0])-cut)
              } else {
                for(j in 1:length(cut)){
                  if(i == 1 && j == 1){
                    branch1 <- c(branch1, cut[1])
                  } else if (j == 1){
                    branch2 <- c(branch2, cut[1]-(min((1:10)[res2[,i]==0])-1) )
                  } else {
                    branch2 <- c(branch2, cut[j]-max(cut[1:(j-1)]))
                  }
                }
                branch3 <- c(branch3, max((1:10)[res2[,i]==0])-max(cut))
              }
            } else {
              branch3 <- c(branch3, sum(res2[,i]==0))
            }
          }
        }
      }
      #branch1
      #branch2
      #branch3
      #plot_tunnel(res[[1]], res[[2]], res[[3]], res[[4]], dual=1)
    }
    
    branch = c( rep(1,length(branch1)), rep(2,length(branch2)), rep(3,length(branch3)))
    length = c(branch1, branch2, branch3)
    d <- data.frame(branch, length)
    posn.d <- position_dodge(width = 0.5)
    ggplot(d, aes(branch, length)) + 
      stat_summary(geom = 'errorbar', position = posn.d, fun.data = mean_se, width=0.25) +
      stat_summary(geom = 'point', position = posn.d, fun.y = mean, size = 3) +
      theme_bw() + theme(aspect.ratio = 1)
    Para_branch <- d
    #save(Para_branch, file="Para_branch.Rdata")
    #write.table(d, "E:Paraneo_branch.txt")
  }
  
  ## to make gif animation
  if(F){
    Parameter <- Paraneo_para
    saveGIF(expr={ 
      RunSimulation(N, Parameter[1], Parameter[2], Parameter[3], Parameter[4],
                    plot_image = 2, dual = 0)
      }, interval = 0.5, movie.name = "E:Paraneo1.gif", loop=1)
  }
}