structural_similarity_measure_function <- function(SS){

KK = nrow(SS)
K_ToShow = KK*0.1
SS_M = SS[1:KK,1:KK]
SS_MUM = SS_M
#str(SS)
#str(SS_M)
index_Max = (which(SS_M == max(SS_M), arr.ind = TRUE))
index_Max = c(index_Max[1],index_Max[2])
index_Max_C = index_Max
# print(index_Max_C)
SS_Diag = SS_M[index_Max]
# str(index_Max)
for (i in 2:KK) {
  SS_M[index_Max_C[1],]=-1
  SS_M[,index_Max_C[2]]=-1
  index_Max_C= which(SS_M == max(SS_M), arr.ind = TRUE)
  index_Max_C = c(index_Max_C[1],index_Max_C[2])
  # print(index_Max_C)
  #Sys.sleep(0.5)
  SS_Diag = cbind(SS_Diag, SS_M[index_Max_C])
  index_Max = rbind(index_Max,index_Max_C)
}

PMR = matrix(0, nrow=KK, ncol=KK)
for(iPM in 1:KK){
  ToSetOne = index_Max[iPM,1]
  PMR[ToSetOne,iPM]=1
}
SS_PR = t(PMR) %*%  SS_MUM 
# levelplot(t(SS_PR))
# Sys.sleep(1)
# levelplot(t(SS_PR))
# Sys.sleep(1)

PMC = matrix(0, nrow=KK, ncol=KK)
for(iPM in 1:KK){
  ToSetOne = index_Max[iPM,2]
  PMC[iPM,ToSetOne]=1
}

SS_PC = SS_MUM %*% t(PMC) 
# levelplot(t(SS_MUM))
# Sys.sleep(1)
# levelplot(t(SS_PC))
# Sys.sleep(1)

SS_P = t(PMR) %*%  SS_MUM %*% t(PMC)
# levelplot(t(SS_MUM))
# Sys.sleep(1)

SS_P_ToShow = SS_P[1:(K_ToShow),1:(K_ToShow)]


# levelplot(t(SS_P_ToShow))
# levelplot(t(SS_P_ToShow), at = seq(0, 1, by=0.01))

# Sys.sleep(1)

# print(sum(diag(SS_P_ToShow))/K_ToShow)
# end_time = Sys.time()
# end_time - start_time
return(sum(diag(SS_P_ToShow))/K_ToShow)

## index_Max = cbind(index_Max,index_Max_C)
  
}


