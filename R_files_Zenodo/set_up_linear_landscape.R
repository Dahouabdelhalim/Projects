





# Set up linear landscape


disp_prob = 0.5*dgeom(0:(num_patches-1),prob=dispersal_probability_param)
#disp_prob = 0.5*ddisclap(0:(num_patches-1),p=dispersal_probability_param)
#disp_prob = c(0.25,0.125,0.0625,0.03,0.02,0.0125)
#disp_prob = c(0.4,0.03,0.02,0.01,0.0001,0.00001)#rep(0.25,3)

connectivity_matrix = matrix(0,num_patches,num_patches)
diag(connectivity_matrix) = 2*disp_prob[1]

for (index in 1:length(disp_prob)){
  
  connectivity_matrix[row(connectivity_matrix) == col(connectivity_matrix)+index] = disp_prob[index+1]
  connectivity_matrix[row(connectivity_matrix) == col(connectivity_matrix)-index] = disp_prob[index+1]
  
}

connectivity_matrix = rbind(matrix(0,nrow=extra_patches,ncol=num_patches),connectivity_matrix)
connectivity_matrix = rbind(connectivity_matrix,matrix(0,nrow=extra_patches,ncol=num_patches))
connectivity_matrix = cbind(matrix(0,ncol=extra_patches,nrow=num_patches+2*extra_patches),connectivity_matrix)
connectivity_matrix = cbind(connectivity_matrix,matrix(0,ncol=extra_patches,nrow=num_patches+2*extra_patches))
#connectivity_matrix[1,1]=0.75
#connectivity_matrix[10,10]=0.75

#image(t(connectivity_matrix))

connectivity=connectivity_matrix

#connectivity_matrix 

avg_disp_dist = sum(disp_prob*(0:(num_patches-1)))

connectivity_matrix[is.na(connectivity_matrix)]=0

