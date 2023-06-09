# Morphometric variation at different spatial scales
# Functions and worked example with adult Homo sapiens

# ---- Packages ----

library(geomorph)  # for worked example only
library(Morpho)  # for functions and worked example

# ---- Functions ----

# Mardia-Dryden distribution for 2D landmark coordinates
md.distri <- function (m_mshape, n, sd = 0.05) {
  
  # Objectives:
  #   Create a matrix of 2D shape coordinates drawn from a Mardia-Dryden distribution
  #   (3D not implemented)
  # Input:
  #   - m_mshape: k x 2 refence matrix (usually the sample mean shape), 
  #   where k is the number of 2D landmarks
  #   - n: number of observations
  #   - sd: standard deviation (default = 0.02)
  # Output:
  #   - x_md: n x 2k matrix of shape coordinates drawn from a Mardia-Dryden distribution
  
  k <- nrow(m_mshape)  # landmark number
  
  x_md <- matrix(NA, nrow = n, ncol = 2 * k)
  colnames(x_md) <- c(paste(1:k,"X", sep = ""), paste(1:k,"Y", sep = ""))
  for (i in 1:k) {
    x_md[, i] <- rnorm(n, mean = m_mshape[i, 1], sd)  # X coordinates
    x_md[, k + i] <- rnorm(n, mean = m_mshape[i, 2], sd)  # Y coordinates
  }
  
  return(x_md)
  
}

# Self-silmilar distribution for 2D landmark coordinates
ssim.distri <- function (m_mshape, n, sd = 0.02, f = 1) {
  
  # Objectives:
  #   Create a matrix of 2D shape coordinates drawn from a self-similar distribution
  #   (3D not implemented)
  # Input:
  #   - m_mshape: k x 2 refence matrix (usually the sample mean shape), 
  #   where k is the number of 2D landmarks
  #   - n: number of observations
  #   - sd: standard deviation (default = 0.02)
  #   - f: scaling factor
  # Output:
  #   - Xdefl: n x 2k matrix of shape coordinates drawn from a self-similar distribution
  
  k <- nrow(m_mshape)  # landmark number
  
  # Create the bending energy matrix as specified in Bookstein (1989)
  m_L <- CreateL(m_mshape)
  m_kxk <- as.matrix(m_L$Lsubk)  # bending energy (BE) matrix
  
  # Take only the (k-3) non-zero eigenvalues and the corresponding eigenvectors
  m_eigen <- eigen(m_kxk)
  m_PW <- m_eigen$vectors[, 1:(k-3)]  # eigenvectors of the BE matrix
  m_be <- m_eigen$values[1:(k-3)]  # eigenvalues of the BE matrix
  
  # Compute a Mardia-Dryden distribution and remove reference shape
  m_x <- matrix(NA, nrow = n, ncol = k)
  m_y <- matrix(NA, nrow = n, ncol = k)
  for (i in 1:k) {
    m_x[, i] <- rnorm(n, mean = 0, sd)
    m_y[, i] <- rnorm(n, mean = 0, sd)
  }
  
  # Self-similar distribution
  m_x_defl <- m_x %*% m_PW %*% diag(m_be ^ (-0.5)) %*% t(m_PW)
  m_y_defl <- m_y %*% m_PW %*% diag(m_be ^ (-0.5)) %*% t(m_PW)
  for (i in 1:n){
    m_x_defl[i, ] <- f * m_x_defl[i, ] + m_mshape[, 1] # add the reference shape
    m_y_defl[i, ] <- f * m_y_defl[i, ] + m_mshape[, 2] # add the reference shape
  }
  colnames(m_x_defl) <- paste(1:k, "X", sep = "")
  colnames(m_y_defl) <- paste(1:k, "Y", sep = "")
  Xdefl <- cbind(m_x_defl, m_y_defl)
  
  return(Xdefl)
  
}

# Function to compute bending energy (BE), principal warps, partial warps scores, etc.
create.pw.be <- function (m_overall, m_mshape) {
  
  # Objectives: 
  #   Computes the partial warps, the partial warps scores and the bending energy 
  #   for the 2D landmark coordinates (3D not implemented)
  # Input:
  #   - m_overall: input k x 2 x n array, where k is the number of 2D landmarks, 
  #   and n is the sample size.
  #   - m_mshape: k x 2 refence matrix (usually the sample mean shape), 
  #   where k is the number of 2D landmarks
  # Output:
  #   - bendingEnergy: bending energy (the (k-3) eigenvalues of the bending energy matrix)
  #   - principalWarps: k x (k-3) matrix of principal warps 
  #   (the k eigenvectors of the bending energy matrix)
  #   - partialWarpScores: n x (2k-6) matrix of partial warp 
  #   (the projection of the vectors of shape coordinates, 
  #   expressed as deviations from the reference shape, onto the principal warps)
  #   - variancePW: variance of the (k-3) partial warps
  #   - Xnonaf: n x 2k matrix of the non-affine component of shape variation
  # Reference:
  #   Bookstein FL. 1989. Principal Warps: Thin-plate splines and the decomposition of deformations. 
  #   IEEE Transactions on pattern analysis and machine intelligence 11(6).
  
  k <- dim(m_overall)[[1]]  # landmark number
  nspec <- dim(m_overall)[[3]]  # specimen number
  
  # Create the bending energy matrix as specified in Bookstein (1989)
  m_L <- CreateL(m_mshape)
  m_kxk <- as.matrix(m_L$Lsubk)  # bending energy matrix
  
  # Take only the (k-3) non-zero eigenvalues and the corresponding eigenvectors
  m_eigen <- eigen(m_kxk)
  m_PW <- m_eigen$vectors[, 1:(k-3)]  # principal warps
  m_be <- m_eigen$values[1:(k-3)]  # bending energy
  
  # Deviation of the landmark coordinates from the reference shape
  m_x <- matrix(NA, nrow = nspec, ncol = k)
  m_y <- matrix(NA, nrow = nspec, ncol = k)
  for (i in 1:nspec){
    m_x[i, ] <- m_overall[, 1, i] - m_mshape[, 1]
    m_y[i, ] <- m_overall[, 2, i] - m_mshape[, 2]
  }
  
  # Partial warps
  m_x_PWscores <- m_x %*% m_PW  # partial warp scores for the X dimension
  m_y_PWscores <- m_y %*% m_PW  # partial warp scores for the Y dimension
  colnames(m_x_PWscores) <- paste(1:(k-3), "X", sep = "")
  colnames(m_y_PWscores) <- paste(1:(k-3), "Y", sep = "")
  m_PWscores <- cbind(m_x_PWscores, m_y_PWscores)
  if (!is.null(dimnames(m_overall)[[3]])) {
    rownames(m_PWscores) <- dimnames(m_overall)[[3]]
  }
  
  # Non-affine component of shape variation
  m_x_nonaf <- m_x_PWscores %*% t(m_PW)
  m_y_nonaf <- m_y_PWscores %*% t(m_PW)
  for (i in 1:nspec){
    m_x_nonaf[i, ] <- m_x_nonaf[i, ] + m_mshape[, 1] # add the reference shape
    m_y_nonaf[i, ] <- m_y_nonaf[i, ] + m_mshape[, 2] # add the reference shape
  }
  m_nonaf <- cbind(m_x_nonaf, m_y_nonaf)
  
  # Variance of each partial warp (mean of the variances for X and Y coordinates)
  varPW <- rep(NA, (k-3))
  for (j in 1:(k-3)) { 
    var_x <- var(m_PWscores[, j])
    var_y <- var(m_PWscores[, (k - 3 + j)])
    varPW[j] <- (var_x + var_y) / 2
  }
  names(varPW) <- paste("PW", 1:(k-3), sep = "")
  
  # Output
  results <- list(bendingEnergy = m_be, 
                  principalWarps = m_PW,
                  partialWarpScores = m_PWscores,
                  variancePW = varPW,
                  Xnonaf = m_nonaf)
  
  return(results)
  
}

# ---- Worked example ----

# Landmark coordinates (curve semi-landmarks are already slid)
homo <- as.matrix(read.table("Landmarks_Homo.txt"))
homo_ar <- arrayspecs(homo, 87, 2)  # 87 two-dimensional landmarks
n_spec <- dim(homo_ar)[[3]]  # number of specimens
k <- dim(homo_ar)[[1]]  # number of landmarks

# Procrustes registration
homo_gpa <- procSym(homo_ar)
m_overall <- homo_gpa$rotated  # Procrustes coordinates
m_mshape <- homo_gpa$mshape  # average shape

# Visualization of the mean shape
plot(m_mshape, asp = 1, main = "Average shape", xlab = "X", ylab = "Y")

# Computation of BW, PW scores etc.
homo_be_pw <- create.pw.be(m_overall, m_mshape)

# Computation of log BE^-1 for the (k-3) partial warps
logInvBE <- log((homo_be_pw$bendingEnergy)^(-1))
# Computation of log PW variance for the (k-3) partial warps
logPWvar <- log(homo_be_pw$variancePW)
# Linear regression of the log PW variance on the log BE^-1
mod <- lm(logPWvar ~ logInvBE)
# Plot log PW variance on log BE^-1 with regression line
plot(logInvBE, logPWvar, col = "white", asp = 1,
     main = "PW variance against inverse BE", xlab = "log 1/BE", ylab = "log PW variance")
text(logInvBE, logPWvar, labels = names(logPWvar), cex = 0.5)
abline(mod, col = "blue")

# Non-affine component of shape distribution
tr_nonaf <- sum(diag(t(homo_be_pw$Xnonaf) %*% homo_be_pw$Xnonaf))  # trace of t(Xnonaf) %*% Xnonaf
plot(homo_be_pw$Xnonaf[, 1:k], homo_be_pw$Xnonaf[, (k+1):(2*k)], asp = 1, las = 1, cex = 0.5, 
     main = "Non-affine component", xlab = "X", ylab = "Y")

# Self-similar distribution
Xdefl <- ssim.distri(m_mshape, n = n_spec, sd = 0.05, f = 1)
tr_defl <- sum(diag(t(Xdefl) %*% Xdefl))  # trace of t(Xdefl) %*% Xdefl
plot(Xdefl[, 1:k], Xdefl[, (k+1):(2*k)], asp = 1, las = 1, cex = 0.5, 
     main = "Self-similar distribution", xlab = "X", ylab = "Y")

# Mardia-Dryden distribution
Xmd <- md.distri(m_mshape, n = n_spec, sd = 0.005)
plot(Xmd[, 1:k], Xmd[, (k+1):(2*k)], asp = 1, las = 1, cex = 0.5, 
     main = "Mardia-Dryden distribution", xlab = "X", ylab = "Y")
