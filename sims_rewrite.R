# Try to rewrite the simulations
rm(list = ls())
packs <- c("pracma","splines2","psych","caret","orthopolynom","parallel","mvtnorm","expm", "matrixStats", 
           "RCAL", "glmnet")
lapply(packs, library, character.only = TRUE)


### Start with the DGP ###
logitPropensity <-function(Z, S){
    dimZ <- dim(Z)[2]
    n <- dim(Z)[1]
    gamma <- c(-1, integer(S-1)+1, integer(dimZ - S))
    probs <- 1/(1 + exp(-1*(Z%*%gamma)))
    D <- rbinom(n,1,probs)
    D
}

probitPropensity <- function(Z, S){
    dimZ <- dim(Z)[2]
    n <- dim(Z)[1]
    gamma <- c(-1.5, integer(S-1)+1, integer(dimZ - S))
    probs <- pnorm(Z%*%gamma)
    D <- rbinom(n,1,probs)
    D
}

# Given Z, generate linear model
LinearModel <- function(D,Z,S,sigma= 1){
    dimZ <- dim(Z)[2]
    n <- dim(Z)[1]
    alpha <- c(integer(S)+1, integer(dimZ - S)) 
    eps <- sigma*rnorm(n)
    Y <- D*(Z%*%alpha) + eps
    Y
}

# Generate Y,D,P,Z 
DGP <- function(n = 500, dimZ = 50, K = 3, rho = 0, sigma = 1, S = 5, mP = 0, mL = 0, deg = 2){
    X <- runif(n-1,min=1,max=2)
    X <- c(X,1.5)
    varMatZ <- rho*toeplitz(2^seq(0,-(dimZ-4),by = -1)) + (1-rho)*diag(dimZ - 3)
    Zinit <- rmvnorm(n = n, sigma = varMatZ)
    ZinitMis <- Zinit + pmax(Zinit,2)**2
    ZinitMis <- (ZinitMis - mean(ZinitMis))/sd(ZinitMis)
    constant <- integer(n) + 1
    
    # Generate Various Mispecifications
    Z <- cbind(constant, X, 0.5*X^2, Zinit) # for if we use 0.5*X^2
    Zmis <- cbind(constant, X, 0.5*X^2, ZinitMis) # for if we use 0.5*X^2
    ZD <- cbind(constant, X, 0.5*X^2, Zinit) # for if we use 0.5*X^2
    ZmisD <- cbind(constant, X, 0.5*X^2, ZinitMis) # for if we use 0.5*X^2
    
    P <- bSpline(X, df = K, degree = deg, intercept = TRUE)
    D <- mP*logitPropensity(Zmis,S)+ (1-mP)*logitPropensity(Z,S)
    Y <- mL*LinearModel(D,Zmis,S,sigma) + (1-mL)*LinearModel(D,Z,S,sigma)
    list("Y" = Y, "Z" = Z, "D" = D, "X" = X, "Pmat" = P, "S" = S, "K" = K, "mP" = mP, "mL" = mL, "dimZ" = dimZ, "rho" = rho, "sigma" = sigma, "n" = n)
}

# Generate the IPW signal 
aIPW <- function(data, probs, means) {
    data$Y * data$D / probs - (data$D / probs - 1)* means
}

# Calibrated Estimation 
gHatBN <- function(data) {
    # Collect information from data
    Y <- data$Y
    D <- data$D
    X <- data$X
    Z <- data$Z
    K <- dim(data$Pmat)[2]
    n <- dim(Z)[1]
    dimZ <- data$dimZ
    ZnoCon <- scale(Z[,2:dimZ])
    
    # Estimate the weighted ATE for each K
    wIPWmat <- cbind()
    aIPWmat <- cbind()
    for (k in 1:K) {
        PW <- data$Pmat[,k] # weights for propensity estimation
        logitModel <- glm.regu.cv(3, 10, y = D, x = ZnoCon, iw = PW, loss = "cal")$sel.fit[,1]
        LW <- PW * D * (1 - logitModel) / logitModel
        linearPenalty <- cv.glmnet(ZnoCon, Y, weights = LW)
        linearModel <- predict(linearPenalty, newx = ZnoCon, s = "lambda.min")
       
        aIPWmat <- cbind(aIPWmat, aIPW(data, logitModel, linearModel)) 
        wIPWmat <- cbind(wIPWmat, PW * aIPW(data, logitModel, linearModel))
    }
    
    # Use this to estimate g
    Q <- t(data$Pmat) %*% data$Pmat / n
    Qinv <- solve(Q)
    beta <- Qinv %*% colMeans(wIPWmat)
    gHat <- data$Pmat %*% beta 
    
    # Also return an estimate of the variance
    epsMat <- cbind()
    for (k in 1:K) {
        epsMat <- cbind(epsMat, aIPWmat[,k] - gHat)
    }
    eps <- Y - gHat
    Omega <- Qinv %*%  t(epsMat * data$Pmat) %*% (epsMat * data$Pmat) %*% Qinv / n
    
    list("beta" = beta, "gHat" = gHat, "Omega" = Omega, "wIPW" = wIPWmat, "Qinv" = Qinv)
    
}

# Standard Estimation using ML loss functions
gHatML <- function(data) {
   # Collect information from data
    Y <- data$Y
    D <- data$D
    X <- data$X
    Z <- data$Z
    K <- data$K
    n <- data$n
    dimZ <- data$dimZ
    ZnoCon <- scale(Z[,2:dimZ])
    
    # Estimate nuisance models
    logitModel <- glm.regu.cv(3, 10, y = D, x = ZnoCon, loss = "ml")$sel.fit[,1]
    linearPenalty <- cv.glmnet(ZnoCon, Y, weights = D / mean(D))
    linearModel <- predict(linearPenalty, newx = ZnoCon, s = "lambda.min")
    
    wIPWmat <- cbind()
    aIPWmat <- cbind()
    for (k in 1:K) {
        PW <- data$Pmat[,k]
        aIPWmat <- cbind(aIPWmat, aIPW(data, logitModel, linearModel)) 
        wIPWmat <- cbind(wIPWmat, PW * aIPW(data, logitModel, linearModel))
    }
    
    # Use this to estimate g
    Q <- t(data$Pmat) %*% data$Pmat / n
    Qinv <- solve(Q)
    beta <- Qinv %*% colMeans(wIPWmat)
    gHat <- data$Pmat %*% beta 
    
    # Also return an estimate of the variance
    epsMat <- cbind()
    for (k in 1:K) {
        epsMat <- cbind(epsMat, aIPWmat[,k] - gHat)
    }
    eps <- Y - gHat
    Omega <- Qinv %*%  t(epsMat * data$Pmat) %*% (epsMat * data$Pmat) %*% Qinv / n
    
    list("beta" = beta, "gHat" = gHat, "Omega" = Omega, "wIPW" = wIPWmat, "Qinv" = Qinv)
}

### Now evaluate the estimates
coverage <- function(alpha, estimation, data) {
    # Collect some object from the data and estimation
    varMatrix <- estimation$Omega
    X <- data$X
    n <- data$n
    
    # Try inference for the psuedo-true parameter
    # g0 <- 1 + X + 0.5*(X**2)
    # pt <- lm(g0 ~ data$Pmat)$fitted.values
    
    # Calculate the standard error
    P1.5 <- (as.matrix(data$Pmat[n,]))
    se <- sqrt(as.numeric(t(P1.5)%*%varMatrix%*%P1.5))/sqrt(n)
    
    # see if the true value falls in the CI
    g1.5 <- 1 + 1.5 + 0.5 * (1.5)**2
    gHat1.5 <- estimation$gHat[n]
    CILower <- gHat1.5 - qnorm(1 - alpha/2)*se
    CIUpper <- gHat1.5 + qnorm(1 - alpha/2)*se
    covered <- 1*(g1.5 <= CIUpper & g1.5 >= CILower)
    
    list("covered" = covered, "se" = se)
}

boot_coverage <- function(alpha, estimation, data) {
  n <- data$n
  K <- data$K 
  X <- data$X 
  
  # Grab wIPWmat and Qinv from estimation
  Qinv = estimation$Qinv
  wIPWmat = estimation$wIPW
  P1.5 <- (as.matrix(data$Pmat[n,]))
  gHat1.5 <- estimation$gHat[n]
 
  # Multiplier Bootstrap the critical value 
  distribution = c()
  for (b in 1:1000) {
    d1 = dim(wIPWmat)[1]
    d2 = dim(wIPWmat)[2]
    mboots = rexp(n)
    mboots = mboots/mean(mboots)
    wIPW.b = sweep(wIPWmat, MARGIN = 1, mboots, "*")
    g.b = t(P1.5) %*% Qinv %*% colMeans(wIPW.b)
    distribution = c(distribution, g.b)
  }
  
  se = sd(distribution - gHat1.5)
  crit = as.numeric(quantile(abs(distribution - gHat1.5), 1 - alpha))
  CIlower = as.numeric(quantile(gHat1.5 - distribution, alpha/2))
  CIupper = as.numeric(quantile(gHat1.5 - distribution, 1 - alpha/2))
  if (CIlower > CIupper) {
    a = CIlower
    b = CIupper
    CIlower = b
    CIupper = a
  }
  g1.5 <- 1 + 1.5 + 0.5 * (1.5)**2
  covered = (g1.5 - gHat1.5 < CIupper)*(g1.5 - gHat1.5 > CIlower)
  
  list("covered" = covered, "se" = se)
}

coverageUniform <- function(estimation, data, alphaVal1 = 0.1, alphaVal2 = 0.05) {
    # Grab some information
    varMatrix <- estimation$Omega
    n <- data$n
    K <- data$K 
    X <- data$X
    seDenomVector <- sqrt(diag(data$Pmat%*%varMatrix%*%t(data$Pmat)))
    seVector <- data$Pmat%*%sqrtm(varMatrix) / seDenomVector
    supportPoints <- 40 #50 seemed to work ok
    distribution = c()
    for (i in 1:1000){
      sampledSeVec =  seVector[sample(1:n, supportPoints),]
      normals <- matrix(rnorm(supportPoints*K), supportPoints, K)
      bootstrap <- abs(rowSums(sampledSeVec*normals))
      bootstrap <- bootstrap[!is.na(bootstrap)]
      distribution = c(max(bootstrap), distribution)
    }
    
    # First Critical Value
    criticalValue1 <- quantile(distribution, 1-alphaVal1/2)
    sigmaX = seDenomVector/sqrt(n)
    CILower1 = estimation$gHat - criticalValue1*sigmaX
    CIUpper1 = estimation$gHat + criticalValue1*sigmaX
    g0 <- (1 + X + 0.5*X**2)[!is.na(CILower1)]
    CILower1 = CILower1[!is.na(CILower1)]
    CIUpper1 = CIUpper1[!is.na(CIUpper1)]
    uniformCoverage1 = 1*(mean( ((g0 >= CILower1) & (g0 <= CIUpper1)) ) > 0.99)
    
    # Second Critical Value
    criticalValue2 <- quantile(distribution, 1-alphaVal2/2)
    sigmaX = seDenomVector/sqrt(n)
    CILower2 = estimation$gHat - criticalValue2*sigmaX
    CIUpper2 = estimation$gHat + criticalValue2*sigmaX
    g0 <- (1 + X + 0.5*X**2)[!is.na(CILower2)]
    CILower2 = CILower2[!is.na(CILower2)]
    CIUpper2 = CIUpper2[!is.na(CIUpper2)]
    uniformCoverage2 = 1*(mean( ((g0 >= CILower2) & (g0 <= CIUpper2)) ) > 0.99)
    list("covered1" = uniformCoverage1, "criticalValue1" = criticalValue1/sqrt(n), "covered2" = uniformCoverage2, "criticalValue2" = criticalValue2/sqrt(n))
}

statsVarious <- function(data){
    # Run Both Estimation Procedures 
    estimationBN <- gHatBN(data)
    print("Finished BN Estimation")
    estimationSC<- gHatML(data)
    print("Finished SC Estimation")
    
    # Calculate some basic statistics
    gHatBN <- estimationBN$gHat
    gHatSC <- estimationSC$gHat
    X <- data$X
    g0 <- 1 + X + 0.5*X**2
    iMSE_BN = mean(gHatBN - g0)
    iVar_BN = mean((gHatBN - g0)^2)
    iMSE_SC = mean(gHatSC - g0)
    iVar_SC = mean((gHatSC - g0)^2)
    
    # Pointwise Coverage and Standard Errors for BN Estimation
    cov90_BN = boot_coverage(0.1, estimationBN, data)$covered
    bn95vec = boot_coverage(0.05, estimationBN, data)
    cov95_BN = bn95vec$covered
    se_BN = bn95vec$se
    
    # Uniform Coverage for BN Estimation
    covUniformBN = coverageUniform(estimationBN, data)
    cov90Uniform_BN = covUniformBN$covered1
    ucv90_BN = covUniformBN$criticalValue1
    cov95Uniform_BN = covUniformBN$covered2
    ucv95_BN = covUniformBN$criticalValue2
    
    # Pointwise Coverage and Standard Errors for SC Estimation 
    cov90_SC = boot_coverage(0.1, estimationSC, data)$covered
    sc95vec = boot_coverage(0.05, estimationSC, data)
    cov95_SC = sc95vec$covered
    se_SC = sc95vec$se
    
    # Uniform Coverage for SC Estimation
    covUniformSC = coverageUniform(estimationSC, data)
    cov90Uniform_SC = covUniformSC$covered1
    ucv90_SC = covUniformSC$criticalValue1
    cov95Uniform_SC = covUniformSC$covered2
    ucv95_SC = covUniformSC$criticalValue2
    
    # Report All Statistics
    list("iMSE_BN" = iMSE_BN, "iVar_BN" = iVar_BN, "iMSE_SC" = iMSE_SC, "iVar_SC" = iVar_SC, "cov90_BN" = cov90_BN, "cov95_BN" = cov95_BN, "se_BN" = se_BN, "cov90_SC" = cov90_SC, "cov95_SC" = cov95_SC, "se_SC"=se_SC, "cov90Uniform_BN" = cov90Uniform_BN, "cov95Uniform_BN" = cov95Uniform_BN, "ucv90_BN" = ucv90_BN, "ucv95_BN" = ucv95_BN, "cov90Uniform_SC" = cov90Uniform_SC, "cov95Uniform_SC" = cov95Uniform_SC, "ucv90_SC" = ucv90_SC, "ucv95_SC" = ucv95_SC)
}

### Finally, run simulations
runSims <- function(numSims, n = 500, dimZ = 50, K = 3, rho = 1, sigma = 1, S = 4, mP = 0, mL = 0, deg = 2){
  iMSE_BN = 0
  iVar_BN = 0
  iMSE_SC = 0
  iVar_SC = 0
  cov90_BN = 0
  cov95_BN = 0
  cov90Uniform_BN = 0
  cov95Uniform_BN = 0
  se_BN = 0
  ucv90_BN = 0
  ucv95_BN = 0
  cov90_SC = 0
  cov95_SC = 0
  cov90Uniform_SC = 0
  cov95Uniform_SC = 0
  ucv90_SC = 0
  ucv95_SC = 0
  se_SC = 0
  simscount = 1
  while (simscount <= numSims){
    print(simscount)
    data <- DGP(n = n, dimZ = dimZ, K = K, rho = rho, sigma = sigma, S = S, mP = mP, mL = mL, deg = deg)
    try({sim = statsVarious(data)
    iMSE_BN = iMSE_BN + sim$iMSE_BN/numSims
    iVar_BN = iVar_BN + sim$iVar_BN/numSims
    iMSE_SC = iMSE_SC + sim$iMSE_SC/numSims
    iVar_SC = iVar_SC + sim$iVar_SC/numSims
    cov90_BN = cov90_BN + sim$cov90_BN/numSims
    cov95_BN = cov95_BN + sim$cov95_BN/numSims
    cov90Uniform_BN = cov90Uniform_BN + sim$cov90Uniform_BN/numSims
    cov95Uniform_BN = cov95Uniform_BN + sim$cov95Uniform_BN/numSims
    ucv90_BN = ucv90_BN + sim$ucv90_BN/numSims
    ucv95_BN = ucv95_BN + sim$ucv95_BN/numSims
    se_BN = se_BN + sim$se_BN/numSims
    cov90_SC = cov90_SC + sim$cov90_SC/numSims
    cov95_SC = cov95_SC + sim$cov95_SC/numSims
    cov90Uniform_SC = cov90Uniform_SC + sim$cov90Uniform_SC/numSims
    cov95Uniform_SC = cov95Uniform_SC + sim$cov95Uniform_SC/numSims
    ucv90_SC = ucv90_SC + sim$ucv90_SC/numSims
    ucv95_SC = ucv95_SC + sim$ucv95_SC/numSims
    se_SC = se_SC + sim$se_SC/numSims
    adjust = numSims/simscount
    print(paste("n:", n, "mP:", mP, "mL:", mL, "K:", K, "S:", S,"dimZ:",dimZ, sep = " "))
    print(c(round(cov90_BN*adjust, 3), round(cov95_BN*adjust, 3), round(cov95Uniform_BN*adjust, 3)))
    print(c(round(cov90_SC*adjust, 3), round(cov95_SC*adjust, 3), round(cov95Uniform_SC*adjust, 3)))
    simscount = simscount + 1
    })
  }
  list("iMSE_BN" = iMSE_BN, "iVar_BN" = iVar_BN, "iMSE_SC" = iMSE_SC, "iVar_SC" = iVar_SC, "cov90_BN" = cov90_BN,
       "cov95_BN" = cov95_BN, "se_BN" = se_BN, "cov90_SC" = cov90_SC, "cov95_SC" = cov95_SC,"se_SC"= se_SC,
       "cov90Uniform_BN" = cov90Uniform_BN, "cov95Uniform_BN" = cov95Uniform_BN, "ucv90_BN" = ucv90_BN,
       "ucv95_BN" = ucv95_BN, "cov90Uniform_SC" = cov90Uniform_SC, "cov95Uniform_SC" = cov95Uniform_SC, 
       "ucv90_SC" = ucv90_SC, "ucv95_SC" = ucv95_SC, "n" = n, "dimZ" = dimZ, "rho" = rho, 
       "K" = K, "sparsity" = S, "mP" = mP, "mL" = mL, "numSims" = numSims)
}

numSims = 500
simsZ = 50
simsS = 7
simsK = 5
simsDeg = 2
simsRho = 1

# Misspecified Propensity, Correctly Specified Linear
# first = runSims(numSims = numSims, n = 400, dimZ = simsZ, K = simsK, rho = simsRho, S = simsS, mP = 1, mL = 0, deg = simsDeg)

# Correctly Specified Propensity, Misspecified Linear
# second = runSims(numSims = numSims, n = 400, dimZ = simsZ, K = simsK, rho = simsRho, S = simsS, mP = 0, mL = 1, deg = simsDeg)

# No misspecification
# third = runSims(numSims = numSims, n = 400, dimZ = simsZ, rho = simsRho, K = simsK, S = simsS, mP = 0, mL = 0, deg = simsDeg)

# Collect all results, write to .csv (if final)
# attempt_400 <- apply(data.frame(t(cbind(first, second, third))), 2, as.character)
# name = paste(400,"n_",simsRho,"rho_",simsDeg,"deg_",simsZ,"Z_",simsS,"S_",simsK,"K.csv",sep="")
# write.csv(attempt_400, file = name)

# Misspecified Propensity, Correctly Specified Linear
# first_800 = runSims(numSims = numSims, n = 800, dimZ = simsZ, K = simsK, rho = simsRho, S = simsS, mP = 1, mL = 0, deg = simsDeg)

# Correctly Specified Propensity, Misspecified Linear
# second_800 = runSims(numSims = numSims, n = 800, dimZ = simsZ, K = simsK, rho = simsRho, S = simsS, mP = 0, mL = 1, deg = simsDeg)

# No misspecification
# third_800 = runSims(numSims = numSims, n = 800, dimZ = simsZ, K = simsK, rho = simsRho, S = simsS, mP = 0, mL = 0, deg = simsDeg)

# Collect all results, write to .csv (if final)
# attempt_800 <- apply(data.frame(t(cbind(first_800, second_800, third_800))), 2, as.character)
# name = paste(800,"n_",simsRho,"rho_",simsDeg,"deg_",simsZ,"Z_",simsS,"S_",simsK,"K.csv",sep="")
# write.csv(attempt_800, file = name)

# Misspecified Propensity, Correctly Specified Linear
# first_2500 = runSims(numSims = numSims, n = 2500, dimZ = simsZ, K = simsK, rho = simsRho, S = simsS, mP = 1, mL = 0, deg = simsDeg)

# Correctly Specified Propensity, Misspecified Linear
# second_2500 = runSims(numSims = numSims, n = 2500, dimZ = simsZ, K = simsK, rho = simsRho, S = simsS, mP = 0, mL = 1, deg = simsDeg)

# No misspecification
# third_2500 = runSims(numSims = numSims, n = 2500, dimZ = simsZ, K = simsK, rho = simsRho, S = simsS, mP = 0, mL = 0, deg = simsDeg)

# Collect all results, write to .csv (if final)
# attempt_2500 <- apply(data.frame(t(cbind(first_2500, second_2500, third_2500))), 2, as.character)
# name = paste(2500,"n_",simsRho,"rho_",simsDeg,"deg_",simsZ,"Z_",simsS,"S_",simsK,"K.csv",sep="")
# write.csv(attempt_1000, file = name)



