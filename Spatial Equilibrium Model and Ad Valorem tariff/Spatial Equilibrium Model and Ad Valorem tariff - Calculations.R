rm(list=ls()) #will remove ALL objects

#******************************************************************************
#   Codes for Problem Set 3 for the Spatial Equilibrium Model
#******************************************************************************

#Parameters
ALPHA <- c(100, 115,125, 110) #Demand function intercepts q <- f(p),
BETA <- c(8,6,6.5,4)          #Demand function slopes (absolute value) q <- f(p) 
OMEGA <- c(-40,-15,-26,-16) #Supply function intercepts q <- f(p)   
GAMMA <- c(6,5,2,1.6)    #Supply function slopes q <- f(p)      

n <- length(ALPHA) # number of regions

# Transport cost matrix
R1 <- c(0,2,4.1,3.2)
R2 <- c(2.2,0,3.2,3.3)      
R3 <- c(3.3,2.9,0,3.4)      
R4 <- c(3.1,3.2,4,0)
#R5 <- c(3.8,2.8,3.1,4.5,0.0)
TC <- matrix(c(R1,R2,R3,R4),nrow = n,byrow = TRUE)

# Converting the intercept and slope parameter in p = f(q) space
A <- ALPHA/BETA   # Demand function intercepts p = f(q)
B <- 1/BETA       # Demand function slopes p = f(q) 
C <- -OMEGA/GAMMA # Supply function intercepts p = f(q)
E <- 1/GAMMA      # Supply function slopes p = f(q)

# Ad valorem tariff matrix
R1 <- c(0.0,0.0,0.0,0)
R2 <- c(0.0,0.0,0.0,0)      
R3 <- c(0.0,0.0,0.0,0)      
R4 <- c(0.0,0.0,0.0,0)
#R5 <- c(0.0,0.0,0.0,0.1,0.0)
AVTh <- matrix(c(R1,R2,R3,R4),nrow = n,byrow = TRUE)
AVTh

# Specific tariff matrix
R1 <- c(0.0,0.0,1.5,0.0)
R2 <- c(0.0,0.0,1.5,0.0)      
R3 <- c(0.0,0.0,0.0,0.0)      
R4 <- c(0.0,0.0,1.5,0.0)
#R5 <- c(0.0,0.0,2,0.0,0.0)
STh <- matrix(c(R1,R2,R3,R4),nrow = n,byrow = TRUE)
STh

###############################################
# NLP
###############################################

# function for objective function for NLP
fcnObj <- function(x){
  
  Qd <- x[1:n]
  Qs <- x[(n+1):(2*n)]
  lamd <- x[(2*n+1):(3*n)]
  lams <- x[(3*n+1):(4*n)]
  Ship <- matrix(x[(4*n+1):(4*n+n**2)],nrow = n,byrow = TRUE)
  
  QUASIWELF <- sum(A*Qd - (1/2)*B*(Qd**2) - C*Qs - (1/2)*E*(Qs**2)) - sum(TC*Ship) -
    sum(ST*Ship) -
    sum(Ship*(lamd - lams)) + 
    sum(Ship*(lamd/(1 + AVT) - lams))
  
    return(signf*QUASIWELF)
}

# function for constraints for NLP
InEqConst <- function(x){
  # InEqConst >= 0
 
  Qd <- x[1:n]
  Qs <- x[(n+1):(2*n)]
  lamd <- x[(2*n+1):(3*n)]
  lams <- x[(3*n+1):(4*n)]
  Ship <- matrix(x[(4*n+1):(4*n+n**2)],nrow = n,byrow = TRUE)
  
  Sconst <- Qs - rowSums(Ship) # rowSums for each row sum across columns
  Dconst <- colSums(Ship) - Qd # colSums for each col sum across rows
  
  Dprice <- lamd - (A - B*Qd)
  Sprice <- C + E*Qs - lams
  
  k <- 0;
  Plink <- rep(1,(n*n))
  for(i in 1:n){
    for(j in 1:n){
      k <- k + 1 # counter
      Plink[k] <- (TC[i,j] + ST[i,j] + lams[i])*(1 + AVT[i,j]) - lamd[j]
    }
  }
  
  return(c(Dconst,Sconst,Dprice,Sprice,Plink,Qd,Qs,c(Ship),lamd,lams))
}


stval <- c(9.97047485, 21.76240207, 44.05548251, 55.20292658, 39.96232243, 42.0954508, 42.2145195, 6.80738731, 34.74989746, 45.08635338, 10.00678127, 11.49793299, 13.10693175, 14.49953391, 10.00762899, 10.00891785, 11.50027079, 13.1074402, 14.49996479, 10.00977651, 9.97047482, 0, 32.124777, 0.00019898, 0, 1e-08, 21.76240204, 9e-08, 20.45211736, 0, 1e-08, 0, 6.80738729, 1e-08, 0, 0, 0, 3e-08, 34.74989743, -1e-08, 0, 3e-08, 5.12331811, 0.0007128, 39.96232244, 0, 0, 5.079, 0.003, 39.985)
AVT <- AVTh*0
ST <- STh*0

signf <- -1 # because 'auglag' is a minimization algorithm
solU <- alabama::auglag(par = stval,  fn = fcnObj, hin = InEqConst)
solnFT <- round(solU$par,digits = 3)

solnFTh <- paste(solnFT, collapse=", ")
solnFTh

QdFT <- solnFT[1:n];QdFT
QsFT <- solnFT[(n+1):(2*n)];QsFT
lamdFT <- solnFT[(2*n+1):(3*n)];lamdFT
lamsFT <- solnFT[(3*n+1):(4*n)];lamsFT
ShipFT <- matrix(solnFT[(4*n+1):(4*n+n**2)],nrow = n,byrow = TRUE);ShipFT
PdFT <- A - B*QdFT;PdFT
PsFT <- C + E*QsFT;PsFT

STRevFT <- colSums(ST*ShipFT);
AVTRevFT <- colSums(AVT*ShipFT*lamsFT);AVTRevFT

# Tariff and Ad Velorem
AVT <- AVTh
ST <- STh

stval <- solnFT
signf <- -1 # because 'auglag' is a minimization algorithm
solU <- alabama::auglag(par = stval,  fn = fcnObj, hin = InEqConst)
solnA1 <- round(solU$par,digits = 3)

QdA1 <- solnA1[1:n];QdA1
QsA1 <- solnA1[(n+1):(2*n)];QsA1
lamdA1 <- solnA1[(2*n+1):(3*n)];lamdA1
lamsA1 <- solnA1[(3*n+1):(4*n)];lamsA1
ShipA1 <- matrix(solnA1[(4*n+1):(4*n+n**2)],nrow = n,byrow = TRUE);ShipA1
PdA1 <- A - B*QdA1;PdA1
PsA1 <- C + E*QsA1;PsA1

STRevA1 <- colSums(ST*ShipA1);STRevA1
AVTRevA1 <- colSums(AVT*ShipA1*lamsA1);AVTRevA1

# Computation of producer surplus, consumer surplus, and gains from trade  

ChPS <- 0.5*(PsA1 - PsFT)*(QsFT + QsA1);ChPS
ChCS <- 0.5*(PdFT - PdA1)*(QdFT + QdA1);ChCS
ChGovRev <- (AVTRevA1 - AVTRevFT) + (STRevA1 - STRevFT);ChGovRev
GFT      <- ChCS + ChPS + ChGovRev;GFT
sum(GFT)

# Exporting Answers in CSV file -------------------------------------------
AnsQPFT <- cbind(QdFT, QsFT, PdFT, PsFT)
AnsShipFT <- ShipFT

AnsQPA1 <- cbind(QdA1, QsA1, PdA1, PsA1)
AnsShipA1 <- ShipA1
AnsWel <- cbind(ChPS,ChCS,ChGovRev,GFT)

PerChShip <- round(100*(ShipA1/ShipFT - 1),2);PerChShip

setwd("Z:/Home work/Homework 5.2")

Ans <- list(FTQP=AnsQPFT, FTShip = AnsShipFT,A1QP=AnsQPA1, A1Ship = AnsShipA1, Welfare = AnsWel)
write.csv(Ans, file = "AnswersTariff.csv", quote = F)

