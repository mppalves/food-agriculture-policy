library("rootSolve")
rm(list=ls())

#parameters
bc <- 0.2
bm <- 0.4
g <- 1.75
a <- 0.4
d <- 2
si <- 0.45
t <- 1
w <- 1
r <- 0.5
pm <- 1.25
theta <- -0.85

#policy variables
s <- 0

ObjFOCs <- function(inputs) {
  xc = inputs[1]
  xm = inputs[2]
  xl = inputs[3]
  L = inputs[4]
  A = inputs[5]
  lamb = inputs[6]
  pc = inputs[7]
  
  #Utilitiy FOCs
  foc1 = bc*(xc^(bc-1))*(xm^bm)*(xl^(1-bc-bm)) - (lamb*pc)
  foc2 = bm*(xc^bc)*(xm^(bm-1))*(xl^(1-bc-bm)) - (lamb*pm)
  foc3 = (1-bc-bm)*(xc^bc)*(xm^bm)*(xl^(-bc-bm)) - (lamb*w)
  
  #Production FOC  
  foc4 = lamb*(1+s)*pc*g*a*(L^(a-1))*(A^si) - lamb*w
  foc5 = lamb*(1+s)*pc*g*si*(L^a)*(A^(si-1)) - lamb*r
  
  
  foc6 = (1+s)*pc*g*(L^a)*(A^si) - w*L + w*t - w*xl - r*a - pc*xc - pm*xm
  
  mcc = xc + g*(pc^theta) - g*(L^a)*(A^si)
  
  return(c(foc1,foc2,foc3,foc4,foc5, foc6, mcc))
}

statv <- c(1,1,1,1,1,1,1)

sol0 <- multiroot(ObjFOCs, statv);round(sol0$root,3)

#Solution a
s <-0.01


sola <- multiroot(ObjFOCs, statv);

#Solution b
s <-0
d <-0.05

solb <- multiroot(ObjFOCs, statv);

XC.base <- round(sol0$root[1])
XL.base <- round(sol0$root[2])
L.base <- round(sol0$root[3])
a.base <- round(sol0$root[4])

XC.A <- round(sola$root[1])
XL.A <- round(sola$root[2])
L.A <- round(sola$root[3])
a.A <- round(sola$root[4])


setwd("Z:/Home work/Homework 6.1")

Ans <- list(XC.base = XC.base, XL.base = XL.base, L.base = L.base, a.base = a.base, XC.A = XC.A ,XL.A = XL.A, L.A = L.A, a.A = a.A ,XC.B = XC.B, XL.B = XL.B, L.B = L.B, a.B = a.B)

write.csv(Ans, file = "FHM.csv", quote = F)