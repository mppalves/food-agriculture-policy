####Parameters####
library("rootSolve")
library('ggplot2')
library('alabama')
library('ploty')

rm(list=ls())

a0a = 1
a1a = 3
a2a = -0.9
b0a = 5
b1a = -0.5
b2a = 1.2
a0b = 11
a1b = 1.2
a2b = -0.8
b0b = 40
b1b = -0.5 
b2b = 0.9
ia = 10
ib = 10
wa = 2
wb = 2


####Functions####

Sa = function(x){
  a = ((a0a*(x)^a1a)*(wa)^a2a)
}
Da = function(x){
  b = ((b0a*(x)^b1a)*(ia)^b2a)
}
Sb = function(x){
  c = ((a0b*(x)^a1b)*(wb)^a2b)
}
Db = function(x){
  d = ((b0b*(x)^b1b)*(ib)^b2b)
}
ES = function(x){
  d = ((a0a*(x)^a1a)*(wa)^a2a) - ((b0a*(x)^b1a)*(ia)^b2a)
}
ED = function(x){
  b = ((b0b*(x)^b1b)*(ib)^b2b) - ((a0b*(x)^a1b)*(wb)^a2b)
}

WCCB = function(x) {
  ib^b2b*b0b*((x^(b1b+1))/(b1b+1))
}

WPCB = function(x) {
  wb^a2b*a0b*((x^(a1b+1))/(a1b+1))
}

####World Equilibrium####

#Question b)

#finding the world price
autarkyPriceWorld = function(inputs) {
  pw=inputs[1]
  qw=inputs[2]
  
  ES = Sa(pw) - Da(pw) - qw
  ED = Db(pw) - Sb(pw) - qw
  
  return(c(ES,ED))
}

solw = multiroot(autarkyPriceWorld, c(1,1))
#quantities/prices
pw = solw$root[1];pw
qw = solw$root[2];qw

#World quantities/prices
pw
qw



#########################################################################
############################### Problem 1 ###############################

# (a) Elasicity ---------------------------------------------------------

ees <- (a0a*a1a*pw^(a1a-1)*wa^a2a - b0a*b1a*pw^(b1a-1)*ia^b2a)*pw/qw; ees
eed <- (b0b*b1b*pw^(b1b-1)*ib^b2b - a0b*a1b*pw^(a1b-1)*wb^a2b)*pw/qw; eed

# answer for 1.(a): 
# elasticities of excess supply = 4.471335 (country A)
# elasticities of excess demand = -1.665381 (country B)

# (b) absolute value of elasticity --------------------------------------

abs(ees) 
abs(eed)

# answer for 1.(b)
# The importing country (country B) will pay higher portion of the tariff. 
# Because country B has less elestic response (in absolute value) than country A.

#########################################################################
############################### Problem 2 ###############################
# NLP

# objective function
fmax <- function(inputs){
  tmax <- inputs[1]
  pamax <- inputs[2]
  
  pbmax<- pamax + tmax
  
  TR <- tmax*ED(pbmax)

  return(-1*TR)
}

# constraints

Eqconst<- function(x){
  tmax <-x[1]
  pamax<-x[2]
  
  pbmax<- pamax+tmax
  
  QED <- ED(pbmax)
  QES <- ES(pamax)
  
  res<- QED-QES #equilibrium
  
  return(res)
  
}

optimum <- auglag(par=c(1,1),fn=fmax,heq=Eqconst)
# Response ------------------------
TmaxR = optimum$par[1];TmaxR
pamaxR = optimum$par[2];pamaxR
pbmaxR = pamaxR + TmaxR;pbmaxR
TR.MaxR = TmaxR*ED(pbmaxR);TR.MaxR

#########################################################################
############################### Problem 3 ###############################

fmaxW <- function(inputs){
  tmax <- inputs[1]
  pamax <- inputs[2]
  
  pbmax <- pamax + tmax
  
  Net.WconsumersB <- WCCB(pw) - WCCB(pbmax)
  Net.WproducersB <- WPCB(pbmax) - WPCB(pw)
  TR = tmax*ED(pbmax);TR
  
  res <- Net.WconsumersB + Net.WproducersB + TR
  
  return(-res)
}

# constraints

EqconstW<- function(x){
  tmax <-x[1]
  pamax<-x[2]
  
  pbmax<- pamax+tmax
  
  QED <- ED(pbmax)
  QES <- ES(pamax)
  
  res<-QED-QES #equilibrium
  
  return(res)
 
}

optimumW <- auglag(par=c(1,1),fn=fmaxW,heq=EqconstW)

# Response ------------------------
TmaxW = optimumW$par[1];TmaxW
pamaxW = optimumW$par[2];pamaxW
pbmaxW = pamaxW + TmaxW;pbmaxW
TR.MaxW = TmaxW*ED(pbmaxW);TR.MaxW

round(data.frame(abs(ees),abs(eed),TmaxR,TR.MaxR,TmaxW,TR.MaxW),3)
