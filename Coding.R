## Replication: The transmission dynamics of gonorrhoea: modelling the reported behaviour of infected patients from Newark, New Jersey

library(deSolve)

## SIS model
# X = susceptible
# Y = symptomatic infected
# A = asymptomatic infected

################################################################

SIS2riskGrs <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    #Males
    dXM1 = mu*NM*phiM1 - mu*XM1 - XM1*betaM*cM1*sum(rhoM[1,]*((c(YF1, YF2, YF3, YF4)+c(AF1, AF2, AF3, AF4))/c(NF1,NF2,NF3,NF4))) + sigmaM*YM1 + gammaM*AM1;
    dYM1 = (1-thetaM)*XM1*betaM*cM1*sum(rhoM[1,]*((c(YF1, YF2, YF3, YF4)+c(AF1, AF2, AF3, AF4))/c(NF1,NF2,NF3,NF4)))-(mu+sigmaM)*YM1; 
    dAM1 = (thetaM)*XM1*betaM*cM1*sum(rhoM[1,]*((c(YF1, YF2, YF3, YF4)+c(AF1, AF2, AF3, AF4))/c(NF1,NF2,NF3,NF4)))-(mu+gammaM)*AM1;
    
    dXM2 = mu*NM*phiM2 - mu*XM2 - XM2*betaM*cM2*sum(rhoM[2,]*((c(YF1, YF2, YF3, YF4)+c(AF1, AF2, AF3, AF4))/c(NF1,NF2,NF3,NF4))) + sigmaM*YM2 + gammaM*AM2;
    dYM2 = (1-thetaM)*XM2*betaM*cM2*sum(rhoM[2,]*((c(YF1, YF2, YF3, YF4)+c(AF1, AF2, AF3, AF4))/c(NF1,NF2,NF3,NF4)))-(mu+sigmaM)*YM2; 
    dAM2 = (thetaM)*XM2*betaM*cM2*sum(rhoM[2,]*((c(YF1, YF2, YF3, YF4)+c(AF1, AF2, AF3, AF4))/c(NF1,NF2,NF3,NF4)))-(mu+gammaM)*AM2;
    
    dXM3 = mu*NM*phiM3 - mu*XM3 - XM3*betaM*cM3*sum(rhoM[3,]*((c(YF1, YF2, YF3, YF4)+c(AF1, AF2, AF3, AF4))/c(NF1,NF2,NF3,NF4))) + sigmaM*YM3 + gammaM*AM3;
    dYM3 = (1-thetaM)*XM3*betaM*cM3*sum(rhoM[3,]*((c(YF1, YF2, YF3, YF4)+c(AF1, AF2, AF3, AF4))/c(NF1,NF2,NF3,NF4)))-(mu+sigmaM)*YM3; 
    dAM3 = (thetaM)*XM3*betaM*cM3*sum(rhoM[3,]*((c(YF1, YF2, YF3, YF4)+c(AF1, AF2, AF3, AF4))/c(NF1,NF2,NF3,NF4)))-(mu+gammaM)*AM3;
    
    dXM4 = mu*NM*phiM4 - mu*XM4 - XM4*betaM*cM4*sum(rhoM[4,]*((c(YF1, YF2, YF3, YF4)+c(AF1, AF2, AF3, AF4))/c(NF1,NF2,NF3,NF4))) + sigmaM*YM4 + gammaM*AM4;
    dYM4 = (1-thetaM)*XM4*betaM*cM4*sum(rhoM[4,]*((c(YF1, YF2, YF3, YF4)+c(AF1, AF2, AF3, AF4))/c(NF1,NF2,NF3,NF4)))-(mu+sigmaM)*YM4; 
    dAM4 = (thetaM)*XM4*betaM*cM4*sum(rhoM[4,]*((c(YF1, YF2, YF3, YF4)+c(AF1, AF2, AF3, AF4))/c(NF1,NF2,NF3,NF4)))-(mu+gammaM)*AM4;
    
    #Females
    dXF1 = mu*NF*phiF1 - mu*XF1 - XF1*betaF*cF1*sum(rhoF[1,]*((c(YM1, YM2, YM3, YM4)+c(AM1, AM2, AM3, AM4))/c(NM1,NM2,NM3,NM4))) + sigmaF*YF1 + gammaF*AF1;
    dYF1 = (1-thetaF)*XF1*betaF*cF1*sum(rhoF[1,]*((c(YM1, YM2, YM3, YM4)+c(AM1, AM2, AM3, AM4))/c(NM1,NM2,NM3,NM4)))-(mu+sigmaF)*YF1; 
    dAF1 = (thetaF)*XF1*betaF*cF1*sum(rhoF[1,]*((c(YM1, YM2, YM3, YM4)+c(AM1, AM2, AM3, AM4))/c(NM1,NM2,NM3,NM4)))-(mu+gammaF)*AF1;
    
    dXF2 = mu*NF*phiF2 - mu*XF2 - XF2*betaF*cF2*sum(rhoF[2,]*((c(YM1, YM2, YM3, YM4)+c(AM1, AM2, AM3, AM4))/c(NM1,NM2,NM3,NM4))) + sigmaF*YF2 + gammaF*AF2;
    dYF2 = (1-thetaF)*XF2*betaF*cF2*sum(rhoF[2,]*((c(YM1, YM2, YM3, YM4)+c(AM1, AM2, AM3, AM4))/c(NM1,NM2,NM3,NM4)))-(mu+sigmaF)*YF2; 
    dAF2 = (thetaF)*XF2*betaF*cF2*sum(rhoF[2,]*((c(YM1, YM2, YM3, YM4)+c(AM1, AM2, AM3, AM4))/c(NM1,NM2,NM3,NM4)))-(mu+gammaF)*AF2;
    
    dXF3 = mu*NF*phiF3 - mu*XF3 - XF3*betaF*cF3*sum(rhoF[3,]*(c(YM1, YM2, YM3, YM4)+c(AM1, AM2, AM3, AM4))/c(NM1,NM2,NM3,NM4)) + sigmaF*YF3 + gammaF*AF3;
    dYF3 = (1-thetaF)*XF3*betaF*cF3*sum(rhoF[3,]*(c(YM1, YM2, YM3, YM4)+c(AM1, AM2, AM3, AM4))/c(NM1,NM2,NM3,NM4))-(mu+sigmaF)*YF3; 
    dAF3 = (thetaF)*XF3*betaF*cF3*sum(rhoF[3,]*(c(YM1, YM2, YM3, YM4)+c(AM1, AM2, AM3, AM4))/c(NM1,NM2,NM3,NM4))-(mu+gammaF)*AF3;
    
    dXF4 = mu*NF*phiF4 - mu*XF4 - XF4*betaF*cF4*sum(rhoF[4,]*(c(YM1, YM2, YM3, YM4)+c(AM1, AM2, AM3, AM4))/c(NM1,NM2,NM3,NM4)) + sigmaF*YF4 + gammaF*AF4;
    dYF4 = (1-thetaF)*XF4*betaF*cF4*sum(rhoF[4,]*(c(YM1, YM2, YM3, YM4)+c(AM1, AM2, AM3, AM4))/c(NM1,NM2,NM3,NM4))-(mu+sigmaF)*YF4; 
    dAF4 = (thetaF)*XF4*betaF*cF4*sum(rhoF[4,]*(c(YM1, YM2, YM3, YM4)+c(AM1, AM2, AM3, AM4))/c(NM1,NM2,NM3,NM4))-(mu+gammaF)*AF4;
   
    return(list(c(dXM1, dYM1, dAM1, dXM2, dYM2, dAM2, dXM3, dYM3, dAM3,
                  dXM4, dYM4, dAM4, dXF1, dYF1, dAF1, dXF2, dYF2, dAF2,
                  dXF3, dYF3, dAF3, dXF4, dYF4, dAF4)))
    
  })
}

#initial condition
mu = 0.2 #net birth and death rate
phiM=c(0.0055,0.0445,0.2,0.75) #proportion in each activity group
phiM1=0.0055 
phiM2=0.0445
phiM3=0.2
phiM4=0.75
phiF=c(0.0055,0.0145,0.1,0.88)
phiF1=0.0055
phiF2=0.0145
phiF3=0.1
phiF4=0.88
cM=c(40,10,2.4,0.5) #mean rate of sex-partner change
cM1=40 
cM2=10
cM3=2.4
cM4=0.5
cF=c(40,30.7,4.8,0.43)
cF1=40
cF2=30.7
cF3=4.8
cF4=0.43
sigmaM = 1/13 #recovery rate for sympt. men
sigmaF = 1/20 #recovery rate for sympt. women
gammaF = 1/185 #recovery rate for asympt. men and women
gammaM = 1/185 #recovery rate for asympt. men and women
betaF = 0.8 #transmission probability men to women
betaM = 0.6 #transmission probability women to men
thetaM = 0.05 #proportion of asympt. men
thetaF = 0.4 #proportion of asympt. women
epsilon = 0.8 #pattern of mixing

NM = 10000
NF = 10000
NM1 = 55
XM1 = 35
YM1 = 10
AM1 = 10
NM2 = 445
XM2 = 425
YM2 = 10
AM2 = 10
NM3 = 2000
XM3 = 1980
YM3 = 10
AM3 = 10
NM4 = 7500
XM4 = 7480
YM4 = 10
AM4 = 10
NF1 = 55
XF1 = 35
YF1 = 10
AF1 = 10
NF2 = 445
XF2 = 425
YF2 = 10
AF2 = 10
NF3 = 2000
XF3 = 1980
YF3 = 10
AF3 = 10
NF4 = 7500
XF4 = 7480
YF4 = 10
AF4 = 10


   
  

#######################################################################
library(matrixcalc)

# Define the function to calculate the rho matrix for male
rhoM_matrix <- matrix(0, nrow = 4, ncol = 4)

for (i in 1:4) {
    for (j in 1:4) {
      if (i == j) {
        rhoM_matrix[i, j] <- epsilon + (1 - epsilon) * ((cF[i]*NF)/sum(cF*NF))
      } else {
        rhoM_matrix[i, j] <- (1 - epsilon) * ((cF[i]*NF)/sum(cF*NF))
      }
    }
}

# Define the function to calculate the rho matrix for female
rhoF_matrix <- matrix(0, nrow = 4, ncol = 4)

for (i in 1:4) {
  for (j in 1:4) {
    if (i == j) {
      rhoF_matrix[i, j] <- epsilon + (1 - epsilon) * ((cM[i]*NM)/sum(cM*NM))
    } else {
      rhoF_matrix[i, j] <- (1 - epsilon) * ((cM[i]*NM)/sum(cM*NM))
    }
  }
}

#rhoF
calculate_rhoF_matrix <- function(epsilon, cM, NM) {
  rhoF_matrix <- matrix(0, nrow = 4, ncol = 4)
  for (i in 1:4) {
    for (j in 1:4) {
      if (i == j) {
        rho_matrix[i, j] <- epsilon + (1 - epsilon) * ((cM*NM)/sum(cM*NM))
      } else {
        rho_matrix[i, j] <- (1 - epsilon) * ((cM*NM)/sum(cM*NM))
      }
    }
  }
  return(rhoF_matrix)
}


################################################################################
state=c(XM1, YM1, AM1, XM2, YM2, AM2, XM3, YM3, AM3,
        XM4, YM4, AM4, XF1, YF1, AF1, XF2, YF2, AF2,
        XF3, YF3, AF3, XF4, YF4, AF4);

parameters=list(mu=mu, rhoM=rhoM_matrix, rhoF=rhoF_matrix, cM=cM, cF=cF, sigmaM=sigmaM, 
            sigmaF=sigmaF, gamma=gamma, betaF=betaF, betaM=betaM, thetaM=thetaM,
            thetaF=thetaF, NF1 = NF1, NF2 = NF2, NF3 = NF3, NF4 = NF4, NM1 = NM1, 
            NM2 = NM2, NM3 = NM3, NM4 = NM4);

times=seq(1,1000*365,by=1);

#Simulation model
sim=ode(y=state,times=times,func=SIS2riskGrs,parms=parameters)
sim=as.data.frame(sim)





