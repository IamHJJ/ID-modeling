## Replication: The transmission dynamics of gonorrhoea: modelling the reported behaviour of infected patients from Newark, New Jersey

library(deSolve)

## SIS model
# X = susceptible
# Y = symptomatic infected
# A = asymptomatic infected

SIS2riskGrs <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    #Males
    dXM1 = mu*NM*phiM1 - mu*XM1 - XM1*betaM*cM1*rhoM1*((YF+AF)/NF) + sigmaM*YM1 + gammaM*AM1;
    dYM1 = thetaM*XM1*betaM*cM1*rhoM*((YF+AF)/NF)-(mu+sigmaM)*YM1; 
    dAM1 = (1-sigma)*XM1*betaM*cM1*rhoM1*((YF+AF)/NF)-(mu+gammaM)*AM1;
    
    dXM2 = mu*NM*phiM2 - mu*XM2 - XM2*betaM*cM2*rhoM2*((YF+AF)/NF) + sigmaM*YM2 + gammaM*AM2;
    dYM2 = thetaM*XM2*betaM*cM2*rhoM*((YF+AF)/NF)-(mu+sigmaM)*YM2; 
    dAM2 = (1-sigma)*XM2*betaM*cM2*rhoM2*((YF+AF)/NF)-(mu+gammaM)*AM2;
    
    dXM3 = mu*NM*phiM3 - mu*XM3 - XM3*betaM*cM3*rhoM3*((YF+AF)/NF) + sigmaM*YM3 + gammaM*AM3;
    dYM3 = thetaM*XM3*betaM*cM3*rhoM*((YF+AF)/NF)-(mu+sigmaM)*YM3; 
    dAM3 = (1-sigma)*XM3*betaM*cM3*rhoM3*((YF+AF)/NF)-(mu+gammaM)*AM3;
    
    dXM4 = mu*NM*phiM4 - mu*XM4 - XM4*betaM*cM4*rhoM4*((YF+AF)/NF) + sigmaM*YM4 + gammaM*AM4;
    dYM4 = thetaM*XM4*betaM*cM4*rhoM*((YF+AF)/NF)-(mu+sigmaM)*YM4; 
    dAM4 = (1-sigma)*XM4*betaM*cM4*rhoM4*((YF+AF)/NF)-(mu+gammaM)*AM4;
    
    #Females
    dXF1 = mu*NF*phiF1 - mu*XF1 - XF1*betaF*cF1*rhoF1*((YM+AM)/NM) + sigmaF*YF1 + gammaF*AF1;
    dYF1 = thetaF*XF1*betaF*cF1*rhoF*((YM+AM)/NM)-(mu+sigmaF)*YF1; 
    dAF1 = (1-sigma)*XF1*betaF*cF1*rhoF1*((YM+AM)/NM)-(mu+gammaF)*AF1;
    
    dXF2 = mu*NF*phiF2 - mu*XF2 - XF2*betaF*cF2*rhoF2*((YM+AM)/NM) + sigmaF*YF2 + gammaF*AF2;
    dYF2 = thetaF*XF2*betaF*cF2*rhoF*((YM+AM)/NM)-(mu+sigmaF)*YF2; 
    dAF2 = (1-sigma)*XF2*betaF*cF2*rhoF2*((YM+AM)/NM)-(mu+gammaF)*AF2;
    
    dXM3 = mu*NF*phiF3 - mu*XF3 - XF3*betaF*cF3*rhoF3*((YM+AM)/NM) + sigmaF*YF3 + gammaF*AF3;
    dYM3 = thetaF*XF3*betaF*cF3*rhoF*((YM+AM)/NM)-(mu+sigmaF)*YF3; 
    dAM3 = (1-sigma)*XF3*betaF*cF3*rhoF3*((YM+AM)/NM)-(mu+gammaF)*AF3;
    
    dXM4 = mu*NF*phiF4 - mu*XF4 - XF4*betaF*cF4*rhoF4*((YM+AM)/NM) + sigmaF*YF4 + gammaF*AF4;
    dYM4 = thetaF*XF4*betaF*cF4*rhoF*((YM+AM)/NM)-(mu+sigmaF)*YF4; 
    dAM4 = (1-sigma)*XF4*betaF*cF4*rhoF4*((YM+AM)/NM)-(mu+gammaF)*AF4;
   
    return(list(c(dXM1, dYM1, dAM1, dXM2, dYM2, dAM2, dXM3, dYM3, dAM3,
                  dXM4, dYM4, dAM4, dXF1, dYF1, dAF1, dXF2, dYF2, dAF2,
                  dXF3, dYF3, dAF3, dXF4, dYF4, dAF4)))
    
  })
}

#initial condition
mu = 0.2 #net birth and death rate
phiM1=0.0055 #proportion in each activity group
phiM2=0.0445
phiM3=0.2
phiM4=0.75
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
gamma = 1/185 #recovery rate for asympt. men and women
betaF = 0.8 #transmission probability men to women
betaM = 0.6 #transmission probability women to men
thetaM = 0.05 #proportion of asympt. men
thetaF = 0.4 #proportion of asympt. women
epsilon = 0.8 #pattern of mixing

NM = 1000 #?
NF = 1000 #?

  
XM1, YM1, AM1, XM2, YM2, AM2, XM3, YM3, AM3,
XM4, YM4, AM4, XF1, YF1, AF1, XF2, YF2, AF2,
XF3, YF3, AF3, XF4, YF4, AF4
  

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

parameters=c(mu=mu, rhoM=rhoM, rhoF=rhoF, cM=cM, cF=cF, sigmaM=sigmaM, 
            sigmaF=sigmaF, gamma=gamma, betaF=betaF, betaM=betaM, thetaM=thetaM,
            thetaF=thetaF);

times=seq(1,1000*365,by=1)

#Simulation model
sim=ode(y=state,times=times,func=SIS2riskGrs,parms=parameters)
sim=as.data.frame(sim)




