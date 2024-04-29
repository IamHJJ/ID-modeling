## Replication: The transmission dynamics of gonorrhoea: modelling the reported behaviour of infected patients from Newark, New Jersey

library(deSolve)
library(tidyverse)

## SIS model
# X = susceptible
# Y = symptomatic infected
# A = asymptomatic infected

#hi
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
mu = 0.2 / 365 #net birth and death rate
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
cM=c(40,10,2.4,0.5)/365 #mean rate of sex-partner change
cM1=40/365
cM2=10/365
cM3=2.4/365
cM4=0.5/365
cF=c(40,30.7,4.8,0.43)/365
cF1=40/365
cF2=30.7/365
cF3=4.8/365
cF4=0.43/365
sigmaM = 1/13 #recovery rate for sympt. men
sigmaF = 1/20 #recovery rate for sympt. women
sigmaM2 = 10/117 #recovery rate for sympt. men after intervention
sigmaF2 = 1/18 #recovery rate for sympt. women after intervention
gammaF = 1/185 #recovery rate for asympt. men and women
gammaM = 1/185 #recovery rate for asympt. men and women
gammaM2 = 10/1665
gammaF2 = 10/1665 #recovery rate for asympt. men and women after intervention
betaM = 0.8 #transmission probability men to women
betaF = 0.6 #transmission probability women to men
betaM2 = 0.72 #transmission probability men to women after intervention
betaF2 = 0.54 #transmission probability women to men after intervention
thetaM = 0.05 #proportion of asympt. men
thetaF = 0.4 #proportion of asympt. women
epsilon = 0.8 #pattern of mixing

NM = 1
NF = 1

NM1 = NM * phiM1
YM1 = 0.00275*0.007
AM1 = 0.00275*0.007*0.05
XM1 = NM1-AM1-YM1

NM2 = NM * phiM2
YM2 = 0.02225*0.007
AM2 = 0.02225*0.007*0.05
XM2 = NM2-AM2-YM2

NM3 = NM * phiM3
YM3 = 0.1*0.007
AM3 = 0.1*0.007*0.05
XM3 = NM3-AM3-YM3

NM4 = NM * phiM4
YM4 = 0.375*0.007
AM4 = 0.375*0.007*0.05
XM4 = NM4-AM4-YM4

NF1 = NF * phiF1
YF1 = 0.00275*0.007
AF1 = 0.00275*0.007*0.4
XF1 = NF1-AF1-YF1

NF2 = NF * phiF2
YF2 = 0.00725*0.007
AF2 = 0.00725*0.007*0.4
XF2 = NF2-AF2-YF2

NF3 = NF * phiF3
YF3 = 0.05*0.007
AF3 = 0.05*0.007*0.4
XF3 = NF3-AF3-YF3

NF4 = NF * phiF4
YF4 = 0.44*0.007
AF4 = 0.44*0.007*0.4
XF4 = NF4-AF4-YF4

NM_matrix=c(NM1,NM2,NM3,NM4)
NF_matrix=c(NF1,NF2,NF3,NF4)


#######################################################################
library(matrixcalc)

# Define the function to calculate the rho matrix for male
rhoM_matrix <- matrix(0, nrow = 4, ncol = 4)

for (i in 1:4) {
  for (j in 1:4) {
    if (i == j) {
      rhoM_matrix[i, j] <- epsilon + (1 - epsilon) * ((cF[j]*NF_matrix[j])/sum(cF*NF_matrix))
    } else {
      rhoM_matrix[i, j] <- (1 - epsilon) * ((cF[j]*NF_matrix[j])/sum(cF*NF_matrix))
    }
  }
}


# Define the function to calculate the rho matrix for female
rhoF_matrix <- matrix(0, nrow = 4, ncol = 4)

for (i in 1:4) {
  for (j in 1:4) {
    if (i == j) {
      rhoF_matrix[i, j] <- epsilon + (1 - epsilon) * ((cM[j]*NM_matrix[j])/sum(cM*NM_matrix))
    } else {
      rhoF_matrix[i, j] <- (1 - epsilon) * ((cM[j]*NM_matrix[j])/sum(cM*NM_matrix))
    }
  }
}


################################################################################
state=c(XM1=XM1, YM1=YM1, AM1=AM1, XM2=XM2, YM2=YM2, AM2=AM2, XM3=XM3, YM3=YM3, AM3=AM3,
        XM4=XM4, YM4=YM4, AM4=AM4, XF1=XF1, YF1=YF1, AF1=AF1, XF2=XF2, YF2=YF2, AF2=AF2,
        XF3=XF3, YF3=YF3, AF3=AF3, XF4=XF4, YF4=YF4, AF4=AF4);

parameters=list(mu=mu, rhoM=rhoM_matrix, rhoF=rhoF_matrix, cM=cM, cF=cF, sigmaM=sigmaM, 
                 sigmaF=sigmaF, gammaF=gammaF, gammaM = gammaM, betaF=betaF, betaM=betaM, thetaM=thetaM,
                 thetaF=thetaF, NF1 = NF1, NF2 = NF2, NF3 = NF3, NF4 = NF4, NM1 = NM1, 
                 NM2 = NM2, NM3 = NM3, NM4 = NM4);

times=seq(1,10*365,by=1);

#Simulation model
sim=ode(y=state,times=times,func=SIS2riskGrs,parms=parameters) %>% data.frame
sim_stat=
  as.data.frame(sim)|>
  mutate(IM1 = YM1+AM1,
         IM2 = YM2+AM2,
         IM3 = YM3+AM3,
         IM4 = YM4+AM4,
         IF1 = YF1+AF1,
         IF2 = YF2+AF2,
         IF3 = YF3+AF3,
         IF4 = YF4+AF4,
         IM = IM1+IM2+IM3+IM4,
         IF = IF1+IF1+IF3+IF4,
         I = IM+IF,
         AM=AM1+AM2+AM3+AM4,
         YM=YM1+YM2+YM3+YM4)

# Prevalence: Overall + Male + Female
plot(x=times, y=sim_stat$YM, xlab = "time", ylab = "Prevalence", type = "l", main = "Prevalence", col = "black");
lines(sim$IM[1:40], type = "l", col = "blue");
lines(sim$IF[1:40], type = "l", col = "red");
legend("topright", legend=c("Overall","Male","Female"), col=c("black","blue","red"), lty=1)

plot(x=times, y=sim$I, xlab = "time", ylab = "Prevalence", type = "l", main = "Prevalence", col = "black");
lines(sim$IM, type = "l", col = "blue");
lines(sim$IF, type = "l", col = "red");
legend("topright", legend=c("Overall","Male","Female"), col=c("black","blue","red"), lty=1)

# Prevalence: Symptomatic Male in each group
plot(x=times[1:40], y=sim$YM1[1:40], xlab = "time", ylab = "Prevalence", type = "l",ylim=c(0, 3.0e-3), main = "Prevalence of Symptomatic Infections among Males", col = "red");
lines(sim$YM2[1:40], type = "l", col = "green");
lines(sim$YM3[1:40], type = "l", col = "blue");
lines(sim$YM4[1:40], type = "l", col = "orange");
legend("topright", legend=c("Group 1", "Group 2", "Group 3", "Group 4"), col=c("red", "green","blue","orange"), lty=1)

# Prevalence: Symptomatic Female
plot(x=times[1:40], y=sim$YF1[1:40], xlab = "time", ylab = "Prevalence", type = "l",ylim=c(0, 3.0e-3), main = "Prevalence of Symptomatic Infections among Females", col = "red");
lines(sim$YF2[1:40], type = "l", col = "green");
lines(sim$YF3[1:40], type = "l", col = "blue");
lines(sim$YF4[1:40], type = "l", col = "orange");
legend("topright", legend=c("Group 1", "Group 2", "Group 3", "Group 4"), col=c("red", "green","blue","orange"), lty=1)


par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(1.5,.5,0),cex=.9)
# Prevalence: Asymptomatic Male
plot(x=times[1:40], y=sim$AM1[1:40], xlab = "time", ylab = "Prevalence", type = "l", ylim=c(0, 1e-3),main = "Prevalence of Asymptomatic Infections among Males", col = "red");
lines(sim$AM2[1:40], type = "l", col = "green");
lines(sim$AM3[1:40], type = "l", col = "blue");
lines(sim$AM4[1:40], type = "l", col = "orange");
legend("topright", legend=c("Group 1", "Group 2", "Group 3", "Group 4"), col=c("red", "green","blue","orange"), lty=1)


# Prevalence: Asymptomatic Female
plot(x=times[1:40], y=sim$AF1[1:40], xlab = "time", ylab = "Prevalence", type = "l", ylim=c(0, 1e-3),main = "Prevalence of Asymptomatic Infections among Females", col = "red");
lines(sim$AF2[1:40], type = "l", col = "green");
lines(sim$AF3[1:40], type = "l", col = "blue");
lines(sim$AF4[1:40], type = "l", col = "orange");
legend("topright", legend=c("Group 1", "Group 2", "Group 3", "Group 4"), col=c("red", "green","blue","orange"), lty=1)



#parameters for changes in sigma for male and female
parameters1=list(mu=mu, rhoM=rhoM_matrix, rhoF=rhoF_matrix, cM=cM, cF=cF, sigmaM=sigmaM2, 
                 sigmaF=sigmaF2, gammaF=gammaF, gammaM = gammaM, betaF=betaF, betaM=betaM, thetaM=thetaM,
                 thetaF=thetaF, NF1 = NF1, NF2 = NF2, NF3 = NF3, NF4 = NF4, NM1 = NM1, 
                 NM2 = NM2, NM3 = NM3, NM4 = NM4);
#Simulation model
sim1=ode(y=state,times=times,func=SIS2riskGrs,parms=parameters1)
sim1=
  as.data.frame(sim1)|>
  mutate(IM1 = YM1+AM1,
         IM2 = YM2+AM2,
         IM3 = YM3+AM3,
         IM4 = YM4+AM4,
         IF1 = YF1+AF1,
         IF2 = YF2+AF2,
         IF3 = YF3+AF3,
         IF4 = YF4+AF4,
         IM = IM1+IM2+IM3+IM4,
         IF = IF1+IF1+IF3+IF4,
         I = IM+IF,
         AM=AM1+AM2+AM3+AM4,
         YM=YM1+YM2+YM3+YM4)

par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(1.5,.5,0),cex=.9)

plot(x=times[1:40], y=sim1$I[1:40]-sim$I[1:40], xlab = "time", ylab = "change in Prevalence", ylim=c(-8.0e-5, 0), 
     type = "l", main = "Change in Prevalence after intervention (sigma)", col = "black");
lines(sim1$IM[1:40]-sim$IM[1:40], type = "l", col = "blue");
lines(sim1$IF[1:40]-sim$IF[1:40], type = "l", col = "red");
legend("topright", legend=c("Overall","Male","Female"), col=c("black","blue","red"), lty=1)

# Prevalence: Symptomatic Male in each group
plot(x=times[1:40], y=sim1$YM1[1:40]-sim$YM1[1:40], xlab = "time", ylab = "change in Prevalence", type = "l",ylim=c(-4.0e-5, 0), main = "Prevalence of Symptomatic Infections among Males (sigma)", col = "red");
lines(sim1$YM2[1:40]-sim$YM2[1:40], type = "l", col = "green");
lines(sim1$YM3[1:40]-sim$YM3[1:40], type = "l", col = "blue");
lines(sim1$YM4[1:40]-sim$YM4[1:40], type = "l", col = "orange");
legend("topright", legend=c("Group 1", "Group 2", "Group 3", "Group 4"), col=c("red", "green","blue","orange"), lty=1)

# Prevalence: Symptomatic Female
plot(x=times[1:40], y=sim1$YF1[1:40]-sim$YF1[1:40], xlab = "time", ylab = "change in Prevalence", type = "l",ylim=c(-4.0e-5, 0), main = "Prevalence of Symptomatic Infections among Females (sigma)", col = "red");
lines(sim1$YF2[1:40]-sim$YF2[1:40], type = "l", col = "green");
lines(sim1$YF3[1:40]-sim$YF3[1:40], type = "l", col = "blue");
lines(sim1$YF4[1:40]-sim$YF4[1:40], type = "l", col = "orange");
legend("topright", legend=c("Group 1", "Group 2", "Group 3", "Group 4"), col=c("red", "green","blue","orange"), lty=1)

# Prevalence: Asymptomatic Male
plot(x=times[1:40], y=sim1$AM1[1:40]-sim$AM1[1:40], xlab = "time", ylab = "change in Prevalence", type = "l", ylim=c(-4.0e-8, 0),
     main = "Prevalence of Asymptomatic Infections among Males (sigma)", col = "red");
lines(sim1$AM2[1:40]-sim$AM2[1:40], type = "l", col = "green");
lines(sim1$AM3[1:40]-sim$AM3[1:40], type = "l", col = "blue");
lines(sim1$AM4[1:40]-sim$AM4[1:40], type = "l", col = "orange");
legend("topright", legend=c("Group 1", "Group 2", "Group 3", "Group 4"), col=c("red", "green","blue","orange"), lty=1)

# Prevalence: Asymptomatic Female
plot(x=times[1:40], y=sim1$AF1[1:40]-sim$AF1[1:40], xlab = "time", ylab = "change in Prevalence", type = "l", ylim=c(-4.0e-8, 0),
     main = "Prevalence of Asymptomatic Infections among Females (sigma)", col = "red");
lines(sim1$AF2[1:40]-sim$AF2[1:40], type = "l", col = "green");
lines(sim1$AF3[1:40]-sim$AF3[1:40], type = "l", col = "blue");
lines(sim1$AF4[1:40]-sim$AF4[1:40], type = "l", col = "orange");
legend("topright", legend=c("Group 1", "Group 2", "Group 3", "Group 4"), col=c("red", "green","blue","orange"), lty=1)

#parameters for changes in gamma for male and female
parameters2=list(mu=mu, rhoM=rhoM_matrix, rhoF=rhoF_matrix, cM=cM, cF=cF, sigmaM=sigmaM, 
                 sigmaF=sigmaF, gammaF=gammaF2, gammaM = gammaF2, betaF=betaF, betaM=betaM, thetaM=thetaM,
                 thetaF=thetaF, NF1 = NF1, NF2 = NF2, NF3 = NF3, NF4 = NF4, NM1 = NM1, 
                 NM2 = NM2, NM3 = NM3, NM4 = NM4);
#Simulation model
sim2=ode(y=state,times=times,func=SIS2riskGrs,parms=parameters2)
sim2=
  as.data.frame(sim2)|>
  mutate(IM1 = YM1+AM1,
         IM2 = YM2+AM2,
         IM3 = YM3+AM3,
         IM4 = YM4+AM4,
         IF1 = YF1+AF1,
         IF2 = YF2+AF2,
         IF3 = YF3+AF3,
         IF4 = YF4+AF4,
         IM = IM1+IM2+IM3+IM4,
         IF = IF1+IF1+IF3+IF4,
         I = IM+IF,
         AM=AM1+AM2+AM3+AM4,
         YM=YM1+YM2+YM3+YM4)

plot(x=times[1:40], y=sim2$I[1:40]-sim$I[1:40], xlab = "time", ylab = "change in Prevalence", ylim=c(-2e-6, 0), 
     type = "l", main = "Change in Prevalence after intervention (gamma)", col = "black");
lines(sim2$IM[1:40]-sim$IM[1:40], type = "l", col = "blue");
lines(sim2$IF[1:40]-sim$IF[1:40], type = "l", col = "red");
legend("topright", legend=c("Overall","Male","Female"), col=c("black","blue","red"), lty=1)

par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(1.5,.5,0),cex=.9)
# Prevalence: Symptomatic Male in each group
plot(x=times[1:40], y=sim2$YM1[1:40]-sim$YM1[1:40], xlab = "time", ylab = "change in Prevalence", type = "l",
     ylim=c(-5.0e-9, 0), main = "Prevalence of Symptomatic Infections among Males (gamma)", col = "red");
lines(sim2$YM2[1:40]-sim$YM2[1:40], type = "l", col = "green");
lines(sim2$YM3[1:40]-sim$YM3[1:40], type = "l", col = "blue");
lines(sim2$YM4[1:40]-sim$YM4[1:40], type = "l", col = "orange");
legend("topright", legend=c("Group 1", "Group 2", "Group 3", "Group 4"), col=c("red", "green","blue","orange"), lty=1)

# Prevalence: Symptomatic Female
plot(x=times[1:40], y=sim2$YF1[1:40]-sim$YF1[1:40], xlab = "time", ylab = "change in Prevalence", type = "l",
     ylim=c(-5.0e-9, 0), main = "Prevalence of Symptomatic Infections among Females (gamma)", col = "red");
lines(sim2$YF2[1:40]-sim$YF2[1:40], type = "l", col = "green");
lines(sim2$YF3[1:40]-sim$YF3[1:40], type = "l", col = "blue");
lines(sim2$YF4[1:40]-sim$YF4[1:40], type = "l", col = "orange");
legend("topright", legend=c("Group 1", "Group 2", "Group 3", "Group 4"), col=c("red", "green","blue","orange"), lty=1)

# Prevalence: Asymptomatic Male
plot(x=times[1:40], y=sim2$AM1[1:40]-sim$AM1[1:40], xlab = "time", ylab = "change in Prevalence", type = "l", ylim=c(-1.5e-6, 0),
     main = "Prevalence of Asymptomatic Infections among Males (gamma)", col = "red");
lines(sim2$AM2[1:40]-sim$AM2[1:40], type = "l", col = "green");
lines(sim2$AM3[1:40]-sim$AM3[1:40], type = "l", col = "blue");
lines(sim2$AM4[1:40]-sim$AM4[1:40], type = "l", col = "orange");
legend("topright", legend=c("Group 1", "Group 2", "Group 3", "Group 4"), col=c("red", "green","blue","orange"), lty=1)

# Prevalence: Asymptomatic Female
plot(x=times[1:40], y=sim2$AF1[1:40]-sim$AF1[1:40], xlab = "time", ylab = "change in Prevalence", type = "l", ylim=c(-1.5e-6, 0),
     main = "Prevalence of Asymptomatic Infections among Females (gamma)", col = "red");
lines(sim2$AF2[1:40]-sim$AF2[1:40], type = "l", col = "green");
lines(sim2$AF3[1:40]-sim$AF3[1:40], type = "l", col = "blue");
lines(sim2$AF4[1:40]-sim$AF4[1:40], type = "l", col = "orange");
legend("topright", legend=c("Group 1", "Group 2", "Group 3", "Group 4"), col=c("red", "green","blue","orange"), lty=1)

#parameters for changes in beta for male and female
parameters3=list(mu=mu, rhoM=rhoM_matrix, rhoF=rhoF_matrix, cM=cM, cF=cF, sigmaM=sigmaM, 
                 sigmaF=sigmaF, gammaF=gammaF, gammaM = gammaM, betaF=betaF2, betaM=betaM2, thetaM=thetaM,
                 thetaF=thetaF, NF1 = NF1, NF2 = NF2, NF3 = NF3, NF4 = NF4, NM1 = NM1, 
                 NM2 = NM2, NM3 = NM3, NM4 = NM4);
#Simulation model
sim3=ode(y=state,times=times,func=SIS2riskGrs,parms=parameters3)
sim3=
  as.data.frame(sim3)|>
  mutate(IM1 = YM1+AM1,
         IM2 = YM2+AM2,
         IM3 = YM3+AM3,
         IM4 = YM4+AM4,
         IF1 = YF1+AF1,
         IF2 = YF2+AF2,
         IF3 = YF3+AF3,
         IF4 = YF4+AF4,
         IM = IM1+IM2+IM3+IM4,
         IF = IF1+IF1+IF3+IF4,
         I = IM+IF,
         AM=AM1+AM2+AM3+AM4,
         YM=YM1+YM2+YM3+YM4)

plot(x=times[1:40], y=sim3$I[1:40]-sim$I[1:40], xlab = "time", ylab = "change in Prevalence", ylim=c(-4.5e-6, 0), 
     type = "l", main = "Change in Prevalence after intervention (beta)", col = "black");
lines(sim3$IM[1:40]-sim$IM[1:40], type = "l", col = "blue");
lines(sim3$IF[1:40]-sim$IF[1:40], type = "l", col = "red");
legend("topright", legend=c("Overall","Male","Female"), col=c("black","blue","red"), lty=1)


par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(1.5,.5,0),cex=.9)
# Prevalence: Symptomatic Male in each group after change in beta
plot(x=times[1:40], y=sim3$YM1[1:40]-sim$YM1[1:40], xlab = "time", ylab = "change in Prevalence", type = "l",
     ylim=c(-8.0e-7, 0), main = "Prevalence of Symptomatic Infections among Males (beta)", col = "red");
lines(sim3$YM2[1:40]-sim$YM2[1:40], type = "l", col = "green");
lines(sim3$YM3[1:40]-sim$YM3[1:40], type = "l", col = "blue");
lines(sim3$YM4[1:40]-sim$YM4[1:40], type = "l", col = "orange");
legend("topright", legend=c("Group 1", "Group 2", "Group 3", "Group 4"), col=c("red", "green","blue","orange"), lty=1)
# Prevalence: Symptomatic Female in each group after change in beta
plot(x=times[1:40], y=sim3$YF1[1:40]-sim$YF1[1:40], xlab = "time", ylab = "change in Prevalence", type = "l",
     ylim=c(-8.0e-7, 0), main = "Prevalence of Symptomatic Infections among Males (beta)", col = "red");
lines(sim3$YF2[1:40]-sim$YF2[1:40], type = "l", col = "green");
lines(sim3$YF3[1:40]-sim$YF3[1:40], type = "l", col = "blue");
lines(sim3$YF4[1:40]-sim$YF4[1:40], type = "l", col = "orange");
legend("topright", legend=c("Group 1", "Group 2", "Group 3", "Group 4"), col=c("red", "green","blue","orange"), lty=1)
# Prevalence: Asymptomatic Male
plot(x=times[1:40], y=sim3$AM1[1:40]-sim$AM1[1:40], xlab = "time", ylab = "change in Prevalence", type = "l", ylim=c(-2e-7, 0),
     main = "Prevalence of Asymptomatic Infections among Males (beta)", col = "red");
lines(sim3$AM2[1:40]-sim$AM2[1:40], type = "l", col = "green");
lines(sim3$AM3[1:40]-sim$AM3[1:40], type = "l", col = "blue");
lines(sim3$AM4[1:40]-sim$AM4[1:40], type = "l", col = "orange");
legend("topright", legend=c("Group 1", "Group 2", "Group 3", "Group 4"), col=c("red", "green","blue","orange"), lty=1)

# Prevalence: Asymptomatic Female
plot(x=times[1:40], y=sim3$AF1[1:40]-sim$AF1[1:40], xlab = "time", ylab = "change in Prevalence", type = "l", ylim=c(-2e-7, 0),
     main = "Prevalence of Asymptomatic Infections among Females (beta)", col = "red");
lines(sim3$AF2[1:40]-sim$AF2[1:40], type = "l", col = "green");
lines(sim3$AF3[1:40]-sim$AF3[1:40], type = "l", col = "blue");
lines(sim3$AF4[1:40]-sim$AF4[1:40], type = "l", col = "orange");
legend("topright", legend=c("Group 1", "Group 2", "Group 3", "Group 4"), col=c("red", "green","blue","orange"), lty=1)


