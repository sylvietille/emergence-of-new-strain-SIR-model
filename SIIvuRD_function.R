#Load following packages 
library(deSolve)
library(reshape2)
library(tidyverse)
library(RColorBrewer)

#define time sequence 
time<-function(start=0, end=50, by=0.01){
  return(seq(start, end, by ))
}

#parameter values for the model 
#alphaI-betaI: effect of vaccine on suscpetible and infected individual
#alphaIv-betaIv: effect of vaccine on vaccinated suscpetible and vaccinated infected individuals
#alphaIu-betaIu: effect of vaccine on mutated infected individuals
#alphaIvu-betaIvu: effect of vaccine on vaccinated mutated infected and vaccinted suscpetible individuals
#vmax : vaccination rate [0,1], proportion of susceptible vaccinated per unit of time
#smax: proportion of suxceptible who will never be vaccinated [0,1]
#d, dv, du, dvu : death rate [0,1]
#r, rv, ru ,rvu : recovery rate [0,1]
#delta: virulence of the mutant 
parameters=c(alphaI=0.,betaI=0., vmax=0.12,smax=0., tsi=1.5, d=0.1, dv=0.01,r=0.1,rv=0.1,du=0.1,ru=0.1,dvu=0.01,rvu=0.1,alphaIv=1,betaIv=1,u=0.0,alphaIu=-0.2,betaIu=-0.2 ,alphaIvu=0.8, betaIvu=0.8,delta=-2)

#define initial state of model groups
#S: Susceptible, Sv: Vaccinated Susceptible, I: Infected, Iv: Vaccinated Infected 
#D: Death, R: Recovered, Iu: Mutated Infected, Ivu : Vaccinated Mutated Infected
state=c(S=0.999, Sv=0, I=0.001, Iv=0.0,D=0, R=0, Iu=0.00, Ivu=0.0)


#function that describe the vaccination 
eff.vac.rate = function(state,para)
{
  S=state[1]
  Smax=para[4]
  vmax=para[3]
  v.eff = (S>0)*(S>Smax)*vmax
  return(v.eff)
}



#define the SIR model
sirvd_mu <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS<-- eff.vac.rate(state,parameters) -S*I*tsi*(1-alphaI)*(1-betaI) - S*Iv*(1-betaIv)*tsi*(1-alphaI) -S*Iu*tsi*(1-alphaI)*(1-delta)*(1-betaIu) -S*Ivu*tsi*(1-delta)*(1-betaIvu)*(1-alphaI)
    dSv<- eff.vac.rate(state,parameters) -Sv*I*(1-alphaIv)*(1-betaI)*tsi - Sv*Iv*(1-alphaIv)*(1-betaIv)*tsi -Sv*Iu*tsi*(1-alphaIvu)*(1-betaIu)*(1-delta) -Sv*Ivu*tsi*(1-alphaIvu)*(1-betaIvu)*(1-delta)
    dI<-S*I*tsi*(1-alphaI)*(1-betaI) + S*Iv*(1-betaIv)*tsi*(1-alphaI) -I*r - I*d -I*u
    dIv<-Sv*I*(1-alphaIv)*(1-betaI)*tsi + Sv*Iv*(1-alphaIv)*(1-betaIv)*tsi - Iv*rv - Iv*dv -Iv*u
    dIu<- S*Iu*tsi*(1-alphaI)*(1-delta)*(1-betaIu) + S*Ivu*tsi*(1-delta)*(1-betaIvu)*(1-alphaI)  -Iu*ru -Iu*du +I*u
    dIvu<- Sv*Iu*tsi*(1-alphaIvu)*(1-betaIu)*(1-delta) + Sv*Ivu*tsi*(1-alphaIvu)*(1-betaIvu)*(1-delta)  -Ivu*rvu - Ivu*dvu +Iv*u
    dD<- I*d + Iv*dv + Iu*du+ Ivu*dvu
    dR<-I*r + Iv*rv + Iu*ru + Ivu*rvu
    return(list(c(dS,dSv,dI,dIv,dD,dR,dIu,dIvu)))
  })
}
  
  
  #plot all groups of a model
  plot_all_groups<-function(out){
    par(mfrow = c(1,2))
    plot(x = out[,1], y = out[,2], type = "l", xlab = "time", ylab="population size",ylim=c(0,1), col="red", lwd=2, main="Population")
    lines(x = out[,1], y = out[,3], type = "l", col="green", lwd=2)
    lines(x = out[,1], y = out[,4], type = "l", col="orange", lwd=2)
    lines(x = out[,1], y = out[,5], type = "l", col="purple", lwd=2)
    lines(x = out[,1], y = out[,6], type = "l", col="black", lwd=2)
    lines(x = out[,1], y = out[,7], type = "l", col="blue", lwd=2)
    lines(x = out[,1], y = out[,8], type = "l", col="yellow", lwd=2)
    lines(x = out[,1], y = out[,9], type = "l", col="pink", lwd=2)  
    
    plot.new()
    legend("center", inset=c(-0.6,-0.30),legend = c("Susceptible ", "Vaccinated Susceptible","Infected","Vaccinated Infected","Death","Recovered",'Mutated infected', 'Mutated Vaccinated infected'), lwd=c(2, 2, 2), col=c("red", "green","orange", "purple","black","blue",'yellow','pink'))

  }

  
  #calculate values of the model in time
  cal_model<-function(times, parameters,state){
    
    out <- ode(y = state, times = times, func = sirvd_mu2, parms = parameters)
    out <- as.data.frame(out)
    
    return( out)
  }
#calculate the probability of invasion and plot it against smax 
  rescu_proba_smax=function(times, state, para,sequ=0.01,N=1000,u=0.01){
    multiplier = seq(0,1,by = sequ) 
    para[4]=1
    
    prob.inv = rep(0,length(multiplier))
    for (i in 1:length(multiplier)){
      
      parameters=c(para[1],para[2], para[3],para[4]*multiplier[i], para[5],para[6], para[7],para[8],para[9],para[10],para[11],para[12],para[13],para[14],para[15],para[16],para[17],para[18],para[19], para[20], para[21],para[22] )
      out <- ode(y = state, times = times, func = sirvd_mu2, parms = parameters)
      out <- as.data.frame(out)
      dt = 0.01
      names(out) = c("Time","Susceptible","Vaccinated Susceptible", "Infected", "Vaccinated Infected", "Death", "Recovered","Mutated infected", "Mutated Vaccinated infected")
      dat.SIR = melt(out,id="Time",measure = c("Susceptible","Vaccinated Susceptible", "Infected", "Vaccinated Infected", "Death", "Recovered","Mutated infected", "Mutated Vaccinated infected"))
      
      
      names(dat.SIR) = c("Time","Compartment","Value")
      
      S = filter(dat.SIR,Compartment == "Susceptible")
      Sv = filter(dat.SIR,Compartment == "Vaccinated Susceptible")
      I = filter(dat.SIR,Compartment == "Infected")
      Iv = filter(dat.SIR,Compartment == "Vaccinated Infected")
      Iu=filter(dat.SIR,Compartment == "Mutated infected")
      Ivu=filter(dat.SIR,Compartment == "Mutated Vaccinated infected")
      
      #u=para[16]
      ru=para[11]
      rvu=para[13]
      du=para[10]
      dvu=para[12]
      tsi=para[5]
      delta=para[21]
      alphaIu=para[17]
      alphaIvu=para[19]
      betaIvu=para[20]
      
      
      Iu_growthrate= S$Value*tsi*(1-delta) + S$Value*tsi*(1-delta)*(1-betaIvu) #recovery rate and death rate ? 
      Ivu_growthrate= Sv$Value*tsi*(1-alphaIu) + Sv$Value*tsi*(1-alphaIvu)*(1-betaIvu)  
      
      prob.inv[i] = 1-prod(1-(I$Value*Iu_growthrate+Iv$Value*Ivu_growthrate)*N*u*dt)
    }  
    
    plot(para[4]*multiplier,prob.inv,type="l",xlab = "smax",main = "",xlim=c(0,1))
    
    return(prob.inv)
  }


#calculate the probability of invasion and plot it against vmax 
rescu_proba_vmax=function(times, state, para,sequ=0.01,N=1000,u=0.01){
  multiplier = seq(0,1,by = sequ) 
  para[3]=1

  prob.inv = rep(0,length(multiplier))
  for (i in 1:length(multiplier)){
    
    parameters=c(para[1],para[2], para[3]*multiplier[i],para[4], para[5],para[6], para[7],para[8],para[9],para[10],para[11],para[12],para[13],para[14],para[15],para[16],para[17],para[18],para[19], para[20], para[21],para[22] )
    out <- ode(y = state, times = times, func = sirvd_mu2, parms = parameters)
    out <- as.data.frame(out)
    dt = 0.01
    names(out) = c("Time","Susceptible","Vaccinated Susceptible", "Infected", "Vaccinated Infected", "Death", "Recovered","Mutated infected", "Mutated Vaccinated infected")
    dat.SIR = melt(out,id="Time",measure = c("Susceptible","Vaccinated Susceptible", "Infected", "Vaccinated Infected", "Death", "Recovered","Mutated infected", "Mutated Vaccinated infected"))

    names(dat.SIR) = c("Time","Compartment","Value")
    
    S = filter(dat.SIR,Compartment == "Susceptible")
    Sv = filter(dat.SIR,Compartment == "Vaccinated Susceptible")
    I = filter(dat.SIR,Compartment == "Infected")
    Iv = filter(dat.SIR,Compartment == "Vaccinated Infected")
    Iu=filter(dat.SIR,Compartment == "Mutated infected")
    Ivu=filter(dat.SIR,Compartment == "Mutated Vaccinated infected")
    
    #u=para[16]
    ru=para[11]
    rvu=para[13]
    du=para[10]
    dvu=para[12]
    tsi=para[5]
    delta=para[21]
    alphaIu=para[17]
    alphaIvu=para[19]
    betaIvu=para[20]
    
    
    Iu_growthrate= S$Value*tsi*(1-delta) + S$Value*tsi*(1-delta)*(1-betaIvu) #recovery rate and death rate ? 
    Ivu_growthrate= Sv$Value*tsi*(1-alphaIu) + Sv$Value*tsi*(1-alphaIvu)*(1-betaIvu)  
    
    prob.inv[i] = 1-prod(1-(I$Value*Iu_growthrate+Iv$Value*Ivu_growthrate)*N*u*dt)
  }  

  plot(para[3]*multiplier,prob.inv,type="l",xlab = "vmax",main = "",xlim=c(0,0.2))
  return(prob.inv)
}

#calculate the number of death compare to vmax and plot it 
death_end_vmax=function(times, para,state, sequ=0.01){
  data<-data.frame(matrix(ncol=2, nrow=0))
  colnames(data)<-c("D", "vmax")
  multiplier = seq(0,1,by = sequ) 
  para[3]=1
  
  for (i in 1:length(multiplier)){
    parameters=c(para[1],para[2],para[3]*multiplier[i],para[4], para[5],para[6], para[7],para[8],para[9],para[10],para[11],para[12],para[13],para[14],para[15],para[16],para[17],para[18],para[19], para[20], para[21],para[22] )
    
    out <- ode(y = state, times = times, func = sirvd_mu2, parms = parameters)
    out <- as.data.frame(out)
    data[nrow(data) +1 ,]<-c(out[dim(out)[1],6], para[3]*multiplier[i])
  }
  plot(data$vmax,data$D,type="l",xlab = "vmax",ylab="Death")
  legend("bottomright", inset=c(-0.30,-0), legend=c(paste("vmax =",para[3]),paste("smax=",para[4]),paste("tsi=",para[5]),paste("tsiu=",para[17]),paste("alpha=",para[1]),paste("alpha2=",para[14]), paste("beta2=",para[15]),paste("alpha3=",para[18]), paste("alpha4=",para[20]), paste("beta4=", para[21])))
  legend("bottomright", inset=c(-0.55,-0), legend=c(paste("u=",para[16]),paste("r=",para[8]),paste("rv=",para[9]),paste("ru=",para[11]),paste("rvu=",para[13]),paste("d=",para[6]),paste("dv=",para[7]),paste("du=",para[10]), paste("dvu=",para[12]),paste("tsiv=",parameters[22])))
  return(data)
}


model=cal_model(times, parameters, state)
plot_all_groups(model)
v=rescu_proba_vmax(times,state,parameters,sequ=0.01,u=10^-6,N=100000)
d=death_end_vmax(times,parameters,state,0.01)
s=rescu_proba_smax(times,state,parameters,sequ=0.01,u=1/1000000,N=100000)
