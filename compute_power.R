library(ggplot2)
library(boot)
library(statmod)
library(numbers)

### Projected sample population sizes
### n_cases: number of case subjects
### n_con: number of control subjects
n_cases <- 5000
n_con <- 9000

# Read in the raw output from the simulations
dat_use <- read.csv("scenario_1.csv", header=T, colClasses = c("numeric", "integer", "integer", "integer", "integer", "integer"))
dat_use$TotCars <- dat_use$CaseCars + dat_use$ControlCars
dat_use$TotPop <- dat_use$TotControls + dat_use$TotCases

dat_use$AF_cases <- dat_use$CaseCars / dat_use$TotCases   
dat_use$AF_controls <- dat_use$ControlCars / dat_use$TotControls

# Average the number of carriers per simulation
carriers <- aggregate(dat_use$TotCars, by=list(dat_use$Penetrance), FUN=mean, na.rm=TRUE)
use <- data.frame(carriers)
agg_pop <- aggregate(dat_use$TotPop, by=list(dat_use$Penetrance), FUN=mean, na.rm=TRUE)

agg_case_cars <- aggregate(dat_use$CaseCars, by=list(dat_use$Penetrance), FUN=mean, na.rm=TRUE)
agg_control_cars <- aggregate(dat_use$ControlCars, by=list(dat_use$Penetrance), FUN=mean, na.rm=TRUE)
agg_cases  <- aggregate(dat_use$TotCases, by=list(dat_use$Penetrance), FUN=mean, na.rm=TRUE)
agg_controls  <- aggregate(dat_use$TotControls, by=list(dat_use$Penetrance), FUN=mean, na.rm=TRUE)


AF_cases <- aggregate(dat_use$AF_cases, by=list(dat_use$Penetrance), FUN=mean, na.rm=TRUE)
AF_controls <- aggregate(dat_use$AF_controls, by=list(dat_use$Penetrance), FUN=mean, na.rm=TRUE)


#################################################################################
#################################################################################
#################################################################################
### POWER CALCULATION 
#################################################################################
#################################################################################
#################################################################################
f_cases <- n_cases/(n_cases+n_con)
f_con <- n_con/(n_cases+n_con)

n_entries <- dim(dat_use)[1]
dat_use$power <- NA
for (curr in 1:n_entries) {
  AF_cases <- dat_use$CaseCars[curr] / dat_use$TotCases[curr]   
  AF_controls <- dat_use$ControlCars[curr] / dat_use$TotControls[curr]
  cum_pow <- power.fisher.test(AF_cases, AF_controls, n_cases, n_con, alpha=5e-07, alternative="two", nsim=100)
  dat_use$power[curr] <- cum_pow
}

powers <- aggregate(dat_use$power, by=list(dat_use$Penetrance), FUN=mean, na.rm=TRUE)
power_CI <- aggregate(dat_use$power, by=list(dat_use$Penetrance), FUN=quantile, na.rm=TRUE)

p<-ggplot() + geom_line(data=powers, aes(x=penetrances, y=powers$x), color="red") +xlab('Penetrance')+ ylab('Power') +
   geom_ribbon(data=power_CI, aes(x= penetrances,ymin=power_CI$x[,2], ymax=power_CI$x[,4]), alpha=0.2) + ylab('Power')
p

write.csv(dat_use, "power_calculation.csv", row.names = F, quote=F)



