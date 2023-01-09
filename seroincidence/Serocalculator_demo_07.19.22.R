

#### Typhoid Seroincidence Tutorial

#Setup
library(tidyverse)




##Load seroincidence package
##Still in development as of 05/12/22
#library(devtools)
#install_github("UCD-SERG/serocalculator")
library(serocalculator)


##Load longitudinal parameter mcmc data and prepare for model
c.hlye.IgG <- readRDS("dmcmc_hlyeigg_09.30.rds") %>%
  mutate(alpha = alpha*365.25, 
         d = r-1) %>%
  select(y1, alpha, d) 


##Load simulated data (lambda = .2) and prepare for model
p.hlye.IgG  <- read_csv("simpophlyeigg.2.csv") %>%
  rename(y= y.smpl,
         a = a.smpl) %>% 
  select(y, a)



## Conditions based on how i simulated the data
cond.hlye.IgG <- data.frame(nu = 2.869835,             # B noise, new noise parameter includes MGH + CA Facts controls 
                            eps = 0.2,            # M noise
                            y.low = 0.0,          # low cutoff
                            y.high = 5e4); 


## Seroincidence estimation
start <- .05

lambda = start # initial estimate: starting value
log.lambda = log(lambda)
log.lmin=log(lambda/10)
log.lmax=log(10*lambda) 



objfunc <- function(llam){
  return(res <- fdev(llam, p.hlye.IgG, c.hlye.IgG, cond.hlye.IgG))
}


fit <- nlm(objfunc,log.lambda,
           hessian=TRUE,print.level=0,stepmax=(log.lmax-log.lmin)/4)


#lambda, lower, upper, LF min
log.lambda.est <- c(exp(fit$estimate),
                    exp(fit$estimate + qnorm(c(0.025))*sqrt(1/fit$hessian)),
                    exp(fit$estimate + qnorm(c(0.975))*sqrt(1/fit$hessian)),
                    fit$minimum)


log.lambda.est
