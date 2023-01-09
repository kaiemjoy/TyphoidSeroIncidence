### KA Parameter Extract for tidy data
library(reshape2)


#output <- "output_under10"
#load(paste(output, "/v9na.pred.rda", sep=""))
load(paste("output/v9na.pred.rda", sep=""))

# parnum: use y0=1; y1=2; t1=3; alpha=4; shape=5


d <- melt(predpar) %>%
  rename(antigen_iso = Var1,
         parameter = Var2,
         iter = Var3) %>% 
  mutate(antigen_iso = factor(antigen_iso),
         parameter = factor(parameter)) %>% 
  mutate(antigen_iso = factor(antigen_iso, labels = c("HlyE_IgA", "HlyE_IgG", "LPS_IgA", "LPS_IgG", "MP_IgA", "MP_IgG", "Vi_IgG")),
         parameter = factor(parameter, labels = c("y0", "y1", "t1", "alpha", "r"))) %>%
  filter(antigen_iso != "MP_IgA" & antigen_iso != "MP_IgG") %>%
  mutate(value = exp(value)) %>%
  mutate(value = ifelse(parameter == "r", value+1, value)) %>%
  pivot_wider(names_from = "parameter", values_from="value") %>%
  rowwise() %>%
  mutate(y1 = y0+y1) 


write_rds(d, "output/mcmc.rds")



  ##SeroCourse
  
  tx2 <- 10^seq(-1,3,0.025)
  
  bt <- function(y0,y1,t1) log(y1/y0)/t1
  
  
  # y0 <-4.501250e-01
  # y1 <- 9.8028190
  # shape <- 0.8518780
  # alpha <- 0.011440335
  # t <- 841.3951
  # t1 <- 6.8954378
  
  ab <- function(t,y0,y1,t1,alpha,shape) {
    beta <- bt(y0,y1,t1);
    yt <- 0;
    if(t <= t1) yt <- y0*exp(beta*t);
    if(t > t1) yt <- (y1^(1-shape)-(1-shape)*alpha*(t-t1))^(1/(1-shape));
    return(yt);
  }
  

  ##Overall SeroCourse
  dT <- data.frame(t=tx2) %>%
    mutate(ID = 1:n()) %>%
    pivot_wider(names_from = ID, values_from = t, names_prefix = "time") %>%
    slice(rep(1:n(), each = nrow(d))) 
  
  
  temp <- cbind(d, dT)  %>% pivot_longer(cols = starts_with("time"), values_to = "t") %>% select(-name)  %>%
    rowwise() %>%
    mutate(res = ab(t,y0,y1,t1,alpha,r)) 

  
  
  Serocourse <- temp %>% group_by(antigen_iso, t) %>%
    summarise(res.med  = quantile(res, 0.5),
              res.low  = quantile(res, 0.025),
              res.high = quantile(res, 0.975)) %>%
    pivot_longer(names_to = "quantile", cols = c("res.med","res.low","res.high"), names_prefix = "res.", values_to = "res") %>%
    mutate(agecat = "all")
  
  

  
  write_rds(Serocourse, "output/serocourse.rds")
  
  