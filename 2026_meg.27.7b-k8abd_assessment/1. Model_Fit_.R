## Run analysis, write model results

library(devtools)
library(FLCore)
library(icesTAF)
library(FLa4a)
library(FLEDA)
library(ggplotFL)
library(gridExtra)
library(icesAdvice)
library(FLFishery)
library(FLasher)
library(msy)
library(a4adiags)
library(devtools)
library(ggplot2)
library(data.table)
library(readr)
library(tidyverse)

setwd("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2026/3.Assessment_october")

set.seed(1234)

# Functions ---
SavePlot <- function(plotname, width = 10, height = 7) {
  file <- paste0(getwd(), "/report/", "meg78_data_", plotname, ".png")
  dev.print(png, file, width = width, height = height, units = "in", res = 300)
}
iterMedians <- function(x) {
  # x es un FLQuant con iter
  apply(x, c(1,2,3,4,5), median, na.rm = TRUE) |>
    FLQuant(dimnames = c(dimnames(x)[1:5], list(iter = "1")),
            units    = units(x))
}
fitMedian <- function(fit) {
  
  # 1. Mediana del stock (FLStock)
  iterMedians <- function(x) {
    apply(x, c(1,2,3,4,5), median, na.rm = TRUE) |>
      FLQuant(dimnames = c(dimnames(x)[1:5], list(iter = "1")),
              units    = units(x))
  }
  
  stock_med <- qapply(fit, iterMedians)
  
  # 2. Predecir todas las iteraciones
  fitted_all <- predict(pars(fit))
  
  # 3. Mediana de cada FLQuant dentro de predict()
  fitted_med <- lapply(fitted_all, function(sublist) {
    lapply(sublist, function(q) {
      arr <- array(q[], dim = dim(q), dimnames = dimnames(q))
      med <- apply(arr, c(1,2,3,4,5), median, na.rm = TRUE)
      dim(med) <- c(dim(q)[1:5], 1)
      dn <- dimnames(q)
      dn$iter <- "1"
      FLQuant(med, dimnames = dn, units = units(q))
    })
  })
  
  # 4. Reconstruir un objeto SCA con stock + fitted_med
  out <- list(
    stock   = stock_med,
    stkmodel = fitted_med$stkmodel,
    qmodel   = fitted_med$qmodel,
    vmodel   = fitted_med$vmodel
  )
  
  return(out)
}
residualsMedian <- function(res) {
  
  # Función interna para colapsar un FLQuant a la mediana
  iterMedians <- function(x) {
    arr <- array(x[], dim = dim(x), dimnames = dimnames(x))
    med <- apply(arr, c(1,2,3,4,5), median, na.rm = TRUE)
    
    # reconstruir FLQuant con iter = 1
    dim(med) <- c(dim(x)[1:5], 1)
    dn <- dimnames(x)
    dn$iter <- "1"
    
    FLQuant(med, dimnames = dn, units = units(x))
  }
  
  # Aplicar a cada FLQuant dentro del objeto a4aFitResiduals
  res_med_list <- lapply(res@.Data, iterMedians)
  
  # Reconstruir un objeto a4aFitResiduals
  out <- new("a4aFitResiduals")
  out@.Data <- res_med_list
  out@names <- res@names
  out@desc  <- paste(res@desc, "(median collapsed)")
  out@lock  <- FALSE
  
  return(out)
}
source('./doc/retro_analysis_f.R') 
source('./doc/retro_analysis_f_mcmc.R')

# Load data ----
load("bootstrap/data/stock/meg78_stock.RData")
load("bootstrap/data/indices/meg78_indices.RData")

# Modify the data ----

## Surveys ----

index<- tun.sel[c("SP_PORC","CPUE.IRLFRsurvey")]

plot(index[[1]]@index[1,])+geom_point()+ggtitle("Porcupine age 1")+ylab("Ind./30min haul")
plot(index[[2]]@index[1,])+geom_point()+ggtitle("IRLFR age 1")+ylab("numbers per 10 hours fished")

index[[2]]@index[ac(1:3),ac(2015:2021)] <- NA # IRLFR
index[[1]]@index[ac(1:3),ac(2015:2021)] <- NA # Porcupine

## Catch ----

plot(stock@catch.n['1',])
stock@catch.n['1',as.character(1984:2000)] <- NA # We do not really believe that the increase in 1-year-olds in the catch is real so we shouldn't formulate a model that treats it as real.

# Define the submodels  ----

qmod <- list(~s(age, k = 4), ~s(age, k = 4))
fmod <- ~s(age, k = 5, by = breakpts(year, c(1990,2013))) + factor(year)
srmod <- ~ factor(year)
n1mod <- ~s(age, k = 3)
vmod <- list(~s(age, k = 3), ~1, ~1)

# Fit the assessment ----

## Unique iter ----
fit01 <- sca(stock, index, fmodel = fmod, qmodel = qmod, srmodel = srmod, vmodel = vmod, n1model = n1mod)
stk01 <- stock + fit01

### Sumamry ----
plot(stk01)
submodels(fit01)

AIC(fit01)
BIC(fit01)

### Retrospective analysis ----
results <- run_retro_analysis(stock, index, fit01, fmod, qmod, srmod, vmod, n1mod)

results$rho_table <- results$rho_table %>% mutate(x = 2025, y = 0)
results$rho_table$qname <- c("F" = "F", "SSB" = "SB", "Recruitment" = "Rec", "Catch" = "C")

new_names <- c("Rec" = "Recruitment", "SB" = "SSB", "C" = "Catch", "F" = "F")
plot(FLStocks(results$retro), col = 1, lwd = 1) +
  facet_wrap(~qname, scales = 'free_y', labeller = labeller(qname = new_names)) +
  geom_text(
    data = results$rho_table,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1, vjust = 0) +
  theme_bw() +
  labs(color = "N years removed")

### Selectivity and catchability ----
fitted <- predict(pars(fit01))
a  <- xyplot(data~age,groups=year,stk01@harvest,type='b',ylim=c(0,1.5),ylab='F',main='Fishing mortality')
a1 <- xyplot(data~age,groups=year,data=fitted$qmodel[1],type='b',ylab='Catchability',main="Porcupine")
a2 <- xyplot(data~age,groups=year,data=fitted$qmodel[2],type='b',ylab='Catchability',main="IRLFR")
grid.arrange(a,a1,a2,ncol=2)

### Residuals plot ----
res <- residuals(fit01, stock, index)
plot(res)

### Observed and predicted catches plot ----
fits01 <- simulate(fit01, 1000)
stks01 <- stock + fits01

pred <- catch(stks01)

pred_median <- apply(pred, 2, median, na.rm = TRUE)
pred_p5     <- apply(pred, 2, quantile, 0.05, na.rm = TRUE)
pred_p95    <- apply(pred, 2, quantile, 0.95, na.rm = TRUE)

df <- tibble(
  Year      = as.numeric(dimnames(pred)$year),
  Observed  = as.numeric(catch(stock)),
  Median    = as.numeric(pred_median),
  P5        = as.numeric(pred_p5),
  P95       = as.numeric(pred_p95)
)

ggplot(df, aes(x = Year)) +
  geom_ribbon(aes(ymin = P5, ymax = P95),
              fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = Median, color = "Predicted"), size = 1.2) +
  geom_line(aes(y = Observed, color = "Observed"), size = 1.2) +
  scale_color_manual(values = c("Observed" = "black",
                                "Predicted" = "steelblue4")) +
  labs(x = "Year",
       y = "Catch (tonnes)",
       color = "Type") +
  theme_bw()

## MCMC ----

mcmc_ctrl <- SCAMCMC(
  mcmc   = 1600000, 
  mcsave = 1000, mcseed = 10) 

fits <- sca(stock, index, fmodel = fmod, qmodel = qmod,
            srmodel = srmod, vmod = vmod, n1mod = n1mod,
            fit = "MCMC", mcmc = mcmc_ctrl)
fits <- burnin(fits, 600)

stks <- stock + fits

### Median ----
fit1 <- fitMedian(fits)
stk1 <- qapply(stks, iterMedians)

### Sumamry ----
plot(stks)
submodels(fits)

### Save ----
save(fmod, qmod, srmod, vmod, n1mod, file = "model/SubModels_Final_october.Rdata")
save(stk1, stk01, stks, stock, index, fit1, fit01, fits, file = "model/MegFit_Final_october.Rdata")
stk_est <- stk1
save(stk_est, file = "model/Meg78_WGMIX_stock_estimated.Rdata")
stk_obs <- stock
save(stk_obs, file = "model/Meg78_WGMIX_stock_input.Rdata")

### Observed and predicted catches  ----

pred <- catch(stks)

pred_median <- apply(pred, 2, median, na.rm = TRUE)
pred_p5     <- apply(pred, 2, quantile, 0.05, na.rm = TRUE)
pred_p95    <- apply(pred, 2, quantile, 0.95, na.rm = TRUE)

df <- tibble(
  Year      = as.numeric(dimnames(pred)$year),
  Observed  = as.numeric(catch(stock)),
  Median    = as.numeric(pred_median),
  P5        = as.numeric(pred_p5),
  P95       = as.numeric(pred_p95)
)

ggplot(df, aes(x = Year)) +
  geom_ribbon(aes(ymin = P5, ymax = P95),
              fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = Median, color = "Predicted"), size = 1.2) +
  geom_line(aes(y = Observed, color = "Observed"), size = 1.2) +
  scale_color_manual(values = c("Observed" = "black",
                                "Predicted" = "steelblue4")) +
  labs(x = "Year",
       y = "Catch (tonnes)",
       color = "Type") +
  theme_bw()

SavePlot('Catch_obs_pred')

write.table(df,'report/Obs_pred_catch.csv', row.names = FALSE, sep = ";")

### Retrospective analysis ----

results <- run_retro_analysis_mcmc(stock, index,
                                   stks, fmod, qmod,
                                   srmod, vmod, n1mod,
                                   mcmc_ctrl,
                                   back = 5)

results$rho_table <- results$rho_table %>% mutate(x = 2025, y = 0)
results$rho_table$qname <- c("F" = "F", "SSB" = "SB", "Recruitment" = "Rec", "Catch" = "C")

new_names <- c("Rec" = "Recruitment", "SB" = "SSB", "C" = "Catch", "F" = "F")

f0 <- as.data.frame(fbar(results$terminal_full)) %>%
  group_by(year) %>%
  summarise(low=quantile(data,.025),
            med=median(data),
            high=quantile(data,.975)) %>%
  mutate(qname="F")

ssb0 <- as.data.frame(ssb(results$terminal_full)) %>%
  group_by(year) %>%
  summarise(low=quantile(data,.025),
            med=median(data),
            high=quantile(data,.975)) %>%
  mutate(qname="SB")

rec0 <- as.data.frame(results$terminal_full@stock.n[1,]) %>%
  group_by(year) %>%
  summarise(low=quantile(data,.025),
            med=median(data),
            high=quantile(data,.975)) %>%
  mutate(qname="Rec")

catch0 <- as.data.frame(catch(results$terminal_full)) %>%
  group_by(year) %>%
  summarise(low=quantile(data,.025),
            med=median(data),
            high=quantile(data,.975)) %>%
  mutate(qname="C")

uncert <- bind_rows(f0, ssb0, rec0, catch0)

plot(FLStocks(results$retro), col = 1, lwd = 1) +
  
  geom_ribbon(
    data = uncert,
    aes(x = year, ymin = low, ymax = high),
    inherit.aes = FALSE,
    fill = "red",
    alpha = 0.15
  ) +
  
  facet_wrap(
    ~qname,
    scales = "free_y",
    labeller = labeller(qname = new_names)
  ) +
  
  # geom_text(
  #   data = results$rho_table,
  #   aes(x = x, y = y, label = label),
  #   inherit.aes = FALSE,
  #   hjust = 1,
  #   vjust = 1
  # ) +
  
  theme_bw() +
  labs(color = "N years removed")

# plot(FLStocks(results$retro), col = 1, lwd = 1) +
#   facet_wrap(~qname, scales = 'free_y', labeller = labeller(qname = new_names)) +
#   geom_text(
#     data = results$rho_table,
#     aes(x = x, y = y, label = label),
#     inherit.aes = FALSE,
#     hjust = 1, vjust = 0) +
#   theme_bw() +
#   labs(color = "N years removed")

save(results, file = "model/MegRetroFit_FINAL.Rdata")

SavePlot('Retro_Fit')

### Selectivity and catchability plot ----

a  <- xyplot(data~age,groups=year,stk1@harvest,type='b',ylim=c(0,1.5),ylab='F',main='Fishing mortality')
a1 <- xyplot(data~age,groups=year,data=fit1$qmodel[1],type='b',ylab='Catchability',main="Porcupine")
a2 <- xyplot(data~age,groups=year,data=fit1$qmodel[2],type='b',ylab='Catchability',main="IRLFR")

a
SavePlot('F')

grid.arrange(a,a1,a2,ncol=2)
SavePlot('F and Q',8,6)

### Residuals plot ----
res <- residuals(fits, stock, index)
res_median <- residualsMedian(res)
plot(res_median)
SavePlot('Residuals1',10,6)

# Res2
res_median$catch.n[is.na(res_median$catch.n)] <- 0 #hack
res_median$SP_PORC[is.na(res_median$SP_PORC)] <- 0 #hack
res_median$CPUE.IRLFRsurvey[is.na(res_median$CPUE.IRLFRsurvey)] <- 0 #hack
bubbles(res_median)
SavePlot('Residuals2',10,6)

# Res3
qqmath(res_median)
SavePlot('Residuals3')

### Summary plot ----
plot(stks) + facet_wrap(~qname,scales='free_y', labeller = labeller(qname = new_names)) +theme_bw()
SavePlot('summary0')

### STF settings ----

# Recruitment intermediate year and advice year:
years <- stk1@range[4]:stk1@range[5]
nyears <- length(years)

Rassumption <- "gmean"
if(Rassumption == "gmean"){
  # GM <- round(exp(mean(log(c(stk1@stock.n[1,1:(nyears-2)])))),0) # thousands. GM time series (without the 	last 2 years)
  GM <- round(exp(mean(log(c(stk1@stock.n[1,ac(2011:2024)])))),0) # thousands. GM time series of the low recruitment period (2011-2024)
}
if(Rassumption == "lastwoyr"){
  GM <- round(exp(mean(log(c(stk1@stock.n[1,(nyears-1):nyears])))),0) # thousands. GM last two years
}

# F status quo (Fsq)

# Extract fishing mortality (F) and plot it
f_values <- harvest(stk1)
f_mean_per_year <- apply(f_values[3:6,], 2, mean)
plot(f_mean_per_year)

# F setting
Fassumption <- "Funscaled"

if (Fassumption == "Fscaled"){
  fsq <- fbar(stk1)[,nyears] # F SCALED to Fbar of final assessment year (if there is a decreasing trend of F in the results of the assessment time-series)
}
if(Fassumption == "Funscaled"){
  fsq <-mean(fbar(stk1)[,nyears-2:0]) # F UNSCALED (Average F last 3 years) 
}

### Summary table SAG lognormal ----
# Summary table using %90 confidence interval
# Lognormal-based confidence intervals (90%) for key variables

stk1@stock.n[1, ac(2025)] <- GM

tsb <- tsb(stk1)
ssb <- ssb(stk1)
catchobs <- stock@catch
discardsobs <- stock@discards
catch <- stk1@catch
recr <- stk1@stock.n[1,]
fbar <- fbar(stk1)

# Parámetro z para 90% CI
z <- 1.6449

# Valores medios
tsb <- tsb(stk1)
ssb <- ssb(stk1)
recr <- stk1@stock.n[1,]
fbar <- fbar(stk1)

# Desviaciones estándar (estimadas a partir de las iteraciones en stks)
tsbse <- apply(tsb(stks), 2, sd)
ssbse <- apply(ssb(stks), 2, sd)
recrse <- apply(stks@stock.n[1,], 2, sd)
recrse[,ac(2025)] <- as.numeric(mean(recrse[1,ac(2011:2024),1,1,1,1]))
fbarse <- apply(fbar(stks), 2, sd)

# Convertir a vectores
tsb_mean <- as.numeric(tsb)
ssb_mean <- as.numeric(ssb)
recr_mean <- as.numeric(recr)
fbar_mean <- as.numeric(fbar)

# Catch, landings y discards
catch <- as.numeric(stk1@catch)
discards <- as.numeric(stock@discards)
landings <- catch - discards

# Función para calcular IC log-normal
calc_lognorm_CI <- function(mean, se, z = 1.6449) {
  log_sd <- sqrt(log(1 + (z * se / mean)^2))
  lower <- exp(log(mean) - log_sd)
  upper <- exp(log(mean) + log_sd)
  list(lower = lower, upper = upper)
}

# Aplicar a los indicadores
tsb_CI <- calc_lognorm_CI(tsb_mean, tsbse)
ssb_CI <- calc_lognorm_CI(ssb_mean, ssbse)
recr_CI <- calc_lognorm_CI(recr_mean, recrse)
fbar_CI <- calc_lognorm_CI(fbar_mean, fbarse)

# Añadir valores para 2025
tsb_mean <- c(tsb_mean, NA)
tsb_CI$lower <- c(tsb_CI$lower, NA)
tsb_CI$upper <- c(tsb_CI$upper, NA)

ssb_mean <- c(ssb_mean, NA)
ssb_CI$lower <- c(ssb_CI$lower, NA)
ssb_CI$upper <- c(ssb_CI$upper, NA)

recr_mean <- c(recr_mean, GM)
recr_CI$lower <- c(recr_CI$lower, NA)
recr_CI$upper <- c(recr_CI$upper, NA)

fbar_mean <- c(fbar_mean, fsq)
fbar_CI$lower <- c(fbar_CI$lower, NA)
fbar_CI$upper <- c(fbar_CI$upper, NA)

catch <- c(catch, NA)
landings <- c(landings, NA)
discards <- c(discards, NA)

# Crear tabla resumen final
summary_df <- data.frame(
  year = c(years, max(years) + 1),
  recr_lo = recr_CI$lower,
  recr = recr_mean,
  recr_hi = recr_CI$upper,
  tsb_lo = tsb_CI$lower,
  tsb = tsb_mean,
  tsb_hi = tsb_CI$upper,
  ssb_lo = ssb_CI$lower,
  ssb = ssb_mean,
  ssb_hi = ssb_CI$upper,
  fbar_lo = fbar_CI$lower,
  fbar = fbar_mean,
  fbar_hi = fbar_CI$upper,
  catch=c(catchobs,NA),
  lan=c(catchobs-discardsobs,NA),
  dis=c(discardsobs,NA)
)
summary_df <- as.data.table(summary_df)
summary_df[, c(2:10,14:16) := lapply(.SD, round, digits=0), .SDcols = c(2:10,14:16)]
summary_df[, c(11:13) := lapply(.SD, icesRound), .SDcols = c(11:13)]

write.csv(summary_df, "report/Summary_SAG_lognorm.csv", row.names = FALSE)

### Sensitivity for template - not needed anymore? ----
ages <- stk1@range[1]:stk1@range[2]
stk3y <- window(stk1, start = max(years) - 2)
p <- apply(landings.n(stk3y) / catch.n(stk3y), 1, mean, na.rm = T)
p <- c(ifelse(is.na(p), 0, p))
sen <- data.frame(
  Age = ages,
  M = c(apply(m(stk3y), 1, mean)),
  Mat = c(apply(mat(stk3y), 1, mean)),
  PF = c(apply(harvest.spwn(stk3y), 1, mean)),
  PM = c(apply(m.spwn(stk3y), 1, mean)),
  Sel = c(apply(harvest(stk3y), 1, mean)) * p,
  WeCa = c(apply(landings.wt(stk3y), 1, mean)),
  Fd = c(apply(harvest(stk3y), 1, mean)) * (1 - p),
  WeCad = c(apply(discards.wt(stk3y), 1, mean)),
  Fi = 0,
  WeCai = 0
)
knitr::kable(sen, row.names = F, digits = c(0, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3))
write.csv(sen, "report/Sen.csv", row.names = F)

### Fitted plots ----
plot(fit01, stock)
SavePlot("Fit1")
plot(fit01, index[1])
SavePlot("Fit2")
plot(fit01, index[2])
SavePlot("Fit3")

fitSumm(fits)

# iters
# 1
# nopar       112.000000
# nlogl               NA
# maxgrad             NA
# nobs        567.000000
# gcv                 NA
# convergence         NA
# accrate       0.331582

### KOBE plot ----

dev.off()

kobe <- data.frame(Btrig=50000,Fmsy=0.25,B=rev(c(ssb(stk1))),F=rev(c(fbar(stk1))))
with(kobe[1,],plot(B/Btrig,F/Fmsy,xlim=c(0,3.5),ylim=c(0,3),main='meg78'))
rect(-1,-1,5,5,col='yellow',border=NA)
rect(-1,1,1,5,col='red',border=NA)
rect(1,-1,5,1,col='green',border=NA)
with(kobe[1,],text(B/Btrig,F/Fmsy,max(years),pos=1))
with(kobe[1,],points(B/Btrig,F/Fmsy,pch=16))
with(kobe,lines(B/Btrig,F/Fmsy))
#plot(kobe)
SavePlot('Kobe',6,6)

dev.off()
