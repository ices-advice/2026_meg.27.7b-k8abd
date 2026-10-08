# model.R - Run forecast
# 2020_sol.27.4_forecast/model.R

# Copyright Iago MOSQUEIRA (WMR), 2020
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


# install.packages("devtools")
# install.packages("FLCore", repo = "http://flr-project.org/R")
# install_github("ices-tools-prod/msy")
# install.packages(c("FLa4a", "FLasher", "ggplotFL"), repos="https://flr-project.org/R")
# install.packages(c("ggplot2", "snpar", "foreach", "data.table"))
# devtools::install_github("flr/a4adiags")

library(FLCore)
library(icesTAF)
library(ggplot2)
library(ggplotFL)
library(FLFishery)
library(FLasher)
library(FLa4a)
library(a4adiags)
library(icesAdvice)
library(FLCore)
library(data.table)

setwd("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2026/3.Assessment_october")
mkdir("model")
set.seed(1234)

# setting  ----
Fassumption <- "Funscaled" #"Fscaled" or "Funscaled"
Rassumption <- "gmean" # "gmean" or "meanblowmean" or "lastwoyr"

# Functions  ----
SavePlot <- function(plotname, width = 10, height = 7) {
  file <- paste0(getwd(), "/model/", "meg78_fwd_", plotname, "_", Fassumption, "_", Rassumption, ".png")
  dev.print(png, file, width = width, height = height, units = "in", res = 300)
}
SavePlotBlim <- function(plotname, width = 10, height = 7) {
  file <- paste0(getwd(), "/report/", "meg78_probBlim_", plotname, ".png")
  dev.print(png, file, width = width, height = height, units = "in", res = 300)
}
iterMedians <- function(x) {
  # x es un FLQuant con iter
  apply(x, c(1,2,3,4,5), median, na.rm = TRUE) |>
    FLQuant(dimnames = c(dimnames(x)[1:5], list(iter = "1")),
            units    = units(x))
}

# Reference points ----

refPts <- data.frame(
  MSYBtrigger =  50000,
  Bpa = 50000,
  Blim = 36000,
  Fpa = 0.25,
  Fp05 = 0.25,
  Fmsy_unconstr = 0.25,
  Fmsy = 0.25,
  Fmsyupper_unconstr = 0.25,
  Fmsyupper = 0.25,
  Fmsylower_unconstr = 0.21,
  Fmsylower = 0.21
)

# Load assessment  ----

load('model/MegFit_Final_october.RData')

# Settings ----
run <- stks 
years <- run@range[4]:run@range[5]
nyears <-length(years)
ages <- run@range[1]:run@range[2]
nages <- length(ages)

# MODEL, ADVICE and FORECAST year

dy <- dims(run)$maxyear #last year of assessment input data
ay <- dy + 1
fy <- ay + 1

# TAC & advice current year

advice <- FLQuant(14869, dimnames=list(age='all', year=ay), units="tonnes") #en WGBIE 2025, advice para el año 2026 son 14869 tonas.
tac <- advice

### RECRUITMENT
if(Rassumption == "meanblowmean"){
# Mean of the R below historical mean (historical mean using all years minus last 2)
  gmean <- exp(mean(log(window(stock.n(run)["1",], end=-2))))
  recsq <- mean(window(stock.n(run)["1",], end=-2)[window(stock.n(run)["1",], end=-2) < gmean])
}
if(Rassumption == "gmean"){
# GEOMEAN but last year
#recsq <- exp(mean(log(window(stock.n(run)["1",], end=-1))))

# GEOMEAN using all years minus last 2 (distintas formas de obtener "recsq")
# specity the number of years removing from the end of the time-series
# recsq <- exp(mean(log(window(stock.n(run)["1",], end=-2))))
## specity the start and final years
# recsq <- exp(mean(log(window(stock.n(run)["1",], start=1984, end=2022))))

# GEOMEAN: para cambiar el reclutamiento del ultimo año del assessment por geomean si consideramos que tenemos mucha incertidumbre.
#plot(run@stock.n['1',])
#run@stock.n['1','2022'] <- recsq #AI: Para la evalucion con datos hasta 2022.

recsq <- round(exp(mean(log(c(stk1@stock.n[1,ac(2011:2024)])))),0) # thousands. GM time series of the low recruitment period (2011-2024)
}
if(Rassumption == "lastwoyr"){
  # GEOMEAN but lasttwo year
  recsq <- exp(mean(log(c(stock.n(run)[1,(nyears-1):nyears]))))
}

# --- SETUP future

run@stock.n[1, ac(2025)] <- recsq

# 3 years, 5 years wts/selex, 3 years discards
fut <- stf(run, nyears=3, wts.nyears = 5, fbar.nyears=5, disc.nyears=3)

### F status quo (Fsq)

if (Fassumption == "Fscaled"){
  Fsq <- expand(fbar(fut)[, ac(dy)], year=ay) # F SCALED to Fbar of final assessment year
}
if(Fassumption == "Funscaled"){
  Fsq <-expand(yearMeans(fbar(fut)[, ac(seq(dy - 2, dy))]), year= ay) # F UNSCALED (Average F last 3 years)
  # Fsq <-FLCore::expand(yearMeans(fbar(fut)[, ac(seq(dy - 2, dy))]), year= ay) # Be careful not to have package (tidyr)instaled. In case of error use: FLCore::expand
  
}

########################################################################################

# SET geomean SRR
gmsrr <- predictModel(model=rec~a, params=FLPar(c(recsq), units="thousands",
  dimnames=list(params="a", year=seq(ay, length=3), iter=1)))

# > gmsrr   # AI: reemplaza el reclutamiento con recsq del año intermedio (2025) y los dos siguientes.
# An object of class "FLQuants": EMPTY
# model:  
#   rec ~ a
# 
# params:  
#   An object of class "FLPar"
# year
# params    2025      2026      2027 
# a        158226      158226     158226 
# units:  thousands 

# GENERATE targets from refpts
refpts <- FLPar(refPts[1,])
targets <- expand(as(refpts, 'FLQuant'), year=fy)

# ##############################################################################
# #AI: added to include a new variable in the output
# # GENERATE targets from refpts
# # need to run the forecast once to see what the SSB at the start of the advice year is: as.numeric(ssb(runs[[1]])[,as.character(fy)])
#
# median(as.numeric(ssb(runs[[1]])[,as.character(fy)]))
# 44895.25 

# # Then enter it manually
# refPts <- cbind(refPts,  44895.25)
# dimnames(refPts)[[2]][ncol(refPts)] <- "constSSB"
# refpts <- FLPar(refPts[1,])
# targets <- expand(as(refpts, 'FLQuant'), year=fy)

# #AI: to complete in the advice the assumptions for interim year and forecast.
# 
# icesRound(median(fbar(runs$Fmsy)[, '2025'])) # 0.162
# icesRound(median(fbar(runs$Fmsy)[, '2026'])) # 0.176
# icesRound(median(fbar(runs$Fmsy)[, '2027'])) # 0.25
# 
# median(ssb(runs$Fmsy)[, '2025']) # 44828.54
# median(ssb(runs$Fmsy)[, '2026']) # 45505.06 # La SSB del año intermedio (2026) hay que incluirla en la tabla del advice donde se da la summary table pero este ultimo a?o no hay valor (NA).
# median(ssb(runs$Fmsy)[, '2027']) # 44895.25
# median(ssb(runs$Fmsy)[, '2028']) # 43058.47
# 
# median(catch(runs$Fmsy)[, '2025']) # 9428.746 # fitted catch for the last assessment year to put in the report Tab 5.1.1.2
# median(catch(runs$Fmsy)[, '2026']) # 11364.65
# median(catch(runs$Fmsy)[, '2027']) # 15420.54
# 
# median(landings(runs$Fmsy)[, '2025']) # 8038.226
# median(landings(runs$Fmsy)[, '2026']) # 9624.963
# median(landings(runs$Fmsy)[, '2027']) # 12916.04
# 
# median(rec(runs$Fmsy)[, '2027']) # 158226
# median(discards(runs$Fmsy)[, '2026']) # 1736.68
# 
# 
# # 
##############################################################################


# catch options table ----

# Targets
C0 <- FLQuant(0, dimnames=list(age='all', year=fy))
Fiy <- FLQuants(fbar=append(Fsq, C0))
Ciy <- FLQuants(catch=append(tac, C0))

# RUN for Fsq
Fsqrun <- fwd(fut, sr=gmsrr, control=as(Fiy, "fwdControl"))

# TEST if Cay < TACay
if(median(catch(Fsqrun)[, ac(ay)]) <= tac) {
  itarget <- Fiy
} else {
  itarget <- Ciy
  Fsqrun <- fwd(fut, sr=gmsrr, control=as(Ciy, "fwdControl"))
}

# DEFINE catch options
catch_options <- list(

  # FMSY
  Fmsy=FLQuants(fbar=targets["Fmsy",]),

  # lowFMSY
  lFmsy=FLQuants(fbar=targets["Fmsylower_unconstr",]),

  # uppFMSY
  uFmsy=FLQuants(fbar=targets["Fmsyupper",]),

  # F0
  F0=FLQuants(fbar=FLQuant(0, dimnames=list(age='all', year=fy))),

  # Fpa
  Fpa=FLQuants(fbar=targets["Fpa",]),

  # # Flim
  # Flim=FLQuants(fbar=targets["Flim",]),

  # Bpa
  Bpa=FLQuants(ssb_flash=targets["Bpa",]),

  # Blim
  Blim=FLQuants(ssb_flash=targets["Blim",]),

  # MSYBtrigger
  MSYBtrigger=FLQuants(ssb_flash=targets["MSYBtrigger",]),

  # F sq
  Fsq=FLQuants(fbar=expand(fbar(Fsqrun)[, ac(ay)], year=fy)),
  
  # SB sq
  SBsq=FLQuants(ssb_flash=expand(ssb(Fsqrun)[, ac(ay+1)], year=fy)),
  
  # TAC
  rotac=FLQuants(catch=expand(tac, year=fy))
)

# C0
F0 <- FLQuants(fbar=FLQuant(0, dimnames=list(age='all', year=fy + 1)))

# Expand the deterministic reference points to 100 iterations
nits <- dim(itarget$fbar)[6]

catch_options$Fmsy$fbar <- FLQuant(
  0.25,
  dimnames = dimnames(window(itarget, end=ay)$fbar)
)
dimnames(catch_options$Fmsy$fbar)$year <- as.character(fy)

F0$fbar <- propagate(F0$fbar, nits)
catch_options$Fmsy$fbar <- propagate(catch_options$Fmsy$fbar, nits)
catch_options$lFmsy$fbar <- propagate(catch_options$lFmsy$fbar, nits)
catch_options$uFmsy$fbar <- propagate(catch_options$uFmsy$fbar, nits)
catch_options$F0$fbar   <- propagate(catch_options$F0$fbar, nits)
catch_options$Fpa$fbar  <- propagate(catch_options$Fpa$fbar, nits)
catch_options$Bpa$ssb_flash <-  propagate(catch_options$Bpa$ssb_flash, nits)
catch_options$Blim$ssb_flash <-  propagate(catch_options$Blim$ssb_flash, nits)
catch_options$MSYBtrigger$ssb_flash <-  propagate(catch_options$MSYBtrigger$ssb_flash, nits)
catch_options$Fsq$fbar <-  propagate(catch_options$Fsq$fbar, nits)
catch_options$SBsq$ssb_flash  <-  propagate(catch_options$SBsq$ssb_flash, nits)
catch_options$rotac$catch <-  propagate(catch_options$rotac$catch, nits)

# CONVERT to fwdControl

make_fc <- function(target, itarget, F0){
  
  qname <- names(target)[1]
  nits <- dim(itarget[[1]])[6]
  
  fc <- fwdControl(
    data.frame(
      year = c(ay, fy, fy + 1),
      quant = c(
        names(itarget)[1],
        names(target)[1],
        names(F0)[1]
      ),
      fishery = NA,
      catch = NA,
      biol = 1
    )
  )
  
  iters <- array(
    NA,
    dim = c(3, nits, 3),
    dimnames = list(
      step = 1:3,
      iter = 1:nits,
      stat = c("min","value","max")
    )
  )
  
  iters[1,, "value"] <- c(itarget[[1]])
  iters[2,, "value"] <- c(target[[1]])
  iters[3,, "value"] <- c(F0[[1]])
  
  fc@iters <- aperm(iters, c(1,3,2))
  
  fc
}

fctls <- lapply(
  catch_options,
  make_fc,
  itarget = window(itarget, end = ay),
  F0 = F0
)


# RUN!

runs <- FLStocks(lapply(fctls, function(x) fwd(fut, sr=gmsrr, control=x)))

# Sustituir SSB por F

flevels <- seq(0.05, 0.80, 0.01)

## Curva F -> P(SSB < Blim) ----
risk_blim <- sapply(flevels, function(ff){
  
  ctrl <- fctls$Fmsy
  
  ctrl@iters["2", "value", ] <- ff
  
  run <- fwd(fut, sr=gmsrr, control=ctrl)
  
  mean(
    c(ssb(run)[, ac(fy+1)]) < refPts$Blim
  )
  
})

plot(
  flevels,
  100 * risk_blim,
  type = "l",
  lwd = 2,
  xlab = "F",
  ylab = "% SSB < Blim"
)

abline(h = 5, col = "red", lwd = 2, lty = 2)

legend(
  "topleft",
  legend = c(
    "5% risk threshold"
  ),
  col = c("red"),
  lty = c(2),
  lwd = 2,
  bty = "n"
)

SavePlotBlim('BlimProb_F')

# # F que produce un 50% de riesgo
# F_Blim50 <- approx(
#   x = risk_blim,
#   y = flevels,
#   xout = 0.5
# )$y

# F que produce un 5% de riesgo
F_Blim05 <- approx(
  x = risk_blim,
  y = flevels,
  xout = 0.05
)$y

## Curva F -> SSB ----

ssb_risk <- sapply(flevels, function(ff){
  
  ctrl <- fctls$Fmsy
  
  ctrl@iters["2","value",] <- ff
  
  run <- fwd(fut, sr=gmsrr, control=ctrl)
  
  median(
    c(ssb(run)[, ac(fy+1)])
  )
})

plot(
  flevels,
  ssb_risk,
  type="l",
  lwd = 2,
  xlab="F",
  ylab="Median SSB"
)

# horizontales
abline(h = refPts$Blim,
       col = "red",
       lty = 2,
       lwd = 2)

abline(h = refPts$MSYBtrigger,
       col = "blue",
       lty = 2,
       lwd = 2)

SSBsq_target <- median(
  c(ssb(runs$Fmsy)[, ac(fy)])
)

abline(h = SSBsq_target,
       col = "darkgreen",
       lty = 2,
       lwd = 2)

legend(
  "bottomleft",
  legend = c(
    paste0("Blim = ", round(refPts$Blim)),
    paste0("MSY Btrigger = Bpa = ", round(refPts$MSYBtrigger)),
    paste0("SSB(", fy, ") = ", round(SSBsq_target))
  ),
  col = c("red", "blue", "darkgreen"),
  lty = 2,
  lwd = 2,
  bty = "n"
)

SavePlotBlim('SSB_F')

# F que da SSB = Blim
F_Blim <- approx(
  x = ssb_risk,
  y = flevels,
  xout = refPts$Blim
)$y

# F que da SSB = MSYBtrigger
F_MSYBtrigger <- approx(
  x = ssb_risk,
  y = flevels,
  xout = refPts$MSYBtrigger
)$y

# F que da SSB(2028)=SSB(2027)
SSBsq_target <- median(
  c(ssb(runs$Fmsy)[, ac(fy)])
)
F_SSBsq <- approx(
  x = ssb_risk,
  y = flevels,
  xout = SSBsq_target
)$y


### F_Blim05, F_Blim, F_MSYBtrigger, F_SSBsq

fctls$Blim05 <- fctls$Blim <- fctls$MSYBtrigger <- fctls$Bpa <- fctls$SBsq <- fctls$Fmsy

fctls$Blim05@iters["2","value",] <- F_Blim05
fctls$Blim@iters["2","value",] <- F_Blim
fctls$MSYBtrigger@iters["2","value",] <- fctls$Bpa@iters["2","value",] <- F_MSYBtrigger
fctls$SBsq@iters["2","value",] <- F_SSBsq

fctls$FmsyNew <- fctls$Fmsy
fctls$FmsyNew@iters["2","value",] <- fctls$FmsyNew@iters["2","value",] * SSBsq_target/refPts$MSYBtrigger

fctls$lFmsyNew <- fctls$lFmsy
fctls$lFmsyNew@iters["2","value",] <- fctls$lFmsyNew@iters["2","value",] * SSBsq_target/refPts$MSYBtrigger

fctls$uFmsyNew <- fctls$uFmsy
fctls$uFmsyNew@iters["2","value",] <- fctls$uFmsyNew@iters["2","value",] * SSBsq_target/refPts$MSYBtrigger

runs <- FLStocks(lapply(fctls, function(x) fwd(fut, sr=gmsrr, control=x)))

# SSB < Blim

prob_blim <- round(
  100 * sapply(names(runs), function(x) {
    mean(c(ssb(runs[[x]])[, ac(fy+1)]) < refPts$Blim)
  }),
  1
)

prob_blim # Probability SSB < Blim considering all the MCMC iterations.

prob_blim_90 <- round(
  100 * sapply(names(runs), function(x) {
    
    ssb_vals <- c(ssb(runs[[x]])[, ac(fy + 1)])
    
    q <- quantile(ssb_vals, c(0.05, 0.95))
    
    mean(ssb_vals[ssb_vals >= q[1] &
                    ssb_vals <= q[2]] < refPts$Blim)
    
  }),
  1
)

prob_blim_90 # Probability SSB < Blim considering the central 90% of MCMC iterations.

# COMPARE

Map(compare, runs, fctls)

# Contribution of each age to the landing and SSB
par(mfrow=c(1,2), mar=c(5,8,4,1), cex=0.8)

# Mediana de landings.n por edad
yield <- apply(
  as.array(runs$Fmsy@landings.n[, nyears+2]),
  1,
  median
)

yield <- drop(yield)

prop <- paste0(round(100 * yield / sum(yield)), "%")
labels <- paste0(
  max(years) - ages + 2,
  rep(c(" (GM)", " (a4a)"), c(2, nages-2))
)

b <- barplot(
  yield,
  horiz = TRUE,
  names = labels,
  las = 1,
  xlab = "Tonnes",
  main = paste("Landings yield", max(years)+2, "(median MCMC)"),
  xlim = c(0, max(yield) * 1.25)
)

text(yield, b, prop, adj = -0.2)
mtext("Cohort", 2, 6)

# ssb
ssb_age <- runs$Fmsy@stock.n *
  runs$Fmsy@stock.wt *
  runs$Fmsy@mat

ssb <- apply(
  as.array(ssb_age[, nyears+3]),
  MARGIN = 1,
  median
)

ssb <- drop(ssb)

prop <- paste0(round(100*ssb/sum(ssb)),'%')
labels <- max(years)-ages+3
b <- barplot(ssb,horiz=T,names=labels,las=1,xlab='Tonnes',main=paste('SSB',max(years)+3),xlim=c(0,max(ssb)*1.25))
text(ssb,b,prop,adj=-0.2)
mtext('Cohort',2,4)

SavePlot('Contrib')

# SAVE 
save(runs, prob_blim, prob_blim_90, risk_blim, ssb_risk, recsq, tac, advice, refPts, file=paste0("model/runs_", Fassumption, "_", Rassumption, ".RData"))
stk_est_stf <- runs[["Fmsy"]]
save(stk_est_stf, file = "model/Meg78_WGMIX_stock_estimated_stf_alliters.Rdata")
stk_est_stf_median <- qapply(stk_est_stf, iterMedians)
save(stk_est_stf_median, file = "model/Meg78_WGMIX_stock_estimated_stf_median.Rdata")

## Summary table  ----
fctls_fmsy <- fctls$Fmsy
fctls_fmsy$value[fctls_fmsy$year == ac(fy+1)] <- fctls_fmsy$value[fctls_fmsy$year == ac(fy)]

run_fmsy <- fwd(fut, sr = gmsrr, control = fctls_fmsy)

# summary plot
new_names <- c("Rec" = "Recruitment", "SB" = "SSB", "C" = "Catch", "F" = "F")
flqs <- FLQuants(sim=iterMedians(stock.n(fits)), det=stock.n(fit01))
keylst <- list(points=FALSE, lines=TRUE, space="right")
plot(run_fmsy) + facet_wrap(~qname,scales='free_y', labeller = labeller(qname = new_names)) +theme_bw()
SavePlot('Blim_forecast')

stks <- run_fmsy

# Summary table using %90 confidence interval
# Lognormal-based confidence intervals (90%) for key variables

# Valores medios
tsb <- apply(tsb(stks), 2, median)
ssb <- apply(ssb(stks), 2, median)
catchobs <- stock@catch
discardsobs <- stock@discards
catch <- apply(stks@catch, 2, median)
recr <- apply(stks@stock.n[1,], 2, median)
fbar <- apply(fbar(stks), 2, median)

# Parámetro z para 90% CI
z <- 1.6449

# Desviaciones estándar (estimadas a partir de las iteraciones en stks)
tsbse <- apply(tsb(stks), 2, sd)
ssbse <- apply(ssb(stks), 2, sd)
recrse <- apply(stks@stock.n[1,], 2, sd)
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

# Crear tabla resumen final
summary_df <- data.frame(
  year = 1984:(max(years)+3),
  recr_lo = c(recr_CI$lower),
  recr = recr_mean,
  recr_hi = c(recr_CI$upper),
  tsb_lo = c(tsb_CI$lower),
  tsb = tsb_mean,
  tsb_hi = c(tsb_CI$upper),
  ssb_lo = c(ssb_CI$lower),
  ssb = ssb_mean,
  ssb_hi = c(ssb_CI$upper),
  fbar_lo = c(fbar_CI$lower),
  fbar = fbar_mean,
  fbar_hi = c(fbar_CI$upper),
  catch=c(catchobs, NA, NA, NA),
  lan=c(catchobs-discardsobs, NA, NA, NA),
  dis=c(discardsobs, NA, NA, NA)
)
summary_df <- as.data.table(summary_df)
summary_df[, c(2:10,14:16) := lapply(.SD, round, digits=0), .SDcols = c(2:10,14:16)]
summary_df[, c(11:13) := lapply(.SD, icesRound), .SDcols = c(11:13)]

# Exportar CSV
write.csv(summary_df, "report/Summary_SAG_lognorm_fwd.csv", row.names = FALSE)
