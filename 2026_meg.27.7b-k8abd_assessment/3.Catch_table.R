#############################################################################################
# catch_table.R - DESC
# /catch_table.R

# Copyright Iago MOSQUEIRA (WMR), 2021
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

### Clean slate
rm(list=ls())

library(icesAdvice)
library(flextable)
library(data.table)
library(writexl)
library(FLFishery)

setwd("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2026/3.Assessment_october")
set.seed(1234)

SavePlot <- function(plotname, width = 10, height = 7) {
  file <- paste0(getwd(), "/report/", "meg78_probBlim_", plotname, ".png")
  dev.print(png, file, width = width, height = height, units = "in", res = 300)
}

# setting
Fassumption <- "Funscaled" #"Fscaled" or "Funscaled"
Rassumption <- "gmean" # "gmean" or "meanblowmean" or "lastwoyr"

load(paste0("model/runs_", Fassumption, "_",Rassumption,".RData"))

forecastYr <- range(runs[[1]])["maxyear"]-1

discardFages <- 1:2  #is this correct?
landFages <- range(runs[[1]])["minfbar"]:range(runs[[1]])["maxfbar"] #It is 3-6.

metrics(window(runs[[1]], start=forecastYr, end=forecastYr),
        list(catch=catch, wanted=landings, unwanted=discards,
             F=fbar))

unlist(lapply(runs, is, "FLStock"))

# TABLE. Advice sheet annual catch options

# tab_options <- lapply(runs, function(x) {
#   
#   data.frame(
#     catch    = median(c(catch(window(x, start=forecastYr, end=forecastYr)))),
#     wanted   = median(c(landings(window(x, start=forecastYr, end=forecastYr)))),
#     unwanted = median(c(discards(window(x, start=forecastYr, end=forecastYr)))),
#     F        = median(c(fbar(window(x, start=forecastYr, end=forecastYr)))),
#     SSB      = median(c(ssb(window(x, start=forecastYr+1, end=forecastYr+1))))
#   )
#   
# })

tab_options <- lapply(runs, function(x) {
  
  ssb_vals <- c(
    ssb(window(x,
               start = forecastYr + 1,
               end   = forecastYr + 1))
  )
  
  data.frame(
    catch    = median(c(catch(window(x, start=forecastYr, end=forecastYr)))),
    wanted   = median(c(landings(window(x, start=forecastYr, end=forecastYr)))),
    unwanted = median(c(discards(window(x, start=forecastYr, end=forecastYr)))),
    F        = median(c(fbar(window(x, start=forecastYr, end=forecastYr)))),
    SSB      = median(ssb_vals),
    
    SSB_LCI  = as.numeric(quantile(ssb_vals, 0.05)),
    SSB_UCI  = as.numeric(quantile(ssb_vals, 0.95))
  )
  
})

tab_options <- rbindlist(tab_options, idcol="basis")

# ADD SSB change, (old - new / new) * 100
ssbiy <- median(
  c(ssb(runs[[1]])[, as.character(forecastYr)])
)
tab_options[, ssbchange:=(SSB - ssbiy) / ssbiy * 100]

# ADD TAC change
tab_options[, tacchange:=(catch - c(tac)) / c(tac) * 100]

# ADD advice change
tab_options[, advicechange:=(catch - c(advice)) / c(advice) * 100]

# ADD SSB > Blim
tab_options[, prob_blim :=prob_blim ]

# FIX zeroes
tab_options[, (2:ncol(tab_options)) := lapply(.SD, function(x) ifelse(abs(x) < 0.0001, 0, x)), .SDcols=2:ncol(tab_options)]

# Remove Roll-over TAC scenario
tab_options <- tab_options[basis != "rotac"]

# SET first column
# tab_options$basis <- c("MSY approach: F[MSY]", "F=MAP F[MSY lower]",
#                        "F=MAP F[MSY upper]", "F=0", "F[pa]", "F[lim]", 
#                        paste0("SSB (",forecastYr+1,")=B[pa]"),
#                        paste0("SSB(",forecastYr+1,")=B[lim]"), 
#                        paste0("SSB(",forecastYr+1,")=MSY B[trigger]"),
#                        paste0("SSB(",forecastYr+1,")=SSB(",forecastYr,")"),
#                        "F[2020]", "Roll-over TAC")

tab_options$basis <- c("MSY approach: F[MSY]", "F=MAP F[MSY lower]",
                       "F=MAP F[MSY upper]", "F=0", "F[pa]", 
                       paste0("SSB (",forecastYr+1,")=B[pa]"),
                       paste0("SSB(",forecastYr+1,")=B[lim]"), 
                       paste0("SSB(",forecastYr+1,")=MSY B[trigger]"),
                       "F[2026]",                                             #"F[2023]" hay que cambiarlo manualmente, en el WGBIE 2023 ponemos F(2023).
                       paste0("SSB(",forecastYr+1,")=SSB(",forecastYr,")"),
                       # "Roll-over TAC")  
                       "F giving 5% risk of SSB < Blim",
                       paste0("F[MSY] × SSB",forecastYr+1,"/MSYB[trigger]"),
                       paste0("F[MSY lower] × SSB",forecastYr+1,"/MSYB[trigger]"),
                       paste0("F[MSY upper] × SSB",forecastYr+1,"/MSYB[trigger]"))
                       

# CALL round and icesRound
tab_options[, c("catch","wanted","unwanted","SSB","SSB_LCI","SSB_UCI", "prob_blim") :=
              lapply(.SD, round, digits = 0),
            .SDcols = c("catch","wanted","unwanted","SSB","SSB_LCI","SSB_UCI", "prob_blim")]

tab_options[, "F" := round(F, 3)]

tab_options[, c("ssbchange","tacchange","advicechange") :=
              lapply(.SD, icesRound),
            .SDcols = c("ssbchange","tacchange","advicechange")]

## ICES Rounding for F
x <- as.logical(tab_options[,"F"]>=0.2)
tab_options[x,"F"] <- round(tab_options[x,"F"],2)

# CREATE table
ft <- flextable(tab_options)

# SET colnames
ftabops <- set_header_labels(ft, 
                             basis="Basis", 
                             catch=paste0("Total catch (",forecastYr,")"), 
                             wanted=paste0("Wanted catch (",forecastYr,")"),
                             unwanted=paste0("Unwanted catch (",forecastYr,")"), 
                             F=paste0("F[total] (ages ",min(landFages),"-",max(landFages),") (",forecastYr,")"),
                             SSB=paste0("SSB (",forecastYr+1,")"), 
                             SSB_LCI="SSB P5",
                             SSB_UCI="SSB P95",
                             ssbchange="% SSB change", 
                             tacchange="% TAC change",
                             advicechange="% Advice change",
                             prob_blim="% SSB < Blim")

dimnames(tab_options)[[2]] <- c("Basis", 
                               paste0("Total catch (",forecastYr,")"), 
                               paste0("Wanted catch (",forecastYr,")"),
                               paste0("Unwanted catch (",forecastYr,")"), 
                               paste0("F[total] (ages ",min(landFages),"-",max(landFages),") (",forecastYr,")"),
                               paste0("SSB (",forecastYr+1,")"), 
                               paste0("SSB P5 (",forecastYr+1,")"),
                               paste0("SSB P95 (",forecastYr+1,")"),
                               "% SSB change", "% TAC change", "% Advice change",
                               "% SSB < Blim")

tab_options$Basis <- ifelse(tab_options$Basis== "SSB (2028)=B[pa]", "SSB (2028)=B[pa]=MSY B[trigger]",  tab_options$Basis)
tab_options <- tab_options[tab_options$Basis != "SSB(2028)=MSY B[trigger]",]

# Order the table
# Define the custom order
custom_order <- c("MSY approach: F[MSY]", "F=MAP F[MSY lower]", "F=MAP F[MSY upper]", "F=0", "F[pa]", "SSB(2028)=B[lim]",
                  "SSB (2028)=B[pa]=MSY B[trigger]", "SSB(2028)=SSB(2027)", "F[2026]", "F giving 5% risk of SSB < Blim",
                  paste0("F[MSY] × SSB",forecastYr+1,"/MSYB[trigger]"),
                  paste0("F[MSY lower] × SSB",forecastYr+1,"/MSYB[trigger]"),
                  paste0("F[MSY upper] × SSB",forecastYr+1,"/MSYB[trigger]")) # "Roll-over TAC")

tab_options$Basis <- factor(tab_options$Basis, levels = custom_order)
tab_options <- tab_options[order(tab_options$Basis), ]

writexl::write_xlsx(tab_options,
                    path=paste0(getwd(),"/report/catch_options_", Fassumption, "_",Rassumption,".xlsx"))


### Blim Plots  ----

## Distribución de SSB para Fmsy

# Las medianas son muy parecidas, la diferencia del riesgo proviene de la cola inferior de la distribución.

par(mfrow=c(1,1))

hist(c(ssb(runs$Fmsy)[, ac(forecastYr+1)]),
     breaks=40,
     freq=FALSE,
     col=rgb(1,0,0,0.4),
     xlab="SSB",
     main="Projected SSB")

# hist(c(ssb(runs$rotac)[, ac(forecastYr+1)]),
#      breaks=40,
#      freq=FALSE,
#      col=rgb(0,0,1,0.4),
#      add=TRUE)

abline(v=refPts$Blim, lwd=2, lty=2)

legend("topright",
       # c("Fmsy","ROTAC","Blim"),
       c("Fmsy","Blim"),
       # fill=c(rgb(1,0,0,0.4),rgb(0,0,1,0.4),NA),
       fill=c(rgb(1,0,0,0.4),NA),
       border=c(NA,NA),
       lty=c(NA,2))

SavePlot("distribution")

# CDF
# ¿Cómo puede ser que Fmsy capture más y tenga menos riesgo?
#   La diferencia está únicamente en esta parte de la cola.

vals_fmsy  <- c(ssb(runs$Fmsy)[, ac(forecastYr+1)])
vals_rotac <- c(ssb(runs$rotac)[, ac(forecastYr+1)])

plot(ecdf(vals_fmsy),
     col="red",
     lwd=2,
     main="Probability distribution of SSB",
     xlab="SSB",
     ylab="Cumulative probability")

lines(ecdf(vals_rotac),
      col="blue",
      lwd=2)

abline(v=refPts$Blim,
       lty=2,
       lwd=2)

legend("topleft",
       c("Fmsy","ROTAC"),
       col=c("red","blue"),
       lwd=2)

SavePlot("CDF")

# Boxplot de escenarios
boxplot(
  list(
    Fmsy  = c(ssb(runs$Fmsy)[, ac(forecastYr+1)]),
    Rotac = c(ssb(runs$rotac)[, ac(forecastYr+1)]),
    F0    = c(ssb(runs$F0)[, ac(forecastYr+1)])
  ),
  ylab="SSB"
)

abline(h=refPts$Blim,
       col="red",
       lwd=2,
       lty=2)

SavePlot("boxplot")

# Risk vs catch
risk <- sapply(names(runs), function(x)
  mean(c(ssb(runs[[x]])[, ac(forecastYr+1)]) < refPts$Blim)
)

catch_med <- sapply(names(runs), function(x)
  median(c(catch(runs[[x]])[, ac(forecastYr)]))
)

catch_med <- catch_med[names(catch_med) != "rotac"]
risk <- risk[names(risk) != "rotac"]

plot(catch_med,
     100 * risk,
     pch = 19,
     xlab = paste0("Median catch (", forecastYr, ", t)"),
     ylab = "Probability SSB < Blim (%)")

labs <- names(catch_med)

x <- catch_med
y <- 100 * risk

# Desplazamientos manuales
x[labs == "Bpa"]         <- x[labs == "Bpa"] - 400
y[labs == "Bpa"]         <- y[labs == "Bpa"] + 1

x[labs == "MSYBtrigger"] <- x[labs == "MSYBtrigger"] + 400
y[labs == "MSYBtrigger"] <- y[labs == "MSYBtrigger"] - 1

x[labs == "Fmsy"]        <- x[labs == "Fmsy"] + 300
y[labs == "Fmsy"]        <- y[labs == "Fmsy"] + 1

x[labs == "uFmsy"]       <- x[labs == "uFmsy"] + 300
y[labs == "uFmsy"]       <- y[labs == "uFmsy"] - 1

x[labs == "Fpa"]         <- x[labs == "Fpa"] - 300
y[labs == "Fpa"]         <- y[labs == "Fpa"] + 2

x[labs == "FmsyNew"] <- x[labs == "FmsyNew"] + 400
y[labs == "FmsyNew"] <- y[labs == "FmsyNew"] - 1

x[labs == "uFmsyNew"]        <- x[labs == "uFmsyNew"] + 300
y[labs == "uFmsyNew"]        <- y[labs == "uFmsyNew"] + 1

text(x, y, labels = labs)

# Líneas de unión para las etiquetas desplazadas
idx <- labs %in% c("Bpa","MSYBtrigger","Fmsy","uFmsy","Fpa")

segments(
  x0 = x[idx],
  y0 = y[idx],
  x1 = catch_med[idx],
  y1 = (100 * risk)[idx],
  col = "grey50",
  lty = 2
)

SavePlot("riskVScatch")

# Observed catch in the last year vs catch option table advices
load("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2026/3.Assessment_october/model/MegFit_Final_october.Rdata")

tmp <- copy(tab_options)
tmp <- tmp[order(`Total catch (2027)`)]

last_catch <- median(
  c(catch(stock)[, "2025"])
)

advice2025 <- 21144    # o la cifra que uses oficialmente
advice2026 <- as.numeric(advice)

op <- par(mar = c(5, 18, 2, 2))

dotchart(
  tmp$`Total catch (2027)`,
  labels = tmp$Basis,
  pch = 19,
  xlab = "Median catch (t)",
  xlim = c(0, max(tmp$`Total catch (2027)`) * 1.1)
)

# Observed catch 2025
abline(
  v = last_catch,
  col = "black",
  lwd = 2,
  lty = 1
)

# Advice for 2025
abline(
  v = advice2025,
  col = "darkgreen",
  lwd = 2,
  lty = 2
)

# Advice for 2026
abline(
  v = advice2026,
  col = "blue",
  lwd = 2,
  lty = 3
)

legend(
  "bottomright",
  legend = c(
    paste0("Observed catch 2025 = ", round(last_catch), " t"),
    paste0("Advice 2025 = ", round(advice2025), " t"),
    paste0("Advice 2026 = ", round(advice2026), " t")
  ),
  col = c("black", "darkgreen", "blue"),
  lty = c(1, 2, 3),
  lwd = 2,
  bty = "n"
)

par(op)

SavePlot("ObsCatchVScatchOpt")