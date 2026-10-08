#' Stock input data
#'
#' Input data related to the stock as an FLStock object
#'
#' @name stock
#' @format RData file
#' @tafOriginator WGBIE
#' @tafYear 2026
#' @tafAccess Public
#' @tafSource script

# AI (03/11/2020): Scrit updated by Ane Iriondo for WKTADSA 2020.

library(icesTAF)
library(xlsx)
library(dplyr)
library(FLCore)
library(openxlsx)
library(tidyr)

setwd("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2026/3.Assessment_october")

end.yr <- 2025

# utility functions
get_sheet <- function(sheet) {
  read.xlsx2(
    taf.boot.path("initial", "data", "Inputdata_MEGRIM_a4a_WGBIE2026_afterChecking.xlsx"),
    sheetName = sheet,
    check.names = FALSE,
    colClasses = numeric()
  )
}

get_age_quant <- function(sheet, units = NA) {
  get_sheet(sheet) %>%
    taf2long(names = c("year", "age", "data")) %>%
    as.FLQuant(units = units)
}

# collate data into an FLStock

# catch.n
catch.n <- get_age_quant("catches.n", "10^3")

# landings.n
landings.n <-
  get_age_quant("landings.n", "10^3") %>%
  window(start = 1984, end = end.yr)


# discards.n
discards.n <-
  get_sheet("discards.n") %>%
  mutate('9' = 0, '10' = 0) %>%
  taf2long(names = c("year", "age", "data")) %>%
  as.FLQuant(units = "10^3") %>%
  window(start = 1984, end= end.yr)


# catches
catch <-
  get_sheet("meg78_caton") %>%
  rename(data = catches) %>%
  select(year, data) %>%
  as.FLQuant(units = "t", quant = "age")

#landings
landings <-
  get_sheet("meg78_caton") %>%
  rename(data = landings) %>%
  select(year, data) %>%
  as.FLQuant(units = "t", quant = "age")

#discards
discards <-
  get_sheet("meg78_caton") %>%
  rename(data = discards) %>%
  select(year, data) %>%
  as.FLQuant(units = "t", quant = "age")

#catches.wt
catch.wt <- get_age_quant("catches.wt", "kg")

#landings.wt
landings.wt <- get_age_quant("landings.wt", "kg")

#discards.wt
discards.wt <- get_age_quant("discards.wt", "kg")


#stock.wt
df <- read.xlsx("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2026/0.Original_country_data/Accession&email/Surveys/Weight/Hans/MegIndexMeanLen.xlsx")

df2 <- df %>%
  mutate(MeanWeight = MeanWeight / 1000) # Pasar a kg

mean_0305 <- df2 %>%
  filter(Year %in% c(2003, 2004, 2005)) %>%
  group_by(Age) %>%
  summarise(MeanWeight = mean(MeanWeight, na.rm = TRUE),
            .groups = "drop") # Media por edad de los tres primeros años

hist_wt <- expand.grid(
  Age = unique(mean_0305$Age),
  Year = 1984:2002
) %>%
  left_join(mean_0305, by = "Age") # Expandir 1984-2002 usando dicha media # Indizea 2003tik aurrera.

wt_all <- bind_rows(
  hist_wt,
  df2 %>% select(Age, Year, MeanWeight) # Juntar históricos + observados
) %>%
  arrange(Age, Year)

wt.mat <- wt_all %>%
  pivot_wider(names_from = Year, values_from = MeanWeight) %>%
  arrange(Age) # Matriz edad x año

ages <- wt.mat$Age
years <- names(wt.mat)[-1]

wt.mat <- as.matrix(wt.mat[, -1])

arr <- array(
  wt.mat,
  dim = c(
    length(ages),
    length(years),
    1, 1, 1, 1
  ),
  dimnames = list(
    age    = as.character(ages),
    year   = as.character(years),
    unit   = "unique",
    season = "all",
    area   = "unique",
    iter   = "1"
  )
)

stock.wt <- FLQuant(arr)

units(stock.wt) <- "kg"

#mortality
m <- replace(stock.wt, TRUE, 0.2)
units(m) <- "m"

#maturity
mat <- get_age_quant("mat", "")

#harvest.spwn, harvest before spawning
harvest.spwn <- replace(mat, TRUE, 0)

#m.spwn, mortality before spawning
m.spwn <- replace(mat, TRUE, 0)

# FLSTOCK
stock <-
  FLStock(
    name = "meg78abd", desc = "run1",
    catch.n = catch.n, landings.n = landings.n, discards.n = discards.n,
    catch.wt = catch.wt, landings.wt = landings.wt, discards.wt = discards.wt, stock.wt = stock.wt,
    catch = catch, landings = landings, discards = discards,
    m = m, mat = mat,
    harvest.spwn = harvest.spwn, m.spwn = m.spwn
  )


range(stock, c("minfbar", "maxfbar")) <- c(3, 6)

# Create plus group at 7 age
stock <- setPlusGroup(stock, 7)

save(stock, file="bootstrap/data/stock/meg78_stock.RData")
