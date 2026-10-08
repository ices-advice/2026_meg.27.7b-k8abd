#' Abundance indices
#'
#' Abundance indices used as an FLIndex object
#'
#' @name indices
#' @format RData file
#' @tafOriginator WGBIE
#' @tafYear 2026
#' @tafAccess Public
#' @tafSource script

library(icesTAF)
library(xlsx)
library(dplyr)
library(FLCore)

setwd("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2026/3.Assessment_october")

# utility functions
get_sheet <- function(sheet) {
  read.xlsx2(
    taf.boot.path("initial", "data", "Inputdata_MEGRIM_a4a_WGBIE2026_afterChecking.xlsx"),
    sheetName = sheet,
    check.names = FALSE,
    colClasses = numeric()
  ) %>%
    mutate_all(function(x) ifelse(x == "NA", NA, x)) %>%
    mutate(across(where(is.character), as.numeric))
}

get_age_quant <- function(sheet, units = "") {
  get_sheet(sheet) %>%
    taf2long(names = c("year", "age", "data")) %>%
    mutate(data = ifelse(data == 0, NA, data)) %>%
    as.FLQuant(units = units)
}

get_index <- function(sheet, name, desc, frange) {
  index <-
    FLIndex(
      name = name,
      desc = desc,
      type = "numbers",
      index = get_age_quant(sheet)
    )
  range(index, c("startf", "endf")) <- frange
  index
}

get_bioindex <- function(sheet, name, desc, frange) {
  index <-
    FLIndexBiomass(
      name = name,
      desc = desc,
      index = get_age_quant(sheet)
    )
  range(index, c("startf", "endf")) <- frange
  index
}


####SURVEYS############################

# EVHOE
evhoe <-
  get_index(
    sheet = "CPUE.EVHOE",
    name = "FR_EVHOE",
    desc = "French IBTS survey index in 7 and 8abd; Catch in numbers per hour;L. whiffiagonnis",
    c(0.75, 1.00) #Q4
  )

# PORCUPINE
sp_porc <-
  get_index(
    sheet = "CPUE.SPGFS",
    name = "SP_PORC",
    desc = "Spanish IBTS Porcupine survey; Cpue in numbers per 30min; L. whiffiagonnis",
    c(0.50, 0.75) #Q3
  )

# IRLFRsurvey
irlfrsurvey <-
  get_index(
    sheet = "CPUE.IRLFRsurvey",
    name = "CPUE.IRLFRsurvey",
    desc = "Combined Irish and French IBTS survey Q4 ; unit Standardised to N*1000/hour; L. whiffiagonnis",
    c(0.75, 1.00)#Q4
  )


# IAM
iam <-
  get_index(
    sheet = "CPUE.IAM",
    name = "CPUE.IAM",
    desc = "Irish Anglerfish and Megrim Survey in Q1 since 2016; unit Standardised to N*1000/hour; L. whiffiagonnis",
    c(0.0, 0.25) #Q1
  )


####COMMERCIAL FLEET############################

# VIGO84-98
vigo8498 <-
  get_index(
    sheet = "CPUE.VIGO84",
    name = "CPUE.Vigo84",
    desc = "Spanish demersal trawlers (Vigo) in subarea 7 from 1984-1998; L. whiffiagonnis",
    c(0.0, 1.00)
  )

# VIGO99-2019
vigo9919 <-
  get_index(
    sheet = "CPUE.VIGO99",
    name = "CPUE.Vigo99",
    desc = "Spanish demersal trawlers (Vigo) in subarea 7 from 1999-2019; L. whiffiagonnis",
    c(0.0, 1.00)
  )

# IRTBB
irtbb <-
  get_index(
    sheet = "LPUE.IRTBB",
    name = "LPUE.IRTBB",
    desc = "Irish beam trawlers; unit Standardised to N0/10SqKm; L. whiffiagonnis",
    c(0.0, 1.00)
  )


# VIgoSampl ### NEW BIOMASS INDEX

# vigosamp <-
#   get_index(
#     sheet = "CPUE.VIGOSAMP",
#     name = "CPUE.VIGOSAMP",
#     desc = "Vigo Standardized CPUE from sampling data; unit Standardised to kg/hour; L. whiffiagonnis",
#     c(0.0, 1.00)
#   )


vigosamp <-
  get_bioindex(
    sheet = "CPUE.VIGOSAMP",
    name = "CPUE.VIGOSAMP",
    desc = "Vigo Standardized CPUE from sampling data; unit Standardised to kg/hour; L. whiffiagonnis",
    c(0.0, 1.00)
  )




######################################


# Create FLIndices
# Select the year and ages of each index
tun.sel <-
  FLIndices(
    evhoe,
    sp_porc,
    trim(iam,age=3:10),
    trim(irlfrsurvey,age=1:10),
    trim(vigo8498, age = 2:9),
    vigo9919,
    trim(irtbb, age = 2:7),
    vigosamp
  )

# Create plus group at 7 age
reduce_pg_index <- function(ix, pg = 7) {
  
  ages <- as.numeric(dimnames(ix@index)$age)
  maxage <- max(ages)
  
  # 1. Sumar edades pg:maxage en la edad pg
  ix@index[ac(pg), ] <- apply(ix@index[ac(pg:maxage), , drop = FALSE], 2, sum, na.rm = TRUE)
  
  # 2. Recortar edades
  ix <- trim(ix, age = 1:pg)
  
  # 3. Actualizar plusgroup
  ix@range["plusgroup"] <- pg
  
  return(ix)
}

# tun.sel[[1]] <- reduce_pg_index(tun.sel[[1]], 7)
tun.sel[[2]] <- reduce_pg_index(tun.sel[[2]], 7)
# tun.sel[[3]] <- reduce_pg_index(tun.sel[[3]], 7)
tun.sel[[4]] <- reduce_pg_index(tun.sel[[4]], 7)
# tun.sel[[5]] <- reduce_pg_index(tun.sel[[5]], 7)
# tun.sel[[6]] <- reduce_pg_index(tun.sel[[6]], 7)
# tun.sel[[7]] <- reduce_pg_index(tun.sel[[7]], 7)
# tun.sel[[8]] <- reduce_pg_index(tun.sel[[8]], 7)
tun.sel[[4]]@index[,"2017"] <- NA # IRLFR no data for 2017

save(tun.sel, file = "bootstrap/data/indices/meg78_indices.RData")

