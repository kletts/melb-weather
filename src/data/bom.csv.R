
library(tidyverse)
library(jsonlite)
library(RCurl)
library(dlm)

stations <- c("Avalon"="94854", 
              "Essendon"="95866", 
              "Frankston"="94876", 
              "Viewbank"="95874")

get_station_data <- function(id) { 
  url <- paste0("http://www.bom.gov.au/fwo/IDV60901/IDV60901.", id, ".json")
  rawdata <- RCurl::getURL(url, ssl.verifypeer = FALSE)
  data <- jsonlite::fromJSON(rawdata)$observations$data |> 
      dplyr::mutate(local_date_time_full = as.POSIXct(local_date_time_full, 
                                                      format = "%Y%m%d%H%M%S", 
                                                      tz="Australia/Melbourne")) } 

data <- purrr::imap_dfr(stations, 
    \(x,y) get_station_data(x) |> dplyr::mutate(name=y)) 

metrics <- c("air_temp"="Temperature (C)", 
          "wind_spd_kmh"="Wind Speed (km/h)", 
          "press"="Pressure (hPa)", 
          "rel_hum"="Relative Humidity (%)")

obstime <- seq(min(data$local_date_time_full), 
               max(data$local_date_time_full), 
               by="30 mins", 
               tz="Australia/Melbourne")

data <- data |> 
  tidyr::pivot_longer(cols=names(metrics), names_to="metric", values_to="value") |> 
  dplyr::reframe(data.frame(obstime=obstime, 
                     valuereg=approx(x= local_date_time_full,
                                     y =value, 
                                     xout = obstime)$y), 
          .by=c("name", "metric")) |> 
  tidyr::pivot_wider(id_cols=c(obstime, metric), names_from="name", values_from="valuereg")

# Multivariate signal -----
fn <- function(parm) { 
  # Optimal parameter update function
  mod <- dlm(
    FF=matrix(c(1,1,1,1,1,1,1,1,0,0,0,0), ncol=3, byrow=FALSE), 
    V=diag(exp(rep(parm[1], 4))), 
    GG=(dlm::dlmModPoly(1) + dlm::dlmModTrig(s=24*2, q=1))$GG, 
    W=diag(exp(parm[2:4])), 
    m0=c(0,0,0), 
    C0=diag(c(1,1,1)))
  return(mod) } 

extract_signal <- function(data, metric) {  
  # Apply multivariate optimisation and trend filter to each metric
  data.m <- data |> 
    data.frame(row.names = "obstime") |>
    data.matrix()
  fit <- dlm::dlmMLE(data.m, parm=log(c(1,1,1,1)), build = fn, hessian = TRUE)
  fitf <- dlm::dlmFilter(data.m, fn(fit$par))
  fitf.conf <- dlm::dlmSvd2var(fitf$U.C, fitf$D.C) 
  data |> 
    dplyr::bind_cols(data.frame(dropFirst(fitf$m))) |> 
    dplyr::rename(Trend=X1, Cycle=X2) |> 
    dplyr::select(-X3) |> 
    dplyr::bind_cols(purrr::map_dfr(fitf.conf[-1], \(x) setNames(diag(x)[1:2], c("TrendVar", "CycleVar"))))
  } 

# Apply and write output to CSV ---- 
data |> 
  dplyr::nest_by(metric) |> 
  dplyr::mutate(trend = list(extract_signal(data, metric))) |> 
  tidyr::unnest(trend) |> 
  dplyr::select(-data) |> 
  dplyr::mutate(metric = metrics[metric]) |> 
  readr::format_csv() |> 
  cat()

