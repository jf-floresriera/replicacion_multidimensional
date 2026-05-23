# analisis.R  –  Script de referencia (Capítulo 21)
# Autor: Jesus Enrique Flores Riera
# Fecha: 22 de mayo de 2026

library(stars)
library(sf)
library(dplyr)
library(ggplot2)
library(readxl)
library(stringr)

# ── PARÁMETROS ──────────────────────────────────────────────────────────────
data_path <- '/home/rstudio/work/01progsig/data/heavy/'
mri_query <- c(10, 20, 50, 100, 250, 500, 700, 1000, 1700, 3000, 7000)

# ── 1. INGESTA DE RASTERS (Huracanes) ───────────────────────────────────────
h_files <- list.files(data_path, pattern = '^h.*\.tif$', full.names = TRUE)
h_files <- stringr::str_sort(h_files, numeric = TRUE)
sth <- read_stars(h_files, quiet = TRUE)
mynames_h <- paste0('h', mri_query)
sth <- setNames(sth, mynames_h)
print(sth)

# ── 2. INGESTA TABULAR (No Huracanes) ───────────────────────────────────────
nhrl <- readxl::read_excel(file.path(data_path, 'nh_returnlevels.xlsx'))
print(head(nhrl))

# ── 3. CONSTRUCCIÓN DEL CUBO VECTORIAL ──────────────────────────────────────
# Convierte malla raster a polígonos sf
sfh <- st_xy2sfc(sth, as.points = FALSE, na.rm = TRUE)

# Une polígonos con tabla de niveles de retorno NH
sf1    <- st_as_sf(sth['h10'])  # malla base
sfnh   <- left_join(sf1, nhrl, by = 'id')
stnh   <- st_as_stars(sfnh)

# ── 4. FUNCIÓN DE EXTRAPOLACIÓN LINEAL ──────────────────────────────────────
extrapola <- function(x, y, xout) {
  if (any(duplicated(x))) {
    datos <- aggregate(y ~ x, data = data.frame(x = x, y = y), FUN = mean)
    x <- datos$x; y <- datos$y
  }
  ord <- order(x); x <- x[ord]; y <- y[ord]; n <- length(x)
  yout <- suppressWarnings(approx(x, y, xout = xout)$y)
  idx_inf <- which(xout < x[1])
  idx_sup <- which(xout > x[n])
  if (length(idx_inf) > 0) {
    slope_low <- (y[2] - y[1]) / (x[2] - x[1])
    yout[idx_inf] <- y[1] + slope_low * (xout[idx_inf] - x[1])
  }
  if (length(idx_sup) > 0) {
    slope_high <- (y[n] - y[n-1]) / (x[n] - x[n-1])
    yout[idx_sup] <- y[n] + slope_high * (xout[idx_sup] - x[n])
  }
  return(yout)
}

# ── 5. FUNCIÓN COMBINACIÓN STAPPLY ──────────────────────────────────────────
# pc = 1 - (1 - ph)(1 - pnh)
combined_columns_stapply <- function(st,
                                     bandsrlh, mrih,
                                     bandsrlnh, mrinh,
                                     mric) {
  rl600 <- seq(1, 600, by = 1)
  # Huracanes
  rlh <- NULL
  for (band in bandsrlh) rlh <- c(rlh, st[band])
  if (any(!is.na(rlh)) && sum(rlh, na.rm=TRUE) > 0) {
    h    <- data.frame(rl = rlh, prob = 1 / mrih)
    hch_prob <- extrapola(x = h$rl, y = h$prob, xout = rl600)
    hch_prob[hch_prob < 0] <- 0; hch_prob[hch_prob > 1] <- 1
  } else {
    hch_prob <- rep(NA, length(rl600))
  }
  # No huracanes
  rlnh <- NULL
  for (band in bandsrlnh) rlnh <- c(rlnh, st[band])
  nh    <- data.frame(rl = rlnh, prob = 1 / mrinh)
  hcnh_prob <- extrapola(x = nh$rl, y = nh$prob, xout = rl600)
  hcnh_prob[hcnh_prob < 0] <- 0; hcnh_prob[hcnh_prob > 1] <- 1
  # Combinación
  if (any(is.na(hch_prob))) {
    hcc_prob <- hcnh_prob
  } else {
    hcc_prob <- 1 - (1 - hch_prob) * (1 - hcnh_prob)
  }
  hcc_prob[hcc_prob < 0] <- 0
  prob_c <- 1 / mric
  hcc  <- data.frame(rl = rl600, prob = hcc_prob)
  rlc_output <- extrapola(x = hcc$prob, y = hcc$rl, xout = prob_c)
  return(rlc_output)
}

# ── 6. VISUALIZACIÓN RESULTADO T=700 ────────────────────────────────────────
# (requiere stnhhcmerged construido previamente en el flujo completo)
# Ejemplo de visualización con el cubo de huracanes
plot(sth['h700'], col = viridis::plasma(100), main = 'Viento huracanado T=700 años',
     reset = FALSE, border = NA)

message('Script analisis.R completado.')
