###############################################################################
# Title: UAP Spatial Analysis in France - Environmental & Sociological Drivers
# Description: Spatial Point Pattern Analysis of Unidentified Aerial Phenomena 
#              in France using Non-Stationary Poisson Processes (PPM).
# Authors (2026) : Michaël Vaillant
# Credits : Earlier prototype by Thibault Laurent, Christine Thomas-Agnan
###############################################################################

# UAP_Analysis.R
# Spatial analysis of GEIPAN UAP cases in Metropolitan France
# Outputs: PDF figures + LaTeX tables
#
# Run (Windows example):
# & "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" .\UAP_Analysis.R

# Packages used
library("png")
library("raster")
library("spdep")
library("spatstat")
library("sp")
library("sf") # Replaces maptools and rgdal
library("classInt")
library("RColorBrewer")
library("vcd")
library("xtable")
library("mnormt")
library("pixmap")

dir.create("cache", showWarnings = FALSE)
###############################################################################
#######           I. Exploratory Analysis
###############################################################################

# Create output directories if they do not exist
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

##############################################
# I.1. Importation of the data with the UAP
cas <- read.csv2("data/geipan/2014-09-19 cas_base_2014.csv", na.strings = "NULL")

# Define the Lambert 2 etendu CRS string
lambert2_proj <- "+proj=lcc +ellps=clrk80 +lat_1=45.898918964419 +lat_2=47.696014502038 +lat_0=46.8 +lon_0=2.337229104484 +x_0=600000 +y_0=2200000 +units=km"

# Importation of the boundary of France using sf
france <- st_read("data/boundaries/FRA_adm0.shp")
st_crs(france) <- 4326 # Set initial CRS to WGS84
france <- st_transform(france, crs = lambert2_proj) # Transform to Lambert 2

# Importation of the boundary of the departments
region <- st_read("data/boundaries/FRA_adm2.shp")
st_crs(region) <- 4326
region <- st_transform(region, crs = lambert2_proj)
# Drop Corsica (assuming row names 94 and 95 correspond to Corsica departments)
region <- region[!(rownames(region) %in% c("94", "95")), ]

# Load the boundary of the communes (saved RData workspace)
load("data/boundaries/FRA_adm5.RData")
communes <- gadm
# Importation of the population data
pop <- read.csv2("data/population/tab.pop.csv")
pop.dep <- aggregate(pop[, c("POP90")], list(dep = pop$C_Dpt), sum)
names(pop.dep) <- c("ID_2", "pop90")
# Merge population with regional spatial data
region <- merge(region, pop.dep, by = "ID_2", sort = FALSE)

##############################################
# I.2. Isolate Metropolitan France (excluding islands)

# Cast multipart polygons to single parts and extract the largest one (Metropolitan France)
france_parts <- st_cast(france, "POLYGON")
mainland_france <- france_parts[which.max(st_area(france_parts)), ]

##############################################
# I.3. Clean the data basis of GEIPAN

# Fix the missing year
cas[cas$num_cas == 852, "AAAA"] = "1993"
cas$AAAA <- as.numeric(as.character(cas$AAAA))

# Exclude large distance errors (centroid vs vertices)
cas <- cas[cas$max_distance < 20, ]

# Correct CRS assignment: 
cas_sf <- st_as_sf(cas, coords = c("Longitude", "Latitude"), crs = lambert2_proj)

# Spatial filter: keep only points falling inside mainland France
cas_fr_sf <- st_intersection(cas_sf, mainland_france)

# Create a new variable period using 'cut'
cas_fr_sf$period <- cut(cas_fr_sf$AAAA,
                        breaks = c(1935, 1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020),
                        labels = c("35-50", "50-60", "60-70", "70-80", "80-90", "90-00", "00-10", "10-20"),
                        include.lowest = TRUE)

# Extract coordinates for spatstat compatibility later in the script
coord.fr <- st_coordinates(cas_fr_sf)

# --- Compatibility Bridge (sf to sp) ---
lambert2 <- CRS(lambert2_proj)
coord.sp <- SpatialPoints(coord.fr, proj4string = lambert2)
region <- as(region, "Spatial")
p3 <- as(mainland_france, "Spatial")
cas.fr <- sf::st_drop_geometry(cas_fr_sf)

##############################################################
#####    II. Spatial exploratory data analysis
##############################################################

tab.cas <- sapply(over(region, coord.sp, returnList = TRUE), length)
region@data$nb_cas <- tab.cas

########## Figure 1
pdf("output/figures/carto1.pdf", width = 5, height = 3.5, pointsize = 8, bg = "white")
plotclr <- brewer.pal(8, "Blues")[c(1, 3, 6, 8)]
op <- par(mfrow = c(1, 2), mar = c(0, 0, 0, 0), mai = rep(0, 4), oma = rep(0, 4))

breaks <- classIntervals(region@data$nb_cas, n = 4, style = "jenks")$brks
interval <- round(breaks, digits = 0)
op2 <- par(mai = c(0, 0.9, 0, 1.5), mar = c(0, 1.7, 0, 2.5))
plot(region, col = plotclr[findInterval(region@data$nb_cas, interval, all.inside = TRUE)], border = 'grey')
decoup <- c(paste("<=", interval[2], sep = ""), paste("]", interval[2], ";", interval[3], "]", sep = ""), 
            paste("]", interval[3], ";", interval[4], "]", sep = ""), paste(">", interval[4], sep = ""))
op3 <- par(bg = "antiquewhite1")
legend("topleft", legend = decoup, title = "UAP Counts", fill = plotclr, cex = 0.8)
par(op3)
par(op2)

breaks <- classIntervals(region@data$nb_cas / region@data$pop90, n = 4, style = "jenks")$brks
interval <- round(breaks, digits = 9)
op2 <- par(mai = c(0, 1.5, 0, 0.5), mar = c(0, 2, 0, 1.7))
plot(region, col = plotclr[findInterval(region@data$nb_cas / region@data$pop90, interval, all.inside = TRUE)], border = 'grey')
interval <- round(interval * 100000, 2)
decoup <- c(paste("<=", interval[2], sep = ""), paste("]", interval[2], ";", interval[3], "]", sep = ""), 
            paste("]", interval[3], ";", interval[4], "]", sep = ""), paste(">", interval[4], sep = ""))
op3 <- par(bg = "antiquewhite1")
legend("topright", legend = decoup, title = "UAP Counts/100,000 inhabitants", fill = plotclr, cex = 0.8)
par(op3)
par(op2)
par(op)
dev.off()

##############################################################
# II.2. Analysis of the variables
##############################################################

PAN_class <- cas.fr$Class_finale
PAN_class <- ifelse(PAN_class == 0.125, "A",
                    ifelse(PAN_class == 0.375, "B",
                           ifelse(PAN_class == -1, "C", "D")))
cas.fr$PAN_class <- PAN_class

crois3 <- table(PAN_class, cas.fr$period)

##############################################################
#####    III. Spatial analysis of the PP
##############################################################

# Create the observation window directly from the sf object
W <- as.owin(mainland_france)

wp <- ppp(x = coord.fr[, 1], y = coord.fr[, 2], window = W)
marks(wp) <- factor(PAN_class)
split.wp <- split(wp)

pdf("output/figures/allPAN.pdf", width = 3.3, height = 3.3, pointsize = 8, bg = "white")
op = par(mar = c(0, 0, 0, 0), oma = rep(0, 4), mai = rep(0, 4))
plot(split.wp, main = "", pch = 3)
par(op)
dev.off()

cas.A <- split.wp[[1]]; cas.B <- split.wp[[2]]
cas.C <- split.wp[[3]]; cas.D <- split.wp[[4]]

# Jitter points to resolve exact duplicates for modeling
cache_casD <- "cache/casD_jittered.rds"
if (file.exists(cache_casD)) {
  cas.D <- readRDS(cache_casD)
} else {
  set.seed(12345)
  cas.D <- rjitter(cas.D, retry=TRUE, nsim=1, drop=TRUE)
  saveRDS(cas.D, cache_casD)
}

M <- quadrat.test(cas.D, nx = 6, ny = 6)

D.cas.D <- density(cas.D, sigma = 20, diggle = TRUE)
D.cas.A <- density(cas.A, sigma = 20, diggle = TRUE)

coords.D <- aggregate(cas.D$x, list(x = as.character(cas.D$x), y = as.character(cas.D$y)), length)
names(coords.D) <- c("x", "y", "weights")
coords.D$x = as.numeric(as.character(coords.D$x))
coords.D$y = as.numeric(as.character(coords.D$y))

pdf("output/figures/intensityD.pdf", width = 5, height = 3.5, pointsize = 6, bg = "white")
op <- par(mfrow = c(1, 2), mar = c(1, 1.5, 1, 1.5), mai = c(0, 0, 0, 0), oma = rep(0, 4))
op2 <- par(mai = c(0, 0.9, 0, 2.5), mar = c(0, 1.7, 0, 2.5))
plot(D.cas.D, main = "", axes = TRUE)
plot(p3, add = TRUE)
points(coords.D$x, coords.D$y, cex = sqrt(coords.D$weights), pch = 16, col = "red")
polygon(c(740, 860, 860, 740), c(1940, 1940, 2060, 2060), lty = 2, lwd = 1.7, border = "white")
legend("topleft", legend = c("1 PAN", "2 PAN", "3 PAN", "4 PAN"), pch = 16, pt.cex = sqrt(1:4))
par(op2)

op2 <- par(mai = c(0.9, 0.9, 0.9, 0.5), mar = c(0.5, 1.7, 0, 1))
plot(D.cas.D[owin(c(740, 860), c(1940, 2060))], 
     main = "", 
     axes = TRUE, 
     xlab = "longitude in kilometer", 
     ylab = "latitude in kilometer")
points(jitter(cas.D[owin(c(740, 860), c(1940, 2060))]$x), 
       jitter(cas.D[owin(c(740, 860), c(1940, 2060))]$y), 
       pch = 16, col = "red")
polygon(c(791.58, 799.07, 799.07, 791.58), c(1992.4, 1992.4, 2000, 2000))
text(795.325, 1996.2, "u*", pch = 15)
x1 <- seq((791.58 + 799.07) / 2 - 10, (791.58 + 799.07) / 2 + 10, 0.001)
lines(x1, 1996.2 + sqrt(100 - (x1 - 795.325)^2), lty = 2)
lines(x1, 1996.2 - sqrt(100 - (x1 - 795.325)^2), lty = 2)
segments(795.325, 1996.2, 800, 1996.2 - sqrt(100 - (795.325 - 800)^2))
text(800, 1996.2 - sqrt(100 - (795.325 - 800)^2), "k=0.000351 (10 km)", cex = 0.8, pos = 1)
x1 <- seq((791.58 + 799.07) / 2 - 20, (791.58 + 799.07) / 2 + 20, 0.001)
lines(x1, 1996.2 + sqrt(400 - (x1 - 795.325)^2), lty = 2)
lines(x1, 1996.2 - sqrt(400 - (x1 - 795.325)^2), lty = 2)
segments(795.325, 1996.2, 810, 1996.2 + sqrt(400 - (795.325 - 810)^2))
text(810, 1996.2 + sqrt(400 - (795.325 - 810)^2), "k=0.000241 (20 km)", cex = 0.8, pos = 4)
x1 <- seq((791.58 + 799.07) / 2 - 50, (791.58 + 799.07) / 2 + 50, 0.001)
lines(x1, 1996.2 + sqrt(2500 - (x1 - 795.325)^2), lty = 2)
lines(x1, 1996.2 - sqrt(2500 - (x1 - 795.325)^2), lty = 2)
segments(795.325, 1996.2, 760, 1996.2 + sqrt(2500 - (795.325 - 760)^2))
text(760, 1996.2 + sqrt(2500 - (795.325 - 760)^2), "k=1.7e-5 (50 km)", cex = 0.8, pos = 3)
par(op2)
par(op)
dev.off()

dmnorm(c(0, 20), varcov = matrix(c(400, 0, 0, 400), 2, 2))

##############################################################
#####    IV. Analysis of the covariates
##############################################################

# IV.1 The density of population
tab <- read.csv("data/population/tab.pop.csv", sep = ';', dec = ',')
tab <- tab[!is.na(tab$Code.INSEE), ]
pop.area <- read.csv2("data/population/pop1990.csv", skip = 0, nrows = 36680)
pop.area$CODEGEO <- as.numeric(as.character(pop.area$CODEGEO))
names(pop.area)[3] <- "pop90"

tab.pop <- merge(tab, pop.area, by.x = "Code.INSEE", by.y = "CODEGEO")

coord.france <- cbind(tab.pop$Longitude, tab.pop$Latitude)
sp.france.0 <- SpatialPoints(coord.france, CRS("+proj=longlat +ellps=WGS84"))
sp.france <- spTransform(sp.france.0, lambert2)

res.france <- point.in.polygon(coordinates(sp.france)[, 1], coordinates(sp.france)[, 2],
                               p3@polygons[[1]]@Polygons[[1]]@coords[, 1], p3@polygons[[1]]@Polygons[[1]]@coords[, 2])
ind.f.NULL <- which(res.france == 0)

wp.france <- ppp(x = coordinates(sp.france)[-ind.f.NULL, 1], y = coordinates(sp.france)[-ind.f.NULL, 2],
                 window = W, unitname = c("metre", "metres"))
marks(wp.france) <- tab.pop$pop90[-ind.f.NULL]

D.pop <- Smooth(wp.france, sigma = 5, dimyx = c(128, 128))
cdf.test(cas.D, D.pop)

fenetre <- owin(c(785, 810), c(1985, 2010))

pdf("output/figures/pop1.pdf", width = 5, height = 3.5, pointsize = 6, bg = "white")
op <- par(mfrow = c(1, 2), mar = c(0, 4, 2, 2), oma = rep(0, 4))
op2 <- par(mai = c(1.1, 0.6, 1.1, 0.6))
plot(fenetre$xrange, fenetre$yrange, type = "n", axes = TRUE, main = "", xlab = "longitude in kilometer", ylab = "latitude in kilometer", asp = 1)
points(wp.france[fenetre]$x, wp.france[fenetre]$y, pch = 16, col = 'royalblue', cex = (marks(wp.france[fenetre]) / 1729)^0.57 * 5)
points(wp.france[fenetre]$x, wp.france[fenetre]$y, pch = 1, cex = (marks(wp.france[fenetre]) / 1729)^0.57 * 5)
plot(communes[communes@data$ID_2 == 26 | communes@data$ID_2 == 7, ], add = TRUE)
par(op2)

op2 <- par(mai = c(0.9, 0.5, 0.9, 0.5))
plot(D.pop[fenetre], main = "", axes = TRUE)
points(wp.france[fenetre]$x, wp.france[fenetre]$y, pch = 16, col = 'red', cex = (marks(wp.france[fenetre]) / 1729)^0.57 * 5)
points(wp.france[fenetre]$x, wp.france[fenetre]$y, pch = 1, cex = (marks(wp.france[fenetre]) / 1729)^0.57 * 5)

polygon(c(791.58, 799.07, 799.07, 791.58), c(1992.4, 1992.4, 2000, 2000), border = "white")
text(795.325, 1996.2, "u*", pch = 15, col = "black", cex = 1.5)

x1 <- seq((791.58 + 799.07) / 2 - 2, (791.58 + 799.07) / 2 + 2, 0.001)
lines(x1, 1996.2 + sqrt(4 - (x1 - 795.325)^2), lwd = 1.5, lty = 2, col = "white")
lines(x1, 1996.2 - sqrt(4 - (x1 - 795.325)^2), lwd = 1.5, lty = 2, col = "white")
segments(795.325, 1996.2, 797, 1996.2 - sqrt(4 - (795.325 - 797)^2), lwd = 1.5, col = "white")
text(797, 1996.2 - sqrt(4 - (795.325 - 797)^2), "k=0.00588 (2 km)", cex = 0.9, pos = 1, col = "black")

x1 <- seq((791.58 + 799.07) / 2 - 5, (791.58 + 799.07) / 2 + 5, 0.001)
lines(x1, 1996.2 + sqrt(25 - (x1 - 795.325)^2), lty = 2, lwd = 1.5, col = "white")
lines(x1, 1996.2 - sqrt(25 - (x1 - 795.325)^2), lty = 2, lwd = 1.5, col = "white")
segments(795.325, 1996.2, 798, 1996.2 + sqrt(25 - (795.325 - 798)^2), lwd = 1.5, col = "white")
text(798, 1996.2 + sqrt(25 - (795.325 - 798)^2), "k=0.00386 (5 km)", cex = 0.9, pos = 4)

x1 <- seq((791.58 + 799.07) / 2 - 10, (791.58 + 799.07) / 2 + 10, 0.001)
lines(x1, 1996.2 + sqrt(100 - (x1 - 795.325)^2), lwd = 1.5, lty = 2, col = "white")
lines(x1, 1996.2 - sqrt(100 - (x1 - 795.325)^2), lwd = 1.5, lty = 2, col = "white")
segments(795.325, 1996.2, 800, 1996.2 - sqrt(100 - (795.325 - 800)^2), lwd = 1.5, col = "white")
text(800, 1996.2 - sqrt(100 - (795.325 - 800)^2), "K=0.000862 (10 km)", cex = 0.9, pos = 3)
legend("topleft", title = "Population density", legend = c("10", "250", "750"), pch = 16, pt.cex = (c(10, 250, 750) / 1757)^0.57 * 5, cex = 0.7)
par(op2)
par(op)

dev.off()

interval <- c(0, 50, 100, 500, 2000, 12000)
tab.stat <- as.table(round(rbind(prop.table(table(findInterval(na.omit(as.vector(D.pop$v)), interval, all.inside = TRUE))),
                                 prop.table(table(findInterval(D.pop[cas.D], interval, all.inside = TRUE)))), 3))
row.names(tab.stat) <- c("Pixels", "PAN D")
colnames(tab.stat) <- c(paste("<=", interval[2], sep = ""), paste("]", interval[2], ";", interval[3], "]", sep = ""), paste("]", interval[3], ";", interval[4], "]", sep = ""),
                        paste("]", interval[4], ";", interval[5], "]", sep = ""), paste(">", interval[5], sep = ""))

print(xtable(tab.stat, digits = 3), file="output/tables/pop_density_table.tex")

D.pop$v[which(D.pop$v < 0)] = 1

tc <- colourmap(terrain.colors(12)[c(1, 3, 5, 9, 12)], breaks = c(0, 50, 100, 500, 2000, 12000))
plot(D.pop, col = tc, main = "Variable 1: Population Density", axes = TRUE)
plot(p3, add = TRUE)
points(coords.D$x, coords.D$y, cex = sqrt(coords.D$weights), pch = 16, col = 'darkblue')

##############################################################
# IV.2 The nuclear sites
##############################################################

# 1. ALWAYS execute the preparation of spatial points (fast)
nucleaire <- read.csv2("data/environment/nuclear.csv")
nucleaire.sp <- SpatialPoints(cbind(nucleaire$Longitude, nucleaire$Latitude), CRS("+proj=longlat +ellps=WGS84"))
nucleaire.sp.m <- spTransform(nucleaire.sp, lambert2) # This object is now found!

coord.nucleaire <- as.data.frame(nucleaire.sp.m)
colnames(coord.nucleaire) <- c("x", "y")

# 2. ONLY cache the heavy computation
cache_nucl <- "cache/D_nucleaire_obj.rds"

if (file.exists(cache_nucl)) {
  D.nucleaire <- readRDS(cache_nucl)
} else {
  grille <- as.data.frame(D.pop)[, 1:2]
  nb.nucleaire <- NULL
  
  for (i in 1:nrow(grille)) {
    x.sp <- SpatialPoints(rbind(grille[i, ], coord.nucleaire))  
    x.dnn <- dnearneigh(x.sp, 0, 20)  
    nb.nucleaire <- c(nb.nucleaire, ifelse(length(x.dnn[[1]]) == 1 && x.dnn[[1]] == 0, 0, length(x.dnn[[1]]))) 
  }

  nucleaire.v <- matrix(NA, 128, 128)
  nucleaire.v[!is.na(D.pop$v)] <- nb.nucleaire
  D.nucleaire <- im(nucleaire.v, D.pop$xcol, D.pop$yrow)
  saveRDS(D.nucleaire, cache_nucl)
}

interval <- c(0, 1, 5, 10, 15)
tab.stat <- as.table(round(rbind(prop.table(table(findInterval(na.omit(as.vector(D.nucleaire$v)), interval, all.inside = TRUE))),
                                prop.table(table(findInterval(D.nucleaire[cas.D], interval, all.inside = TRUE)))), 3))
row.names(tab.stat) <- c("Pixels", "PAN D")
colnames(tab.stat) <- c("0", "1-5", "6-10", "11-14")
print(xtable(tab.stat, digits = 3), file="output/tables/nuclear_sites_table.tex")

res.nucleaire <- point.in.polygon(coordinates(nucleaire.sp.m)[, 1], coordinates(nucleaire.sp.m)[, 2],
                                  p3@polygons[[1]]@Polygons[[1]]@coords[, 1], p3@polygons[[1]]@Polygons[[1]]@coords[, 2])
ind.f.NULL <- which(res.nucleaire == 0)

wp.nucleaire <- ppp(x = coordinates(nucleaire.sp.m)[-ind.f.NULL, 1], y = coordinates(nucleaire.sp.m)[-ind.f.NULL, 2],
                    window = W, unitname = c("metre", "metres"))

D.nucleaire.2 <- density(wp.nucleaire, sigma = 20, diggle = TRUE)
D.nucleaire.2$v[which(D.nucleaire.2$v < 0)] = 0

# Representation of figure 5 in the article
log.nucl <- read.pnm("data/environment/radioactivite.pbm")

pdf("output/figures/nuclear.pdf", width = 5, height = 3.5, pointsize = 7, bg = "white")
op <- par(mfrow = c(1, 2), mar = c(0, 4, 2, 2), oma = rep(0, 4))

plot(p3)

# Extract coordinates beforehand to avoid subsetting errors in the loop
coords_nucl <- coordinates(nucleaire.sp.m)

# Safely loop through the matrix rows
for (i in 1:nrow(coords_nucl)) {
  addlogo(log.nucl, 
          px = c(coords_nucl[i, 1], coords_nucl[i, 1] + 30),
          py = c(coords_nucl[i, 2], coords_nucl[i, 2] + 30), 
          asp = 1)
}

polygon(c(735, 865, 865, 735), c(1935, 1935, 2065, 2065), border = "black", lty = 2)

op2 <- par(mai = c(1, 0.5, 1, 0.5))
plot(D.nucleaire.2[owin(c(735, 865), c(1935, 2065))], main = "", axes = TRUE, xlab = "longitude in kilometer", ylab = "latitude in kilometer")
polygon(c(791.58, 799.07, 799.07, 791.58), c(1992.4, 1992.4, 2000, 2000), border = "white", lwd = 2)
# text(798.325, 1996.2, "u*", pch = 15, add = TRUE, col = "white", adj = 2)
text(798.325, 1996.2, "u*", col="white", adj=2)

ppp.w <- ppp(coord.nucleaire$x, coord.nucleaire$y, window = owin(c(725, 875), c(1925, 2075)))
for (i in 1:length(ppp.w$x)) {
  # eps <- jitter(0, 400)
  # Generate a random spatial shift between -4 and 4 km (since units are km)
  eps <- runif(1, -4, 4)
  addlogo(log.nucl, px = c(ppp.w$x[i], ppp.w$x[i] + 10) + eps,
          py = c(ppp.w$y[i], ppp.w$y[i] + 10), asp = 1)
}
par(op2)
par(op)
dev.off()

cdf.test(cas.D, D.nucleaire.2)

##############################################################
# IV.3 Contaminated Land
##############################################################

# 1. Importation et préparation
pollues <- read.csv2("data/environment/pollues.csv", sep = ";", dec = ",", fileEncoding = "UTF-8")
names(pollues) <- c("Code", "communes", "part_territoires", "sites.pollues_2012")

pollues$Code <- trimws(as.character(pollues$Code))
tab$Code.INSEE <- trimws(as.character(tab$Code.INSEE))

tab.pollues <- merge(tab, pollues, by.x = "Code.INSEE", by.y = "Code", sort = TRUE)

# 2. OPTIMISATION : Filtrage des communes sans pollution (Bonus rapidité)
# On ne garde que les communes ayant au moins 1 site pollué et des coordonnées valides
keep_pollues <- tab.pollues$sites.pollues_2012 > 0 & 
                !is.na(tab.pollues$Longitude) & 
                !is.na(tab.pollues$Latitude)

tab.filtered <- tab.pollues[keep_pollues, ]

# 3. Calcul de la densité (grille 128x128)
df_pop <- as.data.frame(D.pop)
valid_pixels <- !is.na(df_pop$value)
grille_valid <- df_pop[valid_pixels, 1:2]
nb.pollues_valid <- numeric(nrow(grille_valid))

# Préparation des coordonnées des sites pollués pour la boucle
x.sp.sites <- SpatialPoints(cbind(tab.filtered$Longitude, tab.filtered$Latitude), 
                            proj4string = CRS("+proj=longlat +ellps=WGS84"))
x.sp.sites <- spTransform(x.sp.sites, lambert2)
coord.pollues.km <- coordinates(x.sp.sites)

# Boucle optimisée sur les pixels terrestres uniquement
for (i in 1:nrow(grille_valid)) {
  dist_sq <- (grille_valid[i, 1] - coord.pollues.km[, 1])^2 + 
             (grille_valid[i, 2] - coord.pollues.km[, 2])^2
  indices <- which(dist_sq < 25) # Rayon 5km
  nb.pollues_valid[i] <- ifelse(length(indices) == 0, 0, sum(tab.filtered$sites.pollues_2012[indices]))
}

pollues.v <- matrix(NA, 128, 128)
# pollues.v[valid_pixels] <- nb.pollues_valid
idx <- which(valid_pixels)
pollues.v[idx] <- nb.pollues_valid
D.pollues <- im(pollues.v, D.pop$xcol, D.pop$yrow)

# 4. Création de l'objet ppp pour D.pollues2 (Correction crash "marks")
# On transforme tous les points potentiels
sp.pollues.all <- SpatialPoints(cbind(tab.pollues$Longitude, tab.pollues$Latitude), 
                                CRS("+proj=longlat +ellps=WGS84"))
sp.pollues.km <- spTransform(sp.pollues.all, lambert2)
xy.all <- coordinates(sp.pollues.km)

# On filtre ceux qui sont à l'intérieur de la fenêtre W (France métropolitaine)
# ET on filtre les marques (sites.pollues_2012) en même temps
insideW <- inside.owin(xy.all[,1], xy.all[,2], W)

wp.pollues <- ppp(x = xy.all[insideW, 1], 
                  y = xy.all[insideW, 2], 
                  window = W, 
                  marks = tab.pollues$sites.pollues_2012[insideW],
                  unitname = c("kilometre", "kilometres"))

# Calcul de l'intensité lissée (D.pollues2)
D.pollues2 <- density(wp.pollues, sigma = 5, weights = marks(wp.pollues))

# 5. Visualisation
pdf("output/figures/pollues.pdf", width = 5, height = 3.5, pointsize = 8, bg = "white")
op <- par(mfrow = c(1, 2), mar = c(0, 4, 2, 2), oma = rep(0, 4))
plot(W, main = "")

# On n'affiche que les points avec pollution > 0 pour la clarté
wp.visu <- subset(wp.pollues, marks > 0)
points(wp.visu$x, wp.visu$y, cex = (marks(wp.visu) / 36)^0.57 * 3, pch = 16, col = "orange")
points(wp.visu$x, wp.visu$y, cex = (marks(wp.visu) / 36)^0.57 * 3, pch = 1)

polygon(c(735, 865, 865, 735), c(1935, 1935, 2065, 2065), border = "red", lty = 2)
legend("topleft", legend = c("1", "5", "15"), pch = 16, pt.cex = (c(1, 5, 15) / 36)^0.57 * 3, col = "orange", cex = 0.8)

op2 <- par(mai = c(1, 0.5, 1, 0.5))
plot(D.pollues2[owin(c(735, 865), c(1935, 2065))], main = "", axes = TRUE, 
     xlab = "longitude (km)", ylab = "latitude (km)")

# Zoom sur la zone d'intérêt
ppp.zoom <- subset(wp.pollues, window = owin(c(735, 865), c(1935, 2065)))
points(ppp.zoom$x, ppp.zoom$y, pch = 16, cex = (marks(ppp.zoom)/36)^0.57*3, col = 'red')

polygon(c(791.58, 799.07, 799.07, 791.58), c(1992.4, 1992.4, 2000, 2000), border = "white", lwd = 2)
text(805.325, 1996.2, "u*", pch = 15, col = "white", cex = 1.2, adj = 2)
par(op2)
par(op)
dev.off()

# Test de Kolmogorov-Smirnov
cdf.test(cas.D, D.pollues2)

##############################################################
# IV.4 Annual sunshine time (FIXED CRS & MARKS)

if (file.exists("cache/D_soleil.rds")) {
  D.soleil <- readRDS("cache/D_soleil.rds")
} else {
  gr <- read.asciigrid("data/sun/pvgis_g13year00.asc")
  
  # FIX: We must set the CRS on the SpatialGrid object FIRST
  proj4string(gr) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

  # Create the SpatialPointsDataFrame
  gr.sp <- SpatialPointsDataFrame(SpatialPoints(gr), gr@data)
  
  # Ensure the SpatialPoints object also has the CRS defined
  proj4string(gr.sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
  
  # Now the transform will work perfectly
  gr2.sp <- spTransform(gr.sp, lambert2)
  coord.soleil <- coordinates(gr2.sp)

  # Filtering bbox
  keep_bbox <- (coord.soleil[,1] > 72.4 & coord.soleil[,1] < 1031.3) &
               (coord.soleil[,2] > 1703.2 & coord.soleil[,2] < 2677.3)
  gr3.sp <- gr2.sp[keep_bbox, ]

  xy <- coordinates(gr3.sp)
  mark.soleil <- as.numeric(gr3.sp@data[,1])

  # Replace NA with median to avoid rejected points
  mark.soleil[is.na(mark.soleil)] <- median(mark.soleil, na.rm = TRUE)

  # Filter inside W to prevent the "number of points != number of marks" error
  insideW <- inside.owin(xy[,1], xy[,2], W)

  wp.soleil <- ppp(x = xy[insideW,1],
                   y = xy[insideW,2],
                   window = W,
                   marks = mark.soleil[insideW])

  D.soleil <- Smooth(wp.soleil, sigma = 5)
  saveRDS(D.soleil, "cache/D_soleil.rds")
}

# Documentation check: cdf.test will show the correlation with UAP
cdf.test(cas.D, D.soleil)

##############################################################
# IV.5 The percentage of wetlands and water bodies (ROBUST VERSION)

if (file.exists("cache/D_environ.rds")) {
  env_list <- readRDS("cache/D_environ.rds")
  D.forets <- env_list$forets
  D.humidite <- env_list$humidite
} else {
  # We read WITHOUT checking names to avoid the multibyte error
  # And we rename immediately with clean ASCII names
  forets_raw <- read.csv2("data/environment/forets.csv", sep = ';', dec = ',', 
                          fileEncoding = "UTF-8", check.names = FALSE)
  # Based on your snippet: 1=Code, 2=Communes, 3=Artif, 4=Forets
  names(forets_raw) <- c("Code", "Nom", "part_artif", "part_foret")
  
  humidite_raw <- read.csv2("data/environment/humides.csv", sep = ';', dec = ',', 
                            fileEncoding = "UTF-8", check.names = FALSE)
  # Based on your snippet: 1=Code, 2=Communes, 3=Humide
  names(humidite_raw) <- c("Code", "Nom", "part_humide")

  # Alignment: ensure Code is a clean 5-digit string (ex: "01001")
  forets_raw$Code <- trimws(as.character(forets_raw$Code))
  humidite_raw$Code <- trimws(as.character(humidite_raw$Code))
  tab$Code.INSEE <- trimws(as.character(tab$Code.INSEE))

  tab.forets <- merge(tab, forets_raw, by.x = "Code.INSEE", by.y = "Code")
  tab.humidite <- merge(tab, humidite_raw, by.x = "Code.INSEE", by.y = "Code")

  # Security: check if merge succeeded
  if(nrow(tab.forets) == 0) stop("Merge failed for forests. Check INSEE codes.")

  # Spatial Points transformation
  sp_env <- SpatialPoints(cbind(tab.forets$Longitude, tab.forets$Latitude), 
                          proj4string = CRS("+proj=longlat +ellps=WGS84"))
  xy.env <- coordinates(spTransform(sp_env, lambert2))
  
  # Filter inside W to avoid the marks mismatch bug
  insideW_env <- inside.owin(xy.env[,1], xy.env[,2], W)
  
  wp.forets <- ppp(xy.env[insideW_env, 1], xy.env[insideW_env, 2], window = W, 
                   marks = tab.forets$part_foret[insideW_env])
  wp.humidite <- ppp(xy.env[insideW_env, 1], xy.env[insideW_env, 2], window = W, 
                     marks = tab.humidite$part_humide[insideW_env])

  D.forets <- Smooth(wp.forets, sigma = 5, dimyx = c(128, 128))
  D.humidite <- Smooth(wp.humidite, sigma = 5, dimyx = c(128, 128))
  D.humidite$v[D.humidite$v < 0] <- 0
  
  saveRDS(list(forets = D.forets, humidite = D.humidite), "cache/D_environ.rds")
}

##############################################################
# IV.7 Variable airport

if (file.exists("cache/D_aero.rds")) {
  D.aero.2 <- readRDS("cache/D_aero.rds")
} else {
  aeroport <- read.csv2("data/environment/aeroport_clc.csv", fileEncoding = "UTF-8")
  names(aeroport) <- c("Code.INSEE", "Nom_com", "aeroport")
  aeroport <- merge(aeroport, tab.pop, by = "Code.INSEE")

  aero.xy <- coordinates(spTransform(SpatialPoints(cbind(aeroport$Longitude, aeroport$Latitude), 
                                                   proj4string = CRS("+proj=longlat +ellps=WGS84")), lambert2))
  
  insideW_aero <- inside.owin(aero.xy[,1], aero.xy[,2], W)
  wp.aero <- ppp(aero.xy[insideW_aero, 1], aero.xy[insideW_aero, 2], window = W, 
                 marks = aeroport$aeroport[insideW_aero])

  D.aero.2 <- density(wp.aero, sigma = 5, weights = marks(wp.aero), diggle = TRUE)
  saveRDS(D.aero.2, "cache/D_aero.rds")
}


##############################################################
# IV.8 Political Context: FN Vote 2012 (HARD-CLIP + SIGMA SWEEP SAFE)

cache_vote_base <- "cache/voteFN_base_im.rds"
cache_vote      <- "cache/D_voteFN_hardclip.rds"
cache_sigma     <- "cache/voteFN_sigma_sweep.rds"

#  Exact window
W_exact <- spatstat.geom::as.owin(mainland_france)
M_model <- spatstat.geom::as.mask(W_exact, dimyx = dim(D.pop$v))
M_plot  <- spatstat.geom::as.mask(W_exact, dimyx = c(512, 512))

# Raster -> im 
to_im <- function(r) {
  v <- as.matrix(r)
  v <- v[nrow(v):1, ]  # flip Y
  ex <- raster::extent(r)
  x <- seq(ex@xmin, ex@xmax, length.out = ncol(v))
  y <- seq(ex@ymin, ex@ymax, length.out = nrow(v))
  spatstat.geom::im(v, xcol = x, yrow = y)
}

# 1) change base num_im/den_im (on grid, clip BEFORE blur)
if (file.exists(cache_vote_base)) {

  base <- readRDS(cache_vote_base)
  num_im <- base$num_im
  den_im <- base$den_im

} else {

  fn <- png::readPNG("data/environment/vote_FN_2012.png")
  I0 <- 1 - (0.2126*fn[,,1] + 0.7152*fn[,,2] + 0.0722*fn[,,3])

  A <- if (dim(fn)[3] >= 4) fn[,,4] else matrix(1, nrow(I0), ncol(I0))
  A[A < 0.01] <- 0

  r_num <- raster::raster(I0 * A)
  r_den <- raster::raster(A)

  raster::extent(r_num) <- raster::extent(r_den) <- c(-4.9940925, 9.8557315, 40.955078, 51.909617)
  raster::crs(r_num) <- raster::crs(r_den) <- "+proj=longlat +datum=WGS84"

  r_num_l2 <- raster::projectRaster(r_num, crs = lambert2_proj, method = "ngb")
  r_den_l2 <- raster::projectRaster(r_den, crs = lambert2_proj, method = "ngb")

  num_im0 <- to_im(r_num_l2)
  den_im0 <- to_im(r_den_l2)

  # Clip AVANT flou (sur grille modèle)
  num_im <- spatstat.geom::as.im(num_im0, W = M_model)
  den_im <- spatstat.geom::as.im(den_im0, W = M_model)

  saveRDS(list(num_im=num_im, den_im=den_im), cache_vote_base)
}

# 2) Construire/charger D.vote_fn (sigma par défaut = 8 km) + plot hardclip
if (file.exists(cache_vote)) {

  tmp <- readRDS(cache_vote)
  D.vote_fn      <- tmp$D.vote_fn
  D.vote_fn_plot <- tmp$D.vote_fn_plot

} else {

  sigma_km <- 8
  num_bl <- spatstat.explore::blur(num_im, sigma = sigma_km, edge = TRUE, bleed = TRUE)
  den_bl <- spatstat.explore::blur(den_im, sigma = sigma_km, edge = TRUE, bleed = TRUE)

  D.vote_fn_full <- spatstat.geom::eval.im(num_bl / den_bl)

  # Plot : hard clip
  D.vote_fn_plot <- spatstat.geom::as.im(D.vote_fn_full, W = M_plot)
  inside_mat <- outer(D.vote_fn_plot$yrow, D.vote_fn_plot$xcol,
                      function(y, x) spatstat.geom::inside.owin(x, y, W_exact))
  D.vote_fn_plot$v[!inside_mat] <- NA

  den_bl_plot <- spatstat.geom::as.im(den_bl, W = M_plot)
  D.vote_fn_plot$v[den_bl_plot$v < 1e-6] <- NA

  # Model :  grid D.pop, 0 outside France
  D.vote_fn <- spatstat.geom::as.im(D.vote_fn_full, W = M_model)
  D.vote_fn$v[is.na(D.vote_fn$v)] <- 0

  saveRDS(list(D.vote_fn=D.vote_fn, D.vote_fn_plot=D.vote_fn_plot), cache_vote)
}

# Check 
plot(D.vote_fn_plot, main="Vote FN 2012 (Hard-Clipped & Mainland)")
plot(p3, add = TRUE, border="black", lwd=0.7)

###################################################
# SIGMA SWEEP SAFE (utilise num_im/den_im qui existent toujours)

fit_for_sigma <- function(sigma_km, num_im, den_im, M_model) {
  num_bl <- spatstat.explore::blur(num_im, sigma = sigma_km, edge = TRUE, bleed = TRUE)
  den_bl <- spatstat.explore::blur(den_im, sigma = sigma_km, edge = TRUE, bleed = TRUE)

  D_full <- spatstat.geom::eval.im(num_bl / den_bl)
  D_sig  <- spatstat.geom::as.im(D_full, W = M_model)
  D_sig$v[is.na(D_sig$v)] <- 0

  fit <- ppm(cas.D, ~ voteFN, covariates = list(voteFN = D_sig))
  s   <- summary(fit)

  z <- s$coefs.SE.CI["voteFN", "Zval"]
  data.frame(
    sigma_km = sigma_km,
    coef     = coef(fit)[["voteFN"]],
    se       = s$coefs.SE.CI["voteFN", "S.E."],
    z        = z,
    p        = 2 * pnorm(-abs(z)),
    logLik   = as.numeric(logLik(fit)),
    AIC      = AIC(fit)
  )
}

sigmas <- c(4, 6, 8, 10, 12, 15, 18, 20)

if (file.exists(cache_sigma)) {
  sweep <- readRDS(cache_sigma)
} else {
  sweep <- do.call(rbind, lapply(sigmas, fit_for_sigma, num_im=num_im, den_im=den_im, M_model=M_model))
  saveRDS(sweep, cache_sigma)
}

print(sweep[order(sweep$AIC), ])

###################################################
################# All the covariates

data.var <- data.frame(
  pop = log(as.vector(D.pop$v)), 
  pollues = as.vector(D.pollues2$v), 
  nucleaire = as.vector(D.nucleaire$v),
  humidite = as.vector(D.humidite$v), 
  soleil = log(as.vector(D.soleil$v)), 
  aero = as.vector(D.aero.2$v), 
  forets = as.vector(D.forets$v),
  voteFN = as.vector(D.vote_fn$v) # <--- AJOUTÉ
)
cor(na.omit(data.var))

print(xtable(rbind(
  summary(as.vector(D.pop$v))[1:6],
  summary(as.vector(D.pollues2$v))[1:6],
  summary(as.vector(D.nucleaire$v))[1:6],
  summary(as.vector(D.humidite$v))[1:6],
  summary(as.vector(D.soleil$v))[1:6],
  summary(as.vector(D.aero.2$v))[1:6],
  summary(as.vector(D.forets$v))[1:6],
  summary(as.vector(D.vote_fn$v))[1:6] # <--- AJOUTÉ
)), file="output/tables/covariates_summary.tex")

pdf("output/figures/covariates1.pdf", width = 4, height = 4, pointsize = 6, bg = "white")
op = par(mfrow = c(2, 2))
tc <- colourmap(terrain.colors(12)[c(1, 3, 5, 9, 12)], breaks = c(0, 50, 100, 500, 2000, 12000))
plot(D.pop, col = tc, main = "Variable 1: Density of population")
plot(p3, add = TRUE)
points(coords.D$x, coords.D$y, cex = sqrt(coords.D$weights) * 0.6, pch = 3, col = 'darkblue')

tc <- colourmap(terrain.colors(12)[c(1, 3, 5, 9, 12)], breaks = c(0, 0.0002, 0.0005, 0.002, 0.007, 0.009))
plot(D.nucleaire.2, col = tc, main = "Variable 2: Nuclear risk")
plot(p3, add = TRUE)
points(coords.D$x, coords.D$y, cex = sqrt(coords.D$weights) * 0.6, pch = 3, col = 'darkblue')

tc <- colourmap(terrain.colors(12)[c(1, 3, 5, 9, 12)], breaks = c(0, 0.001, 0.02, 0.08, 0.5, 1))
plot(D.pollues2, col = tc, main = "Variable 3: Contaminated sites")
plot(p3, add = TRUE)
points(coords.D$x, coords.D$y, cex = sqrt(coords.D$weights) * 0.6, pch = 3, col = 'darkblue')

tc <- colourmap(brewer.pal(9, "Blues")[c(1, 3, 5, 7, 9)], breaks = c(0, 1, 5, 10, 15, 100))
plot(D.humidite, col = tc, main = "Variable 4: Rate of moisture")
plot(p3, add = TRUE)
points(coords.D$x, coords.D$y, cex = sqrt(coords.D$weights) * 0.6, pch = 3, col = 'red')
par(op)
dev.off()

pdf("output/figures/covariates2.pdf", width = 4, height = 4, pointsize = 6, bg = "white")
op = par(mfrow = c(2, 2))
plot(D.soleil, main = "Variable 5: Yearly sum of global irradiation")
plot(p3, add = TRUE)
points(coords.D$x, coords.D$y, cex = sqrt(coords.D$weights) * 0.6, pch = 3, col = 'red')

plot(D.aero.2, main = "Variable 6: Airport installation")
plot(p3, add = TRUE)
points(coords.D$x, coords.D$y, cex = sqrt(coords.D$weights) * 0.6, pch = 3, col = 'red')

tc <- colourmap(rev(terrain.colors(10)[c(1, 3, 5, 7, 9)]), breaks = c(0, 20, 40, 60, 80, 100))
plot(D.forets, col = tc, main = "Variable 7: Rate of land by forest")
plot(p3, add = TRUE)
points(coords.D$x, coords.D$y, cex = sqrt(coords.D$weights) * 0.6, pch = 3, col = 'darkblue')

plot(D.cas.A, main = "Variable 8: UAP A")
plot(p3, add = TRUE)
points(coords.D$x, coords.D$y, cex = sqrt(coords.D$weights) * 0.6, pch = 3, col = 'red')

plot(D.vote_fn, main = "Variable 9: Vote FN 2012 Intensity")
plot(p3, add = TRUE)
points(coords.D$x, coords.D$y, cex = sqrt(coords.D$weights) * 0.6, pch = 3, col = 'red')

pdf("output/figures/model_breakdown.pdf", width = 8, height = 4)
op <- par(mfrow = c(1, 3))

# Effect pop
plot(rhohat(cas.D, D.pop), main="Effect of Population")

# Effect Nuclear
plot(rhohat(cas.D, D.nucleaire.2), main="Effect of Nuclear")

# Effect FN Vote FN
plot(rhohat(cas.D, D.vote_fn), main="Effect of FN Vote")

# Brute dependancy
plot(rhohat(cas.D, D.vote_fn), main="Relation brute : Intensité PAN D vs Vote FN")

par(op)
dev.off()

fit.D_test <- ppm(cas.D, ~ voteFN, covariates = list(voteFN = D.vote_fn))
print(summary(fit.D_test))


##############################################################
#####    V. Modelisation
##############################################################

# --- 1. (Kolmogorov-Smirnov) ---
cdf.test(cas.D, D.humidite)
cdf.test(cas.D, D.aero.2)
cdf.test(cas.D, D.forets)
cdf.test(cas.D, D.soleil)
cdf.test(cas.D, D.cas.A)
cdf.test(cas.D, D.vote_fn) 

# --- 2. (Full Model) ---
fit.D_full <- ppm(cas.D, ~ log(pop) + pollues + nucleaire + humidite + log(soleil) + aero + forets + cas.A + voteFN,
             covariates = list(pop = D.pop, aero = D.aero.2, nucleaire = D.nucleaire.2, 
                               soleil = D.soleil, forets = D.forets, humidite = D.humidite, 
                               pollues = D.pollues2, cas.A = D.cas.A, voteFN = D.vote_fn))
res = summary(fit.D_full)
print(res)

# Export Full model
t.stat = coefficients(fit.D_full) / res$coefs.SE.CI$S.E.
print(xtable(data.frame(coefficients(fit.D_full), res$coefs.SE.CI$S.E., 2 * (pt(-abs(t.stat), 9000))), 
             digits = 5), file="output/tables/model_full.tex")

# --- 3. Final Model (select variables) ---
fit.D <- ppm(cas.D, ~ log(pop) + pollues + nucleaire + voteFN,
             covariates = list(pop = D.pop, aero = D.aero.2, nucleaire = D.nucleaire.2, 
                               soleil = D.soleil, forets = D.forets, humidite = D.humidite, 
                               pollues = D.pollues2, cas.A = D.cas.A, voteFN = D.vote_fn))
res = summary(fit.D)
print(res)

t.stat = coefficients(fit.D) / res$coefs.SE.CI$S.E.
print(xtable(data.frame(coefficients(fit.D), res$coefs.SE.CI$S.E., 2 * (pt(-abs(t.stat), 9000))), 
             digits = 30), file="output/tables/model_final.tex")

# --- 4. Diagnostics and residuals ---
pred.ppm <- predict(fit.D)
plot(pred.ppm)

res.ppm <- eval.im((D.cas.D - pred.ppm))
tc <- colourmap(cm.colors(80)[c(1, 32:49, 80)], breaks = quantile(as.vector(res.ppm$v), probs = seq(0, 1, 0.05), na.rm = TRUE))

pdf("output/figures/residuals.pdf", width = 3.3, height = 3.3, pointsize = 8, bg = "white")
op = par(mar = c(3, 3, 0, 1))
plot(res.ppm, col = tc, main = "")
plot(p3, add = TRUE)
points(coords.D$x, coords.D$y, cex = sqrt(coords.D$weights), pch = 3, col = 'darkblue')
par(op)
dev.off()

# local predictions
prod(exp((c(1, log(D.pop[owin(c(785, 805), c(1985, 2005))]$v[2, 2]), 
            D.pollues2[owin(c(785, 805), c(1985, 2005))]$v[2, 2], 
            D.nucleaire.2[owin(c(785, 805), c(1985, 2005))]$v[2, 2],
            D.vote_fn[owin(c(785, 805), c(1985, 2005))]$v[2, 2]) * coefficients(fit.D))))

predict(fit.D)[owin(c(785, 805), c(1985, 2005))]$v[2, 2]
D.cas.D[owin(c(785, 805), c(1985, 2005))]$v[2, 2] - predict(fit.D)[owin(c(785, 805), c(1985, 2005))]$v[2, 2]

cat("mean signed =", mean(res.ppm$v, na.rm=TRUE), "\n")
cat("mean abs    =", mean(abs(res.ppm$v), na.rm=TRUE), "\n")
cat("RMS         =", sqrt(mean(res.ppm$v^2, na.rm=TRUE)), "\n")
cat("max abs     =", max(abs(res.ppm$v), na.rm=TRUE), "\n")

anova(update(fit.D_full, . ~ . - voteFN), fit.D_full, test="Chisq")