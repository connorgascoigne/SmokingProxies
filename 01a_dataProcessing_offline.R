# 0. set ----
## 0.1 packages ----

library(tidyverse)
library(sf)
library(patchwork)

## 0.2 directories ---- 

# extract file location of this script
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

# retrieve directories
dir.home <- paste(code.path.splitted[1: (length(code.path.splitted)-2)], collapse = "/")
dir.data <- paste0(dir.home, '/data')
dir.res <- paste0(dir.home, '/results')
dir.data.spatial <- paste0(dir.data, '/spatial')
dir.data.organised <- paste0(dir.data, '/Organised Data')
dir.res.indicators <- paste0(dir.res, '/Indicators')

## 0.3 data load ----

## shape file ----

poly.lsoa <- sf::st_read(dsn = file.path(dir.data.spatial, 'ons11_LSOA'), layer = 'ons11_LSOA')
poly.lsoa.england <- poly.lsoa %>%  dplyr::filter(str_detect(LSOA11CD, '^E'))
nbs.lsoa.england <- spdep::poly2nb(pl = as(poly.lsoa.england, 'Spatial'))

## look up ----

### ons ----

setwd(dir.data)
lookup.oa.lsoa <- read.csv('spatial/ons_lookup/OA11_LSOA11_MSOA11_LAD11_lookUp.csv')
# lookup.postcode.lsoa <- read.csv('spatial/ons_lookup/POSTCODE_OA11_LSOA11_MSOA11_lookUp.csv')
lookup.lsoa.lad <- read.csv('spatial/ons_lookup/LSOA11_LAD11_RGN11_lookUp.csv')

### nhs ----

lookup.gpPractice.lsoa.temp <- 
  read.csv('nhs/GP_LSOA_PATIENT_COUNT_ALL.csv') %>% 
  dplyr::rename(pc.code = PRACTICE_CODE, 
                pc.name = PRACTICE_NAME,
                LSOA11CD = LSOA_CODE,
                lsoa.gp.total.patients = Number.of.Patients) %>% 
  dplyr::filter(grepl('^E', LSOA11CD))

patients.gpPractise.lsoa <-
  lookup.gpPractice.lsoa.temp %>% 
  dplyr::select(LSOA11CD, pc.code, pc.name, lsoa.gp.total.patients) %>% 
  dplyr::mutate(gp.total.patients = sum(lsoa.gp.total.patients),
                .by = 'pc.code') %>% 
  dplyr::mutate(lsoa.gp.average.patients = lsoa.gp.total.patients/gp.total.patients)


## 0.4 plotting variables ----

height <- width <- 10
text.size <- 20

# colour schemes
colours.diverging <- c('#00768D', '#C13E2E')
colours.exceedance <- c('#A8E6CF', '#5CDB95', '#004D3A')
colours.pollutant <- c('#E9E0C8', '#4B0082')
colours.pollutant.discrete <- grDevices::colorRampPalette(colours.pollutant)(10)
colours.covariate <- c("#E6F2F5", "#006D77")
colours.covariate.discrete <- grDevices::colorRampPalette(colours.covariate)(10)

my.map.theme <- function(...){
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 legend.text = ggplot2::element_text(hjust = 0),
                 legend.key = ggplot2::element_rect(fill = NA, colour = NA),
                 panel.background = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 ...)
}

my.theme<-function(...){
  theme(panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour='black'),
        legend.text = ggplot2::element_text(hjust = 0),
        legend.key = ggplot2::element_rect(fill = NA, colour = NA),
        ...)
}

## indicators ----

setwd(dir.data)
data.import.qof <- read.csv('qof-1819-prev-ach-exc-ls-prac.csv')
data.import.caci.2020 <- read.csv('EW_COA Acorn Directory  Tobacco Spend_April 2020.csv')
data.import.caci.2019 <- read.csv('Output_Area Acorn_Directory_Tobacco_Spend_July_2019.csv')
data.import.imd <- readxl::read_xlsx('LSOA11_2019_IMD_SCORES.xlsx', sheet = 'IoD2019 Scores')
data.import.epd <- readRDS('epd.2018.all.rds')
data.import.ahah <- read.csv('AHAH_V3.csv')
data.import.pollution <- read.csv('cmaqurban_20m_2019_mean_by_LSOA.csv')

# finalise data ----

# 1. descriptive function ----

descriptive.function <- function(x){
  
  # x <- rnorm(n = 1000, mean = 0, sd = 1)
  
  mean <- mean(x)
  names(mean) <- 'mean'
  
  qtl <- stats::quantile(x = x, probs = c(0.25, 0.5, 0.75))
  names(qtl) <- c('lower', 'median', 'upper')
  
  rng <- range(x)
  names(rng) <- c('min', 'max')
  
  sd <- stats::sd(x)
  names(sd) <- 'sd'
  
  c(rng[1], qtl[1], qtl[2], mean, qtl[3], rng[2], sd)
  
}

# 2. qof ----

# data sort
## Stage 1: get to LSOA level
data.sort.qof.1 <-
  data.import.qof %>% 
  dplyr::rename(pc.code = pc_code, 
                pc.name = pc_name,
                list.size = list_size) %>% 
  dplyr::select(pc.code, pc.name, list.size, dplyr::starts_with(c('smok02', 'smok04'))) %>% 
  # number smokers per practise
  dplyr::mutate(total.smokers = sum(smok04_denom),
                .by = pc.code) %>% 
  # add number smokers from each practice to gp x lsoa code
  dplyr::right_join(.,
                    patients.gpPractise.lsoa %>% dplyr::select(-pc.name),
                    by = 'pc.code') %>% 
  # select relavent columns
  dplyr::select(LSOA11CD, pc.code, pc.name, 
                lsoa.gp.total.patients, gp.total.patients, lsoa.gp.average.patients,
                total.smokers, list.size, smok02_score, smok02_denom) %>% 
  # calculate smokers per lsoa
  dplyr::mutate(smokers.lsoa = (lsoa.gp.total.patients / list.size) * total.smokers,
                smokers.list = (lsoa.gp.total.patients / gp.total.patients) * total.smokers) %>% 
  # remove gps where there no smokers
  dplyr::filter(total.smokers > 0) %>% 
  # aggregate
  dplyr::summarise(smokers.lsoa = sum(smokers.lsoa),
                   smokers.list = sum(smokers.list),
                   lsoa.population = sum(lsoa.gp.total.patients),
                   .by = 'LSOA11CD')

## Stage 2: complete lsoa structure
data.sort.qof.2 <- 
  data.frame(LSOA11CD = poly.lsoa.england$LSOA11CD) %>% 
  dplyr::left_join(.,
                   data.sort.qof.1,
                   by = 'LSOA11CD')

## Stage 3: Potential impute missing data
# any missing lsoa
missing.data.qof <- 
  data.sort.qof.2 %>% 
  dplyr::filter(rowSums(is.na(.)) > 0)

# impute missing values from neighborhood structure
if(nrow(missing.data.qof) > 0){
  data.sort.qof.2$smokers.lsoa[is.na(data.sort.qof.2$smokers.lsoa)] <-
    sapply(nbs.lsoa.england[which(is.na(data.sort.qof.2$smokers.lsoa))],
           FUN = function(x) { mean(data.sort.qof.2$smokers.lsoa[x], na.rm = TRUE)})
  
  data.sort.qof.2$smokers.list[is.na(data.sort.qof.2$smokers.list)] <-
    sapply(nbs.lsoa.england[which(is.na(data.sort.qof.2$smokers.list))],
           FUN = function(x) { mean(data.sort.qof.2$smokers.list[x], na.rm = TRUE)})
  
  data.sort.qof.2$lsoa.population[is.na(data.sort.qof.2$lsoa.population)] <-
    sapply(nbs.lsoa.england[which(is.na(data.sort.qof.2$lsoa.population))],
           FUN = function(x) { mean(data.sort.qof.2$lsoa.population[x], na.rm = TRUE)})
  
  setwd(dir.res.indicators)
  print(xtable::xtable(missing.data.qof %>% 
                         dplyr::left_join(lookup.lsoa.lad %>% dplyr::select(LSOA11CD, LSOA11NM), by = 'LSOA11CD') %>% 
                         dplyr::relocate(LSOA11NM, .after = LSOA11CD)),
        include.rownames = FALSE,
        file = 'missing_qof.txt')
}

## Stage 4: Complete score
data.qof <-
  data.sort.qof.2 %>% 
  dplyr::mutate(smoking.prev.lsoa = (smokers.lsoa / lsoa.population * 100 ) %>% round(., digits = 2),
                smoking.prev.list = (smokers.list / lsoa.population * 100 ) %>% round(., digits = 2)) %>% 
  dplyr::mutate(value_qof = smoking.prev.lsoa,
                qof = cut(value_qof, unique(quantile(value_qof, seq(from = 0, to = 1, length.out = 11), na.rm = TRUE)), include.lowest=TRUE),
                qof_id = qof %>% factor(., labels = 1:10)) %>% 
  # join LAD names for LDN extraction
  dplyr::left_join(., lookup.lsoa.lad %>% dplyr::select(LSOA11CD, LAD11CD), by = 'LSOA11CD') %>% 
  # join LSOA polygon
  dplyr::left_join(., poly.lsoa.england, by = 'LSOA11CD') %>% 
  sf::st_as_sf()

## Stage 5. descriptive

setwd(dir.res.indicators)
print(xtable::xtable(data.qof %>% 
                       dplyr::pull(value_qof) %>% 
                       descriptive.function() %>%
                       t() %>%
                       as.data.frame() %>% 
                       dplyr::mutate(misisng  = nrow(missing.data.qof)),
                     digits = 2), 
      include.rownames = FALSE,
      file = 'descriptive_qof.txt')

# check
nrow(data.qof) == nrow(poly.lsoa.england)

# london only data
data.qof.london <-
  data.qof %>% 
  dplyr::filter(str_detect(LAD11CD, "^E09"))

# plot
limits.spatial <- c(0, 55)
breaks.spatial <- seq(from = limits.spatial[1], to = limits.spatial[2], by = limits.spatial[2]/4)

# plot
## continuous
map.london.qof <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.qof.london, aes(fill = value_qof), colour = NA) +
  ggplot2::scale_fill_gradient(name = '',
                               low = colours.covariate[1], high = colours.covariate[2],
                               limits = limits.spatial, breaks = breaks.spatial) +
  my.map.theme(legend.position = 'none')
map.qof <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.qof, aes(fill = value_qof), colour = NA) +
  ggplot2::scale_fill_gradient(name = 'Smokers per 1000 people',
                               guide = guide_colourbar(title.position = 'top'),
                               low = colours.covariate[1], high = colours.covariate[2],
                               limits = limits.spatial, breaks = breaks.spatial) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.qof), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5))
## discrete
map.london.qof.decile <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.qof.london, aes(fill = qof_id), colour = NA) +
  ggplot2::scale_fill_manual(name = '',
                             values = colours.covariate.discrete) +
  my.map.theme(legend.position = 'none')
map.qof.decile <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.qof, aes(fill = qof_id), colour = NA) +
  ggplot2::scale_fill_manual(name = 'Smokers per 1000 people',
                             guide = ggplot2::guide_legend(title.position = 'top'),
                             labels = c('1 (Least)', 2:9, '10 (Most)'),
                             values = colours.covariate.discrete) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.qof.decile), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  # ggplot2::labs(title = 'Pollution: 03') +
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5)); (map.qof/map.qof.decile)
## historgram
breaks.qof <-
  data.qof %>% 
  dplyr::pull(value_qof) %>% 
  quantile(x = ., probs = seq(from = 0, to = 1, length.out = 11), na.rm = TRUE) %>% 
  unique()
histogram.qof <-
  ggplot2::ggplot(data = data.qof %>% sf::st_drop_geometry(), aes(x = value_qof)) +
  ggplot2::geom_histogram(colour = 'black', fill = 'white', bins = 250) +
  ggplot2::geom_vline(xintercept = breaks.qof, colour = 'red3', linetype = 'dashed', linewidth = 0.75) +
  ggplot2::labs(x = 'Smokers per 1000 people', y = 'Count') +
  my.theme(text = element_text(size = text.size))

# save
setwd(dir.res.indicators)
ggplot2::ggsave(plot = map.qof,
                filename = 'map_qof.png',
                width = width, 
                height = height)
ggplot2::ggsave(plot = map.qof.decile,
                filename = 'map_qof_decile.png',
                width = width, 
                height = height)
ggplot2::ggsave(plot = histogram.qof,
                filename = 'histogram_qof.png',
                width = width, 
                height = height)

# 3. epd ----

# data sort
## Stage 1: get to LSOA level
data.sort.epd.1 <- 
  data.import.epd %>% 
  dplyr::rename(pc.code = PRACTICE_CODE) %>% 
  dplyr::summarise(total = sum(TOTAL_QUANTITY),
                   .by = pc.code) %>% 
  # add total from each practice to gp x lsoa code
  dplyr::right_join(.,
                    patients.gpPractise.lsoa %>% dplyr::select(-pc.name),
                    by = 'pc.code') %>% 
  dplyr::mutate(total = dplyr::if_else(is.na(total), 0, total)) %>% 
  # rate per lsoa
  dplyr::mutate(rate.lsoa = lsoa.gp.total.patients / gp.total.patients * total) %>% 
  # check should be same as total
  dplyr::mutate(check =  sum(rate.lsoa),
                .by = 'pc.code') %>% 
  dplyr::summarise(rate.lsoa = sum(rate.lsoa),
                   lsoa.gp.total.patients = sum(lsoa.gp.total.patients),
                   .by = 'LSOA11CD')

## Stage 2: complete lsoa structure
data.sort.epd.2 <- 
  data.frame(LSOA11CD = poly.lsoa.england$LSOA11CD) %>% 
  dplyr::left_join(.,
                   data.sort.epd.1,
                   by = 'LSOA11CD')

## Stage 3: Potential impute missing data
# any missing lsoa
missing.data.epd <- 
  data.sort.epd.2 %>% 
  dplyr::filter(rowSums(is.na(.)) > 0)

# impute missing values from neighborhood structure
if(nrow(missing.data.epd) > 0){
  data.sort.epd.2$rate.lsoa[is.na(data.sort.epd.2$rate.lsoa)] <-
    sapply(nbs.lsoa.england[which(is.na(data.sort.epd.2$rate.lsoa))],
           FUN = function(x) { mean(data.sort.epd.2$rate.lsoa[x], na.rm = TRUE)})
  
  
  data.sort.epd.2$lsoa.gp.total.patients[is.na(data.sort.epd.2$lsoa.gp.total.patients)] <-
    sapply(nbs.lsoa.england[which(is.na(data.sort.epd.2$lsoa.gp.total.patients))],
           FUN = function(x) { mean(data.sort.epd.2$lsoa.gp.total.patients[x], na.rm = TRUE)})
  
  setwd(dir.res.indicators)
  print(xtable::xtable(missing.data.epd %>% 
                         dplyr::left_join(lookup.lsoa.lad %>% dplyr::select(LSOA11CD, LSOA11NM), by = 'LSOA11CD') %>% 
                         dplyr::relocate(LSOA11NM, .after = LSOA11CD)),
        include.rownames = FALSE,
        file = 'missing_epd.txt')
}

## Stage 4: Complete score
data.epd <-
  data.sort.epd.2 %>% 
  dplyr::mutate(rate.person.lsoa = rate.lsoa / lsoa.gp.total.patients) %>% 
  dplyr::mutate(value_epd = rate.person.lsoa,
                epd = cut(value_epd, unique(quantile(value_epd, seq(from = 0, to = 1, length.out = 11), na.rm = TRUE)), include.lowest=TRUE),
                epd_id = epd %>% factor(., labels = 1:10)) %>% 
  # join LAD names for LDN extraction
  dplyr::left_join(., lookup.lsoa.lad %>% dplyr::select(LSOA11CD, LAD11CD), by = 'LSOA11CD') %>% 
  # join LSOA polygon
  dplyr::left_join(., poly.lsoa.england, by = 'LSOA11CD') %>% 
  sf::st_as_sf()

## Stage 5. descriptive

setwd(dir.res.indicators)
print(xtable::xtable(data.epd %>% 
                       dplyr::pull(value_epd) %>% 
                       descriptive.function() %>%
                       t() %>%
                       as.data.frame() %>% 
                       dplyr::mutate(misisng  = nrow(missing.data.epd)),
                     digits = 2), 
      include.rownames = FALSE,
      file = 'descriptive_epd.txt')

# check
nrow(data.epd) == nrow(poly.lsoa.england)

# london only data
data.epd.london <-
  data.epd %>% 
  dplyr::filter(str_detect(LAD11CD, "^E09"))

# plot
limits.spatial <- c(0, 7)
breaks.spatial <- seq(from = limits.spatial[1], to = limits.spatial[2], by = limits.spatial[2]/4)

# plot
## continuous
map.london.epd <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.epd.london, aes(fill = value_epd), colour = NA) +
  ggplot2::scale_fill_gradient(name = '',
                               low = colours.covariate[1], high = colours.covariate[2],
                               limits = limits.spatial, breaks = breaks.spatial) +
  my.map.theme(legend.position = 'none')
map.epd <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.epd, aes(fill = value_epd), colour = NA) +
  ggplot2::scale_fill_gradient(name = 'Nicotine replacement therapy\nprescriptions per person',
                               guide = guide_colourbar(title.position = 'top'),
                               low = colours.covariate[1], high = colours.covariate[2],
                               limits = limits.spatial, breaks = breaks.spatial) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.epd), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5))
## discrete
map.london.epd.decile <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.epd.london, aes(fill = epd_id), colour = NA) +
  ggplot2::scale_fill_manual(name = '',
                             values = colours.covariate.discrete) +
  my.map.theme(legend.position = 'none')
map.epd.decile <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.epd, aes(fill = epd_id), colour = NA) +
  ggplot2::scale_fill_manual(name = 'Nicotine replacement therapy\nprescriptions per person',
                             guide = ggplot2::guide_legend(title.position = 'top'),
                             labels = c('1 (Least)', 2:9, '10 (Most)'),
                             values = colours.covariate.discrete) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.epd.decile), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  # ggplot2::labs(title = 'Pollution: 03') +
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5)); (map.epd/map.epd.decile)
## historgram
breaks.epd <-
  data.epd %>% 
  dplyr::pull(value_epd) %>% 
  quantile(x = ., probs = seq(from = 0, to = 1, length.out = 11), na.rm = TRUE) %>% 
  unique()
histogram.epd <-
  ggplot2::ggplot(data = data.epd %>% sf::st_drop_geometry(), aes(x = value_epd)) +
  ggplot2::geom_histogram(colour = 'black', fill = 'white', bins = 250) +
  ggplot2::geom_vline(xintercept = breaks.epd, colour = 'red3', linetype = 'dashed', linewidth = 0.75) +
  ggplot2::labs(x = 'Nicotine replacement therapy\nprescriptions per person', y = 'Count') +
  my.theme(text = element_text(size = text.size))

# save
setwd(dir.res.indicators)
ggplot2::ggsave(plot = map.epd,
                filename = 'map_epd.png',
                width = width, 
                height = height)
ggplot2::ggsave(plot = map.epd.decile,
                filename = 'map_epd_decile.png',
                width = width, 
                height = height)
ggplot2::ggsave(plot = histogram.epd,
                filename = 'histogram_epd.png',
                width = width, 
                height = height)

# 4. caci ----

# data sort
## Stage 1: get to LSOA level
data.sort.caci.1 <-
  data.import.caci.2019 %>% 
  dplyr::rename(OA11CD = Output.Area.ID,
                tobacco.spend = Tobacco.Expenditure..pw.) %>% 
  dplyr::select(OA11CD, tobacco.spend, Total.Population) %>%
  # filter so there are no 0 spend and 0 population
  dplyr::filter(tobacco.spend > 0, Total.Population > 0) %>% 
  # join LSOA
  dplyr::left_join(.,
                   lookup.oa.lsoa %>% dplyr::select(OA11CD, LSOA11CD),
                   by = 'OA11CD') %>% 
  # england only
  dplyr::filter(grepl('^E', LSOA11CD)) %>% 
  # LSOA level analysis  
  dplyr::summarise(tobacco.spend = sum(tobacco.spend),
                   population = sum(Total.Population),
                   .by = 'LSOA11CD')

## Stage 2: complete lsoa structure
data.sort.caci.2 <- 
  data.frame(LSOA11CD = poly.lsoa.england$LSOA11CD) %>% 
  dplyr::left_join(.,
                   data.sort.caci.1,
                   by = 'LSOA11CD')

## Stage 3: Potential impute missing data
# any missing lsoa
missing.data.caci <- 
  data.sort.caci.2 %>% 
  dplyr::filter(rowSums(is.na(.)) > 0)

# impute missing values from neighborhood structure
if(nrow(missing.data.caci) > 0){
  data.sort.caci.2$tobacco.spend[is.na(data.sort.caci.2$tobacco.spend)] <-
    sapply(nbs.lsoa.england[which(is.na(data.sort.caci.2$tobacco.spend))],
           FUN = function(x) { mean(data.sort.caci.2$tobacco.spend[x], na.rm = TRUE)})
  
  data.sort.caci.2$population[is.na(data.sort.caci.2$population)] <-
    sapply(nbs.lsoa.england[which(is.na(data.sort.caci.2$population))],
           FUN = function(x) { mean(data.sort.caci.2$population[x], na.rm = TRUE)})
  
  setwd(dir.res.indicators)
  print(xtable::xtable(missing.data.caci %>% 
                         dplyr::left_join(lookup.lsoa.lad %>% dplyr::select(LSOA11CD, LSOA11NM), by = 'LSOA11CD') %>% 
                         dplyr::relocate(LSOA11NM, .after = LSOA11CD)),
        include.rownames = FALSE,
        file = 'missing_caci.txt')
}

## Stage 4: Complete score
data.caci <-
  data.sort.caci.2 %>% 
  # spending per pop with decile
  dplyr::mutate(spending.rate = tobacco.spend / population) %>% 
  dplyr::mutate(value_caci = spending.rate,
                caci = cut(value_caci, unique(quantile(value_caci, seq(from = 0, to = 1, length.out = 11), na.rm = TRUE)), include.lowest=TRUE),
                caci_id = caci %>% factor(., labels = 1:10)) %>% 
  # join LAD names for LDN extraction
  dplyr::left_join(., lookup.lsoa.lad %>% dplyr::select(LSOA11CD, LAD11CD), by = 'LSOA11CD') %>% 
  # join LSOA polygon
  dplyr::left_join(., poly.lsoa.england, by = 'LSOA11CD') %>% 
  sf::st_as_sf()

## Stage 5. descriptive

setwd(dir.res.indicators)
print(xtable::xtable(data.caci %>% 
                       dplyr::pull(value_caci) %>% 
                       descriptive.function() %>%
                       t() %>%
                       as.data.frame() %>% 
                       dplyr::mutate(misisng  = nrow(missing.data.caci)),
                     digits = 2), 
      include.rownames = FALSE,
      file = 'descriptive_caci.txt')

# check
nrow(data.caci) == nrow(poly.lsoa.england)

# london only data
data.caci.london <-
  data.caci %>% 
  dplyr::filter(str_detect(LAD11CD, "^E09"))

# plot
limits.spatial <- c(0, 15)
breaks.spatial <- seq(from = limits.spatial[1], to = limits.spatial[2], by = limits.spatial[2]/4)

# plot
## continuous
map.london.caci <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.caci.london, aes(fill = value_caci), colour = NA) +
  ggplot2::scale_fill_gradient(name = '',
                               low = colours.covariate[1], high = colours.covariate[2],
                               limits = limits.spatial, breaks = breaks.spatial) +
  my.map.theme(legend.position = 'none')
map.caci <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.caci, aes(fill = value_caci), colour = NA) +
  ggplot2::scale_fill_gradient(name = 'Tobacco spending per person per week (£)',
                               guide = guide_colourbar(title.position = 'top'),
                               low = colours.covariate[1], high = colours.covariate[2],
                               limits = limits.spatial, breaks = breaks.spatial) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.caci), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5))
## discrete
map.london.caci.decile <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.caci.london, aes(fill = caci_id), colour = NA) +
  ggplot2::scale_fill_manual(name = '',
                             values = colours.covariate.discrete) +
  my.map.theme(legend.position = 'none')
map.caci.decile <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.caci, aes(fill = caci_id), colour = NA) +
  ggplot2::scale_fill_manual(name = 'Tobacco spending per person per week (£)',
                             guide = ggplot2::guide_legend(title.position = 'top'),
                             labels = c('1 (Least)', 2:9, '10 (Most)'),
                             values = colours.covariate.discrete) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.caci.decile), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  # ggplot2::labs(title = 'Pollution: 03') +
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5)); (map.caci/map.caci.decile)
## historgram
breaks.caci <-
  data.caci %>% 
  dplyr::pull(value_caci) %>% 
  quantile(x = ., probs = seq(from = 0, to = 1, length.out = 11), na.rm = TRUE) %>% 
  unique()
histogram.caci <-
  ggplot2::ggplot(data = data.caci %>% sf::st_drop_geometry(), aes(x = value_caci)) +
  ggplot2::geom_histogram(colour = 'black', fill = 'white', bins = 250) +
  ggplot2::geom_vline(xintercept = breaks.caci, colour = 'red3', linetype = 'dashed', linewidth = 0.75) +
  ggplot2::labs(x = 'Tobacco spending per person per week (£)', y = 'Count') +
  my.theme(text = element_text(size = text.size))

# save
setwd(dir.res.indicators)
ggplot2::ggsave(plot = map.caci,
                filename = 'map_caci.png',
                width = width, 
                height = height)
ggplot2::ggsave(plot = map.caci.decile,
                filename = 'map_caci_decile.png',
                width = width, 
                height = height)
ggplot2::ggsave(plot = histogram.caci,
                filename = 'histogram_caci.png',
                width = width, 
                height = height)

# 5. ahah ----

# data sort
## Stage 1: get to LSOA level
data.sort.ahah.1 <- 
  data.import.ahah %>% 
  dplyr::rename(LSOA11CD = lsoa11) %>% 
  dplyr::select(LSOA11CD, ah3tob, ah3tob_rnk, ah3tob_pct) %>% 
  # filter for england only
  dplyr::filter(grepl('^E', LSOA11CD))

## Stage 2: complete lsoa structure
data.sort.ahah.2 <- 
  data.frame(LSOA11CD = poly.lsoa.england$LSOA11CD) %>% 
  dplyr::left_join(.,
                   data.sort.ahah.1,
                   by = 'LSOA11CD')

## Stage 3: Potential impute missing data
# any missing lsoa
missing.data.ahah <- 
  data.sort.ahah.2 %>% 
  dplyr::filter(rowSums(is.na(.)) > 0)

# impute missing values from neighborhood structure
if(nrow(missing.data.ahah) > 0){
  data.sort.ahah.2$ah3tob[is.na(data.sort.ahah.2$ah3tob)] <-
    sapply(nbs.lsoa.england[which(is.na(data.sort.ahah.2$ah3tob))],
           FUN = function(x) { mean(data.sort.ahah.2$ah3tob[x], na.rm = TRUE)})
  
  setwd(dir.res.indicators)
  print(xtable::xtable(missing.data.ahah %>% 
                         dplyr::left_join(lookup.lsoa.lad %>% dplyr::select(LSOA11CD, LSOA11NM), by = 'LSOA11CD') %>% 
                         dplyr::relocate(LSOA11NM, .after = LSOA11CD)),
        include.rownames = FALSE,
        file = 'missing_ahah.txt')
}

## Stage 4: Complete score
data.ahah <-
  data.sort.ahah.2 %>% 
  dplyr::mutate(value_ahah = ah3tob,
                ahah = cut(value_ahah, unique(quantile(value_ahah, seq(from = 0, to = 1, length.out = 11), na.rm = TRUE)), include.lowest=TRUE),
                ahah_id = ahah %>% factor(., labels = 1:10)) %>% 
  # join LAD names for LDN extraction
  dplyr::left_join(., lookup.lsoa.lad %>% dplyr::select(LSOA11CD, LAD11CD), by = 'LSOA11CD') %>% 
  # join LSOA polygon
  dplyr::left_join(., poly.lsoa.england, by = 'LSOA11CD') %>% 
  sf::st_as_sf()

## Stage 5. descriptive

setwd(dir.res.indicators)
print(xtable::xtable(data.ahah %>% 
                       dplyr::pull(value_ahah) %>% 
                       descriptive.function() %>%
                       t() %>%
                       as.data.frame() %>% 
                       dplyr::mutate(misisng  = nrow(missing.data.ahah)),
                     digits = 2), 
      include.rownames = FALSE,
      file = 'descriptive_ahah.txt')

# check
nrow(data.ahah) == nrow(poly.lsoa.england)

# london only data
data.ahah.london <-
  data.ahah %>% 
  dplyr::filter(str_detect(LAD11CD, "^E09"))

# plot
limits.spatial <- c(0, 170)
breaks.spatial <- seq(from = limits.spatial[1], to = limits.spatial[2], by = limits.spatial[2]/4)

# plot
## continuous
map.london.ahah <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.ahah.london, aes(fill = value_ahah), colour = NA) +
  ggplot2::scale_fill_gradient(name = '',
                               low = colours.covariate[1], high = colours.covariate[2],
                               limits = limits.spatial, breaks = breaks.spatial) +
  my.map.theme(legend.position = 'none')
map.ahah <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.ahah, aes(fill = value_ahah), colour = NA) +
  ggplot2::scale_fill_gradient(name = 'Distance to tobacconist (km)',
                               guide = guide_colourbar(title.position = 'top'),
                               low = colours.covariate[1], high = colours.covariate[2],
                               limits = limits.spatial, breaks = breaks.spatial) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.ahah), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5))
## discrete
map.london.ahah.decile <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.ahah.london, aes(fill = ahah_id), colour = NA) +
  ggplot2::scale_fill_manual(name = '',
                             values = colours.covariate.discrete) +
  my.map.theme(legend.position = 'none')
map.ahah.decile <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.ahah, aes(fill = ahah_id), colour = NA) +
  ggplot2::scale_fill_manual(name = 'Distance to tobacconist (km)',
                             guide = ggplot2::guide_legend(title.position = 'top'),
                             labels = c('1 (Nearest)', 2:9, '10 (Furthest)'),
                             values = colours.covariate.discrete) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.ahah.decile), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  # ggplot2::labs(title = 'Pollution: 03') +
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5)); (map.ahah/map.ahah.decile)
## historgram
breaks.ahah <-
  data.ahah %>% 
  dplyr::pull(value_ahah) %>% 
  quantile(x = ., probs = seq(from = 0, to = 1, length.out = 11), na.rm = TRUE) %>% 
  unique()
histogram.ahah <-
  ggplot2::ggplot(data = data.ahah %>% sf::st_drop_geometry(), aes(x = value_ahah)) +
  ggplot2::geom_histogram(colour = 'black', fill = 'white', bins = 250) +
  ggplot2::geom_vline(xintercept = breaks.ahah, colour = 'red3', linetype = 'dashed', linewidth = 0.75) +
  ggplot2::labs(x = 'Distance to tobacconist (km)', y = 'Count') +
  my.theme(text = element_text(size = text.size))

# save
setwd(dir.res.indicators)
ggplot2::ggsave(plot = map.ahah,
                filename = 'map_ahah.png',
                width = width, 
                height = height)
ggplot2::ggsave(plot = map.ahah.decile,
                filename = 'map_ahah_decile.png',
                width = width, 
                height = height)
ggplot2::ggsave(plot = histogram.ahah,
                filename = 'histogram_ahah.png',
                width = width, 
                height = height)

# 6. imd ----

# data sort
## Stage 1: get to LSOA level
data.sort.imd.1 <- 
  data.import.imd %>% 
  dplyr::mutate(LSOA11CD = `LSOA code (2011)`,
                imdScore = `Index of Multiple Deprivation (IMD) Score`) %>% 
  dplyr::select(LSOA11CD, imdScore)

## Stage 2: complete lsoa structure
data.sort.imd.2 <- 
  data.frame(LSOA11CD = poly.lsoa.england$LSOA11CD) %>% 
  dplyr::left_join(.,
                   data.sort.imd.1,
                   by = 'LSOA11CD')

## Stage 3: Potential impute missing data
# any missing lsoa
missing.data.imd <- 
  data.sort.imd.2 %>% 
  dplyr::filter(rowSums(is.na(.)) > 0)

# imput missing values from neighborhood structure
if(nrow(missing.data.imd) > 0){
  data.sort.imd.2$imdScore[is.na(data.sort.imd.2$imdScore)] <-
    sapply(nbs.lsoa.england[which(is.na(data.sort.imd.2$imdScore))],
           FUN = function(x) { mean(data.sort.imd.2$imdScore[x], na.rm = TRUE)})
  
  setwd(dir.res.indicators)
  print(xtable::xtable(missing.imd %>% 
                         dplyr::left_join(lookup.lsoa.lad %>% dplyr::select(LSOA11CD, LSOA11NM), by = 'LSOA11CD') %>% 
                         dplyr::relocate(LSOA11NM, .after = LSOA11CD)),
        include.rownames = FALSE,
        file = 'missing_imd.txt')
}

## Stage 4: Complete score
data.imd <-
  data.sort.imd.2 %>% 
  dplyr::mutate(value_imd = imdScore,
                imd = cut(value_imd, unique(quantile(value_imd, seq(from = 0, to = 1, length.out = 11), na.rm = TRUE)), include.lowest=TRUE),
                imd_id = imd %>% factor(., labels = 1:10)) %>% 
  # join LAD names for LDN extraction
  dplyr::left_join(., lookup.lsoa.lad %>% dplyr::select(LSOA11CD, LAD11CD), by = 'LSOA11CD') %>% 
  # join LSOA polygon
  dplyr::left_join(., poly.lsoa.england, by = 'LSOA11CD') %>% 
  sf::st_as_sf()

## Stage 5. descriptive

setwd(dir.res.indicators)
print(xtable::xtable(data.imd %>% 
                       dplyr::pull(value_imd) %>% 
                       descriptive.function() %>%
                       t() %>%
                       as.data.frame() %>% 
                       dplyr::mutate(misisng  = nrow(missing.data.imd)),
                     digits = 2), 
      include.rownames = FALSE,
      file = 'descriptive_imd.txt')
# check
nrow(data.imd) == nrow(poly.lsoa.england)

# london only data
data.imd.london <-
  data.imd %>% 
  dplyr::filter(str_detect(LAD11CD, "^E09"))

# plot
limits.spatial <- c(0, 100)
breaks.spatial <- seq(from = limits.spatial[1], to = limits.spatial[2], by = limits.spatial[2]/4)

# plot
## continuous
map.london.imd <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.imd.london, aes(fill = value_imd), colour = NA) +
  ggplot2::scale_fill_gradient(name = '',
                               low = colours.covariate[1], high = colours.covariate[2],
                               limits = limits.spatial, breaks = breaks.spatial) +
  my.map.theme(legend.position = 'none')
map.imd <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.imd, aes(fill = value_imd), colour = NA) +
  ggplot2::scale_fill_gradient(name = 'Index of multiple deprivation score',
                               guide = guide_colourbar(title.position = 'top'),
                               low = colours.covariate[1], high = colours.covariate[2],
                               limits = limits.spatial, breaks = breaks.spatial) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.imd), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5))
## discrete
map.london.imd.decile <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.imd.london, aes(fill = imd_id), colour = NA) +
  ggplot2::scale_fill_manual(name = '',
                             values = colours.covariate.discrete) +
  my.map.theme(legend.position = 'none')
map.imd.decile <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.imd, aes(fill = imd_id), colour = NA) +
  ggplot2::scale_fill_manual(name = 'Index of multiple deprivation score',
                             guide = ggplot2::guide_legend(title.position = 'top'),
                             labels = c('1 (Lowest)', 2:9, '10 (Highest)'),
                             values = colours.covariate.discrete) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.imd.decile), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  # ggplot2::labs(title = 'Pollution: 03') +
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5)); (map.imd/map.imd.decile)
## historgram
breaks.imd <-
  data.imd %>% 
  dplyr::pull(value_imd) %>% 
  quantile(x = ., probs = seq(from = 0, to = 1, length.out = 11), na.rm = TRUE) %>% 
  unique()
histogram.imd <-
  ggplot2::ggplot(data = data.imd %>% sf::st_drop_geometry(), aes(x = value_imd)) +
  ggplot2::geom_histogram(colour = 'black', fill = 'white', bins = 250) +
  ggplot2::geom_vline(xintercept = breaks.imd, colour = 'red3', linetype = 'dashed', linewidth = 0.75) +
  ggplot2::labs(x = 'Index of multiple deprivation score', y = 'Count') +
  my.theme(text = element_text(size = text.size))

# save
setwd(dir.res.indicators)
ggplot2::ggsave(plot = map.imd,
                filename = 'map_imd.png',
                width = width, 
                height = height)
ggplot2::ggsave(plot = map.imd.decile,
                filename = 'map_imd_decile.png',
                width = width, 
                height = height)
ggplot2::ggsave(plot = histogram.imd,
                filename = 'histogram_imd.png',
                width = width, 
                height = height)

# 7. pollution -----

# data sort
## Stage 1: get to LSOA level
data.sort.pollution.1 <- 
  data.import.pollution %>% 
  # remove name as from poly
  dplyr::select(-LSOA11NM, -GlobalID)

## Stage 2: complete lsoa structure
data.sort.pollution.2 <- 
  data.frame(LSOA11CD = poly.lsoa.england$LSOA11CD) %>% 
  dplyr::left_join(.,
                   data.sort.pollution.1,
                   by = 'LSOA11CD')

## Stage 3: Potential impute missing data
# any missing lsoa
missing.data.pollution <- 
  data.sort.pollution.2 %>% 
  dplyr::filter(rowSums(is.na(.)) > 0)

# impute missing values from neighborhood structure
if(nrow(missing.data.pollution) > 0){
  data.sort.pollution.2$PM25_ugm3[is.na(data.sort.pollution.2$PM25_ugm3)] <-
    sapply(nbs.lsoa.england[which(is.na(data.sort.pollution.2$PM25_ugm3))],
           FUN = function(x) { mean(data.sort.pollution.2$PM25_ugm3[x], na.rm = TRUE)})
  
  data.sort.pollution.2$PM10_ugm3[is.na(data.sort.pollution.2$PM10_ugm3)] <-
    sapply(nbs.lsoa.england[which(is.na(data.sort.pollution.2$PM10_ugm3))],
           FUN = function(x) { mean(data.sort.pollution.2$PM10_ugm3[x], na.rm = TRUE)})
  
  data.sort.pollution.2$NO2_ugm3[is.na(data.sort.pollution.2$NO2_ugm3)] <-
    sapply(nbs.lsoa.england[which(is.na(data.sort.pollution.2$NO2_ugm3))],
           FUN = function(x) { mean(data.sort.pollution.2$NO2_ugm3[x], na.rm = TRUE)})
  
  data.sort.pollution.2$O3_ugm3[is.na(data.sort.pollution.2$O3_ugm3)] <-
    sapply(nbs.lsoa.england[which(is.na(data.sort.pollution.2$O3_ugm3))],
           FUN = function(x) { mean(data.sort.pollution.2$O3_ugm3[x], na.rm = TRUE)})
  
  setwd(dir.res.indicators)
  print(xtable::xtable(missing.data.pollution %>% 
                         dplyr::left_join(lookup.lsoa.lad %>% dplyr::select(LSOA11CD, LSOA11NM), by = 'LSOA11CD') %>% 
                         dplyr::relocate(LSOA11NM, .after = LSOA11CD)),
        include.rownames = FALSE,
        file = 'missing_pollution.txt')
}

## Stage 4: Complete score
data.pollution <-
  data.sort.pollution.2 %>%
  # nicer naming
  dplyr::rename(pm25 = PM25_ugm3,
                pm10 = PM10_ugm3,
                no2 = NO2_ugm3,
                o3 = O3_ugm3) %>% 
  dplyr::mutate(pm25_decile = dplyr::ntile(x = pm25, n = 10) %>% factor(x = ., levels = 1:10, labels = c('1 (Lowest)', 2:9, '10 (Highest)')),
                pm10_decile = dplyr::ntile(x = pm10, n = 10) %>% factor(x = ., levels = 1:10, labels = c('1 (Lowest)', 2:9, '10 (Highest)')),
                no2_decile = dplyr::ntile(x = no2, n = 10) %>% factor(x = ., levels = 1:10, labels = c('1 (Lowest)', 2:9, '10 (Highest)')),
                o3_decile = dplyr::ntile(x = o3, n = 10) %>% factor(x = ., levels = 1:10, labels = c('1 (Lowest)', 2:9, '10 (Highest)'))) %>% 
  # join LAD names for LDN extraction
  dplyr::left_join(., lookup.lsoa.lad %>% dplyr::select(LSOA11CD, LAD11CD), by = 'LSOA11CD') %>% 
  # join LSOA polygon
  dplyr::left_join(., poly.lsoa.england, by = 'LSOA11CD') %>% 
  sf::st_as_sf()

## Stage 5. descriptive

setwd(dir.res.indicators)
print(xtable::xtable(data.pollution %>% 
                       dplyr::select(pm25, pm10, no2, o3) %>%
                       sf::st_drop_geometry() %>% 
                       purrr::map_df(.f = ~ descriptive.function(.x), .id = 'variable') %>%
                       dplyr::mutate(missing = 0),
                     digits = 2), 
      include.rownames = FALSE,
      file = 'descriptive_pollution.txt')

# check
nrow(data.pollution) == nrow(poly.lsoa.england)
  
# london only data
data.pollution.london <-
  data.pollution %>% 
  dplyr::filter(str_detect(LAD11CD, "^E09"))

## 7.1. pm25 ----

# limits
limits.pm25 <- c(0, 20)

# breaks
breaks.pm25 <- seq(from = limits.pm25[1], to = limits.pm25[2], by = limits.pm25[2]/4)

# continous
map.london.pm25 <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.pollution.london, aes(fill = pm25), colour = NA) +
  ggplot2::scale_fill_gradient(name = '',
                               low = colours.pollutant[1], high = colours.pollutant[2],
                               limits = limits.pm25, breaks = breaks.pm25) +
  my.map.theme(legend.position = 'none')
map.pm25 <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.pollution, aes(fill = pm25), colour = NA) +
  ggplot2::scale_fill_gradient(name = unname(latex2exp ::TeX(c('$PM_{2.5}(\\mu{g}/m^{3}$)'))),
                               guide = guide_colourbar(title.position = 'top'),
                               low = colours.pollutant[1], high = colours.pollutant[2],
                               limits = limits.pm25, breaks = breaks.pm25) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.pm25), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  # ggplot2::labs(title = 'Pollution: PM2.5') +
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5))

# discrete
map.london.pm25.decile <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.pollution.london, aes(fill = pm25_decile), colour = NA) +
  ggplot2::scale_fill_manual(name = '',
                             values = colours.pollutant.discrete) +
  my.map.theme(legend.position = 'none')
map.pm25.decile <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.pollution, aes(fill = pm25_decile), colour = NA) +
  ggplot2::scale_fill_manual(name = unname(latex2exp ::TeX(c('$PM_{2.5}(\\mu{g}/m^{3}$) Decile'))),
                             guide = ggplot2::guide_legend(title.position = 'top'),
                             values = colours.pollutant.discrete) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.pm25.decile), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  # ggplot2::labs(title = 'Pollution: PM2.5') +
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5)); 

# save
setwd(dir.res.indicators)
ggplot2::ggsave(plot = map.pm25,
                filename = 'map_pm25.png',
                width = width, 
                height = height)
ggplot2::ggsave(plot = map.pm25.decile,
                filename = 'map_pm25_decile.png',
                width = width, 
                height = height)

## 7.2. pm10 ----

# limits
limits.pm10 <- c(0, 35)

# breaks
breaks.pm10 <- seq(from = limits.pm10[1], to = limits.pm10[2], by = limits.pm10[2]/4)

# continous
map.london.pm10 <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.pollution.london, aes(fill = pm10), colour = NA) +
  ggplot2::scale_fill_gradient(name = '',
                               low = colours.pollutant[1], high = colours.pollutant[2],
                               limits = limits.pm10, breaks = breaks.pm10) +
  my.map.theme(legend.position = 'none')
map.pm10 <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.pollution, aes(fill = pm10), colour = NA) +
  ggplot2::scale_fill_gradient(name = unname(latex2exp ::TeX(c('$PM_{10}(\\mu{g}/m^{3}$)'))),
                               guide = guide_colourbar(title.position = 'top'),
                               low = colours.pollutant[1], high = colours.pollutant[2],
                               limits = limits.pm10, breaks = breaks.pm10) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.pm10), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  # ggplot2::labs(title = 'Pollution: PM10') +
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5))

# discrete
map.london.pm10.decile <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.pollution.london, aes(fill = pm10_decile), colour = NA) +
  ggplot2::scale_fill_manual(name = '',
                             values = colours.pollutant.discrete) +
  my.map.theme(legend.position = 'none')
map.pm10.decile <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.pollution, aes(fill = pm10_decile), colour = NA) +
  ggplot2::scale_fill_manual(name = unname(latex2exp ::TeX(c('$PM_{10}(\\mu{g}/m^{3}$) Decile'))),
                             guide = ggplot2::guide_legend(title.position = 'top'),
                             values = colours.pollutant.discrete) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.pm10.decile), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  # ggplot2::labs(title = 'Pollution: PM10') +
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5)); 

# save
setwd(dir.res.indicators)
ggplot2::ggsave(plot = map.pm10,
                filename = 'map_pm10.png',
                width = width, 
                height = height)
ggplot2::ggsave(plot = map.pm10.decile,
                filename = 'map_pm10_decile.png',
                width = width, 
                height = height)

## 7.3. no2 ----

# limits
limits.no2 <- c(0, 60)

# breaks
breaks.no2 <- seq(from = limits.no2[1], to = limits.no2[2], by = limits.no2[2]/4)

# continous
map.london.no2 <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.pollution.london, aes(fill = no2), colour = NA) +
  ggplot2::scale_fill_gradient(name = '',
                               low = colours.pollutant[1], high = colours.pollutant[2],
                               limits = limits.no2, breaks = breaks.no2) +
  my.map.theme(legend.position = 'none')
map.no2 <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.pollution, aes(fill = no2), colour = NA) +
  ggplot2::scale_fill_gradient(name = unname(latex2exp ::TeX(c('$NO_{2}(\\mu{g}/m^{3}$)'))),
                               guide = guide_colourbar(title.position = 'top'),
                               low = colours.pollutant[1], high = colours.pollutant[2],
                               limits = limits.no2, breaks = breaks.no2) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.no2), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  # ggplot2::labs(title = 'Pollution: NO2') +
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5))

# discrete
map.london.no2.decile <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.pollution.london, aes(fill = no2_decile), colour = NA) +
  ggplot2::scale_fill_manual(name = '',
                             values = colours.pollutant.discrete) +
  my.map.theme(legend.position = 'none')
map.no2.decile <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.pollution, aes(fill = no2_decile), colour = NA) +
  ggplot2::scale_fill_manual(name = unname(latex2exp ::TeX(c('$NO_{2}(\\mu{g}/m^{3}$) Decile'))),
                             guide = ggplot2::guide_legend(title.position = 'top'),
                             values = colours.pollutant.discrete) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.no2.decile), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  # ggplot2::labs(title = 'Pollution: NO2') +
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5)); 

# save
setwd(dir.res.indicators)
ggplot2::ggsave(plot = map.no2,
                filename = 'map_no2.png',
                width = width, 
                height = height)
ggplot2::ggsave(plot = map.no2.decile,
                filename = 'map_no2_decile.png',
                width = width, 
                height = height)

## 7.4. o3 ----

# limits
limits.o3 <- c(0, 65)

# breaks
breaks.o3 <- seq(from = limits.o3[1], to = limits.o3[2], by = limits.o3[2]/4)

# continous
map.london.o3 <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.pollution.london, aes(fill = o3), colour = NA) +
  ggplot2::scale_fill_gradient(name = '',
                               low = colours.pollutant[1], high = colours.pollutant[2],
                               limits = limits.o3, breaks = breaks.o3) +
  my.map.theme(legend.position = 'none')
map.o3 <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.pollution, aes(fill = o3), colour = NA) +
  ggplot2::scale_fill_gradient(name = unname(latex2exp ::TeX(c('$O_{3}(\\mu{g}/m^{3}$)'))),
                               guide = guide_colourbar(title.position = 'top'),
                               low = colours.pollutant[1], high = colours.pollutant[2],
                               limits = limits.o3, breaks = breaks.o3) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.o3), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  # ggplot2::labs(title = 'Pollution: 03') +
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5))

# discrete
map.london.o3.decile <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.pollution.london, aes(fill = o3_decile), colour = NA) +
  ggplot2::scale_fill_manual(name = '',
                             values = colours.pollutant.discrete) +
  my.map.theme(legend.position = 'none')
map.o3.decile <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.pollution, aes(fill = o3_decile), colour = NA) +
  ggplot2::scale_fill_manual(name = unname(latex2exp ::TeX(c('$O_{3}(\\mu{g}/m^{3}$) Decile'))),
                             guide = ggplot2::guide_legend(title.position = 'top'),
                             values = colours.pollutant.discrete) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.o3.decile), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  # ggplot2::labs(title = 'Pollution: 03') +
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5)); 

# save
setwd(dir.res.indicators)
ggplot2::ggsave(plot = map.o3,
                filename = 'map_o3.png',
                width = width, 
                height = height)
ggplot2::ggsave(plot = map.o3.decile,
                filename = 'map_o3_decile.png',
                width = width, 
                height = height)

# 99. collect all ----

data.indicators <-
  list(data.qof = data.qof %>% sf::st_drop_geometry() %>% dplyr::select(LSOA11CD, value_qof, qof, qof_id),
       data.epd = data.epd %>% sf::st_drop_geometry() %>% dplyr::select(LSOA11CD, value_epd, epd, epd_id),
       data.caci = data.caci %>% sf::st_drop_geometry() %>% dplyr::select(LSOA11CD, value_caci, caci, caci_id),
       data.ahah = data.ahah %>% sf::st_drop_geometry() %>% dplyr::select(LSOA11CD, value_ahah, ahah, ahah_id),
       data.imd = data.imd %>% sf::st_drop_geometry() %>% dplyr::select(LSOA11CD, value_imd, imd, imd_id),
       data.pollution = data.pollution %>% sf::st_drop_geometry() %>% dplyr::select(LSOA11CD, pm25, pm10, no2, o3))

setwd(dir.data.organised)
save(data.indicators, file = 'LSOA11_INDICATORS.rda')

