# packages ----

if(!require('tidyverse')) {
  install.packages('tidyverse', dependencies = TRUE)
}
if(!require('patchwork')) {
  install.packages('patchwork', dependencies = TRUE)
}

# directories ----

# retrieve directories
code.path <- rstudioapi::getActiveDocumentContext()$path
dir.home <- sub('(NIH_Smoking_Proxies).*', '\\1', code.path)
dir.code <- paste0(dir.home, '/Code')
dir.data <- paste0(dir.home, '/Data')
dir.res <- paste0(dir.home, '/Results')
dir.res.indicators <- paste0(dir.res, '/Indicators')
dir.data.spatial <- paste0(dir.data, '/shapeFiles')

# import ----

## functions ----

setwd(dir.code)
source('00_functions.R')

## spatial information ----

setwd(dir.data)
load('Organised Data/ENGLAND_SPATIAL_INFORMATION_LSOA_MSOA_LAD.rda')
poly.lsoa <- sf::st_read(dsn = dir.data.spatial, layer = 'ONS11_LSOA')
poly.lsoa.england <- poly.lsoa %>% dplyr::filter(stringr::str_detect(LSOA11CD, '^E'))
nbs.lsoa.england <- spdep::poly2nb(pl = as(poly.lsoa.england, 'Spatial'))
lookup.lsoa.lad <- read.csv('ONS Code Look Ups/LSOA11_LAD11_RGN11_lookUp.csv')

## data ----

setwd(dir.data)
data.import.extract <- read.csv('ConnorExtract.csv')
data.import.lcrr <- read.csv('LCRate2017NoCount.csv')
data.import.react <- read.csv('react_smoking_rates_export.csv')
load('Organised Data/LSOA11_POPULATION_TOTAL.rda')
load('Organised Data/LSOA11_INDICATORS.rda')

# plot arguments ----

# plot arguments
text.size <- 20
width <- height <- 10

# colour schemes
colours.diverging <- c('#007D5E', '#E9E0C8', '#C13E2E')
colours.exceedance <- c('#A0D6D0', '#009B92', '#00768D')
colours.covariate <- c('#E6F2F5', '#006D77')
colours.covariate.discrete <- grDevices::colorRampPalette(colours.covariate)(10)

# sort SAHSU indicators ----

## react ----

# Stage 1: get to lsoa level
data.sort.react.1 <- 
  data.import.react %>% 
  dplyr::rename(LSOA11CD = U_LSOA) %>% 
  dplyr::select(LSOA11CD, rate_current)

# Stage 2: complete lsoa structure
data.sort.react.2 <- 
  data.frame(LSOA11CD = lsoa.info$LSOA11CD) %>% 
  dplyr::left_join(., 
                   data.sort.react.1,
                   by = 'LSOA11CD')

# Stage 3: potential impute missing data
missing.data.react <-
  data.sort.react.2 %>% 
  dplyr::filter(rowSums(is.na(.)) > 0)

# impute missing values from neighbourhood structure
if(nrow(missing.data.react) > 0){
  data.sort.react.2$rate_current[is.na(data.sort.react.2$rate_current)] <-
    sapply(nbs.lsoa.england[which(is.na(data.sort.react.2$rate_current))],
           FUN = function(x) {mean(data.sort.react.2$rate_current[x], na.rm = TRUE)})
  
  setwd(dir.res.indicators)
  print(xtable::xtable(missing.data.react %>% 
                         dplyr::left_join(., lookup.lsoa.lad %>% dplyr::select(LSOA11CD, LSOA11NM), by = 'LSOA11CD') %>% 
                         dplyr::relocate(LSOA11NM, .after = LSOA11CD)),
        include.rownames = FALSE,
        file = 'missing_react.txt')
}

# Stage 4: complete score
data.react <-
  data.sort.react.2 %>% 
  dplyr::mutate(value_react = rate_current,
                react = cut(value_react, unique(quantile(value_react, seq(from = 0, to = 1, length.out = 11), na.rm = TRUE)), include.lowest = TRUE),
                react_id = react %>% factor(., labels = 1:10)) %>% 
  dplyr::left_join(., lookup.lsoa.lad %>% dplyr::select(LSOA11CD, LAD11CD), by = 'LSOA11CD') %>% 
  dplyr::left_join(.,poly.lsoa.england, by = 'LSOA11CD') %>% 
  sf::st_as_sf()

# stage 5: descriptive

setwd(dir.res.indicators)
print(xtable::xtable(data.react %>% 
                       dplyr::pull(value_react) %>%
                       {. * 100} %>% 
                       descriptive.function() %>% 
                       t() %>% 
                       as.data.frame() %>% 
                       dplyr::mutate(missing = nrow(missing.data.react)),
                     digits = 2),
      include.rownames = FALSE,
      file = 'descriptive_react.txt')

data.react.london <-
  data.react %>% 
  dplyr::filter(stringr::str_detect(LAD11CD, '^E09'))

limits.spatial <- 100*c(0, 0.5)
breaks.spatial <- seq(from = limits.spatial[1], to = limits.spatial[2], by = limits.spatial[2]/4)

# plot
## continuous
map.london.react <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.react.london, aes(fill = 100*value_react), colour = NA) +
  ggplot2::scale_fill_gradient(name = '',
                               low = colours.covariate[1], high = colours.covariate[2],
                               limits = limits.spatial, breaks = breaks.spatial) +
  my.map.theme(legend.position = 'none')
map.react <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.react, aes(fill = 100*value_react), colour = NA) +
  ggplot2::scale_fill_gradient(name = 'Proportion of smokers (%)',
                               guide = guide_colourbar(title.position = 'top'),
                               low = colours.covariate[1], high = colours.covariate[2],
                               limits = limits.spatial, breaks = breaks.spatial) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.react), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5))
## discrete
map.london.react.decile <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.react.london, aes(fill = react_id), colour = NA) +
  ggplot2::scale_fill_manual(name = '',
                             values = colours.covariate.discrete) +
  my.map.theme(legend.position = 'none')
map.react.decile <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.react, aes(fill = react_id), colour = NA) +
  ggplot2::scale_fill_manual(name = 'Proportion of smokers (%)',
                             guide = ggplot2::guide_legend(title.position = 'top'),
                             labels = c('1 (Least)', 2:9, '10 (Most)'),
                             values = colours.covariate.discrete) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.react.decile), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  # ggplot2::labs(title = 'Pollution: 03') +
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5)); (map.react/map.react.decile)
## historgram
breaks.react <-
  data.react %>% 
  dplyr::pull(value_react) %>%
  {. * 100} %>% 
  stats::quantile(x = ., probs = seq(from = 0, to = 1, length.out = 11), na.rm = TRUE) %>% 
  unique()
histogram.react <- 
  ggplot2::ggplot(data = data.react %>% sf::st_drop_geometry(), aes(x = 100*value_react)) +
  ggplot2::geom_histogram(colour = 'black', fill = 'white', bins = 250) +
  ggplot2::geom_vline(xintercept = breaks.react, colour = 'red3', linetype = 'dashed', linewidth = 0.75) +
  ggplot2::labs(x = 'Proportion of smokers (%)', y = 'Count') +
  my.theme(text = element_text(size = text.size))
# save
setwd(dir.res.indicators)
ggplot2::ggsave(plot = map.react,
                filename = 'map_react.png',
                width = width, 
                height = height)
ggplot2::ggsave(plot = map.react.decile,
                filename = 'map_react_decile.png',
                width = width, 
                height = height)
ggplot2::ggsave(plot = histogram.react,
                filename = 'histogram_react.png',
                width = width, 
                height = height)

## lcrr ----

# Stage 1: get to lsoa level
data.sort.lcrr.1 <- 
  data.import.lcrr %>% 
  dplyr::rename(LSOA11CD = LSOA11) %>% 
  dplyr::select(LSOA11CD, SmoothedRR)

# Stage 2: complete lsoa structure
data.sort.lcrr.2 <- 
  data.frame(LSOA11CD = lsoa.info$LSOA11CD) %>% 
  dplyr::left_join(., 
                   data.sort.lcrr.1,
                   by = 'LSOA11CD') %>% 
  ### Manually replace Sutton 021A which has smoothedRR as 47.524978 ###
  dplyr::mutate(SmoothedRR = dplyr::if_else(SmoothedRR > 47, NA, SmoothedRR))

# Stage 3: potential impute missing data
missing.data.lcrr <-
  data.sort.lcrr.2 %>% 
  dplyr::filter(rowSums(is.na(.)) > 0)

# impute missing values from neighbourhood structure
if(nrow(missing.data.lcrr) > 0){
  data.sort.lcrr.2$SmoothedRR[is.na(data.sort.lcrr.2$SmoothedRR)] <-
    sapply(nbs.lsoa.england[which(is.na(data.sort.lcrr.2$SmoothedRR))],
           FUN = function(x) {mean(data.sort.lcrr.2$SmoothedRR[x], na.rm = TRUE)})
  
  setwd(dir.res.indicators)
  print(xtable::xtable(missing.data.lcrr %>% 
                         dplyr::left_join(., lookup.lsoa.lad %>% dplyr::select(LSOA11CD, LSOA11NM), by = 'LSOA11CD') %>% 
                         dplyr::relocate(LSOA11NM, .after = LSOA11CD)),
        include.rownames = FALSE,
        file = 'missing_lcrr.txt')
  
}

# Stage 4: complete score
data.lcrr <-
  data.sort.lcrr.2 %>% 
  dplyr::mutate(value_lcrr = SmoothedRR,
                lcrr = cut(value_lcrr, unique(quantile(value_lcrr, seq(from = 0, to = 1, length.out = 11), na.rm = TRUE)), include.lowest = TRUE),
                lcrr_id = lcrr %>% factor(., labels = 1:10)) %>% 
  dplyr::left_join(., lookup.lsoa.lad %>% dplyr::select(LSOA11CD, LAD11CD), by = 'LSOA11CD') %>% 
  dplyr::left_join(.,poly.lsoa.england, by = 'LSOA11CD') %>% 
  sf::st_as_sf()

# stage 5: descriptive

setwd(dir.res.indicators)
print(xtable::xtable(data.lcrr %>% 
                       dplyr::pull(value_lcrr) %>% 
                       descriptive.function %>% 
                       t() %>% 
                       as.data.frame() %>% 
                       dplyr::mutate(missing = nrow(missing.data.lcrr)),
                     digits = 2),
      include.rownames = FALSE,
      file = 'descriptive_lcrr.txt')

data.lcrr.london <-
  data.lcrr %>% 
  dplyr::filter(stringr::str_detect(LAD11CD, '^E09'))

# Smoothed RR for E01004090 is 47.524978. Next highest is E01006558 at 4.350886
limits.spatial <- c(0, 5)
breaks.spatial <- seq(from = limits.spatial[1], to = limits.spatial[2], by = limits.spatial[2]/4)

# plot
## continuous
map.london.lcrr <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.lcrr.london, aes(fill = value_lcrr), colour = NA) +
  ggplot2::scale_fill_gradient(name = '',
                               low = colours.covariate[1], high = colours.covariate[2],
                               limits = limits.spatial, breaks = breaks.spatial) +
  my.map.theme(legend.position = 'none')
map.lcrr <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.lcrr, aes(fill = value_lcrr), colour = NA) +
  ggplot2::scale_fill_gradient(name = 'Relative risk to national average',
                               guide = guide_colourbar(title.position = 'top'),
                               low = colours.covariate[1], high = colours.covariate[2],
                               limits = limits.spatial, breaks = breaks.spatial) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.lcrr), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5))
## discrete
map.london.lcrr.decile <-
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.lcrr.london, aes(fill = lcrr_id), colour = NA) +
  ggplot2::scale_fill_manual(name = '',
                             values = colours.covariate.discrete) +
  my.map.theme(legend.position = 'none')
map.lcrr.decile <- 
  ggplot2::ggplot() +
  ggplot2::geom_sf(data = data.lcrr, aes(fill = lcrr_id), colour = NA) +
  ggplot2::scale_fill_manual(name = 'Relative risk to national average',
                             guide = ggplot2::guide_legend(title.position = 'top'),
                             labels = c('1 (Lowest)', 2:9, '10 (Highest)'),
                             values = colours.covariate.discrete) +
  # add london in seperately
  ggplot2::geom_rect(aes(xmin = 490000, xmax = 560000, ymin = 150000, ymax = 210000), color = 'red3', linewidth = 1, fill = NA) +
  ggplot2::annotation_custom(ggplot2::ggplotGrob(map.london.lcrr.decile), xmin = 82700, xmax = 300000, ymin = 160000, ymax = 450000) +
  # title 
  # ggplot2::labs(title = 'Pollution: 03') +
  my.map.theme(plot.title = element_text(hjust = 0.5),
               text = element_text(size = text.size),
               legend.position = 'bottom',
               legend.key.width = unit(2, 'cm'),
               legend.title = element_text(hjust = 0.5)); (map.lcrr/map.lcrr.decile)
## historgram
breaks.lcrr <-
  data.lcrr %>% 
  dplyr::pull(value_lcrr) %>% 
  quantile(x = ., probs = seq(from = 0, to = 1, length.out = 11), na.rm = TRUE) %>% 
  unique()
histogram.lcrr <-
  ggplot2::ggplot(data = data.lcrr %>% sf::st_drop_geometry(), aes(x = value_lcrr)) +
  ggplot2::geom_histogram(colour = 'black', fill = 'white', bins = 250) +
  ggplot2::geom_vline(xintercept = breaks.lcrr, colour = 'red3', linetype = 'dashed', linewidth = 0.75) +
  ggplot2::labs(x = 'Relative risk to national average', y = 'Count') +
  my.theme(text = element_text(size = text.size))

# save
setwd(dir.res.indicators)
ggplot2::ggsave(plot = map.lcrr,
                filename = 'map_lcrr.png',
                width = width, 
                height = height)
ggplot2::ggsave(plot = map.lcrr.decile,
                filename = 'map_lcrr_decile.png',
                width = width, 
                height = height)
ggplot2::ggsave(plot = histogram.lcrr,
                filename = 'histogram_lcrr.png',
                width = width, 
                height = height)

# finalise data ----

## admission data ----

data.admissions <-
  data.import.extract %>% 
  # reformat dates
  dplyr::mutate(ADMIDATE = ADMIDATE %>% as.Date(., format = '%d/%m/%Y'),
                DISDATE = DISDATE %>% as.Date(., format = '%d/%m/%Y')) %>% 
  dplyr::mutate(YEAR = ADMIDATE %>% lubridate::year(),
                AGE = ADMIAGE,
                ADMISSIONS = 1) %>% 
  # remove missing 
  dplyr::filter(!is.na(AGE)) %>%
  # new age label 
  dplyr::mutate(AGE.NEW = dplyr::if_else(AGE < 90, as.character(AGE), '90+')) %>% 
  # summarise by LSOA, year, age and sex
  dplyr::summarise(ADMISSIONS = sum(ADMISSIONS),
                   .by = c('LSOA11CD', 'YEAR', 'AGE.NEW', 'SEX')) %>% 
  # rename age
  dplyr::rename(AGE = AGE.NEW) %>% 
  # age grouped into 5 year groups
  dplyr::mutate(AGE.GROUP = dplyr::case_when(AGE %in% c(0:4) ~ '0-4',
                                             AGE %in% c(5:9) ~ '5-9',
                                             AGE %in% c(10:14) ~ '10-14',
                                             AGE %in% c(15:19) ~ '15-19',
                                             AGE %in% c(20:24) ~ '20-24',
                                             AGE %in% c(25:29) ~ '25-29',
                                             AGE %in% c(30:34) ~ '30-34',
                                             AGE %in% c(35:39) ~ '35-39',
                                             AGE %in% c(40:44) ~ '40-44',
                                             AGE %in% c(45:49) ~ '45-49',
                                             AGE %in% c(50:54) ~ '50-54',
                                             AGE %in% c(55:59) ~ '55-59',
                                             AGE %in% c(60:64) ~ '60-64',
                                             AGE %in% c(65:69) ~ '65-69',
                                             AGE %in% c(70:74) ~ '70-74',
                                             AGE %in% c(75:79) ~ '75-79',
                                             AGE %in% c(80:84) ~ '80-84',
                                             AGE %in% c(85:89) ~ '85-89',
                                             .default = '90+'),
                SEX = dplyr::if_else(SEX == 1, 'MALES', 'FEMALES')) %>% 
  dplyr::summarise(Y = sum(ADMISSIONS),
                   .by = c('LSOA11CD', 'YEAR', 'AGE.GROUP', 'SEX')) %>% 
  dplyr::arrange(LSOA11CD, YEAR, AGE.GROUP, SEX)

## population data ----

data.population <- 
  population.lsoa %>% 
  # age grouped into 5 year groups
  dplyr::mutate(AGE.GROUP = dplyr::case_when(AGE %in% c(0:4) ~ '0-4',
                                             AGE %in% c(5:9) ~ '5-9',
                                             AGE %in% c(10:14) ~ '10-14',
                                             AGE %in% c(15:19) ~ '15-19',
                                             AGE %in% c(20:24) ~ '20-24',
                                             AGE %in% c(25:29) ~ '25-29',
                                             AGE %in% c(30:34) ~ '30-34',
                                             AGE %in% c(35:39) ~ '35-39',
                                             AGE %in% c(40:44) ~ '40-44',
                                             AGE %in% c(45:49) ~ '45-49',
                                             AGE %in% c(50:54) ~ '50-54',
                                             AGE %in% c(55:59) ~ '55-59',
                                             AGE %in% c(60:64) ~ '60-64',
                                             AGE %in% c(65:69) ~ '65-69',
                                             AGE %in% c(70:74) ~ '70-74',
                                             AGE %in% c(75:79) ~ '75-79',
                                             AGE %in% c(80:84) ~ '80-84',
                                             AGE %in% c(85:89) ~ '85-89',
                                             .default = '90+')) %>% 
  dplyr::summarise(N = sum(POPULATION),
                   .by = c('LSOA11CD', 'YEAR', 'AGE.GROUP', 'SEX')) %>% 
  dplyr::arrange(LSOA11CD, YEAR, AGE.GROUP, SEX)

## indirect standardisation ----

### balanced data ----

age.label <- c('0-4', '5-9', '10-14', '15-19', '20-24', '25-29', '30-34', '35-39',
               '40-44', '45-49', '50-54', '55-59', '60-64', '65-69', '70-74',
               '75-79', '80-84', '85-89', '90+')
sex.label <- c('MALES', 'FEMALES')

data.balanced <-
  expand.grid(LSOA11CD = lsoa.info$LSOA11CD,
              YEAR = 2017:2019,
              AGE.GROUP = age.label,
              SEX = sex.label)

### reference data ----

# want AGE-SEX specific rates per LSOA
data.reference <-
  data.balanced %>%
  # join HES data to balanced data
  dplyr::left_join(., data.admissions, 
                   by = c('LSOA11CD', 'YEAR', 'AGE.GROUP', 'SEX')) %>% 
  ## replace na's with 0
  dplyr::mutate(Y = Y %>% tidyr::replace_na(., replace = 0)) %>% 
  # join (yearly) population totals
  dplyr::left_join(., data.population,
                   by = c('LSOA11CD', 'YEAR', 'AGE.GROUP', 'SEX')) %>% 
  # sum per LSOA-AGE-SEX
  dplyr::summarise(Y = sum(Y),
                   N = sum(N),
                   .by = c('AGE.GROUP', 'SEX')) %>% 
  # AGE-SEX specific rates
  dplyr::mutate(R = Y/N)

### final data ----

data.final <- 
  data.balanced %>%
  # join HES data to balanced data
  ## replace na's with 0
  dplyr::left_join(., data.admissions, 
                   by = c('LSOA11CD', 'YEAR', 'AGE.GROUP', 'SEX')) %>% 
  dplyr::mutate(Y = Y %>% tidyr::replace_na(., replace = 0)) %>%
  # join population total
  dplyr::left_join(., data.population,
                   by = c('LSOA11CD', 'YEAR', 'AGE.GROUP', 'SEX')) %>% 
  # join age-sex rates
  dplyr::left_join(., data.reference %>% dplyr::select(AGE.GROUP, SEX, R),
                   by = c('AGE.GROUP', 'SEX')) %>% 
  # expected = population * age-sex specific rate
  dplyr::mutate(E = N * R) %>% 
  # summarise over age-sex as adjusted for
  dplyr::summarise(Y = sum(Y),
                   N = sum(N),
                   E = sum(E),
                   .by = c('LSOA11CD')) %>% 
  # join indicators
  ## qof
  dplyr::left_join(.,
                   data.indicators[['data.qof']] %>% 
                     # select only needed
                     dplyr::select(LSOA11CD, qof, qof_id),
                   by = 'LSOA11CD') %>%
  ## react
  dplyr::left_join(.,
                   data.react %>% 
                     dplyr::select(LSOA11CD, react, react_id),
                   by = 'LSOA11CD') %>% 
  ## lcrr
  dplyr::left_join(.,
                   data.lcrr %>%  
                     sf::st_drop_geometry() %>% 
                     dplyr::select(LSOA11CD, lcrr, lcrr_id),
                   by = 'LSOA11CD') %>%
  ## epd
  dplyr::left_join(.,
                   data.indicators[['data.epd']] %>% 
                     dplyr::select(LSOA11CD, epd, epd_id),
                   by = 'LSOA11CD') %>%
  ## caci
  dplyr::left_join(.,
                   data.indicators[['data.caci']] %>%
                     dplyr::select(LSOA11CD, caci, caci_id),
                   by = 'LSOA11CD') %>%
  ## ahah
  dplyr::left_join(.,
                   data.indicators[['data.ahah']] %>%  
                     dplyr::select(LSOA11CD, ahah, ahah_id),
                   by = 'LSOA11CD') %>%
  ## imd 
  dplyr::left_join(.,
                   data.indicators[['data.imd']] %>% 
                     dplyr::select(LSOA11CD, imd, imd_id),
                   by = 'LSOA11CD') %>%
  ## pollution
  dplyr::left_join(.,
                   data.indicators[['data.pollution']] %>% 
                     dplyr::select(LSOA11CD, pm25, pm10, no2, o3) %>% 
                     sf::st_drop_geometry(),
                   by = 'LSOA11CD') %>% 
  # labels
  ## correct spatial labels
  dplyr::left_join(., 
                   lsoa.info %>% dplyr::mutate(space_id = 1:n()) %>% dplyr::select(LSOA11CD, space_id),
                   by = 'LSOA11CD') %>% 
  dplyr::relocate(c('qof_id', 'react_id', 'lcrr_id', 'caci_id', 'epd_id', 'ahah_id', 'imd_id'), .after = 'space_id')

### check ----

data.final$Y %>% sum() == data.final$E %>% sum()
data.final %>% is.na() %>% sum()

### save ----

setwd(dir.data)
save(data.final, file = 'Organised Data/LSOA11_HES_FINAL.rda')
