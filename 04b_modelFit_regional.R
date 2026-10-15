# 0.1. packages ----

if(!require('tidyverse')) {
  install.packages('tidyverse', dependencies = TRUE)
}

# 0.2. directories ----

# retrieve directories
code.path <- rstudioapi::getActiveDocumentContext()$path
dir.home <- sub('(NIH_Smoking_Proxies).*', '\\1', code.path)
dir.code <- paste0(dir.home, '/Code')
dir.data <- paste0(dir.home, '/Data')
dir.res <- paste0(dir.home, '/Results')
dir.data.organised <- paste0(dir.data, '/Organised Data')
dir.data.spatial <- paste0(dir.data, '/shapeFiles')
dir.res.model.fit <- paste0(dir.res, '/Model Fits')
dir.res.main.analysis <- paste0(dir.res, '/Main Analysis')

# 0.3. import ----

## 0.3.1. functions ----

setwd(dir.code)
source('00_functions.R')

## 0.3.2. data ----

setwd(dir.data.organised)
load('LSOA11_HES_FINAL.rda')

## 0.3.3. spatial information ----

setwd(dir.data)

lookup.lsoa.rgn <- readr::read_csv('ONS Code Look Ups/LSOA11_LAD11_RGN11_lookUp.csv')
poly.lsoa <- sf::st_read(dsn = dir.data.spatial, layer = 'ONS11_LSOA')
names.region <- 
  lookup.lsoa.rgn %>%
  dplyr::distinct(RGN11NM) %>%
  dplyr::filter(RGN11NM != 'Wales') %>% 
  dplyr::arrange() %>% 
  dplyr::pull()


## 0.4. output arguments ----

# plot 
height <- width <- 10
text.size <- 20

# colour schemes
colours.diverging <- c('#007D5E', '#E9E0C8', '#C13E2E')
colours.exceedance <- c('#A0D6D0', '#009B92', '#00768D')
colours.sequential <- c('#E9E0C8', '#A175A7')

# table
digits.table <- 4

# 1. regional adjacency matrices ----

setwd(paste0(dir.data, '/Organised Data'))
if(file.exists('ENGLAND_AMAT_LSOA_BY_REGION.rda')){
  load('ENGLAND_AMAT_LSOA_BY_REGION.rda')
} else {
  
  # function to turn poly into regional lsoa amat
  region.amat <- function(poly, lookup, region){
    
    # poly <- poly.lsoa
    # lookup <- lookup.lsoa.rgn
    # region <- 'South East'
    
    poly.filter <-
      poly %>%
      dplyr::filter(stringr::str_detect(LSOA11CD, '^E')) %>%
      dplyr::left_join(.,
                       lookup %>% dplyr::select(LSOA11CD, RGN11NM),
                       by = 'LSOA11CD') %>%
      dplyr::filter(RGN11NM == !!region)
    
    a.mat <- get.amat(poly.filter)
    rownames(a.mat) <- colnames(a.mat) <- paste0(gsub(' ', '', region), '_', 1:nrow(poly.filter))
    
    info <-
      data.frame(LSOA11CD = poly.filter$LSOA11CD,
                 LSOA11NM = poly.filter$LSOA11NM,
                 RGN11NM = poly.filter$RGN11NM) %>% 
      dplyr::mutate(space_id = 1:n(),
                    internal = paste0(gsub(' ', '', region), '_', 1:n()))
    
    return(list(a.mat = a.mat, info = info))
    
  }
  
  # region name as list 
  names.region.list <-
    names.region %>% 
    purrr::set_names() %>% 
    as.list()
  
  # run function
  lsoa.region.amat.information <- 
    lapply(X = names.region.list,
           FUN = region.amat,
           poly = poly.lsoa,
           lookup = lookup.lsoa.rgn)
  
  # save
  save(lsoa.region.amat.information, file = 'ENGLAND_AMAT_LSOA_BY_REGION.rda')
}

# 2. data ----

## 2.1. components of exposure, region, model ----

# exposures 
components.exposure <- c('pm25', 'pm10', 'no2', 'o3', 'o3plus')

# regions
components.region <- 
  lookup.lsoa.rgn %>%
  dplyr::distinct(RGN11NM) %>%
  dplyr::filter(RGN11NM != 'Wales') %>% 
  dplyr::arrange() %>% 
  dplyr::pull()

# models and formulas
components.model <-
  data.frame(model = 
               c(# model 0
                 'M0',
                 # model 1
                 'M1-QOF', 'M1-REACT', 'M1-LCRR', 'M1-EPD', 'M1-CACI', 'M1-AHAH', 'M1-IMD',
                 # model 2
                 'M2-QOF', 'M2-REACT', 'M2-LCRR', 'M2-EPD', 'M2-CACI', 'M2-AHAH',
                 # model 3
                 'M3'),
             formula.update = 
               c(# model 0
                 NA, 
                 # model 1
                 'qof_id', 'react_id', 'lcrr_id', 'epd_id', 'caci_id', 'ahah_id', 'imd_id',
                 # model 2
                 'imd_id + qof_id', 'imd_id + react_id', 'imd_id + lcrr_id', 'imd_id + epd_id', 'imd_id + caci_id', 'imd_id + ahah_id',
                 # model 3
                 'qof_id + react_id + lcrr_id + epd_id + caci_id + ahah_id + imd_id'))

parameter.names.data.frame <-
  data.frame(parameter.name.inla = c('(Intercept)', 'exposure',
                                     paste0('qof_id', 2:10), 
                                     paste0('react_id', 2:10),  
                                     paste0('lcrr_id', 2:10), 
                                     paste0('epd_id', 2:10), 
                                     paste0('caci_id', 2:10), 
                                     paste0('ahah_id', 2:10), 
                                     paste0('imd_id', 2:10), 
                                     'Precision for space_id', 'Phi for space_id'),
             parameter.name.label = c('Intercept', 'Exposure',
                                      rep('QOF', times = 9),  
                                      rep('REACT', times = 9),  
                                      rep('LCRR', times = 9), 
                                      rep('EPD', times = 9),
                                      rep('CACI', times = 9), 
                                      rep('AHAH', times = 9),
                                      rep('IMD', times = 9), 
                                      'Precision', 'Phi'))

## 2.2. complete lists ----

### 2.3.0. all ---- 

components.all <-
  expand.grid(exposure = components.exposure,
              region = components.region) %>% 
  dplyr::cross_join(., 
                    components.model) %>%
  dplyr::mutate(name = paste0(exposure, '_', region, '_', model) %>% gsub(pattern = ' ', replacement = '', x = .))

### 2.3.1. pm25 ----

components.list.pm25 <-
  components.all %>% 
  dplyr::filter(exposure == 'pm25') %>% 
  split(., f = .$name)

### 2.3.2. pm10 ----

components.list.pm10 <-
  components.all %>% 
  dplyr::filter(exposure == 'pm10') %>% 
  split(., f = .$name)

### 2.3.3. no2 ----

components.list.no2 <-
  components.all %>% 
  dplyr::filter(exposure == 'no2') %>% 
  split(., f = .$name)

### 2.3.4. o3 ----

components.list.o3 <-
  components.all %>% 
  dplyr::filter(exposure == 'o3') %>% 
  split(., f = .$name)

### 2.3.5. o3plus ----

components.list.o3plus <-
  expand.grid(exposure = components.exposure,
              region = components.region) %>% 
  dplyr::cross_join(., 
                    data.frame(model = 
                                 c(# model 0
                                   'M0',
                                   # model 1
                                   'M1-QOF', 'M1-REACT', 'M1-LCRR', 'M1-EPD', 'M1-CACI', 'M1-AHAH', 'M1-IMD',
                                   # model 2
                                   'M2-QOF', 'M2-REACT', 'M2-LCRR', 'M2-EPD', 'M2-CACI', 'M2-AHAH',
                                   # model 3
                                   'M3'),
                               formula.update = 
                                 c(# model 0
                                   'no2', 
                                   # model 1
                                   'no2 + qof_id', 'no2 + react_id', 'no2 + lcrr_id', 'no2 + epd_id', 'no2 + caci_id', 'no2 + ahah_id', 'no2 + imd_id',
                                   # model 2
                                   'no2 + imd_id + qof_id', 'no2 + imd_id + react_id', 'no2 + imd_id + lcrr_id', 'no2 + imd_id + epd_id', 'no2 + imd_id + caci_id', 'no2 + imd_id + ahah_id',
                                   # model 3
                                   'no2 + qof_id + react_id + lcrr_id + epd_id + caci_id + ahah_id + imd_id'))) %>%
  dplyr::mutate(name = paste0(exposure, 'plus_', region, '_', model) %>% gsub(pattern = ' ', replacement = '', x = .)) %>% 
  dplyr::filter(exposure == 'o3') %>% 
  split(., f = .$name)

# 3. run all models ----

## 3.0. function ----

regional.model.fit <- function(X, data, info.region.lsoa){
  
  # 0. function arguments ----
  
  # X <- components.list.pm25[[1]]
  # data <- data.final
  # info.region.lsoa <- lsoa.region.amat.information
  
  # 1. extract arguments ----
  
  selected.exposure <- X$exposure %>% as.character()
  selected.region <- X$region %>% as.character()
  selected.model <- X$model %>% as.character()
  selected.formula.update <- X$formula.update %>% as.character()
  
  # 2. spatial information ----
  
  info.all <- info.region.lsoa[[selected.region]]
  info.amat <- info.all$a.mat
  info.space <- info.all$info
  
  # 3. organise data ----
  
  data.2 <-
    data %>% 
    dplyr::select(-space_id) %>% 
    dplyr::left_join(., info.space, by = 'LSOA11CD') %>% 
    # only selected region will NOT be NA
    dplyr::filter(!is.na(RGN11NM)) %>% 
    # new column 'exposure' is the pollutant selected
    dplyr::mutate(exposure = .data[[selected.exposure]])
  
  # 4. fit model ----
  
  # do not assign this to anything
  # will return model.fit value
  model.fit(formula.update = selected.formula.update,
            data = data.2,
            a.mat = info.amat)
  
  # return(list(selected.exposure = X$exposure %>% as.character(),
  #             selected.region = X$region %>% as.character(),
  #             selected.model = X$model %>% as.character(),
  #             selected.formula.update = X$formula.update %>% as.character()))
  
  # selected.formula.update
  
  # selected.formula.update == 'NA'
  # is.na(selected.formula.update)
  
}

## 3.1. pm25 ----

setwd(paste0(dir.res.model.fit, '/Regional'))
if(file.exists('parameters_regional_pm25.csv')){
  
  # load
  parameters.regional.pm25 <- readr::read_csv(file = 'parameters_regional_pm25.csv')
  
} else {
  
  # fit all models
  fits.regional.pm25 <-
    lapply(X = components.list.pm25,
           FUN = regional.model.fit,
           data = data.final,
           info.region.lsoa = lsoa.region.amat.information)
  
  # parameters
  parameters.regional.pm25 <- 
    purrr::map(.x = names(fits.regional.pm25),
               .f = function(x) { fits.regional.pm25[[x]]$parameters %>% dplyr::mutate(name = x) } ) %>% 
    dplyr::bind_rows() %>%
    # better parameter labels
    dplyr::rename(parameter.name.inla = rowname) %>% 
    dplyr::left_join(., parameter.names.data.frame, by = 'parameter.name.inla') %>% 
    # join for speicific exposure/region/model names
    dplyr::left_join(., components.all %>% dplyr::select(-formula.update),
                     by = 'name') %>% 
    # nicer labels
    dplyr::mutate(exposure = exposure %>% factor(., levels = components.exposure),
                  region = region %>% factor(., levels = components.region %>% sort()),
                  model = model %>% factor(., levels = components.model$model %>% unique()),
                  model.type = model %>% gsub(pattern = '(-.*)', replacement = '', x = .))
  
  # save
  save(fits.regional.pm25,  file = 'fits_regional_pm25.rda')
  readr::write_csv(x = parameters.regional.pm25,
                   file = 'parameters_regional_pm25.csv')
  
  # clear
  rm(fits.regional.pm25,
     parameters.regional.pm25)
  gc()
  
}

## 3.2. pm10 ----

setwd(paste0(dir.res.model.fit, '/Regional'))
if(file.exists('parameters_regional_pm10.csv')){
  
  # load
  parameters.regional.pm10 <- readr::read_csv(file = 'parameters_regional_pm10.csv')
  
} else {
  
  # fit all models
  fits.regional.pm10 <-
    lapply(X = components.list.pm10,
           FUN = regional.model.fit,
           data = data.final,
           info.region.lsoa = lsoa.region.amat.information)
  
  # parameters
  parameters.regional.pm10 <- 
    purrr::map(.x = names(fits.regional.pm10),
               .f = function(x) { fits.regional.pm10[[x]]$parameters %>% dplyr::mutate(name = x) } ) %>% 
    dplyr::bind_rows() %>%
    # better parameter labels
    dplyr::rename(parameter.name.inla = rowname) %>% 
    dplyr::left_join(., parameter.names.data.frame, by = 'parameter.name.inla') %>% 
    # join for speicific exposure/region/model names
    dplyr::left_join(., components.all %>% dplyr::select(-formula.update),
                     by = 'name') %>% 
    # nicer labels
    dplyr::mutate(exposure = exposure %>% factor(., levels = components.exposure),
                  region = region %>% factor(., levels = components.region %>% sort()),
                  model = model %>% factor(., levels = components.model$model %>% unique()),
                  model.type = model %>% gsub(pattern = '(-.*)', replacement = '', x = .))
  
  # save
  save(fits.regional.pm10,  file = 'fits_regional_pm10.rda')
  readr::write_csv(x = parameters.regional.pm10,
                   file = 'parameters_regional_pm10.csv')
  
  # clear
  rm(fits.regional.pm10,
     parameters.regional.pm10)
  gc()
  
}

## 3.3. no2 ----

setwd(paste0(dir.res.model.fit, '/Regional'))
if(file.exists('parameters_regional_no2.csv')){
  
  # load
  parameters.regional.no2 <- readr::read_csv(file = 'parameters_regional_no2.csv')
  
} else {
  
  # fit all models
  fits.regional.no2 <-
    lapply(X = components.list.no2,
           FUN = regional.model.fit,
           data = data.final,
           info.region.lsoa = lsoa.region.amat.information)
  
  # parameters
  parameters.regional.no2 <- 
    purrr::map(.x = names(fits.regional.no2),
               .f = function(x) { fits.regional.no2[[x]]$parameters %>% dplyr::mutate(name = x) } ) %>% 
    dplyr::bind_rows() %>%
    # better parameter labels
    dplyr::rename(parameter.name.inla = rowname) %>% 
    dplyr::left_join(., parameter.names.data.frame, by = 'parameter.name.inla') %>% 
    # join for speicific exposure/region/model names
    dplyr::left_join(., components.all %>% dplyr::select(-formula.update),
                     by = 'name') %>% 
    # nicer labels
    dplyr::mutate(exposure = exposure %>% factor(., levels = components.exposure),
                  region = region %>% factor(., levels = components.region %>% sort()),
                  model = model %>% factor(., levels = components.model$model %>% unique()),
                  model.type = model %>% gsub(pattern = '(-.*)', replacement = '', x = .))
  
  # save
  save(fits.regional.no2,  file = 'fits_regional_no2.rda')
  readr::write_csv(x = parameters.regional.no2,
                   file = 'parameters_regional_no2.csv')
  
  # clear
  rm(fits.regional.no2,
     parameters.regional.no2)
  gc()
  
}

## 3.4. o3 ----

setwd(paste0(dir.res.model.fit, '/Regional'))
if(file.exists('parameters_regional_o3.csv')){
  
  # load
  parameters.regional.o3 <- readr::read_csv(file = 'parameters_regional_o3.csv')
  
} else {
  
  # fit all models
  fits.regional.o3 <-
    lapply(X = components.list.o3,
           FUN = regional.model.fit,
           data = data.final,
           info.region.lsoa = lsoa.region.amat.information)
  
  # parameters
  parameters.regional.o3 <- 
    purrr::map(.x = names(fits.regional.o3),
               .f = function(x) { fits.regional.o3[[x]]$parameters %>% dplyr::mutate(name = x) } ) %>% 
    dplyr::bind_rows() %>%
    # better parameter labels
    dplyr::rename(parameter.name.inla = rowname) %>% 
    dplyr::left_join(., parameter.names.data.frame, by = 'parameter.name.inla') %>% 
    # join for speicific exposure/region/model names
    dplyr::left_join(., components.all %>% dplyr::select(-formula.update),
                     by = 'name') %>% 
    # nicer labels
    dplyr::mutate(exposure = exposure %>% factor(., levels = components.exposure),
                  region = region %>% factor(., levels = components.region %>% sort()),
                  model = model %>% factor(., levels = components.model$model %>% unique()),
                  model.type = model %>% gsub(pattern = '(-.*)', replacement = '', x = .)) %>% 
    dplyr::mutate(parameter.name.label = dplyr::if_else(exposure == 'o3plus' & parameter.name.inla == 'no2', 'Adjustment', parameter.name.label))
  
  # save
  save(fits.regional.o3,  file = 'fits_regional_o3.rda')
  readr::write_csv(x = parameters.regional.o3,
                   file = 'parameters_regional_o3.csv')
  
  # clear
  rm(fits.regional.o3,
     parameters.regional.o3)
  gc()
  
}

## 3.5. o3plus ----

setwd(paste0(dir.res.model.fit, '/Regional'))
if(file.exists('parameters_regional_o3plus.csv')){
  
  # load
  parameters.regional.o3plus <- readr::read_csv(file = 'parameters_regional_o3plus.csv')
  
} else {
  
  # fit all models
  fits.regional.o3plus <-
    lapply(X = components.list.o3plus,
           FUN = regional.model.fit,
           data = data.final,
           info.region.lsoa = lsoa.region.amat.information)
  
  # parameters
  parameters.regional.o3plus <- 
    purrr::map(.x = names(fits.regional.o3plus),
               .f = function(x) { fits.regional.o3plus[[x]]$parameters %>% dplyr::mutate(name = x) } ) %>% 
    dplyr::bind_rows() %>%
    # better parameter labels
    dplyr::rename(parameter.name.inla = rowname) %>% 
    dplyr::left_join(., parameter.names.data.frame, by = 'parameter.name.inla') %>% 
    # join for speicific exposure/region/model names
    dplyr::left_join(., components.all %>% dplyr::select(-formula.update),
                     by = 'name') %>% 
    # nicer labels
    dplyr::mutate(exposure = exposure %>% factor(., levels = components.exposure),
                  region = region %>% factor(., levels = components.region %>% sort()),
                  model = model %>% factor(., levels = components.model$model %>% unique()),
                  model.type = model %>% gsub(pattern = '(-.*)', replacement = '', x = .))
  
  # save
  save(fits.regional.o3plus,  file = 'fits_regional_o3plus.rda')
  readr::write_csv(x = parameters.regional.o3plus,
                   file = 'parameters_regional_o3plus.csv')
  
  # clear
  rm(fits.regional.o3plus,
     parameters.regional.o3plus)
  gc()
  
}

## 3.6. combine all ----

parameters.regional.all <- 
  dplyr::bind_rows(parameters.regional.pm25,
                   parameters.regional.pm10,
                   parameters.regional.no2,
                   parameters.regional.o3,
                   parameters.regional.o3plus)

# 4. results ----

## 4.0. set up ----

## 4.0.1. functions ----

exposure.parameter.plot <- function(data, exposure, term = c('exposure', 'Phi for space_id')[1], text.size = 20){
  
  # 0. function arguments ----
  
  # data <- parameters.regional.all
  # exposure <- 'pm25'
  # term <- 'expsoure'
  
  # 1. data sorting ----
  
  ## 1.1. plot data ----
  
  data.2 <- 
    data %>% 
    dplyr::filter(parameter.name.inla == {{term}} & exposure == {{exposure}}) %>% 
    dplyr::arrange(region) %>% 
    dplyr::mutate(index = 1:n(),
                  .by = c('region'))
  
  ## 1.2. axis names ----
  
  axis.names <-
    data.2 %>% 
    dplyr::select(model, model.type, index) %>% 
    dplyr::distinct() %>% 
    dplyr::arrange(index)
  
  axis.lines <-
    axis.names %>% 
    dplyr::group_by(model.type) %>% 
    dplyr::slice_max(index) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(!(model.type == 'M3')) %>% 
    dplyr::mutate(x = index + 0.5,
                  xend = index + 0.5,
                  y = -Inf,
                  yend = Inf)
  
  # 2. plot ----
  
  plot <- 
    ggplot2::ggplot(data = data.2, aes(x = index, y = median)) +
    # point and error bar
    ggplot2::geom_point(aes(group = region), position = ggplot2::position_dodge(0.5)) +
    ggplot2::geom_errorbar(aes(ymin = lower, ymax = upper, group = region), width = 0, position = ggplot2::position_dodge(0.5)) +
    # labelling and axis sort
    ## axis labels
    ggplot2::labs(x = '', y = '') +
    ggplot2::scale_x_continuous(labels = axis.names$model,
                                breaks = 1:nrow(axis.names),
                                trans = 'reverse') +
    ## add divisors for models
    ggplot2::geom_segment(data = axis.lines, 
                          aes(x = x, xend = xend, y = y, yend = yend),
                          colour = 'grey') +
    ## flip axis
    ggplot2::coord_flip() +
    # facet
    ggplot2::facet_wrap(~ region, nrow = 1) +
    # theme
    my.theme(text = element_text(size = text.size),
             legend.title = ggplot2::element_blank(),
             legend.position = 'none')
  
  if(term == 'Phi for space_id'){
    plot 
  } else {
    plot +
    # line at y = 0
    ggplot2::geom_hline(yintercept = 0, colour = 'red3', linetype = 'dashed')
  }
  
}

exposure.parameter.heatmap <- function(data, exposure, term = c('exposure', 'Phi for space_id')[1], text.size = 20){
  
  # 0. function arguments ----
  
  # data <- parameters.all
  # exposure <- 'pm25'
  # term <- 'exposure'
  
  # 1. data sorting ----
  
  data.2 <- 
    data %>% 
    dplyr::filter(parameter.name.inla == {{term}} & exposure == {{exposure}}) 
  
  # 2. plot ----
  
  if(term == 'Phi for space_id') {
    midpoint <- 0.5
  } else {
    midpoint <- 0
  }
  
  ggplot2::ggplot(data = data.2, aes(x = region, y = model, fill = median)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(name = '',
                                  midpoint = midpoint,
                                  low = colours.diverging[1], 
                                  mid = colours.diverging[2],
                                  high = colours.diverging[3]) +
    ggplot2::scale_y_discrete(name = '',
                              limits = rev(levels(data.2$model))) +
    ggplot2::scale_x_discrete(name = '') +
    my.theme(text = ggplot2::element_text(size = text.size),
             legend.key.width = unit(2, 'cm'),
             legend.position = 'bottom',
             legend.title = ggplot2::element_text(hjust = 0.5),
             axis.text.x = ggplot2::element_text(angle = 90/5, vjust = 1, hjust = 1))
  
}

## 4.1. exposure ----

exposure.region.parameters.plot.pm25 <- exposure.parameter.plot(data = parameters.regional.all, exposure = 'pm25', term = 'exposure')
exposure.region.parameters.plot.pm10 <- exposure.parameter.plot(data = parameters.regional.all, exposure = 'pm10', term = 'exposure')
exposure.region.parameters.plot.no2 <- exposure.parameter.plot(data = parameters.regional.all, exposure = 'no2', term = 'exposure')
exposure.region.parameters.plot.o3 <- exposure.parameter.plot(data = parameters.regional.all, exposure = 'o3', term = 'exposure')

setwd(paste0(dir.res.main.analysis, '/Regional'))
ggplot2::ggsave(plot = exposure.region.parameters.plot.pm25,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'exposure.region.parameters.plot.pm25'),
                                  '.png'),
                width = width, 
                height = height)
ggplot2::ggsave(plot = exposure.region.parameters.plot.pm10,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'exposure.region.parameters.plot.pm10'),
                                  '.png'),
                width = width, 
                height = height)
ggplot2::ggsave(plot = exposure.region.parameters.plot.no2,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'exposure.region.parameters.plot.no2'),
                                  '.png'),
                width = width, 
                height = height)
ggplot2::ggsave(plot = exposure.region.parameters.plot.o3,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'exposure.region.parameters.plot.o3'),
                                  '.png'),
                width = width, 
                height = height)

## 4.2 spatial mixing ----

spatial.mixing.region.parameters.plot.pm25 <- exposure.parameter.plot(data = parameters.regional.all, exposure = 'pm25', term = 'Phi for space_id')
spatial.mixing.region.parameters.plot.pm10 <- exposure.parameter.plot(data = parameters.regional.all, exposure = 'pm10', term = 'Phi for space_id')
spatial.mixing.region.parameters.plot.no2 <- exposure.parameter.plot(data = parameters.regional.all, exposure = 'no2', term = 'Phi for space_id')
spatial.mixing.region.parameters.plot.o3 <- exposure.parameter.plot(data = parameters.regional.all, exposure = 'o3', term = 'Phi for space_id')

setwd(paste0(dir.res.main.analysis, '/Regional'))
ggplot2::ggsave(plot = spatial.mixing.region.parameters.plot.pm25,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'spatial.mixing.region.parameters.plot.pm25'),
                                  '.png'),
                width = width, 
                height = height)
ggplot2::ggsave(plot = spatial.mixing.region.parameters.plot.pm10,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'spatial.mixing.region.parameters.plot.pm10'),
                                  '.png'),
                width = width, 
                height = height)
ggplot2::ggsave(plot = spatial.mixing.region.parameters.plot.no2,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'spatial.mixing.region.parameters.plot.no2'),
                                  '.png'),
                width = width, 
                height = height)
ggplot2::ggsave(plot = spatial.mixing.region.parameters.plot.o3,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'spatial.mixing.region.parameters.plot.o3'),
                                  '.png'),
                width = width, 
                height = height)

## 4.3 all parameters table ----

table.parameters <-
  parameters.regional.all %>% 
  dplyr::select(exposure, region, model, parameter.name.label , estimate) %>% 
  dplyr::arrange(exposure, region)

setwd(paste0(dir.res.main.analysis, '/Regional'))
print(xtable::xtable(table.parameters),
      include.rownames = FALSE,
      file = 'model_parameters_region.txt')


