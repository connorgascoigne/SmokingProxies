# 0.1. packages ----

if(!require('tidyverse')) {
  install.packages('tidyverse', dependencies = TRUE)
}
if(!require('INLA')) {
  install.packages('INLA', dependencies = TRUE)
}

# 0.2. directories ----

# retrieve directories
code.path <- rstudioapi::getActiveDocumentContext()$path
dir.home <- sub('(NIH_Smoking_Proxies).*', '\\1', code.path)
dir.code <- paste0(dir.home, '/Code')
dir.data <- paste0(dir.home, '/Data')
dir.res <- paste0(dir.home, '/Results')
dir.data.organised <- paste0(dir.data, '/Organised Data')
dir.res.model.fit <- paste0(dir.res, '/Model Fits')
dir.res.main.analysis <- paste0(dir.res, '/Main Analysis')

# 0.3. imports ----

## 0.3.1. functions ----

setwd(dir.code)
source('00_functions.R')

## 0.3.2. organised data ----

setwd(dir.data.organised)
load('LSOA11_HES_FINAL.rda')
load('ENGLAND_SPATIAL_INFORMATION_LSOA_MSOA_LAD.rda')
load('ENGLAND_AMAT_LSOA_MSOA_LAD.rda')
rm(msoa.info, msoa.mat, lad.info, lad.mat) # remove unneeded variables (space)

# 0.4. output arguments ----

# plot 
height <- width <- 10
text.size <- 20

# table
digits.table <- 4

# 1. data sort ----

## 1.1. components of exposure, model ----

# exposures 
components.exposure <- c('pm25', 'pm10', 'no2', 'o3', 'o3plus')

# models and formulas
components.model <-
  data.frame(model = 
               c(# model 0
                 'M0',
                 # model 1
                 'M1-QOF', 'M1-REACT', 'M1-LCRR', 'M1-EPD', 'M1-CACI', 'M1-AHAH', 'M1-IMD',
                 # model 2
                 'M2-QOF', 'M2-REACT', 'M2-LCRR', 'M2-EPD', 'M2-CACI', 'M2-AHAH',
                 # # model additional
                 # 'M2a-REACT', 'M2a-LCRR', 'M2a-EPD','M2a-CACI', 'M2a-AHAH',
                 # model 3
                 'M3'),
             formula.update = 
               c(# model 0
                 NA, 
                 # model 1
                 'qof_id', 'react_id', 'lcrr_id', 'epd_id', 'caci_id', 'ahah_id', 'imd_id',
                 # model 2
                 'imd_id + qof_id', 'imd_id + react_id', 'imd_id + lcrr_id', 'imd_id + epd_id', 'imd_id + caci_id', 'imd_id + ahah_id',
                 # # model additional
                 # 'qof_id + react_id', 'qof_id + lcrr_id', 'qof_id + epd_id', 'qof_id + caci_id', 'qof_id + ahah_id',
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

## 1.2. complete lists ----

### 1.2.0. all ----

components.all <-
  data.frame(exposure = components.exposure) %>% 
  dplyr::cross_join(., 
                    components.model) %>%
  dplyr::mutate(name = paste0(exposure, '_', model) %>% gsub(pattern = ' ', replacement = '', x = .))

### 1.2.1. pm25 ----

components.list.pm25 <-
  components.all %>% 
  dplyr::filter(exposure == 'pm25') %>% 
  split(., f = .$name)

### 1.2.2. pm10 ----

components.list.pm10 <-
  components.all %>% 
  dplyr::filter(exposure == 'pm10') %>% 
  split(., f = .$name)

### 1.2.3. no2 ----

components.list.no2 <-
  components.all %>% 
  dplyr::filter(exposure == 'no2') %>% 
  split(., f = .$name)

### 1.2.4. o3 ----

components.list.o3 <-
  components.all %>% 
  dplyr::filter(exposure == 'o3') %>% 
  split(., f = .$name)


### 1.2.5. o3plus ----

components.list.o3plus <-
  data.frame(exposure = 'o3') %>% 
  dplyr::cross_join(., 
                    data.frame(model = 
                                 c(# model 0
                                   'M0',
                                   # model 1
                                   'M1-QOF', 'M1-REACT', 'M1-LCRR', 'M1-EPD', 'M1-CACI', 'M1-AHAH', 'M1-IMD',
                                   # model 2
                                   'M2-QOF', 'M2-REACT', 'M2-LCRR', 'M2-EPD', 'M2-CACI', 'M2-AHAH',
                                   # # model additional
                                   # 'M2a-REACT', 'M2a-LCRR', 'M2a-EPD','M2a-CACI', 'M2a-AHAH',
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
  dplyr::mutate(name = paste0(exposure, 'plus_', model) %>% gsub(pattern = ' ', replacement = '', x = .)) %>%
  dplyr::filter(exposure == 'o3') %>%  
  split(., f = .$name)

# 2. fitting models ----

## 2.0. function ----

national.model.fit <- function(X, data, a.mat){
  
  # 0. function arguments ----
  
  # X <- components.all[1,]
  # data <- data.final.2
  
  # 1. extract arguments ----
  
  selected.exposure <- X$exposure %>% as.character()
  selected.model <- X$model %>% as.character()
  selected.formula.update <- X$formula.update %>% as.character()
  
  # 2. organise data ----
  
  data.2 <-
    data %>% 
    # new column 'exposure' is the pollutant selected
    dplyr::mutate(exposure = .data[[selected.exposure]])
  
  # 4. fit model ----
  
  # do not assign this to anything
  # will return model.fit value
  model.fit(formula.update = selected.formula.update,
            data = data.2,
            a.mat = a.mat)
  
}

## 2.1. pm25 ----

setwd(paste0(dir.res.model.fit, '/National'))
if(file.exists('parameters_national_pm25.csv')){
  
  # load
  parameters.pm25 <- readr::read_csv(file = 'parameters_national_pm25.csv')
  
} else {
  
  # fit
  fits.pm25 <-
    lapply(X = components.list.pm25,
           FUN = national.model.fit,
           data = data.final,
           a.mat = lsoa.mat)
  
  # parameters
  parameters.pm25 <- 
    purrr::map(.x = names(fits.pm25),
               .f = function(x) { fits.pm25[[x]]$parameters %>% dplyr::mutate(name = x) } ) %>% 
    dplyr::bind_rows() %>%
    # better parameter labels
    dplyr::rename(parameter.name.inla = rowname) %>% 
    dplyr::left_join(., parameter.names.data.frame, by = 'parameter.name.inla') %>% 
    # join for specific exposure/region/model names
    dplyr::left_join(., components.all %>% dplyr::select(-formula.update),
                     by = 'name') %>% 
    # nicer labels
    dplyr::mutate(exposure = exposure %>% factor(., levels = components.exposure),
                  region = 'England',
                  model = model %>% factor(., levels = components.model$model %>% unique()),
                  model.type = model %>% gsub(pattern = '(-.*)', replacement = '', x = .))
  
  # save 
  save(fits.pm25, file = 'fits_national_pm25.rda')
  readr::write_csv(x = parameters.pm25,
                   file = 'parameters_national_pm25.csv')
  
  # clear
  rm(fits.pm25, 
     parameters.pm25)
  gc()
  
  # # load
  # readr::read_csv(file = 'parameters_national_pm25.csv')
  
}

## 2.2. pm10 ----

setwd(paste0(dir.res.model.fit, '/National'))
if(file.exists('parameters_national_pm10.csv')){
  
  # load
  parameters.pm10 <- readr::read_csv(file = 'parameters_national_pm10.csv')
  
} else {
  
  # fit
  fits.pm10 <-
    lapply(X = components.list.pm10,
           FUN = national.model.fit,
           data = data.final,
           a.mat = lsoa.mat)
  
  # parameters
  parameters.pm10 <- 
    purrr::map(.x = names(fits.pm10),
               .f = function(x) { fits.pm10[[x]]$parameters %>% dplyr::mutate(name = x) } ) %>% 
    dplyr::bind_rows() %>%
    # better parameter labels
    dplyr::rename(parameter.name.inla = rowname) %>% 
    dplyr::left_join(., parameter.names.data.frame, by = 'parameter.name.inla') %>% 
    # join for specific exposure/region/model names
    dplyr::left_join(., components.all %>% dplyr::select(-formula.update),
                     by = 'name') %>% 
    # nicer labels
    dplyr::mutate(exposure = exposure %>% factor(., levels = components.exposure),
                  region = 'England',
                  model = model %>% factor(., levels = components.model$model %>% unique()),
                  model.type = model %>% gsub(pattern = '(-.*)', replacement = '', x = .))
  
  # save 
  save(fits.pm10, file = 'fits_national_pm10.rda')
  readr::write_csv(x = parameters.pm10,
                   file = 'parameters_national_pm10.csv')
  
  # clear
  rm(fits.pm10,
     parameters.pm10)
  gc()
  
}

## 2.3. no2 ----

setwd(paste0(dir.res.model.fit, '/National'))
if(file.exists('parameters_national_no2.csv')){
  
  # load
  parameters.no2 <- readr::read_csv(file = 'parameters_national_no2.csv')
  
} else {
  
  # fit
  fits.no2 <-
    lapply(X = components.list.no2,
           FUN = national.model.fit,
           data = data.final,
           a.mat = lsoa.mat)
  
  # parameters
  parameters.no2 <- 
    purrr::map(.x = names(fits.no2),
               .f = function(x) { fits.no2[[x]]$parameters %>% dplyr::mutate(name = x) } ) %>% 
    dplyr::bind_rows() %>%
    # better parameter labels
    dplyr::rename(parameter.name.inla = rowname) %>% 
    dplyr::left_join(., parameter.names.data.frame, by = 'parameter.name.inla') %>% 
    # join for specific exposure/region/model names
    dplyr::left_join(., components.all %>% dplyr::select(-formula.update),
                     by = 'name') %>% 
    # nicer labels
    dplyr::mutate(exposure = exposure %>% factor(., levels = components.exposure),
                  region = 'England',
                  model = model %>% factor(., levels = components.model$model %>% unique()),
                  model.type = model %>% gsub(pattern = '(-.*)', replacement = '', x = .))
  
  # save 
  save(fits.no2, file = 'fits_national_no2.rda')
  readr::write_csv(x = parameters.no2,
                   file = 'parameters_national_no2.csv')
  
  # clear
  rm(fits.no2, 
     parameters.no2)
  gc()
  
}

## 2.4. o3 ----

setwd(paste0(dir.res.model.fit, '/National'))
if(file.exists('parameters_national_o3.csv')){
  
  # load
  parameters.o3 <- readr::read_csv(file = 'parameters_national_o3.csv')
  
} else {
  
  # fit
  fits.o3 <-
    lapply(X = components.list.o3,
           FUN = national.model.fit,
           data = data.final,
           a.mat = lsoa.mat)
  
  # parameters
  parameters.o3 <- 
    purrr::map(.x = names(fits.o3),
               .f = function(x) { fits.o3[[x]]$parameters %>% dplyr::mutate(name = x) } ) %>% 
    dplyr::bind_rows() %>%
    # better parameter labels
    dplyr::rename(parameter.name.inla = rowname) %>% 
    dplyr::left_join(., parameter.names.data.frame, by = 'parameter.name.inla') %>% 
    # join for specific exposure/region/model names
    dplyr::left_join(., components.all %>% dplyr::select(-formula.update),
                     by = 'name') %>% 
    # nicer labels
    dplyr::mutate(exposure = exposure %>% factor(., levels = components.exposure),
                  region = 'England',
                  model = model %>% factor(., levels = components.model$model %>% unique()),
                  model.type = model %>% gsub(pattern = '(-.*)', replacement = '', x = .))
  
  # save 
  save(fits.o3, file = 'fits_national_o3.rda')
  readr::write_csv(x = parameters.o3,
                   file = 'parameters_national_o3.csv')
  
  # clear
  rm(fits.o3, 
     parameters.o3)
  gc()
  
}

## 2.5. o3plus ----

setwd(paste0(dir.res.model.fit, '/National'))
if(file.exists('parameters_national_o3plus.csv')){
  
  # load
  parameters.o3plus <- readr::read_csv(file = 'parameters_national_o3plus.csv')
  
} else {
  
  # fit
  fits.o3plus <-
    lapply(X = components.list.o3plus,
           FUN = national.model.fit,
           data = data.final,
           a.mat = lsoa.mat)
  
  # parameters
  parameters.o3plus <- 
    purrr::map(.x = names(fits.o3plus),
               .f = function(x) { fits.o3plus[[x]]$parameters %>% dplyr::mutate(name = x) } ) %>% 
    dplyr::bind_rows() %>%
    # better parameter labels
    dplyr::rename(parameter.name.inla = rowname) %>% 
    dplyr::left_join(., parameter.names.data.frame, by = 'parameter.name.inla') %>% 
    # join for specific exposure/region/model names
    dplyr::left_join(., components.all %>% dplyr::select(-formula.update),
                     by = 'name') %>% 
    # nicer labels
    dplyr::mutate(exposure = exposure %>% factor(., levels = components.exposure),
                  region = 'England',
                  model = model %>% factor(., levels = components.model$model %>% unique()),
                  model.type = model %>% gsub(pattern = '(-.*)', replacement = '', x = .)) %>% 
    dplyr::mutate(parameter.name.label = dplyr::if_else(exposure == 'o3plus' & parameter.name.inla == 'no2', 'Adjustment', parameter.name.label))
  
  # save 
  save(fits.o3plus, file = 'fits_national_o3plus.rda')
  readr::write_csv(x = parameters.o3plus,
                   file = 'parameters_national_o3plus.csv')
  
  # clear
  rm(fits.o3plus, 
     parameters.o3plus)
  gc()
  
}

## 2.6. combine all ----

parameters.all <- 
  dplyr::bind_rows(parameters.pm25,
                   parameters.pm10,
                   parameters.no2,
                   parameters.o3,
                   parameters.o3plus)

# 3. results ----


## 3.0. set up ----

### 3.0.1. functions ----

exposure.parameter.plot <- function(data, exposure, term = c('exposure', 'Phi for space_id')[1], text.size = 20){
  
  # 0. function arguments ----
  
  # data <- parameters.all
  # exposure <- 'pm25'
  # term <- 'exposure'
  
  # 1. data sorting ----
  
  ## 1.1. plot data ----
  
  data.2 <- 
    data %>% 
    dplyr::filter(parameter.name.inla == {{term}} & exposure == {{exposure}}) %>% 
    dplyr::arrange(region) %>% 
    dplyr::mutate(index = 1:n(),
                  .by = 'exposure')
  
  ## 1.2. axis names ----
  
  axis.names <-
    data.2 %>% 
    dplyr::select(region, model, model.type, index) %>% 
    dplyr::distinct() %>% 
    dplyr::arrange(region, index)
  
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
    ggplot2::geom_point(position = ggplot2::position_dodge(0.5)) +
    ggplot2::geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, position = ggplot2::position_dodge(0.5)) +
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

## 3.1. exposure ----

exposure.nation.parameters.plot.pm25 <- exposure.parameter.plot(data = parameters.all, exposure = 'pm25')
exposure.nation.parameters.plot.pm10 <- exposure.parameter.plot(data = parameters.all, exposure = 'pm10')
exposure.nation.parameters.plot.no2 <- exposure.parameter.plot(data = parameters.all, exposure = 'no2')
exposure.nation.parameters.plot.o3 <- exposure.parameter.plot(data = parameters.all, exposure = 'o3')

setwd(paste0(dir.res.main.analysis, '/National'))
ggplot2::ggsave(plot = exposure.nation.parameters.plot.pm25,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'exposure.nation.parameters.plot.pm25'),
                                  '.png'),
                width = width, 
                height = height)
ggplot2::ggsave(plot = exposure.nation.parameters.plot.pm10,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'exposure.nation.parameters.plot.pm10'),
                                  '.png'),
                width = width, 
                height = height)
ggplot2::ggsave(plot = exposure.nation.parameters.plot.no2,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'exposure.nation.parameters.plot.no2'),
                                  '.png'),
                width = width, 
                height = height)
ggplot2::ggsave(plot = exposure.nation.parameters.plot.o3,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'exposure.nation.parameters.plot.o3'),
                                  '.png'),
                width = width, 
                height = height)

## 3.2. spatial  mixing ----

exposure.nation.mixing.plot.pm25 <- exposure.parameter.plot(data = parameters.all, exposure = 'pm25')
exposure.nation.mixing.plot.pm10 <- exposure.parameter.plot(data = parameters.all, exposure = 'pm10')
exposure.nation.mixing.plot.no2 <- exposure.parameter.plot(data = parameters.all, exposure = 'no2')
exposure.nation.mixing.plot.o3 <- exposure.parameter.plot(data = parameters.all, exposure = 'o3')

setwd(paste0(dir.res.main.analysis, '/National'))
ggplot2::ggsave(plot = exposure.nation.mixing.plot.pm25,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'exposure.nation.mixing.plot.pm25'),
                                  '.png'),
                width = width, 
                height = height)
ggplot2::ggsave(plot = exposure.nation.mixing.plot.pm10,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'exposure.nation.mixing.plot.pm10'),
                                  '.png'),
                width = width, 
                height = height)
ggplot2::ggsave(plot = exposure.nation.mixing.plot.no2,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'exposure.nation.mixing.plot.no2'),
                                  '.png'),
                width = width, 
                height = height)
ggplot2::ggsave(plot = exposure.nation.mixing.plot.o3,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'exposure.nation.mixing.plot.o3'),
                                  '.png'),
                width = width, 
                height = height)

## 3.3. table ----

table.parameters <-
  parameters.all %>% 
  dplyr::select(exposure, region, model, parameter.name.label , estimate) %>% 
  dplyr::arrange(exposure, region)

setwd(paste0(dir.res.main.analysis, '/National'))
print(xtable::xtable(table.parameters),
      include.rownames = FALSE,
      file = 'model_parameters_nation.txt')
