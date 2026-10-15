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

## 0.4. output arguments ----

# plot 
height <- width <- 10
text.size <- 10

# colour schemes
colours.diverging <- c('#007D5E', '#E9E0C8', '#C13E2E')
colours.exceedance <- c('#A0D6D0', '#009B92', '#00768D')
colours.sequential <- c('#E9E0C8', '#A175A7')

# table
digits.table <- 4

## 0.5. model results ----

### 0.5.1. national ----

setwd(paste0(dir.res.model.fit, '/National'))
parameters.national.pm25 <- readr::read_csv(file = 'parameters_national_pm25.csv', show_col_types = FALSE)
parameters.national.pm10 <- readr::read_csv(file = 'parameters_national_pm10.csv', show_col_types = FALSE)
parameters.national.no2 <- readr::read_csv(file = 'parameters_national_no2.csv', show_col_types = FALSE)
parameters.national.o3 <- readr::read_csv(file = 'parameters_national_o3.csv', show_col_types = FALSE)
parameters.national.o3plus <- readr::read_csv(file = 'parameters_national_o3plus.csv', show_col_types = FALSE)

### 0.5.2. regional -----

setwd(paste0(dir.res.model.fit, '/Regional'))
parameters.regional.pm25 <- readr::read_csv(file = 'parameters_regional_pm25.csv', show_col_types = FALSE)
parameters.regional.pm10 <- readr::read_csv(file = 'parameters_regional_pm10.csv', show_col_types = FALSE)
parameters.regional.no2 <- readr::read_csv(file = 'parameters_regional_no2.csv', show_col_types = FALSE)
parameters.regional.o3 <- readr::read_csv(file = 'parameters_regional_o3.csv', show_col_types = FALSE)
parameters.regional.o3plus <- readr::read_csv(file = 'parameters_regional_o3plus.csv', show_col_types = FALSE)

## 0.6. labelling ----

label.exposure <- 
  data.frame(level = c('pm25', 'pm10', 'no2', 'o3', 'o3plus'),
             label = c('beta[PM[2.5]]', 'beta[PM[10]]', 'beta[NO[2]]', 'beta[O[3]]', 'beta[O[3] + NO[2]]'))
label.model <- c('M0',
                 'M1-IMD', 'M1-QOF', 'M1-REACT', 'M1-LCRR', 'M1-EPD', 'M1-CACI', 'M1-AHAH',
                 'M2-QOF', 'M2-REACT', 'M2-LCRR', 'M2-EPD', 'M2-CACI', 'M2-AHAH',
                 'M3')
label.model.type <- c('M0', 'M1', 'M2', 'M3')
label.region <- 
  lookup.lsoa.rgn %>%
  dplyr::distinct(RGN11NM) %>%
  dplyr::filter(RGN11NM != 'Wales') %>% 
  dplyr::arrange() %>% 
  dplyr::pull() %>% 
  sort() %>% 
  c('England', .)

# 1. results ----

## 1.0. set up ----

### 1.0.1 functions ----

generate.colour.values <- function(n.groups) {
  # Fallback colors for 1 or 2 groups
  colours.default <- c('#E41A1C', '#377EB8') # Adjust these colors as needed
  
  if (n.groups == 1) {
    # Return the first color if only one group is requested
    return(colours.default[1])
  } else if (n.groups == 2) {
    # Return the first two colors for two groups
    return(colours.default[1:2])
  } else if (n.groups <= max(RColorBrewer::brewer.pal.info$maxcolors)) {
    # Use RColorBrewer for 3 to max colors
    return(RColorBrewer::brewer.pal(n.groups, 'Set3'))
  } else {
    # For more than the maximum supported colors in RColorBrewer, use a color ramp
    return(colorRampPalette(RColorBrewer::brewer.pal(8, 'Set3'))(n.groups))
  }
}

exposure.parameter.plot <- function(data, term = c('exposure', 'Phi for space_id')[1], text.size = text.size){
  
  # 0. function arguments ----
  
  # data <- parameters.all %>% dplyr::filter(exposure == 'PM[2.5]') %>% dplyr::filter(region == 'England')
  # term <- 'exposure'
  # text.size <- 10
  
  # 1. data sorting ----
  
  ## 1.1. plot data ----
  
  data.2 <- 
    data %>% 
    dplyr::filter(parameter.name.inla == {{term}}) %>% 
    dplyr::arrange(exposure, region, model) %>% 
    dplyr::mutate(index = 1:n(),
                  .by = c('exposure', 'region'))
  
  ## 1.2. x-axis names ----
  
  x.axis.breaks <-
    data.2 %>% 
    dplyr::select(model, model.type, index) %>% 
    dplyr::distinct() %>% 
    dplyr::arrange(index)
  
  x.axis.lines <-
    x.axis.breaks %>% 
    dplyr::group_by(model.type) %>% 
    dplyr::slice_max(index) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(!(model.type == 'M3')) %>% 
    dplyr::mutate(x = index + 0.5,
                  xend = index + 0.5,
                  y = -Inf,
                  yend = Inf)
  
  ## 1.3. y-axis names ----
  
  # y.axis.breaks <-
  #   c(floor(min(data.2$lower)*10)/10,
  #     0,
  #     ceiling(max(data.2$upper)*10)/10)
  
  ## 1.4. region colours ----
  
  number.region <- 
    data.2 %>% 
    dplyr::filter(region != 'England') %>% 
    dplyr::pull(region) %>% 
    unique() %>% 
    length()
  
  if(number.region == 0){
    colours.region <- 'black'
    legend.position <- 'none'
  } else if (number.region == 1){
    colours.region <- c('black', 'red3')
    legend.position <- 'bottom'
  } else if(number.region == 2){
    colours.region <- c('black', 'red3', 'blue3')
    legend.position <- 'bottom'
  } else if(number.region >= 3){
    colours.region <-
      number.region %>% 
      generate.colour.values() %>% 
      colorspace::darken(., amount  = 0.1) %>% 
      c('black', .)
    legend.position <- 'bottom'
  }
  
  ## 1.5. M1-IMD shaded region ----
  # Extract the lower and upper 95% CrI for M1-IMD for each exposure
  # This creates a data frame with one row per exposure, giving the 
  # lower and upper bounds of the M1-IMD shaded rectangle
  
  imd.ribbon <-
    data.2 %>%
    dplyr::filter(model == 'M1-IMD', region == 'England') %>%
    dplyr::select(exposure, lower, upper) %>%
    dplyr::rename(imd.lower = lower,
                  imd.upper = upper)
  
  # 2. plot ----
  
  plot <- 
    ggplot2::ggplot(data = data.2, aes(x = index, y = median)) +
    # M1-IMD shaded region (added first so it sits behind all other elements)
    ggplot2::geom_rect(data = imd.ribbon,
                       inherit.aes = FALSE,
                       aes(xmin = -Inf,
                           xmax = Inf,
                           ymin = imd.lower,
                           ymax = imd.upper),
                       fill = 'skyblue',
                       alpha = 0.3) +
    # point and error bar
    ggplot2::geom_point(aes(colour = region, group = interaction(region, exposure)), position = ggplot2::position_dodge(0.5)) +
    ggplot2::geom_errorbar(aes(ymin = lower, ymax = upper, colour = region, group = interaction(region, exposure)), width = 0, position = ggplot2::position_dodge(0.5)) +
    # colour
    ggplot2::scale_color_manual(name = '',
                                values = colours.region) +
    # axis sort
    ## x-axis (to be flipped to y-axis)
    ### labels
    ggplot2::scale_x_continuous(name = '',
                                labels = x.axis.breaks$model,
                                breaks = 1:nrow(x.axis.breaks),
                                trans = 'reverse') +
    
    ### lines
    ggplot2::geom_segment(data = x.axis.lines, 
                          aes(x = x, xend = xend, y = y, yend = yend),
                          colour = 'grey') +
    ## y-axis (to be flipped to x-axis)
    ggplot2::scale_y_continuous(name = '') +
    ## flip
    ggplot2::coord_flip() +
    # facet
    ggplot2::facet_wrap(~ exposure, 
                        nrow = 1,
                        scales = 'free_x',
                        labeller = ggplot2::label_parsed) +
    # theme
    my.theme(text = element_text(size = text.size),
             strip.text = element_text(size = 1.25 * text.size),
             legend.title = ggplot2::element_blank(),
             legend.position = legend.position)
  
  if(term == 'Phi for space_id'){
    plot 
  } else {
    plot +
      # line at y = 0
      ggplot2::geom_hline(yintercept = 0, colour = 'black', linetype = 'dashed')
  }
  
}

### 1.0.2. data ----

parameters.all <- 
  dplyr::bind_rows(parameters.national.pm25,
                   parameters.national.pm10,
                   parameters.national.no2,
                   parameters.national.o3,
                   parameters.national.o3plus,
                   parameters.regional.pm25,
                   parameters.regional.pm10,
                   parameters.regional.no2,
                   parameters.regional.o3,
                   parameters.regional.o3plus) %>% 
  dplyr::mutate(exposure = exposure %>% factor(., labels = label.exposure$label, levels = label.exposure$level),
                region = region %>% factor(., levels = label.region),
                model = model %>% factor(., levels = label.model),
                model.type = model.type %>% factor(., levels = label.model.type))

## 1.1. exposure ----

### 1.1.0. pm25 only ----

exposure.plot.regional <- exposure.parameter.plot(data = parameters.all %>% dplyr::filter(exposure == 'beta[PM[2.5]]'), term = 'exposure', text.size = text.size)
exposure.plot.national <- exposure.parameter.plot(data = parameters.all %>% dplyr::filter(exposure == 'beta[PM[2.5]]') %>% dplyr::filter(region == 'England'), term = 'exposure', text.size = text.size)

setwd(dir.res.main.analysis)
ggplot2::ggsave(plot = exposure.plot.regional,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'pm25.plot.regional'),
                                  '.png'),
                width = width, 
                height = height)
ggplot2::ggsave(plot = exposure.plot.national,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'pm25.plot.national'),
                                  '.png'),
                width = width, 
                height = height)

### 1.1.1. without o3plus ----

exposure.plot.regional <- exposure.parameter.plot(data = parameters.all %>% dplyr::filter(exposure != 'beta[O[3] + NO[2]]'), term = 'exposure', text.size = text.size)
exposure.plot.national <- exposure.parameter.plot(data = parameters.all %>% dplyr::filter(exposure != 'beta[O[3] + NO[2]]') %>% dplyr::filter(region == 'England'), term = 'exposure', text.size = text.size)

setwd(dir.res.main.analysis)
ggplot2::ggsave(plot = exposure.plot.regional,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'exposure.plot.regional'),
                                  '.png'),
                width = width, 
                height = height)
ggplot2::ggsave(plot = exposure.plot.national,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'exposure.plot.national'),
                                  '.png'),
                width = width, 
                height = height)

### 1.1.2. with o3plus ----

exposure.plot.regional <- exposure.parameter.plot(data = parameters.all, term = 'exposure', text.size = text.size)
exposure.plot.national <- exposure.parameter.plot(data = parameters.all %>% dplyr::filter(region == 'England'), term = 'exposure', text.size = text.size)

setwd(dir.res.main.analysis)
ggplot2::ggsave(plot = exposure.plot.regional,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'exposure.plot.regional'),
                                  '_o3plus.png'),
                width = width, 
                height = height)
ggplot2::ggsave(plot = exposure.plot.national,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'exposure.plot.national'),
                                  '_o3plus.png'),
                width = width, 
                height = height)

## 1.2 spatial mixing ----

### 1.2.1 without o3plus ----

spatial.mixing.plot.regional <- exposure.parameter.plot(data = parameters.all %>% dplyr::filter(exposure != 'beta[O[3] + NO[2]]'), term = 'Phi for space_id', text.size = text.size)
spatial.mixing.plot.national <- exposure.parameter.plot(data = parameters.all %>% dplyr::filter(exposure != 'beta[O[3] + NO[2]]') %>% dplyr::filter(region == 'England'), term = 'Phi for space_id', text.size = text.size)

setwd(dir.res.main.analysis)
ggplot2::ggsave(plot = spatial.mixing.plot.regional,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'spatial.mixing.plot.regional'),
                                  '.png'),
                width = width, 
                height = height)
ggplot2::ggsave(plot = spatial.mixing.plot.national,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'spatial.mixing.plot.national'),
                                  '.png'),
                width = width, 
                height = height)

### 1.2.1 with o3plus ----

spatial.mixing.plot.regional <- exposure.parameter.plot(data = parameters.all, term = 'Phi for space_id', text.size = text.size)
spatial.mixing.plot.national <- exposure.parameter.plot(data = parameters.all %>% dplyr::filter(region == 'England'), term = 'Phi for space_id', text.size = text.size)

setwd(dir.res.main.analysis)
ggplot2::ggsave(plot = spatial.mixing.plot.regional,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'spatial.mixing.plot.regional'),
                                  '_o3plus.png'),
                width = width, 
                height = height)
ggplot2::ggsave(plot = spatial.mixing.plot.national,
                filename = paste0('model_parameters_', 
                                  gsub(pattern = '\\.', replacement = '_', x = 'spatial.mixing.plot.national'),
                                  '_o3plus.png'),
                width = width, 
                height = height)

## 1.3 all parameters table ----

### 1.3.0. pm25 only ---- 

table.parameters <-
  parameters.all %>% 
  dplyr::filter(exposure == 'beta[PM[2.5]]') %>% 
  dplyr::select(exposure, region, model, parameter.name.label , estimate) %>% 
  dplyr::arrange(exposure, region)

setwd(dir.res.main.analysis)
print(xtable::xtable(table.parameters),
      include.rownames = FALSE,
      file = 'model_parameters_region_pm25.txt')


table.parameters.all.csv <-
  parameters.all %>% 
  dplyr::filter(exposure != 'beta[PM[2.5]]') %>% 
  dplyr::select(exposure, region, model, parameter.name.label , median, lower, upper, estimate) %>% 
  dplyr::arrange(exposure, region)

table.parameters.exposure.csv <-
  parameters.all %>% 
  dplyr::filter(exposure == 'beta[PM[2.5]]') %>% 
  dplyr::filter(parameter.name.label == 'Exposure') %>% 
  dplyr::select(exposure, region, model, parameter.name.label , median, lower, upper, estimate) %>% 
  dplyr::arrange(exposure, region)

setwd(dir.res.main.analysis)
readr::write_csv(x = table.parameters.all.csv,
                 file = 'model_parameters_region_all_pm25.csv')
readr::write_csv(x = table.parameters.exposure.csv,
                 file = 'model_parameters_region_exposure_pm25.csv')

### 1.3.1. without o3plus ----

table.parameters <-
  parameters.all %>% 
  dplyr::filter(exposure != 'beta[O[3] + NO[2]]') %>% 
  dplyr::select(exposure, region, model, parameter.name.label , estimate) %>% 
  dplyr::arrange(exposure, region)

setwd(dir.res.main.analysis)
print(xtable::xtable(table.parameters),
      include.rownames = FALSE,
      file = 'model_parameters_region.txt')


table.parameters.all.csv <-
  parameters.all %>% 
  dplyr::filter(exposure != 'beta[O[3] + NO[2]]') %>% 
  dplyr::select(exposure, region, model, parameter.name.label , median, lower, upper, estimate) %>% 
  dplyr::arrange(exposure, region)

table.parameters.exposure.csv <-
  parameters.all %>% 
  dplyr::filter(exposure != 'beta[O[3] + NO[2]]') %>% 
  dplyr::filter(parameter.name.label == 'Exposure') %>% 
  dplyr::select(exposure, region, model, parameter.name.label , median, lower, upper, estimate) %>% 
  dplyr::arrange(exposure, region)

setwd(dir.res.main.analysis)
readr::write_csv(x = table.parameters.all.csv,
                 file = 'model_parameters_region_all.csv')
readr::write_csv(x = table.parameters.exposure.csv,
                 file = 'model_parameters_region_exposure.csv')

### 1.3.1. with o3plus ----

table.parameters <-
  parameters.all %>% 
  dplyr::select(exposure, region, model, parameter.name.label , estimate) %>% 
  dplyr::arrange(exposure, region)

setwd(dir.res.main.analysis)
print(xtable::xtable(table.parameters),
      include.rownames = FALSE,
      file = 'model_parameters_region_o3plus.txt')


table.parameters.all.csv <-
  parameters.all %>% 
  dplyr::select(exposure, region, model, parameter.name.label , median, lower, upper, estimate) %>% 
  dplyr::arrange(exposure, region)

table.parameters.exposure.csv <-
  parameters.all %>% 
  dplyr::filter(parameter.name.label == 'Exposure') %>% 
  dplyr::select(exposure, region, model, parameter.name.label , median, lower, upper, estimate) %>% 
  dplyr::arrange(exposure, region)

setwd(dir.res.main.analysis)
readr::write_csv(x = table.parameters.all.csv,
                 file = 'model_parameters_region_all_o3plus.csv')
readr::write_csv(x = table.parameters.exposure.csv,
                 file = 'model_parameters_region_exposure_o3plus.csv')
