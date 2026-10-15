# packages ----

if(!require('tidyverse')) {
  install.packages('tidyverse', dependencies = TRUE)
}
if(!require('INLA')) {
  install.packages('INLA', dependencies = TRUE)
}

# directories ----

# retrieve directories
code.path <- rstudioapi::getActiveDocumentContext()$path
dir.home <- sub('(NIH_Smoking_Proxies).*', '\\1', code.path)
dir.code <- paste0(dir.home, '/Code')
dir.data <- paste0(dir.home, '/Data')
dir.res <- paste0(dir.home, '/Results')
dir.data.organised <- paste0(dir.data, '/Organised Data')
dir.res.exposure.selection <- paste0(dir.res, '/Exposure Selection')

# imports ----

## functions ----

setwd(dir.code)
source('00_functions.R')

## organised data ----

setwd(dir.data.organised)
load('LSOA11_HES_FINAL.rda')
load('ENGLAND_SPATIAL_INFORMATION_LSOA_MSOA_LAD.rda')
load('ENGLAND_AMAT_LSOA_MSOA_LAD.rda')
rm(msoa.info, msoa.mat, lad.info, lad.mat) # remove unneeded variables (space)

# output arguments ----

# plot
height <- width <- 10

# table
digits.table <- 4

# model fitting ----

## inla arguments ----

pc.u <- 1
pc.alpha <- 0.01
pu.u.phi <- 0.5
pc.alpha.phi <- 2/3
formula <- 
  Y ~ 
  1 + 
  pollution +
  f(space_id, 
    model = 'bym2', 
    graph = lsoa.mat, 
    constr = TRUE,
    scale.model = TRUE,
    hyper = list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)),
                 phi = list(prior = 'pc', param = c(pu.u.phi, pc.alpha.phi))))
family <- 'poisson'
control.compute <- list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE)
control.predictor <- list(compute = TRUE, link = 1)
control.inla <- list(strategy = 'auto', int.strategy = 'auto')
control.mode <- list(restart = TRUE)
inla.mode <- c('classic', 'twostage', 'experimental')[3]
verbose <- FALSE

## PM 2.5 ----

data.final.pm25 <-
  data.final %>% 
  dplyr::mutate(pollution = pm25)

fit.pm25 <- 
  INLA::inla(formula = formula,
             family = family,
             data = data.final.pm25,
             E = E,
             control.compute = control.compute,
             control.predictor = control.predictor,
             control.inla = control.inla,
             control.mode = control.mode,
             inla.mode = inla.mode,
             num.threads = round(parallel::detectCores()*0.8),
             verbose = verbose)

## PM 10 ----

data.final.pm10 <-
  data.final %>% 
  dplyr::mutate(pollution = pm10)

fit.pm10 <- 
  INLA::inla(formula = formula,
             family = family,
             data = data.final.pm10,
             E = E,
             control.compute = control.compute,
             control.predictor = control.predictor,
             control.inla = control.inla,
             control.mode = control.mode,
             inla.mode = inla.mode,
             num.threads = round(parallel::detectCores()*0.8),
             verbose = verbose)

## NO2 ----

data.final.no2 <-
  data.final %>% 
  dplyr::mutate(pollution = no2)

fit.no2 <- 
  INLA::inla(formula = formula,
             family = family,
             data = data.final.no2,
             E = E,
             control.compute = control.compute,
             control.predictor = control.predictor,
             control.inla = control.inla,
             control.mode = control.mode,
             inla.mode = inla.mode,
             num.threads = round(parallel::detectCores()*0.8),
             verbose = verbose)

## O3 ----

data.final.o3 <-
  data.final %>% 
  dplyr::mutate(pollution = o3)

fit.o3 <- 
  INLA::inla(formula = formula,
             family = family,
             data = data.final.o3,
             E = E,
             control.compute = control.compute,
             control.predictor = control.predictor,
             control.inla = control.inla,
             control.mode = control.mode,
             inla.mode = inla.mode,
             num.threads = round(parallel::detectCores()*0.8),
             verbose = verbose)

# results ----

## pm2.5 ----

parameters.pm25 <- 
  dplyr::bind_rows(
    fit.pm25$summary.fixed[,c(3:5)] %>% 
      tibble::rownames_to_column(),
    fit.pm25$summary.hyperpar[,c(3:5)] %>% 
      tibble::rownames_to_column()) %>% 
  dplyr::mutate(median = `0.5quant`,
                lower = `0.025quant`,
                upper = `0.975quant`,
                estimate = paste0(median %>% round(., digits = digits.table), ' (',
                                  lower %>% round(., digits = digits.table), ', ',
                                  upper %>% round(., digits = digits.table), ')')) %>% 
  dplyr::select(rowname, median, lower, upper, estimate)

## pm10 ----

parameters.pm10 <- 
  dplyr::bind_rows(
    fit.pm10$summary.fixed[,c(3:5)] %>% 
      tibble::rownames_to_column(),
    fit.pm10$summary.hyperpar[,c(3:5)] %>% 
      tibble::rownames_to_column()) %>% 
  dplyr::mutate(median = `0.5quant`,
                lower = `0.025quant`,
                upper = `0.975quant`,
                estimate = paste0(median %>% round(., digits = digits.table), ' (',
                                  lower %>% round(., digits = digits.table), ', ',
                                  upper %>% round(., digits = digits.table), ')')) %>% 
  dplyr::select(rowname, median, lower, upper, estimate)

## no2 ----

parameters.no2 <- 
  dplyr::bind_rows(
    fit.no2$summary.fixed[,c(3:5)] %>% 
      tibble::rownames_to_column(),
    fit.no2$summary.hyperpar[,c(3:5)] %>% 
      tibble::rownames_to_column()) %>% 
  dplyr::mutate(median = `0.5quant`,
                lower = `0.025quant`,
                upper = `0.975quant`,
                estimate = paste0(median %>% round(., digits = digits.table), ' (',
                                  lower %>% round(., digits = digits.table), ', ',
                                  upper %>% round(., digits = digits.table), ')')) %>% 
  dplyr::select(rowname, median, lower, upper, estimate)

## o3 ----

parameters.o3 <- 
  dplyr::bind_rows(
    fit.o3$summary.fixed[,c(3:5)] %>% 
      tibble::rownames_to_column(),
    fit.o3$summary.hyperpar[,c(3:5)] %>% 
      tibble::rownames_to_column()) %>% 
  dplyr::mutate(median = `0.5quant`,
                lower = `0.025quant`,
                upper = `0.975quant`,
                estimate = paste0(median %>% round(., digits = digits.table), ' (',
                                  lower %>% round(., digits = digits.table), ', ',
                                  upper %>% round(., digits = digits.table), ')')) %>% 
  dplyr::select(rowname, median, lower, upper, estimate)

## all pollutants ----

parameter.names.inla <- c('(Intercept)', 'pollution', 'Precision for space_id', 'Phi for space_id')
parameter.names.label <- c('Intercept', 'Pollution', 'Precision', 'Phi')
# parameter.names.inla <- c('pollution')
# parameter.names.label <- c('Pollution')

parameters.all.template <-
  data.frame(rowname = parameter.names.inla,
             parameter = parameter.names.label %>% factor(., levels = parameter.names.label, labels = parameter.names.label))

parameters.all.full <-
  parameters.all.template %>% 
  dplyr::left_join(.,
                   parameters.pm25,
                   by = 'rowname') %>% 
  dplyr::left_join(.,
                   parameters.pm10,
                   suffix = c('.PM25', ''),
                   by = 'rowname') %>% 
  dplyr::left_join(.,
                   parameters.no2,
                   suffix = c('.PM10', ''),
                   by = 'rowname') %>% 
  dplyr::left_join(.,
                   parameters.o3,
                   suffix = c('.NO2', '.O3'),
                   by = 'rowname')

#### plot ----

parameters.plot.data <- 
  parameters.all.full %>% 
  dplyr::filter(parameter %in% c('Pollution', 'Phi')) %>% 
  dplyr::select(-dplyr::starts_with('estimate'), -rowname) %>% 
  dplyr::mutate(dplyr::across(-c(parameter), ~ dplyr::if_else(is.na(.x), 0, .x))) %>%
  tidyr::pivot_longer(cols = -c(parameter),
                      names_to = c('metric', 'model'),
                      names_pattern = '(.*)\\.(.*)',
                      values_to = 'value') %>% 
  tidyr::pivot_wider(names_from = 'metric',
                     values_from = 'value') %>% 
  dplyr::mutate(model = model %>% factor(., levels = c('PM25', 'PM10', 'NO2', 'O3'), labels = c('PM2.5', 'PM10', 'NO2', 'O3'))) %>% 
  # ORDER FIRST THEN INDEX
  dplyr::arrange(parameter, model) %>%
  dplyr::group_by(model) %>% 
  dplyr::mutate(index = 1:n()) %>% 
  dplyr::ungroup()

parameters.axis.names <-
  parameters.plot.data %>% 
  dplyr::select(parameter) %>% 
  dplyr::distinct() %>% 
  dplyr::arrange(parameter) %>% 
  dplyr::mutate(level = dplyr::case_when(parameter == 'Intercept' ~ '',
                                         parameter == 'Pollution'~ '',
                                         parameter == 'Precision'~ '',
                                         parameter == 'Phi'~ '',
                                         .default = parameter %>% as.character()))

parameters.plot <-
  ggplot2::ggplot(data = parameters.plot.data %>% dplyr::filter(parameter != 'Phi'), aes(x = index, y = median)) +
  ggplot2::geom_hline(yintercept = 0, colour = 'black', linetype = 'dashed') +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5), aes(group = model, colour = model)) +
  ggplot2::geom_errorbar(aes(ymin = lower, ymax = upper, group = model, colour = model), width = 0, position = ggplot2::position_dodge(0.5)) +
  ggplot2::scale_colour_manual(values = c('red3', 'blue3', 'green4', 'purple3')) +
  ggplot2::labs(x = '', y = '') +
  # plot annotation
  ## axis labels
  ggplot2::scale_x_continuous(labels = parameters.axis.names$level,
                              breaks = 1:nrow(parameters.axis.names),
                              trans = 'reverse') +
  # ## splitting by exposure
  # ### lines
  # ggplot2::geom_segment(data = data.frame(y = rep(-Inf, times = 1), yend = rep(Inf, times = 1),
  #                                         x = c(1.5),
  #                                         xend = c(1.5)),
  #                       aes(x = x, xend = xend, y = y, yend = yend),
  #                       colour = 'grey') +
  # ### labels
  # ggplot2::annotate('text',
  #                   y = rep(-0.5, times = 2),
  #                   x = c(1, 2),
  #                   label = parameter.names.label,
  #                   fontface = 'bold',
  #                   colour = 'grey',
  #                   size = 4,
  #                   hjust = 0) +
  # make sure coor_flip at the end
  ggplot2::coord_flip() +
  my.theme(legend.title = ggplot2::element_blank(),
           # text = ggplot2::element_text(size = text.size),
           legend.position = 'bottom'); parameters.plot

phi.plot <-
  ggplot2::ggplot(data = parameters.plot.data %>% dplyr::filter(parameter == 'Phi'), aes(x = index, y = median)) +
  # ggplot2::geom_hline(yintercept = 0.5, colour = 'black', linetype = 'dashed') +
  ggplot2::geom_point(position = ggplot2::position_dodge(0.5), aes(group = model, colour = model)) +
  ggplot2::geom_errorbar(aes(ymin = lower, ymax = upper, group = model, colour = model), width = 0, position = ggplot2::position_dodge(0.5)) +
  ggplot2::scale_colour_manual(values = c('red3', 'blue3', 'green4', 'purple3')) +
  ggplot2::labs(x = '', y = '') +
  # plot annotation
  ## axis labels
  ggplot2::scale_x_continuous(labels = parameters.axis.names$level,
                              breaks = 1:nrow(parameters.axis.names),
                              trans = 'reverse') +
  # ## splitting by exposure
  # ### lines
  # ggplot2::geom_segment(data = data.frame(y = rep(-Inf, times = 1), yend = rep(Inf, times = 1),
  #                                         x = c(1.5),
  #                                         xend = c(1.5)),
  #                       aes(x = x, xend = xend, y = y, yend = yend),
  #                       colour = 'grey') +
  # ### labels
  # ggplot2::annotate('text',
  #                   y = rep(-0.5, times = 2),
  #                   x = c(1, 2),
  #                   label = parameter.names.label,
  #                   fontface = 'bold',
  #                   colour = 'grey',
  #                   size = 4,
  #                   hjust = 0) +
  # make sure coor_flip at the end
  ggplot2::coord_flip() +
  # ggplot2::facet_grid(~ parameter, scale = 'free') +
  my.theme(legend.title = ggplot2::element_blank(),
           # text = ggplot2::element_text(size = text.size),
           legend.position = 'bottom'); phi.plot

setwd(dir.res.exposure.selection)
ggplot2::ggsave(plot = parameters.plot,
                filename = 'exposure_parameters_regression.png',
                width = width, 
                height = height)
ggplot2::ggsave(plot = phi.plot,
                filename = 'exposure_parameters_spatial_mixing.png',
                width = width, 
                height = height)

#### table ----

parameters.table <-
  parameters.all.full %>% 
  dplyr::select(parameter, estimate.PM25, estimate.PM10, estimate.NO2, estimate.O3) %>% 
  dplyr::mutate(dplyr::across(-c(parameter), ~ dplyr::if_else(is.na(.x), '-', .x)))

setwd(dir.res.exposure.selection)
print(xtable::xtable(parameters.table),
      include.rownames = FALSE,
      file = 'exposure_parameters.txt')



