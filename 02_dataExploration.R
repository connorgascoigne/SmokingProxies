# packages ----

if(!require('tidyverse')) {
  install.packages('tidyverse', dependencies = TRUE)
}
if(!require('reshape2')) {
  install.packages('reshape2', dependencies = TRUE)
}

# directories ----

# retrieve directories
code.path <- rstudioapi::getActiveDocumentContext()$path
dir.home <- sub('(NIH_Smoking_Proxies).*', '\\1', code.path)
dir.code <- paste0(dir.home, '/Code')
dir.data <- paste0(dir.home, '/Data')
dir.res <- paste0(dir.home, '/Results')
dir.data.spatial <- paste0(dir.data, '/shapeFiles')
dir.data.organised <- paste0(dir.data, '/Organised Data')
dir.res.data.exploration <- paste0(dir.res, '/Data Exploration')

# imports ----

## functions ----

setwd(dir.code)
source('00_functions.R')

## organised data ----

setwd(dir.data.organised)
load('LSOA11_HES_FINAL.rda')

# output arguments ----

# plot 
height <- width <- 10
text.size <- 20

# colour schemes
colours.diverging <- c('#007D5E', '#E9E0C8', '#C13E2E')
colours.exceedance <- c('#A0D6D0', '#009B92', '#00768D')
colours.sequential <- c('#E9E0C8', '#A175A7')

# table
digits.table <- 4

# correlation ----

parameter.code.name <- c('imd', 'qof', 'react', 'lcrr', 'epd', 'caci', 'ahah')
parameter.code.name.id <- c(paste0(parameter.code.name, '_id'), 'pm25', 'pm10', 'no2', 'o3')
parameter.plot.name <- c(
  'imd_id' = 'IMD', 
  'qof_id' = 'QOF',
  'react_id' = 'REACT', 
  'lcrr_id' = 'LCRR', 
  'epd_id' = 'EPD', 
  'caci_id' = 'CACI', 
  'ahah_id' = 'AHAH',
  'pm25' = 'PM[2.5]', 
  'pm10' = 'PM[10]', 
  'no2' = 'NO[2]', 
  'o3' = 'O[3]'
)

correlation.matrix.data.full <-
  data.final %>% 
  dplyr::select(dplyr::all_of(parameter.code.name.id)) %>% 
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.numeric(.x))) %>% 
  cor(., method = 'spearman')

get.lower.triangle <- function(matrix){
  matrix[lower.tri(matrix)] <- NA
  return(matrix)
}

correlation.matrix.data.lower <-
  get.lower.triangle(correlation.matrix.data.full)

correlation.matrix.plot.data <-
  correlation.matrix.data.lower %>% 
  reshape2::melt(., varnames = c('row', 'column'), value.name = 'correlation', na.rm = TRUE) %>% 
  tibble::as_tibble() # %>% 
# dplyr::mutate(row = row %>% factor(., levels = parameter.code.name.id, labels = parameter.plot.name),
#               column = column %>% factor(., levels = parameter.code.name.id, labels = parameter.plot.name))

setwd(dir.res.data.exploration)
save(correlation.matrix.data.full, file = 'correlation_data_full.rda')
save(correlation.matrix.plot.data, file = 'correlation_data_plot.rda')

correlation.matrix.plot <-
  ggplot2::ggplot(data = correlation.matrix.plot.data, aes(x = row, y = column, fill = correlation, label = round(correlation, digits = 2))) +
  ggplot2::geom_tile(colour = 'white') +
  ggplot2::scale_fill_gradient2(name = '',
                                midpoint = 0,
                                low = colours.diverging[3], 
                                mid = colours.diverging[2],
                                high = colours.diverging[1],
                                limits = c(-1, 1)) +
  ggplot2::geom_text(colour = 'black', size = 4) +
  ggplot2::coord_flip() +
  ggplot2::scale_x_discrete(
    limits = parameter.code.name.id %>% rev(),
    labels = function(x) parse(text = parameter.plot.name[x])
  ) +
  ggplot2::scale_y_discrete(
    # limits = parameter.plot.name %>% rev(),
    labels = function(y) parse(text = parameter.plot.name[y])
  ) +
  ggplot2::labs(y = '', x = '') +
  ggplot2::theme(text = ggplot2::element_text(size = text.size),
                 legend.position = 'bottom',
                 legend.key.width = unit(2, 'cm'),
                 axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
                 axis.ticks = element_blank(),
                 panel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 axis.line = ggplot2::element_blank(),
                 legend.text = ggplot2::element_text(hjust = 0),
                 legend.key = ggplot2::element_rect(fill = NA, colour = NA),
                 legend.title = ggplot2::element_blank()); correlation.matrix.plot

setwd(dir.res.data.exploration)
ggplot2::ggsave(plot = correlation.matrix.plot,
                filename = 'correlation_matrix.png',
                width = width, height = height)
