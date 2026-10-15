# descriptive statistics ----

descriptive.function <- function(x){
  
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

# link and inv-link ----

logit <- function(x){
  
  log(x / (1 - x))
  
}

expit <- function(x){
  
  exp(x) / (1 + exp(x))
  
}

# get the adjacency matrix ----

get.amat <- function (poly){
  
  mat.temp <- spdep::poly2nb(as(poly, 'Spatial'))
  mat <- spdep::nb2mat(mat.temp, zero.policy = TRUE)
  
  return(mat)
}

# themes for plots ----

my.theme<-function(...){
  theme(panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour = 'black'),
        # legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(hjust = 0),
        legend.key = ggplot2::element_rect(fill = NA, colour = NA),
        ...)
}

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

# model fit ----

model.fit <- function(formula.update,
                      data,
                      family = 'poisson',
                      a.mat = NULL,
                      pc.u = 1, pc.alpha = 0.01,
                      pu.u.phi = 0.5, pc.alpha.phi = 2/3,
                      control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE),
                      control.predictor = list(compute = TRUE, link = 1),
                      control.inla = list(strategy = 'auto', int.strategy = 'auto'),
                      control.mode = list(restart = TRUE),
                      inla.mode = c('classic', 'twostage', 'experimental')[3],
                      verbose = FALSE){
  
  # 0. function arguments ----
  
  # formula <- 'Y ~ 1'
  # family <- 'poisson'
  # pc.u <- 1; pc.alpha <- 0.01; pu.u.phi <- 0.5; pc.alpha.phi <- 2/3
  # control.compute <- list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE)
  # control.predictor <- list(compute = TRUE, link = 1)
  # control.inla <- list(strategy = 'auto', int.strategy = 'auto')
  # control.mode <- list(restart = TRUE)
  # inla.mode <- c('classic', 'twostage', 'experimental')[3]
  # verbose <- FALSE
  
  # 1. model parameters ----
  
  hyper.pc <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)))
  hyper.pc.space <- list(prec = list(prior = 'pc.prec', param = c(pc.u, pc.alpha)), phi = list(prior = 'pc', param = c(pu.u.phi, pc.alpha.phi)))
  
  # 2. model formula ----
  
  formula.base <- Y ~ 1 + exposure + f(space_id, model = 'bym2', graph = a.mat, constr = TRUE, scale.model = TRUE, hyper = hyper.pc.space)
  # formula.base <- Y ~ 1 + exposure + f(space_id, model = 'iid', hyper = hyper.pc)
  
  # update formula if NOT NA
  if(is.na(formula.update)){
    formula.final <- formula.base
  } else {
    formula.final <- update(formula.base, paste('~ . +', formula.update))
  }
  
  # 3. model fit ----
  
  fit <- 
    INLA::inla(formula = formula.final,
               family = family,
               data = data,
               E = E,
               control.compute = control.compute,
               control.predictor = control.predictor,
               control.inla = control.inla,
               control.mode = control.mode,
               inla.mode = inla.mode,
               num.threads = round(parallel::detectCores()*0.8),
               verbose = verbose)
  
  # 4. extract parameters ----
  
  parameters <- 
    dplyr::bind_rows(
      fit$summary.fixed[,c(3:5)] %>% 
        tibble::rownames_to_column(),
      fit$summary.hyperpar[,c(3:5)] %>% 
        tibble::rownames_to_column()) %>% 
    dplyr::mutate(median = `0.5quant`,
                  lower = `0.025quant`,
                  upper = `0.975quant`,
                  estimate = paste0(median %>% round(., digits = 4), ' (',
                                    lower %>% round(., digits = 4), ', ',
                                    upper %>% round(., digits = 4), ')')) %>% 
    dplyr::select(rowname, median, lower, upper, estimate)
  
  # 5. extract scores ----
  
  scores <-
    data.frame(dic = fit$dic$dic,
               waic = fit$waic$waic,
               cpo = fit$cpo$cpo %>% log() %>% mean() %>% {. * -1})
  
  # 6. return ----
  
  return(list(formula = formula, fit = fit, parameters = parameters, scores = scores))
  
}

# posterior distribution summary ----

my.summary <- function(x, CI = 0.95) {
  
  # x = theta1[,1]; CI = 0.95; 
  
  # lower and upper CI calculators
  lowerCI <- (1 - CI)/2
  upperCI <- 1 - lowerCI
  
  qntl <- quantile(x, probs = c(lowerCI, 0.5, upperCI))
  data.frame(mean = mean(x), variance = var(x), lower = qntl[1], median = qntl[2], upper = qntl[3])
  
}