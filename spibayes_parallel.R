
# *------------------------------------------------------------------
# | PROGRAM NAME: Bayesian GAM model for non-linear trend analysis
# | FILE NAME: spibayes_parallel.R
# | DATE: Oct.22.2020
# | CREATED BY:  Kay Sung     
# *----------------------------------------------------------------
# | PURPOSE: Bayesian analysis for modeling seasonal and long-term
# |         trends of each gauges. It runs parallel using cores in
# |         individual computer.
# | 
# *------------------------------------------------------------------


### For custom MLE functions
require(here)
require(tidyverse)
require(here)

require(rnoaa)
require(tidyverse)
require(curl)
require(here)

#require(spibayes)
require(rstan)
require(mgcv)
require(viridis)
require(ggthemes)

### Packages for spi
require(fitdistrplus)
require(lubridate)
require(zoo)
require(cmdstanr)
select <- dplyr::select

theme_set(theme_classic(8))


############################  
########Define path
###############################  
here_path <- here::here()

data_path <- file.path(here_path, "./data")
output_path <- file.path(here_path, "./output")

### Set up figure folder
write_figures_path <- file.path(output_path, "figures/spline_synth")
dir.create(write_figures_path, recursive=TRUE, showWarnings = FALSE)  


##########################  
sapply(list.files(pattern="[.]R$", path="./functions", full.names=TRUE), source)

find_z <- function(orig_basis){
  C <- rep(1, nrow(orig_basis)) %*% orig_basis
  qrc <- qr(t(C))
  Z <- qr.Q(qrc,complete=TRUE)[,(nrow(C)+1):ncol(C)]
  return(Z)
}

bayes_tensor <- function(id) {
  #download ghcnd data
  show(id)
  options(noaakey = "GpOcGoGckLcTOEsaNABdZoRxXrOCPTkY") 
  station_data <- ghcnd_search(stationid = id, var = "PRCP")
  
  #Skip the gagues which have less then 85% of entire timeseries
  
  
  ###########################################################################
  ###  Create knots
  ###########################################################################
  ### Create knots
  n_knots_jdate <- 18
  n_knots_year <- 8
  
  knots_jdate <- seq(1,365,length=n_knots_jdate)
  knots_year <- seq(1925,2018,length=n_knots_year)
  
  knot_loc <- list(jdate = knots_jdate, year = knots_year)
  
  ###########################################################################
  ###  extract data
  ###########################################################################
  ### For now only work with gamma model (no hurdle model)
  n_roll <- 92
  n_roll_min <- n_roll - 14
  
  accum_df <- station_data$prcp %>%
    complete(date = seq.Date(min(date), max(date), by="1 day")) %>%
    arrange(date) %>%
    mutate(jdate = yday(date)) %>%
    mutate(month_day = paste0(month(date),"-",day(date))) %>%
    mutate(year = year(date))

# drop the rollmean value if number of NA precip is included more than 14
  
  fitting_df <- accum_df %>%
    filter(year >= 1925 | year  <= 2019)  %>%
    filter(jdate <= 365) %>% 
    mutate(roll_mean_3 = rollmeanr(x=prcp, k=n_roll, fill=NA, na.rm=TRUE)) %>%
    mutate(numdays_notna = rollsumr(x=!is.na(prcp), k=n_roll, fill=NA, na.rm=TRUE)) %>%
    mutate(roll_mean_3 = case_when(numdays_notna > n_roll_min ~ roll_mean_3,
                                   TRUE ~ NA_real_)) %>%
    mutate(precip = roll_mean_3/10) %>%
    mutate(zero = (precip == 0)) %>%
    select(id, date, year, jdate, zero,precip)
  
  data_all <- fitting_df %>%
    drop_na(precip)
  
  data_pos <- fitting_df %>%
    filter(precip > 0)
  
  num_data <- nrow(fitting_df)
  num_nona <- nrow(data_all)
  
  if (num_nona/ num_data < 0.85){
    return(id)
  }
  
  ###########################################################################
  ###  Create spline basis using full time series
  ###########################################################################
  
  ### Extract the number of knots
  n_knots <- sapply(knot_loc, FUN = length)
  
  ### Create tensor spline basis for precip
  ti_model_pos <- gam(list(precip 
                           ~ ti(jdate, bs=c("cc"), k = n_knots[1]) + ti(year, bs=c("cr"), k = n_knots[2]) + ti(jdate,year, bs=c("cc", "cr"), k = n_knots),
                           ~ ti(jdate, bs=c("cc"), k = n_knots[1]) + ti(year, bs=c("cr"), k = n_knots[2]) + ti(jdate,year, bs=c("cc", "cr"), k = n_knots)),
                      data = data_pos, 
                      family=gammals(link=list("identity","log")),
                      fit = FALSE, 
                      select = TRUE)
  
  ### Extract basis matrix for mean
  X_mean_jdate <- PredictMat(ti_model_pos$smooth[[1]],data_pos)
  X_mean_year <- PredictMat(ti_model_pos$smooth[[2]],data_pos)
  X_mean_tensor <- PredictMat(ti_model_pos$smooth[[3]],data_pos)
  
  ### Extract basis matrix for dispersion
  X_disp_jdate <- PredictMat(ti_model_pos$smooth[[4]],data_pos)
  X_disp_year <- PredictMat(ti_model_pos$smooth[[5]],data_pos)
  X_disp_tensor <- PredictMat(ti_model_pos$smooth[[6]],data_pos)
  
  ### Extract S Penalty matrices for positive
  s_mean_jdate <- ti_model_pos$S[[1]]
  s_mean_year <- ti_model_pos$S[[2]]
  s_mean_year_double <- ti_model_pos$S[[3]]
  s_mean_tensor_jdate <- ti_model_pos$S[[4]]
  s_mean_tensor_year <- ti_model_pos$S[[5]]
  
  s_disp_jdate <- ti_model_pos$S[[6]]
  s_disp_year <- ti_model_pos$S[[7]]
  s_disp_year_double <- ti_model_pos$S[[8]]
  s_disp_tensor_jdate <- ti_model_pos$S[[9]]
  s_disp_tensor_year <- ti_model_pos$S[[10]]
  
  ### Extract S scaling function
  s_scale_pos <- sapply(ti_model_pos$smooth, function(x){x$S.scale})
  #s_scale_zero <- sapply(ti_model_zero$smooth, function(x){x$S.scale})
  
  ### Extract basis dimensions
  basis_dim_mean <- c(dim(X_mean_jdate)[2], dim(X_mean_year)[2], dim(X_mean_tensor)[2])
  basis_dim_disp <- c(dim(X_disp_jdate)[2], dim(X_disp_year)[2], dim(X_disp_tensor)[2])
  
  
  ###########################################################################
  ###  Calculate prior and initial values from MLE fit
  ###########################################################################
  n_years <- length(unique(lubridate::year(data_all$date)))
  
  ### Fit a cyclic spline model using mgcv
  ti_init <- gam(list(precip 
                      ~ ti(jdate, bs=c("cc"), k = n_knots[1]),
                      ~ ti(jdate, bs=c("cc"), k = n_knots[1])),
                 data = data_pos, 
                 family=gammals(link=list("identity","log")),
                 fit = TRUE, 
                 select = TRUE)
  
  theta_model <- mgcv::gam(zero ~ 1, data=data_all, knots = list(jdate = knot_loc$jdate), select=TRUE,family=binomial)
  
  coef_init_df <- data.frame(t(coef(ti_init)))
  
  ### Process the beta intercept
  b_0_mean_init <- select(coef_init_df,contains("Intercept"))[1]
  b_0_disp_init <- select(coef_init_df,contains("Intercept"))[2]
  b_0_theta_init <- coef(theta_model)[1]
  
  ### Process the beta init
  b_init_mean_jdate <- c(unlist(select(coef_init_df,contains("ti.jdate"))))
  b_init_mean_year <- rep(0,dim(X_mean_year)[2])
  b_init_mean_tensor <- rep(0,dim(X_mean_tensor)[2])
  
  b_init_disp_jdate <- c(unlist(select(coef_init_df,contains("ti.1.jdate"))))
  b_init_disp_year <- rep(0,dim(X_disp_year)[2])
  b_init_disp_tensor <- rep(0,dim(X_disp_tensor)[2])
  
  ### Pull the estimated penalty from a jdate only model
  lambda_mean_jdate_est <- ti_init$sp[1]
  lambda_disp_jdate_est <- ti_init$sp[2]
  
  ### Penalize annual an order of magnitude higher
  year_pen_ratio <- 100
  
  lambda_mean_init <- c(lambda_mean_jdate_est, lambda_mean_jdate_est*year_pen_ratio, lambda_mean_jdate_est*year_pen_ratio^2, lambda_mean_jdate_est, lambda_mean_jdate_est*year_pen_ratio)
  lambda_disp_init <- c(lambda_disp_jdate_est, lambda_disp_jdate_est*year_pen_ratio, lambda_disp_jdate_est*year_pen_ratio^2, lambda_disp_jdate_est, lambda_disp_jdate_est*year_pen_ratio)
  
  ### Create matrices for gamma distributed priors on lambda
  lambda_mean_prior <- matrix(NA, 5, 2)
  ### Add the shape parameter for lambda. Tight near normal for everything except the double penalty
  lambda_mean_prior[,1] <- c(500, 5000, 0.5, 5000, 5000)
  
  ### Replicate for other lambdas
  lambda_disp_prior <-lambda_mean_prior
  
  ### Calculate the rate parameter
  lambda_mean_prior[,2] <-lambda_mean_prior[,1]/lambda_mean_init
  lambda_disp_prior[,2] <-lambda_disp_prior[,1]/lambda_disp_init
  
  
  ### Create output list and return
  init_vals_output <- list(
    b_0_mean = b_0_mean_init, 
    b_0_disp = b_0_disp_init, 
    b_0_theta = b_0_theta_init, 
    b_mean_jdate = b_init_mean_jdate, 
    b_mean_year = b_init_mean_year, 
    b_mean_tensor = b_init_mean_tensor, 
    b_disp_jdate = b_init_disp_jdate, 
    b_disp_year = b_init_disp_year, 
    b_disp_tensor = b_init_disp_tensor, 
    lambda_mean = lambda_mean_init, 
    lambda_disp = lambda_disp_init
  )
  
  ### b_0 prior
  b_0_prior <- as.matrix(summary(ti_init)$p.table[,1:2])
  b_0_prior <- rbind(b_0_prior, matrix(summary(theta_model)$p.table[1:2], 1, 2))
  b_0_prior[,2] <- b_0_prior[,2] * 10
  

  ###########################################################################
  ###  Set up to Run the model
  ###########################################################################
  
  ### Extract precipitation
  y_pos <-  data_pos$precip
  y_zero <- data_all$zero
  
  ### Calculate some necesary dimensions
  N <- length(y_zero)
  N_pos <- length(y_pos)
  
  ### Create the data to send to stan model
  data_fitting <- list(
    N = N, 
    N_pos = N_pos,
    
    basis_dim_mean = basis_dim_mean, 
    basis_dim_disp = basis_dim_disp,
    #basis_dim_theta = basis_dim_theta,
    
    y_pos= y_pos, 
    y_zero = y_zero,
    
    X_mean_jdate = X_mean_jdate, 
    X_mean_year = X_mean_year, 
    X_mean_tensor = X_mean_tensor, 
    
    X_disp_jdate = X_disp_jdate, 
    X_disp_year = X_disp_year, 
    X_disp_tensor = X_disp_tensor, 
    
    b_0_prior = b_0_prior,
    
    s_mean_jdate = s_mean_jdate,
    s_mean_year = s_mean_year,
    s_mean_year_double = s_mean_year_double,
    s_mean_tensor_jdate = s_mean_tensor_jdate,
    s_mean_tensor_year = s_mean_tensor_year,
    
    s_disp_jdate = s_disp_jdate,
    s_disp_year = s_disp_year,
    s_disp_year_double = s_disp_year_double,
    s_disp_tensor_jdate = s_disp_tensor_jdate,
    s_disp_tensor_year = s_disp_tensor_year,
    
    lambda_mean_prior = lambda_mean_prior,
    lambda_disp_prior = lambda_disp_prior
  )
  
  ### Make the penalties positive definite
  penalty_index <- which(startsWith(names(data_fitting), "s_"))
  for (j in penalty_index){
    data_fitting[[j]] <- as.matrix(Matrix::nearPD(data_fitting[[j]])$mat)
  }
  
  #str(data_fitting)
  
  ### Reparse for init vals
  init_vals <- list(init_vals_output)
  
  ### Set up for multiple chains
  n_chains <- 3
  
  for (j in seq(2,n_chains)){
    init_vals[[j]] <- init_vals[[1]]
    
    ### Only tweak the year penalties
    init_vals[[j]]$lambda_mean[c(2,3,5)] <- init_vals[[j]]$lambda_mean[c(2,3,5)] * exp(runif(3, -3,3))
    init_vals[[j]]$lambda_disp[c(2,3,5)] <- init_vals[[j]]$lambda_disp[c(2,3,5)] * exp(runif(3, -3,3))
  }
  
 # str(init_vals)
  
  
  ###########################################################################
  ###  Compile the model
  ###########################################################################
  
  #mod <- cmdstan_model("./models/tensor_ti_second.stan")
  mod <- cmdstan_model(file.path(here_path,"tensor_ti_six_nozero.stan"))
  #stan("./models/tensor_penalized.stan", data = data_fitting, iter = 20)
  
  ###########################################################################
  ###  Quick Run not full Bayesian
  ###########################################################################
  ### Good for checking the lambda priors
  
  # Start the clock!
  #ptm <- proc.time()
  
  ### Quick model fit without full Bayesian
  ### Use Stan's LBFGS algorithm
  fit_quick <- mod$optimize(
    data = data_fitting,
    seed = 123,
    refresh = 10,
    init = list(init_vals[[1]])
  )
  
  # Stop the clock
  #proc.time() - ptm
  
  ### Quick plot
  b_0_mean <- fit_quick$mle("b_0_mean")
  
  b_mean_jdate <- fit_quick$mle("b_mean_jdate")
  b_mean_year <- fit_quick$mle("b_mean_year")
  b_mean_tensor <- fit_quick$mle("b_mean_tensor")
  
  b_disp_jate <-fit_quick$mle("b_disp_jdate")
  b_disp_year <-fit_quick$mle("b_disp_year")
  b_disp_tensor <-fit_quick$mle("b_disp_tensor")
  
  b_para <- c(b_0_mean, b_mean_jdate,b_mean_year,b_mean_tensor,b_disp_jate, b_disp_year,b_disp_tensor)           
  
  filename <- paste0(id,"_para.rda")
  save(b_para, file = filename, path = output_path)
  
  ##### Check lambda
  #fit_quick$mle("lambda_mean")
  #fit_quick$mle("lambda_disp")
  
  ### No spline for zeros
  # b_0_theta <- fit_quick$mle("b_0_theta")
  # exp(b_0_theta) / (1 + exp(b_0_theta))
  # median(data_pos$theta)
  
  #filename <- paste0(id,"_tensor.rda")
  #save(data_pos, file = filename, path = write_figures_path)
  ###########################################################################
  ###  Run the full model
  ###########################################################################
  # Start the clock!
  # ptm <- proc.time()
  
  #n_chains
  
  ### Sample from the model
  # fit_pen <- mod$sample(
  #   data = data_fitting,
  #   init = init_vals,
  #   seed = 123,
  #   chains = n_chains,
  #   output_dir = file.path(here_path,"/output"), 
  #   iter_warmup = 400, 
  #   iter_sampling = 400, 
  #   refresh = 100,
  #   save_warmup = TRUE
  # )
  
  # Stop the clock
  # proc.time() - ptm
  
  ###########################################################################
  ###  Plotting
  ###########################################################################
   plot_df <- data_pos
   plot_df$mean_est <- X_mean_jdate %*% b_mean_jdate + b_0_mean
   plot_df$mean_init <- X_mean_jdate %*% init_vals[[1]]$b_mean_jdate + init_vals[[1]]$b_0_mean[[1]]
   plot_df$mean_est_full <- X_mean_jdate %*% b_mean_jdate + X_mean_year %*% b_mean_year + X_mean_tensor %*% b_mean_tensor + b_0_mean
   plot_df$year_partial <- X_mean_year %*% b_mean_year + b_0_mean
   plot_df$tensor_partial <- X_mean_tensor %*% b_mean_tensor# + b_0_mean

  filename <- paste0(id,"_plot.rda")
  save(plot_df, file = filename, path = output_path)
  
 # load("USW00012916_plot.rda")

  #ggplot(plot_df, aes(x=jdate)) + geom_line(aes(y=precip, group = year)) + geom_line(aes(y=exp(mean_est)), colour = "red")+ geom_line(aes(y=exp(mean_init)), colour = "blue")
  #ggsave("./output/bayesian_meanprecip.png", p, width =12.5, height = 8, dpi = 300)

 #p <- ggplot(plot_df, aes(x=jdate)) + geom_line(aes(y=precip, group = year),color = 'grey') + 
#   geom_line(aes(y=exp(mean_est_full), group = year)) + scale_color_viridis()
#   theme_classic()
 
# p
 
# filename <- paste0(id,"mean.png")
# ggsave(filename = filename, plot = p, path = write_figures_path, width =12.5, height = 8, dpi = 300)

# p <- ggplot(plot_df, aes(x=jdate, y = year)) + geom_raster(aes(fill=exp(mean_est_full))) + scale_fill_viridis()
# filename <- paste0(id,"mean_raster.png")
# ggsave(filename = filename, plot = p, path = write_figures_path, width =12.5, height = 8, dpi = 300)
 #print(filename)

}

##############################################################
#plan for multi process
##############################################################
require(furrr)
station <- read.csv("GHCNDgauges.csv")

station_list <- station %>%
  filter(variable == "PRCP") %>% 
  filter(X < 1920 & X.1 > 2018) %>%
  mutate(id = as.character(ID)) %>% 
  select(id, Latitude,Longitude)

#source(file.path(here_path,"/functions/bayes_tensor.R"))
numCores <- parallel::detectCores() - 1
plan(multicore, workers = numCores)


map(station_list$id, bayes_tensor)

