

### Trelliscope
### ToxPi
# Reduce to chemical class rather than abbrev, 
# aggregate concentrations?  mols can be added



#logistic:  detected or not, plot prob of prediction, correlated ot others?
# 



library(data.table)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(sf)
library(maps)
library(lubridate)
library(scales)
library(raster)
library(rgdal)



##### Plot over the State ####
PFAS <- read.csv("/Users/zilberds/Desktop/PFAS/nc_pfas_joined.csv")
states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))


counties <- st_as_sf(map("county", plot = FALSE, fill = TRUE))
counties <- subset(counties, grepl("north carolina", counties$ID))


ggplot() +
  geom_sf(data = counties, fill = NA, color = gray(.5))+
  geom_point(data = PFAS,aes(nc_public_water_Longitude,
                             nc_public_water_Latitude))

#### View Gamma vs NB vs Exp dist  ######
y_obs = PFAS$Conc_ppt[!is.na(PFAS$Conc_ppt)]
d_thresh = 1
plot(density(y_obs, bw = .5), xlim = c(0, 100))
mu = mean(y_obs)
varnc = var(y_obs)
param_p = mu/varnc
param_r = mu^2/(varnc-mu)
param_s = varnc - mu # = param_p*mu/(1-param_p)
hist(y_obs, breaks = 2500, xlim = c(0,50), probability = T)
scale_factor_ex = 1/(1-pexp(1, rate = 1/mu))
lines(1:100,dexp(1:100, rate = 1/mu)*scale_factor_ex, col=2, lwd=2)
scale_factor = 1/(1-pnbinom(1, mu = mu, size =param_r))
lines(1:100,dnbinom(1:100, mu = mu, size =param_r)*scale_factor, col=3, lwd=2)
scale_factor_g = 1/(1- pgamma(1,shape =mu^2/varnc, scale = 1/param_p))
lines(1:100,dgamma(1:100, shape =mu^2/varnc, scale = 1/param_p)*scale_factor_g,col=4)

####  check what a GP does without link #####

input_abbr="PFBS"
mdf = PFAS[which(PFAS$Abbreviation==input_abbr),]
cc_mdf = mdf[!is.na(mdf$Conc_ppt),]



x_pred = seq(min(PFAS$nc_public_water_Longitude), max(PFAS$nc_public_water_Longitude), length.out=100)
y_pred = seq(min(PFAS$nc_public_water_Latitude), max(PFAS$nc_public_water_Latitude), length.out=100)
pred_loc = expand.grid(x_pred, y_pred)
obs_val = cc_mdf$Conc_ppt
obs_loc = cbind(cc_mdf$nc_public_water_Longitude, cc_mdf$nc_public_water_Latitude)
K_nn = Matern(rdist(obs_loc), range = 1, smoothness = 1)+diag(rep(.001, length(obs_val)))
K_pn = Matern(rdist(pred_loc, obs_loc), range = 6, smoothness = 1)
mean_pred = K_pn%*%solve(K_nn, obs_val)

NC_state = states[states$ID=="north carolina",]


NC_pred_df = cbind(pred_loc, mean_pred)
pred_raster = rasterFromXYZ(NC_pred_df)
cropped_pred = mask(pred_raster, st_as_sf(NC_state$geom))
plot(cropped_pred)





#### Fit with Link ####

# fit logistic model for NA:
# first create new column, NA=0 or Obs=1
# write log likelihood, gradient (score), hessian 
# latent field to predict via GP Vecchia?  Gradient descent
# Consider alternative kernels

# choose r, obs threshold delta

test_df = PFAS[PFAS$Abbreviation=="PFNA",]
test_df = test_df[,c('Conc_ppt', "nc_public_water_Latitude","nc_public_water_Longitude")]
locs = cbind(test_df$nc_public_water_Latitude, test_df$nc_public_water_Longitude)
test_df = test_df[!duplicated(locs),]

llh_gamma_cens = function(y,z,a,delta, K){
  total_llh = 0
  for(i in 1:length(y)){
    if(!is.na(y[i])){
      point_llh = a*log(a)-z[i]*a-log(gamma(a))+(a-1)*log(y[i])-a*y[i]*exp(-z[i])
      total_llh = total_llh + point_llh
    }else{
      cumul_dense =  pgamma(delta, shape=a, scale = 1/a*exp(z[i]))
      total_llh = total_llh + cumul_dense
    }
  }
  total_llh - z%*%solve(K, z)/2 - length(z)/2*log(2*pi)-1/2* log(det(K))
  return(total_llh)
}

gradient_llh_gamma_cens = function(y,z,a,delta, K){
  grad_z = rep(0, length(z))
  for(i in 1:length(grad_z)){
    if(is.na(y[i])){
      F_aB = pgamma(delta, shape=a, scale =  1/a*exp(z[i]))
      F_a1B = pgamma(delta, shape=a+1, scale = 1/a*exp(z[i]))
      grad_z[i] = a*(F_a1B-F_aB)/F_aB
    }else{
      grad_z[i] = a*(exp(-z[i])*y[i]-1)
    }
  }
  gp_component = -solve(K, z)
  return(grad_z + gp_component)
}


hessian_llh_gamma_cens=  function(y,z,a,delta, K){
  hess_z = rep(0, length(z))
  for(i in 1:length(z)){
    if(is.na(y[i])){
      F_aB = pgamma(delta, shape=a, scale =  1/a*exp(z[i]))
      F_a1B = pgamma(delta, shape=a+1, scale = 1/a*exp(z[i]))
      F_a2B = pgamma(delta, shape=a+2, scale = 1/a*exp(z[i]))
      F_aB_prime = a*(F_a1B-F_aB)
      F_a1B_prime = (a+1)*(F_a2B-F_a1B)
      quotient = (F_a1B_prime*F_aB-F_aB_prime*F_a1B)/F_aB_prime^2
      #negate for consistency with real case: negative folded into update
      hess_z[i]=-(a*quotient-1)
    }else{
      hess_z[i] = a*exp(-z[i])*y[i]
    }
  }
  gp_component = solve(K)
  return(diag(hess_z) + gp_component)
}



y_obs = test_df$Conc_ppt
#y_obs[is.na(y_obs)]=1
z_est = rep(-1, length(y_obs))
locs = cbind(test_df$nc_public_water_Latitude, test_df$nc_public_water_Longitude)
K = Matern(rdist(locs), range = 2, smoothness = 1.5)
K_inv = solve(K)

for(i in 1:15){
  hess_mat = hessian_llh_gamma_cens(y_obs, z_est, a=2, delta = 1, K=K)
  grad_vec = gradient_llh_gamma_cens(y_obs, z_est, a=2, delta = 1, K=K)
  z_est = z_est + solve(hess_mat,grad_vec)
  #z_est = z_est - 1e-6 *grad_vec
  print(llh_gamma_cens(y_obs, z_est, a=2, delta = 1, K=K))
  print(z_est[1:5])
}


x_pred = seq(min(PFAS$nc_public_water_Longitude), max(PFAS$nc_public_water_Longitude), length.out=100)
y_pred = seq(min(PFAS$nc_public_water_Latitude), max(PFAS$nc_public_water_Latitude), length.out=100)
pred_loc = expand.grid(x_pred, y_pred)

K_pn = Matern(rdist(pred_loc, locs), range = 6, smoothness = 1)
mean_pred = K_pn%*%solve(K, z_est)

NC_state = states[states$ID=="north carolina",]


NC_pred_df = cbind(pred_loc, mean_pred)
pred_raster = rasterFromXYZ(NC_pred_df)
cropped_pred = mask(pred_raster, st_as_sf(NC_state$geom))
plot(cropped_pred)





#### Testing  ####
locs=seq(0, 10, by=.2)
Ktest = Matern(rdist(locs), range = 2, smoothness = 1.5)
z_latent = t(chol(Ktest))%*%rnorm(length(locs))
y = rgamma(length(locs), shape = 2, scale = 1/2*exp(z_latent))
par(mfrow=c(1,2))
plot(locs, z_latent, main="Latent GP")
y[y<.5] = NA
plot(locs, y, main="Observed + Censored")
points(locs[is.na(y)], rep(.5, sum(is.na(y))),col=2)
K=Ktest
z_est = rep(0, length(locs))
for(i in 1:25){
  hess_mat = hessian_llh_gamma_cens(y, z_est, a=2, delta = .5, K=K)
  grad_vec = gradient_llh_gamma_cens(y, z_est, a=2, delta = .5, K=K)
  z_est = z_est + solve(hess_mat,grad_vec)
}
plot(locs, z_latent, main = "Latent GP, True + Estim")
points(locs, z_est, col=3)
plot(locs, y, main = "Observed + Predicted")
points(locs[is.na(y)], rep(.5, sum(is.na(y))),col=2)
lines(locs, exp(z_est), col=3)