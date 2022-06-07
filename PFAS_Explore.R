

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
library(fields)



##### Plot over the State ####
PFAS <- read.csv("/Users/zilberds/Desktop/PFAS/nc_pfas_joined.csv")
states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))


counties <- st_as_sf(map("county", plot = FALSE, fill = TRUE))
counties <- subset(counties, grepl("north carolina", counties$ID))


ggplot() +
  geom_sf(data = counties, fill = NA, color = gray(.5))+
  geom_point(data = PFAS,aes(nc_public_water_Longitude,
                             nc_public_water_Latitude))




#### Fit with Link ####

# fit logistic model for NA:
# first create new column, NA=0 or Obs=1
# write log likelihood, gradient (score), hessian 
# latent field to predict via GP Vecchia?  Gradient descent
# Consider alternative kernels

# Prior for beta:  Jeffry, solve(t(X)%*%X)

# choose r, obs threshold delta
llh_gamma_cens = function(y,z,a,delta, K, eta, beta, V){
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
  
  GP_llh =  - t(eta)%*%solve(K, eta)/2 - length(eta)/2*log(2*pi)-1/2* log(det(K))
  beta_llh = 0
  if(!is.null(dim(V))){
    beta_llh = - t(beta)%*%solve(V, beta)/2 - length(eta)/2*log(2*pi)-1/2* log(det(V))
  }
  return(total_llh+beta_llh+GP_llh)
}

gradient_llh_gamma_cens = function(y,z,a,delta, K, eta = NA){
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
  if(!all(is.na(eta))) gp_component = -solve(K, eta)
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
      quotient = (F_a1B_prime*F_aB-F_aB_prime^2)/F_aB^2
      #negate for consistency with real case: negative folded into update
      hess_z[i]=-(a*quotient-1)
    }else{
      hess_z[i] = a*exp(-z[i])*y[i]
    }
  }
  gp_component = solve(K)
  return(diag(hess_z) + gp_component)
}




# setup constants
NC_state = states[states$ID=="north carolina",]
x_axis_pred = seq(min(PFAS$nc_public_water_Longitude)-.5, max(PFAS$nc_public_water_Longitude), length.out=100)
y_axis_pred = seq(min(PFAS$nc_public_water_Latitude), max(PFAS$nc_public_water_Latitude)+.1, length.out=100)
pred_loc = expand.grid(x_axis_pred, y_axis_pred)


# predict field 
test_df = PFAS[PFAS$Abbreviation=="PFMOAA",]
test_df = test_df[,c('Conc_ppt', "nc_public_water_Latitude",
                     "nc_public_water_Longitude",
                     "nc_public_water_Depth",
                     "nc_public_water_RECHRATE")]

locs = cbind(test_df$nc_public_water_Longitude,test_df$nc_public_water_Latitude)
test_df = test_df[!duplicated(locs),]
locs = locs[!duplicated(locs),]
y_obs = test_df$Conc_ppt
X = cbind(test_df$nc_public_water_Depth, test_df$nc_public_water_RECHRATE)
X[which(is.na(X[,1])),1]=0
#y_obs[is.na(y_obs)]=1
eta_est = rep(-1, length(y_obs))
B_est = c(0,0)
V = diag(rep(1, length(B_est)))
K = Matern(rdist(locs), range = 1, smoothness = 1.5)
K_inv = solve(K)
z_est = X%*%B_est + eta_est
for(i in 1:50){
  hess_mat_z = hessian_llh_gamma_cens(y_obs, z_est, a=2, delta = 1, 
                                      K=K)
  grad_vec_z = gradient_llh_gamma_cens(y_obs, z_est, a=2, delta = 1, 
                                       K=K, eta = eta_est)
  
  if(!is.null(dim(V))){
    hess_mat_beta = t(X)%*%solve(K)%*%X+solve(V)
    grad_vec_beta = t(X)%*%solve(K, eta_est) -solve(V, B_est)
    B_est = B_est +solve(hess_mat_beta,grad_vec_beta)
  }
  
  z_est = z_est + solve(hess_mat_z,grad_vec_z)
  eta_est = z_est - X%*%B_est 
  
  print(norm(solve(hess_mat_z,grad_vec_z)))
  print(eta_est[1:5])
  print(B_est)
}

##LR test:  2(L(Beta)-L(beta = 0)) = 2*log p(beta) ~ Chisq 2dof
1-pchisq(t(B_est)%*%solve(V, B_est), df = 2)
# beta: not significant!

K_pn = Matern(rdist(pred_loc, locs), range = 1, smoothness = 1.5)
mean_pred = K_pn%*%solve(K, z_est)
#mean_pred = K_pn%*%solve(K,  X%*%B_est)



NC_pred_df = cbind(pred_loc, "value"=exp(mean_pred))
pred_raster = rasterFromXYZ(NC_pred_df)
cropped_pred = mask(pred_raster, st_as_sf(NC_state$geom))
plot(cropped_pred, main = "Predicted PFNA, <1 = Censored ")

crop_pred_df = as.data.frame(cropped_pred, xy=T)
crop_pred_df = crop_pred_df[complete.cases(crop_pred_df),]
ggplot() +  
  geom_tile(data=crop_pred_df, aes(x=x, y=y, fill=value), alpha=0.8) +
  geom_sf(data = counties, fill = NA, color = gray(.5))+
  geom_point(data = PFAS,aes(nc_public_water_Longitude,
                             nc_public_water_Latitude), size=1)





#### River modeling ####
# goal:  given obs in the data set, make predictions along rivers
#  river predictions can be extended to non-river predictions with additive kernel

library(SSN)
library(broom)

file_dir = "/Users/zilberds/Desktop/PFAS/river_shapefiles/Flowline_SA03N_NSI/"
shp_file_dir = "/Users/zilberds/Desktop/PFAS/river_shapefiles/Flowline_SA03N_NSI/Flowline_SA03N_NSI.shp"
riv_shape <- readOGR( dsn=file_dir,
                      verbose=TRUE)
riv_shape@data
riv_shape@lines
summary(riv_shape)
length(riv_shape)
class(riv_shape)
crs(riv_shape)
extent(riv_shape)

plot(riv_shape[1:5000,])
riv_shape[1:5,]



riv_shp_tidy = tidy(riv_shape[100:3000,])
ggplot() +
  geom_path(data = riv_shp_tidy, aes( x = long, y = lat, group = group)) 

riv_mask= mask(riv_shp_tidy, st_as_sf(NC_state$geom))
# Projected coordinate system: NAD_1983_Albers
# TODO: try to project onto lat long or vice versa
# mask and overlay on NC
# 


## Check public water sources shape
file_dir_pubwater = "/Users/zilberds/Desktop/PFAS/Public_Water_Supply_Water_Sources/"
pubwat_shape <- readOGR( dsn=file_dir_pubwater,
                         verbose=TRUE)
crs(pubwat_shape)
plot(pubwat_shape)

#### Exploratory data distirbutions#####
######Gamma vs NB vs Exp dist fits ######
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


#### Verification of NR  ####
locs=seq(0, 50, by=.2)
Ktest = Matern(rdist(locs), range = 2, smoothness = 2.5)
eta_latent = t(chol(Ktest))%*%rnorm(length(locs))
X = sin(locs)
beta_coef = 1.3
z_latent=eta_latent+beta_coef*X
# generate data
y = rgamma(length(locs), shape = 2, scale = 1/2*exp(z_latent))
y_obs=y
par(mfrow=c(1,2))
plot(locs, z_latent, main="Latent GP")
# censor small values
y[y<.5] = NA
plot(locs, y, main="Observed + Censored")
points(locs[is.na(y)], rep(.5, sum(is.na(y))),col=2)
K=Ktest
B_est = c(0)
V = diag(rep(1, length(B_est)))
eta_est = rep(-1, length(y_obs))
z_est = X%*%B_est + eta_est
# Newton Raphson iterations to find the MAP
z_est = rep(0, length(locs))
for(i in 1:25){
  # hess_mat = hessian_llh_gamma_cens(y, z_est, a=2, delta = .5, K=K)
  # grad_vec = gradient_llh_gamma_cens(y, z_est, a=2, delta = .5, K=K)
  # print( norm(solve(hess_mat,grad_vec)))
  # z_est = z_est + solve(hess_mat,grad_vec)


  hess_mat_z = hessian_llh_gamma_cens(y_obs, z_est, a=2, delta = 1, 
                                      K=K)
  grad_vec_z = gradient_llh_gamma_cens(y_obs, z_est, a=2, delta = 1, 
                                       K=K, eta = eta_est)
  
  if(!is.null(dim(V))){
    hess_mat_beta = t(X)%*%solve(K)%*%X+solve(V)
    grad_vec_beta = t(X)%*%solve(K, eta_est) -solve(V, B_est)
    B_est = B_est +solve(hess_mat_beta,grad_vec_beta)
  }
  
  z_est = z_est + solve(hess_mat_z,grad_vec_z)
  eta_est = z_est - X%*%B_est 
  
  print( norm(solve(hess_mat_z,grad_vec_z)))
}
# plot the latent field, compare to truth
plot(locs, z_latent, main = "Latent GP, True + Estim")
points(locs, z_est, col=3)
# plot the observed field, compare to censored
plot(locs, y, main = "Observed + Predicted")
points(locs[is.na(y)], rep(.5, sum(is.na(y))),col=2)
lines(locs, exp(z_est), col=3)

