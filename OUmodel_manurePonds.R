#####################################################################
######### define pdf and random generator for univariate OU #########
#####################################################################

rOU <- function(n, x0 = 0, t = 1, sigma2 = 1, alpha = 0.1, theta = 0){
  mean <- theta + (x0-theta)*exp(-1*alpha*t)
  var <- sigma2 / (2*alpha) * (1-exp(-2*alpha*t))
  rnorm(n = n, mean = mean, sd = sqrt(var))
}

library(latex2exp)

nlines <- 100 #number of realizations of the process to plot
x0 <- 1 #starting state of the process
theta <- 10 #location of the optimum
sigma2 = 1 #brownian motion rate parameter
alpha = .75 #strength of the rubber band / returning force
t <- 0:4000/500 #discrete time steps over which to simulate the process
x <- matrix(0,nrow = nlines, ncol = length(t)); x[,1] <- x0 #initialize the matrix of realizations
for(j in 1:nlines){
  for(i in 2:length(t)){
    x[j,i] <- rOU(n=1, x0=x[j,i-1], t = t[i]-t[i-1], sigma2=sigma2, alpha=alpha, theta=theta)
  }
}
plot(t,x[1,], type = "l", ylim = range(x), col=rgb(0,0,0,0.2), ylab = "x")
abline(h = theta, col=2, lwd = 2) #plot the optimum as a horizontal line
for(j in 2:nlines){lines(t,x[j,],col=rgb(0,0,0,0.2))}
title(main = TeX(paste0("Realizations from OU: $\\alpha = ", alpha, ",\\, \\theta = ", theta, ",\\, \\sigma^2 = ", sigma2, "$")))

rOU(10)

time <- (1:1000)/500


dOU <- function(x0, x1, t, sigma2, alpha, theta, log = T){
  mean <- theta + (x0-theta)*exp(-1*alpha*t)
  var <- sigma2 / (2*alpha) * (1-exp(-2*alpha*t))
  dnorm(x = x1, mean = mean, sd = sqrt(var), log = log)
}

dOU(x0 = 0,x1 = rOU(10), t = 1, sigma2 = 1, alpha = 0.1, theta = 0)

#############################################
######### simulate model parameters #########
#############################################

nPonds <- 70
#how long before the present did ponds start
samplePeriodPerPond <- sample(100:1000, size = nPonds, replace = T) 
#how often do ponds sample N per day? (i.e. with what rate?)
samplingRatePerPond <- runif(min = 0.01, max = 0.03, n = nPonds)
#let's assume sampling is a poisson process
sampleSojournsPerPond <- as.list(rep(0,nPonds))
for(i in 1:nPonds){
  time <- rexp(samplingRatePerPond[i], n = 100)
  sampleSojournsPerPond[[i]] <-  time[1:which(cumsum(time) > samplePeriodPerPond[i])[1]]
}

# #simulate pond-specific optima
# perPondTheta <- rnorm(nPonds, mean = 10, sd = 1)
# #simulate pond-specific returning forces, needs to be >0
# perPondAlpha <- rgamma(nPonds, 0.1, 0.1)
# #simulate pond-specific magnitude of noise, needs to be >0
# perPondSigma2 <- rgamma(nPonds, 0.1, 0.1)

#or, alternatively, simulate from the prior specified in the hierarchical model below
thetaMU <- rnorm(1,10,2)
thetaSigma <- rexp(1,1)
perPondTheta <- rnorm(nPonds, thetaMU, thetaSigma) 

alphaA <- rexp(1,1)
alphaB <- rexp(1,1)
perPondAlpha <- rgamma(nPonds, shape = alphaA, rate = alphaB)

sigmaBMA <- rexp(1,1)
sigmaBMB <- rexp(1,1)
perPondSigma2 <- rgamma(nPonds, shape = sigmaBMA, rate = sigmaBMB)

#################################
######### simulate data #########
#################################

#this collapses variability throughout the pond and through time, since temporally adjascent samples will also be spatially adjascent
#assuming observations are obtained from the same location

d <- as.data.frame(matrix(0, 1e5, 3)) #initialize empty data frame
colnames(d) <- c("pond", "day", "nitrogen")

foo <- 1 #index counter
for(i in 1:nPonds){
  sojourns <- sampleSojournsPerPond[[i]]
  ndata <- length(sojourns) + 1
  theta <- perPondTheta[i]
  alpha <- perPondAlpha[i]
  sigma2 <- perPondSigma2[i]
  Nitr <- rep(0, ndata)
  #generate starting value assuming process is at stationarity 
  #just simulate after a million days with arbitrary x0
  Nitr[1] <- rOU(n = 1, x0 = 0, t = 1e6, sigma2 = sigma2, alpha = alpha, theta = theta)
  for(j in 2:ndata){
    Nitr[j] <- rOU(n = 1, x0 = Nitr[j-1], t = sojourns[j-1], sigma2 = sigma2, alpha = alpha, theta = theta)
  }
  d[foo:(foo + ndata - 1),] <- cbind(i, c(0, round(cumsum(sojourns))), Nitr)
  foo <- foo + ndata
}
#trim excess
d <- d[d$pond != 0,]
d$pond <- as.integer(d$pond)

#sample random starting dates
startingDates <- sample(seq(as.Date('2012/01/01'), as.Date('2014/01/01'), by="day"), nPonds)
d$day <- d$day + startingDates[d$pond]

# d is what a normal data table would look like
# but we need to process it into sequential (i.e. independent) pairs of observations
# (we effectively simulated these from the start, but for clarity it is better to produce it post-hoc)
n_ponds <- length(unique(d$pond))
obs_per_pond <- as.vector(table(d$pond))
data <- as.data.frame(matrix(0, sum(as.vector(obs_per_pond-1)), 5))
colnames(data) <- c("pond", "interval", "nitroStart", "nitroEnd", "diffNitro")
foo <- 1
for(i in 1:n_ponds){
  pond <- i
  nobs <- obs_per_pond[i]
  pd <- d[d$pond == pond,]
  pd <- pd[order(pd$day),] #order by day
  diffNitro <- sapply(1:(nobs-1), function(x) pd$nitrogen[x+1] - pd$nitrogen[x])
  diffDays <- sapply(1:(nobs-1), function(x) pd$day[x+1] - pd$day[x])
  data[foo:(foo+nobs-2),] <- cbind(pond, diffDays, pd$nitrogen[1:(nobs-1)], pd$nitrogen[2:(nobs)], diffNitro)
  foo <- foo + nobs - 1
}
data$nitroEnd - data$nitroStart  - data$diffNitro #quick sanity check

#####################################
######### fit model in Stan #########
#####################################

#we're gonna use RStan later, but imo it's a lot easier to quickly read models specified in the "Rethinking" syntax
library(rethinking)

data <- data[data$interval != 0,] #model breaks if interval is 0, even though I rounded and 0 is really <0.5
data$pond <- as.integer(data$pond)

#treat all ponds as the same; i.e. ignore variation in process between ponds, disregard hierarchical structure
m0 <- map2stan(
  alist(
    nitroEnd ~ dnorm(mu, sigma),
    mu <- theta + (nitroStart - theta) * exp(-1 * alpha * interval),
    sigma <- 0 + sqrt(sigma2),
    sigma2 <- sigmaBM / (2 * alpha) * (1 - exp(-2 * alpha * interval)),
    
    #these can all be informed by talking to farmers
    theta ~ dnorm(theta_MU, theta_SIGMA), 
    theta_MU ~ dnorm(10,2),
    theta_SIGMA ~ dexp(1),
    
    alpha ~ dgamma(alpha_A, alpha_B),
    alpha_A ~ dexp(1),
    alpha_B ~ dexp(1),
    
    sigmaBM ~ dgamma(sigmaBM_A, sigmaBM_B),
    sigmaBM_A ~ dexp(1),
    sigmaBM_B ~ dexp(1)
    
  ) ,
  data= data,
  iter = 5e3, warmup = 5e3, chains = 1, cores = 1, sample = T,
  start = list(alpha_A=0.1, alpha_B=0.1, sigmaBM_A = 0.1, sigmaBM_B = 0.1, theta_MU = 10, theta_SIGMA = 1)
)
cat(m0$model)

summary(m0)
plot(m0)
pairs(m0)

m1 <- map2stan(
  alist(
    nitroEnd ~ dnorm(mu, sigma),
    mu <- theta[pond] + ((nitroStart - theta[pond]) * exp(-1 * alpha[pond] * interval)),
    sigma <- 0 + sqrt(sigma2),
    sigma2 <- (sigmaBM[pond] * (1 - exp(-2 * alpha[pond] * interval))) / (2 * alpha[pond]),
    
    #these can all be informed by talking to farmers
    theta[pond] ~ dnorm(theta_MU, theta_SIGMA), 
    theta_MU ~ dnorm(10,2),
    theta_SIGMA ~ dexp(1),
    
    alpha[pond] ~ dgamma(alpha_A, alpha_B),
    alpha_A ~ dexp(1),
    alpha_B ~ dexp(1),
    
    sigmaBM[pond] ~ dgamma(sigmaBM_A, sigmaBM_B),
    sigmaBM_A ~ dexp(1),
    sigmaBM_B ~ dexp(1)

  ) ,
  data= data, constraints=list(alpha_A="lower=0", alpha_B="lower=0", sigmaBM_A="lower=0", sigmaBM_B="lower=0", theta_SIGMA="lower=0", theta="lower=0", alpha="lower=0", sigmaBM="lower=0"),
  iter = 2e3, warmup = 1e3, chains = 1, cores = 1, sample = T,
  start = list(alpha_A=0.1, alpha_B=0.1, sigmaBM_A = 0.1, sigmaBM_B = 0.1, theta_MU = 10, theta_SIGMA = 1)
)

#hmm it looks like this fails to initialize -- maybe it's struggling to respect the lower bounds?
#let's just implement the model in raw Stan and pass that to RStan

m1t <- 
"data{
  int<lower=1> N;
  int<lower=1> N_pond;
  real nitroEnd[N];
  int pond[N];
  real interval[N];
  real nitroStart[N];
}
parameters{
  real<lower=0> alpha_A;
  real<lower=0> alpha_B;
  real<lower=0> sigmaBM_A;
  real<lower=0> sigmaBM_B;
  real theta_MU;
  real<lower=0> theta_SIGMA;
  vector<lower=0>[N_pond] theta;
  vector<lower=0>[N_pond] alpha;
  vector<lower=0>[N_pond] sigmaBM;
}
model{
  vector[N] sigma2;
  vector[N] sigma;
  vector[N] mu;
  sigmaBM_B ~ exponential( 1 );
  sigmaBM_A ~ exponential( 1 );
  for ( j in 1:N_pond ) sigmaBM[j] ~ gamma( sigmaBM_A , sigmaBM_B );
  alpha_B ~ exponential( 1 );
  alpha_A ~ exponential( 1 );
  for ( j in 1:N_pond ) alpha[j] ~ gamma( alpha_A , alpha_B );
  theta_SIGMA ~ cauchy( 0 , 1 );
  theta_MU ~ normal( 10 , 2 );
  theta ~ normal( theta_MU , theta_SIGMA );
  for ( i in 1:N ) {
    sigma2[i] = (sigmaBM[pond[i]] * (1 - exp(-2 * alpha[pond[i]] * interval[i])))/(2 *      alpha[pond[i]]);
  }
  for ( i in 1:N ) {
    sigma[i] = 0 + sqrt(sigma2[i]);
  }
  for ( i in 1:N ) {
    mu[i] = theta[pond[i]] + ((nitroStart[i] - theta[pond[i]]) * exp(-1 * alpha[pond[i]] *      interval[i]));
  }
  nitroEnd ~ normal( mu , sigma );
}
generated quantities{
  vector[N] sigma2;
  vector[N] sigma;
  vector[N] mu;
  real dev;
  dev = 0;
  for ( i in 1:N ) {
    sigma2[i] = (sigmaBM[pond[i]] * (1 - exp(-2 * alpha[pond[i]] * interval[i])))/(2 *      alpha[pond[i]]);
  }
  for ( i in 1:N ) {
    sigma[i] = 0 + sqrt(sigma2[i]);
  }
  for ( i in 1:N ) {
    mu[i] = theta[pond[i]] + ((nitroStart[i] - theta[pond[i]]) * exp(-1 * alpha[pond[i]] *      interval[i]));
  }
  dev = dev + (-2)*normal_lpdf( nitroEnd | mu , sigma );
}"

m1 <- stan(model_code = m1t, verbose = T, init_r = 0.01, chains = 1, iter = 5e5, control = list(adapt_delta = 0.9925), thin = 5,
           data=list(pond = data$pond, nitroEnd = data$nitroEnd, nitroStart = data$nitroStart, 
           interval = data$interval, N = length(data[,1]), N_pond = length(unique(data$pond))))

save(m1, file = "OU_manurePonds_hierarchicalModel")
load(file = "OU_manurePonds_hierarchicalModel")

summary(m1)
m1s <- extract.samples(m1)
infThetas <- m1s$theta
infAlphas <- m1s$alpha
infSigmaBM <- m1s$sigmaBM

#let's look at how well we retrieve pond-specific parameters
dev.off()
plot(apply(infThetas, 2, mean), perPondTheta, xlab = "posterior mean", ylab = "true value"); abline(0,1,lwd=2,col=2); title(("theta"))
plot(apply(infAlphas, 2, mean), perPondAlpha, xlab = "posterior mean", ylab = "true value"); abline(0,1,lwd=2,col=2); title(("alpha"))
plot(apply(infSigmaBM, 2, mean), perPondSigma2, xlab = "posterior mean", ylab = "true value"); abline(0,1,lwd=2,col=2); title(("sigma2_BM"))

#often in OU models the 'returning force' and 'brownian noise' parameters are non-identifiable, 
#but with so much data here we find ourselves in possession of decent ability to tease them apart 
plot(m1s$sigmaBM_A/m1s$sigmaBM_B, m1s$alpha_A/m1s$alpha_B)
cor(m1s$sigmaBM_A/m1s$sigmaBM_B, m1s$alpha_A/m1s$alpha_B)

#what about hyperparameters?
par(mfrow = c(3,2))
hist(m1s$theta_MU, breaks = 100); abline(v = thetaMU, col = 2, lwd = 2)
hist(m1s$theta_SIGMA, breaks = 100); abline(v = thetaSigma, col = 2, lwd = 2); legend("true value", x = "topright", lwd = 2, col = 2)
hist(m1s$alpha_A, breaks = 100); abline(v = alphaA, col = 2, lwd = 2)
hist(m1s$alpha_B, breaks = 100); abline(v = alphaB, col = 2, lwd = 2)
hist(m1s$sigmaBM_A, breaks = 100); abline(v = sigmaBMA, col = 2, lwd = 2)
hist(m1s$sigmaBM_B, breaks = 100); abline(v = sigmaBMB, col = 2, lwd = 2)

#if we wanted to do guess future values, we could could generate posterior predictive distributions easily
#for example, say the owner of pond #37 tells us they took a nitrogen measurement of 10.1 two weeks ago
#our prediction for its current value would be:
pond_id <- 17
n_days <- 14
starting_value <- 11.1
postTheta <- m1s$theta[,pond_id]; postTheta <- postTheta[round(seq(1, length(postTheta), length.out = 1000))]
postAlpha <- m1s$alpha[,pond_id]; postAlpha <- postAlpha[round(seq(1, length(postAlpha), length.out = 1000))]
postSigma2BM <- m1s$sigmaBM[,pond_id]; postSigma2BM <- postSigma2BM[round(seq(1, length(postSigma2BM), length.out = 1000))]

samp_per_samp <- 100
preds <- as.vector(sapply(1:length(postTheta), function(x) rOU(n = samp_per_samp, x0 = starting_value, t = n_days, sigma2 = postSigma2BM[x], alpha = postAlpha[x], theta = postTheta[x])))
par(mfrow = c(1,1))
hist(preds, breaks = 100, xlim = c(5,15), xlab = "predicted nitrogen value", 
     main = paste0("starting value = ", starting_value, ", pond id = ", pond_id, ", time since last measurement = ", round(n_days, 3), " days")); 
abline(v = mean(preds), col = 2, lwd = 2); legend("mean prediction", x = "topright", lwd = 2, col = 2)

#conversely, if he measured that pond an hour ago, we'd tell him the current value is
n_days <- 2 / 24 
preds <- as.vector(sapply(1:length(postTheta), function(x) rOU(n = samp_per_samp, x0 = starting_value, t = n_days, sigma2 = postSigma2BM[x], alpha = postAlpha[x], theta = postTheta[x])))
hist(preds, breaks = 100, xlim = c(5,15), xlab = "predicted nitrogen value", freq = FALSE,
     main = paste0("starting value = ", starting_value, ", pond_id = ", pond_id, ", time since last measurement = ", round(n_days, 3), " days")); 
abline(v = mean(preds), col = 2, lwd = 2); legend("mean prediction", x = "topright", lwd = 2, col = 2)
#much greater certainty!

#conversely, if he measured that pond 2 days ago, we'd tell him the current value is
n_days <- 0.75
preds <- as.vector(sapply(1:length(postTheta), function(x) rOU(n = samp_per_samp, x0 = starting_value, t = n_days, sigma2 = postSigma2BM[x], alpha = postAlpha[x], theta = postTheta[x])))
hist(preds, breaks = 100, xlim = c(5,15), xlab = "predicted nitrogen value", freq = FALSE,
     main = paste0("starting value = ", starting_value, ", pond_id = ", pond_id, ", time since last measurement = ", round(n_days, 3), " days")); 
abline(v = mean(preds), col = 2, lwd = 2); legend("mean prediction", x = "topright", lwd = 2, col = 2)

#conversely, if he measured that pond 2 days ago, we'd tell him the current value is
n_days <- 7
preds <- as.vector(sapply(1:length(postTheta), function(x) rOU(n = samp_per_samp, x0 = starting_value, t = n_days, sigma2 = postSigma2BM[x], alpha = postAlpha[x], theta = postTheta[x])))
hist(preds, breaks = 100, xlim = c(5,15), xlab = "predicted nitrogen value", freq = FALSE,
     main = paste0("starting value = ", starting_value, ", pond_id = ", pond_id, ", time since last measurement = ", round(n_days, 3), " days")); 
abline(v = mean(preds), col = 2, lwd = 2); legend("mean prediction", x = "topright", lwd = 2, col = 2)

#conversely, if he measured that pond 2,000,000 days ago, we'd tell him the current value is
n_days <- 2e5
preds <- as.vector(sapply(1:length(postTheta), function(x) rOU(n = samp_per_samp, x0 = starting_value, t = n_days, sigma2 = postSigma2BM[x], alpha = postAlpha[x], theta = postTheta[x])))
hist(preds, breaks = 100, xlim = c(5,15), xlab = "predicted nitrogen value", freq = FALSE,
     main = paste0("starting value = ", starting_value, ", pond_id = ", pond_id, ", time since last measurement = ", round(n_days, 3), " days")); 
abline(v = mean(preds), col = 2, lwd = 2); legend("mean prediction", x = "topright", lwd = 2, col = 2)
#it has returned to the limiting distribution of the process


#(real data may have more or less certainty, relative to how much noise is actually in the measured dataset)
