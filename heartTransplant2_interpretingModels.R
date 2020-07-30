#### INTERPRETING AND COMPARING THE SUCCESFULLY FITTED MODELS ####
setwd("/Volumes/1TB/heart_transplant/")
library(rethinking)
load("m8_data_noP9")
load(file = "m8_noncentered_noP9") #WAIC 100959.5, gets 100% of WAIC weight
# load(file = "m10_noncentered_noP9") #WAIC 143984.8
# load(file = "m12_noncentered_noP9") #WAIC 103238.7
load(file = "m8_pat_key_noP9")
pat_key <- pat_key[order(pat_key[,2]),]
load(file = "m8_prot_key_noP9")
prot_key <- prot_key[order(prot_key[,2]),]

# compare(m8, m12)

m8s <- extract.samples(m8, n = 5000)
precis_m8 <- precis(m8)
m8@formula

#starting intercept values
dens(exp(m8s$A0b + m8s$A0j[,1]*m8s$sigA0j), add = F)
for(i in 2:ncol(m8s$A0j)){
  dens(exp(m8s$A0b + m8s$A0j[,i]*m8s$sigA0j), add = T, col = rgb(0, 0, 0, 0.15))
}
quantile(sapply(1:ncol(m8s$A0j), function(prot) mean(exp(m8s$A0b + m8s$A0j[,prot]*m8s$sigA0j))), probs = 0:100/100)
intercept_prot_inds_decr <- order(sapply(1:ncol(m8s$A0j), function(prot) mean(exp(m8s$A0b + m8s$A0j[,prot]*m8s$sigA0j))), decreasing = T)
highest_intercepts_order_names <- prot_key[intercept_prot_inds_decr,1]
highest_intercepts_order_means <- sort(sapply(1:ncol(m8s$A0j), function(prot) mean(exp(m8s$A0b + m8s$A0j[,prot]*m8s$sigA0j))), decreasing = T)
highest_intercepts_order_means_in_nonrejectors <- sort(sapply(1:ncol(m8s$A0j), function(prot) mean(exp(m8s$A0b + m8s$A0j[,prot]*m8s$sigA0j))), decreasing = T)
highest_intercepts_order_means_in_rejectors <- (sapply(1:ncol(m8s$A0j), function(prot) mean(exp(m8s$A0b + m8s$A0j[,prot]*m8s$sigA0j + m8s$a_rej))))[intercept_prot_inds_decr]
pre_transplant_intercept <- cbind(protein_name = as.character(highest_intercepts_order_names), 
                                  mean_spectral_count_in_nonrejectors = round(highest_intercepts_order_means_in_nonrejectors, 3),
                                  mean_spectral_count_in_rejectors = round(highest_intercepts_order_means_in_rejectors, 3))
write.csv(pre_transplant_intercept, file = "pre-transplant_intercept_values_m8_noP9.csv")

m8@formula
#relative change through time per unit day values, post-transplant / pre-rejection
dens((exp((m8s$bR0b*2 + m8s$bR0j[,1]*m8s$sigbR0j)*0.03)-1)*100, add = F)
for(i in 2:ncol(m8s$A0j)){
  dens((exp((m8s$bR0b*2 + m8s$bR0j[,i]*m8s$sigbR0j)*0.03)-1)*100, add = T)
  
}
hist(sapply(1:ncol(m8s$bR0j), function(prot) mean((exp((m8s$bR0b*2 + m8s$bR0j[,prot]*m8s$sigbR0j)*0.03)-1)*100)))
highest_slope_order_names <- prot_key[order(sapply(1:ncol(m8s$bR0j), function(prot) mean((exp((m8s$bR0b*2 + m8s$bR0j[,prot]*m8s$sigbR0j)*0.03)-1)*100)), decreasing = T),1]

highest_slope_order_means_in_nonrejectors <- sort(sapply(1:ncol(m8s$bR0j), function(prot) mean((exp((m8s$bR0b*2 + m8s$bR0j[,prot]*m8s$sigbR0j)*0.03)-1)*100)), decreasing = T)
highest_slope_order_probsGreatZero_in_nonrejectors <- sapply(1:ncol(m8s$bR0j), function(prot) sum(((exp((m8s$bR0b*2 + m8s$bR0j[,prot]*m8s$sigbR0j)*0.03)-1)*100) > 0) / nrow(m8s$bR0j))[
  order(sapply(1:ncol(m8s$bR0j), function(prot) mean((exp((m8s$bR0b*2 + m8s$bR0j[,prot]*m8s$sigbR0j)*0.03)-1)*100)), decreasing = T)]

highest_slope_order_means_in_rejectors <- (sapply(1:ncol(m8s$bR0j), function(prot) mean((exp((m8s$bR0b*2 + m8s$bR0j[,prot]*m8s$sigbR0j + m8s$b_rej)*0.03)-1)*100)))[
  order(sapply(1:ncol(m8s$bR0j), function(prot) mean((exp((m8s$bR0b*2 + m8s$bR0j[,prot]*m8s$sigbR0j)*0.03)-1)*100)), decreasing = T)]
highest_slope_order_probsGreatZero_in_rejectors <- sapply(1:ncol(m8s$bR0j), function(prot) sum(((exp((m8s$bR0b*2 + m8s$bR0j[,prot]*m8s$sigbR0j + m8s$b_rej)*0.03)-1)*100) > 0) / nrow(m8s$bR0j))[
  order(sapply(1:ncol(m8s$bR0j), function(prot) mean((exp((m8s$bR0b*2 + m8s$bR0j[,prot]*m8s$sigbR0j)*0.03)-1)*100)), decreasing = T)]


post_transplant_slope_per_day <- cbind(protein_name = as.character(highest_slope_order_names), 
                                       mean_spectral_count_proportional_percentage_change_per_day_in_nonrejectors = round(highest_slope_order_means_in_nonrejectors, 3),
                                       probability_change_per_day_is_positive_in_nonrejectors = round(highest_slope_order_probsGreatZero_in_nonrejectors, 4), 
                                       mean_spectral_count_proportional_percentage_change_per_day_in_rejectors = round(highest_slope_order_means_in_rejectors, 3),
                                       probability_change_per_day_is_positive_in_rejectors = round(highest_slope_order_probsGreatZero_in_rejectors, 4))
write.csv(post_transplant_slope_per_day, file = "post-transplant_pre-rejection_slope_values_m8_noP9.csv")


m8@formula
#relative change through time per unit day values, post-transplant / post-rejection
dens((exp((m8s$bR1b*2 + m8s$bR1j[,1]*m8s$sigbR1j)*0.03)-1)*100, add = F)
for(i in 2:ncol(m8s$A0j)){
  dens((exp((m8s$bR1b*2 + m8s$bR1j[,i]*m8s$sigbR1j)*0.03)-1)*100, add = T)
  
}
hist(sapply(1:ncol(m8s$bR1j), function(prot) mean((exp((m8s$bR1b*2 + m8s$bR1j[,prot]*m8s$sigbR1j)*0.03)-1)*100)))
highest_slope_order_names <- prot_key[order(sapply(1:ncol(m8s$bR1j), function(prot) mean((exp((m8s$bR1b*2 + m8s$bR1j[,prot]*m8s$sigbR1j)*0.03)-1)*100)), decreasing = T),1]

highest_slope_order_means <- sort(sapply(1:ncol(m8s$bR1j), function(prot) mean((exp((m8s$bR1b*2 + m8s$bR1j[,prot]*m8s$sigbR1j)*0.03)-1)*100)), decreasing = T)
highest_slope_order_probsGreatZero <- sapply(1:ncol(m8s$bR1j), function(prot) sum(((exp((m8s$bR1b*2 + m8s$bR1j[,prot]*m8s$sigbR1j)*0.03)-1)*100) > 0) / nrow(m8s$bR1j))[
  order(sapply(1:ncol(m8s$bR1j), function(prot) mean((exp((m8s$bR1b*2 + m8s$bR1j[,prot]*m8s$sigbR1j)*0.03)-1)*100)), decreasing = T)]


post_rejection_slope_per_day <- cbind(protein_name = as.character(highest_slope_order_names), 
                                       mean_spectral_count_proportional_percentage_change_per_day = round(highest_slope_order_means_in_nonrejectors, 3),
                                       probability_change_per_day_is_positive = round(highest_slope_order_probsGreatZero_in_nonrejectors, 4))
write.csv(post_transplant_slope_per_day, file = "post-transplant_post-rejection_slope_values_m8_noP9.csv")

m8@formula
#discrete drop at rejection event
dens(exp(m8s$A1b + m8s$A1j[,1]*m8s$sigA1j)*100, add = F)
for(i in 2:ncol(m8s$A1j)){
  dens(exp(m8s$A1b + m8s$A1j[,i]*m8s$sigA1j)*100, add = T, col = rgb(0, 0, 0, 0.15))
}
quantile(sapply(1:ncol(m8s$A1j), function(prot) mean(exp(m8s$A1b + m8s$A1j[,prot]*m8s$sigA1j))), probs = 0:100/100)
change_prot_inds_incr <- order(sapply(1:ncol(m8s$A1j), function(prot) mean(exp(m8s$A1b + m8s$A1j[,prot]*m8s$sigA1j))), decreasing = F)
highest_changes_order_names <- prot_key[change_prot_inds_incr,1]
highest_changes_order_means <- 100 * sort(sapply(1:ncol(m8s$A1j), function(prot) mean(exp(m8s$A1b + m8s$A1j[,prot]*m8s$sigA1j))), decreasing = F)
prob_counts_went_down <- (sapply(1:ncol(m8s$A1j), function(prot) sum(exp(m8s$A1b + m8s$A1j[,prot]*m8s$sigA1j) < 1) / nrow(m8s$A1j)))[change_prot_inds_incr]

at_rejection_change <- cbind(protein_name = as.character(highest_changes_order_names), 
                                  mean_percent_spectral_count_change_at_rejection_event = round(highest_changes_order_means, 3),
                                  probability_spectral_count_went_down = round(prob_counts_went_down, 4))
write.csv(at_rejection_change, file = "discrete_rejection_event_percent_change_m8_noP9.csv")
### ## ## ## ## ## ###
## now let's do m11 ## 
### ## ## ## ## ## ### 

load("m8_data_noP9")
load(file = "m8_pat_key_noP9")
pat_key <- pat_key[order(pat_key[,2]),]
load(file = "m8_prot_key_noP9")
prot_key <- prot_key[order(prot_key[,2]),]

load(file = "m11_noncentered_noP9") #WAIC 100959.5, gets 100% of WAIC weight
# load(file = "m10_noncentered_noP9") #WAIC 143984.8
# load(file = "m12_noncentered_noP9") #WAIC 103238.7


# compare(m11, m12)

m11s <- extract.samples(m11, n = 10000)
precis_m11 <- precis(m11)
m11@formula

#starting intercept values
dens(exp(m11s$A0b + m11s$A0j[,1]*m11s$sigA0j), add = F)
for(i in 2:ncol(m11s$A0j)){
  dens(exp(m11s$A0b + m11s$A0j[,i]*m11s$sigA0j), add = T, col = rgb(0, 0, 0, 0.15))
}
quantile(sapply(1:ncol(m11s$A0j), function(prot) mean(exp(m11s$A0b + m11s$A0j[,prot]*m11s$sigA0j))), probs = 0:100/100)
intercept_prot_inds_decr <- order(sapply(1:ncol(m11s$A0j), function(prot) mean(exp(m11s$A0b + m11s$A0j[,prot]*m11s$sigA0j))), decreasing = T)
highest_intercepts_order_names <- prot_key[intercept_prot_inds_decr,1]
highest_intercepts_order_means <- sort(sapply(1:ncol(m11s$A0j), function(prot) mean(exp(m11s$A0b + m11s$A0j[,prot]*m11s$sigA0j))), decreasing = T)
highest_intercepts_order_means_in_nonrejectors <- sort(sapply(1:ncol(m11s$A0j), function(prot) mean(exp(m11s$A0b + m11s$A0j[,prot]*m11s$sigA0j))), decreasing = T)
highest_intercepts_order_means_in_rejectors <- (sapply(1:ncol(m11s$A0j), function(prot) mean(exp(m11s$A0b + m11s$A0j[,prot]*m11s$sigA0j + m11s$a_rej + 
                                                                                                   m11s$A0jr[,prot] * m11s$sigA0jr))))[intercept_prot_inds_decr]
pre_transplant_intercept <- cbind(protein_name = as.character(highest_intercepts_order_names), 
                                  mean_spectral_count_in_nonrejectors = round(highest_intercepts_order_means_in_nonrejectors, 3),
                                  mean_spectral_count_in_rejectors = round(highest_intercepts_order_means_in_rejectors, 3))
write.csv(pre_transplant_intercept, file = "pre-transplant_intercept_values_m11_noP9.csv")

m11@formula
#relative change through time per unit day values, post-transplant / pre-rejection
dens((exp((m11s$bR0b*2 + m11s$bR0j[,1]*m11s$sigbR0j)*0.03)-1)*100, add = F)
for(i in 2:ncol(m11s$A0j)){
  dens((exp((m11s$bR0b*2 + m11s$bR0j[,i]*m11s$sigbR0j)*0.03)-1)*100, add = T)
  
}
hist(sapply(1:ncol(m11s$bR0j), function(prot) mean((exp((m11s$bR0b*2 + m11s$bR0j[,prot]*m11s$sigbR0j)*0.03)-1)*100)))
highest_slope_order_names <- prot_key[order(sapply(1:ncol(m11s$bR0j), function(prot) mean((exp((m11s$bR0b*2 + m11s$bR0j[,prot]*m11s$sigbR0j)*0.03)-1)*100)), decreasing = T),1]

highest_slope_order_means_in_nonrejectors <- sort(sapply(1:ncol(m11s$bR0j), function(prot) mean((exp((m11s$bR0b*2 + m11s$bR0j[,prot]*m11s$sigbR0j)*0.03)-1)*100)), decreasing = T)
highest_slope_order_probsGreatZero_in_nonrejectors <- sapply(1:ncol(m11s$bR0j), function(prot) sum(((exp((m11s$bR0b*2 + m11s$bR0j[,prot]*m11s$sigbR0j)*0.03)-1)*100) > 0) / nrow(m11s$bR0j))[
  order(sapply(1:ncol(m11s$bR0j), function(prot) mean((exp((m11s$bR0b*2 + m11s$bR0j[,prot]*m11s$sigbR0j)*0.03)-1)*100)), decreasing = T)]

highest_slope_order_means_in_rejectors <- (sapply(1:ncol(m11s$bR0j), function(prot) mean((exp((m11s$bR0b*2 + m11s$bR0j[,prot]*m11s$sigbR0j + 
                                                                                                 m11s$b_rej + m11s$bR0jr[,prot] * m11s$sigbR0jr)*0.03)-1)*100)))[
                        order(sapply(1:ncol(m11s$bR0j), function(prot) mean((exp((m11s$bR0b*2 + m11s$bR0j[,prot]*m11s$sigbR0j)*0.03)-1)*100)), decreasing = T)]
highest_slope_order_probsGreatZero_in_rejectors <- sapply(1:ncol(m11s$bR0j), function(prot) sum(((exp((m11s$bR0b*2 + m11s$bR0j[,prot]*m11s$sigbR0j + m11s$b_rej + m11s$bR0jr[,prot] * m11s$sigbR0jr)*0.03)-1)*100) > 0) / nrow(m11s$bR0j))[
  order(sapply(1:ncol(m11s$bR0j), function(prot) mean((exp((m11s$bR0b*2 + m11s$bR0j[,prot]*m11s$sigbR0j)*0.03)-1)*100)), decreasing = T)]


post_transplant_slope_per_day <- cbind(protein_name = as.character(highest_slope_order_names), 
                                       mean_spectral_count_proportional_percentage_change_per_day_in_nonrejectors = round(highest_slope_order_means_in_nonrejectors, 3),
                                       probability_change_per_day_is_positive_in_nonrejectors = round(highest_slope_order_probsGreatZero_in_nonrejectors, 4), 
                                       mean_spectral_count_proportional_percentage_change_per_day_in_rejectors = round(highest_slope_order_means_in_rejectors, 3),
                                       probability_change_per_day_is_positive_in_rejectors = round(highest_slope_order_probsGreatZero_in_rejectors, 4))
write.csv(post_transplant_slope_per_day, file = "post-transplant_pre-rejection_slope_values_m11_noP9.csv")


m11@formula
#relative change through time per unit day values, post-transplant / post-rejection
dens((exp((m11s$bR1b*2 + m11s$bR1j[,1]*m11s$sigbR1j)*0.03)-1)*100, add = F)
for(i in 2:ncol(m11s$A0j)){
  dens((exp((m11s$bR1b*2 + m11s$bR1j[,i]*m11s$sigbR1j)*0.03)-1)*100, add = T)
  
}
hist(sapply(1:ncol(m11s$bR1j), function(prot) mean((exp((m11s$bR1b*2 + m11s$bR1j[,prot]*m11s$sigbR1j)*0.03)-1)*100)))
highest_slope_order_names <- prot_key[order(sapply(1:ncol(m11s$bR1j), function(prot) mean((exp((m11s$bR1b*2 + m11s$bR1j[,prot]*m11s$sigbR1j)*0.03)-1)*100)), decreasing = T),1]

highest_slope_order_means <- sort(sapply(1:ncol(m11s$bR1j), function(prot) mean((exp((m11s$bR1b*2 + m11s$bR1j[,prot]*m11s$sigbR1j)*0.03)-1)*100)), decreasing = T)
highest_slope_order_probsGreatZero <- sapply(1:ncol(m11s$bR1j), function(prot) sum(((exp((m11s$bR1b*2 + m11s$bR1j[,prot]*m11s$sigbR1j)*0.03)-1)*100) > 0) / nrow(m11s$bR1j))[
  order(sapply(1:ncol(m11s$bR1j), function(prot) mean((exp((m11s$bR1b*2 + m11s$bR1j[,prot]*m11s$sigbR1j)*0.03)-1)*100)), decreasing = T)]


post_rejection_slope_per_day <- cbind(protein_name = as.character(highest_slope_order_names), 
                                      mean_spectral_count_proportional_percentage_change_per_day = round(highest_slope_order_means_in_nonrejectors, 3),
                                      probability_change_per_day_is_positive = round(highest_slope_order_probsGreatZero_in_nonrejectors, 4))
write.csv(post_transplant_slope_per_day, file = "post-transplant_post-rejection_slope_values_m11_noP9.csv")

m11@formula
#discrete drop at rejection event
dens(exp(m11s$A1b + m11s$A1j[,1]*m11s$sigA1j)*100, add = F)
for(i in 2:ncol(m11s$A1j)){
  dens(exp(m11s$A1b + m11s$A1j[,i]*m11s$sigA1j)*100, add = T, col = rgb(0, 0, 0, 0.15))
}
quantile(sapply(1:ncol(m11s$A1j), function(prot) mean(exp(m11s$A1b + m11s$A1j[,prot]*m11s$sigA1j))), probs = 0:100/100)
change_prot_inds_incr <- order(sapply(1:ncol(m11s$A1j), function(prot) mean(exp(m11s$A1b + m11s$A1j[,prot]*m11s$sigA1j))), decreasing = F)
highest_changes_order_names <- prot_key[change_prot_inds_incr,1]
highest_changes_order_means <- 100 * sort(sapply(1:ncol(m11s$A1j), function(prot) mean(exp(m11s$A1b + m11s$A1j[,prot]*m11s$sigA1j))), decreasing = F)
prob_counts_went_down <- (sapply(1:ncol(m11s$A1j), function(prot) sum(exp(m11s$A1b + m11s$A1j[,prot]*m11s$sigA1j) < 1) / nrow(m11s$A1j)))[change_prot_inds_incr]

at_rejection_change <- cbind(protein_name = as.character(highest_changes_order_names), 
                             mean_percent_spectral_count_change_at_rejection_event = round(highest_changes_order_means, 3),
                             probability_spectral_count_went_down = round(prob_counts_went_down, 4))
write.csv(at_rejection_change, file = "discrete_rejection_event_percent_change_m11_noP9.csv")

###############################
### interpreting new models ###
###############################

setwd("/Volumes/1TB/heart_transplant/")
library(rethinking)

#WAICs, lower is better
load("mf1_np_g_l") #WAIC 53253.89, PSIS 53387.26, wins unambiguously, WAIC_weight = 1
load("mf1_np_g_nl") #WAIC 53096.38, PSIS 53283.19
load("mf1_np_ng_nl") #WAIC 69634.83, PSIS 69929.92
load("mf1_np_ng_l") #WAIC 69851.41, PSIS 70072.67

WAIC(mf1_g_l)
PSIS(mf1_g_l)
WAIC(mf1_g_nl)
PSIS(mf1_g_nl)
WAIC(mf1_ng_nl)
PSIS(mf1_ng_nl)
WAIC(mf1_ng_l)
PSIS(mf1_ng_l)
compare(mf1_g_l, mf1_g_nl)

load("mf2_op_g_l_noP9") #WAIC 41759.4, WAIC weight = 1, PSIS weight = 1
load("mf2_op_g_nl_noP9") #WAIC 41853.1, 
compare(mf2_g_l, mf2_g_nl,func= PSIS)
precis(mf2_g_nl)

#interpreting the stairstep models
setwd("/Volumes/1TB/heart_transplant/")
library(rethinking)
load("m12_noncentered")
load("m12ng_noncentered")
load("m12_noncentered_noP9")

####################################################################
### interpreting simple non-gamma model of just the new patients ###
####################################################################
load("mf1_np_ng_nl") #WAIC 69634.83, PSIS 69929.92
precis_mf1_ng_nl <- precis(mf1_ng_nl)
lppd_mf1_ng_nl <- lppd(mf1_ng_nl)


mf1_ng_nls <- extract.samples(mf1_ng_nl, n = 10000)
precis(mf1_ng_nl)

load(file = "m8_pat_key_newpatients")
pat_key <- pat_key[order(pat_key[,2]),]
load(file = "m8_prot_key_newpatients")
prot_key <- prot_key[order(prot_key[,2]),]

mf1_ng_nl@formula
#starting intercept values
dens(exp(mf1_ng_nls$A0b + mf1_ng_nls$A0j[,1]), add = F)
for(i in 2:ncol(mf1_ng_nls$A0j)){
  dens(exp(mf1_ng_nls$A0b + mf1_ng_nls$A0j[,i]*mf1_ng_nls$sigA0j), add = T)
}
hist(sapply(1:ncol(mf1_ng_nls$A0j), function(prot) mean(exp(mf1_ng_nls$A0b + mf1_ng_nls$A0j[,prot]*mf1_ng_nls$sigA0j))))
highest_intercepts_order_names <- prot_key[order(sapply(1:ncol(mf1_ng_nls$A0j), function(prot) mean(exp(mf1_ng_nls$A0b + mf1_ng_nls$A0j[,prot]*mf1_ng_nls$sigA0j))), decreasing = T),1]
highest_intercepts_order_means <- sort(sapply(1:ncol(mf1_ng_nls$A0j), function(prot) mean(exp(mf1_ng_nls$A0b + mf1_ng_nls$A0j[,prot]*mf1_ng_nls$sigA0j))), decreasing = T)
pre_transplant_intercept <- cbind(protein_name = as.character(highest_intercepts_order_names), mean_spectral_count = round(highest_intercepts_order_means, 3))
write.csv(pre_transplant_intercept, file = "pre-transplant_intercept_values_mf1_NewPatients_noGamma_noLog.csv")

mf1_ng_nl@formula
#change through time per unit day values
dens((exp((mf1_ng_nls$bR0b*0.01 + mf1_ng_nls$bR0j[,1]*mf1_ng_nls$sigbR0j))-1)*100, add = F)
for(i in 2:ncol(mf1_ng_nls$A0j)){
  dens((exp((mf1_ng_nls$bR0b*0.01 + mf1_ng_nls$bR0j[,i]*mf1_ng_nls$sigbR0j))-1)*100, add = T)
  
}
hist(sapply(1:ncol(mf1_ng_nls$bR0j), function(prot) mean((exp((mf1_ng_nls$bR0b*2 + mf1_ng_nls$bR0j[,prot]*mf1_ng_nls$sigbR0j)*0.03)-1)*100)))
highest_slope_order_names <- prot_key[order(sapply(1:ncol(mf1_ng_nls$bR0j), function(prot) mean((exp((mf1_ng_nls$bR0b*0.01 + mf1_ng_nls$bR0j[,prot]*mf1_ng_nls$sigbR0j))-1)*100)), decreasing = T),1]
highest_slope_order_means <- sort(sapply(1:ncol(mf1_ng_nls$bR0j), function(prot) mean((exp((mf1_ng_nls$bR0b*0.01 + mf1_ng_nls$bR0j[,prot]*mf1_ng_nls$sigbR0j))-1)*100)), decreasing = T)
highest_slope_order_propsGreatZero <- sapply(1:ncol(mf1_ng_nls$bR0j), function(prot) sum(((exp((mf1_ng_nls$bR0b*0.01 + mf1_ng_nls$bR0j[,prot]*mf1_ng_nls$sigbR0j))-1)*100) > 0) / nrow(mf1_ng_nls$bR0j))[
  order(sapply(1:ncol(mf1_ng_nls$bR0j), function(prot) mean((exp((mf1_ng_nls$bR0b*0.01 + mf1_ng_nls$bR0j[,prot]*mf1_ng_nls$sigbR0j))-1)*100)), decreasing = T)]
post_transplant_slope_per_day <- cbind(protein_name = as.character(highest_slope_order_names), 
                                       mean_spectral_count_proportional_percentage_change_per_day = round(highest_slope_order_means, 3),
                                       probability_change_per_day_is_positive = round(highest_slope_order_propsGreatZero, 3))
write.csv(post_transplant_slope_per_day, file = "post-transplant_slope_values_mf1_NewPatients_noGamma_noLog.csv")


order(sapply(1:length(mf1_ng_nls$bR0j[1,]), function(prot) mean(mf1_ng_nls$bR0j[,prot])), decreasing = T)

#############################################
##same analysis, just using the old priors###
#############################################

load("mf1_np_ng_nl_oldPriors")
load("m8_data_newpatients")
summary(mf1_ng_nl)
precis_mf1_ng_nl <- precis(mf1_ng_nl, depth = 2)
lppd_mf1_ng_nl <- lppd(mf1_ng_nl)
mf1_ng_nls <- extract.samples(mf1_ng_nl, n = 10000)


load("mf1_ng_nls_newpats_oldpriors_samples")
load("mf1_ng_nls_newpats_oldpriors_precis")


load(file = "m8_pat_key_newpatients")
pat_key <- pat_key[order(pat_key[,2]),]
load(file = "m8_prot_key_newpatients")
prot_key <- prot_key[order(prot_key[,2]),]

mf1_ng_nl@formula
#starting intercept values
dens(exp(mf1_ng_nls$A0b + mf1_ng_nls$A0j[,1]), add = F)
for(i in 2:ncol(mf1_ng_nls$A0j)){
  dens(exp(mf1_ng_nls$A0b + mf1_ng_nls$A0j[,i]*mf1_ng_nls$sigA0j), add = T)
}
hist(sapply(1:ncol(mf1_ng_nls$A0j), function(prot) mean(exp(mf1_ng_nls$A0b + mf1_ng_nls$A0j[,prot]*mf1_ng_nls$sigA0j))))
highest_intercepts_order_names <- prot_key[order(sapply(1:ncol(mf1_ng_nls$A0j), function(prot) mean(exp(mf1_ng_nls$A0b + mf1_ng_nls$A0j[,prot]*mf1_ng_nls$sigA0j))), decreasing = T),1]
highest_intercepts_order_means <- sort(sapply(1:ncol(mf1_ng_nls$A0j), function(prot) mean(exp(mf1_ng_nls$A0b + mf1_ng_nls$A0j[,prot]*mf1_ng_nls$sigA0j))), decreasing = T)
pre_transplant_intercept <- cbind(protein_name = as.character(highest_intercepts_order_names), mean_spectral_count = round(highest_intercepts_order_means, 3))
write.csv(pre_transplant_intercept, file = "pre-transplant_intercept_values_mf1_NewPatients_noGamma_noLog_oldPriors.csv")

mf1_ng_nl@formula
#change through time per unit day values
dens((exp((mf1_ng_nls$bR0b*0.01 + mf1_ng_nls$bR0j[,1]*mf1_ng_nls$sigbR0j))-1)*100, add = F, xlim = c(-0.5, 0.75), ylim = c(0,10))
for(i in 2:ncol(mf1_ng_nls$A0j)){
  dens((exp((mf1_ng_nls$bR0b*0.01 + mf1_ng_nls$bR0j[,i]*mf1_ng_nls$sigbR0j))-1)*100, add = T)
}
hist(sapply(1:ncol(mf1_ng_nls$bR0j), function(prot) mean((exp((mf1_ng_nls$bR0b*0.01 + mf1_ng_nls$bR0j[,prot]*mf1_ng_nls$sigbR0j))-1)*100)))
highest_slope_order_names <- prot_key[order(sapply(1:ncol(mf1_ng_nls$bR0j), function(prot) mean((exp((mf1_ng_nls$bR0b*0.01 + mf1_ng_nls$bR0j[,prot]*mf1_ng_nls$sigbR0j))-1)*100)), decreasing = T),1]
highest_slope_order_means <- sort(sapply(1:ncol(mf1_ng_nls$bR0j), function(prot) mean((exp((mf1_ng_nls$bR0b*0.01 + mf1_ng_nls$bR0j[,prot]*mf1_ng_nls$sigbR0j))-1)*100)), decreasing = T)
highest_slope_order_propsGreatZero <- sapply(1:ncol(mf1_ng_nls$bR0j), function(prot) sum(((exp((mf1_ng_nls$bR0b*0.01 + mf1_ng_nls$bR0j[,prot]*mf1_ng_nls$sigbR0j))-1)*100) > 0) / nrow(mf1_ng_nls$bR0j))[
  order(sapply(1:ncol(mf1_ng_nls$bR0j), function(prot) mean((exp((mf1_ng_nls$bR0b*0.01 + mf1_ng_nls$bR0j[,prot]*mf1_ng_nls$sigbR0j))-1)*100)), decreasing = T)]
post_transplant_slope_per_day <- cbind(protein_name = as.character(highest_slope_order_names), 
                                       mean_spectral_count_proportional_percentage_change_per_day = round(highest_slope_order_means, 3),
                                       probability_change_per_day_is_positive = round(highest_slope_order_propsGreatZero, 3))
write.csv(post_transplant_slope_per_day, file = "post-transplant_slope_values_mf1_NewPatients_noGamma_noLog_oldPriors.csv")


order(sapply(1:length(mf1_ng_nls$bR0j[1,]), function(prot) mean(mf1_ng_nls$bR0j[,prot])), decreasing = T)


#interpreting the stairstep models
setwd("/Volumes/1TB/heart_transplant/")
library(rethinking)
load("m16_noncentered_noP9")
load("m12_noncentered_noP9")
load("m13_noncentered_noP9")
WAIC(m16) #98030.76, weight = 1
WAIC(m12) #103238.7
WAIC(m13) #103240.9
compare(m13, m16)
lppd_m16 <- lppd(m16)
precis_m16 <- precis(m16, depth = 2)
precis_m16_simple <- precis(m16, depth = 1)
precis(m13)


## generating excel files for model 16
m16s <- extract.samples(m16, n = 10000)
load(file = "m8_pat_key_noP9")
pat_key <- pat_key[order(pat_key[,2]),]
load(file = "m8_prot_key_noP9")
prot_key <- prot_key[order(prot_key[,2]),]


m16@formula
#starting intercept values
dens(exp(m16s$A0b + m16s$A0j[,1]*m16s$sigA0j), add = F, ylim = c(0,10), xlim = c(0,5), col = rgb(0,0,0,0.5))
for(i in 2:ncol(m16s$A0j)){
  dens(exp(m16s$A0b + m16s$A0j[,i]*m16s$sigA0j), add = T, col = rgb(0,0,0,0.05))
}
hist(sapply(1:ncol(m16s$A0j), function(prot) mean(exp(m16s$A0b + m16s$A0j[,prot]*m16s$sigA0j)))) #new patients
hist(sapply(1:ncol(m16s$A0j), function(prot) mean(exp(m16s$A0b + m16s$A0j[,prot]*m16s$sigA0j + m16s$a0_rej + m16s$A0jr[,prot]*m16s$sigA0jr))), add = T) 

order_names <- order(sapply(1:ncol(m16s$A0j), function(prot) mean(exp(m16s$A0b + m16s$A0j[,prot]*m16s$sigA0j))), decreasing = T)
highest_intercepts_order_names <- prot_key[order_names,1]
highest_intercepts_order_means <- sort(sapply(1:ncol(m16s$A0j), function(prot) mean(exp(m16s$A0b + m16s$A0j[,prot]*m16s$sigA0j))), decreasing = T)
pre_transplant_intercept <- cbind(protein_name = as.character(highest_intercepts_order_names), mean_spectral_count = round(highest_intercepts_order_means, 3))
write.csv(pre_transplant_intercept, file = "pre-transplant_intercept_values_stairstepM16_non-rejecting.csv")

order_names <- order(sapply(1:ncol(m16s$A0j), function(prot) mean(exp(m16s$A0b + m16s$A0j[,prot]*m16s$sigA0j 
                                                                      + m16s$a0_rej + m16s$A0jr[,prot]*m16s$sigA0jr))), decreasing = T)
highest_intercepts_order_names <- prot_key[order_names,1]
highest_intercepts_order_means <- sort(sapply(1:ncol(m16s$A0j), function(prot) mean(exp(m16s$A0b + m16s$A0j[,prot]*m16s$sigA0j
                                                                                                  + m16s$a0_rej + m16s$A0jr[,prot]*m16s$sigA0jr))), decreasing = T)
pre_transplant_intercept <- cbind(protein_name = as.character(highest_intercepts_order_names), mean_spectral_count = round(highest_intercepts_order_means, 3))
write.csv(pre_transplant_intercept, file = "pre-transplant_intercept_values_stairstepM16_rejecting.csv")

#change when entering phase 2, post-transplant / pre-rejection
hist(sapply(1:ncol(m16s$A1j), function(prot) mean(exp(m16s$A1b + m16s$A1j[,prot]*m16s$sigA1j)))) 
hist(sapply(1:ncol(m16s$A1j), function(prot) mean(exp(m16s$A1b + m16s$A1j[,prot]*m16s$sigA1j + m16s$a1_rej + m16s$A1jr[,prot]*m16s$sigA1jr))), add = F) 

order_names <- order(sapply(1:ncol(m16s$A1j), function(prot) mean(exp(m16s$A1b + m16s$A1j[,prot]*m16s$sigA1j))), decreasing = T)
highest_intercepts_order_names <- prot_key[order_names,1]
highest_intercepts_order_means <- sort(sapply(1:ncol(m16s$A1j), function(prot) mean(exp(m16s$A1b + m16s$A1j[,prot]*m16s$sigA1j)))*100, decreasing = T)
prob_counts_went_down <- (sapply(1:ncol(m16s$A1j), function(prot) sum(exp(m16s$A1b + m16s$A1j[,prot]*m16s$sigA1j) < 1) / nrow(m16s$A1j)))[order_names]

post_transplant_change <- cbind(protein_name = as.character(highest_intercepts_order_names), 
                                mean_percent_spectral_count_change_at_transplant_event = round(highest_intercepts_order_means, 3),
                                probability_spectral_count_went_down = prob_counts_went_down)
write.csv(post_transplant_change, file = "postTransplant_preRejection_Change_stairstepM16_non-rejecting.csv")


order_names <- order(sapply(1:ncol(m16s$A1j), function(prot) mean(exp(m16s$A1b + m16s$A1j[,prot]*m16s$sigA1j
                                                                      + m16s$a1_rej + m16s$A1jr[,prot]*m16s$sigA1jr))), decreasing = T)
highest_intercepts_order_names <- prot_key[order_names,1]
highest_intercepts_order_means <- sort(sapply(1:ncol(m16s$A1j), function(prot) mean(exp(m16s$A1b + m16s$A1j[,prot]*m16s$sigA1j
                                                                                        + m16s$a1_rej + m16s$A1jr[,prot]*m16s$sigA1jr)))*100, decreasing = T)
prob_counts_went_down <- (sapply(1:ncol(m16s$A1j), function(prot) sum(exp(m16s$A1b + m16s$A1j[,prot]*m16s$sigA1j
                                                                          + m16s$a1_rej + m16s$A1jr[,prot]*m16s$sigA1jr) < 1) / nrow(m16s$A1j)))[order_names]

post_transplant_change <- cbind(protein_name = as.character(highest_intercepts_order_names), 
                                mean_percent_spectral_count_change_at_transplant_event = round(highest_intercepts_order_means, 3),
                                probability_spectral_count_went_down = prob_counts_went_down)
write.csv(post_transplant_change, file = "postTransplant_preRejection_Change_stairstepM16_rejecting.csv")

#discrete drop at rejection event
change_prot_inds_incr <- order(sapply(1:ncol(m16s$A2j), function(prot) mean(exp(m16s$A2b + m16s$A2j[,prot]*m16s$sigA2j))), decreasing = F)
highest_changes_order_names <- prot_key[change_prot_inds_incr,1]
highest_changes_order_means <- 100 * sort(sapply(1:ncol(m16s$A2j), function(prot) mean(exp(m16s$A2b + m16s$A2j[,prot]*m16s$sigA2j))), decreasing = F)
prob_counts_went_down <- (sapply(1:ncol(m16s$A2j), function(prot) sum(exp(m16s$A2b + m16s$A2j[,prot]*m16s$sigA2j) < 1) / nrow(m16s$A2j)))[change_prot_inds_incr]
at_rejection_change <- cbind(protein_name = as.character(highest_changes_order_names), 
                             mean_percent_spectral_count_change_at_rejection_event = round(highest_changes_order_means, 3),
                             probability_spectral_count_went_down = round(prob_counts_went_down, 4))
write.csv(at_rejection_change, file = "discrete_rejection_event_percent_change__stairstepM16.csv")

load("m8_data_noP9")
sapply(1:length(unique(d_sub$patient_id)), function(pat) 
  d_sub$spectral_count[d_sub$patient_id == pat & d_sub$afterTransplant == 0] / d_sub$spectral_count[d_sub$patient_id == pat & d_sub$afterTransplant == 1 & d_sub$afterRejection == 0]
)


#look at pearson residuals for new patient model with old priors
setwd("/Volumes/1TB/heart_transplant/")
library(rethinking)
load("mf1_np_ng_nl_oldPriors")
load("m8_data_newpatients")
summary(mf1_ng_nl)
lppd_mf1_ng_nl <- lppd(mf1_ng_nl)


load(file = "m8_pat_key_newpatients")
pat_key <- pat_key[order(pat_key[,2]),]
load(file = "m8_prot_key_newpatients")
prot_key <- prot_key[order(prot_key[,2]),]
sort(lppd_mf1_ng_nl)[1:10]
d_sub[order(lppd_mf1_ng_nl)[1:10],]
pat_key
prot_key[660,] #model does not like the super high spectral count

mf1_ng_nlp <- link(mf1_ng_nl, n = 2000)
hist(d_sub$spectral_count, xlim = c(0,50), breaks = 200, col = rgb(0,0,1,alpha = 0.5), ylim = c(0,5000))
hist(rpois(n = length(mf1_ng_nlp$lambda[1,]), lambda = mf1_ng_nlp$lambda[1,]), add = T, col = rgb(1,0,0,alpha = 0.5), breaks = 200)
hist(sapply(1:length(d_sub$spectral_count), function(x) ppois(d_sub$spectral_count[x], mf1_ng_nlp$lambda[,x])))
# hist(sapply(1:length(d_sub$spectral_count), function(x) sum(rpois(2000,mf1_ng_nlp$lambda[,x]) >= d_sub$spectral_count[x])/2000))

#plot raw residuals from mean expectation
plot(sapply(1:length(mf1_ng_nlp$lambda[1,]), function(x) mean(mf1_ng_nlp$lambda[,x])), 
     sapply(1:length(mf1_ng_nlp$lambda[1,]), function(x) mean(mf1_ng_nlp$lambda[,x])) - d_sub$spectral_count,
     col = rgb(0,0,0,0.3)) 

#plot pearson residuals -- mean lambda
plot(sapply(1:length(mf1_ng_nlp$lambda[1,]), function(x) mean(mf1_ng_nlp$lambda[,x])), 
     (sapply(1:length(mf1_ng_nlp$lambda[1,]), function(x) mean(mf1_ng_nlp$lambda[,x])) - d_sub$spectral_count) / sqrt(sapply(1:length(mf1_ng_nlp$lambda[1,]), function(x) mean(mf1_ng_nlp$lambda[,x]))),
     xlab = "Predicted Rate", ylab = "Standardized Residual from Predicted Rate")
hist((sapply(1:length(mf1_ng_nlp$lambda[1,]), function(x) mean(mf1_ng_nlp$lambda[,x])) - d_sub$spectral_count) / sqrt(sapply(1:length(mf1_ng_nlp$lambda[1,]), function(x) mean(mf1_ng_nlp$lambda[,x]))),
     breaks = 500, xlim = c(-5, 5), freq = F)
lines(seq(-5, 5, length.out = 200), dnorm(x = seq(-5, 5, length.out = 200), mean = 0, sd = 1))

#plot pearson residuals -- sampled lambda
numSamps <- 400
indices <- floor(seq(from = 1, to = length(mf1_ng_nlp$lambda[,1]), length.out = numSamps))
resids <- rep(0, length(d_sub$spectral_count)*length(indices))
for(x in 1:length(d_sub$spectral_count)){
  if (x %% 100 == 0) {cat(paste(x, " "))}
  resids[((x-1)*numSamps+1):(x*numSamps)] <- (mf1_ng_nlp$lambda[,x][indices] - d_sub$spectral_count[x]) / sqrt(mf1_ng_nlp$lambda[,x][indices])
}
hist(resids, breaks = 500, xlim = c(-5, 5), freq = F)
lines(seq(-5, 5, length.out = 200), dnorm(x = seq(-5, 5, length.out = 200), mean = 0, sd = 1))
