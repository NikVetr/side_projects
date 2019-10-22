

#NOTE: This assumes that consumption in different food categories is independent, which is definitely not the case
iter <- 1e4
percDiff <- vector(length = iter)
for(i in 1:iter){
  if(i %% 1000 == 0) {print(i)}
  
  numIndiv_control <- 864
  prob_veg_control <- .074 #probability with which a given person is veg*ns
  abstain_veg_control <- rbinom( numIndiv_control , 1 , prob_veg_control ) #simulating veg*ns
  
  #pork
  prob_abstain_completely_control1 <- .20
  #these are for non-veg*ns who don't eat food from this category for various reasons (e.g. religious, cultural, taste, etc.)
  #For simplicity, I'm assuming these are independent from individuals who don't eat this food due to having a veg*n diet    
  rate_consume_control1 <- 2 #the rate at which non-abstainers eat from food category 1 -- i.e. pork
  abstain_control1 <- rbinom( numIndiv_control , 1 , prob_abstain_completely_control1 ) #simulating non-veg*n abstainers
  abstain_control1[abstain_veg_control == 1] <- 1 #ensuring veg*n abstainers are represented
  cons_control1 <- (1-abstain_control1)*rpois( numIndiv_control , rate_consume_control1 ) #simulating servings of pork consumed
  
  #beef
  prob_abstain_completely_control2 <- .15 
  rate_consume_control2 <- 2
  abstain_control2 <- rbinom( numIndiv_control , 1 , prob_abstain_completely_control2 )
  abstain_control2[abstain_veg_control == 1] <- 1
  cons_control2 <- (1-abstain_control2)*rpois( numIndiv_control , rate_consume_control2 )
  
  #bird
  prob_abstain_completely_control3 <- .05 
  rate_consume_control3 <- 3
  abstain_control3 <- rbinom( numIndiv_control , 1 , prob_abstain_completely_control3 )
  abstain_control3[abstain_veg_control == 1] <- 1
  cons_control3 <- (1-abstain_control3)*rpois( numIndiv_control , rate_consume_control3 )
  
  #sea
  prob_abstain_completely_control4 <- .1 
  rate_consume_control4 <- 1
  abstain_control4 <- rbinom( numIndiv_control , 1 , prob_abstain_completely_control4 )
  abstain_control4[abstain_veg_control == 1] <- 1
  cons_control4 <- (1-abstain_control4)*rpois( numIndiv_control , rate_consume_control4 )
  
  cons_control <- cons_control1 + cons_control2 + cons_control3 + cons_control4
  
  perc_zeros_control <- sum(cons_control == 0)/length(cons_control) * 100
  
  ##### Treatment Group #####
  
  numIndiv_treat <- 934
  prob_veg_treat <- .074
  abstain_veg_treat <- rbinom( numIndiv_treat , 1 , prob_veg_treat )
  
  #pork
  prob_abstain_completely_treat1 <- .20 
  rate_consume_treat1 <- 2
  abstain_treat1 <- rbinom( numIndiv_treat , 1 , prob_abstain_completely_treat1 )
  abstain_treat1[abstain_veg_treat == 1] <- 1
  cons_treat1 <- (1-abstain_treat1)*rpois( numIndiv_treat , rate_consume_treat1 )

  #beef
  prob_abstain_completely_treat2 <- .15
  rate_consume_treat2 <- 2
  abstain_treat2 <- rbinom( numIndiv_treat , 1 , prob_abstain_completely_treat2 )
  abstain_treat2[abstain_veg_treat == 1] <- 1
  cons_treat2 <- (1-abstain_treat2)*rpois( numIndiv_treat , rate_consume_treat2 )

  #bird
  prob_abstain_completely_treat3 <- .05
  rate_consume_treat3 <- 3
  abstain_treat3 <- rbinom( numIndiv_treat , 1 , prob_abstain_completely_treat3 )
  abstain_treat3[abstain_veg_treat == 1] <- 1
  cons_treat3 <- (1-abstain_treat3)*rpois( numIndiv_treat , rate_consume_treat3 )

  #sea
  prob_abstain_completely_treat4 <- .1
  rate_consume_treat4 <- 1
  abstain_treat4 <- rbinom( numIndiv_treat , 1 , prob_abstain_completely_treat4 )
  abstain_treat4[abstain_veg_treat == 1] <- 1
  cons_treat4 <- (1-abstain_treat4)*rpois( numIndiv_treat , rate_consume_treat4 )
  
  cons_treat <- cons_treat1 + cons_treat2 + cons_treat3 + cons_treat4
  
  perc_zeros_treat <- sum(cons_treat == 0)/length(cons_treat) * 100
  
  percDiff[i] <- perc_zeros_treat - perc_zeros_control
}
dens(percDiff)
sum(percDiff > 2)/iter #"p-value" for a "one-tailed" test

#simple test of perc diff as stat
probVeg <- .074
diffs <- sapply(1:1e5, function(x) sum(rbinom(934, 1, prob = probVeg))/934 * 100 - sum(rbinom(864, 1, prob = probVeg))/864 * 100)
sum(diffs > 2) / length(diffs)

#ratio
probVeg <- .074
rats <- sapply(1:1e5, function(x) (sum(rbinom(934, 1, prob = probVeg))/934 * 100) / (sum(rbinom(864, 1, prob = probVeg))/864 * 100))
sum(rats > 1.3125) / length(rats)
