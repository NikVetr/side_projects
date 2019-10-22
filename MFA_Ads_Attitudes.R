#cite rethinking package and richard's book
#load packages
library(rethinking)
library(mcmcplots)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)

#who uses .xls any more lol? Let's save them as CSVs in excel and load them into R
# control <- read.csv(choose.files())
# experimental <- read.csv(choose.files()) #weird errors, hmmm...
# str(control)
# str(experimental)
# control$treat <- 0 #create a treatment variable
# experimental$treat <- 1
# d <- rbind(control, experimental) #combine the two data frames
# head(d) #check to make sure everything looks ok...
# str(d) 
#...actually, hot damn, let's not do it this way. Let's combine the data in excel and clean it up a bit there
#It's not as sexy or officially replicable as doing everything in R, but you can access my combined 
#and cleaned CSV file here:


# d <- read.csv(choose.files()) #read the CSV with the combined data in
d <- read.csv("C:\\Users\\Nikolai\\Dropbox\\MFA_Survey\\CombinedData.csv")
str(d) #much better!!
#let's reorder some levels
d$NowVs4MonPast <- factor(d$NowVs4MonPast, levels = c("less", "same", "more"))
d$NowVs4MonFut <- factor(d$NowVs4MonFut, levels = c("none", "less", "same", "more")) #why have a none in there?? none distinguishes an absolute amount,
#all the other ones are relative amounts! Someone eating none can be eating less or same! And the 4MonPast doesn't have a none option! Argh. Whatever.
#let's say "none" means "super less" since we get rid of the people who claimed to be vegetarians when they took the survey anyway
d$MeatUniqueSmart <- factor(d$MeatUniqueSmart, levels = c("Strongly disagree", "Somewhat disagree", "Neither agree nor disagree", "Somewhat agree",
                                                          "Strongly agree"))
d$LessMeatGood <- factor(d$LessMeatGood, levels = c("Strongly disagree", "Somewhat disagree", "Neither agree nor disagree", "Somewhat agree",
                                                    "Strongly agree"))
d$MeatReplace <- factor(d$MeatReplace, levels = c("Strongly disagree", "Somewhat disagree", "Neither agree nor disagree", "Somewhat agree",
                                                  "Strongly agree"))

#now let's look at and get a sense of the data

d$anyChangeP4M <- as.factor(paste0(d$P4MonIncr, d$P4MonNC, d$P4MonDecr))
summary(d$anyChangeP4M)
#since this asks about any change in the past 4 months, we get some sorta flip-floppy answers, where some people increased, decreased, 
#and did not change their behaviors, all at the same time. I dunno that I really want to look at this variable, especially since something 
#like "NowVs4MonPast" has a much more clear-cut interpretation. Here, we don't know if it's a lasting or super ephemeral change, if they 
#did it purely by chance or if they did it for unrelated reasons (like going backpacking or something)

#Let's look at the other data and clean it up a bit

summary(d$Gender) #there's an "Other_NS" option; I dunno what that means, but I'll just recode it as "Other"
d$Gender[which(d$Gender == "Other_NS")] <- "Other"
d$Gender <- factor(d$Gender)

summary(d$Age) #Why is there a "12 and under" category? Age is a perfectly good numeric variable!
d$c.Age <- as.character(d$Age)
d$c.Age[which(d$c.Age == "12 or younger" | d$c.Age == "12")] <- "12"
d$c.Age <- as.numeric(d$c.Age)
summary(d$c.Age)
#can you not ask people younger than 12 how old they are? Don't you need to be at least 13 to be on facebook anyway? 
#Whatever, there are only 85 pf them, for now I'll recode the 12 and unders as 12, which seems like a reasonable assumption 
#(there aren't going to be 5 year olds
#taking surveys on the interner lol) and I want age to be numeric and not categorical. I can later try coding it as missing or accomodating fuzziness
#in age assignment statistically (by modeling) or dumping all the data from 12 and unders, since it might have resulted from a different person taking 
#the survey than the one who watched the video from facebook. If it doesn't take too long I can rerun the analysis under each approach
#to assess sensitivity to assumptions
#let's also standardize this variable to ease later interpretation and improve MCMC performance
d$Age.s <- (d$c.Age - mean(d$c.Age)) / sd(d$c.Age) #also called "z-scores"

#not looking at NAs for now
consumption <- d[,c("Treat", "X2DaysPork", "X2DaysBeef", "X2DaysDairy", "X2DaysEggs", "X2DaysBird", "X2DaysSea")] %>% 
  gather(Type, Amount, -Treat)
ggplot(data = consumption, aes(x=Type, y=Amount, color=Treat)) +
  geom_jitter() +
  theme_bw() #do note -- this removes over half the data, which we'll impute later
#hmm, interesting. I guess when the survey asked "In the past two days, how many servings have 
# you had of.Port, Beef, Dairy, Eggs, Chicken, Fish?" it didn't allow responses above 9? I don't
#think it coded them as NAs, since all the variables have 228 NAs This might be failing to capture 
#some of the effect (I've known a lot of people to eat several dozen servings of dairy, beef, eggs, chicken
#etc. over 2 days), though maybe young women on facebook are less likely to do so than 250 lb male powerlifters
#Also, eyeballing the amounts consumed in this quick n' dirty exploratory data analysis, I don't see any effect
#of treatment. It is clear that people eat a lot of servings of dairy and relatively few servings of fish, though

#To standardize or not to standardize? Since these are outcomes and I sorta care about the units and I want to impute all those NAs, I won't for now

consumption <- d[,c("Gender", "X2DaysPork", "X2DaysBeef", "X2DaysDairy", "X2DaysEggs", "X2DaysBird", "X2DaysSea")] %>% 
  gather(Type, Amount, -Gender)
consumption <- consumption[consumption[,1] != "",]
ggplot(data = consumption, aes(x=Type, y=Amount, color=Gender)) +
  geom_jitter() +
  theme_bw()


consumption <- d[,c("c.Age", "X2DaysPork")] %>% gather(Type, Amount, -c.Age)
porkAge <- ggplot(data = consumption, aes(x=c.Age, y=Amount)) + geom_jitter(alpha = .2) + theme_bw() + ggtitle("Pork")

consumption <- d[,c("c.Age", "X2DaysBeef")] %>% gather(Type, Amount, -c.Age)
beefAge <- ggplot(data = consumption, aes(x=c.Age, y=Amount)) + geom_jitter(alpha = .2) + theme_bw() + ggtitle("Beef")

consumption <- d[,c("c.Age", "X2DaysDairy")] %>% gather(Type, Amount, -c.Age)
dairyAge <- ggplot(data = consumption, aes(x=c.Age, y=Amount)) + geom_jitter(alpha = .2) + theme_bw() + ggtitle("Dairy")

consumption <- d[,c("c.Age", "X2DaysEggs")] %>% gather(Type, Amount, -c.Age)
eggsAge <- ggplot(data = consumption, aes(x=c.Age, y=Amount)) + geom_jitter(alpha = .2) + theme_bw() + ggtitle("Eggs")

consumption <- d[,c("c.Age", "X2DaysBird")] %>% gather(Type, Amount, -c.Age)
birdAge <- ggplot(data = consumption, aes(x=c.Age, y=Amount)) + geom_jitter(alpha = .2) + theme_bw() + ggtitle("Bird")

consumption <- d[,c("c.Age", "X2DaysSea")] %>% gather(Type, Amount, -c.Age)
seaAge <- ggplot(data = consumption, aes(x=c.Age, y=Amount)) + geom_jitter(alpha = .2) + theme_bw() + ggtitle("Sea")

grid.arrange(porkAge, beefAge, dairyAge, eggsAge, birdAge, seaAge, ncol = 3)

summary(d$MeatUniqueSmart) 
summary(d$LessMeatGood)
summary(d$MeatReplace) #lotsa missing values here, and since 368 are missing from each I guess these people just didn't get to that page in the
#survey? Can we assume these are MCAR?
attitude <- d[,c("Treat", "MeatUniqueSmart", "LessMeatGood", "MeatReplace")] %>% 
  gather(Question, LikertAnswer, -Treat)
attitude$LikertAnswer <- factor(attitude$LikertAnswer, levels = c("Strongly disagree", "Somewhat disagree", "Neither agree nor disagree", 
                                                                  "Somewhat agree", "Strongly agree"))
ggplot(data = na.omit(attitude), aes(x=Question, y=LikertAnswer, color=Treat)) +
  geom_jitter() +
  theme_bw() #again, we omit NAs
#again, no obvious relationship with treatment

summary(d$Country) #eh

#Do we care about people who already didn't eat meat and still don't eat meat? I guess we might in case watching the ad persuaded them to 
#consider trying it in the future or changed their opinions on the  moral value of animals. Let's look at these people briefly
veg <- which(d$NoMeatAlready == "DidNotEatAlready")
d$AlrVeg <- ifelse(d$NoMeatAlready == "DidNotEatAlready",1,0) #let's create another column of this too
d$AlrVeg[which(d$anyChangeP4M == "" & d$AlrVeg == 0)] <- NA #preserve missing data
#such a shame this couldn't have asked "Where you a vegetarian 4 months ago?" without asking if they don't eat meat now. Oh well! As such, we'll underestimate
#the amount of vegetarians 4 months ago, since some of them will have started eating meat again, especially if recidivism rates are high. There also might be confusion 
#over what "meat" entails, since some people don't think fish are made of meat, and others don't think chickens are.
attitudeVeg <- d[veg ,c("Treat", "MeatUniqueSmart", "LessMeatGood", "MeatReplace")] %>% 
  gather(Question, LikertAnswer, -Treat)
attitudeVeg$LikertAnswer <- factor(attitudeVeg$LikertAnswer, levels = c("Strongly disagree", "Somewhat disagree", "Neither agree nor disagree", 
                                                                        "Somewhat agree", "Strongly agree"))
ggplot(data = na.omit(attitudeVeg), aes(x=Question, y=LikertAnswer)) +
  geom_jitter() +
  theme_bw()
summary(d$NowVs4MonPast[veg]) #so despite not eating meat then and now, some people are eating more meat and many less meat. Weird!
summary(d$NowVs4MonFut[veg]) #it looks like some of our non-meat eaters will be eating no meat 4 months from now, but most will be eating none and some less. That's good
#at least!

#how much meat are the folks who don't "eat meat, chicken, and fish" actually eating in the last few days?
par(mfrow = c(2,3))
hist(d$X2DaysPork[veg]); hist(d$X2DaysBeef[veg]); hist(d$X2DaysDairy[veg]); hist(d$X2DaysEggs[veg]); hist(d$X2DaysBird[veg]); hist(d$X2DaysSea[veg])
consumptionVeg <- d[veg,c("Treat", "X2DaysPork", "X2DaysBeef", "X2DaysDairy", "X2DaysEggs", "X2DaysBird", "X2DaysSea")] %>% 
  gather(Type, Amount, -Treat)
ggplot(data = consumptionVeg, aes(x=Type, y=Amount)) +
  geom_jitter() +
  theme_bw()
#most but not all consume no servings, though many are consuming dairy or eggs.
#I'm not super sure how to deal with them, especially the ones who are actually eating meat
length(veg)
#there are only 151 of them -- enough for each to get a first gen, Kanto Pokemon
#the data is suspect, but let's try to still use it in our model
#and let's only include the variables of interest, too -- the ones we'll actually use in our subsequent analysis
d$Age_s <- d$Age.s

d2 <- d[, c("Treat", "X2DaysPork", "X2DaysBeef", "X2DaysDairy", "X2DaysEggs", "X2DaysBird", "X2DaysSea", "NowVs4MonPast",
            "NowVs4MonFut", "MeatUniqueSmart", "LessMeatGood", "MeatReplace", "Gender", "Country", "Age_s")]

d2$Country[d2$Country == ""] <- NA
d2$country_id <- coerce_index(d2$Country)
d2$Gender[d2$Gender == ""] <- NA
d2$gender_id <- coerce_index(d2$Gender)
d2$Treat <- as.integer(d2$Treat) #change from factor to integer
#d2$AlrVeg <- as.integer(d2$AlrVeg)
str(d2)

# There might be some concomitant variable bias or post-treatment bias if we include the belief questions
#as predictor variables, unfortunately. Since the survey was administered after the treatment, including belief questions
#might mask the effects of the video intervention if the video people to lean toward vegetarianism by changing 
#their beliefs (since multiple regression asks what the value of knowing some bit of information is after you know
#other information, like knowing treatment status after belief). It would have been nice to know people's opinions 
#before seeing the video (cos then we could
#tease out if maybe the video was only effective for those who already believed animals had moral value or something)
#, but that would have been tricky given the model design See Rethinking: Model comparison doesn't help

#don't just compare a bajillion models because you might stumble upon a peculiarly overfit one by chance that fools
#WAIC and similar metrics via the curse of tippecanoe because, by chance, it's really well fit to the data
#can be usefully labeled as "data dredging" but still, exercise judgment when constructing models

set.seed(95616) #UC-Davis zip code

# What'll our candidate models be? Outcome variables are the X2Days... outcomes, NowVs4MonPast, NowVs4MonFut, and the attitude outcomes
#Predictor variables are Treat, Age.s, Gender, and Country
# 
# #the simplest model of interest
# ~Treat
# 
# ##two predictiors
# ~Treat + Age.s
# ~Treat + Treat*Age.S + Age.s
# #Age seems important, because as people age I'd anticipate them eating more due to higher metabolic requirements, not 
# # being as susceptible to marketing, and having greater control of their diets.
# 
# ~Treat + Gender 
# ~Treat + Treat*Gender + Gender
# #boys (and other?) stereotypically eat more meat than girls (since eating meat is the enterprise of manly men), and since the video
# #seems to be about why "Ariana Grande leaves meat off her plate", there might be some gendered response to the video due to in-group bias
# #or similarity effects or whatever
# 
# ~Treat + Country
# ~Treat + Treat*Country + Country
# #I'll happily buy that different countries can vary in their attitudes and meat consumption, and furthermore that Country might affect
# #the effectiveness of the video (since maybe Australians don't know who Ariana Grande is or something??)
# 
# ##Three predictors, one interaction
# ~Treat + Age.s + Gender
# ~Treat + Treat*Age.s + Age.s + Gender
# ~Treat + Treat*Gender + Gender + Age.s
# 
# ~Treat + Age.s + Country
# ~Treat + Treat*Age.s + Age.s + Country
# ~Treat + Treat*Country + Country + Age.s
# 
# ~Treat + Gender + Country
# ~Treat + Treat*Gender + Gender + Country
# ~Treat + Treat*Country + Country + Gender
# 
# #Four Predictors
# ~Treat + Gender + Country + Age.s
# ~Treat + Treat*Gender + Gender + Country + Age.s
# ~Treat + Treat*Country + Gender + Country + Age.s
# ~Treat + Treat*Age.s + Gender + Country + Age.s

#obviously there are a ton more models we can try, but let's look at how these do first.


dcc <- d2[ complete.cases(d2) , ] #see below -- we also need our observations to be equal for model comparison purposes


#recode responses as numbers 
dcc$LMG_i <- as.numeric(dcc$LessMeatGood)

#doing some plotting
par(mfrow = c(1,3))
simplehist(dcc$LMG_i, xlim = c(1,5), xlab = "response")
title("Histogram of Responses")

# discrete proportion of each response value
pr_k <- table( dcc$LMG_i ) / nrow(dcc)
# cumsum converts to cumulative proportions
cum_pr_k <- cumsum( pr_k )
# plot
plot( 1:5 , cum_pr_k , type="b" , xlab="response" ,
      ylab="cumulative proportion" , ylim=c(0,1) )
title("Cumulative Proportion of Responses")

logit <- function(x) log(x/(1-x)) # convenience function 11.4
lco <- logit( cum_pr_k ) 

plot( 1:5 , lco , type="b" , xlab="response" ,
      ylab="log cumulative odds")
title("Log Cumulative Odds of Responses")

# Fitting Models #

#simplest model, just intercepts
m1 <- map2stan(
  alist(
    LMG_i ~ dordlogit( phi , cutpoints ),
    phi <- 0, #note: link function already embedded in dordlogit likelihood
    cutpoints ~ dnorm(0,5)
  ) ,
  data= list(LMG_i = dcc$LMG_i) ,
  iter = 25000, warmup = 5000,
  start=list(cutpoints = c(-1.88, -.76, .62, 1.75)) ) #start the chain off somewhere decent
precis(m1, depth = 2)

#effect of treatment on phi
m2 <- map2stan(
  alist(
    LMG_i ~ dordlogit( phi , cutpoints ),
    phi <- Treat * pT, 
    cutpoints ~ dnorm(0,5),
    pT ~ dnorm(0,5)
  ) ,
  data= list(LMG_i = dcc$LMG_i, Treat = dcc$Treat) ,
  iter = 25000, warmup = 5000,
  start=list(cutpoints = c(-1.88, -.76, .62, 1.75)) ) 
precis(m2, depth = 2)

#effect of treatment on phi
m3 <- map2stan(
  alist(
    LMG_i ~ dordlogit( phi , cutpoints ),
    phi <- Treat * pT + Age_s * pA, 
    cutpoints ~ dnorm(0,5),
    c(pT, pA) ~ dnorm(0,5)
  ) ,
  data= list(LMG_i = dcc$LMG_i, Treat = dcc$Treat, Age_s = dcc$Age_s) ,
  iter = 25000, warmup = 5000,
  start=list(cutpoints = c(-1.88, -.76, .62, 1.75)) ) 
precis(m3, depth = 2)

m4 <- map2stan(
  alist(
    LMG_i ~ dordlogit( phi , cutpoints ),
    phi <- Treat * pT + Age_s * pA + pG[gender_id], 
    cutpoints ~ dnorm(0,5),
    c(pT, pA) ~ dnorm(0,5),
    pG[gender_id] ~ dnorm(0,sigGp),
    sigGp ~ dcauchy(0,1)
  ) ,
  data= list(LMG_i = dcc$LMG_i, Treat = dcc$Treat, Age_s = dcc$Age_s, gender_id = dcc$gender_id) ,
  iter = 25000, warmup = 5000,
  start=list(cutpoints = c(-1.88, -.76, .62, 1.75)) ) 
precis(m4, depth = 2)

m5 <- map2stan(
  alist(
    LMG_i ~ dordlogit( phi , cutpoints ),
    phi <- Treat * pT + Age_s * pA + pG[gender_id] + pC[country_id], 
    cutpoints ~ dnorm(0,5),
    c(pT, pA) ~ dnorm(0,5),
    pG[gender_id] ~ dnorm(0,sigGp),
    sigGp ~ dcauchy(0,1),
    pC[country_id] ~ dnorm(0,sigCp),
    sigCp ~ dcauchy(0,1)
  ) ,
  data= list(LMG_i = dcc$LMG_i, Treat = dcc$Treat, Age_s = dcc$Age_s, gender_id = dcc$gender_id, country_id = dcc$country_id) ,
  iter = 25000, warmup = 5000,
  start=list(cutpoints = c(-1.88, -.76, .62, 1.75)) ) 
precis(m5, depth = 2)

m6 <- map2stan(
  alist(
    LMG_i ~ dordlogit( phi , cutpoints ),
    phi <- Treat * pT + pG[gender_id] + pC[country_id], 
    cutpoints ~ dnorm(0,5),
    pT ~ dnorm(0,5),
    pG[gender_id] ~ dnorm(0,sigGp),
    sigGp ~ dcauchy(0,1),
    pC[country_id] ~ dnorm(0,sigCp),
    sigCp ~ dcauchy(0,1)
  ) ,
  data= list(LMG_i = dcc$LMG_i, Treat = dcc$Treat, gender_id = dcc$gender_id, country_id = dcc$country_id) ,
  iter = 25000, warmup = 5000,
  start=list(cutpoints = c(-1.88, -.76, .62, 1.75)) ) 
precis(m6, depth = 2)

m7 <- map2stan(
  alist(
    LMG_i ~ dordlogit( phi , cutpoints ),
    phi <- Treat * pT + Treat * Age_s * pTA + Age_s * pA + pG[gender_id] + pC[country_id], 
    cutpoints ~ dnorm(0,5),
    c(pT, pA, pTA) ~ dnorm(0,5),
    pG[gender_id] ~ dnorm(0,sigGp),
    sigGp ~ dcauchy(0,1),
    pC[country_id] ~ dnorm(0,sigCp),
    sigCp ~ dcauchy(0,1)
  ) ,
  data= list(LMG_i = dcc$LMG_i, Treat = dcc$Treat, Age_s = dcc$Age_s, gender_id = dcc$gender_id, country_id = dcc$country_id) ,
  iter = 25000, warmup = 5000,
  start=list(cutpoints = c(-1.88, -.76, .62, 1.75)) ) 
precis(m7, depth = 2)

m8 <- map2stan(
  alist(
    LMG_i ~ dordlogit( phi , cutpoints ),
    phi <- Treat * pT + Treat * Age_s * pTA + Age_s * pA, 
    cutpoints ~ dnorm(0,5),
    c(pT, pA, pTA) ~ dnorm(0,5)
  ) ,
  data= list(LMG_i = dcc$LMG_i, Treat = dcc$Treat, Age_s = dcc$Age_s) ,
  iter = 25000, warmup = 5000,
  start=list(cutpoints = c(-1.88, -.76, .62, 1.75)) ) 
precis(m8, depth = 2)

## No effect of treatment ##
#effect of treatment on phi

m9n <- map2stan(
  alist(
    LMG_i ~ dordlogit( phi , cutpoints ),
    phi <- Age_s * pA, 
    cutpoints ~ dnorm(0,5),
    pA ~ dnorm(0,5)
  ) ,
  data= list(LMG_i = dcc$LMG_i, Age_s = dcc$Age_s) ,
  iter = 25000, warmup = 5000,
  start=list(cutpoints = c(-1.88, -.76, .62, 1.75)) ) 
precis(m9n, depth = 2)

m10n <- map2stan(
  alist(
    LMG_i ~ dordlogit( phi , cutpoints ),
    phi <- Age_s * pA + pG[gender_id], 
    cutpoints ~ dnorm(0,5),
    pA ~ dnorm(0,5),
    pG[gender_id] ~ dnorm(0,sigGp),
    sigGp ~ dcauchy(0,1)
  ) ,
  data= list(LMG_i = dcc$LMG_i, Age_s = dcc$Age_s, gender_id = dcc$gender_id) ,
  iter = 25000, warmup = 5000,
  start=list(cutpoints = c(-1.88, -.76, .62, 1.75)) ) 
precis(m10n, depth = 2)

m11n <- map2stan(
  alist(
    LMG_i ~ dordlogit( phi , cutpoints ),
    phi <- Age_s * pA + pG[gender_id] + pC[country_id], 
    cutpoints ~ dnorm(0,5),
    pA ~ dnorm(0,5),
    pG[gender_id] ~ dnorm(0,sigGp),
    sigGp ~ dcauchy(0,1),
    pC[country_id] ~ dnorm(0,sigCp),
    sigCp ~ dcauchy(0,1)
  ) ,
  data= list(LMG_i = dcc$LMG_i, Age_s = dcc$Age_s, gender_id = dcc$gender_id, country_id = dcc$country_id) ,
  iter = 25000, warmup = 5000,
  start=list(cutpoints = c(-1.88, -.76, .62, 1.75)) ) 
precis(m11n, depth = 2)

m12n <- map2stan(
  alist(
    LMG_i ~ dordlogit( phi , cutpoints ),
    phi <- pG[gender_id] + pC[country_id], 
    cutpoints ~ dnorm(0,5),
    pG[gender_id] ~ dnorm(0,sigGp),
    sigGp ~ dcauchy(0,1),
    pC[country_id] ~ dnorm(0,sigCp),
    sigCp ~ dcauchy(0,1)
  ) ,
  data= list(LMG_i = dcc$LMG_i, gender_id = dcc$gender_id, country_id = dcc$country_id) ,
  iter = 25000, warmup = 5000,
  start=list(cutpoints = c(-1.88, -.76, .62, 1.75)) ) 
precis(m12n, depth = 2)

m13 <- map2stan(
  alist(
    LMG_i ~ dordlogit( phi , cutpoints ),
    phi <- Treat * pT + Treat * Age_s * pTA + Age_s * pA + pG[gender_id] + pC[country_id] + Treat*pGT[gender_id], 
    cutpoints ~ dnorm(0,5),
    c(pT, pA, pTA) ~ dnorm(0,5),
    pG[gender_id] ~ dnorm(0,sigGp),
    sigGp ~ dcauchy(0,1),
    pGT[gender_id] ~ dnorm(0,sigGTp),
    sigGTp ~ dcauchy(0,1),
    pC[country_id] ~ dnorm(0,sigCp),
    sigCp ~ dcauchy(0,1)
  ) ,
  data= list(LMG_i = dcc$LMG_i, Treat = dcc$Treat, Age_s = dcc$Age_s, gender_id = dcc$gender_id, country_id = dcc$country_id) ,
  iter = 25000, warmup = 5000,
  start=list(cutpoints = c(-1.88, -.76, .62, 1.75)) ) 
precis(m13, depth = 2)

m14 <- map2stan(
  alist(
    LMG_i ~ dordlogit( phi , cutpoints ),
    phi <- Treat * pT + pG[gender_id] + Treat*pGT[gender_id], 
    cutpoints ~ dnorm(0,5),
    pT ~ dnorm(0,5),
    pG[gender_id] ~ dnorm(0,sigGp),
    sigGp ~ dcauchy(0,1),
    pGT[gender_id] ~ dnorm(0,sigGTp),
    sigGTp ~ dcauchy(0,1)
  ) ,
  data= list(LMG_i = dcc$LMG_i, Treat = dcc$Treat, Age_s = dcc$Age_s, gender_id = dcc$gender_id, country_id = dcc$country_id) ,
  iter = 25000, warmup = 5000,
  start=list(cutpoints = c(-1.88, -.76, .62, 1.75)) ) 
precis(m14, depth = 2)

compare(m1, m2, m3, m4, m5, m6, m7, m8, m9n, m10n, m11n, m12n, m13, m14, refresh = 0.1)

coeftab(m12n, m11n, m6, m5, m10n, m7, m13, m14, m4)

d.pred <- data.frame(
  Treat = 1, 
  gender_id = 1,
  country_id = 6,
  Age_s = 0
)

sims <- ensemble(m12n, m11n, m6, m5, m10n, m7, m13, m14, m4, data = d.pred, n = 1e3)
a <- sim(m6, data = d.pred)



#posterior predictions from model 6

compare(m1, m2, m3, m4, m5, m6, m7, m8, m13, m14)

#extract samples from the posterior

nsamp <- 200 #number of samples to extract from the posterior and subsequently plot
post <- extract.samples( m6 , n = nsamp) 

plot( 1 , 1 , type="n" , xlab="Treatment" , ylab="probability" , 
      xlim=c(0,1) , ylim=c(0,1) , xaxp=c(0,1,1) , yaxp=c(0,1,2) )

Treat <- 0:1 #range of treatment values to iterate over
gender_id <- 1 #value for gender
country_id <- 6 #value for country (i.e. USA)
for ( s in 1:nsamp ) {
  p <- lapply(1:length(post), function (x) if(length(dim(post[[x]])) == 2) post[[x]][s,] else  post[[x]][s])
  names(p) <- names(post)
  ak <- as.numeric(p$cutpoints[1:4])
  phi <- p$pT*Treat + p$pG[gender_id] + p$pC[country_id]
  pk <- pordlogit( 1:4 , a=ak , phi=phi )
  for ( i in 1:4 )
    lines( Treat , pk[,i] , col=col.alpha(rangi2,0.1) )
}
title("Effect of Treatment on Cumulative Probability of Response")
abline(h = 1, lty = 2)
abline(h = 0, lty = 2)

text(x = .8, y= .05, labels = "Strongly Disagree")
text(x = .8, y= .2, labels = "Somewhat Disagree")
text(x = .8, y= .5, labels = "Neither Agree nor Disagree")
text(x = .8, y= .75, labels = "Somewhat Agree")
text(x = .8, y= .95, labels = "Strongly Agree")


#computing difference between control and treatment

nsamp <- 5000 #number of samples to extract from the posterior and subsequently plot predictions of
post <- extract.samples( m6 , n = nsamp) 
Treat <- 0:1 #range of treatment values to iterate over
gender_id <- 1 #value for gender
country_id <- 6 #value for country (i.e. USA)
diffs <- matrix(nrow = 4, ncol = nsamp)
for ( s in 1:nsamp ) {
  p <- lapply(1:length(post), function (x) if(length(dim(post[[x]])) == 2) post[[x]][s,] else  post[[x]][s])
  names(p) <- names(post)
  ak <- as.numeric(p$cutpoints[1:4])
  phi <- p$pT*Treat + p$pG[gender_id] + p$pC[country_id]
  pk <- pordlogit( 1:4 , a=ak , phi=phi )
  for ( i in 1:4 )
    diffs[i,s] <- pk[2,i] - pk[1,i]
}

dens(diffs[1,], lwd = 2)
text(x = mean(diffs[1,]), y = 2, labels = paste0("1 (", sum(diffs[1,] < 0) / length(diffs[1,]), ")"), lwd = 2)
for(i in 2:4){
  dens(diffs[i,], add = T, col = i, lwd = 2)
  text(x = mean(diffs[i,]), y = i*2, labels = paste0(i, " (", sum(diffs[i,] < 0) / length(diffs[i,]), ")"), col = i, lwd = 2)
}


#testing the effects of age given treatment


nsamp <- 200 #number of samples to extract from the posterior and subsequently plot
post <- extract.samples( m5 , n = nsamp) 

age.seq <- seq(from = -1.5, to = 2, by = .2)
par(mfrow = c(6,3))

for(ages in 1:length(age.seq)){
plot( 1 , 1 , type="n" , xlab="Treatment" , ylab="probability" , 
      xlim=c(0,1) , ylim=c(0,1) , xaxp=c(0,1,1) , yaxp=c(0,1,2) )

Treat <- 0:1 #range of treatment values to iterate over
gender_id <- 1 #value for gender (i.e. woman)
country_id <- 6 #value for country (i.e. USA)
Age_s <- age.seq[ages]
for ( s in 1:nsamp ) {
  p <- lapply(1:length(post), function (x) if(length(dim(post[[x]])) == 2) post[[x]][s,] else  post[[x]][s])
  names(p) <- names(post)
  ak <- as.numeric(p$cutpoints[1:4])
  phi <- p$pT*Treat + p$pG[gender_id] + p$pC[country_id] + p$pA*Age_s
  pk <- pordlogit( 1:4 , a=ak , phi=phi )
  for ( i in 1:4 )
    lines( Treat , pk[,i] , col=col.alpha(rangi2,0.1) )
}
title(paste0("Age_s = ", Age_s))
abline(h = 1, lty = 2)
abline(h = 0, lty = 2)

text(x = .8, y= .05, labels = "StD")
text(x = .8, y= .2, labels = "SoD")
text(x = .8, y= .5, labels = "NAND")
text(x = .8, y= .75, labels = "SoA")
text(x = .8, y= .95, labels = "StA")
}

# age diffs distribution


nsamp <- 5000 #number of samples to extract from the posterior and subsequently plot
post <- extract.samples( m5 , n = nsamp) 

age.seq <- seq(from = -1.5, to = 2, by = .2)
age.seq[8] <- -0.1
par(mfrow = c(6,3))

for(ages in 1:length(age.seq)){
  Treat <- 0:1 #range of treatment values to iterate over
  gender_id <- 1 #value for gender (i.e. woman)
  country_id <- 6 #value for country (i.e. USA)
  Age_s <- age.seq[ages]
  diffs <- matrix(nrow = 4, ncol = nsamp)
  for ( s in 1:nsamp ) {
    p <- lapply(1:length(post), function (x) if(length(dim(post[[x]])) == 2) post[[x]][s,] else  post[[x]][s])
    names(p) <- names(post)
    ak <- as.numeric(p$cutpoints[1:4])
    phi <- p$pT*Treat + p$pG[gender_id] + p$pC[country_id] + p$pA*Age_s
    pk <- pordlogit( 1:4 , a=ak , phi=phi )
    for ( i in 1:4 )
      diffs[i,s] <- pk[2,i] - pk[1,i]
  }

  dens(diffs[1,], lwd = 2)
  title(paste0("Age_s = ", Age_s))
  text(x = mean(diffs[1,]), y = 2, labels = paste0("1 (", sum(diffs[1,] < 0) / length(diffs[1,]), ")"), lwd = 2)
  for(i in 2:4){
    dens(diffs[i,], add = T, col = i, lwd = 2)
    text(x = mean(diffs[i,]), y = i*2, labels = paste0(i, " (", sum(diffs[i,] < 0) / length(diffs[i,]), ")"), col = i, lwd = 2)
}
}



#probing m7



#testing the effects of age given treatment


nsamp <- 200 #number of samples to extract from the posterior and subsequently plot
post <- extract.samples( m7 , n = nsamp) 

age.seq <- seq(from = -1.5, to = 2, by = .2)
age.seq[8] <- -0.1
par(mfrow = c(6,3))

for(ages in 1:length(age.seq)){
  plot( 1 , 1 , type="n" , xlab="Treatment" , ylab="probability" , 
        xlim=c(0,1) , ylim=c(0,1) , xaxp=c(0,1,1) , yaxp=c(0,1,2) )
  
  Treat <- 0:1 #range of treatment values to iterate over
  gender_id <- 1 #value for gender (i.e. woman)
  country_id <- 6 #value for country (i.e. USA)
  Age_s <- age.seq[ages]
  for ( s in 1:nsamp ) {
    p <- lapply(1:length(post), function (x) if(length(dim(post[[x]])) == 2) post[[x]][s,] else  post[[x]][s])
    names(p) <- names(post)
    ak <- as.numeric(p$cutpoints[1:4])
    phi <- p$pT*Treat + p$pG[gender_id] + p$pC[country_id] + p$pA*Age_s + p$pTA * Treat * Age_s
    pk <- pordlogit( 1:4 , a=ak , phi=phi )
    for ( i in 1:4 )
      lines( Treat , pk[,i] , col=col.alpha(rangi2,0.1) )
  }
  title(paste0("Age_s = ", Age_s))
  abline(h = 1, lty = 2)
  abline(h = 0, lty = 2)
  
  text(x = .8, y= .05, labels = "StD")
  text(x = .8, y= .2, labels = "SoD")
  text(x = .8, y= .5, labels = "NAND")
  text(x = .8, y= .75, labels = "SoA")
  text(x = .8, y= .95, labels = "StA")
}

# age diffs distribution


nsamp <- 5000 #number of samples to extract from the posterior and subsequently plot
post <- extract.samples( m7 , n = nsamp) 

age.seq <- seq(from = -1.5, to = 2, by = .2)
age.seq[8] <- -0.1
par(mfrow = c(6,3))


for(ages in 1:length(age.seq)){
  Treat <- 0:1 #range of treatment values to iterate over
  gender_id <- 1 #value for gender (i.e. woman)
  country_id <- 6 #value for country (i.e. USA)
  Age_s <- age.seq[ages]
  diffs <- matrix(nrow = 4, ncol = nsamp)
  for ( s in 1:nsamp ) {
    p <- lapply(1:length(post), function (x) if(length(dim(post[[x]])) == 2) post[[x]][s,] else  post[[x]][s])
    names(p) <- names(post)
    ak <- as.numeric(p$cutpoints[1:4])
    phi <- p$pT*Treat + p$pG[gender_id] + p$pC[country_id] + p$pA*Age_s + p$pTA * Treat * Age_s
    pk <- pordlogit( 1:4 , a=ak , phi=phi )
    for ( i in 1:4 )
      diffs[i,s] <- pk[2,i] - pk[1,i]
  }
  
  dens(diffs[1,], lwd = 2)
  title(paste0("Age_s = ", Age_s))
  text(x = mean(diffs[1,]), y = 2, labels = paste0("1 (", sum(diffs[1,] < 0) / length(diffs[1,]), ")"), lwd = 2)
  for(i in 2:4){
    dens(diffs[i,], add = T, col = i, lwd = 2)
    text(x = mean(diffs[i,]), y = i*2, labels = paste0(i, " (", sum(diffs[i,] < 0) / length(diffs[i,]), ")"), col = i, lwd = 2)
  }
}


## MAKING ANIMATIONS FROM THE PREVIOUS CODE ##
## MAKING ANIMATIONS FROM THE PREVIOUS CODE ##
## MAKING ANIMATIONS FROM THE PREVIOUS CODE ##
## MAKING ANIMATIONS FROM THE PREVIOUS CODE ##
## MAKING ANIMATIONS FROM THE PREVIOUS CODE ##



png(paste0("testplot", i, ".png"))
dev.off()


#testing the effects of age given treatment


nsamp <- 500 #number of samples to extract from the posterior and subsequently plot
post <- extract.samples( m5 , n = nsamp) 

age.seq <- seq(from = -1.5, to = 2, by = .1)
age.seq[15] <- -0.1
par(mfrow = c(1,1))
setwd("C:/Users/Nikolai/Desktop/AttitudeAnimations")
for(ages in 1:length(age.seq)){
  print(ages)
  png(paste0("m5Changes", ages, ".png"))
  plot( 1 , 1 , type="n" , xlab="Treatment" , ylab="probability" , 
        xlim=c(0,1) , ylim=c(0,1) , xaxp=c(0,1,1) , yaxp=c(0,1,2) )
  
  Treat <- 0:1 #range of treatment values to iterate over
  gender_id <- 1 #value for gender (i.e. woman)
  country_id <- 6 #value for country (i.e. USA)
  Age_s <- age.seq[ages]
  for ( s in 1:nsamp ) {
    p <- lapply(1:length(post), function (x) if(length(dim(post[[x]])) == 2) post[[x]][s,] else  post[[x]][s])
    names(p) <- names(post)
    ak <- as.numeric(p$cutpoints[1:4])
    phi <- p$pT*Treat + p$pG[gender_id] + p$pC[country_id] + p$pA*Age_s
    pk <- pordlogit( 1:4 , a=ak , phi=phi )
    for ( i in 1:4 )
      lines( Treat , pk[,i] , col=col.alpha(rangi2,0.05) )
  }
  title(paste0("Age_s = ", Age_s))
  abline(h = 1, lty = 2)
  abline(h = 0, lty = 2)
  
  text(x = .8, y= .05, labels = "StD")
  text(x = .8, y= .2, labels = "SoD")
  text(x = .8, y= .5, labels = "NAND")
  text(x = .8, y= .75, labels = "SoA")
  text(x = .8, y= .95, labels = "StA")
  dev.off()
}

# age diffs distribution


nsamp <- 10000 #number of samples to extract from the posterior and subsequently plot
post <- extract.samples( m5 , n = nsamp) 

age.seq <- seq(from = -1.5, to = 2, by = .1)
age.seq[15] <- -0.1
par(mfrow = c(1,1))

for(ages in 1:length(age.seq)){
  print(ages)
  png(paste0("m5Dist", ages, ".png"))
  Treat <- 0:1 #range of treatment values to iterate over
  gender_id <- 1 #value for gender (i.e. woman)
  country_id <- 6 #value for country (i.e. USA)
  Age_s <- age.seq[ages]
  diffs <- matrix(nrow = 4, ncol = nsamp)
  for ( s in 1:nsamp ) {
    p <- lapply(1:length(post), function (x) if(length(dim(post[[x]])) == 2) post[[x]][s,] else  post[[x]][s])
    names(p) <- names(post)
    ak <- as.numeric(p$cutpoints[1:4])
    phi <- p$pT*Treat + p$pG[gender_id] + p$pC[country_id] + p$pA*Age_s
    pk <- pordlogit( 1:4 , a=ak , phi=phi )
    for ( i in 1:4 )
      diffs[i,s] <- pk[2,i] - pk[1,i]
  }
  
  dens(diffs[1,], lwd = 2, xlim = c(-.1, .1), ylim = c(0,45))
  title(paste0("Age_s = ", Age_s))
  text(x = mean(diffs[1,]), y = 2, labels = paste0("1 (", sum(diffs[1,] < 0) / length(diffs[1,]), ")"), lwd = 2)
  for(i in 2:4){
    dens(diffs[i,], add = T, col = i, lwd = 2)
    text(x = mean(diffs[i,]), y = i*2, labels = paste0(i, " (", sum(diffs[i,] < 0) / length(diffs[i,]), ")"), col = i, lwd = 2)
  }
  dev.off()
}



#probing m7



#testing the effects of age given treatment


nsamp <- 500 #number of samples to extract from the posterior and subsequently plot
post <- extract.samples( m7 , n = nsamp) 

age.seq <- seq(from = -1.5, to = 2, by = .1)
age.seq[15] <- -0.1
par(mfrow = c(1,1))

for(ages in 1:length(age.seq)){
  print(ages)
  png(paste0("m7Changes", ages, ".png"))
  plot( 1 , 1 , type="n" , xlab="Treatment" , ylab="probability" , 
        xlim=c(0,1) , ylim=c(0,1) , xaxp=c(0,1,1) , yaxp=c(0,1,2) )
  
  Treat <- 0:1 #range of treatment values to iterate over
  gender_id <- 1 #value for gender (i.e. woman)
  country_id <- 6 #value for country (i.e. USA)
  Age_s <- age.seq[ages]
  for ( s in 1:nsamp ) {
    p <- lapply(1:length(post), function (x) if(length(dim(post[[x]])) == 2) post[[x]][s,] else  post[[x]][s])
    names(p) <- names(post)
    ak <- as.numeric(p$cutpoints[1:4])
    phi <- p$pT*Treat + p$pG[gender_id] + p$pC[country_id] + p$pA*Age_s + p$pTA * Treat * Age_s
    pk <- pordlogit( 1:4 , a=ak , phi=phi )
    for ( i in 1:4 )
      lines( Treat , pk[,i] , col=col.alpha(rangi2,0.05) )
  }
  title(paste0("Age_s = ", Age_s))
  abline(h = 1, lty = 2)
  abline(h = 0, lty = 2)
  
  text(x = .8, y= .05, labels = "StD")
  text(x = .8, y= .2, labels = "SoD")
  text(x = .8, y= .5, labels = "NAND")
  text(x = .8, y= .75, labels = "SoA")
  text(x = .8, y= .95, labels = "StA")
  dev.off()
}

# age diffs distribution


nsamp <- 10000 #number of samples to extract from the posterior and subsequently plot
post <- extract.samples( m7 , n = nsamp) 

age.seq <- seq(from = -1.5, to = 2, by = .1)
age.seq[15] <- -0.1
par(mfrow = c(1,1))


for(ages in 1:length(age.seq)){
  print(ages)
  png(paste0("m7Dist", ages, ".png"))
  Treat <- 0:1 #range of treatment values to iterate over
  gender_id <- 1 #value for gender (i.e. woman)
  country_id <- 6 #value for country (i.e. USA)
  Age_s <- age.seq[ages]
  diffs <- matrix(nrow = 4, ncol = nsamp)
  for ( s in 1:nsamp ) {
    p <- lapply(1:length(post), function (x) if(length(dim(post[[x]])) == 2) post[[x]][s,] else  post[[x]][s])
    names(p) <- names(post)
    ak <- as.numeric(p$cutpoints[1:4])
    phi <- p$pT*Treat + p$pG[gender_id] + p$pC[country_id] + p$pA*Age_s + p$pTA * Treat * Age_s
    pk <- pordlogit( 1:4 , a=ak , phi=phi )
    for ( i in 1:4 )
      diffs[i,s] <- pk[2,i] - pk[1,i]
  }
  
  dens(diffs[1,], lwd = 2, xlim = c(-.1, .1), ylim = c(0,45))
  title(paste0("Age_s = ", Age_s))
  text(x = mean(diffs[1,]), y = 2, labels = paste0("1 (", sum(diffs[1,] < 0) / length(diffs[1,]), ")"), lwd = 2)
  for(i in 2:4){
    dens(diffs[i,], add = T, col = i, lwd = 2)
    text(x = mean(diffs[i,]), y = i*2, labels = paste0(i, " (", sum(diffs[i,] < 0) / length(diffs[i,]), ")"), col = i, lwd = 2)
  }
  dev.off()
}

#compare m7 pTA prior and posterior


dens(post$pTA, xlim = c(-10, 10), col=col.alpha(rangi2,1)); dens(rnorm(n=1e5, mean = 0, sd = 5), add = T)
legend(7,4.5, c("Prior","Posterior"), lty=c(1,1), lwd=c(2.5,2.5),col=c(1,col.alpha(rangi2,1)))
title("Prior and Posterior Distributions for the Interaction pTA between Age and Treatment (m7)")


save.image("C:/Users/Nikolai/Dropbox/R Files/attitude_MG.RData")



