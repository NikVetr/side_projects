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
d <- read.csv("Data/CombinedData.csv")
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

d2 <- d[, c("Treat", "X2DaysPork", "X2DaysBeef", "X2DaysDairy", "X2DaysEggs", "X2DaysBird", "X2DaysSea", "NowVs4MonPast",
            "NowVs4MonFut", "MeatUniqueSmart", "LessMeatGood", "MeatReplace", "Gender", "Country", "Age.s")]

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


# --------- #Chicken Data# ---------- #


m0cc.c <- map2stan(
  alist(
    X2DaysBird ~ dzipois( p , lambda ),
    logit(p) <- pI, 
    log(lambda) <- lI, 
    c(pI) ~ dnorm(-2,0.5),
    c(lI) ~ dnorm(1,1)
  ) ,
  data= list(X2DaysBird = dcc$X2DaysBird),
  iter = 25000, warmup = 5000
)
precis(m0cc.c)

m1cc.c <- map2stan(
  alist(
    X2DaysBird ~ dzipois( p , lambda ),
    logit(p) <- pI + pT*Treat,
    log(lambda) <- lI + lT*Treat,
    pI ~ dnorm(-2, 0.5),
    pT ~ dnorm(0,0.5),
    lI ~ dnorm(1,1),
    lT ~ dnorm(0,1)
  ) ,
  data= list(Treat = dcc$Treat, X2DaysBird = dcc$X2DaysBird),
  iter = 25000, warmup = 5000
)
precis(m1cc.c)

m2cc.c <- map2stan( #no interactions, just main effects
  alist(
    X2DaysBird ~ dzipois( p , lambda ),
    logit(p) <- pI + pT*Treat + pG[gender_id] + pA*Age_s + pC[country_id],
    log(lambda) <- lI + lT*Treat + lG[gender_id] + lA*Age_s + lC[country_id],
    
    pI ~ dnorm(-2, 0.5),
    c(pT, pA) ~ dnorm(0,0.5),
    pC[country_id] ~ dnorm(0,sigCp),
    sigCp ~ dcauchy(0,0.5),
    pG[gender_id] ~ dnorm(0,sigGp),
    sigGp ~ dcauchy(0,0.5),
    
    lI ~ dnorm(1,1),
    c(lT, lA) ~ dnorm(0,1),
    lG[gender_id] ~ dnorm(0,sigGl),
    sigGl ~ dcauchy(0,0.5),
    lC[country_id] ~ dnorm(0,sigCl),
    sigCl ~ dcauchy(0,0.5)
  ) ,
  data= list(Treat = dcc$Treat, X2DaysBird = dcc$X2DaysBird, gender_id = dcc$gender_id,
             Age_s = dcc$Age.s, country_id = dcc$country_id),
  iter = 25000, warmup = 5000
)
precis(m2cc.c, depth = 2)

m3cc.c <- map2stan( #no interactions, just Treatment and Gender effects
  alist(
    X2DaysBird ~ dzipois( p , lambda ),
    logit(p) <- pI + pT*Treat + pG[gender_id],
    log(lambda) <- lI + lT*Treat + lG[gender_id],
    
    pI ~ dnorm(-2, 0.5),
    pT ~ dnorm(0,0.5),
    pG[gender_id] ~ dnorm(0,sigGp),
    sigGp ~ dcauchy(0,0.5),
    
    lI ~ dnorm(1,1),
    lT ~ dnorm(0,1),
    lG[gender_id] ~ dnorm(0,sigGl),
    sigGl ~ dcauchy(0,0.5)
  ) ,
  data= list(Treat = dcc$Treat, X2DaysBird = dcc$X2DaysBird, gender_id = dcc$gender_id),
  iter = 25000, warmup = 5000
)
precis(m3cc.c, depth = 2)

m4cc.c <- map2stan( #no interactions, just Treatment and Country effects
  alist(
    X2DaysBird ~ dzipois( p , lambda ),
    logit(p) <- pI + pT*Treat + pC[country_id],
    log(lambda) <- lI + lT*Treat + lC[country_id],
    
    pI ~ dnorm(-2, 0.5),
    pT ~ dnorm(0,0.5),
    pC[country_id] ~ dnorm(0,sigCp),
    sigCp ~ dcauchy(0,0.5),
    
    lI ~ dnorm(1,1),
    lT ~ dnorm(0,1),
    lC[country_id] ~ dnorm(0,sigCl),
    sigCl ~ dcauchy(0,0.5)
  ) ,
  data= list(Treat = dcc$Treat, X2DaysBird = dcc$X2DaysBird, country_id = dcc$country_id),
  iter = 25000, warmup = 5000
)
precis(m4cc.c, depth = 2)

m5cc.c <- map2stan( #no interactions, just Treatment and Age effects
  alist(
    X2DaysBird ~ dzipois( p , lambda ),
    logit(p) <- pI + pT*Treat + pA*Age_s,
    log(lambda) <- lI + lT*Treat + lA*Age_s,
    
    pI ~ dnorm(-2,0.5),
    c(pT, pA) ~ dnorm(0,0.5),
    
    lI ~ dnorm(1,1),
    c(lT, lA) ~ dnorm(0,1)
  ) ,
  data= list(Treat = dcc$Treat, X2DaysBird = dcc$X2DaysBird, Age_s = dcc$Age.s),
  iter = 25000, warmup = 5000
)
precis(m5cc.c, depth = 2)

m6cc.c <- map2stan( #A-T interaction and main effects
  alist(
    X2DaysBird ~ dzipois( p , lambda ),
    logit(p) <- pI + pT*Treat + pAT*Treat*Age_s + pG[gender_id] + pA*Age_s + pC[country_id],
    log(lambda) <- lI + lT*Treat + lAT*Age_s*Treat + lG[gender_id] + lA*Age_s + lC[country_id],
    
    
    pI ~ dnorm(-2, 0.5),
    c(pT, pA, pAT) ~ dnorm(0,0.5),
    pC[country_id] ~ dnorm(0,sigCp),
    sigCp ~ dcauchy(0,0.5),
    pG[gender_id] ~ dnorm(0,sigGp),
    sigGp ~ dcauchy(0,0.5),
    
    lI ~ dnorm(1,1),
    c(lT, lA, lAT) ~ dnorm(0,1),
    lG[gender_id] ~ dnorm(0,sigGl),
    sigGl ~ dcauchy(0,0.5),
    lC[country_id] ~ dnorm(0,sigCl),
    sigCl ~ dcauchy(0,0.5)
  ) ,
  data= list(Treat = dcc$Treat, X2DaysBird = dcc$X2DaysBird, gender_id = dcc$gender_id,
             Age_s = dcc$Age.s, country_id = dcc$country_id),
  iter = 30000, warmup = 6000
)
precis(m6cc.c, depth = 2)

m7cc.c <- map2stan( #A-T-G interaction and main effects
  alist(
    X2DaysBird ~ dzipois( p , lambda ),
    logit(p) <- pI + pATG*Treat + pG[gender_id] + pA*Age_s + pC[country_id],
    pATG <- pT + pAT * Age_s + pGT[gender_id],
    log(lambda) <- lI + lATG*Treat + lG[gender_id] + lA*Age_s + lC[country_id],
    lATG <- lT + lAT * Age_s + lGT[gender_id],
    
    
    pI ~ dnorm(-2, 0.5),
    c(pT, pA, pAT) ~ dnorm(0,0.5),
    pC[country_id] ~ dnorm(0,sigCp),
    sigCp ~ dcauchy(0,0.5),
    pG[gender_id] ~ dnorm(0,sigGp),
    sigGp ~ dcauchy(0,0.5),
    pGT[gender_id] ~ dnorm(0,sigGTp),
    sigGTp ~ dcauchy(0,0.5),
    
    lI ~ dnorm(1,1),
    c(lT, lA, lAT) ~ dnorm(0,1),
    lG[gender_id] ~ dnorm(0,sigGl),
    sigGl ~ dcauchy(0,0.5),
    lGT[gender_id] ~ dnorm(0,sigGTl),
    sigGTl ~ dcauchy(0,0.5),
    lC[country_id] ~ dnorm(0,sigCl),
    sigCl ~ dcauchy(0,0.5)
  ) ,
  data= list(Treat = dcc$Treat, X2DaysBird = dcc$X2DaysBird, gender_id = dcc$gender_id,
             Age_s = dcc$Age.s, country_id = dcc$country_id),
  iter = 25000, warmup = 5000
)
precis(m7cc.c, depth = 2)

m8cc.c <- map2stan( #A-T-G-C interaction and main effects
  alist(
    X2DaysBird ~ dzipois( p , lambda ),
    logit(p) <- pI + pATGC*Treat + pG[gender_id] + pA*Age_s + pC[country_id],
    pATGC <- pT + pAT * Age_s + pGT[gender_id] + pCT[country_id],
    log(lambda) <- lI + lATGC*Treat + lG[gender_id] + lA*Age_s + lC[country_id],
    lATGC <- lT + lAT * Age_s + lGT[gender_id] + lCT[country_id],
    
    pI ~ dnorm(-2, 0.5),
    c(pT, pA, pAT) ~ dnorm(0,0.5),
    pC[country_id] ~ dnorm(0,sigCp),
    sigCp ~ dcauchy(0,0.5),
    pG[gender_id] ~ dnorm(0,sigGp),
    sigGp ~ dcauchy(0,0.5),
    pGT[gender_id] ~ dnorm(0,sigGTp),
    sigGTp ~ dcauchy(0,0.5),
    pCT[country_id] ~ dnorm(0,sigCTp),
    sigCTp ~ dcauchy(0,0.5),
    
    lI ~ dnorm(1,1),
    c(lT, lA, lAT) ~ dnorm(0,1),
    lG[gender_id] ~ dnorm(0,sigGl),
    sigGl ~ dcauchy(0,0.5),
    lGT[gender_id] ~ dnorm(0,sigGTl),
    sigGTl ~ dcauchy(0,0.5),
    lCT[country_id] ~ dnorm(0,sigCTl),
    sigCTl ~ dcauchy(0,0.5),
    lC[country_id] ~ dnorm(0,sigCl),
    sigCl ~ dcauchy(0,0.5)
  ) ,
  data= list(Treat = dcc$Treat, X2DaysBird = dcc$X2DaysBird, gender_id = dcc$gender_id,
             Age_s = dcc$Age.s, country_id = dcc$country_id),
  iter = 25000, warmup = 5000
)
precis(m8cc.c, depth = 2)

m9cc.c <- map2stan( 
  alist(
    X2DaysBird ~ dzipois( p , lambda ),
    logit(p) <- pI + pT*Treat,
    log(lambda) <- lI + lATG*Treat + lG[gender_id] + lA*Age_s + lC[country_id],
    lATG <- lT + lAT * Age_s + lGT[gender_id],
    
    pI ~ dnorm(-2, 0.5),
    pT ~ dnorm(0,0.5),
    
    lI ~ dnorm(1,1),
    c(lT, lA, lAT) ~ dnorm(0,1),
    lG[gender_id] ~ dnorm(0,sigGl),
    sigGl ~ dcauchy(0,0.5),
    lGT[gender_id] ~ dnorm(0,sigGTl),
    sigGTl ~ dcauchy(0,0.5),
    lC[country_id] ~ dnorm(0,sigCl),
    sigCl ~ dcauchy(0,0.5)
  ) ,
  data= list(Treat = dcc$Treat, X2DaysBird = dcc$X2DaysBird, gender_id = dcc$gender_id,
             Age_s = dcc$Age.s, country_id = dcc$country_id),
  iter = 25000, warmup = 5000
)
precis(m9cc.c, depth = 2)

m10cc.c <- map2stan( #A-T-G interaction on lambda and main effects; no effect on p of treatment
  alist(
    X2DaysBird ~ dzipois( p , lambda ),
    logit(p) <- pI,
    log(lambda) <- lI + lATG*Treat + lG[gender_id] + lA*Age_s + lC[country_id],
    lATG <- lT + lAT * Age_s + lGT[gender_id],
    
    pI ~ dnorm(-2, 0.5),
    
    lI ~ dnorm(1,1),
    c(lT, lA, lAT) ~ dnorm(0,1),
    lG[gender_id] ~ dnorm(0,sigGl),
    sigGl ~ dcauchy(0,0.5),
    lGT[gender_id] ~ dnorm(0,sigGTl),
    sigGTl ~ dcauchy(0,0.5),
    lC[country_id] ~ dnorm(0,sigCl),
    sigCl ~ dcauchy(0,0.5)
  ) ,
  data= list(Treat = dcc$Treat, X2DaysBird = dcc$X2DaysBird, gender_id = dcc$gender_id,
             Age_s = dcc$Age.s, country_id = dcc$country_id),
  iter = 25000, warmup = 5000
)
precis(m10cc.c, depth = 2)

m12cc.c <- map2stan( #A-T-G-C interaction on lambda and main effects; no effect on p of treatment
  alist(
    X2DaysBird ~ dzipois( p , lambda ),
    logit(p) <- pI,
    log(lambda) <- lI + lATGC*Treat + lG[gender_id] + lA*Age_s + lC[country_id],
    lATGC <- lT + lAT * Age_s + lGT[gender_id] + lCT[country_id],
    
    pI ~ dnorm(-2, 0.5),
    lI ~ dnorm(1,1),
    c(lT, lA, lAT) ~ dnorm(0,1),
    lG[gender_id] ~ dnorm(0,sigGl),
    sigGl ~ dcauchy(0,0.5),
    lGT[gender_id] ~ dnorm(0,sigGTl),
    sigGTl ~ dcauchy(0,0.5),
    lCT[country_id] ~ dnorm(0,sigCTl),
    sigCTl ~ dcauchy(0,0.5),
    lC[country_id] ~ dnorm(0,sigCl),
    sigCl ~ dcauchy(0,0.5)
  ) ,
  data= list(Treat = dcc$Treat, X2DaysBird = dcc$X2DaysBird, gender_id = dcc$gender_id,
             Age_s = dcc$Age.s, country_id = dcc$country_id),
  iter = 25000, warmup = 5000
)
precis(m12cc.c, depth = 2)

m13cc.c <- map2stan( #A-T interaction on lambda and main effects; no effect on p of treatment
  alist(
    X2DaysBird ~ dzipois( p , lambda ),
    logit(p) <- pI,
    log(lambda) <- lI + lT*Treat + lAT*Age_s*Treat + lG[gender_id] + lA*Age_s + lC[country_id],
    
    pI ~ dnorm(-2, 0.5),
    lI ~ dnorm(1,1),
    c(lT, lA, lAT) ~ dnorm(0,1),
    lG[gender_id] ~ dnorm(0,sigGl),
    sigGl ~ dcauchy(0,0.5),
    lC[country_id] ~ dnorm(0,sigCl),
    sigCl ~ dcauchy(0,0.5)
  ) ,
  data= list(Treat = dcc$Treat, X2DaysBird = dcc$X2DaysBird, gender_id = dcc$gender_id,
             Age_s = dcc$Age.s, country_id = dcc$country_id),
  iter = 25000, warmup = 5000
)
precis(m13cc.c, depth = 2)

m14cc.c <- map2stan( #A-G interaction on lambda and main effects; no effect on p of treatment
  alist(
    X2DaysBird ~ dzipois( p , lambda ),
    logit(p) <- pI,
    log(lambda) <- lI + lTG*Treat + lG[gender_id] + lA*Age_s + lC[country_id],
    lTG <- lT + lGT[gender_id],
    
    pI ~ dnorm(-2, 0.5),
    lI ~ dnorm(1,1),    
    c(lT, lA) ~ dnorm(0,1),
    lG[gender_id] ~ dnorm(0,sigGl),
    sigGl ~ dcauchy(0,0.5),
    lGT[gender_id] ~ dnorm(0,sigGTl),
    sigGTl ~ dcauchy(0,0.5),
    lC[country_id] ~ dnorm(0,sigCl),
    sigCl ~ dcauchy(0,0.5)
  ) ,
  data= list(Treat = dcc$Treat, X2DaysBird = dcc$X2DaysBird, gender_id = dcc$gender_id,
             Age_s = dcc$Age.s, country_id = dcc$country_id),
  iter = 25000, warmup = 5000
)
precis(m14cc.c, depth = 2)


# ----- #NO TREATMENT EFFECTS ----- #

m2cc.c.nt <- map2stan( #no interactions, just main effects
  alist(
    X2DaysBird ~ dzipois( p , lambda ),
    logit(p) <- pI + pG[gender_id] + pA*Age_s + pC[country_id],
    log(lambda) <- lI + lG[gender_id] + lA*Age_s + lC[country_id],
    
    pI ~ dnorm(-2, 0.5),
    pA ~ dnorm(0,0.5),
    pC[country_id] ~ dnorm(0,sigCp),
    sigCp ~ dcauchy(0,0.5),
    pG[gender_id] ~ dnorm(0,sigGp),
    sigGp ~ dcauchy(0,0.5),
    
    lI ~ dnorm(1,1),
    lA ~ dnorm(0,1),
    lG[gender_id] ~ dnorm(0,sigGl),
    sigGl ~ dcauchy(0,0.5),
    lC[country_id] ~ dnorm(0,sigCl),
    sigCl ~ dcauchy(0,0.5)
  ) ,
  data= list(X2DaysBird = dcc$X2DaysBird, gender_id = dcc$gender_id,
             Age_s = dcc$Age.s, country_id = dcc$country_id),
  iter = 25000, warmup = 5000
)
precis(m2cc.c.nt, depth = 2)

m3cc.c.nt <- map2stan( #no interactions, just AlrVeg, Treatment, and Gender effects
  alist(
    X2DaysBird ~ dzipois( p , lambda ),
    logit(p) <- pI + pG[gender_id],
    log(lambda) <- lI + lG[gender_id],
    
    pI ~ dnorm(-2, 0.5),
    pG[gender_id] ~ dnorm(0,sigGp),
    sigGp ~ dcauchy(0,0.5),
    
    lI ~ dnorm(1,1),
    lG[gender_id] ~ dnorm(0,sigGl),
    sigGl ~ dcauchy(0,0.5)
  ) ,
  data= list(X2DaysBird = dcc$X2DaysBird, gender_id = dcc$gender_id),
  iter = 30000, warmup = 6000
)
precis(m3cc.c.nt, depth = 2)

m4cc.c.nt <- map2stan( #no interactions, just AlrVeg, Treatment, and Country effects
  alist(
    X2DaysBird ~ dzipois( p , lambda ),
    logit(p) <- pI + pC[country_id],
    log(lambda) <- lI + lC[country_id],
    
    pI ~ dnorm(-2, 0.5),
    pC[country_id] ~ dnorm(0,sigCp),
    sigCp ~ dcauchy(0,0.5),
    
    lI ~ dnorm(1,1),
    lC[country_id] ~ dnorm(0,sigCl),
    sigCl ~ dcauchy(0,0.5)
  ) ,
  data= list(X2DaysBird = dcc$X2DaysBird, country_id = dcc$country_id),
  iter = 25000, warmup = 5000
)
precis(m4cc.c.nt, depth = 2)

m5cc.c.nt <- map2stan( #no interactions, just AlrVeg, Treatment, and Age effects
  alist(
    X2DaysBird ~ dzipois( p , lambda ),
    logit(p) <- pI + pA*Age_s,
    log(lambda) <- lI + lA*Age_s,
    
    pI ~ dnorm(-2, 0.5),
    lI ~ dnorm(1,1),
    pA ~ dnorm(0,0.5),
    lA ~ dnorm(0,1)
  ) ,
  data= list(X2DaysBird = dcc$X2DaysBird, Age_s = dcc$Age.s),
  iter = 25000, warmup = 5000
)
precis(m5cc.c.nt, depth = 2)

m7cc.c.nt <- map2stan( #A-T-G interaction and main effects
  alist(
    X2DaysBird ~ dzipois( p , lambda ),
    logit(p) <- pI + pG[gender_id] + pA*Age_s + pC[country_id],
    log(lambda) <- lI + lG[gender_id] + lA*Age_s + lC[country_id],
    
    pI ~ dnorm(-2, 0.5),
    lI ~ dnorm(1,1),
    pA ~ dnorm(0,0.5),
    pC[country_id] ~ dnorm(0,sigCp),
    sigCp ~ dcauchy(0,0.5),
    pG[gender_id] ~ dnorm(0,sigGp),
    sigGp ~ dcauchy(0,0.5),
    
    lI ~ dnorm(1,1),
    lA ~ dnorm(0,1),
    lG[gender_id] ~ dnorm(0,sigGl),
    sigGl ~ dcauchy(0,0.5),
    lC[country_id] ~ dnorm(0,sigCl),
    sigCl ~ dcauchy(0,0.5)
  ) ,
  data= list(X2DaysBird = dcc$X2DaysBird, gender_id = dcc$gender_id,
             Age_s = dcc$Age.s, country_id = dcc$country_id),
  iter = 25000, warmup = 5000
)
precis(m7cc.c.nt, depth = 2)

m11cc.c.nt <- map2stan( #no treatment effects
  alist(
    X2DaysBird ~ dzipois( p , lambda ),
    logit(p) <- pI,
    log(lambda) <- lI  + lG[gender_id] + lA*Age_s + lC[country_id],
    
    pI ~ dnorm(-2, 0.5),
    lI ~ dnorm(1,1),
    lA ~ dnorm(0,1),
    lG[gender_id] ~ dnorm(0,sigGl),
    sigGl ~ dcauchy(0,0.5),
    lC[country_id] ~ dnorm(0,sigCl),
    sigCl ~ dcauchy(0,0.5)
  ) ,
  data= list(X2DaysBird = dcc$X2DaysBird, gender_id = dcc$gender_id,
             Age_s = dcc$Age.s, country_id = dcc$country_id),
  iter = 25000, warmup = 5000
)
precis(m11cc.c.nt, depth = 2)

save.image("C:/Users/Nikolai/Dropbox/R Files/ChickensAndPork.R.RData")

compare(m1cc.c, m2cc.c, m3cc.c, m4cc.c, m5cc.c, m6cc.c, m7cc.c, m8cc.c, m9cc.c, m10cc.c, m11cc.c.nt, 
        m12cc.c, m13cc.c, m14cc.c, m2cc.c.nt, m3cc.c.nt, m4cc.c.nt, m5cc.c.nt, m7cc.c.nt)

compare(m11cc.c.nt, m14cc.c, m2cc.c.nt, m7cc.c.nt, m13cc.c, m10cc.c, m9cc.c, m12cc.c, m2cc.c, m4cc.c.nt)

coeftab(m11cc.c.nt, m14cc.c, m2cc.c.nt, m7cc.c.nt, m13cc.c, m10cc.c, m9cc.c, m12cc.c, m2cc.c, m4cc.c.nt)

sims <- ensemble(m11cc.c.nt, m14cc.c, m2cc.c.nt, m7cc.c.nt, m13cc.c, m10cc.c, m9cc.c, m12cc.c, m2cc.c, m4cc.c.nt, 
         data = d.pred, do_sim = T, do_link = F)$sim

a <- link(m11cc.c.nt, data = d.pred)
str(a)

d.pred <- data.frame(
  Treat = c(0,1),
  gender_id = 1,
  country_id = 6,
  Age_s = 0
)

sims <- ensemble(m11cc.c.nt, m14cc.c, m2cc.c.nt, m7cc.c.nt, m13cc.c, m10cc.c, m9cc.c, m12cc.c, m2cc.c, m4cc.c.nt, 
                 data = d.pred, do_sim = T, do_link = F, n = 5e3)$sim



# Histogram Colored (blue and red)
hist(sims[,1], col=rgb(1,0,0,0.3), main = "Simulated Bird Consumption", xlab = "Number of Servings")
hist(sims[,2], col=rgb(0,0,1,0.3), add=T)
legend(7,1800, c("Control","Treatment"), lty=c(1,1), lwd=c(2.5,2.5),col=c(rgb(1,0,0,0.3),rgb(0,0,1,0.3))) 


#simulate values for p and lambda
a <- list()
a[[1]] <- as.data.frame(link(data = d.pred, m11cc.c.nt, n = 1e4)) #.41
a[[2]] <- as.data.frame(link(data = d.pred, m14cc.c, n = 1e4)) #.25
a[[3]] <- as.data.frame(link(data = d.pred, m2cc.c.nt, n = 1e4)) # .10
a[[4]] <- as.data.frame(link(data = d.pred, m7cc.c.nt, n = 1e4)) # .07
a[[5]] <- as.data.frame(link(data = d.pred, m13cc.c, n = 1e4)) #.06
a[[6]] <- as.data.frame(link(data = d.pred, m10cc.c, n = 1e4)) #.05
a[[7]] <- as.data.frame(link(data = d.pred, m9cc.c, n = 1e4)) #.03

#sample proportional to WAIC weights
b <- list()
b[[1]] <- sample_n(as.data.frame(a[[1]]), size = 4100)
b[[2]] <- sample_n(as.data.frame(a[[2]]), size = 3100)
b[[3]] <- sample_n(as.data.frame(a[[3]]), size = 1200)
b[[4]] <- sample_n(as.data.frame(a[[4]]), size = 700)
b[[5]] <- sample_n(as.data.frame(a[[5]]), size = 600)
b[[6]] <- sample_n(as.data.frame(a[[6]]), size = 500)
b[[7]] <- sample_n(as.data.frame(a[[7]]), size = 300)

#put it all in one data frame
values <- b[[1]][,c("p.1", "p.2", "lambda.1", "lambda.2")]
for (i in 2:7) { values <- rbind(values, b[[i]][,c("p.1", "p.2", "lambda.1", "lambda.2")])}

#plot the distributions
par(mfrow = c(2,1))
#ps
dens(values$p.1)
dens(values$p.2, add = T, col = 2)
title("Model-Averaged Estimates for p for Bird Consumption")
legend(0.208,31, c("Control","Treatment"), lty=c(1,1), lwd=c(2.5,2.5),col=c(1,2)) 

#lambdas
dens(values$lambda.1)
dens(values$lambda.2, add = T, col = 2)
title("Model-Averaged Estimates for Lambda for Bird Consumption")
legend(2.65,5.2, c("Control","Treatment"), lty=c(1,1), lwd=c(2.5,2.5),col=c(1,2)) 

#difference in lambda between control and treatment
par(mfrow = c(1,1))
dens(values$lambda.2 - values$lambda.1)
title("Estimated Difference in Lambda Between Control and Treatment")

diff <- values$lambda.2 - values$lambda.1
sum( diff == 0 ) / length( diff )


#effects of age on american women with and without treatment

age.seq <- seq(from = -1.5, to = 2, by = .1)

d.pred.NT <- data.frame(
  Treat = 0,
  gender_id = 1,
  country_id = 6,
  Age_s = age.seq
)

d.pred.T <- data.frame(
  Treat = 1,
  gender_id = 1,
  country_id = 6,
  Age_s = age.seq
)


sims.a.T <- ensemble(m11cc.c.nt, m14cc.c, m2cc.c.nt, m7cc.c.nt, m13cc.c, m10cc.c, m9cc.c, m12cc.c, m2cc.c, m4cc.c.nt, 
                     data = d.pred.T, do_sim = T, do_link = F, n = 5e3)$sim

consMean.T <- apply(sims.a.T, 2, mean)
consPI.T <- apply(sims.a.T, 2, PI, prob = .89)

sims.a.NT <- ensemble(m11cc.c.nt, m14cc.c, m2cc.c.nt, m7cc.c.nt, m13cc.c, m10cc.c, m9cc.c, m12cc.c, m2cc.c, m4cc.c.nt, 
                      data = d.pred.NT, do_sim = T, do_link = F, n = 5e3)$sim

consMean.NT <- apply(sims.a.NT, 2, mean)
consPI.NT <- apply(sims.a.NT, 2, PI, prob = .89)


plot( jitter(X2DaysBird) ~ jitter(Age.s) , dcc , col=col.alpha(1,0.1) ) 

lines( age.seq , consMean.T , col = 2, lwd = 2)
shade( consPI.T , age.seq, col = col.alpha(2, 0.1) )

lines( age.seq , consMean.NT , col = 1, lwd = 2)
shade( consPI.NT , age.seq)

legend(1.3,9.4, c("Control","Treatment"), lty=c(1,1), lwd=c(2.5,2.5),col=c(1, 2)) 
title("Predicted Effects of Age on Consumption")






#simulate values for p and lambda
a <- list()
a[[1]] <- as.data.frame(link(data = d.pred, m11cc.c.nt, n = 1e4)) #.41
a[[2]] <- as.data.frame(link(data = d.pred, m14cc.c, n = 1e4)) #.25
a[[3]] <- as.data.frame(link(data = d.pred, m2cc.c.nt, n = 1e4)) # .10
a[[4]] <- as.data.frame(link(data = d.pred, m7cc.c.nt, n = 1e4)) # .07
a[[5]] <- as.data.frame(link(data = d.pred, m13cc.c, n = 1e4)) #.06
a[[6]] <- as.data.frame(link(data = d.pred, m10cc.c, n = 1e4)) #.05
a[[7]] <- as.data.frame(link(data = d.pred, m9cc.c, n = 1e4)) #.03

#sample proportional to WAIC weights
b <- list()
b[[1]] <- sample_n(as.data.frame(a[[1]]), size = 4100)
b[[2]] <- sample_n(as.data.frame(a[[2]]), size = 3100)
b[[3]] <- sample_n(as.data.frame(a[[3]]), size = 1200)
b[[4]] <- sample_n(as.data.frame(a[[4]]), size = 700)
b[[5]] <- sample_n(as.data.frame(a[[5]]), size = 600)
b[[6]] <- sample_n(as.data.frame(a[[6]]), size = 500)
b[[7]] <- sample_n(as.data.frame(a[[7]]), size = 300)








## Estimate lambda and p for a variety of ages ##

age.seq <- seq(from = -1.5, to = 2, by = .1)

d.pred.NT <- data.frame(
  Treat = 0,
  gender_id = 1,
  country_id = 6,
  Age_s = age.seq
)

d.pred.T <- data.frame(
  Treat = 1,
  gender_id = 1,
  country_id = 6,
  Age_s = age.seq
)

#simulate values for p and lambda for NT

aNT <- list()
aNT[[1]] <- as.data.frame(link(data = d.pred.NT, m11cc.c.nt, n = 1e4)) #.41
aNT[[2]] <- as.data.frame(link(data = d.pred.NT, m14cc.c, n = 1e4)) #.25
aNT[[3]] <- as.data.frame(link(data = d.pred.NT, m2cc.c.nt, n = 1e4)) # .10
aNT[[4]] <- as.data.frame(link(data = d.pred.NT, m7cc.c.nt, n = 1e4)) # .07
aNT[[5]] <- as.data.frame(link(data = d.pred.NT, m13cc.c, n = 1e4)) #.06
aNT[[6]] <- as.data.frame(link(data = d.pred.NT, m10cc.c, n = 1e4)) #.05
aNT[[7]] <- as.data.frame(link(data = d.pred.NT, m9cc.c, n = 1e4)) #.03


#simulate values for p and lambda for T
aT <- list()
aT[[1]] <- as.data.frame(link(data = d.pred.T, m11cc.c.nt, n = 1e4)) #.41
aT[[2]] <- as.data.frame(link(data = d.pred.T, m14cc.c, n = 1e4)) #.25
aT[[3]] <- as.data.frame(link(data = d.pred.T, m2cc.c.nt, n = 1e4)) # .10
aT[[4]] <- as.data.frame(link(data = d.pred.T, m7cc.c.nt, n = 1e4)) # .07
aT[[5]] <- as.data.frame(link(data = d.pred.T, m13cc.c, n = 1e4)) #.06
aT[[6]] <- as.data.frame(link(data = d.pred.T, m10cc.c, n = 1e4)) #.05
aT[[7]] <- as.data.frame(link(data = d.pred.T, m9cc.c, n = 1e4)) #.03

#sample proportional to WAIC weights for NT
bNT <- list()
bNT[[1]] <- sample_n(as.data.frame(aNT[[1]]), size = 4100)
bNT[[2]] <- sample_n(as.data.frame(aNT[[2]]), size = 2500)
bNT[[3]] <- sample_n(as.data.frame(aNT[[3]]), size = 1000)
bNT[[4]] <- sample_n(as.data.frame(aNT[[4]]), size = 700)
bNT[[5]] <- sample_n(as.data.frame(aNT[[5]]), size = 600)
bNT[[6]] <- sample_n(as.data.frame(aNT[[4]]), size = 500)
bNT[[7]] <- sample_n(as.data.frame(aNT[[5]]), size = 300)

#sample proportional to WAIC weights for T
bT <- list()
bT[[1]] <- sample_n(as.data.frame(aT[[1]]), size = 4100)
bT[[2]] <- sample_n(as.data.frame(aT[[2]]), size = 2500)
bT[[3]] <- sample_n(as.data.frame(aT[[3]]), size = 1000)
bT[[4]] <- sample_n(as.data.frame(aT[[4]]), size = 700)
bT[[5]] <- sample_n(as.data.frame(aT[[5]]), size = 600)
bT[[6]] <- sample_n(as.data.frame(aT[[4]]), size = 500)
bT[[7]] <- sample_n(as.data.frame(aT[[5]]), size = 300)

#put it all in one data frame for NT
valuesNT <- bNT[[1]][,c(paste0("p.", seq(1:36)), paste0("lambda.", seq(1:36)))]
for (i in 2:7) { valuesNT <- rbind(valuesNT, bNT[[i]][,c(paste0("p.", seq(1:36)), paste0("lambda.", seq(1:36)))])}

#put it all in one data frame for T
valuesT <- bT[[1]][,c(paste0("p.", seq(1:36)), paste0("lambda.", seq(1:36)))]
for (i in 2:7) { valuesT <- rbind(valuesT, bT[[i]][,c(paste0("p.", seq(1:36)), paste0("lambda.", seq(1:36)))])}

valuesNT.mean <- apply(valuesNT, 2, mean)
valuesT.mean <- apply(valuesT, 2, mean)

valuesNT.HPDI <- apply(valuesNT, 2, HPDI)
valuesT.HPDI <- apply(valuesT, 2, HPDI)

#convert back to original age values
trueAgeSeq <- age.seq * sd(d$c.Age) + mean(d$c.Age)

#plotting change in p with age
plot(trueAgeSeq, valuesNT.mean[1:36], type = "l", ylim = c(0.13, 0.2), lwd = 2, xlab = "Age", ylab = "Predicted p")
lines(trueAgeSeq, valuesT.mean[1:36], col = 2)
shade(valuesNT.HPDI[,1:36], trueAgeSeq)
shade(valuesT.HPDI[,1:36], trueAgeSeq, col = col.alpha(2, alpha = .1))
legend(22.8,.2, c("Control","Treatment"), lty=c(1,1), lwd=c(2.5,2.5),col=c(1,2), bty = "n") 
title("Estimated Effects of Exposure to Ad on Probability of Bird Abstention")

#plotting change in lambda with age
plot(trueAgeSeq, valuesNT.mean[37:72], type = "l", ylim = c(2, 3), lwd = 2, xlab = "Age", ylab = "Predicted ??")
lines(trueAgeSeq, valuesT.mean[37:72], col = 2, lwd = 2)
shade(valuesNT.HPDI[,37:72], trueAgeSeq)
shade(valuesT.HPDI[,37:72], trueAgeSeq, col = col.alpha(2, alpha = .1))
legend(23.1,3, c("Control","Treatment"), lty=c(1,1), lwd=c(2.5,2.5),col=c(1,2)) 
title("Estimated Effects of Exposure to Ad on Rate of Bird Consumption")

#plot the distributions of change in lambda as individuals age
par(mfrow = c(1,1))
for (i in 37:72){
  print(i)
  png(paste0("lambdaChicken", i, ".png"), width = 900)
  diffs <- valuesT[,i] - valuesNT[,i]
  dens(diffs, xlim = c(-.75, .75), ylim = c(0, 3.5))
  title(paste0("Posterior Distribution for Estimated Differences between Treatment and Control on Bird Consumption (??) for Age: ~", round(trueAgeSeq[i-36])))
  abline(v = 0, lty = 2)
  text(x = mean(diffs), y = .1, labels = paste0("% Neg: ", round(sum( diffs < 0 ) / length( diffs ) * 100, digits = 2), "%"))
  dev.off()
}

diffs <- valuesT[,52] - valuesNT[,52]
dens(diffs, xlim = c(-.5, .5), ylim = c(0, 3.5))
title("Posterior Distribution for Estimated Differences between Treatment and Control on Bird Consumption Rate (??)")
abline(v = 0, lty = 2)

