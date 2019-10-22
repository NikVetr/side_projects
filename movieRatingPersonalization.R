smallData <- F
approxProbit <- F
logitNormal <- T

sigmoid <- function(x) {1/(1+exp(-x))}; logit <- function(p) {log(p/(1-p))}

if(smallData){
  movies <- read.csv("ml-latest-small/movies.csv")
  ratings <- read.csv(file = "ml-latest-small/ratings.csv")[,-4]
} else {
  ratings <- readLines("ml-10M100K/ratings.dat")
  ratings <- as.data.frame(t(sapply(1:length(ratings), function(x) 
    as.numeric(strsplit(ratings[x], "::")[[1]])))[,-4], col.names = c("userId","movieId","rating"))
  colnames(ratings) <- c("userId","movieId","rating")
  movies <- readLines("ml-10M100K/movies.dat")
  movies <- as.data.frame(t(sapply(1:length(movies), function(x) strsplit(movies[x], split = "::")[[1]])))
  colnames(movies) <- c("movieId","title","genres")
}

imdb <- read.csv(file = "IDMB_top250.csv")

#convert ratings to approx. probit
if(approxProbit){
  ratings$rating <- qnorm(replace(ratings$rating, list = which(ratings$rating == 5), 4.875)/5, 0, 1)
}

#convert ratings ~normal~ by assuming logit-normal
if(logitNormal){
  ratings$rating <-logit(replace(ratings$rating, list = which(ratings$rating == 5), 4.875)/5)
}
#find subset of top 100 most commonly rated movies
ntop <- 300
top <- as.numeric(attr(sort(table(ratings$movieId), decreasing = T)[1:ntop], which = "dimnames")[[1]])
topMovs <- movies$title[match(top, movies$movieId)]
ratings <- ratings[sapply(1:length(ratings$movieId), function(x) any(ratings$movieId[x] == top)),]
meanRatings <- sapply(1:ntop, function(x) mean(ratings[ratings$movieId == top[x],3]))
users <- unique(ratings$userId)
ratmoves <- sapply(1:length(users), function(x) ratings[ratings$userId == users[x],2])
rats <- sapply(1:length(users), function(x) ratings[ratings$userId == users[x],3])
ratmoves <- sapply(1:length(ratmoves), function(x) match(ratmoves[[x]], top))
sparseVecNA <- function(length, indices, data){
  foo <- rep(NA, length) 
  foo[indices] <- data
  foo
}
topRatingsMatrix <- t(sapply(1:length(users), function(x) sparseVecNA(length = ntop, indices = ratmoves[[x]], data = rats[[x]])))
topCovs <- cov(topRatingsMatrix, us = "pairw"); topCovs[is.na(topCovs)] <- 0
library(Matrix)
topCovsPSD <- as.matrix(nearPD(topCovs)$mat)

nRats <- 10
if(logitNormal){
  rawRatings <- sample(1:10, nRats, replace = T)/2
  persRatings <- sparseVecNA(length = ntop, 
                             data = logit(replace(rawRatings, list = which(rawRatings == 5), 4.75)/5), indices = sample(1:ntop, nRats, replace = F))
} else if(approxProbit) {
  
} else {
  persRatings <- sparseVecNA(length = ntop, data = sample(1:10, nRats, replace = T)/2, indices = sample(1:ntop, nRats, replace = F))
}
# specific movie mode
persRatings <- sparseVecNA(length = ntop, data = rep(logit(4.9/5), sum(startsWith(as.character(topMovs), "Harry"))),
                           indices = which(startsWith(as.character(topMovs), "Harry")))
persMovs <- which(!is.na(persRatings))
persRats <- persRatings[persMovs]
condMVN <- function(means, cov, obs, inds){
  A <- cov[-inds, -inds]
  C <- cov[inds, inds]
  B <- cov[-inds, inds]
  condCov <- A - B%*%solve(C)%*%t(B)
  condMean <- as.vector(means[-inds] + B%*%solve(C)%*%(obs - means[inds]))
  return(list(means = condMean, cov = condCov))
}
preds <- condMVN(means = meanRatings, cov = topCovsPSD, obs = persRats, inds = persMovs)
predMovs <- (1:length(top))[-persMovs]
if(logitNormal){
  head(cbind(as.character(topMovs[persMovs[order(persRats, decreasing = T)]]), sort(sigmoid(persRats)*5, decreasing = T)), 20)
} 
if(logitNormal){
    head(cbind(as.character(topMovs[predMovs[order(preds$means, decreasing = T)]]), sort(sigmoid(preds$means)*5, decreasing = T)), 20)
}
if(!logitNormal & !approxProbit) {
  head(cbind(as.character(topMovs[persMovs[order(persRats, decreasing = T)]]), sort(persRats, decreasing = T)), 20)
  head(cbind(as.character(topMovs[predMovs[order(preds$means, decreasing = T)]]), sort(preds$means, decreasing = T)), 20)
}

#find maximally surprising movies
if(logitNormal){

  head(cbind(as.character(topMovs[predMovs[order(sigmoid(preds$means) - sigmoid(meanRatings[-persMovs]), decreasing = T)]]), 
             sort((sigmoid(preds$means) - sigmoid(meanRatings[-persMovs]))*5, decreasing = T)), 20)
}