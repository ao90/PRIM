# PRIM package
Implementation of the Patient Rule Induction Method (PRIM) as suggested by Friedman and Fisher (1999) in R


 
 
## Installation
 
You can install the latest version from
[Github](https://github.com/ao90/PRIM)
 
```s
# install.packages("devtools")
library(devtools)
install_github("ao90/PRIM") 
```
 
## Usage
 
```s
library(PRIM)
 
# generating random data:
set.seed(123)
n <- 500
x1 <- runif(n = n, min = -1)
x2 <- runif(n = n, min = -1)
x3 <- runif(n = n, min = -1)
cat <- as.factor(sample(c("a","b","c", "d"), size = n, replace = TRUE))
wsk <- (1-sqrt(x1^2+x2^2)/sqrt(2))
y <- as.logical(rbinom(n = n, prob = wsk, size = 1))
dat <- cbind.data.frame(y, x1, x2, x3, cat)
plot(dat$x1, dat$x2, col=dat$y+1, pch=16)
remove(x1, x2, x3, y, wsk, cat, n)

# run (multiple) peeling and show trajectory:
prim <- PRIM_peel_bs(formula=y ~ ., data=dat, beta_min = .01)
plot(prim) # multiple trajectory
head(prim$box) # box definitions

# apply the PRIM function to find the best boxes with a support of at least 0.1:
p <- PRIM(y~., data=dat, beta_min = 0.1, max_boxes = 3)
p
 
```
 
## License
 
This package is free and open source software, licensed under GPL (>= 3).
