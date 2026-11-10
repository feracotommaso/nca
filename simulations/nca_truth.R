
library(NCA)

?NCA::nca_random()

r <- c(NA_real_)
for (i in 1:10) {
  print(i)
  d<-nca_random(
    n = 1000,
    intercepts = 0.3,
    slopes = 0.5,
    distribution.x = "normal",
    distribution.y = "normal",
    mean.x = .5,
    mean.y = .5,
    sd.x = .2,
    sd.y = .2
  )
  r[i] <- cor(d$X,d$Y)
}

hist(r)
mean(r)
median(r)

plot(d$X,d$Y)

nca_analysis(d,1,2)
cor(d$X,d$Y)

hist(d$Y, breaks = 30)

# N = 1000
> mean(r)
[1] 0.0787068
> median(r)
[1] 0.07706909