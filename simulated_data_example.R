
dat1 = read.csv("simulated_data/4x4.B.change.1.csv")
dat2 = read.csv("simulated_data/4x4.B.change.2.csv")
dat3 = read.csv("simulated_data/4x4.B.change.3.csv")
dat4 = read.csv("simulated_data/4x4.B.change.4.csv")

y = array(NA, dim=c(4, nrow(dat1), 4))
y[1,,] = as.matrix(dat1[,c("g1","g2","P1","P2")])
y[2,,] = as.matrix(dat2[,c("g1","g2","P1","P2")])
y[3,,] = as.matrix(dat3[,c("g1","g2","P1","P2")])
y[4,,] = as.matrix(dat4[,c("g1","g2","P1","P2")])


# fit the tvvarss model
model = tvvarss(y=y, include_trend=TRUE, de_mean=TRUE, mcmc_chain=3, mcmc_iter=700, mcmc_warmup = 300)

# fit tvvarss model without de-meaning
model = tvvarss(y=y, de_mean=FALSE)
# spp names are "Grazer 1","Grazer 2","Predator 1","Predator 2"
B = extract(model, pars=c("B"))
B_mean = apply(B[[1]], c(2,3,4), mean)
par(mfrow = c(4,4), mai=c(0.1,0.1,0.1,0.1))
for(i in 1:4) {
  for(j in 1:4) {
    if(i==j) plot(B_mean[,j,i], ylim=c(0,1))
    if(i!=j) plot(B_mean[,j,i], ylim=c(-1,1))
  }
}

