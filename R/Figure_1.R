library("tvvarss")
library("broom")
library("rstan")

load("output_linear_chain_staticB.Rdata")

n_simulations = 100
sim_config = expand.grid("B_diag" = c("low", "med", "high"),
  "sim" = 1:n_simulations, "site" = c(1,2,4))

output = matrix(NA, nrow(sim_config), 16)
# just look at B[1,1]
for(i in 101:length(saved_output)) {

trueB = c(saved_output[[i]]$sim_output$B_mat[,,1])
estB = saved_output[[i]]$estimate$estimate
output[i,] = trueB - estB

}

output = as.data.frame(output)
names(output) = saved_output[[101]]$estimate$term
sim_config = cbind(sim_config, output)

# reorder columns
sim_config = sim_config[, c(1:3, 4, 8, 12, 16, 5, 9, 13, 17, 6, 10, 14, 18, 7, 11, 15, 19)]
# make boxplot of bias as function of site
par(mfrow = c(4,4), mgp = c(2,1,0), mai = c(0.5,0.5,0.1,0.1))
for(i in 1:16) {
  boxplot(sim_config[,c(3+i)] ~ sim_config$site, outline = F, xlab = "Site", main=names(sim_config)[c(3+i)])
  abline(0, 0, col="red")
}

