# install.packages(c("parallel", "mgcv", "dlnm", "mvmeta", "doParallel", "data.table"))
library(parallel)
library(mgcv)
library(dlnm)
library(mvmeta)
library(splines)
library(doParallel)
library(data.table)

## To run the analysis for Zone2, Zone3, or Zone1, simply replace all instances of "Zone1"
## Please change the following paths to match your own file locations.
data <- fread("/Users/huasong/Desktop/data/3.Summary data/Zone1.csv")

name_list <- unique(data$Point_num) 
variables <- c("PO4", "NO3", "SST", "SSW") 
lag <- c(0, 15)
output_dir <- "/Users/huasong/Desktop/"

## The number of cores should be selected based on the computer’s configuration.
num_cores <- min(5, detectCores() - 2)
cl <- makeCluster(num_cores)
clusterEvalQ(cl, {
  library(mgcv)
  library(dlnm)
  library(mvmeta)
  library(splines)
})


calculate_lag_effect <- function(data_subset, variable, lag) {

  ylag_25 <- matrix(NA, nrow = length(unique(data_subset$Point_num)), ncol = 3)
  ylag_75 <- matrix(NA, nrow = length(unique(data_subset$Point_num)), ncol = 3)
  Slag_25 <- list()
  Slag_75 <- list()

  
  site_names <- unique(data_subset$Point_num)
  

  for (i in seq_along(site_names)) {
    site_data <- subset(data_subset, Point_num == site_names[i])
    

    argvar <- list(fun = "ns", df = 3)
    arglag <- list(fun = "ns", df = 3)
    cb <- crossbasis(site_data[[variable]], lag = lag, argvar = argvar, arglag = arglag)
    

    other_vars <- c("PO4", "NO3", "SST", "SSW")
    other_vars <- setdiff(other_vars, variable)
    smooth_terms <- paste0("ns(site_data$", other_vars, ", 4)", collapse = " + ")
    formula_str <- paste("Chla ~ cb +", smooth_terms, "+ ns(site_data$Time, 4*30)")
    formula <- as.formula(formula_str)
    

    model <- gam(formula, family = Gamma("log"), data = site_data)
    
    suppressWarnings({
      crlag_25 <- crossreduce(cb, model, type = "var", value = quantile(site_data[[variable]], 0.25), cen = quantile(site_data[[variable]], 0))
      crlag_75 <- crossreduce(cb, model, type = "var", value = quantile(site_data[[variable]], 0.75), cen = quantile(site_data[[variable]], 0))
    })
    

    ylag_25[i, ] <- coef(crlag_25)
    Slag_25[[i]] <- vcov(crlag_25)
    ylag_75[i, ] <- coef(crlag_75)
    Slag_75[[i]] <- vcov(crlag_75)
  }
  

  return(list(ylag_25 = ylag_25, ylag_75 = ylag_75, Slag_25 = Slag_25, Slag_75 = Slag_75))
}


batch_size <- 1000  
num_batches <- ceiling(length(name_list) / batch_size)

all_results <- list()  

for (batch in 1:num_batches) {
  cat("Processing batch", batch, "of", num_batches, "\n")
  

  batch_names <- name_list[((batch - 1) * batch_size + 1):min(batch * batch_size, length(name_list))]
  batch_data <- subset(data, Point_num %in% batch_names)
  

  clusterExport(cl, varlist = c("calculate_lag_effect", "batch_data", "lag", "variables"))
  batch_results <- parLapply(cl, variables, function(var) {
    calculate_lag_effect(batch_data, var, lag)
  })
  

  all_results <- append(all_results, batch_results)
}
save(all_results, file = "/Users/huasong/Desktop/all_results_Zone1.RData")


ylag_25_po4 <- do.call(rbind, lapply(all_results[1], function(res) res$ylag_25))
ylag_75_po4 <- do.call(rbind, lapply(all_results[1], function(res) res$ylag_75))


ylag_25_no3 <- do.call(rbind, lapply(all_results[2], function(res) res$ylag_25))
ylag_75_no3 <- do.call(rbind, lapply(all_results[2], function(res) res$ylag_75))


ylag_25_sst <- do.call(rbind, lapply(all_results[3], function(res) res$ylag_25))
ylag_75_sst <- do.call(rbind, lapply(all_results[3], function(res) res$ylag_75))


ylag_25_ssw <- do.call(rbind, lapply(all_results[4], function(res) res$ylag_25))
ylag_75_ssw <- do.call(rbind, lapply(all_results[4], function(res) res$ylag_75))


Slag_25_po4 <- do.call(rbind,lapply(all_results[1], function(res) res$Slag_25))
Slag_75_po4 <- do.call(rbind,lapply(all_results[1], function(res) res$Slag_75))


Slag_25_no3 <- do.call(rbind,lapply(all_results[2], function(res) res$Slag_25))
Slag_75_no3 <- do.call(rbind,lapply(all_results[2], function(res) res$Slag_75))


Slag_25_sst <- do.call(rbind,lapply(all_results[3], function(res) res$Slag_25))
Slag_75_sst <- do.call(rbind,lapply(all_results[3], function(res) res$Slag_75))


Slag_25_ssw <- do.call(rbind,lapply(all_results[4], function(res) res$Slag_25))
Slag_75_ssw <- do.call(rbind,lapply(all_results[4], function(res) res$Slag_75))


mvlag_25_po4 <- mvmeta(ylag_25_po4 ~ 1, Slag_25_po4, method = "reml")
mvlag_75_po4 <- mvmeta(ylag_75_po4 ~ 1, Slag_75_po4, method = "reml")


mvlag_25_no3 <- mvmeta(ylag_25_no3 ~ 1, Slag_25_no3, method = "reml")
mvlag_75_no3 <- mvmeta(ylag_75_no3 ~ 1, Slag_75_no3, method = "reml")


mvlag_25_sst <- mvmeta(ylag_25_sst ~ 1, Slag_25_sst, method = "reml")
mvlag_75_sst <- mvmeta(ylag_75_sst ~ 1, Slag_75_sst, method = "reml")


mvlag_25_ssw <- mvmeta(ylag_25_ssw ~ 1, Slag_25_ssw, method = "reml")
mvlag_75_ssw <- mvmeta(ylag_75_ssw ~ 1, Slag_75_ssw, method = "reml")


mvlag_25_all <- list(mvlag_25_po4, mvlag_25_no3, mvlag_25_sst, mvlag_25_ssw)
mvlag_75_all <- list(mvlag_75_po4, mvlag_75_no3, mvlag_75_sst, mvlag_75_ssw)




xlag <- 0:15 
blag <-lapply(1:length(variables), function(i) {
  do.call("onebasis", c(list(x = xlag), attr(crossbasis(data[[variables[i]]], lag = lag, argvar = list(fun = "ns", df = 3), arglag = list(fun = "ns", df = 3)), "arglag")))
})

lag_25_all <- lapply(1:length(variables), function(i) {
  crosspred(blag[[i]], coef = coef(mvlag_25_all[[i]]), vcov = vcov(mvlag_25_all[[i]]), model.link = "log", at = 0:15, cen = quantile(data[[variables[i]]], 0))
})
lag_75_all <- lapply(1:length(variables), function(i) {
  crosspred(blag[[i]], coef = coef(mvlag_75_all[[i]]), vcov = vcov(mvlag_75_all[[i]]), model.link = "log", at = 0:15, cen = quantile(data[[variables[i]]], 0))
})


for (i in 1:length(variables)) {
  write.csv(data.frame(variable = variables[i], lag = lag_25_all[[i]]$predvar, RRfit = lag_25_all[[i]]$allRRfit, RRlow = lag_25_all[[i]]$allRRlow, RRhigh = lag_25_all[[i]]$allRRhigh),
            paste0(output_dir, variables[i], "_lag25_all_Zone1.csv"), row.names = FALSE)
  write.csv(data.frame(variable = variables[i], lag = lag_75_all[[i]]$predvar, RRfit = lag_75_all[[i]]$allRRfit, RRlow = lag_75_all[[i]]$allRRlow, RRhigh = lag_75_all[[i]]$allRRhigh),
            paste0(output_dir, variables[i], "_lag75_all_Zone1.csv"), row.names = FALSE)

}


stopCluster(cl)
cat("All computations are complete! The data has been saved.\n")


for (i in 1:length(variables)) {
  assign(paste0(variables[i], "_lag_25_all"), read.csv(paste0(output_dir, variables[i], "_lag25_all_Zone1.csv")))
  assign(paste0(variables[i], "_lag_75_all"), read.csv(paste0(output_dir, variables[i], "_lag75_all_Zone1.csv")))
 
}

## Please change the following paths to match your own file locations.
png("/Users/huasong/Desktop/lag_Zone1_all.png", width=18, height=17, units="cm", res=300)
par(mfrow=c(2,2))


plot(PO4_lag_25_all$lag, PO4_lag_25_all$RRfit, type = "l", col = "red", lwd = 2,
     ylab = "CRR", xlab = "Lag", cex.lab = 1.5,
     ylim = range(c(PO4_lag_25_all$RRlow, PO4_lag_25_all$RRhigh, PO4_lag_75_all$RRlow, PO4_lag_75_all$RRhigh)))
polygon(c(PO4_lag_25_all$lag, rev(PO4_lag_25_all$lag)),
        c(PO4_lag_25_all$RRlow, rev(PO4_lag_25_all$RRhigh)),
        col = adjustcolor("red", alpha.f = 0.2), border = NA)
lines(PO4_lag_75_all$lag, PO4_lag_75_all$RRfit, col = "blue", lwd = 2)
polygon(c(PO4_lag_75_all$lag, rev(PO4_lag_75_all$lag)),
        c(PO4_lag_75_all$RRlow, rev(PO4_lag_75_all$RRhigh)),
        col = adjustcolor("blue", alpha.f = 0.2), border = NA)
abline(h = 1, lty = 2, col = "black")
legend("topright", legend = c(expression(PO[4]*" 25th"), expression(PO[4]*" 75th")),
       col = c("red", "blue"), lty = 1, lwd = 2, bty="n", cex = 1)


plot(NO3_lag_25_all$lag, NO3_lag_25_all$RRfit, type = "l", col = "red", lwd = 2,
     ylab = "CRR", xlab = "Lag", cex.lab = 1.5,
     ylim = range(c(NO3_lag_25_all$RRlow, NO3_lag_25_all$RRhigh, NO3_lag_75_all$RRlow, NO3_lag_75_all$RRhigh)))
polygon(c(NO3_lag_25_all$lag, rev(NO3_lag_25_all$lag)),
        c(NO3_lag_25_all$RRlow, rev(NO3_lag_25_all$RRhigh)),
        col = adjustcolor("red", alpha.f = 0.2), border = NA)
lines(NO3_lag_75_all$lag, NO3_lag_75_all$RRfit, col = "blue", lwd = 2)
polygon(c(NO3_lag_75_all$lag, rev(NO3_lag_75_all$lag)),
        c(NO3_lag_75_all$RRlow, rev(NO3_lag_75_all$RRhigh)),
        col = adjustcolor("blue", alpha.f = 0.2), border = NA)
abline(h = 1, lty = 2, col = "black")
legend("topright", legend = c(expression(NO[3]*" 25th"), expression(NO[3]*" 75th")),
       col = c("red", "blue"), lty = 1, lwd = 2, bty="n", cex = 1)


plot(SST_lag_25_all$lag, SST_lag_25_all$RRfit, type = "l", col = "red", lwd = 2,
     ylab = "CRR", xlab = "Lag", cex.lab = 1.5,
     ylim = range(c(SST_lag_25_all$RRlow, SST_lag_25_all$RRhigh+0.03, SST_lag_75_all$RRlow, SST_lag_75_all$RRhigh+0.03)))
polygon(c(SST_lag_25_all$lag, rev(SST_lag_25_all$lag)),
        c(SST_lag_25_all$RRlow, rev(SST_lag_25_all$RRhigh)),
        col = adjustcolor("red", alpha.f = 0.2), border = NA)
lines(SST_lag_75_all$lag, SST_lag_75_all$RRfit, col = "blue", lwd = 2)
polygon(c(SST_lag_75_all$lag, rev(SST_lag_75_all$lag)),
        c(SST_lag_75_all$RRlow, rev(SST_lag_75_all$RRhigh)),
        col = adjustcolor("blue", alpha.f = 0.2), border = NA)
abline(h = 1, lty = 2, col = "black")
legend("topright", legend = c(expression(SST*" 25th"), expression(SST*" 75th")),
       col = c("red", "blue"), lty = 1, lwd = 2, bty="n", cex = 1)


plot(SSW_lag_25_all$lag, SSW_lag_25_all$RRfit, type = "l", col = "red", lwd = 2,
     ylab = "CRR", xlab = "Lag", cex.lab = 1.5,
     ylim = range(c(SSW_lag_25_all$RRlow, SSW_lag_25_all$RRhigh, SSW_lag_75_all$RRlow, SSW_lag_75_all$RRhigh)))
polygon(c(SSW_lag_25_all$lag, rev(SSW_lag_25_all$lag)),
        c(SSW_lag_25_all$RRlow, rev(SSW_lag_25_all$RRhigh)),
        col = adjustcolor("red", alpha.f = 0.2), border = NA)
lines(SSW_lag_75_all$lag, SSW_lag_75_all$RRfit, col = "blue", lwd = 2)
polygon(c(SSW_lag_75_all$lag, rev(SSW_lag_75_all$lag)),
        c(SSW_lag_75_all$RRlow, rev(SSW_lag_75_all$RRhigh)),
        col = adjustcolor("blue", alpha.f = 0.2), border = NA)
abline(h = 1, lty = 2, col = "black")
legend("topright", legend = c(expression(SSW*" 25th"), expression(SSW*" 75th")),
       col = c("red", "blue"), lty = 1, lwd = 2, bty="n", cex = 1)
mtext("Zone1",side=3,line=-2,outer=TRUE,cex=2,font=2,adj=0,col="black")
dev.off()