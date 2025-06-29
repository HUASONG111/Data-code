# install.packages(c("dlnm", "mgcv", "dplyr", "ggplot2", "data.table", "cowplot"))
library(dlnm) 
library(mgcv)  
library(splines) 
library(dplyr)   
library(ggplot2) 
library(parallel) 
library(data.table)
library(cowplot)

## Set number of cores based on computer configuration 
## Detect total number of logical CPU cores (num_cores < total_cores)
# total_cores <- parallel::detectCores(logical = TRUE)
num_cores <- 5 
cl <- makeCluster(num_cores) 

clusterEvalQ(cl, {
  library(dlnm)
  library(mgcv)
  library(splines)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(parallel)
  library(cowplot)
})

## Please change the following paths to match your own file locations.
## Set 'data_dir' to the folder where your input data files are stored,
## and 'output_dir' to the folder where you want to save the results.
data_dir <- "/Users/huasong/Desktop/data/1.Original data of each single point/"
output_dir <- "/Users/huasong/Desktop/"

clusterExport(cl, c('data_dir', 'output_dir', 'fread'))

process_file <- function(i) {
  results_list <- list()
  mat_po4_list <- list()
  mat_no3_list <- list()
  mat_sst_list <- list()
  mat_ssw_list <- list()
  cum_po4_list <- list()
  cum_no3_list <- list()
  cum_sst_list <- list()
  cum_ssw_list <- list()
  final_results <- list()
  lag_results <- list()
  
  cat(sprintf("Processing file %d of 100\n", i))
  filename <- sprintf("%schl_%02d.csv", data_dir, i)
  data <- fread(filename) 
  if (any(is.na(data))) {
    warning(sprintf("skip file: %s", filename))
    return(NULL)
  }
  argvar <- list(fun = "ns", df = 3) 
  arglag <- list(fun = "ns", df = 3) 
  
  cb.po4 <- crossbasis(data$PO4, lag = 15, argvar = argvar, arglag = arglag) 
  mode_dlnm_po4 <- gam(Chla ~ cb.po4 + ns(NO3, 4) + ns(SST, 4) + ns(SSW, 4) + ns(Day,120), family = Gamma(link="log"), data = data)
  pred.po4 <- crosspred(cb.po4, mode_dlnm_po4, by = 0.01, cen = min(data$PO4), cumul = TRUE)
  po4delta <- which.min(abs(pred.po4$predvar - (min(data$PO4) + 10)))
  mat_po4 <- as.data.frame(round(cbind(pred.po4$matRRfit[po4delta,], pred.po4$matRRlow[po4delta,], pred.po4$matRRhigh[po4delta,]), digits = 3))
  mat_po4 <- mat_po4 %>% mutate(erf = (V1 - 1) * 100, erl = (V2 - 1) * 100, erh = (V3 - 1) * 100, lag = c(0:15))
  cum_po4 <- as.data.frame(round(cbind(pred.po4$cumRRfit[po4delta,], pred.po4$cumRRlow[po4delta,], pred.po4$cumRRhigh[po4delta,]), digits = 3))
  cum_po4 <- cum_po4 %>% mutate(erf = (V1 - 1) * 100, erl = (V2 - 1) * 100, erh = (V3 - 1) * 100, lag = c(0:15))
  
  
  cb.no3 <- crossbasis(data$NO3, lag = 15, argvar = argvar, arglag = arglag)
  mode_dlnm_no3 <- gam(Chla ~ cb.no3 + ns(PO4, 4) + ns(SST, 4) + ns(SSW, 4) +ns(Day,120), family = Gamma(link="log"), data = data)
  pred.no3 <- crosspred(cb.no3, mode_dlnm_no3, by = 0.01, cen = min(data$NO3), cumul = TRUE)
  no3delta <- which.min(abs(pred.no3$predvar - (min(data$NO3) + 10)))
  mat_no3 <- as.data.frame(round(cbind(pred.no3$matRRfit[no3delta,], pred.no3$matRRlow[no3delta,], pred.no3$matRRhigh[no3delta,]), digits = 3))
  mat_no3 <- mat_no3 %>% mutate(erf = (V1 - 1) * 100, erl = (V2 - 1) * 100, erh = (V3 - 1) * 100, lag = c(0:15))
  cum_no3 <- as.data.frame(round(cbind(pred.no3$cumRRfit[no3delta,], pred.no3$cumRRlow[no3delta,], pred.no3$cumRRhigh[no3delta,]), digits = 3))
  cum_no3 <- cum_no3 %>% mutate(erf = (V1 - 1) * 100, erl = (V2 - 1) * 100, erh = (V3 - 1) * 100, lag = c(0:15))
  
  cb.sst <- crossbasis(data$SST, lag = 15, argvar = argvar, arglag = arglag)
  mode_dlnm_sst <- gam(Chla ~ cb.sst + ns(SSW, 4) + ns(PO4, 4) + ns(NO3, 4) +ns(Day,120), family = Gamma(link="log"), data = data)
  pred.sst <- crosspred(cb.sst, mode_dlnm_sst, by = 0.01, cen = min(data$SST), cumul = TRUE)
  sstdelta <- which.min(abs(pred.sst$predvar - (min(data$SST) + 1)))
  mat_sst <- as.data.frame(round(cbind(pred.sst$matRRfit[sstdelta,], pred.sst$matRRlow[sstdelta,], pred.sst$matRRhigh[sstdelta,]), digits = 3))
  mat_sst <- mat_sst %>% mutate(erf = (V1 - 1) * 100, erl = (V2 - 1) * 100, erh = (V3 - 1) * 100, lag = c(0:15))
  cum_sst <- as.data.frame(round(cbind(pred.sst$cumRRfit[sstdelta,], pred.sst$cumRRlow[sstdelta,], pred.sst$cumRRhigh[sstdelta,]), digits = 3))
  cum_sst <- cum_sst %>% mutate(erf = (V1 - 1) * 100, erl = (V2 - 1) * 100, erh = (V3 - 1) * 100, lag = c(0:15))
  
  cb.ssw <- crossbasis(data$SSW, lag = 15, argvar = argvar, arglag = arglag)
  mode_dlnm_ssw <- gam(Chla ~ cb.ssw + ns(SST, 4) + ns(PO4, 4) + ns(NO3, 4) +ns(Day,120), family = Gamma(link="log"), data = data)
  pred.ssw <- crosspred(cb.ssw, mode_dlnm_ssw, by = 0.01, cen = min(data$SSW), cumul = TRUE)
  sswdelta <- which.min(abs(pred.ssw$predvar - (min(data$SSW) + 1)))
  mat_ssw <- as.data.frame(round(cbind(pred.ssw$matRRfit[sswdelta,], pred.ssw$matRRlow[sswdelta,], pred.ssw$matRRhigh[sswdelta,]), digits = 3))
  mat_ssw <- mat_ssw %>% mutate(erf = (V1 - 1) * 100, erl = (V2 - 1) * 100, erh = (V3 - 1) * 100, lag = c(0:15))
  cum_ssw <- as.data.frame(round(cbind(pred.ssw$cumRRfit[sswdelta,], pred.ssw$cumRRlow[sswdelta,], pred.ssw$cumRRhigh[sswdelta,]), digits = 3))
  cum_ssw <- cum_ssw %>% mutate(erf = (V1 - 1) * 100, erl = (V2 - 1) * 100, erh = (V3 - 1) * 100, lag = c(0:15))
  
  mat_po4_list[[i]] <- cbind(Filename = i, mat_po4)
  mat_no3_list[[i]] <- cbind(Filename = i, mat_no3)
  mat_sst_list[[i]] <- cbind(Filename = i, mat_sst)
  mat_ssw_list[[i]] <- cbind(Filename = i, mat_ssw)
  cum_po4_list[[i]] <- cbind(Filename = i, cum_po4)
  cum_no3_list[[i]] <- cbind(Filename = i, cum_no3)
  cum_sst_list[[i]] <- cbind(Filename = i, cum_sst)
  cum_ssw_list[[i]] <- cbind(Filename = i, cum_ssw)
  
  
  mat_list <- list(mat_po4, mat_no3, mat_sst, mat_ssw)
  cum_list <- list(cum_po4, cum_no3, cum_sst, cum_ssw)
  variables <- c("PO4", "NO3", "SST", "SSW")
  
  for (j in 1:4) {
    mat <- mat_list[[j]]
    cum <- cum_list[[j]]
    
    if (all(cum$erh * cum$erl < 0)) {
      lag_start <- -999
      duration <- 0
      effect_sign <- NA
    } else if (all(mat$erh * mat$erl > 0) && all(cum$erh * cum$erl > 0)) {
      lag_start <- 0
      duration <- 16
      effect_sign <- ifelse(mat$erh[1] > 0, "Positive", "Negative")
    } else if (any(cum$erh * cum$erl > 0) && all(mat$erh * mat$erl <= 0)) {
      lag_start <- -999
      duration <- 0
      effect_sign <- NA
    } else {
      valid_days <- which(cum$erh * cum$erl > 0)	
      if (length(valid_days) > 0) {
        if (all(cum$erh[valid_days] > 0) || all(cum$erl[valid_days] < 0)) {
          effect_day <- which(mat$erh * mat$erl > 0)[1]
          lag_start <- effect_day-1
          end_day <- which(mat$erh * mat$erl <= 0)[which(mat$erh * mat$erl <= 0) > effect_day]
          
          if (length(end_day) == 0) {
            duration <- 17 - effect_day  # 
          } else {
            duration <- end_day[1] - effect_day   #
          }
          effect_sign <- ifelse(mat$erh[effect_day] > 0, "Positive", "Negative")
        } else {
          
          negative_days <- valid_days[cum$erh[valid_days] < 0 & cum$erl[valid_days] < 0]
          positive_days <- valid_days[cum$erh[valid_days] > 0 & cum$erl[valid_days] > 0]
          
          if (length(negative_days) > 0) {
            max_negative_erh <- max(abs(cum$erh[negative_days]))
          } else {
            max_negative_erh <- -Inf 
          }
          
          
          if (length(positive_days) > 0) {
            max_positive_erl <- max(cum$erl[positive_days])
          } else {
            max_positive_erl <- -Inf 
          }
          
          
          if ( max_positive_erl > max_negative_erh ) {
            effect_day <- which(mat$erh > 0 & mat$erl > 0)[1]
            lag_start <- effect_day-1
            end_day <- which(mat$erh * mat$erl <= 0)[which(mat$erh * mat$erl <= 0) > effect_day][1]
            
            if (length(end_day) == 0) {
              duration <- 17 - effect_day  
            } else if (is.na(end_day)) {
              end_day <- which(mat$erh < 0 & mat$erl < 0)[which(mat$erh<0 & mat$erl < 0) > effect_day][1]
              if (is.na(end_day) || length(end_day) == 0) {
                duration <- 17 - effect_day 
              } else {
                duration <- end_day[1] - effect_day
              }
            } else {
              duration <- end_day[1] - effect_day   
            }
            effect_sign <- "Positive"
          } else if (max_positive_erl < max_negative_erh) {
            effect_day <- which(mat$erh < 0 & mat$erl < 0)[1]
            lag_start <- effect_day - 1
            end_day <- which(mat$erh * mat$erl <= 0)[which(mat$erh * mat$erl <= 0) > effect_day]
            
            if (length(end_day) == 0) {
              duration <- 17 - effect_day  
            } else if (is.na(end_day)) {
              end_day <- which(mat$erh > 0 & mat$erl > 0)[which(mat$erh>0 & mat$erl > 0) > effect_day][1]
              if (is.na(end_day) || length(end_day) == 0) {
                duration <- 17 - effect_day  
              } else {
                duration <- end_day[1] - effect_day
              }
            } else {
              duration <- end_day[1] - effect_day  
            }
            effect_sign <- "Negative"
          } else {  
            first_positive_day <- which(cum$erh * cum$erl > 0)[1]
            if (cum$erh[first_positive_day] > 0) {
              effect_day <- which(mat$erh > 0 & mat$erl > 0)[1]
              lag_start <- effect_day - 1
              end_day <- which(mat$erh * mat$erl <= 0)[which(mat$erh * mat$erl <= 0) > effect_day]
              
              if (length(end_day) == 0) {
                duration <- 17 - effect_day 
              } else {
                duration <- end_day[1] - effect_day  
              }
              effect_sign <- "Positive"
            } else {
              effect_day <- which(mat$erh < 0 & mat$erl < 0)[1]
              lag_start <- effect_day - 1
              end_day <- which(mat$erh * mat$erl <= 0)[which(mat$erh * mat$erl <= 0) > effect_day]
              
              if (length(end_day) == 0) {
                duration <- 17 - effect_day  
              } else {
                duration <- end_day[1] - effect_day  
              }
              effect_sign <- "Negative"
            }
          }
        }
      } else {
        lag_start <- -999
        duration <- 0
        effect_sign <- NA
      }
    }
    
    lag_results[[length(lag_results) + 1]] <- data.frame(Filename = i, Variable = variables[j], Start_Lag = lag_start, Duration = duration, Effect_Sign = effect_sign)
  }
  
  final_results[[length(final_results) + 1]] <- do.call(rbind, lag_results)
  

  mat_combined <- bind_rows(
    mat_po4 %>% mutate(Label = "PO[4]"),
    mat_no3 %>% mutate(Label = "NO[3]"),
    mat_sst %>% mutate(Label = "SST"),
    mat_ssw %>% mutate(Label = "SSW")
  )
  
  cum_combined <- bind_rows(
    cum_po4 %>% mutate(Label = "PO[4]"),
    cum_no3 %>% mutate(Label = "NO[3]"),
    cum_sst %>% mutate(Label = "SST"),
    cum_ssw %>% mutate(Label = "SSW")
  )
  
  p1 <- ggplot(mat_combined, aes(x = lag, y = erf, color = Label)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = erl, ymax = erh), width = 0.2, size = 1) +
    geom_hline(yintercept = 0, color = "grey20", linetype = "dashed") +
    scale_colour_discrete(name = "Factor") +
    xlab("Lag time (day)") +
    ylab("CER (95% CI)") +
    ggtitle(sprintf("(a) Effect of a single factor in specific lag times", i)) +
    facet_grid(Label ~ ., scales = "free_y", labeller = label_parsed) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.position = 'none',
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          text = element_text(family = "serif"),
          plot.margin = margin(10, 10, 10, 10),
          strip.text = element_text(size = 16, face = "bold"))
  
  p2 <- ggplot(cum_combined, aes(x = lag, color = Label, fill = Label)) +
    geom_ribbon(aes(ymin = erl, ymax = erh), alpha = 0.3, linetype = 0) +
    geom_smooth(aes(y = erf), method = "loess", span = 0.3, se = FALSE, size = 1.1) +
    geom_hline(yintercept = 0, color = "grey20", linetype = "dashed") +
    guides(fill = guide_legend(title = NULL), color = guide_legend(title = NULL)) +
    xlab("Lag time (day)") +
    ylab("CER (95% CI)") +
    ggtitle(sprintf("(b) Effect of a single factor in cumulative lag times", i)) +
    facet_grid(Label ~ ., scales = "free_y", labeller = label_parsed) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.position = 'none',
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          text = element_text(family = "serif"),
          plot.margin = margin(10, 10, 10, 10),
          strip.text = element_text(size = 16, face = "bold"))
  
  
  plot_combined <- plot_grid(p1, p2, ncol = 2)
  
  
  output_file <- sprintf("%splot_combined_%02d.png", output_dir, i)
  ggsave(output_file, plot = plot_combined, width = 12, height = 6, dpi = 300)
  cat(sprintf("Plot saved: %s\n", output_file))
  
  
  return(list(mat_po4 = mat_po4, mat_no3 = mat_no3, mat_sst = mat_sst, mat_ssw = mat_ssw,
              cum_po4 = cum_po4, cum_no3 = cum_no3, cum_sst = cum_sst, cum_ssw = cum_ssw,
              final_results = final_results))
}
## Run the 'process_file' function in parallel on file index 1 to 25.
## Choose a small spatial subset (latitude: -41° to -40°, longitude: 2° to 3°) 
## with a spatial resolution of 0.25°, to demonstrate the structure and output format
## Adjust the numeric range (e.g., 1:10 or 21:40) based on the number and names of files to be processed.
results <- parLapply(cl, 1:25, process_file)

final_results_df <- do.call(rbind, lapply(results, function(x) {
  
  if (length(x$final_results) > 0) {
    return(do.call(rbind, x$final_results))
  } else {
    return(NULL)
  }
}))

head(final_results_df)


lat_seq <- seq(-41, -40, by = 0.25)  
lon_seq <- seq(2, 3, by = 0.25)      

grid <- expand.grid(Lat = lat_seq,Lon = lon_seq) 
grid <- grid[order(grid$Lon, grid$Lat), ] 

stopifnot(nrow(grid) == length(unique(final_results_df$Filename)))

for (var in variables) {
  df_sub <- final_results_df[final_results_df$Variable == var, ]
  df_sub <- df_sub[order(df_sub$Filename), ] 
  # Lagday
  lag_df <- data.frame(Filename = df_sub$Filename,
                       Lat = grid$Lat,
                       Lon = grid$Lon,
                       Start_Lag = df_sub$Start_Lag)
  lag_filename <- sprintf("%slagday_%s.csv", output_dir, tolower(var))
  write.table(lag_df, lag_filename, sep = ",", row.names = FALSE, col.names = FALSE)
  
  # Duration
  dur_df <- data.frame(Filename = df_sub$Filename,
                       Lat = grid$Lat,
                       Lon = grid$Lon,
                       Duration = df_sub$Duration)
  dur_filename <- sprintf("%sduration_%s.csv", output_dir, tolower(var))
  write.table(dur_df, dur_filename, sep = ",", row.names = FALSE, col.names = FALSE)
}