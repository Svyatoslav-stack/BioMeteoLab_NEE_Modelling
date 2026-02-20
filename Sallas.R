> ## -- NEE -- ##
  > `#fluxdata: The chamber measurement data with time stamps, measurement ids and needed environmental variables`
> 
  > # `TmStamp: POSIXct datetime. For example the midpoint of the chamber measurement`
  > # `collar_id: character. The identifier of the measurement point`
  > # `group_id: character. The identifier of collar group (e.g. hummock, hollow)`
  > # `meas_id: character. An identifier for the individual measurement. E.g. TmStamp + collar_id`
  > # `PAR_tr: level of shading: 0 is ambient, 100 is full shade.`
  > # `CO2_C: numeric. The CO2 flux.`
  > # `Reco: The CO2 flux of the full shaded measurement. NEE is measured in sets of 4 with different` 
  > #       `levels of shading. This is the Reco for each set from the same collar (same date and collar_id, subsequent measurements).` 
  > #       `I.e there is a Reco for each of the CO2_C measurement.`
  > # `PAR: PAR measured during the chamber closure, e.g. mean PAR during closure`
  > # `T_soil: Soil temperature measured during the chamber closure`
  > 
  > 
  > # =======NEE========== #
  > # ==================== #
  > `##Light response model##`
> 
  > # `---- Setup ----`
  > `library(dplyr)`
> `library(purrr)`
> `library(stringr)  # for safe filenames, optional`
> 
  > # `---- Plot layout ----`
  > `par(mfrow = c(4, 3), mar = c(4, 4, 2, 1))`
> 
  > `#Set the data groups####`
> `df <- fluxdata`
> `df$date <- as.Date(df$TmStamp)`
> 
  > # `Split by group`
  > `groups <- split(df, df$group_id)` 
> 
  > # ---- Function to fit the model ---- #
  > 
  > `starts_from_data <- function(PAR, NEE) {`
>   # `Low light to estimate R (respiration)`
  >   `q    <- stats::quantile(PAR, probs = 0.2, na.rm = TRUE)`
>   `iLow <- which(PAR <= q)`
>   `if (length(iLow) < 3) iLow <- order(PAR)[seq_len(min(5, length(PAR)))]`
>   
  >   `R0 <- stats::median(NEE[iLow], na.rm = TRUE)`
>   
  >   # `Pseudo-GPP >= 0`
  >   `GPP_proxy <- pmax(R0 - NEE, 0)`
>   
  >   # `Initial slope (through origin) from low PAR`
  >   `slope0 <- try(stats::coef(stats::lm(GPP_proxy[iLow] ~ 0 + PAR[iLow]))[1], silent = TRUE)`
>   `if (!is.finite(slope0)) slope0 <- 0.05`
>   
  >   # `Pmax from high tail of pseudo-GPP (robust to outliers)`
  >   `Pmax0 <- max(stats::quantile(GPP_proxy, 0.95, na.rm = TRUE), na.rm = TRUE)`
>   `if (!is.finite(Pmax0) || Pmax0 <= 0) Pmax0 <- max(GPP_proxy, na.rm = TRUE)`
>   
  >   `list(R = as.numeric(R0),`
>        `alpha = max(as.numeric(slope0), 1e-6),`
>        `Pmax = max(as.numeric(Pmax0), 1e-6))`
> `}`
> 
  > `rmse <- function(obs, fit) sqrt(mean((obs - fit)^2, na.rm = TRUE))`
> `mae  <- function(obs, fit) mean(abs(obs - fit), na.rm = TRUE)`
> 
  > # `---- Fitting the model ----`
  > 
  > `results_list <- list() #Parameters will be collecter here`
> 
  > `par(mfrow=c(3,3))`
> `for (g in names(groups)) {`
>   `dd <- groups[[g]]`
>   
  >   # `Keep only needed cols and drop NA/Inf rows`
  >   `needed <- c("mean_PAR_egm5", "CO2_C_f")`
>   `if (!all(needed %in% names(dd))) next`
>   `dd <- dd[is.finite(dd$mean_PAR_egm5) & is.finite(dd$CO2_C_f), ]`
>   `if (nrow(dd) < 6 || length(unique(dd$mean_PAR_egm5)) < 4) next`
>   
  >   `x <- dd$PAR`
>   `y <- dd$CO2_C_f`
>   
  >   # `Plot the data`
  >   `plot(x, y,`
>        `main = g,`
>        `xlab = "PAR (mean_closure_PAR)",`
>        `ylab = "NEE (CO2_C_f)",`
>        `pch = 19, col = "darkgreen", cex = 0.8)`
>   `abline(h=0, col='grey', lty=2)`
>   
  >   # `Data-driven starting values`
  >   `sv <- starts_from_data(x, y)`
>   
  >   # `Try nls(), capture warnings`
  >   `warn_msgs <- character(0)`
>   `fit <- try(`
>     `withCallingHandlers(`
>       `stats::nls(`
>         `,`
>         `data = dd,`
>         `start = sv,`
>         `algorithm = "port",`
>         `lower = c(R = -Inf, alpha = 1e-8, Pmax = 1e-8),`
>         `upper = c(R =  Inf, alpha = 1e+3, Pmax = 1e+4),`
>         `control = stats::nls.control(maxiter = 300, tol = 1e-8, warnOnly = TRUE, minFactor = 1/2^20)`
>       `),`
>       `warning = function(w) {`
>         `warn_msgs <<- c(warn_msgs, conditionMessage(w))`
>         `invokeRestart("muffleWarning")`
>       `}`
>     `),`
>     `silent = TRUE`
>   `)`
>   
  >   `method_used <- "nls"`
>   
  >   `fits_by_group[[g]] <- fit`
>   
  >   # `Metrics and statistics`
  >   `co   <- stats::coef(fit)`
>   `yhat <- stats::fitted(fit)`
>   `rs   <- y - yhat`
>   `rss  <- sum(rs^2, na.rm = TRUE)`
>   `tss  <- sum((y - mean(y))^2, na.rm = TRUE)`
>   `r2   <- 1 - rss / tss`
>   
  >   `vc  <- safe_vcov(fit)`
>   `seR <- if (!is.null(vc) && "R"     %in% rownames(vc)) sqrt(vc["R",     "R"])     else NA_real_`
>   `seA <- if (!is.null(vc) && "alpha" %in% rownames(vc)) sqrt(vc["alpha", "alpha"]) else NA_real_`
>   `seP <- if (!is.null(vc) && "Pmax"  %in% rownames(vc)) sqrt(vc["Pmax",  "Pmax"])  else NA_real_`
>   
  > 
  >   # `Derived ecophysiological quantities`
  >   `Ic <- NA_real_`
>   `if (is.finite(co["alpha"]) && co["alpha"] > 0 && is.finite(co["Pmax"]) && is.finite(co["R"]) &&`
>       `(co["Pmax"] > co["R"])) {`
>     `Ic <- (co["R"] * co["Pmax"]) / (co["alpha"] * (co["Pmax"] - co["R"]))`
>   `}`
>   `I_half <- if (is.finite(co["alpha"]) && co["alpha"] > 0) co["Pmax"] / co["alpha"] else NA_real_`
>   
  >   `AICv <- try(stats::AIC(fit), silent = TRUE); if (inherits(AICv, "try-error")) AICv <- NA_real_`
>   `BICv <- try(stats::BIC(fit), silent = TRUE); if (inherits(BICv, "try-error")) BICv <- NA_real_`
>   
  >   `s <- summary(fit)`
>   `param_tab <- s$parameters`
>   `p_R     <- if ("R"     %in% rownames(param_tab)) param_tab["R",     "Pr(>|t|)"] else NA_real_`
>   `p_alpha <- if ("alpha" %in% rownames(param_tab)) param_tab["alpha", "Pr(>|t|)"] else NA_real_`
>   `p_Pmax  <- if ("Pmax"  %in% rownames(param_tab)) param_tab["Pmax",  "Pr(>|t|)"] else NA_real_`
>   
  >   # `Extra Sum of Squares`
  >   `n <- length(y)`
>   `p_full <- 3L; p_null <- 1L`
>   `df1 <- p_full - p_null`
>   `df2 <- n - p_full`
>   `RSS1 <- sum((y - yhat)^2, na.rm = TRUE)`
>   `RSS0 <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)`
>   `model_p <- NA_real_`
>   `if (is.finite(RSS0) && is.finite(RSS1) && RSS0 > RSS1 && df2 > 0) {`
>     `Fstat <- ((RSS0 - RSS1) / df1) / (RSS1 / df2)`
>     `if (is.finite(Fstat) && Fstat >= 0) model_p <- stats::pf(Fstat, df1, df2, lower.tail = FALSE)`
>   `}`
>   
  >   # `Plot fitted curve`
  >   `xseq  <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 200)`
>   `ypred <- predict(fit, newdata = data.frame(mean_PAR_egm5 = xseq))`
>   `lines(xseq, ypred, col = "blue", lwd = 2)`
> 
  >   # `---- Store results ----`
  >   `results_list[[g]] <- data.frame(`
>     `group = g,`
>     `plot_id = if ("plot_id" %in% names(dd)) paste(unique(dd$plot_id), collapse = ";") else NA_character_,`
>     `date    = if ("date"    %in% names(dd)) paste(unique(dd$date),    collapse = ";") else NA_character_,`
>     `converged = TRUE,`
>     `method = method_used,`
>     `warning = paste(unique(warn_msgs), collapse = " | "),`
>     `n = nrow(dd),`
>     `par_min = min(x), par_max = max(x), par_span = diff(range(x)),`
>     `R = unname(co["R"]), alpha = unname(co["alpha"]), Pmax = unname(co["Pmax"]),`
>     `se_R = seR, se_alpha = seA, se_Pmax = seP,`
>     `rmse = rmse(y, yhat), mae = mae(y, yhat), r2 = r2,`
>     `AIC = as.numeric(AICv), BIC = as.numeric(BICv),`
>     `Ic = as.numeric(Ic), I_half = as.numeric(I_half),`
>     `p_R = as.numeric(p_R), p_alpha = as.numeric(p_alpha), p_Pmax = as.numeric(p_Pmax),`
>     `model_p = as.numeric(model_p),`
>     `stringsAsFactors = FALSE`
>   `)`
> `}`
> 
  > `fit_summary <- if (length(results_list)) do.call(rbind, results_list) else`
>   `data.frame()`
> 
  > # `Peek at results`
  > `print(utils::head(fit_summary))`
> 
  > 
  > # ==================== #
  > 
  > ## `Each group (e.g. hummocks, hollows) gets a set of model parameters for each measurement day. They should be in the results_list` 
  > ## `These parameters are then interpolated between measurement days to get one parameter value for each of the day of the growing season.`
  > 
  > # =======Reco========== #
  > # ===================== #
  > ## ==Lloyd & Taylor##== #
  > 
  > 
  > # `---- Lloyd & Taylor function ----`
  > `lloyd_taylor_resp <- function(T_C, Rref, E0, Tref_C = 10, T0_K = 227.13) {`
>   `Tk     <- T_C + 273.15`
>   `Tref_K <- Tref_C + 273.15`
>   `Rref * exp(E0 * (1 / (Tref_K - T0_K) - 1 / (Tk - T0_K)))`
> `}`
> 
  > 
  > # `---- Data preparation ----`
  > `df_Reco <- fluxdata[fluxdata$PAR_tr == '100', ]`
> `df_Reco$date <- as.Date(df_Reco$TmStamp)`
> 
  > # `Split by group`
  > `groups_Reco <- split(df_Reco, df_Reco$group_id)`
> 
  > # `---- Helpers ----`
  > `rmse <- function(obs, fit) sqrt(mean((obs - fit)^2, na.rm = TRUE))`
> `mae  <- function(obs, fit) mean(abs(obs - fit), na.rm = TRUE)`
> 
  > `#This is for calculating SE` 
> `bootstrap_se <- function(model, data, nboot = 200) {`
>   `coefs <- replicate(nboot, {`
>     `idx <- sample(seq_len(nrow(data)), replace = TRUE)`
>     `boot_data <- data[idx, ]`
>     `fit_boot <- try(nls(`
>       `Reco ~ lloyd_taylor_resp(mean_soilT_comp, Rref, E0, Tref_C = 10),`
>       `data = boot_data,`
>       `start = as.list(coef(model)),`
>       `control = nls.control(maxiter = 200, tol = 1e-6, warnOnly = TRUE)`
>     `), silent = TRUE)`
>     `if (inherits(fit_boot, "try-error")) return(c(Rref = NA, E0 = NA))`
>     `coef(fit_boot)`
>   `})`
>   `apply(coefs, 1, sd, na.rm = TRUE)`
> `}`
> 
  > # `---- Containers ----`
  > `results_list <- list()`
> `fits_by_group <- list()`
> `par(mfrow = c(3, 3))`
> 
  > # `---- Main loop ----`
  > `for (g in names(groups_Reco)) {`
>   `dd <- groups_Reco[[g]]`
>   
  >   # `Filter valid rows`
  >   `needed <- c("T_soil", "Reco")`
>   `if (!all(needed %in% names(dd))) next`
>   `dd <- dd[is.finite(dd$T_soil) & is.finite(dd$Reco), ]`
>   `if (nrow(dd) < 6 || length(unique(dd$T_soil)) < 4) next`
>   
  >   `x <- dd$T_soil`
>   `y <- dd$Reco`
>   
  >   # `Plot raw data`
  >   `plot(x, y,`
>        `main = g,`
>        `xlab = "Soil Temperature (°C)",`
>        `ylab = "Reco",`
>        `pch = 19, col = "brown", cex = 0.8,` 
>        `ylim=c(0, 450))`
>   
  >   # `Starting values`
  >   `Tref_C <- 10`
>   `near <- which(abs(dd$T_soil - Tref_C) <= 1)`
>   `if (length(near) == 0) near <- which.min(abs(dd$T_soil - Tref_C))`
>   
  >   `Rref_start <- mean(dd$Reco[near], na.rm = TRUE)`
>   `E0_start   <- 200`
>   
  >   # `Fit nls model`
  >   `warn_msgs <- character(0)`
>   `fit <- try(`
>     `withCallingHandlers(`
>       `nls(`
>         `Reco ~ lloyd_taylor_resp(T_soil, Rref, E0, Tref_C = Tref_C),`
>         `data = dd,`
>         `start = list(Rref = Rref_start, E0 = E0_start),`
>         `control = nls.control(maxiter = 200, tol = 1e-6, warnOnly = TRUE)`
>       `),`
>       `warning = function(w) {`
>         `warn_msgs <<- c(warn_msgs, conditionMessage(w))`
>         `invokeRestart("muffleWarning")`
>       `}`
>     `),`
>     `silent = TRUE`
>   `)`
>   
  >   `if (inherits(fit, "try-error")) next`
>   
  >   `fits_by_group[[g]] <- fit`
>   
  >   # `Extract coefficients`
  >   `co <- coef(fit)`
>   `yhat <- fitted(fit)`
>   `rs <- y - yhat`
>   `rss <- sum(rs^2, na.rm = TRUE)`
>   `tss <- sum((y - mean(y))^2, na.rm = TRUE)`
>   `r2 <- 1 - rss / tss`
>   
  >   # `Bootstrap SE`
  >   `se_boot <- bootstrap_se(fit, dd, nboot = 200)`
> 
  >   # `Derived quantities`
  >   `AICv <- try(AIC(fit), silent = TRUE); if (inherits(AICv, "try-error")) AICv <- NA_real_`
>   `BICv <- try(BIC(fit), silent = TRUE); if (inherits(BICv, "try-error")) BICv <- NA_real_`
>   
  >   # `---- Parameter p-values ----`
  >   `s <- summary(fit)`
>   `param_tab <- s$parameters`
>   `p_Rref <- if ("Rref" %in% rownames(param_tab)) param_tab["Rref", "Pr(>|t|)"] else NA_real_`
>   `p_E0   <- if ("E0"   %in% rownames(param_tab)) param_tab["E0",   "Pr(>|t|)"] else NA_real_`
>   
  >   # `---- Model-level p-value (F-test) ----`
  >   `n <- length(y)`
>   `p_full <- 2L; p_null <- 1L`
>   `df1 <- p_full - p_null`
>   `df2 <- n - p_full`
>   `RSS1 <- sum((y - yhat)^2, na.rm = TRUE)`
>   `RSS0 <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)`
>   `model_p <- NA_real_`
>   `if (is.finite(RSS0) && is.finite(RSS1) && RSS0 > RSS1 && df2 > 0) {`
>     `Fstat <- ((RSS0 - RSS1) / df1) / (RSS1 / df2)`
>     `if (is.finite(Fstat) && Fstat >= 0) model_p <- stats::pf(Fstat, df1, df2, lower.tail = FALSE)`
>   `}`
>   
  >   # `Plot fitted curve`
  >   `xseq <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 200)`
>   `ypred <- lloyd_taylor_resp(xseq, co["Rref"], co["E0"], Tref_C = Tref_C)`
>   `lines(xseq, ypred, col = "blue", lwd = 2)`
>   
  >   # `Store results`
  >   `results_list[[g]] <- data.frame(`
>     `group = g,`
>     `plot_id = if ("plot_id" %in% names(dd)) paste(unique(dd$plot_id), collapse = ";") else NA_character_,`
>     `date = if ("date" %in% names(dd)) paste(unique(dd$date), collapse = ";") else NA_character_,`
>     `converged = TRUE,`
>     `method = "nls",`
>     `warning = paste(unique(warn_msgs), collapse = " | "),`
>     `n = nrow(dd),`
>     `temp_min = min(x), temp_max = max(x), temp_span = diff(range(x)),`
>     `Rref = unname(co["Rref"]), E0 = unname(co["E0"]),`
>     `se_Rref_boot = se_boot["Rref"], se_E0_boot = se_boot["E0"],`
>     `rmse = rmse(y, yhat), mae = mae(y, yhat), r2 = r2,`
>     `AIC = as.numeric(AICv), BIC = as.numeric(BICv),`
>     `p_Rref = as.numeric(p_Rref), p_E0 = as.numeric(p_E0),`
>     `model_p = as.numeric(model_p),`
>     `stringsAsFactors = FALSE`
>   `)`
> `}`
> 
  > `fit_summary_Reco <- if (length(results_list)) do.call(rbind, results_list) else data.frame() #this collects the parameters etc.`
> 
  > ##
  > # ===NEE predictions=== #
  > # ===================== #
  > 
  > # `To predict NEE between the measurements, you will need a continuous time series of PAR and T_soil (preferably 30 min data). It is critical that the continuous measurements are harmonized to the PAR and T_soil used in the fitting of the models above,` 
  > # `because these models are not linear, and can be very sensitive. If they are not harmonized, the predictions may under-or overertimate the components of NEE (Reco, GPP).`
  > 
  > ## `data: hhourly PAR and T_soil, Reco parameters (Rref and E0), NEE parameters (alpha, Pmax)  interpolated between measurement days, for each of the collar groups separately`
  > 
  > `alpha = data$alpha_ip #interpolated alpha between measurement days`
> `Pmax = data$$Pmax_ip #interpolated Pmax between measurement days`
> `PAR = data$PAR`
> `Tsoil_K = data$T_soil+273.15` 
> `Rref = data$Rref`
> `E0 = data$E0`
> 
  > `data$NEE_pred <- (Rref * exp(E0 * (1 / (Tref_K - T0_K) - 1 / (Tsoil_K - T0_K)))) - ((alpha * Pmax * PAR) / (Pmax + alpha * PAR))`
> `data$photos_pred <- -1*((alpha * Pmax * PAR) / (Pmax + alpha * PAR)) #atmospheric convention: negative is uptake by photosynthesis`
> `data$Reco_pred <- (Rref * exp(E0 * (1 / (Tref_K - T0_K) - 1 / (Tsoil_K - T0_K))))` 
> 
  > `##Make sure your sign conventions make sense, make sure you know the structure of your data, think through how to group your data. Make sure you know what the unit of your flux is, and what the temporal resolution of your` 
> `# continuous time series is, and take those into acsount when compiling growing season balances.`
> 
  > 
  > 