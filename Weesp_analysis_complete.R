# Angstel/Vecht bog-oak chronology (1283 BCE - 156 CE) ####
# Analysis code for Arts and Bazelmans, The Holocene
# The script runs from top to bottom and writes the figures to pdf_figures/ and png_figures/
# and the tables to tables/, in the order in which they appear in the paper:
#   Main text: Figure 1, Figure 2, Table 1, Figure 3, Figure 4, persistence test (Section 4.3),
#              Figure 5, Figure 6, Figure 7
#   Supplement: Figure S1, Table S5, Figure S2, Figure S3, Table S11, Table S12,
#               fitted delay (S7), Figure S4, Table S14, Figure S5, Table S15
# In RStudio each figure and table is a section in the document outline (Ctrl+Shift+O)
# Each section is wrapped in { } so it can be folded, and run as one block by placing
# the cursor on the opening { and pressing Ctrl+Enter. Run the settings, functions and data first.

# Files needed in D:/Phd/documents/Weesp (data_dir below):
#   10_Angstel_Vecht_driechronos_std050_vsc101_constant040_16082026.xlsx
#   (or AV_diagnostic_data.csv and AV_constant40.csv), the growth index and its diagnostics
#   stal_scot.csv, Uamh an Tartair growth rate (Baker et al., 2015)
#   sawtooth2020.txt, Sawtooth Lake AMV index (Lapointe et al., 2020)
#   palgeo.jpg, the map of Figure 1
# The two comparison records are downloaded if they are missing.
# Completed surrogate runs are saved in saved_results/ and are not rerun (see settings).


# Packages and settings ####
{
  library(readxl)
  library(astrochron)     # taner(), linterp(), sortNave(), mtmML96(), mtmPL()
  have_multitaper <- requireNamespace("multitaper", quietly = TRUE)
  
  set.seed(1337)
  
  # Folder with the data; figures and tables are written to subfolders of it
  data_dir <- "D:/Phd/documents/Weesp"
  if (dir.exists(data_dir)) setwd(data_dir)
  for (d in c("pdf_figures", "png_figures", "tables", "saved_results"))
    if (!dir.exists(d)) dir.create(d)
  
  data_file  <- "10_Angstel_Vecht_driechronos_std050_vsc101_constant040_16082026.xlsx"
  chronology <- "std050"        # "std050" | "vsc101" | "constant40_median"
  year_min   <- -1283           # analysis window
  year_max   <-  156
  
  # Surrogate counts used for the published results (draft runs used 1000, 500 and 5000)
  n_surrogates          <- 10000   # Table 1, Figure 3
  n_surrogate_fields    <- 2000    # persistence test (Section 4.3)
  n_surrogates_bandpass <- 20000   # band-pass correlations (Figures S3, S4, Tables S12, S14)
  n_surrogates_step     <- 500    # change near 572 BCE (Figure S5)
  
  # Surrogate runs that are completed are saved in saved_results/ and loaded on the next run,
  # so they are not rerun. A saved run is only reused for the same chronology and surrogate count;
  # delete the file (or set use_saved_results <- FALSE) to rerun it.
  use_saved_results <- TRUE
  
  # Wavelet settings
  omega     <- 10
  dj_figure <- 1 / 100
  dj_test   <- 1 / 20             # coarser grid for the surrogates
  upper_period <- 256
  
  # Candidate periods tested, with a +/- 10% window
  candidate_periods <- c(2.2, 3.6, 5.9, 8, 12, 15.4, 24, 38, 60, 80, 110, 200)
  candidate_window  <- 0.10
  
  # Bands removed before estimating the signal-corrected AR(1) coefficient.
  # The 2.2 and 200 yr bands are left out: removing the 2.2 yr band raises phi
  notch_periods <- c(3.6, 5.9, 8, 12, 15.4, 24, 38, 60, 80, 110)
  
  # Octave bands for the persistence test
  octave_bands <- data.frame(
    Band = c("2-4", "4-8", "8-16", "16-32", "32-64", "64-128"),
    Period_min = c(2, 4, 8, 16, 32, 64),
    Period_max = c(4, 8, 16, 32, 64, 128),
    stringsAsFactors = FALSE)
}


# Functions ####

## Progress bars and saved results for the surrogate loops ####
{
  # Run expr once and save the result, or load it when it was saved before.
  # The random number state after the run is saved with it, so later sections
  # give the same result whether a run is loaded or recomputed
  run_or_load <- function(name, expr) {
    file <- file.path("saved_results", paste0(name, ".rds"))
    if (use_saved_results && file.exists(file)) {
      z <- readRDS(file)
      assign(".Random.seed", z$seed, envir = globalenv())
      cat("\nLoaded", file, "(delete it to rerun)\n")
      return(z$value)
    }
    value <- expr
    saveRDS(list(value = value, seed = get(".Random.seed", envir = globalenv())), file)
    value
  }
  
  # Start a progress bar with a label
  progress_start <- function(n, label) {
    cat("\n", label, " (", n, " surrogates)\n", sep = "")
    flush.console()
    txtProgressBar(min = 0, max = n, style = 3)
  }
  
  # replicate() with a progress bar
  replicate_progress <- function(n, label, expr) {
    f <- eval.parent(substitute(function(...) expr))
    pb <- progress_start(n, label); k <- 0
    out <- sapply(integer(n), function(...) {
      r <- f(); k <<- k + 1; setTxtProgressBar(pb, k); r
    }, simplify = "array")
    close(pb)
    out
  }
}

## Plotting: colours, colour key, cone of influence, BCE/CE axis ####
# Colours: scico batlowK with quantile colour breaks, as in WaverideR
{
  palette_name <- "batlowK"
  
  batlow_colours <- function(n = 100) {
    if (requireNamespace("scico", quietly = TRUE))
      scico::scico(n = n, palette = palette_name)
    else hcl.colors(n, "YlOrRd", rev = TRUE)
  }
  
  quantile_breaks <- function(z, n = 100)
    quantile(z, probs = seq(0, 1, length.out = n + 1), na.rm = TRUE)
  
  # Colour key for the scalograms
  draw_colour_key <- function(breaks, cols = batlow_colours(100), lab = "Power",
                              n_lab = 4, cex = 0.55, frac = 0.42,
                              bar_h = 0.22, y0 = 0.52, align = "left",
                              lab_pos = "right") {
    nb <- length(breaks) - 1
    par(xaxs = "i", yaxs = "i")
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    x0 <- switch(align, left = 0.02, centre = (1 - frac) / 2, right = 1 - frac - 0.02)
    xs <- seq(x0, x0 + frac, length.out = nb + 1)
    rect(xs[-(nb + 1)], y0, xs[-1], y0 + bar_h, col = cols, border = NA)
    rect(x0, y0, x0 + frac, y0 + bar_h, border = "grey25", lwd = 0.7)
    at <- round(seq(1, nb + 1, length.out = n_lab))
    text(x = xs[at], y = y0 - 0.16, labels = signif(breaks[at], 2),
         cex = cex, xpd = NA)
    if (lab_pos == "right") {
      text(x = x0 + frac + 0.015, y = y0 + bar_h / 2, labels = lab,
           cex = cex, adj = 0, xpd = NA)
    } else {
      text(x = x0 + frac / 2, y = y0 - 0.30, labels = lab, cex = cex, xpd = NA)
    }
  }
  
  shade_coi <- function(xs, coi, period_min, period_max) {
    y <- log2(pmax(pmin(coi, period_max), period_min))
    polygon(c(xs, rev(xs)),
            c(y, rep(log2(period_max), length(xs))),
            col = rgb(1, 1, 1, 0.55), border = NA)
    lines(xs, y, lwd = 1.2, lty = 2, col = "grey20")
  }
  
  # BCE/CE axis, there is no year 0 so BCE years are plotted at year + 1
  add_bce_ce_axis <- function(by = 250) {
    labs <- seq(-1250, 250, by = by)
    labs <- labs[!(labs %in% c(-1, 0, 1))]
    at <- ifelse(labs < 0, labs + 1, labs)
    minor <- seq(-1300, 200, by = 50)
    axis(1, at = ifelse(minor < 0, minor + 1, minor), labels = FALSE, tck = -0.02)
    axis(1, at = at, labels = labs, tck = -0.05)
    axis(1, at = 1, labels = "1", font = 2, tck = -0.05)
  }
  
  # Panel label in the upper left corner
  add_label <- function(lab) {
    usr <- par("usr")
    text(usr[1] + 0.02 * (usr[2] - usr[1]),
         usr[4] - 0.075 * (usr[4] - usr[3]), lab, font = 2, cex = 1.2)
  }
}


## Reading the chronology and its diagnostics ####
{
  read_chronology <- function(proxy = chronology) {
    
    if (file.exists(data_file)) {
      d <- as.data.frame(read_excel(data_file, sheet = "Diagnostic_data"))
      names(d)[names(d) == "samp.depth"] <- "samp_depth"
      if (proxy == "constant40_median") {
        d2 <- as.data.frame(read_excel(data_file, sheet = "Constant40"))
        return(data.frame(Year = d2$Year, GI = d2$constant40_median,
                          samp_depth = d2$Original_n))
      }
      return(data.frame(Year = d$Year, GI = d[[proxy]],
                        samp_depth = d$samp_depth))
    }
    
    d <- read.csv("AV_diagnostic_data.csv")
    names(d)[2] <- "samp_depth"
    data.frame(Year = d$Year, GI = d[[proxy]], samp_depth = d$samp_depth)
  }
  
  read_diagnostics <- function() {
    if (file.exists(data_file)) {
      d  <- as.data.frame(read_excel(data_file, sheet = "Diagnostic_data"))
      names(d)[names(d) == "samp.depth"] <- "samp_depth"
      c40 <- as.data.frame(read_excel(data_file, sheet = "Constant40"))
      d$constant40_median <- c40$constant40_median[match(d$Year, c40$Year)]
    } else {
      d <- read.csv("AV_diagnostic_data.csv")
      names(d)[names(d) == "samp.depth"] <- "samp_depth"
    }
    d
  }
  
  running_variance <- function(x, k = 101) {
    w  <- rep(1 / k, k)
    m1 <- as.numeric(stats::filter(x, w, sides = 2))
    m2 <- as.numeric(stats::filter(x^2, w, sides = 2))
    (m2 - m1^2) * k / (k - 1)
  }
}


## Wavelet transform and average spectral power ####
# Morlet transform, used for the data and for every surrogate
{
  morlet_wavelet <- function(x, dt = 1, dj = dj_test, omega0 = omega,
                             period_min = 2, period_max = upper_period) {
    
    n <- length(x)
    x <- x - mean(x)
    sdx <- sd(x)
    
    fac <- (4 * pi) / (omega0 + sqrt(2 + omega0^2))
    s0 <- period_min * dt / fac
    smax <- period_max * dt / fac
    scales <- s0 * 2^((0:floor(log2(smax / s0) / dj)) * dj)
    
    N <- 2^ceiling(log2(n))
    xpad <- c(x, rep(0, N - n))
    k <- (1:floor(N / 2)) * (2 * pi) / (N * dt)
    k <- c(0, k, -k[(floor((N - 1) / 2)):1])
    f <- fft(xpad)
    
    W <- matrix(complex(real = 0, imaginary = 0),
                nrow = length(scales), ncol = n)
    
    for (i in seq_along(scales)) {
      s <- scales[i]
      daughter <- sqrt(2 * pi * s / dt) * pi^(-0.25) *
        exp(-((s * k - omega0)^2) / 2 * (k > 0)) * (k > 0)
      W[i, ] <- fft(f * daughter, inverse = TRUE)[1:n] / N
    }
    
    period <- scales * fac
    edge <- pmin(seq_len(n) - 1, n - seq_len(n))
    edge[edge == 0] <- 1e-5
    coi <- fac * dt * edge
    
    # power is normalised by the variance and rectified by scale (Liu et al., 2007)
    list(W = W, power = (Mod(W)^2) / (sdx^2 * scales), period = period,
         scales = scales, coi = coi, valid = outer(period, coi, "<="), n = n)
  }
  
  average_spectral_power <- function(w)
    vapply(seq_along(w$period), function(i) {
      v <- w$valid[i, ]
      if (!any(v)) NA_real_ else mean(w$power[i, v])
    }, numeric(1))
}


## Background (noise) models ####
{
  ar1_coefficient <- function(x)
    as.numeric(arima(x, order = c(1, 0, 0), include.mean = TRUE,
                     method = "ML")$coef["ar1"])
  
  remove_bands <- function(x, periods, tol = candidate_window) {
    n <- length(x); mu <- mean(x)
    f <- fft(x - mu); freq <- (0:(n - 1)) / n
    drop <- rep(FALSE, n)
    for (p in periods) {
      lo <- 1 / (p * (1 + tol)); hi <- 1 / (p * (1 - tol))
      drop <- drop | (freq >= lo & freq <= hi) |
        (freq >= 1 - hi & freq <= 1 - lo)
    }
    f[drop] <- 0
    Re(fft(f, inverse = TRUE) / n) + mu
  }
  
  simulate_white <- function(n, par) rnorm(n, par$mean, par$sd)
  
  simulate_ar1 <- function(n, par)
    as.numeric(arima.sim(list(ar = par$phi), n = n,
                         sd = par$sd * sqrt(1 - par$phi^2))) + par$mean
  
  simulate_powerlaw <- function(n, par) {
    nn <- 2^ceiling(log2(n))
    fr <- (1:(nn / 2)) / nn
    sp <- complex(modulus = fr^(-par$beta / 2),
                  argument = runif(nn / 2, 0, 2 * pi))
    y <- Re(fft(c(0, sp, Conj(rev(sp))[-1]), inverse = TRUE))[1:n]
    par$mean + par$sd * (y - mean(y)) / sd(y)
  }
  
  # The power-law exponent is fitted to an unrectified spectrum (the multitaper spectrum),
  # fitted to the rectified spectrum it becomes negative
  fit_powerlaw_exponent <- function(period, power_rect, scales, x) {
    
    if (have_multitaper) {
      nfft <- 2^ceiling(log2(length(x)))
      fit <- multitaper::spec.mtm(timeSeries = as.numeric(x), nw = 4, k = 7,
                                  nFFT = nfft, returnZeroFreq = FALSE, Ftest = FALSE, jackknife = FALSE,
                                  adaptiveWeighting = TRUE, plot = FALSE, deltat = 1, dtUnits = "year")
      per <- 1 / as.numeric(fit$freq)
      pw  <- as.numeric(fit$spec)
    } else {
      per <- period
      pw  <- power_rect * scales
    }
    
    keep <- is.finite(pw) & pw > 0 & per >= 2 & per <= 128
    for (p in notch_periods)
      keep <- keep & !(per >= p * (1 - candidate_window) & per <= p * (1 + candidate_window))
    
    freq <- 1 / per
    -as.numeric(coef(lm(log(pw[keep]) ~ log(freq[keep])))[2])
  }
  
  make_backgrounds <- function(x, period, gws, scales) {
    
    phi_obs  <- ar1_coefficient(x)
    phi_corr <- ar1_coefficient(remove_bands(x, notch_periods))
    beta     <- fit_powerlaw_exponent(period, gws, scales, x)
    
    cat("\nBackground models\n")
    cat("  AR(1) phi, observed series      :", round(phi_obs, 3), "\n")
    cat("  AR(1) phi, candidate bands out  :", round(phi_corr, 3), "\n")
    cat("  power-law beta                  :", round(beta, 3), "\n")
    
    list(
      white = list(sim = simulate_white,
                   par = list(mean = mean(x), sd = sd(x)),
                   label = "white noise"),
      ar1_observed = list(sim = simulate_ar1,
                          par = list(mean = mean(x), sd = sd(x), phi = phi_obs),
                          label = sprintf("AR(1) phi=%.3f (observed)", phi_obs)),
      ar1_corrected = list(sim = simulate_ar1,
                           par = list(mean = mean(x), sd = sd(x), phi = phi_corr),
                           label = sprintf("AR(1) phi=%.3f (signal-corrected)",
                                           phi_corr)),
      powerlaw = list(sim = simulate_powerlaw,
                      par = list(mean = mean(x), sd = sd(x), beta = beta),
                      label = sprintf("power law beta=%.3f", beta))
    )
  }
}


## Significance tests ####
# 95% level per scale and maximum-in-band p-values of the average spectral power
{
  surrogate_significance <- function(bg, nsim = n_surrogates) {
    
    nper <- length(w_test$period)
    crit <- matrix(NA_real_, nrow = nper, ncol = nsim)
    gws_sim <- matrix(NA_real_, nrow = nper, ncol = nsim)
    
    pb <- progress_start(nsim, paste("Table 1, background:", bg$label))
    for (s in seq_len(nsim)) {
      ws <- morlet_wavelet(bg$sim(n, bg$par), dj = dj_test)
      crit[, s] <- vapply(seq_len(nper), function(i) {
        v <- ws$valid[i, ]
        if (!any(v)) NA_real_ else quantile(ws$power[i, v], 0.95, names = FALSE)
      }, numeric(1))
      gws_sim[, s] <- average_spectral_power(ws)
      setTxtProgressBar(pb, s)
    }
    close(pb)
    
    cand <- do.call(rbind, lapply(candidate_periods, function(p) {
      rr <- which(w_test$period >= p * (1 - candidate_window) &
                    w_test$period <= p * (1 + candidate_window))
      if (!length(rr)) return(NULL)
      obs_max <- max(gws_obs[rr], na.rm = TRUE)
      sim_max <- apply(gws_sim[rr, , drop = FALSE], 2, max, na.rm = TRUE)
      data.frame(Candidate = p,
                 Peak_period = w_test$period[rr][which.max(gws_obs[rr])],
                 Peak_power = obs_max,
                 p_max_in_window = mean(sim_max >= obs_max),
                 Background = bg$label,
                 stringsAsFactors = FALSE)
    }))
    
    list(crit_scale = rowMeans(crit, na.rm = TRUE),
         gws95 = apply(gws_sim, 1, quantile, probs = 0.95, na.rm = TRUE),
         cand = cand)
  }
  
  # Longest run of years above the 95% level, for the persistence test
  longest_run_above <- function(v, crit) {
    r <- rle(v > crit); i <- which(r$values)
    if (!length(i)) 0L else max(r$lengths[i])
  }
}


## Cycle extraction and share of spectral power ####
{
  power_share <- function(w, cycle, up = 1.2, down = 0.8) {
    
    rr <- which(w$period >= cycle * down & w$period <= cycle * up)
    band <- rep(NA_real_, w$n)
    total <- rep(NA_real_, w$n)
    
    for (j in seq_len(w$n)) {
      vb <- rr[w$valid[rr, j]]
      va <- which(w$valid[, j])
      if (length(vb)) band[j] <- mean(w$power[vb, j])
      if (length(va)) total[j] <- sum(w$power[va, j])
    }
    
    frac <- band * length(rr) / total
    frac[!is.finite(frac)] <- 0
    # mask years where fewer than half the periods lie outside the cone of influence
    nvalid <- colSums(w$valid)
    frac[nvalid < 0.5 * nrow(w$valid)] <- NA
    cbind(x = growth$PlotYear, band = band, total = total, fraction = frac)
  }
  
  taner_cycle <- function(centre, lo, hi, roll = 10^20) {
    taner(VA25mean, flow = 1 / hi, fhigh = 1 / lo,
          roll = roll, xmax = 1 / (centre / 3),
          detrend = FALSE, genplot = FALSE)
  }
  
  # Number of backgrounds under which a cycle is significant (from Table 1)
  n_backgrounds_passed <- function(centre) {
    if (!exists("cycle_table")) return("")
    z <- cycle_table[abs(cycle_table$Candidate - centre) < 1e-6, ]
    if (!nrow(z)) return("")
    npass <- sum(z$p_max_in_window < 0.05)
    paste0(" [", npass, "/4 backgrounds]")
  }
}


## Bicoherence and combination tones ####
# Bicoherence after Kim and Powers (1979), based on the astrochron function,
# conditioning factor at the 75th percentile of the denominator
{
  bicoherence <- function(dat, segments = 8, overlap = 50, padfac = 50,
                          CL = 90, maxF = 1 / 5,
                          demean = TRUE, detrend = TRUE, taper = TRUE) {
    
    d <- data.frame(dat)
    npts <- nrow(d)
    dt <- abs(d[2, 1] - d[1, 1])
    
    segpts <- floor(npts / (segments - (segments - 1) * overlap / 100))
    padpts <- segpts * padfac
    if (padpts %% 2 != 0) padpts <- padpts + 1
    df <- 1 / (padpts * dt)
    freq <- df * (seq_len(padpts) - 1)
    
    ft <- matrix(NA_complex_, nrow = padpts, ncol = segments)
    ihold <- 1
    for (i in seq_len(segments)) {
      dd <- d[ihold:(ihold + segpts - 1), ]
      if (demean) dd[, 2] <- dd[, 2] - mean(dd[, 2])
      if (detrend) dd[, 2] <- residuals(lm(dd[, 2] ~ dd[, 1]))
      if (taper) dd[, 2] <- dd[, 2] *
          0.5 * (1 - cos(2 * pi * seq_len(segpts) / (segpts + 1)))
      ft[, i] <- fft(c(dd[, 2], rep(0, padpts - segpts))) / segpts
      ihold <- ihold + floor(segpts * (100 - overlap) / 100)
    }
    
    pwr <- rowMeans(Mod(ft)^2)
    nfreq <- floor(maxF / df)
    
    B <- matrix(0 + 0i, nfreq, nfreq)
    f1f2 <- matrix(0, nfreq, nfreq)
    f3 <- matrix(0, nfreq, nfreq)
    
    for (ii in seq_len(segments)) {
      for (i in seq_len(nfreq)) {
        B[i, ] <- B[i, ] + ft[i, ii] * ft[seq_len(nfreq), ii] *
          Conj(ft[i + seq_len(nfreq) - 1, ii])
        f1f2[i, ] <- f1f2[i, ] + Mod(ft[i, ii] * ft[seq_len(nfreq), ii])^2
        f3[i, ] <- f3[i, ] + Mod(ft[i + seq_len(nfreq) - 1, ii])^2
      }
    }
    
    B <- B / segments; f1f2 <- f1f2 / segments; f3 <- f3 / segments
    den <- f1f2 * f3
    CF <- sort(den)[floor(length(den) * 0.75)]
    bic2 <- Mod(B)^2 / (den + CF)
    
    list(freq = freq[seq_len(nfreq)], bic2 = bic2, pwr = pwr[seq_len(nfreq)],
         CLset = qchisq(CL / 100, df = 2) / (2 * segments),
         CF = CF, CL = CL, segments = segments, overlap = overlap)
  }
  
  # Band-pass the two parents, multiply, band-pass the product at the sum frequency,
  # and compare with the cycle extracted directly from the record
  combination_tone <- function(P1, P2, bw = 0.10) {
    
    # 1/60 + 1/100 -> 37.5 yr and 1/24 + 1/40 -> 15.0 yr
    f1 <- 1 / P1; f2 <- 1 / P2
    f_sum <- f1 + f2
    Nyq <- 0.5
    
    bp <- function(dat, f, bw) {
      taner(dat, padfac = 2, flow = f - bw * f, fhigh = f + bw * f,
            roll = 10^3, demean = TRUE, detrend = FALSE, addmean = FALSE,
            output = 1, xmin = 0, xmax = Nyq, genplot = FALSE,
            check = TRUE, verbose = FALSE)
    }
    
    b1 <- bp(VA25mean, f1, bw)          # parent 1
    b2 <- bp(VA25mean, f2, bw)          # parent 2
    direct <- bp(VA25mean, f_sum, bw)   # daughter, straight from the record
    
    prod <- cbind(b1[, 1], b1[, 2] * b2[, 2])   # quadratic interaction
    synth <- bp(prod, f_sum, bw)                # daughter, synthesised
    
    r <- cor(direct[, 2], synth[, 2])
    
    cat(sprintf("\n%g + %g -> %.1f yr  |  r(direct, synthetic) = %+.3f\n",
                P1, P2, 1 / f_sum, r))
    
    list(parent1 = b1, parent2 = b2, direct = direct, synthetic = synth,
         period = 1 / f_sum, r = r)
  }
}


## Cross-wavelet transform ####
# After analyze_Xwavelet() in WaverideR (Arts, 2023), cross-power is rectified by scale
{
  cross_wavelet <- function(data_1, data_2, dj = dj_figure,
                            lowerPeriod = 2, upperPeriod = 256,
                            omega_nr = omega) {
    
    data_1 <- data_1[order(data_1[, 1]), ]
    data_2 <- data_2[order(data_2[, 1]), ]
    
    xmin <- max(min(data_1[, 1]), min(data_2[, 1]))
    xmax <- min(max(data_1[, 1]), max(data_2[, 1]))
    dx <- max(median(diff(data_1[, 1])), median(diff(data_2[, 1])))
    grid <- seq(xmin, xmax, by = dx)
    
    y1 <- approx(data_1[, 1], data_1[, 2], grid, rule = 1)$y
    y2 <- approx(data_2[, 1], data_2[, 2], grid, rule = 1)$y
    
    wt1 <- morlet_wavelet(y1, dt = dx, dj = dj, omega0 = omega_nr,
                          period_min = lowerPeriod, period_max = upperPeriod)
    wt2 <- morlet_wavelet(y2, dt = dx, dj = dj, omega0 = omega_nr,
                          period_min = lowerPeriod, period_max = upperPeriod)
    
    nr <- length(wt1$period); nc <- length(grid)
    
    # rectify by scale
    Wave.xy <- (wt1$W * Conj(wt2$W)) / wt1$scales
    Power.xy <- Mod(Wave.xy)
    Phase.xy <- Arg(Wave.xy)
    
    list(Wave = t(Wave.xy), Power = t(Power.xy), Phase = t(Phase.xy),
         Power.avg = rowMeans(as.matrix(Power.xy)),
         dt = dx, dj = dj, Scale = wt1$scales, Period = wt1$period,
         nc = nc, nr = nr,
         axis.1 = grid, axis.2 = log2(wt1$period),
         omega_nr = omega_nr,
         x1 = grid, y1 = y1, x2 = grid, y2 = y2,
         coi = pmin(wt1$coi, wt2$coi),
         valid = wt1$valid & wt2$valid)
  }
  
  # Cross-wavelet for plotting, computed and drawn in two steps so the colour key uses the same breaks,
  # phase arrows are placed on ridges of cross-power
  cross_wavelet_for_plot <- function(A, B, lo = 2, hi = 192, q = 0.85) {
    X <- cross_wavelet(cbind(A$x, A$y), cbind(B$x, B$y),
                       dj = dj_figure, lowerPeriod = 1.5,
                       upperPeriod = hi, omega_nr = 15)
    rows <- which(X$Period >= lo & X$Period <= hi)
    z <- t(t(X$Power)[rows, ])
    list(X = X, rows = rows, z = z, brk = quantile_breaks(z), lo = lo, hi = hi, q = q)
  }
  
  plot_cross_wavelet <- function(K, ylab = "Period (years)") {
    X <- K$X; rows <- K$rows
    xs <- X$axis.1; ys <- log2(X$Period[rows])
    image(xs, ys, K$z, col = batlow_colours(100), breaks = K$brk, useRaster = TRUE,
          axes = FALSE, xlab = "", ylab = ylab, cex.lab = 0.9,
          ylim = range(ys), xaxs = "i")
    shade_coi(xs, X$coi, K$lo, K$hi)
    at <- 2^(1:8); at <- at[at >= K$lo & at <= K$hi]
    axis(2, at = log2(at), labels = at, las = 1, cex.axis = 0.8)
    add_bce_ce_axis(250)
    box()
    # one arrow per cell of a time-period grid, at the strongest ridge
    Pw <- t(X$Power); nr <- X$nr
    up <- rbind(Pw[-1, , drop = FALSE], -Inf)
    dn <- rbind(Inf, Pw[-nr, , drop = FALSE])
    md <- Pw > up & Pw > dn; md[1, ] <- FALSE; md[nr, ] <- FALSE
    sel <- md & (Pw >= quantile(X$Power, K$q, na.rm = TRUE)) & X$valid
    Ph <- t(X$Phase)
    nt <- 42; np <- 16
    tbin <- cut(seq_len(X$nc), breaks = nt, labels = FALSE)
    pbin <- cut(log2(X$Period), breaks = np, labels = FALSE)
    for (tb in seq_len(nt)) for (pb in seq_len(np)) {
      cand <- which(sel & outer(pbin == pb, tbin == tb, "&"), arr.ind = TRUE)
      cand <- cand[cand[, 1] %in% rows, , drop = FALSE]
      if (!nrow(cand)) next
      k <- cand[which.max(Pw[cand]), , drop = FALSE]
      r <- k[1]; cc <- k[2]
      dx <- diff(range(xs)) * 0.014 * cos(Ph[r, cc])
      dy <- diff(range(ys)) * 0.028 * sin(Ph[r, cc])
      arrows(xs[cc] - dx, log2(X$Period[r]) - dy,
             xs[cc] + dx, log2(X$Period[r]) + dy,
             length = 0.020, lwd = 1.0, col = "black")
    }
  }
  
  # Phase by period: power-weighted mean outside the cone of influence,
  # positive when the comparison record leads
  edges <- c(28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 70)
  phase_by_period <- function(A, B) {
    X <- cross_wavelet(cbind(A$x, A$y), cbind(B$x, B$y), dj = dj_figure,
                       lowerPeriod = 1.5, upperPeriod = 192, omega_nr = 15)
    Pw <- t(X$Power); Ph <- t(X$Phase)
    t(vapply(seq_len(length(edges) - 1), function(k) {
      rows <- which(X$Period >= edges[k] & X$Period < edges[k + 1])
      v <- X$valid[rows, ]; w <- Pw[rows, ][v]
      wu <- (Pw * X$Period)[rows, ][v]    # unrectified
      z  <- sum(w  * exp(1i * Ph[rows, ][v])) / sum(w)
      zu <- sum(wu * exp(1i * Ph[rows, ][v])) / sum(wu)
      c(phase = Arg(z) * 180 / pi, R = Mod(z), power = mean(w),
        phase_unrect = Arg(zu) * 180 / pi, power_unrect = mean(wu))
    }, numeric(5)))
  }
}


## Band-pass comparison of two records ####
{
  download_if_missing <- function(url, file) {
    if (!file.exists(file))
      try(download.file(url, file, quiet = TRUE), silent = TRUE)
    file.exists(file)
  }
  
  bandpass <- function(dat, lo, hi)
    taner(dat, padfac = 2, flow = 1 / hi, fhigh = 1 / lo, roll = 10^3,
          demean = TRUE, detrend = FALSE, addmean = FALSE, output = 1,
          xmin = 0, xmax = 0.5, genplot = FALSE, check = TRUE, verbose = FALSE)
  
  # Phase randomisation, keeps the spectrum (and autocorrelation) of the series
  phase_randomise <- function(v) {
    n <- length(v); f <- fft(v); k <- 2:floor(n / 2)
    f[k] <- Mod(f[k]) * exp(1i * runif(length(k), 0, 2 * pi))
    f[n - k + 2] <- Conj(f[k])
    Re(fft(f, inverse = TRUE)) / n
  }
  
  # Band-pass correlation at zero lag and over lags, tested against phase-randomised
  # surrogates of the first record, and the 201-yr running correlation
  bands <- list(c(8, 16), c(16, 32), c(32, 64))
  win <- 201; h <- (win - 1) / 2
  
  bandpass_comparison <- function(A, B, lab) {
    tab <- do.call(rbind, lapply(bands, function(bd) {
      a <- bandpass(A, bd[1], bd[2]); b <- bandpass(B, bd[1], bd[2])
      r <- cor(a[, 2], b[, 2])
      cc <- ccf(a[, 2], b[, 2], lag.max = round(bd[2] / 2), plot = FALSE)
      i <- which.max(abs(cc$acf))
      data.frame(Pair = lab, Band = sprintf("%d-%d", bd[1], bd[2]),
                 r0 = round(r, 3), ccf_max = round(cc$acf[i], 3),
                 lag_at_max = cc$lag[i])
    }))
    a <- bandpass(A, 32, 64); b <- bandpass(B, 32, 64)
    robs <- cor(a[, 2], b[, 2])
    rs <- replicate_progress(n_surrogates_bandpass, paste("Band-pass correlation 32-64 yr,", lab), {
      z <- data.frame(x = A$x, y = phase_randomise(A$y))
      cor(bandpass(z, 32, 64)[, 2], b[, 2])
    })
    pv <- mean(abs(rs) >= abs(robs))
    n <- nrow(a); rr <- rep(NA_real_, n)
    for (i in seq_len(n)) {
      lo <- i - h; hi <- i + h
      if (lo >= 1 && hi <= n) rr[i] <- cor(a[lo:hi, 2], b[lo:hi, 2])
    }
    cat(sprintf("%s\n", lab))
    print(tab, row.names = FALSE)
    cat(sprintf("  32-64 yr: r = %+0.3f, two-sided p = %.4f (5000 surrogates); running %d-yr median %+0.2f, same sign in %.0f%% of windows\n\n",
                robs, pv, win, median(rr, na.rm = TRUE),
                100 * mean(sign(rr) == sign(robs), na.rm = TRUE)))
    tab$p_32_64 <- c(NA, NA, round(pv, 3))
    list(a = a, b = b, r = robs, p = pv, run = rr, tab = tab)
  }
  
  max_abs_ccf <- function(a, g) max(abs(ccf(a, g, lag.max = 20, plot = FALSE)$acf))
  
  # Panels for Figure S4
  plot_sawtooth_panels <- function(S, K, lab1, col1, ttl, show_stats = FALSE) {
    # colour key and title
    par(mar = c(0.4, 5.0, 1.6, 5.0))
    draw_colour_key(K$brk, lab = "Cross-wavelet power", frac = 0.45,
                    bar_h = 0.30, y0 = 0.46)
    mtext(ttl, side = 3, line = 0.2, adj = 0, cex = 0.85, font = 2)
    # scalogram
    par(mar = c(2.8, 5.0, 0.6, 5.0))
    plot_cross_wavelet(K)
    mtext("A", side = 3, line = 0.1, adj = -0.09, font = 2, cex = 0.9)
    # band-passed records
    par(mar = c(2.8, 5.0, 2.0, 5.0))
    plot(S$a[, 1], S$a[, 2], type = "l", lwd = 1.7, col = col1, xaxt = "n",
         xlab = "", ylab = paste0(lab1, " (32-64 yr)"), xaxs = "i",
         cex.lab = 0.9)
    add_bce_ce_axis(250)
    mtext("B", side = 3, line = 0.4, adj = -0.09, font = 2, cex = 0.9)
    par(new = TRUE)
    plot(S$b[, 1], S$b[, 2], type = "l", lwd = 1.7, col = "#B2182B",
         axes = FALSE, xlab = "", ylab = "", xaxs = "i")
    axis(4, col = "#B2182B", col.axis = "#B2182B")
    mtext("Angstel/Vecht GI (32-64 yr)", side = 4, line = 2.7,
          col = "#B2182B", cex = 0.72)
    if (show_stats)
      mtext(sprintf("r = %+0.2f, p = %.3f", S$r, S$p), side = 3, line = 0.3,
            adj = 1, cex = 0.72, col = "grey25")
    # running correlation
    par(mar = c(2.8, 5.0, 2.0, 5.0))
    plot(S$a[, 1], S$run, type = "l", lwd = 1.7, ylim = c(-1, 1), xaxt = "n",
         xlab = "", ylab = "201-yr running correlation", xaxs = "i",
         cex.lab = 0.9)
    abline(h = 0, lty = 2)
    add_bce_ce_axis(250)
    mtext("C", side = 3, line = 0.4, adj = -0.09, font = 2, cex = 0.9)
    mtext("Year (BCE/CE)", side = 1, line = 2.2, cex = 0.75)
  }
}


## Change in spectral organisation ####
# Log ratio of power at 48-80 yr and 14-17 yr, smoothed over 31 yr
{
  band_power_ratio <- function(x, yr, smooth = 31) {
    w <- morlet_wavelet(x, dj = dj_test)
    b1 <- w$period >= 48 & w$period <= 80
    b2 <- w$period >= 14 & w$period <= 17
    p1 <- colMeans(ifelse(w$valid[b1, , drop = FALSE], w$power[b1, , drop = FALSE], NA), na.rm = TRUE)
    p2 <- colMeans(ifelse(w$valid[b2, , drop = FALSE], w$power[b2, , drop = FALSE], NA), na.rm = TRUE)
    lr <- log(p1 / p2)
    k  <- min(smooth, length(lr))
    sm <- stats::filter(lr, rep(1 / k, k), sides = 2)
    as.numeric(sm)
  }
  
  # Largest step in the ratio, with at least 200 valid years on both sides
  largest_step <- function(sm, yr, edge = 200) {
    n <- length(sm)
    ok <- which(!is.na(sm))
    idx <- ok[ok > min(ok) + edge & ok < max(ok) - edge]
    d <- vapply(idx, function(i)
      mean(sm[(i + 1):max(ok)], na.rm = TRUE) - mean(sm[min(ok):i], na.rm = TRUE), 0)
    list(year = yr[idx[which.max(abs(d))]], stat = max(abs(d)),
         idx = idx, d = d)
  }
  
  # Plotting year back to historical year (no year 0)
  historical_year <- function(py) ifelse(py <= 0, py - 1, py)
}


# Load data ####

## Angstel/Vecht growth index and diagnostics ####
{
  growth <- read_chronology()
  growth <- growth[growth$Year >= year_min & growth$Year <= year_max, ]
  growth$PlotYear <- ifelse(growth$Year < 0, growth$Year + 1, growth$Year)
  
  cat("\nProxy:", chronology,
      "\nYears:", nrow(growth),
      "(", min(growth$Year), "to", max(growth$Year), ")",
      "\nSample depth:", min(growth$samp_depth), "-", max(growth$samp_depth),
      " median", median(growth$samp_depth), "\n")
  
  # Growth index on the plotting time axis
  VA25mean <- data.frame(x = growth$PlotYear, y = growth$GI)
  VA25mean <- linterp(VA25mean, genplot = FALSE, verbose = FALSE)
  VA25mean <- sortNave(VA25mean, genplot = FALSE, verbose = FALSE)
  
  # Sample depth, running Rbar, EPS and running variance (Figure S1)
  # EPS is computed from sample depth and the 101 yr running Rbar (Wigley et al., 1984)
  DG <- read_diagnostics()
  DG <- DG[DG$Year >= year_min & DG$Year <= year_max, ]
  DG$PlotYear <- ifelse(DG$Year < 0, DG$Year + 1, DG$Year)
  DG$EPS <- DG$samp_depth * DG$Rbar101 / (1 + (DG$samp_depth - 1) * DG$Rbar101)
  DG$var_std <- running_variance(DG$std050)
  DG$var_vsc <- running_variance(DG$vsc101)
}

## Comparison records on the common window 889 BCE - 156 CE ####
{
  stal_url <- "https://raw.githubusercontent.com/stratigraphy/Weesp-dendro-cycles/main/stal_scot.csv"
  saw_url  <- paste0("https://www.ncei.noaa.gov/pub/data/paleo/paleolimnology/",
                     "northamerica/canada/ellesmere/sawtooth2020-noaa.txt")
  
  have_stal <- download_if_missing(stal_url, "stal_scot.csv")
  have_saw  <- download_if_missing(saw_url,  "sawtooth2020.txt")
  if (!have_stal || !have_saw)
    stop("The comparison needs stal_scot.csv and sawtooth2020.txt; ",
         "download them from the URLs given in the header and rerun.")
  
  stal <- read.csv("stal_scot.csv")[, c(2, 3)]
  stal <- stal[stats::complete.cases(stal), ]
  stal <- linterp(sortNave(stal, genplot = FALSE, verbose = FALSE),
                  genplot = FALSE, verbose = FALSE)
  names(stal) <- c("x", "y")
  
  raw <- readLines("sawtooth2020.txt", warn = FALSE)
  raw <- raw[!grepl("^#", raw) & !grepl("^depth", raw)]
  sp  <- strsplit(trimws(raw), "[ \t]+"); sp <- sp[vapply(sp, length, 1L) >= 5]
  saw <- data.frame(x = as.numeric(vapply(sp, `[`, "", 2)),
                    y = as.numeric(vapply(sp, `[`, "", 5)))
  saw <- saw[saw$y > -999, ]; saw <- saw[order(saw$x), ]
  
  GIraw <- VA25mean; names(GIraw) <- c("x", "y")
  
  gmin <- max(min(stal$x), min(saw$x), min(GIraw$x))
  gmax <- min(max(stal$x), max(saw$x), max(GIraw$x))
  grid <- seq(gmin, gmax, by = 1)
  cat(sprintf("\nCommon window for all comparisons: %d to %d = %d years\n",
              gmin, gmax, length(grid)))
  on_grid <- function(d) data.frame(
    x = grid, y = as.numeric(scale(approx(d$x, d$y, grid)$y)))
  ST <- on_grid(stal); SW <- on_grid(saw); GI <- on_grid(GIraw)
}


# MAIN TEXT ####

# Figure 1: palaeogeographic map ####
# The map is drawn from Vos et al. (2018) and is not produced by the analysis;
# it is only copied to the output folder when it can be downloaded
# Output: pdf_figures/Figure1.pdf, png_figures/Figure1.png
{
  # palgeo.jpg is the map modified after Vos et al. (2018), stored in the data folder
  if (file.exists("palgeo.jpg") && requireNamespace("jpeg", quietly = TRUE)) {
    img <- jpeg::readJPEG("palgeo.jpg")
    
    plot_figure1 <- function() {
      par(mar = c(0, 0, 0, 0))
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
           xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
      rasterImage(img, 0, 0, 1, 1)
    }
    
    asp <- dim(img)[1] / dim(img)[2]
    pdf("pdf_figures/Figure1.pdf", width = 8, height = 8 * asp, onefile = FALSE)
    plot_figure1()
    dev.off()
    
    png("png_figures/Figure1.png", width = 2000, height = round(2000 * asp), res = 250)
    plot_figure1()
    dev.off()
  } else cat("Figure 1: palgeo.jpg or the jpeg package not found, figure skipped\n")
}


# Figure 2: growth index and sample depth ####
# Output: pdf_figures/Figure2.pdf, png_figures/Figure2.png
{
  plot_figure2 <- function() {
    
    par(mfrow = c(2, 1), mar = c(4, 4, 2, 4), oma = c(0, 0, 1, 0))
    
    plot(growth$PlotYear, growth$GI, type = "l", col = "grey25", lwd = 0.7, xaxt = "n",
         xaxs = "i", xlab = "", ylab = "Growth index",
         main = "Angstel/Vecht growth index (std050)", cex.main = 1)
    add_bce_ce_axis()
    
    plot(growth$PlotYear, growth$samp_depth, type = "l", col = "#2C7BB6", lwd = 1.2,
         xaxt = "n", xaxs = "i", xlab = "Year", ylab = "Sample depth",
         main = "Replication", cex.main = 1)
    abline(h = c(10, 40), lty = c(3, 2))
    add_bce_ce_axis()
    
  }
  
  pdf("pdf_figures/Figure2.pdf", width = 12, height = 7.2, onefile = FALSE)
  plot_figure2()
  dev.off()
  
  png("png_figures/Figure2.png", width = 3000, height = 1800, res = 250)
  plot_figure2()
  dev.off()
}


# Table 1: significance of the spectral peaks ####
# Maximum-in-band test against four backgrounds
# Output: tables/Table_1_cycle_significance.csv
{
  x <- VA25mean[, 2]
  n <- length(x)
  
  # Transform used for the statistics
  w_test <- morlet_wavelet(x, dj = dj_test)
  gws_obs <- average_spectral_power(w_test)
  
  backgrounds <- make_backgrounds(x, w_test$period, gws_obs, w_test$scales)
  
  cat("\nTable 1: surrogate significance\n")
  sig <- run_or_load(paste0("Table1_", chronology, "_n", n_surrogates),
                     lapply(backgrounds, surrogate_significance))
  
  cycle_table <- do.call(rbind, lapply(sig, function(z) z$cand))
  rownames(cycle_table) <- NULL
  write.csv(cycle_table,
            file.path("tables", "Table_1_cycle_significance.csv"), row.names = FALSE)
  
  print(reshape(cycle_table[, c("Candidate", "Background", "p_max_in_window")],
                idvar = "Candidate", timevar = "Background", direction = "wide"),
        row.names = FALSE)
}


# Figure 3: wavelet scalogram with significance ####
# Output: pdf_figures/Figure3.pdf, png_figures/Figure3.png
{
  plot_figure3 <- function() {
    
    layout(matrix(c(1, 4, 2, 3), nrow = 2, byrow = TRUE),
           widths = c(8, 2), heights = c(1, 3))
    
    # Growth index
    par(mar = c(0, 5, 2, 0))
    plot(growth$PlotYear, growth$GI, type = "l", xaxt = "n", xaxs = "i",
         xlab = "", ylab = "GI", main = "Angstel/Vecht std050")
    
    # Scalogram
    par(mar = c(4, 5, 0, 0))
    brk <- quantile_breaks(w_test$power)
    image(x = growth$PlotYear, y = log2(w_test$period), z = t(w_test$power),
          col = batlow_colours(100), breaks = brk, useRaster = TRUE,
          xaxs = "i", axes = FALSE,
          xlab = "Year (BCE/CE)", ylab = "Period (years)")
    add_bce_ce_axis()
    pt <- c(2, 4, 8, 16, 32, 64, 128, 256)
    axis(2, at = log2(pt), labels = pt, las = 2)
    box()
    
    # 95% contours for the two AR(1) backgrounds
    for (bn in c("ar1_observed", "ar1_corrected")) {
      mask <- (w_test$power > sig[[bn]]$crit_scale) & w_test$valid
      contour(x = growth$PlotYear, y = log2(w_test$period), z = t(mask * 1),
              levels = 0.5, add = TRUE, drawlabels = FALSE,
              lwd = if (bn == "ar1_observed") 1.6 else 1.1,
              lty = if (bn == "ar1_observed") 1 else 2,
              col = if (bn == "ar1_observed") "black" else "grey25")
    }
    
    # Cone of influence
    shade_coi(growth$PlotYear, w_test$coi, 2, max(w_test$period))
    
    # Average spectral power with the 95% level of each background
    par(mar = c(4, 0, 0, 1))
    plot(gws_obs, log2(w_test$period), type = "l", lwd = 2, yaxs = "i",
         yaxt = "n", xlab = "Avg. power", ylab = "")
    cols <- c(white = "grey50", ar1_observed = "firebrick",
              ar1_corrected = "darkorange", powerlaw = "steelblue")
    for (bn in names(cols))
      lines(sig[[bn]]$gws95, log2(w_test$period), col = cols[bn], lty = 2)
    abline(h = log2(candidate_periods), lty = 3, col = "grey80")
    
    par(mar = c(0, 0.5, 2, 1))
    plot.new()
    legend("top", bty = "n", cex = 0.78, lwd = c(2.2, 1.2, 1.2, 1.2, 1.2),
           lty = c(1, 2, 2, 2, 2), col = c("black", cols),
           legend = c("observed", "white noise 95%", expression(paste("AR(1) ", phi, " = 0.662, 95%")),
                      expression(paste("AR(1) ", phi, " = 0.510, 95%")), "power law 95%"))
    # Colour key
    par(new = TRUE, mar = c(0.5, 0.8, 9.0, 1.2))
    draw_colour_key(brk, lab = "Scalogram power", frac = 0.86, bar_h = 0.26,
                    y0 = 0.50, align = "centre", lab_pos = "below", n_lab = 3)
    
  }
  
  pdf("pdf_figures/Figure3.pdf", width = 12.8, height = 8, onefile = FALSE)
  plot_figure3()
  dev.off()
  
  png("png_figures/Figure3.png", width = 3200, height = 2000, res = 250)
  plot_figure3()
  dev.off()
}


# Figure 4: extracted cycles and their share of spectral power ####
# Cycles supported by at least three of the four backgrounds (Table 1), and the borderline ~24 yr cycle
# Output: pdf_figures/Figure4.pdf, png_figures/Figure4.png
{
  w_plot <- morlet_wavelet(x, dj = dj_figure)
  
  cyc <- list(
    "60"   = list(f = taner_cycle(60, 50, 67),       lo = 50,   hi = 67),
    "38"   = list(f = taner_cycle(38, 31.5, 44.5),   lo = 31.5, hi = 44.5),
    "24"   = list(f = taner_cycle(24, 20, 28),       lo = 20,   hi = 28),
    "15.4" = list(f = taner_cycle(15.4, 12, 19),     lo = 12,   hi = 19)
  )
  panels <- names(cyc)
  
  # Share of total wavelet power carried by each cycle (period window 0.8-1.2 x the cycle),
  # computed for every year, including cells inside the cone of influence
  xp    <- growth$PlotYear
  xlims <- range(xp)
  share <- sapply(panels, function(k) {
    cc <- as.numeric(k)
    rr <- which(w_plot$period >= cc * 0.8 & w_plot$period <= cc * 1.2)
    colSums(w_plot$power[rr, , drop = FALSE]) / colSums(w_plot$power)
  })
  
  cols4 <- setNames(
    if (requireNamespace("scico", quietly = TRUE))
      scico::scico(4, palette = "batlow", begin = 0.1, end = 0.85)
    else c("#6baed6", "#74c476", "#9e9ac8", "#fd8d3c"),
    panels)
  
  plot_figure4 <- function() {
    
    layout(matrix(1:5, nrow = 5), heights = c(1, 1, 1, 1, 3.4))
    par(mar = c(0, 5, 1.6, 1))
    
    for (i in seq_along(panels)) {
      nmv <- panels[i]
      ttl <- sprintf("~%s yr cycle (passband %.0f-%.0f yr)%s",
                     nmv, cyc[[nmv]]$lo, cyc[[nmv]]$hi,
                     n_backgrounds_passed(as.numeric(nmv)))
      plot(VA25mean, type = "l", col = "grey65", lwd = 0.6, xaxt = "n",
           xlab = "", ylab = "GI", main = ttl, cex.main = 0.95,
           xlim = xlims, xaxs = "i")
      lines(cyc[[nmv]]$f, col = cols4[nmv], lwd = 2.2)
      add_label(LETTERS[i])
    }
    
    par(mar = c(4, 5, 1.6, 1))
    plot(0, 0, type = "n", xlim = xlims, ylim = c(0, 100), xaxt = "n",
         xlab = "Year (BCE/CE)", ylab = "Percentage spectral power (%)",
         main = "Share of total wavelet power carried by each cycle", cex.main = 0.95,
         xaxs = "i", yaxs = "i")
    add_bce_ce_axis()
    
    for (i in seq_along(panels)) {
      keys <- panels[i:length(panels)]
      top  <- rowSums(share[, keys, drop = FALSE]) * 100
      polygon(c(xp, rev(xp)), c(top, rep(0, length(xp))),
              col = cols4[panels[i]], border = NA)
    }
    
    # Cone of influence of the shortest cycle, beyond which the shares are edge-affected
    coi_years <- range(xp[w_plot$coi >= 15.4 * 1.2])
    abline(v = coi_years, lty = 2, col = "grey20")
    box()
    
    add_label("E")
    legend("topright", legend = paste0("~", panels, " yr"), fill = cols4[panels],
           border = NA, bty = "n", cex = 0.9, title = "Cycles")
    
  }
  
  pdf("pdf_figures/Figure4.pdf", width = 12, height = 9.2, onefile = FALSE)
  plot_figure4()
  dev.off()
  
  png("png_figures/Figure4.png", width = 3000, height = 2300, res = 250)
  plot_figure4()
  dev.off()
}

# Section 4.3: persistence (area-wise test) ####
# Longest run above the 95% level per octave band, against the two AR(1) backgrounds
# Output: tables/Section_4.3_persistence_test.csv
{
  band_def <- lapply(seq_len(nrow(octave_bands)), function(b) {
    last <- b == nrow(octave_bands)
    rr <- if (last) which(w_test$period >= octave_bands$Period_min[b] &
                            w_test$period <= octave_bands$Period_max[b])
    else which(w_test$period >= octave_bands$Period_min[b] &
                 w_test$period < octave_bands$Period_max[b])
    if (!length(rr)) return(NULL)
    cc <- which(apply(w_test$valid[rr, , drop = FALSE], 2, all))
    if (!length(cc)) return(NULL)
    list(rows = rr, cols = cc,
         value = colMeans(w_test$power[rr, cc, drop = FALSE]))
  })
  
  areawise_table <- run_or_load(paste0("Persistence_", chronology, "_n", n_surrogate_fields), {
    areawise_table <- list(); ai <- 0
    
    for (bn in c("ar1_observed", "ar1_corrected")) {
      
      sim_bands <- vector("list", n_surrogate_fields)
      
      pb <- progress_start(n_surrogate_fields, paste("Persistence test, background:", backgrounds[[bn]]$label))
      for (s in seq_len(n_surrogate_fields)) {
        ws <- morlet_wavelet(backgrounds[[bn]]$sim(n, backgrounds[[bn]]$par),
                             dj = dj_test)
        sim_bands[[s]] <- lapply(band_def, function(bd)
          if (is.null(bd)) NULL
          else colMeans(ws$power[bd$rows, bd$cols, drop = FALSE]))
        setTxtProgressBar(pb, s)
      }
      close(pb)
      
      for (b in seq_len(nrow(octave_bands))) {
        bd <- band_def[[b]]
        if (is.null(bd)) next
        crit <- quantile(unlist(lapply(sim_bands, function(z) z[[b]])),
                         0.95, names = FALSE)
        obs_run <- longest_run_above(bd$value, crit)
        sim_runs <- vapply(sim_bands,
                           function(z) longest_run_above(z[[b]], crit), numeric(1))
        ai <- ai + 1
        areawise_table[[ai]] <- data.frame(
          Band = octave_bands$Band[b],
          Background = backgrounds[[bn]]$label,
          N_years_outside_COI = length(bd$cols),
          Percent_years_above = 100 * mean(bd$value > crit),
          Longest_run_years = obs_run,
          Null_p95_run = quantile(sim_runs, 0.95, names = FALSE),
          p_areawise = mean(sim_runs >= obs_run),
          stringsAsFactors = FALSE)
      }
    }
    
    do.call(rbind, areawise_table)
  })
  write.csv(areawise_table, file.path("tables", "Section_4.3_persistence_test.csv"),
            row.names = FALSE)
  print(areawise_table, row.names = FALSE)
}


# Figure 5: bicoherence ####
# Output: pdf_figures/Figure5.pdf, png_figures/Figure5.png
{
  bic <- bicoherence(VA25mean)
  
  cat("\nBicoherence settings: segments =", bic$segments,
      " overlap =", bic$overlap, "%  CL =", bic$CL,
      "%  conditioning factor = 75th percentile of the denominator\n")
  
  plot_figure5 <- function() {
    
    layout(matrix(c(2, 4, 1, 3), nrow = 2, byrow = TRUE),
           heights = c(1, 2), widths = c(2, 1))
    
    per_lab <- c(100, 60, 40, 24, 16, 12, 9.5, 8, 6)
    per_lines <- 1 / per_lab
    
    # --- Main panel ---
    par(mar = c(5, 5, 0, 0))
    image(bic$freq, bic$freq, bic$bic2, col = batlow_colours(100),
          breaks = quantile_breaks(bic$bic2), useRaster = TRUE,
          xlab = "Frequency (cycles/year)", ylab = "Frequency (cycles/year)")
    abline(v = per_lines, col = "#FFFFFF28", lwd = 4)
    abline(h = per_lines, col = "#FFFFFF28", lwd = 4)
    contour(bic$freq, bic$freq, bic$bic2, add = TRUE, levels = bic$CLset,
            drawlabels = FALSE, lwd = 0.9, col = "grey15")
    abline(a = 0, b = 1, lwd = 1.2, lty = 3, col = "grey20")
    
    # Sum-frequency interactions shown in Figure 6
    for (ij in list(c(24, 40), c(60, 100))) {
      x0 <- 1 / ij[1]; y0 <- 1 / ij[2]
      xe <- max(bic$freq) * 1.2
      segments(x0, y0, xe, y0 - (xe - x0), col = "grey10", lwd = 2, lty = 2)
      points(x0, y0, pch = 21, bg = "white", col = "grey10", cex = 1.2)
      points(y0, x0, pch = 21, bg = "white", col = "grey10", cex = 1.2)
    }
    box()
    
    # --- Top panel, power spectrum with period labels ---
    par(mar = c(0, 5, 2, 0))
    plot(bic$freq, bic$pwr, type = "l", lwd = 2, log = "y", xaxs = "i",
         xaxt = "n", xlab = "", ylab = "Power")
    abline(v = per_lines, col = "#80808038", lwd = 5)
    for (i in seq_along(per_lab))
      mtext(per_lab[i], side = 3, at = per_lines[i], line = 0.15, cex = 0.7)
    box()
    
    # --- Right panel, power spectrum rotated ---
    par(mar = c(5, 0, 0, 2))
    plot(bic$pwr, bic$freq, type = "l", lwd = 2, log = "x", yaxs = "i",
         yaxt = "n", xlab = "Power", ylab = "")
    abline(h = per_lines, col = "#80808038", lwd = 5)
    for (i in seq_along(per_lab))
      mtext(per_lab[i], side = 4, at = per_lines[i], line = 0.3, las = 1,
            cex = 0.7)
    box()
    
    # --- Colour key ---
    par(mar = c(2, 2, 3, 2))
    draw_colour_key(quantile_breaks(bic$bic2), lab = "Magnitude-squared bicoherence",
                    frac = 0.60, bar_h = 0.16, y0 = 0.55, align = "centre",
                    lab_pos = "below", n_lab = 4)
    
  }
  
  pdf("pdf_figures/Figure5.pdf", width = 11.2, height = 10.4, onefile = FALSE)
  plot_figure5()
  dev.off()
  
  png("png_figures/Figure5.png", width = 2800, height = 2600, res = 250)
  plot_figure5()
  dev.off()
}



# Figure 6: combination tones ####
# Output: pdf_figures/Figure6.pdf, png_figures/Figure6.png
{
  comb_37 <- combination_tone(60, 100)   # -> 37.5 yr
  comb_15 <- combination_tone(24, 40)    # -> 15.0 yr
  
  # All curves span the full interval; every panel and its second axis use the same x-limits
  xlims <- range(VA25mean[, 1])
  
  plot_figure6 <- function() {
    
    par(mfrow = c(3, 1), mar = c(4, 4, 2, 4))
    
    plot(VA25mean, type = "l", xaxt = "n", xlab = "", ylab = "GI", col = "grey65",
         lwd = 0.6, xlim = xlims, xaxs = "i", cex.main = 1,
         main = "A. Growth index with the 15.0- and 37.5-yr components")
    lines(comb_15$direct[, 1], comb_15$direct[, 2] + mean(x), col = "#2C7BB6", lwd = 1.6)
    lines(comb_37$direct[, 1], comb_37$direct[, 2] + mean(x), col = "#B2182B", lwd = 1.6)
    add_bce_ce_axis()
    
    plot(comb_37$direct[, 1], comb_37$direct[, 2] + mean(x), type = "l", lwd = 1.8,
         col = "#B2182B", xaxt = "n", xlab = "", ylab = "GI", xlim = xlims, xaxs = "i",
         cex.main = 1,
         main = sprintf("B. ~%.0f yr: extracted (red) vs synthesised from 60 and 100 yr (sum tone) (black), r = %+.2f",
                        comb_37$period, comb_37$r))
    add_bce_ce_axis()
    par(new = TRUE)
    plot(comb_37$synthetic, type = "l", axes = FALSE, xlab = "", ylab = "",
         xlim = xlims, xaxs = "i")
    axis(4); mtext("Nonlinear combination", side = 4, line = 2.5, cex = 0.8)
    
    plot(comb_15$direct[, 1], comb_15$direct[, 2] + mean(x), type = "l", lwd = 1.8,
         col = "#2C7BB6", xaxt = "n", xlab = "Year (BCE/CE)", ylab = "GI",
         xlim = xlims, xaxs = "i", cex.main = 1,
         main = sprintf("C. ~%.0f yr: extracted (blue) vs synthesised from 24 and 40 yr (black), r = %+.2f",
                        comb_15$period, comb_15$r))
    add_bce_ce_axis()
    par(new = TRUE)
    plot(comb_15$synthetic, type = "l", axes = FALSE, xlab = "", ylab = "",
         xlim = xlims, xaxs = "i")
    axis(4); mtext("Nonlinear combination", side = 4, line = 2.5, cex = 0.8)
    
  }
  
  pdf("pdf_figures/Figure6.pdf", width = 12, height = 8.8, onefile = FALSE)
  plot_figure6()
  dev.off()
  
  png("png_figures/Figure6.png", width = 3000, height = 2200, res = 250)
  plot_figure6()
  dev.off()
}


# Figure 7: cross-wavelet of Uamh an Tartair and the growth index ####
# Output: pdf_figures/Figure7.pdf, png_figures/Figure7.png
{
  K1 <- cross_wavelet_for_plot(ST, GI)
  
  plot_figure7 <- function() {
    layout(matrix(1:2, ncol = 1), heights = c(0.16, 1))
    par(oma = c(3.0, 0, 0.2, 0))
    par(mar = c(0.4, 5.0, 1.6, 2.0))
    draw_colour_key(K1$brk, lab = "Cross-wavelet power", frac = 0.34,
                    bar_h = 0.30, y0 = 0.46)
    mtext("Uamh an Tartair stalagmite growth rate x Angstel/Vecht GI, 889 BCE-156 CE",
          side = 3, line = 0.2, adj = 0, cex = 0.8, font = 2)
    par(mar = c(0.6, 5.0, 0.6, 2.0))
    plot_cross_wavelet(K1)
    mtext("Year (BCE/CE)", side = 1, outer = TRUE, line = 1.6, cex = 0.8)
  }
  
  pdf("pdf_figures/Figure7.pdf", width = 12, height = 6.8, onefile = FALSE)
  plot_figure7()
  dev.off()
  
  png("png_figures/Figure7.png", width = 3000, height = 1700, res = 250)
  plot_figure7()
  dev.off()
}


# SUPPLEMENTARY MATERIAL ####

# Figure S1: sample depth, running variance, EPS and 32-64 yr power (S1) ####
# Output: pdf_figures/FigureS1.pdf, png_figures/FigureS1.png
{
  w_s1 <- morlet_wavelet(DG$std050, dj = dj_test)
  rows <- which(w_s1$period >= 32 & w_s1$period <= 64)
  DG$bp3264 <- vapply(seq_len(nrow(DG)), function(j) {
    v <- w_s1$valid[rows, j]
    if (!any(v)) NA_real_ else mean(w_s1$power[rows[v], j])
  }, numeric(1))
  
  cat("\nFigure S1 diagnostics (1283 BCE - 156 CE)\n")
  cat(sprintf("  EPS >= 0.85 in %.1f%% of years, median EPS %.2f\n",
              100 * mean(DG$EPS >= 0.85, na.rm = TRUE), median(DG$EPS, na.rm = TRUE)))
  cat(sprintf("  r(sample depth, running variance)     = %+.2f\n",
              cor(DG$samp_depth, DG$var_std, use = "complete.obs")))
  cat(sprintf("  r(sample depth, running Rbar)         = %+.2f\n",
              cor(DG$samp_depth, DG$Rbar101, use = "complete.obs")))
  cat(sprintf("  r(sample depth, 32-64-yr band power)  = %+.2f\n",
              cor(DG$samp_depth, DG$bp3264, use = "complete.obs")))
  cat(sprintf("  r(running variance, 32-64-yr power)   = %+.2f\n",
              cor(DG$var_std, DG$bp3264, use = "complete.obs")))
  cat(sprintf("  r(running variance std050, vsc101)    = %+.2f\n",
              cor(DG$var_std, DG$var_vsc, use = "complete.obs")))
  
  plot_figures1 <- function() {
    par(mfrow = c(4, 1), mar = c(1.5, 4.5, 1.2, 1), oma = c(3, 0, 0, 0))
    plot(DG$PlotYear, DG$samp_depth, type = "l", lwd = 1.2, col = "#2C7BB6",
         xaxt = "n", xaxs = "i", xlab = "", ylab = "Sample depth")
    abline(h = c(10, 40), lty = 3); mtext("A", side = 3, adj = 0, font = 2, cex = 0.8)
    add_bce_ce_axis()
    yl <- range(c(DG$var_std, DG$var_vsc), na.rm = TRUE)
    plot(DG$PlotYear, DG$var_std, type = "l", lwd = 1.2, col = "black", ylim = yl,
         xaxt = "n", xaxs = "i", xlab = "", ylab = "Running variance (101 yr)")
    lines(DG$PlotYear, DG$var_vsc, lwd = 1.2, col = "#B2182B")
    legend("topright", c("std050", "vsc101"), col = c("black", "#B2182B"),
           lwd = 1.2, bty = "n", cex = 0.8)
    mtext("B", side = 3, adj = 0, font = 2, cex = 0.8)
    add_bce_ce_axis()
    plot(DG$PlotYear, DG$EPS, type = "l", lwd = 1.2, col = "grey25", ylim = c(0.6, 1),
         xaxt = "n", xaxs = "i", xlab = "", ylab = "EPS")
    abline(h = 0.85, lty = 3); mtext("C", side = 3, adj = 0, font = 2, cex = 0.8)
    add_bce_ce_axis()
    plot(DG$PlotYear, DG$bp3264, type = "l", lwd = 1.2, col = "#1B7837",
         xaxt = "n", xaxs = "i", xlab = "", ylab = "Power, 32-64 yr")
    mtext("D", side = 3, adj = 0, font = 2, cex = 0.8)
    add_bce_ce_axis()
    mtext("Year", side = 1, outer = TRUE, line = 1.5, cex = 0.8)
  }
  
  pdf("pdf_figures/FigureS1.pdf", width = 11.2, height = 9.6, onefile = FALSE)
  plot_figures1()
  dev.off()
  
  png("png_figures/FigureS1.png", width = 2800, height = 2400, res = 250)
  plot_figures1()
  dev.off()
}


# Table S5: multitaper cross-check (S2) ####
# Full analysis window, robust red noise (Mann and Lees, 1996) and power law;
# the harmonic F-test is listed but not used
# Output: tables/Table_S5_multitaper_crosscheck.csv
{
  GIfull <- data.frame(x = growth$PlotYear, y = growth$GI)
  GIfull <- sortNave(linterp(GIfull, genplot = FALSE, verbose = FALSE),
                     genplot = FALSE, verbose = FALSE)
  
  cand <- c(2.2, 3.6, 5.9, 8, 12, 15.4, 24, 38, 60, 80, 110)
  
  cat("--- MTM with a robust (Mann and Lees, 1996) red-noise background ---\n")
  ml <- mtmML96(GIfull, tbw = 3, padfac = 5, output = 1, genplot = FALSE,
                verbose = FALSE)
  cat("--- MTM with a fitted power-law background ---\n")
  pl <- mtmPL(GIfull, tbw = 3, padfac = 5, output = 1, genplot = FALSE,
              verbose = FALSE)
  
  cat("\ncolumns ML96:", paste(names(ml), collapse = " | "), "\n")
  cat("columns PL  :", paste(names(pl), collapse = " | "), "\n\n")
  
  # Robust red noise (AR1_CL) and power law (PowerLaw_CL); the harmonic F-test is not used
  tab <- do.call(rbind, lapply(cand, function(P) {
    fr <- ml$Frequency
    sel <- which(1 / fr >= P * 0.9 & 1 / fr <= P * 1.1)
    if (!length(sel)) return(NULL)
    i <- sel[which.max(ml$Power[sel])]
    data.frame(Candidate = P,
               MTM_peak  = round(1 / fr[i], 2),
               ML96_redCL = round(max(ml$AR1_CL[sel], na.rm = TRUE), 1),
               PL_CL      = round(max(pl$PowerLaw_CL[sel], na.rm = TRUE), 1),
               HarmF_CL   = round(max(ml$Harmonic_CL[sel], na.rm = TRUE), 1))
  }))
  print(tab, row.names = FALSE)
  write.csv(tab, file.path("tables", "Table_S5_multitaper_crosscheck.csv"), row.names = FALSE)
}


# Figure S2: the three chronologies over 1209-40 BCE (S5) ####
# Output: pdf_figures/FigureS2.pdf, png_figures/FigureS2.png
{
  core <- DG[DG$Year >= -1209 & DG$Year <= -40, ]
  chron <- c(std050 = "std050", vsc101 = "vsc101",
             constant40 = "constant40_median")
  gws3 <- lapply(chron, function(cn) {
    w <- morlet_wavelet(core[[cn]], dj = dj_figure)
    list(period = w$period, gws = average_spectral_power(w))
  })
  cat("\nFigure S2: correlation of CWT average spectral power, 1209-40 BCE\n")
  cat(sprintf("  std050 vs constant40 : %.4f\n", cor(gws3$std050$gws, gws3$constant40$gws, use = "complete.obs")))
  cat(sprintf("  std050 vs vsc101     : %.4f\n", cor(gws3$std050$gws, gws3$vsc101$gws, use = "complete.obs")))
  cat(sprintf("  vsc101 vs constant40 : %.4f\n", cor(gws3$vsc101$gws, gws3$constant40$gws, use = "complete.obs")))
  
  plot_figures2 <- function() {
    par(mar = c(4.5, 4.5, 1, 1))
    cols3 <- c("black", "#B2182B", "#2166AC")
    yl <- range(unlist(lapply(gws3, `[[`, "gws")), na.rm = TRUE)
    plot(gws3$std050$period, gws3$std050$gws, type = "n", log = "x", ylim = yl,
         xlab = "Period (years)", ylab = "CWT average spectral power",
         xaxt = "n")
    axis(1, at = c(2, 4, 8, 16, 32, 64, 128, 256))
    # constant40 is drawn last, dashed, because it lies on top of std050
    for (i in 1:3) lines(gws3[[i]]$period, gws3[[i]]$gws, col = cols3[i],
                         lwd = c(2.6, 1.4, 1.6)[i], lty = c(1, 1, 2)[i])
    legend("topright", c("std050 (original)", "vsc101 (variance-stabilised)",
                         "constant40 (40 series per year)"),
           col = cols3, lwd = c(2.6, 1.4, 1.6), lty = c(1, 1, 2), bty = "n", cex = 0.85)
  }
  
  pdf("pdf_figures/FigureS2.pdf", width = 11.2, height = 7.2, onefile = FALSE)
  plot_figures2()
  dev.off()
  
  png("png_figures/FigureS2.png", width = 2800, height = 1800, res = 250)
  plot_figures2()
  dev.off()
}


# Figure S3: band-passed stalagmite record and growth index (S7) ####
# Band-pass correlations for all three pairings, also used for Figure S4 and Table S14
# Output: pdf_figures/FigureS3.pdf, png_figures/FigureS3.png
{
  # Stalagmite vs Sawtooth (S3) is not plotted, it is reported in Table S14
  bp_all <- run_or_load(paste0("FigureS3_bandpass_", chronology, "_n", n_surrogates_bandpass), list(
    S1 = bandpass_comparison(ST, GI, "Uamh an Tartair stalagmite growth rate vs Angstel/Vecht GI"),
    S2 = bandpass_comparison(SW, GI, "Sawtooth Lake AMVI vs Angstel/Vecht GI"),
    S3 = bandpass_comparison(ST, SW, "Uamh an Tartair stalagmite growth rate vs Sawtooth Lake AMVI")))
  S1 <- bp_all$S1; S2 <- bp_all$S2; S3 <- bp_all$S3
  
  plot_figures3 <- function() {
    par(mfrow = c(2, 1), oma = c(1.6, 0, 0, 0))
    par(mar = c(2.8, 5.0, 2.0, 5.0))
    plot(S1$a[, 1], S1$a[, 2], type = "l", lwd = 1.7, col = "grey25", xaxt = "n",
         xlab = "", ylab = "Stalagmite growth rate (32-64 yr)", xaxs = "i", cex.lab = 0.9)
    add_bce_ce_axis(250)
    mtext("A", side = 3, line = 0.4, adj = -0.06, font = 2)
    par(new = TRUE)
    plot(S1$b[, 1], S1$b[, 2], type = "l", lwd = 1.7, col = "#B2182B",
         axes = FALSE, xlab = "", ylab = "", xaxs = "i")
    axis(4, col = "#B2182B", col.axis = "#B2182B")
    mtext("Angstel/Vecht GI (32-64 yr)", side = 4, line = 2.7, col = "#B2182B", cex = 0.72)
    plot(S1$a[, 1], S1$run, type = "l", lwd = 1.7, ylim = c(-1, 1), xaxt = "n",
         xlab = "", ylab = "201-yr running correlation", xaxs = "i", cex.lab = 0.9)
    abline(h = 0, lty = 2)
    add_bce_ce_axis(250)
    mtext("B", side = 3, line = 0.4, adj = -0.06, font = 2)
    mtext("Year (BCE/CE)", side = 1, outer = TRUE, line = 0.4, cex = 0.8)
  }
  
  pdf("pdf_figures/FigureS3.pdf", width = 12, height = 8.4, onefile = FALSE)
  plot_figures3()
  dev.off()
  
  png("png_figures/FigureS3.png", width = 3000, height = 2100, res = 250)
  plot_figures3()
  dev.off()
}


# Table S11: cross-wavelet phase by period (S7) ####
# Output: tables/Table_S11_phase_by_period.csv
{
  pST <- phase_by_period(ST, GI); pSW <- phase_by_period(SW, GI)
  tabS11 <- data.frame(Period = sprintf("%d-%d", head(edges, -1), edges[-1]),
                       ST_GI_phase = round(pST[, "phase"]), ST_GI_power = signif(pST[, "power"], 3),
                       ST_GI_phase_unrect = round(pST[, "phase_unrect"]),
                       ST_GI_power_unrect = signif(pST[, "power_unrect"], 3),
                       SW_GI_phase = round(pSW[, "phase"]), SW_GI_power = signif(pSW[, "power"], 3),
                       SW_GI_phase_unrect = round(pSW[, "phase_unrect"]),
                       SW_GI_power_unrect = signif(pSW[, "power_unrect"], 3))
  write.csv(tabS11, file.path("tables", "Table_S11_phase_by_period.csv"), row.names = FALSE)
  print(tabS11, row.names = FALSE)
}


# Table S12: band-pass correlation in two sub-bands (S7) ####
# The split at 48 yr was made after looking at the phase, so these p-values are descriptive
# Output: tables/Table_S12_subband_correlation.csv
{
  tabS12 <- run_or_load(paste0("TableS12_", chronology, "_n", n_surrogates_bandpass),
                        do.call(rbind, lapply(list(list("Stalagmite", ST), list("Sawtooth", SW)), function(pp)
                          do.call(rbind, lapply(list(c(32, 48), c(48, 80)), function(b) {
                            g  <- bandpass(GI, b[1], b[2])[, 2]
                            a  <- bandpass(pp[[2]], b[1], b[2])[, 2]
                            cc <- ccf(a, g, lag.max = 40, plot = FALSE)
                            r0 <- cor(a, g); m0 <- max_abs_ccf(a, g)
                            nul <- replicate_progress(n_surrogates_bandpass,
                                                      sprintf("Table S12, %s, %d-%d yr", pp[[1]], b[1], b[2]), {
                                                        as <- bandpass(data.frame(x = pp[[2]]$x, y = phase_randomise(pp[[2]]$y)), b[1], b[2])[, 2]
                                                        c(cor(as, g), max_abs_ccf(as, g)) })
                            data.frame(Pair = pp[[1]], Band = sprintf("%d-%d", b[1], b[2]),
                                       r0 = round(r0, 3), p_r0 = mean(abs(nul[1, ]) >= abs(r0)),
                                       min_r = round(min(cc$acf), 3), lag_min = cc$lag[which.min(cc$acf)],
                                       max_r = round(max(cc$acf), 3), lag_max = cc$lag[which.max(cc$acf)],
                                       maxabs_20 = round(m0, 3), p_maxabs = mean(nul[2, ] >= m0))
                          })))))
  write.csv(tabS12, file.path("tables", "Table_S12_subband_correlation.csv"), row.names = FALSE)
  print(tabS12, row.names = FALSE)
}


# S7: single delay fitted to the phase and lagged correlation ####
# A constant delay D rotates the phase with period as 360 D / P, fitted over 28-70 yr;
# the largest lagged correlation (32-80 yr) is tested with the lag search included
# Output: tables/S7_fitted_delay.csv
{
  Xl <- cross_wavelet(cbind(ST$x, ST$y), cbind(GI$x, GI$y), dj = dj_figure,
                      lowerPeriod = 1.5, upperPeriod = 192, omega_nr = 15)
  Pl <- matrix(Xl$Period, Xl$nr, Xl$nc); ok <- Xl$valid & Pl >= 28 & Pl <= 70
  wl <- t(Xl$Power)[ok]; phl <- t(Xl$Phase)[ok]; pel <- Pl[ok]
  fitD <- sapply(0:30, function(D) sum(wl * cos(phl - 2 * pi * D / pel)) / sum(wl))
  fitC <- max(sapply(seq(-pi, pi, length.out = 361),
                     function(a) sum(wl * cos(phl - a)) / sum(wl)))
  cat(sprintf("Fitted delay: stalagmite leads by %d yr (mean cos %.2f; best constant phase %.2f)\n",
              (0:30)[which.max(fitD)], max(fitD), fitC))
  a80 <- bandpass(ST, 32, 80)[, 2]; g80 <- bandpass(GI, 32, 80)[, 2]
  cc80 <- ccf(a80, g80, lag.max = 30, plot = FALSE)
  m80 <- max(cc80$acf)
  n80 <- run_or_load(paste0("S7_lag_", chronology, "_n", n_surrogates_bandpass),
                     replicate_progress(n_surrogates_bandpass, "S7, lagged correlation 32-80 yr", {
                       as <- bandpass(data.frame(x = ST$x, y = phase_randomise(ST$y)), 32, 80)[, 2]
                       max(ccf(as, g80, lag.max = 30, plot = FALSE)$acf) }))
  cat(sprintf("32-80 yr: largest lagged r = %+.2f at lag %d (negative = stalagmite leads), p = %.3f\n",
              m80, cc80$lag[which.max(cc80$acf)], mean(n80 >= m80)))
  write.csv(data.frame(fitted_delay_yr = (0:30)[which.max(fitD)], mean_cos_delay = round(max(fitD), 2),
                       mean_cos_constant_phase = round(fitC, 2), r_lagged_32_80 = round(m80, 2),
                       lag_yr = cc80$lag[which.max(cc80$acf)], p = mean(n80 >= m80)),
            file.path("tables", "S7_fitted_delay.csv"), row.names = FALSE)
}


# Figure S4: comparison with the Sawtooth Lake AMV index (S8) ####
# Output: pdf_figures/FigureS4.pdf, png_figures/FigureS4.png
{
  K2 <- cross_wavelet_for_plot(SW, GI)
  
  plot_figures4 <- function() {
    layout(matrix(1:4, ncol = 1), heights = c(0.17, 1, 1, 1))
    par(oma = c(0.6, 0, 0.4, 0))
    plot_sawtooth_panels(S2, K2, "Sawtooth Lake AMVI", "#1B4F72",
                         "Sawtooth Lake, Ellesmere Island, 889 BCE-156 CE", show_stats = TRUE)
  }
  
  pdf("pdf_figures/FigureS4.pdf", width = 8, height = 11.6, onefile = FALSE)
  plot_figures4()
  dev.off()
  
  png("png_figures/FigureS4.png", width = 2000, height = 2900, res = 250)
  plot_figures4()
  dev.off()
}


# Table S14: band-pass correlations for the three pairings (S8) ####
# Output: tables/Table_S14_bandpass_correlations.csv
{
  write.csv(rbind(S1$tab, S2$tab, S3$tab),
            file.path("tables", "Table_S14_bandpass_correlations.csv"),
            row.names = FALSE)
  print(rbind(S1$tab, S2$tab, S3$tab), row.names = FALSE)
}


# Figure S5: change in spectral organisation near 572 BCE (S9) ####
# Full analysis window
# Output: pdf_figures/FigureS5.pdf, png_figures/FigureS5.png, tables/FigS5_band_power_ratio.csv, tables/FigS5_step_test.csv
{
  x  <- growth$GI
  yr <- growth$PlotYear
  
  # Ratio used for the step test (cells outside the cone of influence)
  sm  <- band_power_ratio(x, yr)
  obs <- largest_step(sm, yr)
  obs_hist <- historical_year(obs$year)
  cat(sprintf("\nObserved strongest transition in the 48-80 / 14-17 yr power ratio: %d BCE (plot-axis %d), step = %.3f\n",
              -round(obs_hist), round(obs$year), obs$stat))
  
  # Years where the ratio is defined (outside the cone of influence); the figure shows only these
  xlim_s5 <- range(yr[is.finite(sm)])
  
  # AR(1) surrogates
  phi <- as.numeric(arima(x, order = c(1, 0, 0), method = "ML")$coef["ar1"])
  cat(sprintf("AR(1) phi of the std050 chronology: %.3f\n", phi))
  nsim <- n_surrogates_step
  step_null <- run_or_load(paste0("FigureS5_", chronology, "_n", nsim), {
    null_stat <- numeric(nsim); null_year <- numeric(nsim)
    pb <- progress_start(nsim, "Figure S5, AR(1) surrogates")
    for (i in seq_len(nsim)) {
      xs <- as.numeric(arima.sim(list(ar = phi), n = length(x),
                                 sd = sd(x) * sqrt(1 - phi^2))) + mean(x)
      s  <- largest_step(band_power_ratio(xs, yr), yr)
      null_stat[i] <- s$stat; null_year[i] <- s$year
      setTxtProgressBar(pb, i)
    }
    close(pb)
    list(stat = null_stat, year = null_year)
  })
  null_stat <- step_null$stat; null_year <- step_null$year
  p <- mean(null_stat >= obs$stat)
  cat(sprintf("AR(1) surrogates (n = %d): median step %.3f, 95th pct %.3f, p = %.3f\n",
              nsim, median(null_stat), quantile(null_stat, 0.95), p))
  cat(sprintf("Surrogate transition years spread over %d to %d BCE (sd %.0f yr)\n",
              -round(historical_year(min(null_year))), -round(historical_year(max(null_year))), sd(null_year)))
  
  write.csv(data.frame(Year = yr, log_ratio_48_80_over_14_17 = sm),
            file.path("tables", "FigS5_band_power_ratio.csv"), row.names = FALSE)
  write.csv(data.frame(observed_year_BCE = -round(obs_hist), observed_step = round(obs$stat, 3),
                       phi = round(phi, 3), n_sim = nsim,
                       null_median = round(median(null_stat), 3),
                       null_p95 = round(as.numeric(quantile(null_stat, 0.95)), 3),
                       p_value = round(p, 3),
                       null_year_min_BCE = -round(historical_year(min(null_year))),
                       null_year_max_BCE = -round(historical_year(max(null_year))),
                       null_year_sd = round(sd(null_year))),
            file.path("tables", "FigS5_step_test.csv"), row.names = FALSE)
  
  plot_figures5 <- function() {
    par(mar = c(4.5, 4.5, 2.5, 1))
    plot(yr, sm, type = "l", lwd = 2, col = "grey20", xaxt = "n",
         xlim = xlim_s5, xaxs = "i",
         xlab = "Year (BCE/CE)", ylab = "log( power 48-80 yr / power 14-17 yr )",
         main = sprintf("31-yr smoothed band-power ratio (step at %d BCE, p = %.2f against AR(1))",
                        -round(obs_hist), p))
    add_bce_ce_axis(250)
    abline(h = 0, lty = 3)
    abline(v = obs$year, col = "#B2182B", lwd = 2)
    legend("topleft", bty = "n", cex = 0.85,
           legend = sprintf("strongest step (%d BCE)", -round(obs_hist)),
           col = "#B2182B", lwd = 2)
  }
  
  pdf("pdf_figures/FigureS5.pdf", width = 10.4, height = 6, onefile = FALSE)
  plot_figures5()
  dev.off()
  
  png("png_figures/FigureS5.png", width = 2600, height = 1500, res = 250)
  plot_figures5()
  dev.off()
}


# Table S15: sensitivity of the AR(1) coefficient (S10) ####
# Full analysis window (x from Figure S5)
# Output: tables/Table_S15_ar1_sensitivity.csv
{
  cat(sprintf("n = %d\n", length(x)))
  cat(sprintf("%-58s phi = %.3f\n", "observed series (no notching)", ar1_coefficient(x)))
  cat(sprintf("%-58s phi = %.3f\n", "lag-1 sample autocorrelation (for comparison)", acf(x, plot=FALSE, lag.max=1)$acf[2]))
  sets <- list(
    "manuscript notch list: 3.6-110 yr, 10 bands"        = c(3.6,5.9,8,12,15.4,24,38,60,80,110),
    "  + 2.2 yr (near-Nyquist) added"                    = c(2.2,3.6,5.9,8,12,15.4,24,38,60,80,110),
    "  + 200 yr added"                                   = c(3.6,5.9,8,12,15.4,24,38,60,80,110,200),
    "only the four supported components (15.4,24,38,60)" = c(15.4,24,38,60),
    "only the multidecadal band (38, 60)"                = c(38,60),
    "only decadal and shorter (3.6-24)"                  = c(3.6,5.9,8,12,15.4,24))
  for (nm in names(sets)) cat(sprintf("%-58s phi = %.3f\n", nm, ar1_coefficient(remove_bands(x, sets[[nm]]))))
  cat("\nNote: adding the 2.2-year band RAISES phi. Removing variance close to\n")
  cat("the Nyquist period lengthens the effective decorrelation time, so the\n")
  cat("near-Nyquist band is left in the series when the background is fitted.\n")
  cat("\nSensitivity to the +/- tolerance around each notched period\n")
  for (tol in c(0.05,0.075,0.10,0.15,0.20))
    cat(sprintf("  tol = +/-%4.1f%%  (fraction of spectrum removed ~%.0f%%)   phi = %.3f\n",
                100*tol, 100*min(1, 2*tol*10), ar1_coefficient(remove_bands(x, c(3.6,5.9,8,12,15.4,24,38,60,80,110), tol = tol))))
  cat("\nphi falls as the notch widens, because more of the series is removed.\n")
  cat("The second AR(1) background is therefore a bound on the background, not\n")
  cat("an estimate of it; the two AR(1) variants together bracket the answer.\n")
  
  tols <- c(0.05, 0.075, 0.10, 0.15, 0.20)
  tabS15 <- data.frame(
    Estimated_on = c("observed series, no notching", names(sets),
                     sprintf("notch 3.6-110 yr, +/- %g%%", 100 * tols)),
    phi = round(c(ar1_coefficient(x),
                  sapply(sets, function(b) ar1_coefficient(remove_bands(x, b))),
                  sapply(tols, function(t) ar1_coefficient(remove_bands(x, notch_periods, tol = t)))), 3))
  write.csv(tabS15, file.path("tables", "Table_S15_ar1_sensitivity.csv"), row.names = FALSE)
}


# Session information ####
# Output: tables/sessionInfo_Weesp.txt
{
  capture.output(sessionInfo(),
                 file = file.path("tables", "sessionInfo_Weesp.txt"))
  
  cat("\nComplete. Figures are in pdf_figures and png_figures, tables in tables\n")
  
  if (n_surrogates < 10000 || n_surrogate_fields < 2000 || n_surrogates_bandpass < 20000) {
    cat("\n")
    cat("  DRAFT SURROGATE COUNTS WERE USED\n")
    cat(sprintf("    n_surrogates       = %6d   (publication: >= 10000)\n", n_surrogates))
    cat(sprintf("    n_surrogate_fields = %6d   (publication: >=  2000)\n", n_surrogate_fields))
    cat(sprintf("    n_surrogates_bandpass    = %6d   (publication: >= 20000)\n", n_surrogates_bandpass))
    cat("  The verdicts are the same at these counts, but the 95% lines of\n")
    cat("  Figure 3 carry visible Monte Carlo jitter and the p-values move\n")
    cat("  in the third decimal. Raise the counts for the final run and\n")
    cat("  refresh the p-values quoted in the manuscript from the output.\n")
  }
}


