# Angstel/Vecht bog-oak chronology (1283 BCE - 156 CE) ####
# Analysis code for Arts and Bazelmans, The Holocene
# This one script produces every figure, table and p-value of the paper and its supplement.
# Open Weesp.Rproj in RStudio, then source this file (or run it section by section).
# Input data are read from data/; results are written to output/:
#   output/figures_pdf/  vector PDF (Calibri)       output/figures_png/  PNG, 300 dpi
#   output/tables/       CSV tables                 output/saved_results/  completed surrogate runs
# The sections follow the order of the paper:
#   Main text: Figure 1, Figure 2, Table 1, Figure 3, Figure 4, persistence test (Section 4.3),
#              Figure 5, Figure 6, Figure 7
#   Supplement: Figure S1, Table S5, Figure S2, Figure S3, Table S9, Table S10, Table S11,
#               Table S12, fitted delay (S7), Table S13, Figure S4, Table S14, Figure S5,
#               Table S15, Tables S1-S4 and S6-S8 (supplementary suite at omega0 = 6),
#               detrending variants of S6 (only when their input file is present)
# In RStudio each figure and table is a section in the document outline (Ctrl+Shift+O)
# Each section is wrapped in { } so it can be folded, and run as one block by placing
# the cursor on the opening { and pressing Ctrl+Enter. Run the settings, functions and data first.

# Files needed in data/:
#   10_Angstel_Vecht_driechronos_std050_vsc101_constant040_16082026.xlsx
#   (or AV_diagnostic_data.csv and AV_constant40.csv), the growth index and its diagnostics
#   stal_scot.csv, Uamh an Tartair growth rate (Baker et al., 2015)
#   sawtooth2020.txt, Sawtooth Lake AMV index (Lapointe et al., 2020)
#   palgeo.jpg, the map of Figure 1
#   AV_objects_07.rds (optional), J. Bazelmans' detrending variants for S6
# The two comparison records are downloaded if they are missing.
# Completed surrogate runs are saved in output/saved_results/ and are not rerun (see settings).


# Packages and settings ####
{
  library(readxl)
  library(astrochron)     # taner(), linterp(), sortNave(), mtmML96(), mtmPL()
  library(multitaper)     # spec.mtm(), multitaper spectra of the supplementary suite
  library(systemfonts)    # glyph widths of Calibri for the PDF figures
  have_multitaper <- TRUE

  set.seed(1337)

  # Paths are relative to the project folder (the folder of Weesp.Rproj and this script).
  # The working directory is not changed: files are addressed through in_data() and in_out()
  data_dir <- "data"
  out_dir <- "output"
  if (!dir.exists(data_dir)) {
    stop("Folder data/ not found. Open Weesp.Rproj, or set the working directory to the ",
        "folder that holds this script, and run it again.")
  }
  in_data <- function(...) file.path(data_dir, ...)
  in_out <- function(...) file.path(out_dir, ...)
  for (d in c("figures_pdf", "figures_png", "tables", "saved_results")) {
    if (!dir.exists(in_out(d))) {
      dir.create(in_out(d), recursive = TRUE)
    }
  }

  data_file  <- in_data("10_Angstel_Vecht_driechronos_std050_vsc101_constant040_16082026.xlsx")
  chronology <- "std050"        # "std050" | "vsc101" | "constant40_median"
  year_min   <- -1283           # analysis window
  year_max   <-  156

  # Text size of the figures: the figures are drawn at their full width (8-13 in), so the
  # point size is raised from R's default of 12 to keep text legible at print size.
  # Margins keep their size in inches (see record_figure)
  figure_pointsize <- 15
  figure_font <- "Calibri"

  # Surrogate counts used for the published results (draft runs used 1000, 500 and 5000)
  n_surrogates          <- 10000   # Table 1, Figure 3
  n_surrogate_fields    <- 2000    # persistence test (Section 4.3)
  n_surrogates_bandpass <- 20000   # band-pass correlations (Figures S3, S4, Tables S12, S14)
  n_surrogates_step     <- 500    # change near 572 BCE (Figure S5)
  n_surrogates_s13      <- 1000   # significance of the Sawtooth record itself (Table S13)
  n_surrogates_suite    <- 2000   # supplementary suite at omega0 = 6 (Tables S2-S4, S7, S8)
  n_fields_suite        <- 1000   # its duration-based persistence test (Table S6)
  n_surrogates_variants <- 5000   # detrending variants (S6), optional

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
    file <- in_out("saved_results", paste0(name, ".rds"))
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

  # Evaluate expr n times with a progress bar. A scalar result gives a vector of length n,
  # a vector result a matrix with one column per evaluation (the shapes that replicate gives)
  replicate_progress <- function(n, label, expr) {
    f <- eval.parent(substitute(function() expr))
    pb <- progress_start(n, label)
    res <- vector("list", n)
    for (k in seq_len(n)) {
      res[[k]] <- f()
      setTxtProgressBar(pb, k)
    }
    close(pb)
    if (length(res[[1]]) == 1L) {
      out <- unlist(res, use.names = FALSE)
    } else {
      out <- do.call(cbind, res)
    }
    out
  }

  # Evaluate expr with its own random number stream, and restore the main stream afterwards,
  # so that a section added to the script does not change the results of later sections
  with_seed <- function(seed, expr) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
    expr
  }
}

## Writing the figures ####
# Each figure is drawn once on a null device and recorded. The recording is written as a
# vector PDF (Calibri, editable in Illustrator) and as a 300 dpi PNG. No cairo devices are
# used. In the PDF a hyphen after a letter, digit, ) or ~ is set as a soft hyphen, because
# pdf() sets "-" as a minus sign; numbers keep their minus sign. Text is set at
# figure_pointsize, and plot titles are set in the regular (not bold) face.
{
  # Calibri is not one of the standard PDF fonts, so pdf() is given its glyph widths in an
  # AFM file written from the installed font. The PDF names the font Calibri and does not
  # embed it: Illustrator and PDF viewers use the installed Calibri, and the text stays
  # editable
  write_afm <- function(path, index, afm_file) {
    enc_file <- system.file("enc", "WinAnsi.enc", package = "grDevices")
    enc <- readLines(enc_file)
    enc <- sub("%.*$", "", enc)
    tok <- unlist(strsplit(enc, "[[:space:]]+"))
    tok <- tok[grepl("^/", tok)]
    glyph_names <- sub("^/", "", tok[-1])
    chars <- character(length(glyph_names))
    for (k in seq_along(glyph_names)) {
      code <- k - 1
      chars[k] <- if (code >= 32) { iconv(rawToChar(as.raw(code)), "CP1252", "UTF-8") } else { NA }
    }
    keep <- !is.na(chars) & glyph_names != ".notdef"
    glyph_names <- c(glyph_names[keep], "minus")
    chars <- c(chars[keep], "\u2212")
    chars[glyph_names == "hyphen"] <- "-"
    first <- !duplicated(glyph_names)
    glyph_names <- glyph_names[first]
    chars <- enc2utf8(chars[first])

    g <- systemfonts::glyph_info(chars, path = path, index = index, size = 1000, res = 72)
    info <- systemfonts::font_info(path = path, index = index)
    bb <- matrix(0, length(chars), 4)
    for (k in seq_along(chars)) {
      bb[k, ] <- round(unlist(g$bbox[k]))
    }
    cap <- bb[glyph_names == "H", 4]
    xh <- bb[glyph_names == "x", 4]
    lines <- c("StartFontMetrics 2.0",
        paste("FontName", info$name),
        paste("FullName", info$name),
        paste("FamilyName", info$family),
        paste("Weight", if (info$bold) { "Bold" } else { "Regular" }),
        paste("ItalicAngle", if (info$italic) { -11 } else { 0 }),
        "IsFixedPitch false",
        paste("FontBBox", min(bb[, 1]), min(bb[, 3]), max(bb[, 2]), max(bb[, 4])),
        "EncodingScheme FontSpecific",
        paste("CapHeight", cap),
        paste("XHeight", xh),
        paste("Ascender", max(bb[, 4])),
        paste("Descender", min(bb[, 3])),
        paste("StartCharMetrics", length(chars)),
        sprintf("C -1 ; WX %d ; N %s ; B %d %d %d %d ;", as.integer(round(g$x_advance)),
            glyph_names, as.integer(bb[, 1]), as.integer(bb[, 3]), as.integer(bb[, 2]),
            as.integer(bb[, 4])),
        "EndCharMetrics",
        "EndFontMetrics")
    writeLines(lines, afm_file)
    afm_file
  }

  # Calibri is registered once per R session; pdf() refuses to register a font that is in use
  if (!(figure_font %in% names(grDevices::pdfFonts()))) {
    afm <- character(4)
    styles <- data.frame(weight = c("normal", "bold", "normal", "bold"),
        italic = c(FALSE, FALSE, TRUE, TRUE), stringsAsFactors = FALSE)
    for (k in 1:4) {
      f <- systemfonts::match_fonts(figure_font, italic = styles$italic[k], weight = styles$weight[k])
      afm[k] <- write_afm(f$path, f$index, file.path(tempdir(), paste0(figure_font, "_", k, ".afm")))
    }
    pdf_font_list <- list(grDevices::Type1Font(figure_font, afm))
    names(pdf_font_list) <- figure_font
    do.call(grDevices::pdfFonts, pdf_font_list)
  }
  if (.Platform$OS.type == "windows") {
    font_list <- list(grDevices::windowsFont(figure_font))
    names(font_list) <- figure_font
    do.call(grDevices::windowsFonts, font_list)
  }
  pdf_font <- figure_font
  png_font <- figure_font

  record_figure <- function(plot_fun, width, height) {
    pdf(NULL, width = width, height = height, family = pdf_font, encoding = "WinAnsi",
        pointsize = figure_pointsize)
    dev.control("enable")
    # layout() and par(mfrow) reset mex, cex and font.main, so they are set again before
    # every new plot: mex keeps the margins at their size for 12-point text, cex = 1 stops
    # multi-panel layouts from shrinking the text, and titles are regular
    old_hook <- getHook("before.plot.new")
    on.exit(setHook("before.plot.new", old_hook, "replace"))
    setHook("before.plot.new", function() {
      par(mex = 12 / figure_pointsize, cex = 1, font.main = 1)
    })
    plot_fun()
    fig <- recordPlot()
    invisible(dev.off())
    fig
  }

  # Remove an old output file first; on Windows a file that is open in a viewer is locked
  fig_file <- function(dir, name, ext) {
    f <- in_out(dir, paste0(name, ext))
    if (file.exists(f) && !file.remove(f)) {
      stop(f, " is open in another program; close it and run this section again")
    }
    f
  }

  soft_hyphen <- function(v) {
    if (inherits(v, "NativeSymbolInfo")) {
      return(v)
    }
    if (is.character(v)) {
      return(gsub("(?<=[[:alnum:])~])-", "\u00ad", v, perl = TRUE))
    }
    if (is.list(v) || is.call(v) || is.expression(v)) {
      for (i in seq_along(v)) {
        if (!is.null(v[[i]]) && !is.symbol(v[[i]])) {
          v[[i]] <- soft_hyphen(v[[i]])
        }
      }
    }
    v
  }

  save_figure <- function(fig, name, width, height) {
    pdf(fig_file("figures_pdf", name, ".pdf"), width = width, height = height,
        onefile = FALSE, bg = "white", family = pdf_font, encoding = "WinAnsi",
        useDingbats = FALSE, pointsize = figure_pointsize)
    fig_pdf <- fig
    fig_pdf[[1]] <- soft_hyphen(fig[[1]])
    replayPlot(fig_pdf)
    invisible(dev.off())
    png(fig_file("figures_png", name, ".png"), width = width, height = height, units = "in",
        res = 300, family = png_font, bg = "white", pointsize = figure_pointsize)
    replayPlot(fig)
    invisible(dev.off())
  }
}

## Plotting: colours, colour key, cone of influence, BCE/CE axis ####
# Colours: scico batlowK with quantile colour breaks, as in WaverideR
{
  palette_name <- "batlowK"

  batlow_colours <- function(n = 100) {
    if (requireNamespace("scico", quietly = TRUE)) {
      scico::scico(n = n, palette = palette_name)
    } else {
      hcl.colors(n, "YlOrRd", rev = TRUE)
    }
  }

  quantile_breaks <- function(z, n = 100)
      quantile(z, probs = seq(0, 1, length.out = n + 1), na.rm = TRUE)

  # Colour key for the scalograms
  draw_colour_key <- function(breaks, cols = batlow_colours(100), lab = "Power",
      n_lab = 4, cex = 0.8, frac = 0.42,
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

  # BCE/CE axis. There is no year 0: 1 BCE (-1) is followed by 1 CE, so BCE years are plotted at year + 1
  add_bce_ce_axis <- function(by = 250) {
    labs <- seq(-1250, 250, by = by)
    labs <- labs[!(labs %in% c(-1, 0, 1))]
    at <- ifelse(labs < 0, labs + 1, labs)
    minor <- c(seq(-1300, -50, by = 50), seq(50, 200, by = 50))   # no tick at year 0
    axis(1, at = ifelse(minor < 0, minor + 1, minor), labels = FALSE, tck = -0.02)
    axis(1, at = at, labels = labs, tck = -0.05)
    # only year 1 (1 CE) is marked at the BCE/CE boundary; there is no tick at year 0
    axis(1, at = 1, labels = "1", font = 2, tck = -0.05)
  }

  # Panel label outside the upper left corner of the current panel
  add_label <- function(lab) {
    mtext(lab, side = 3, line = 0.2, at = grconvertX(0.01, "nfc", "user"), adj = 0,
        font = 2, cex = 1.2)
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

    d <- read.csv(in_data("AV_diagnostic_data.csv"))
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
      d <- read.csv(in_data("AV_diagnostic_data.csv"))
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

  average_spectral_power <- function(w) {
    out <- numeric(length(w$period))
    for (i in seq_along(w$period)) {
      v <- w$valid[i, ]
      out[i] <- if (!any(v)) { NA_real_ } else { mean(w$power[i, v]) }
    }
    out
  }
}


## Background (noise) models ####
{
  ar1_coefficient <- function(x)
      as.numeric(arima(x, order = c(1, 0, 0), include.mean = TRUE,
          method = "ML")$coef["ar1"])

  remove_bands <- function(x, periods, tol = candidate_window) {
    n <- length(x)
    mu <- mean(x)
    f <- fft(x - mu)
    freq <- (0:(n - 1)) / n
    drop <- rep(FALSE, n)
    for (p in periods) {
      lo <- 1 / (p * (1 + tol))
      hi <- 1 / (p * (1 - tol))
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
    for (p in notch_periods) {
      keep <- keep & !(per >= p * (1 - candidate_window) & per <= p * (1 + candidate_window))
    }

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
  # Row-wise quantile and column-wise maximum of a matrix, ignoring NA
  row_quantile <- function(M, prob) {
    out <- numeric(nrow(M))
    for (i in seq_len(nrow(M))) {
      out[i] <- quantile(M[i, ], probs = prob, na.rm = TRUE, names = FALSE)
    }
    out
  }
  column_max <- function(M) {
    out <- numeric(ncol(M))
    for (j in seq_len(ncol(M))) {
      out[j] <- max(M[, j], na.rm = TRUE)
    }
    out
  }

  surrogate_significance <- function(bg, nsim = n_surrogates) {

    nper <- length(w_test$period)
    crit <- matrix(NA_real_, nrow = nper, ncol = nsim)
    gws_sim <- matrix(NA_real_, nrow = nper, ncol = nsim)

    pb <- progress_start(nsim, paste("Table 1, background:", bg$label))
    for (s in seq_len(nsim)) {
      ws <- morlet_wavelet(bg$sim(n, bg$par), dj = dj_test)
      for (i in seq_len(nper)) {
        v <- ws$valid[i, ]
        crit[i, s] <- if (!any(v)) { NA_real_ } else { quantile(ws$power[i, v], 0.95, names = FALSE) }
      }
      gws_sim[, s] <- average_spectral_power(ws)
      setTxtProgressBar(pb, s)
    }
    close(pb)

    cand <- vector("list", length(candidate_periods))
    for (j in seq_along(candidate_periods)) {
      p <- candidate_periods[j]
      rr <- which(w_test$period >= p * (1 - candidate_window) &
          w_test$period <= p * (1 + candidate_window))
      if (!length(rr)) {
        next
      }
      obs_max <- max(gws_obs[rr], na.rm = TRUE)
      sim_max <- column_max(gws_sim[rr, , drop = FALSE])
      cand[[j]] <- data.frame(Candidate = p,
          Peak_period = w_test$period[rr][which.max(gws_obs[rr])],
          Peak_power = obs_max,
          p_max_in_window = mean(sim_max >= obs_max),
          Background = bg$label,
          stringsAsFactors = FALSE)
    }
    cand <- do.call(rbind, cand)

    list(crit_scale = rowMeans(crit, na.rm = TRUE),
        gws95 = row_quantile(gws_sim, 0.95),
        cand = cand)
  }

  # Longest run of years above the 95% level, for the persistence test
  longest_run_above <- function(v, crit) {
    r <- rle(v > crit)
    i <- which(r$values)
    if (!length(i)) {
      0L
    } else {
      max(r$lengths[i])
    }
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
      if (length(vb)) {
        band[j] <- mean(w$power[vb, j])
      }
      if (length(va)) {
        total[j] <- sum(w$power[va, j])
      }
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
    if (!exists("cycle_table")) {
      return("")
    }
    z <- cycle_table[abs(cycle_table$Candidate - centre) < 1e-6, ]
    if (!nrow(z)) {
      return("")
    }
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
    if (padpts %% 2 != 0) {
      padpts <- padpts + 1
    }
    df <- 1 / (padpts * dt)
    freq <- df * (seq_len(padpts) - 1)

    ft <- matrix(NA_complex_, nrow = padpts, ncol = segments)
    ihold <- 1
    for (i in seq_len(segments)) {
      dd <- d[ihold:(ihold + segpts - 1), ]
      if (demean) {
        dd[, 2] <- dd[, 2] - mean(dd[, 2])
      }
      if (detrend) {
        dd[, 2] <- residuals(lm(dd[, 2] ~ dd[, 1]))
      }
      if (taper) {
        dd[, 2] <- dd[, 2] *
            0.5 * (1 - cos(2 * pi * seq_len(segpts) / (segpts + 1)))
      }
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

    B <- B / segments
    f1f2 <- f1f2 / segments
    f3 <- f3 / segments
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
    f1 <- 1 / P1
    f2 <- 1 / P2
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

    nr <- length(wt1$period)
    nc <- length(grid)

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
    X <- K$X
    rows <- K$rows
    xs <- X$axis.1
    ys <- log2(X$Period[rows])
    image(xs, ys, K$z, col = batlow_colours(100), breaks = K$brk, useRaster = TRUE,
        axes = FALSE, xlab = "", ylab = ylab, cex.lab = 0.9,
        ylim = range(ys), xaxs = "i")
    shade_coi(xs, X$coi, K$lo, K$hi)
    at <- 2^(1:8)
    at <- at[at >= K$lo & at <= K$hi]
    axis(2, at = log2(at), labels = at, las = 1, cex.axis = 0.8)
    add_bce_ce_axis(250)
    box()
    # one arrow per cell of a time-period grid, at the strongest ridge
    Pw <- t(X$Power)
    nr <- X$nr
    up <- rbind(Pw[-1, , drop = FALSE], -Inf)
    dn <- rbind(Inf, Pw[-nr, , drop = FALSE])
    md <- Pw > up & Pw > dn
    md[1, ] <- FALSE
    md[nr, ] <- FALSE
    sel <- md & (Pw >= quantile(X$Power, K$q, na.rm = TRUE)) & X$valid
    Ph <- t(X$Phase)
    nt <- 42
    np <- 16
    tbin <- cut(seq_len(X$nc), breaks = nt, labels = FALSE)
    pbin <- cut(log2(X$Period), breaks = np, labels = FALSE)
    for (tb in seq_len(nt)) {
      for (pb in seq_len(np)) {
        cand <- which(sel & outer(pbin == pb, tbin == tb, "&"), arr.ind = TRUE)
        cand <- cand[cand[, 1] %in% rows, , drop = FALSE]
        if (!nrow(cand)) {
          next
        }
        k <- cand[which.max(Pw[cand]), , drop = FALSE]
        r <- k[1]
        cc <- k[2]
        dx <- diff(range(xs)) * 0.014 * cos(Ph[r, cc])
        dy <- diff(range(ys)) * 0.028 * sin(Ph[r, cc])
        arrows(xs[cc] - dx, log2(X$Period[r]) - dy,
            xs[cc] + dx, log2(X$Period[r]) + dy,
            length = 0.020, lwd = 1.0, col = "black")
      }
    }
  }

  # Phase by period: power-weighted mean outside the cone of influence,
  # positive when the comparison record leads
  edges <- c(28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 70)
  phase_by_period <- function(A, B) {
    X <- cross_wavelet(cbind(A$x, A$y), cbind(B$x, B$y), dj = dj_figure,
        lowerPeriod = 1.5, upperPeriod = 192, omega_nr = 15)
    Pw <- t(X$Power)
    Ph <- t(X$Phase)
    nb <- length(edges) - 1
    out <- matrix(NA_real_, nrow = nb, ncol = 5,
        dimnames = list(NULL, c("phase", "R", "power", "phase_unrect", "power_unrect")))
    for (k in seq_len(nb)) {
      rows <- which(X$Period >= edges[k] & X$Period < edges[k + 1])
      v <- X$valid[rows, ]
      w <- Pw[rows, ][v]
      wu <- (Pw * X$Period)[rows, ][v]    # unrectified
      z  <- sum(w  * exp(1i * Ph[rows, ][v])) / sum(w)
      zu <- sum(wu * exp(1i * Ph[rows, ][v])) / sum(wu)
      out[k, ] <- c(Arg(z) * 180 / pi, Mod(z), mean(w), Arg(zu) * 180 / pi, mean(wu))
    }
    out
  }
}


## Band-pass comparison of two records ####
{
  download_if_missing <- function(url, file) {
    if (!file.exists(file)) {
      try(download.file(url, file, quiet = TRUE), silent = TRUE)
    }
    file.exists(file)
  }

  bandpass <- function(dat, lo, hi)
      taner(dat, padfac = 2, flow = 1 / hi, fhigh = 1 / lo, roll = 10^3,
          demean = TRUE, detrend = FALSE, addmean = FALSE, output = 1,
          xmin = 0, xmax = 0.5, genplot = FALSE, check = TRUE, verbose = FALSE)

  # Phase randomisation, keeps the spectrum (and autocorrelation) of the series
  phase_randomise <- function(v) {
    n <- length(v)
    f <- fft(v)
    k <- 2:floor(n / 2)
    f[k] <- Mod(f[k]) * exp(1i * runif(length(k), 0, 2 * pi))
    f[n - k + 2] <- Conj(f[k])
    Re(fft(f, inverse = TRUE)) / n
  }

  # Band-pass correlation at zero lag and over lags, tested against phase-randomised
  # surrogates of the first record, and the 201-yr running correlation
  bands <- list(c(8, 16), c(16, 32), c(32, 64))
  win <- 201
  h <- (win - 1) / 2

  bandpass_comparison <- function(A, B, lab) {
    tab <- vector("list", length(bands))
    for (j in seq_along(bands)) {
      bd <- bands[[j]]
      a <- bandpass(A, bd[1], bd[2])
      b <- bandpass(B, bd[1], bd[2])
      r <- cor(a[, 2], b[, 2])
      cc <- ccf(a[, 2], b[, 2], lag.max = round(bd[2] / 2), plot = FALSE)
      i <- which.max(abs(cc$acf))
      tab[[j]] <- data.frame(Pair = lab, Band = sprintf("%d-%d", bd[1], bd[2]),
          r0 = round(r, 3), ccf_max = round(cc$acf[i], 3),
          lag_at_max = cc$lag[i])
    }
    tab <- do.call(rbind, tab)
    a <- bandpass(A, 32, 64)
    b <- bandpass(B, 32, 64)
    robs <- cor(a[, 2], b[, 2])
    rs <- replicate_progress(n_surrogates_bandpass, paste("Band-pass correlation 32-64 yr,", lab), {
      z <- data.frame(x = A$x, y = phase_randomise(A$y))
      cor(bandpass(z, 32, 64)[, 2], b[, 2])
    })
    pv <- mean(abs(rs) >= abs(robs))
    n <- nrow(a)
    rr <- rep(NA_real_, n)
    for (i in seq_len(n)) {
      lo <- i - h
      hi <- i + h
      if (lo >= 1 && hi <= n) {
        rr[i] <- cor(a[lo:hi, 2], b[lo:hi, 2])
      }
    }
    cat(sprintf("%s\n", lab))
    print(tab, row.names = FALSE)
    cat(sprintf("  32-64 yr: r = %+0.3f, two-sided p = %.4f (%d surrogates); running %d-yr median %+0.2f, same sign in %.0f%% of windows\n\n",
        robs, pv, n_surrogates_bandpass, win, median(rr, na.rm = TRUE),
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
        bar_h = 0.30, y0 = 0.46, cex = 0.65)
    mtext(ttl, side = 3, line = 0.2, adj = 0, cex = 1)
    # scalogram
    par(mar = c(2.8, 5.0, 0.6, 5.0))
    plot_cross_wavelet(K)
    add_label("A")
    # band-passed records
    par(mar = c(2.8, 5.0, 2.0, 5.0))
    plot(S$a[, 1], S$a[, 2], type = "l", lwd = 1.7, col = col1, xaxt = "n",
        xlab = "", ylab = paste0(lab1, " (32-64 yr)"), xaxs = "i",
        cex.lab = 0.9)
    add_bce_ce_axis(250)
    add_label("B")
    par(new = TRUE)
    plot(S$b[, 1], S$b[, 2], type = "l", lwd = 1.7, col = "#B2182B",
        axes = FALSE, xlab = "", ylab = "", xaxs = "i")
    axis(4, col = "#B2182B", col.axis = "#B2182B")
    mtext("Angstel/Vecht GI (32-64 yr)", side = 4, line = 2.7,
        col = "#B2182B", cex = 0.72)
    if (show_stats) {
      mtext(sprintf("r = %+0.2f, p = %.3f", S$r, S$p), side = 3, line = 0.3,
          adj = 1, cex = 0.72, col = "grey25")
    }
    # running correlation
    par(mar = c(2.8, 5.0, 2.0, 5.0))
    plot(S$a[, 1], S$run, type = "l", lwd = 1.7, ylim = c(-1, 1), xaxt = "n",
        xlab = "", ylab = "201-yr running correlation", xaxs = "i",
        cex.lab = 0.9)
    abline(h = 0, lty = 2)
    add_bce_ce_axis(250)
    add_label("C")
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
    d <- numeric(length(idx))
    for (j in seq_along(idx)) {
      i <- idx[j]
      d[j] <- mean(sm[(i + 1):max(ok)], na.rm = TRUE) - mean(sm[min(ok):i], na.rm = TRUE)
    }
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

## Comparison records on the common window 890 BCE - 156 CE ####
{
  stal_url <- "https://raw.githubusercontent.com/stratigraphy/Weesp-dendro-cycles/main/stal_scot.csv"
  saw_url  <- paste0("https://www.ncei.noaa.gov/pub/data/paleo/paleolimnology/",
      "northamerica/canada/ellesmere/sawtooth2020-noaa.txt")

  have_stal <- download_if_missing(stal_url, in_data("stal_scot.csv"))
  have_saw  <- download_if_missing(saw_url,  in_data("sawtooth2020.txt"))
  if (!have_stal || !have_saw) {
    stop("The comparison needs stal_scot.csv and sawtooth2020.txt; ",
        "download them from the URLs given in the header and rerun.")
  }

  stal <- read.csv(in_data("stal_scot.csv"))[, c(2, 3)]
  stal <- stal[stats::complete.cases(stal), ]
  # The record is in historical years (no year 0); BCE years are moved up by one so that
  # 1 BCE and 1 CE are consecutive, as for the growth index, and no year 0 is interpolated
  stal[, 1] <- ifelse(stal[, 1] < 0, stal[, 1] + 1, stal[, 1])
  stal <- linterp(sortNave(stal, genplot = FALSE, verbose = FALSE),
      genplot = FALSE, verbose = FALSE)
  names(stal) <- c("x", "y")

  raw <- readLines(in_data("sawtooth2020.txt"), warn = FALSE)
  raw <- raw[!grepl("^#", raw) & !grepl("^depth", raw)]
  sp  <- strsplit(trimws(raw), "[ \t]+")
  sx <- rep(NA_character_, length(sp))
  sy <- rep(NA_character_, length(sp))
  for (i_sp in seq_along(sp)) {
    if (length(sp[[i_sp]]) >= 5) {
      sx[i_sp] <- sp[[i_sp]][2]
      sy[i_sp] <- sp[[i_sp]][5]
    }
  }
  keep_rows <- !is.na(sx)
  saw <- data.frame(x = as.numeric(sx[keep_rows]), y = as.numeric(sy[keep_rows]))
  saw <- saw[saw$y > -999, ]
  saw <- saw[order(saw$x), ]
  # The Sawtooth years run through a year 0 (yearCE = 1950 - ageBP), so they are already
  # plot years: -889 is 890 BCE, the first year of the common window

  GIraw <- VA25mean
  names(GIraw) <- c("x", "y")

  gmin <- max(min(stal$x), min(saw$x), min(GIraw$x))
  gmax <- min(max(stal$x), max(saw$x), max(GIraw$x))
  grid <- seq(gmin, gmax, by = 1)
  cat(sprintf("\nCommon window for all comparisons: %d to %d = %d years\n",
      gmin, gmax, length(grid)))
  on_grid <- function(d) data.frame(
      x = grid, y = as.numeric(scale(approx(d$x, d$y, grid)$y)))
  ST <- on_grid(stal)
  SW <- on_grid(saw)
  GI <- on_grid(GIraw)
}


# MAIN TEXT ####

# Figure 2: growth index and sample depth ####
# Output: figures_pdf/Figure2.pdf, figures_png/Figure2.png
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

  save_figure(record_figure(plot_figure2, 12, 7.2), "Figure2", 12, 7.2)
}


# Table 1: significance of the spectral peaks ####
# Maximum-in-band test against four backgrounds
# Output: output/tables/Table_1_cycle_significance.csv
{
  x <- VA25mean[, 2]
  n <- length(x)

  # Transform used for the statistics
  w_test <- morlet_wavelet(x, dj = dj_test)
  gws_obs <- average_spectral_power(w_test)

  backgrounds <- make_backgrounds(x, w_test$period, gws_obs, w_test$scales)

  cat("\nTable 1: surrogate significance\n")
  sig <- run_or_load(paste0("Table1_", chronology, "_n", n_surrogates), {
    sig <- vector("list", length(backgrounds))
    names(sig) <- names(backgrounds)
    for (bn in names(backgrounds)) {
      sig[[bn]] <- surrogate_significance(backgrounds[[bn]])
    }
    sig
  })

  cycle_table <- vector("list", length(sig))
  for (j_sig in seq_along(sig)) {
    cycle_table[[j_sig]] <- sig[[j_sig]]$cand
  }
  cycle_table <- do.call(rbind, cycle_table)
  rownames(cycle_table) <- NULL
  write.csv(cycle_table,
      in_out("tables", "Table_1_cycle_significance.csv"), row.names = FALSE)

  print(reshape(cycle_table[, c("Candidate", "Background", "p_max_in_window")],
      idvar = "Candidate", timevar = "Background", direction = "wide"),
      row.names = FALSE)
}


# Figure 3: wavelet scalogram with significance ####
# Output: figures_pdf/Figure3.pdf, figures_png/Figure3.png
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
          lwd = if (bn == "ar1_observed") { 1.6 } else { 1.1 },
          lty = if (bn == "ar1_observed") { 1 } else { 2 },
          col = if (bn == "ar1_observed") { "black" } else { "grey25" })
    }

    # Cone of influence
    shade_coi(growth$PlotYear, w_test$coi, 2, max(w_test$period))

    # Average spectral power with the 95% level of each background
    par(mar = c(4, 0, 0, 1))
    plot(gws_obs, log2(w_test$period), type = "l", lwd = 2, yaxs = "i",
        yaxt = "n", xlab = "Avg. power", ylab = "")
    cols <- c(white = "grey50", ar1_observed = "firebrick",
        ar1_corrected = "darkorange", powerlaw = "steelblue")
    for (bn in names(cols)) {
      lines(sig[[bn]]$gws95, log2(w_test$period), col = cols[bn], lty = 2)
    }
    abline(h = log2(candidate_periods), lty = 3, col = "grey80")

    par(mar = c(0, 0.5, 2, 1))
    plot.new()
    legend("top", bty = "n", cex = 0.7, lwd = c(2.2, 1.2, 1.2, 1.2, 1.2),
        lty = c(1, 2, 2, 2, 2), col = c("black", cols),
        legend = c("observed", "white noise 95%", expression(paste("AR(1) ", phi, " = 0.662, 95%")),
            expression(paste("AR(1) ", phi, " = 0.510, 95%")), "power law 95%"))
    # Colour key
    par(new = TRUE, mar = c(0.2, 0.8, 7.6, 1.2))
    draw_colour_key(brk, lab = "Scalogram power", frac = 0.86, bar_h = 0.26, cex = 0.6,
        y0 = 0.50, align = "centre", lab_pos = "below", n_lab = 2)

  }

  save_figure(record_figure(plot_figure3, 12.8, 8), "Figure3", 12.8, 8)
}


# Figure 4: extracted cycles and their share of spectral power ####
# Cycles supported by at least three of the four backgrounds (Table 1), and the borderline ~24 yr cycle
# Output: figures_pdf/Figure4.pdf, figures_png/Figure4.png
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
  share <- matrix(NA_real_, nrow = ncol(w_plot$power), ncol = length(panels),
      dimnames = list(NULL, panels))
  for (panel_k in panels) {
    rr_k <- which(w_plot$period >= as.numeric(panel_k) * 0.8 & w_plot$period <= as.numeric(panel_k) * 1.2)
    share[, panel_k] <- colSums(w_plot$power[rr_k, , drop = FALSE]) / colSums(w_plot$power)
  }

  cols4 <- setNames(
      if (requireNamespace("scico", quietly = TRUE)) {
        scico::scico(4, palette = "batlow", begin = 0.1, end = 0.85)
      } else { c("#6baed6", "#74c476", "#9e9ac8", "#fd8d3c") },
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

  save_figure(record_figure(plot_figure4, 12, 9.2), "Figure4", 12, 9.2)
}

# Section 4.3: persistence (area-wise test) ####
# Longest run above the 95% level per octave band, against the two AR(1) backgrounds
# Output: output/tables/Section_4.3_persistence_test.csv
{
  band_def <- vector("list", nrow(octave_bands))
  for (b in seq_len(nrow(octave_bands))) {
    last <- b == nrow(octave_bands)
    rr <- if (last) {
      which(w_test$period >= octave_bands$Period_min[b] &
          w_test$period <= octave_bands$Period_max[b])
    } else {
      which(w_test$period >= octave_bands$Period_min[b] &
          w_test$period < octave_bands$Period_max[b])
    }
    if (!length(rr)) {
      next
    }
    # years in which the whole band lies outside the cone of influence
    cc <- which(colSums(!w_test$valid[rr, , drop = FALSE]) == 0)
    if (!length(cc)) {
      next
    }
    band_def[b] <- list(list(rows = rr, cols = cc,
        value = colMeans(w_test$power[rr, cc, drop = FALSE])))
  }

  areawise_table <- run_or_load(paste0("Persistence_", chronology, "_n", n_surrogate_fields), {
    areawise_table <- list()
    ai <- 0

    for (bn in c("ar1_observed", "ar1_corrected")) {

      sim_bands <- vector("list", n_surrogate_fields)

      pb <- progress_start(n_surrogate_fields, paste("Persistence test, background:", backgrounds[[bn]]$label))
      for (s in seq_len(n_surrogate_fields)) {
        ws <- morlet_wavelet(backgrounds[[bn]]$sim(n, backgrounds[[bn]]$par),
            dj = dj_test)
        sb <- vector("list", length(band_def))
        for (b in seq_along(band_def)) {
          bd <- band_def[[b]]
          if (!is.null(bd)) {
            sb[b] <- list(colMeans(ws$power[bd$rows, bd$cols, drop = FALSE]))
          }
        }
        sim_bands[[s]] <- sb
        setTxtProgressBar(pb, s)
      }
      close(pb)

      for (b in seq_len(nrow(octave_bands))) {
        bd <- band_def[[b]]
        if (is.null(bd)) {
          next
        }
        pooled <- vector("list", length(sim_bands))
        for (s in seq_along(sim_bands)) {
          pooled[[s]] <- sim_bands[[s]][[b]]
        }
        crit <- quantile(unlist(pooled), 0.95, names = FALSE)
        obs_run <- longest_run_above(bd$value, crit)
        sim_runs <- numeric(length(sim_bands))
        for (s in seq_along(sim_bands)) {
          sim_runs[s] <- longest_run_above(sim_bands[[s]][[b]], crit)
        }
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
  write.csv(areawise_table, in_out("tables", "Section_4.3_persistence_test.csv"),
      row.names = FALSE)
  print(areawise_table, row.names = FALSE)
}


# Figure 5: bicoherence ####
# Output: figures_pdf/Figure5.pdf, figures_png/Figure5.png
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

    # Main panel
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
      x0 <- 1 / ij[1]
      y0 <- 1 / ij[2]
      xe <- max(bic$freq) * 1.2
      segments(x0, y0, xe, y0 - (xe - x0), col = "grey10", lwd = 2, lty = 2)
      points(x0, y0, pch = 21, bg = "white", col = "grey10", cex = 1.2)
      points(y0, x0, pch = 21, bg = "white", col = "grey10", cex = 1.2)
    }
    box()

    # Top panel, power spectrum with period labels
    par(mar = c(0, 5, 2, 0))
    plot(bic$freq, bic$pwr, type = "l", lwd = 2, log = "y", xaxs = "i",
        xaxt = "n", xlab = "", ylab = "Power")
    abline(v = per_lines, col = "#80808038", lwd = 5)
    for (i in seq_along(per_lab)) {
      mtext(per_lab[i], side = 3, at = per_lines[i], line = 0.15, cex = 0.7)
    }
    box()

    # Right panel, power spectrum rotated
    par(mar = c(5, 0, 0, 2))
    plot(bic$pwr, bic$freq, type = "l", lwd = 2, log = "x", yaxs = "i",
        yaxt = "n", xlab = "Power", ylab = "")
    abline(h = per_lines, col = "#80808038", lwd = 5)
    for (i in seq_along(per_lab)) {
      mtext(per_lab[i], side = 4, at = per_lines[i], line = 0.3, las = 1,
          cex = 0.7)
    }
    box()

    # Colour key
    par(mar = c(2, 2, 3, 2))
    draw_colour_key(quantile_breaks(bic$bic2), lab = "Magnitude-squared bicoherence",
        frac = 0.60, bar_h = 0.16, y0 = 0.55, align = "centre",
        lab_pos = "below", n_lab = 4)

  }

  save_figure(record_figure(plot_figure5, 11.2, 10.4), "Figure5", 11.2, 10.4)
}


# Figure 6: combination tones ####
# Output: figures_pdf/Figure6.pdf, figures_png/Figure6.png
{
  comb_37 <- combination_tone(60, 100)   # -> 37.5 yr
  comb_15 <- combination_tone(24, 40)    # -> 15.0 yr

  # All curves span the full interval; every panel and its second axis use the same x-limits
  xlims <- range(VA25mean[, 1])

  plot_figure6 <- function() {

    par(mfrow = c(3, 1), mar = c(4, 4, 2, 4))

    plot(VA25mean, type = "l", xaxt = "n", xlab = "", ylab = "GI", col = "grey65",
        lwd = 0.6, xlim = xlims, xaxs = "i", cex.main = 1,
        main = "Growth index with the 15.0- and 37.5-yr components")
    add_label("A")
    lines(comb_15$direct[, 1], comb_15$direct[, 2] + mean(x), col = "#2C7BB6", lwd = 1.6)
    lines(comb_37$direct[, 1], comb_37$direct[, 2] + mean(x), col = "#B2182B", lwd = 1.6)
    add_bce_ce_axis()

    plot(comb_37$direct[, 1], comb_37$direct[, 2] + mean(x), type = "l", lwd = 1.8,
        col = "#B2182B", xaxt = "n", xlab = "", ylab = "GI", xlim = xlims, xaxs = "i",
        cex.main = 1,
        main = sprintf("~%.0f yr: extracted (red) vs synthesised from 60 and 100 yr (sum tone) (black)",
            comb_37$period))
    add_label("B")
    add_bce_ce_axis()
    par(new = TRUE)
    plot(comb_37$synthetic, type = "l", axes = FALSE, xlab = "", ylab = "",
        xlim = xlims, xaxs = "i")
    axis(4)
    mtext("Nonlinear combination", side = 4, line = 2.5, cex = 0.8)

    plot(comb_15$direct[, 1], comb_15$direct[, 2] + mean(x), type = "l", lwd = 1.8,
        col = "#2C7BB6", xaxt = "n", xlab = "Year (BCE/CE)", ylab = "GI",
        xlim = xlims, xaxs = "i", cex.main = 1,
        main = sprintf("~%.0f yr: extracted (blue) vs synthesised from 24 and 40 yr (black)",
            comb_15$period))
    add_label("C")
    add_bce_ce_axis()
    par(new = TRUE)
    plot(comb_15$synthetic, type = "l", axes = FALSE, xlab = "", ylab = "",
        xlim = xlims, xaxs = "i")
    axis(4)
    mtext("Nonlinear combination", side = 4, line = 2.5, cex = 0.8)

  }

  save_figure(record_figure(plot_figure6, 12, 8.8), "Figure6", 12, 8.8)
}


# Figure 7: cross-wavelet of Uamh an Tartair and the growth index ####
# Output: figures_pdf/Figure7.pdf, figures_png/Figure7.png
{
  K1 <- cross_wavelet_for_plot(ST, GI)

  plot_figure7 <- function() {
    layout(matrix(1:2, ncol = 1), heights = c(0.16, 1))
    par(oma = c(3.0, 0, 0.2, 0))
    par(mar = c(0.4, 5.0, 1.6, 2.0))
    draw_colour_key(K1$brk, lab = "Cross-wavelet power", frac = 0.34,
        bar_h = 0.30, y0 = 0.46)
    mtext("Uamh an Tartair stalagmite growth rate vs Angstel/Vecht GI, 890 BCE-156 CE",
        side = 3, line = 0.2, adj = 0, cex = 1)
    par(mar = c(0.6, 5.0, 0.6, 2.0))
    plot_cross_wavelet(K1)
    mtext("Year (BCE/CE)", side = 1, outer = TRUE, line = 1.6, cex = 0.8)
  }

  save_figure(record_figure(plot_figure7, 12, 6.8), "Figure7", 12, 6.8)
}


# SUPPLEMENTARY MATERIAL ####

# Figure S1: sample depth, running variance, EPS and 32-64 yr power (S1) ####
# Output: figures_pdf/FigureS1.pdf, figures_png/FigureS1.png
{
  w_s1 <- morlet_wavelet(DG$std050, dj = dj_test)
  rows <- which(w_s1$period >= 32 & w_s1$period <= 64)
  DG$bp3264 <- NA_real_
  for (j_dg in seq_len(nrow(DG))) {
    v <- w_s1$valid[rows, j_dg]
    if (any(v)) {
      DG$bp3264[j_dg] <- mean(w_s1$power[rows[v], j_dg])
    }
  }

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
    abline(h = c(10, 40), lty = 3)
    add_label("A")
    add_bce_ce_axis()
    yl <- range(c(DG$var_std, DG$var_vsc), na.rm = TRUE)
    plot(DG$PlotYear, DG$var_std, type = "l", lwd = 1.2, col = "black", ylim = yl,
        xaxt = "n", xaxs = "i", xlab = "", ylab = "Variance, 101 yr")
    lines(DG$PlotYear, DG$var_vsc, lwd = 1.2, col = "#B2182B")
    legend("topright", c("std050", "vsc101"), col = c("black", "#B2182B"),
        lwd = 1.2, bty = "n", cex = 0.8)
    add_label("B")
    add_bce_ce_axis()
    plot(DG$PlotYear, DG$EPS, type = "l", lwd = 1.2, col = "grey25", ylim = c(0.6, 1),
        xaxt = "n", xaxs = "i", xlab = "", ylab = "EPS")
    abline(h = 0.85, lty = 3)
    add_label("C")
    add_bce_ce_axis()
    plot(DG$PlotYear, DG$bp3264, type = "l", lwd = 1.2, col = "#1B7837",
        xaxt = "n", xaxs = "i", xlab = "", ylab = "Power, 32-64 yr")
    add_label("D")
    add_bce_ce_axis()
    mtext("Year", side = 1, outer = TRUE, line = 1.5, cex = 0.8)
  }

  save_figure(record_figure(plot_figures1, 11.2, 9.6), "FigureS1", 11.2, 9.6)
}


# Table S5: multitaper cross-check (S2) ####
# Full analysis window, robust red noise (Mann and Lees, 1996) and power law;
# the harmonic F-test is listed but not used
# Output: output/tables/Table_S5_multitaper_crosscheck.csv
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
  tab <- vector("list", length(cand))
  fr <- ml$Frequency
  for (j_c in seq_along(cand)) {
    P_c <- cand[j_c]
    sel <- which(1 / fr >= P_c * 0.9 & 1 / fr <= P_c * 1.1)
    if (!length(sel)) {
      next
    }
    i_c <- sel[which.max(ml$Power[sel])]
    tab[[j_c]] <- data.frame(Candidate = P_c,
        MTM_peak  = round(1 / fr[i_c], 2),
        ML96_redCL = round(max(ml$AR1_CL[sel], na.rm = TRUE), 1),
        PL_CL      = round(max(pl$PowerLaw_CL[sel], na.rm = TRUE), 1),
        HarmF_CL   = round(max(ml$Harmonic_CL[sel], na.rm = TRUE), 1))
  }
  tab <- do.call(rbind, tab)
  print(tab, row.names = FALSE)
  write.csv(tab, in_out("tables", "Table_S5_multitaper_crosscheck.csv"), row.names = FALSE)
}


# Figure S2: the three chronologies over 1209-40 BCE (S5) ####
# Output: figures_pdf/FigureS2.pdf, figures_png/FigureS2.png
{
  core <- DG[DG$Year >= -1209 & DG$Year <= -40, ]
  chron <- c(std050 = "std050", vsc101 = "vsc101",
      constant40 = "constant40_median")
  gws3 <- vector("list", length(chron))
  names(gws3) <- names(chron)
  for (cn3 in names(chron)) {
    w3 <- morlet_wavelet(core[[chron[[cn3]]]], dj = dj_figure)
    gws3[[cn3]] <- list(period = w3$period, gws = average_spectral_power(w3))
  }
  cat("\nFigure S2: correlation of CWT average spectral power, 1209-40 BCE\n")
  cat(sprintf("  std050 vs constant40 : %.4f\n", cor(gws3$std050$gws, gws3$constant40$gws, use = "complete.obs")))
  cat(sprintf("  std050 vs vsc101     : %.4f\n", cor(gws3$std050$gws, gws3$vsc101$gws, use = "complete.obs")))
  cat(sprintf("  vsc101 vs constant40 : %.4f\n", cor(gws3$vsc101$gws, gws3$constant40$gws, use = "complete.obs")))

  plot_figures2 <- function() {
    par(mar = c(4.5, 4.5, 1, 1))
    cols3 <- c("black", "#B2182B", "#2166AC")
    yl <- range(c(gws3$std050$gws, gws3$vsc101$gws, gws3$constant40$gws), na.rm = TRUE)
    plot(gws3$std050$period, gws3$std050$gws, type = "n", log = "x", ylim = yl,
        xlab = "Period (years)", ylab = "CWT average spectral power",
        xaxt = "n")
    axis(1, at = c(2, 4, 8, 16, 32, 64, 128, 256))
    # constant40 is drawn last, dashed, because it lies on top of std050
    for (i in 1:3) {
      lines(gws3[[i]]$period, gws3[[i]]$gws, col = cols3[i],
          lwd = c(2.6, 1.4, 1.6)[i], lty = c(1, 1, 2)[i])
    }
    legend("topright", c("std050 (original)", "vsc101 (variance-stabilised)",
        "constant40 (40 series per year)"),
        col = cols3, lwd = c(2.6, 1.4, 1.6), lty = c(1, 1, 2), bty = "n", cex = 0.85)
  }

  save_figure(record_figure(plot_figures2, 11.2, 7.2), "FigureS2", 11.2, 7.2)
}


# Figure S3: band-passed stalagmite record and growth index (S7) ####
# Band-pass correlations for all three pairings, also used for Figure S4 and Table S14
# Output: figures_pdf/FigureS3.pdf, figures_png/FigureS3.png
{
  # Stalagmite vs Sawtooth (S3) is not plotted, it is reported in Table S14
  bp_all <- run_or_load(paste0("FigureS3_bandpass_", chronology, "_n", n_surrogates_bandpass), list(
      S1 = bandpass_comparison(ST, GI, "Uamh an Tartair stalagmite growth rate vs Angstel/Vecht GI"),
      S2 = bandpass_comparison(SW, GI, "Sawtooth Lake AMVI vs Angstel/Vecht GI"),
      S3 = bandpass_comparison(ST, SW, "Uamh an Tartair stalagmite growth rate vs Sawtooth Lake AMVI")))
  S1 <- bp_all$S1
  S2 <- bp_all$S2
  S3 <- bp_all$S3

  plot_figures3 <- function() {
    par(mfrow = c(2, 1), oma = c(1.6, 0, 0, 0))
    par(mar = c(2.8, 5.0, 2.0, 5.0))
    plot(S1$a[, 1], S1$a[, 2], type = "l", lwd = 1.7, col = "grey25", xaxt = "n",
        xlab = "", ylab = "Stalagmite growth rate (32-64 yr)", xaxs = "i", cex.lab = 0.9)
    add_bce_ce_axis(250)
    add_label("A")
    par(new = TRUE)
    plot(S1$b[, 1], S1$b[, 2], type = "l", lwd = 1.7, col = "#B2182B",
        axes = FALSE, xlab = "", ylab = "", xaxs = "i")
    axis(4, col = "#B2182B", col.axis = "#B2182B")
    mtext("Angstel/Vecht GI (32-64 yr)", side = 4, line = 2.7, col = "#B2182B", cex = 0.72)
    plot(S1$a[, 1], S1$run, type = "l", lwd = 1.7, ylim = c(-1, 1), xaxt = "n",
        xlab = "", ylab = "201-yr running correlation", xaxs = "i", cex.lab = 0.9)
    abline(h = 0, lty = 2)
    add_bce_ce_axis(250)
    add_label("B")
    mtext("Year (BCE/CE)", side = 1, outer = TRUE, line = 0.4, cex = 0.8)
  }

  save_figure(record_figure(plot_figures3, 12, 8.4), "FigureS3", 12, 8.4)
}


# Table S9: phase on the cross-power ridges, by band (S7) ####
# Ridges are local maxima of cross-power along the period axis, above the 80th percentile
# of cross-power and outside the cone of influence; R is the mean resultant length
# Output: output/tables/Table_S9_ridge_phase.csv
{
  X9 <- cross_wavelet(cbind(ST$x, ST$y), cbind(GI$x, GI$y), dj = dj_figure,
      lowerPeriod = 1.5, upperPeriod = 192, omega_nr = 15)
  Pw9 <- t(X9$Power)
  Ph9 <- t(X9$Phase)
  up9 <- rbind(Pw9[-1, , drop = FALSE], -Inf)
  dn9 <- rbind(Inf, Pw9[-X9$nr, , drop = FALSE])
  ridge9 <- Pw9 > up9 & Pw9 > dn9
  ridge9[1, ] <- FALSE
  ridge9[X9$nr, ] <- FALSE
  ridge9 <- ridge9 & (Pw9 >= quantile(X9$Power, 0.8, na.rm = TRUE)) & X9$valid
  bands9 <- list(c(8, 16), c(16, 32), c(32, 64), c(64, 128))
  tabS9 <- data.frame(Band = character(length(bands9)), Ridge_cells = integer(length(bands9)),
      Mean_phase_deg = NA_real_, R = NA_real_, stringsAsFactors = FALSE)
  for (k9 in seq_along(bands9)) {
    rows9 <- which(X9$Period >= bands9[[k9]][1] & X9$Period <= bands9[[k9]][2])
    ph9 <- Ph9[rows9, ][ridge9[rows9, ]]
    z9 <- mean(exp(1i * ph9))
    tabS9$Band[k9] <- sprintf("%d-%d", bands9[[k9]][1], bands9[[k9]][2])
    tabS9$Ridge_cells[k9] <- length(ph9)
    tabS9$Mean_phase_deg[k9] <- round(Arg(z9) * 180 / pi)
    tabS9$R[k9] <- round(Mod(z9), 2)
  }
  cat("\nTable S9: phase on the cross-power ridges, stalagmite x GI (positive: stalagmite leads)\n")
  print(tabS9, row.names = FALSE)
  write.csv(tabS9, in_out("tables", "Table_S9_ridge_phase.csv"), row.names = FALSE)
}


# Table S10: band-pass correlation over the whole stalagmite overlap (S7) ####
# The stalagmite record and the growth index over their full overlap (991 BCE - 156 CE),
# lags within +/- 32 yr, positive when the growth index leads
# Output: output/tables/Table_S10_bandpass_full_overlap.csv
{
  g10 <- seq(max(min(stal$x), min(GIraw$x)), min(max(stal$x), max(GIraw$x)), by = 1)
  ST10 <- data.frame(x = g10, y = as.numeric(scale(approx(stal$x, stal$y, g10)$y)))
  GI10 <- data.frame(x = g10, y = as.numeric(scale(approx(GIraw$x, GIraw$y, g10)$y)))
  cat(sprintf("\nTable S10: stalagmite and growth index overlap %d to %d = %d years\n",
      min(g10), max(g10), length(g10)))
  tabS10 <- vector("list", length(bands))
  for (k10 in seq_along(bands)) {
    a10 <- bandpass(ST10, bands[[k10]][1], bands[[k10]][2])[, 2]
    b10 <- bandpass(GI10, bands[[k10]][1], bands[[k10]][2])[, 2]
    cc10 <- ccf(a10, b10, lag.max = 32, plot = FALSE)
    tabS10[[k10]] <- data.frame(
        Band = sprintf("%d-%d", bands[[k10]][1], bands[[k10]][2]),
        r0 = round(cor(a10, b10), 2),
        min_r = round(min(cc10$acf), 2), lag_min = cc10$lag[which.min(cc10$acf)],
        max_r = round(max(cc10$acf), 2), lag_max = cc10$lag[which.max(cc10$acf)])
  }
  tabS10 <- do.call(rbind, tabS10)
  print(tabS10, row.names = FALSE)
  write.csv(tabS10, in_out("tables", "Table_S10_bandpass_full_overlap.csv"), row.names = FALSE)

  # 201-yr running correlation of the 32-64 yr band over the full overlap
  a10 <- bandpass(ST10, 32, 64)[, 2]
  b10 <- bandpass(GI10, 32, 64)[, 2]
  rr10 <- rep(NA_real_, length(a10))
  for (i10 in seq_along(a10)) {
    if (i10 - h >= 1 && i10 + h <= length(a10)) {
      rr10[i10] <- cor(a10[(i10 - h):(i10 + h)], b10[(i10 - h):(i10 + h)])
    }
  }
  r10 <- cor(a10, b10)
  cat(sprintf("  32-64 yr over the full overlap: r = %+0.3f; running %d-yr median %+0.2f, same sign in %.0f%% of windows\n",
      r10, win, median(rr10, na.rm = TRUE), 100 * mean(sign(rr10) == sign(r10), na.rm = TRUE)))
}


# S7: common window restricted to years with at least 40 series ####
# Zero-lag band-pass correlations over the common window up to 40 BCE (plot year -39),
# which leaves out the poorly replicated final segment
# Output: output/tables/S7_restricted_window.csv
{
  ST40 <- ST[ST$x <= -39, ]
  GI40 <- GI[GI$x <= -39, ]
  bands40 <- list(c(32, 64), c(32, 48), c(48, 80))
  tab40 <- data.frame(Band = character(length(bands40)), r0_restricted = NA_real_,
      r0_common_window = NA_real_, stringsAsFactors = FALSE)
  for (k40 in seq_along(bands40)) {
    lo40 <- bands40[[k40]][1]
    hi40 <- bands40[[k40]][2]
    tab40$Band[k40] <- sprintf("%d-%d", lo40, hi40)
    tab40$r0_restricted[k40] <- round(cor(bandpass(ST40, lo40, hi40)[, 2],
        bandpass(GI40, lo40, hi40)[, 2]), 3)
    tab40$r0_common_window[k40] <- round(cor(bandpass(ST, lo40, hi40)[, 2],
        bandpass(GI, lo40, hi40)[, 2]), 3)
  }
  cat(sprintf("\nStalagmite vs GI, restricted window %d to %d (plot years)\n", min(ST40$x), max(ST40$x)))
  print(tab40, row.names = FALSE)
  write.csv(tab40, in_out("tables", "S7_restricted_window.csv"), row.names = FALSE)
}


# Table S11: cross-wavelet phase by period (S7) ####
# Output: output/tables/Table_S11_phase_by_period.csv
{
  pST <- phase_by_period(ST, GI)
  pSW <- phase_by_period(SW, GI)
  tabS11 <- data.frame(Period = sprintf("%d-%d", head(edges, -1), edges[-1]),
      ST_GI_phase = round(pST[, "phase"]), ST_GI_power = signif(pST[, "power"], 3),
      ST_GI_phase_unrect = round(pST[, "phase_unrect"]),
      ST_GI_power_unrect = signif(pST[, "power_unrect"], 3),
      SW_GI_phase = round(pSW[, "phase"]), SW_GI_power = signif(pSW[, "power"], 3),
      SW_GI_phase_unrect = round(pSW[, "phase_unrect"]),
      SW_GI_power_unrect = signif(pSW[, "power_unrect"], 3))
  write.csv(tabS11, in_out("tables", "Table_S11_phase_by_period.csv"), row.names = FALSE)
  print(tabS11, row.names = FALSE)
}


# Table S12: band-pass correlation in two sub-bands (S7) ####
# The split at 48 yr was made after looking at the phase, so these p-values are descriptive
# Output: output/tables/Table_S12_subband_correlation.csv
{
  tabS12 <- run_or_load(paste0("TableS12_", chronology, "_n", n_surrogates_bandpass), local({
    pairs12 <- list(list("Stalagmite", ST), list("Sawtooth", SW))
    bands12 <- list(c(32, 48), c(48, 80))
    rows12 <- list()
    for (pp in pairs12) {
      for (b in bands12) {
        g  <- bandpass(GI, b[1], b[2])[, 2]
        a  <- bandpass(pp[[2]], b[1], b[2])[, 2]
        cc <- ccf(a, g, lag.max = 40, plot = FALSE)
        r0 <- cor(a, g)
        m0 <- max_abs_ccf(a, g)
        nul <- replicate_progress(n_surrogates_bandpass,
            sprintf("Table S12, %s, %d-%d yr", pp[[1]], b[1], b[2]), {
              as <- bandpass(data.frame(x = pp[[2]]$x, y = phase_randomise(pp[[2]]$y)), b[1], b[2])[, 2]
              c(cor(as, g), max_abs_ccf(as, g))
            })
        rows12[[length(rows12) + 1]] <- data.frame(
            Pair = pp[[1]], Band = sprintf("%d-%d", b[1], b[2]),
            r0 = round(r0, 3), p_r0 = mean(abs(nul[1, ]) >= abs(r0)),
            min_r = round(min(cc$acf), 3), lag_min = cc$lag[which.min(cc$acf)],
            max_r = round(max(cc$acf), 3), lag_max = cc$lag[which.max(cc$acf)],
            maxabs_20 = round(m0, 3), p_maxabs = mean(nul[2, ] >= m0))
      }
    }
    do.call(rbind, rows12)
  }))
  write.csv(tabS12, in_out("tables", "Table_S12_subband_correlation.csv"), row.names = FALSE)
  print(tabS12, row.names = FALSE)
}


# S7: single delay fitted to the phase and lagged correlation ####
# A constant delay D rotates the phase with period as 360 D / P, fitted over 28-70 yr;
# the largest lagged correlation (32-80 yr) is tested with the lag search included
# Output: output/tables/S7_fitted_delay.csv
{
  Xl <- cross_wavelet(cbind(ST$x, ST$y), cbind(GI$x, GI$y), dj = dj_figure,
      lowerPeriod = 1.5, upperPeriod = 192, omega_nr = 15)
  Pl <- matrix(Xl$Period, Xl$nr, Xl$nc)
  ok <- Xl$valid & Pl >= 28 & Pl <= 70
  wl <- t(Xl$Power)[ok]
  phl <- t(Xl$Phase)[ok]
  pel <- Pl[ok]
  delays <- 0:30
  fitD <- numeric(length(delays))
  for (j_d in seq_along(delays)) {
    fitD[j_d] <- sum(wl * cos(phl - 2 * pi * delays[j_d] / pel)) / sum(wl)
  }
  angles <- seq(-pi, pi, length.out = 361)
  fitA <- numeric(length(angles))
  for (j_a in seq_along(angles)) {
    fitA[j_a] <- sum(wl * cos(phl - angles[j_a])) / sum(wl)
  }
  fitC <- max(fitA)
  cat(sprintf("Fitted delay: stalagmite leads by %d yr (mean cos %.2f; best constant phase %.2f)\n",
      (0:30)[which.max(fitD)], max(fitD), fitC))
  a80 <- bandpass(ST, 32, 80)[, 2]
  g80 <- bandpass(GI, 32, 80)[, 2]
  cc80 <- ccf(a80, g80, lag.max = 30, plot = FALSE)
  m80 <- max(cc80$acf)
  n80 <- run_or_load(paste0("S7_lag_", chronology, "_n", n_surrogates_bandpass),
      replicate_progress(n_surrogates_bandpass, "S7, lagged correlation 32-80 yr", {
        as <- bandpass(data.frame(x = ST$x, y = phase_randomise(ST$y)), 32, 80)[, 2]
        max(ccf(as, g80, lag.max = 30, plot = FALSE)$acf)
      }))
  cat(sprintf("32-80 yr: largest lagged r = %+.2f at lag %d (negative = stalagmite leads), p = %.3f\n",
      m80, cc80$lag[which.max(cc80$acf)], mean(n80 >= m80)))
  write.csv(data.frame(fitted_delay_yr = (0:30)[which.max(fitD)], mean_cos_delay = round(max(fitD), 2),
      mean_cos_constant_phase = round(fitC, 2), r_lagged_32_80 = round(m80, 2),
      lag_yr = cc80$lag[which.max(cc80$acf)], p = mean(n80 >= m80)),
      in_out("tables", "S7_fitted_delay.csv"), row.names = FALSE)
}


# Table S13: significance of the Sawtooth Lake AMV index itself (S8) ####
# The maximum-in-band test of Table 1, applied to the Sawtooth record over the common
# window with the four backgrounds fitted to that record. The surrogates use their own
# random number stream (with_seed), so the later sections give the same results.
# Output: output/tables/Table_S13_sawtooth_band_significance.csv
{
  tabS13 <- run_or_load(paste0("TableS13_sawtooth_n", n_surrogates_s13), with_seed(2013, local({
    xs13 <- SW$y
    w13 <- morlet_wavelet(xs13, dj = dj_test)
    gws13 <- average_spectral_power(w13)
    cat("\nTable S13: Sawtooth Lake AMVI, 890 BCE - 156 CE")
    bg13 <- make_backgrounds(xs13, w13$period, gws13, w13$scales)
    bands13 <- list(c(8, 16), c(16, 32), c(32, 64))
    rows13 <- list()
    for (bn13 in names(bg13)) {
      sims13 <- matrix(NA_real_, nrow = length(w13$period), ncol = n_surrogates_s13)
      pb <- progress_start(n_surrogates_s13, paste("Table S13, background:", bg13[[bn13]]$label))
      for (s13 in seq_len(n_surrogates_s13)) {
        ws13 <- morlet_wavelet(bg13[[bn13]]$sim(length(xs13), bg13[[bn13]]$par), dj = dj_test)
        sims13[, s13] <- average_spectral_power(ws13)
        setTxtProgressBar(pb, s13)
      }
      close(pb)
      for (b13 in bands13) {
        rr13 <- which(w13$period >= b13[1] & w13$period <= b13[2])
        obs13 <- max(gws13[rr13], na.rm = TRUE)
        rows13[[length(rows13) + 1]] <- data.frame(
            Band = sprintf("%d-%d", b13[1], b13[2]),
            Peak_period = round(w13$period[rr13][which.max(gws13[rr13])], 1),
            Background = bg13[[bn13]]$label,
            p_max_in_band = mean(column_max(sims13[rr13, , drop = FALSE]) >= obs13),
            stringsAsFactors = FALSE)
      }
    }
    do.call(rbind, rows13)
  })))
  print(reshape(tabS13[, c("Band", "Background", "p_max_in_band")],
      idvar = "Band", timevar = "Background", direction = "wide"), row.names = FALSE)
  write.csv(tabS13, in_out("tables", "Table_S13_sawtooth_band_significance.csv"), row.names = FALSE)
}


# Figure S4: comparison with the Sawtooth Lake AMV index (S8) ####
# Output: figures_pdf/FigureS4.pdf, figures_png/FigureS4.png
{
  K2 <- cross_wavelet_for_plot(SW, GI)

  plot_figures4 <- function() {
    layout(matrix(1:4, ncol = 1), heights = c(0.22, 1, 1, 1))
    par(oma = c(0.6, 0, 0.4, 0))
    plot_sawtooth_panels(S2, K2, "Sawtooth Lake AMVI", "#1B4F72",
        "Sawtooth Lake, Ellesmere Island, vs Angstel/Vecht GI, 890 BCE-156 CE", show_stats = TRUE)
  }

  save_figure(record_figure(plot_figures4, 8, 11.6), "FigureS4", 8, 11.6)
}


# Table S14: band-pass correlations for the three pairings (S8) ####
# Output: output/tables/Table_S14_bandpass_correlations.csv
{
  write.csv(rbind(S1$tab, S2$tab, S3$tab),
      in_out("tables", "Table_S14_bandpass_correlations.csv"),
      row.names = FALSE)
  print(rbind(S1$tab, S2$tab, S3$tab), row.names = FALSE)
}


# Figure S5: change in spectral organisation near 572 BCE (S9) ####
# Full analysis window
# Output: figures_pdf/FigureS5.pdf, figures_png/FigureS5.png, tables/FigS5_band_power_ratio.csv, tables/FigS5_step_test.csv
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
    null_stat <- numeric(nsim)
    null_year <- numeric(nsim)
    pb <- progress_start(nsim, "Figure S5, AR(1) surrogates")
    for (i in seq_len(nsim)) {
      xs <- as.numeric(arima.sim(list(ar = phi), n = length(x),
          sd = sd(x) * sqrt(1 - phi^2))) + mean(x)
      s  <- largest_step(band_power_ratio(xs, yr), yr)
      null_stat[i] <- s$stat
      null_year[i] <- s$year
      setTxtProgressBar(pb, i)
    }
    close(pb)
    list(stat = null_stat, year = null_year)
  })
  null_stat <- step_null$stat
  null_year <- step_null$year
  p <- mean(null_stat >= obs$stat)
  cat(sprintf("AR(1) surrogates (n = %d): median step %.3f, 95th pct %.3f, p = %.3f\n",
      nsim, median(null_stat), quantile(null_stat, 0.95), p))
  cat(sprintf("Surrogate transition years spread over %d to %d BCE (sd %.0f yr)\n",
      -round(historical_year(min(null_year))), -round(historical_year(max(null_year))), sd(null_year)))

  write.csv(data.frame(Year = yr, log_ratio_48_80_over_14_17 = sm),
      in_out("tables", "FigS5_band_power_ratio.csv"), row.names = FALSE)
  write.csv(data.frame(observed_year_BCE = -round(obs_hist), observed_step = round(obs$stat, 3),
      phi = round(phi, 3), n_sim = nsim,
      null_median = round(median(null_stat), 3),
      null_p95 = round(as.numeric(quantile(null_stat, 0.95)), 3),
      p_value = round(p, 3),
      null_year_min_BCE = -round(historical_year(min(null_year))),
      null_year_max_BCE = -round(historical_year(max(null_year))),
      null_year_sd = round(sd(null_year))),
      in_out("tables", "FigS5_step_test.csv"), row.names = FALSE)

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

  save_figure(record_figure(plot_figures5, 10.4, 6), "FigureS5", 10.4, 6)
}


# Table S15: sensitivity of the AR(1) coefficient (S10) ####
# Full analysis window (x from Figure S5)
# Output: output/tables/Table_S15_ar1_sensitivity.csv
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
  for (nm in names(sets)) {
    cat(sprintf("%-58s phi = %.3f\n", nm, ar1_coefficient(remove_bands(x, sets[[nm]]))))
  }
  cat("\nNote: adding the 2.2-year band RAISES phi. Removing variance close to\n")
  cat("the Nyquist period lengthens the effective decorrelation time, so the\n")
  cat("near-Nyquist band is left in the series when the background is fitted.\n")
  cat("\nSensitivity to the +/- tolerance around each notched period\n")
  for (tol in c(0.05,0.075,0.10,0.15,0.20)) {
    cat(sprintf("  tol = +/-%4.1f%%  (fraction of spectrum removed ~%.0f%%)   phi = %.3f\n",
        100*tol, 100*min(1, 2*tol*10), ar1_coefficient(remove_bands(x, c(3.6,5.9,8,12,15.4,24,38,60,80,110), tol = tol))))
  }
  cat("\nphi falls as the notch widens, because more of the series is removed.\n")
  cat("The candidate-band-excluded AR(1) background is therefore a sensitivity\n")
  cat("scenario, not an estimate of the true background.\n")

  tols <- c(0.05, 0.075, 0.10, 0.15, 0.20)
  phi_sets <- numeric(length(sets))
  for (j_s in seq_along(sets)) {
    phi_sets[j_s] <- ar1_coefficient(remove_bands(x, sets[[j_s]]))
  }
  phi_tols <- numeric(length(tols))
  for (j_t in seq_along(tols)) {
    phi_tols[j_t] <- ar1_coefficient(remove_bands(x, notch_periods, tol = tols[j_t]))
  }
  tabS15 <- data.frame(
      Estimated_on = c("observed series, no notching", names(sets),
          sprintf("notch 3.6-110 yr, +/- %g%%", 100 * tols)),
      phi = round(c(ar1_coefficient(x), phi_sets, phi_tols), 3))
  write.csv(tabS15, in_out("tables", "Table_S15_ar1_sensitivity.csv"), row.names = FALSE)
}


# SUPPLEMENTARY SUITE AT OMEGA0 = 6 ####

# Tables S1-S4 and S6-S8: four backgrounds, three chronologies and a prewhitened run (S2-S5) ####
# The full factorial of three chronologies, two windows, a prewhitened series and four
# backgrounds, at omega0 = 6 and dj = 1/20. It runs on its own random number stream
# (seed 20092026), so it does not change the results of the other sections.
# Output: output/tables/Table_S1_background_parameters.csv,
#         Tables_S2_S3_S7_S8_candidate_window_tests.csv, Table_S4_octave_band_tests.csv,
#         Table_S6_duration_test_omega6.csv and three diagnostic tables
{
  suite <- run_or_load(paste0("SupplementSuite_omega6_n", n_surrogates_suite, "_", n_fields_suite),
      with_seed(20092026, local({

        AV_nsim_spectral <- n_surrogates_suite   # surrogates for MTM / global wavelet
        AV_nsim_local    <- n_fields_suite       # surrogates for the local (2-D) wavelet field

        AV_mtm_nw <- 4
        AV_mtm_k  <- 7

        AV_dj          <- 1 / 20   # 20 wavelet scales per octave
        AV_period_min  <- 2
        AV_period_max  <- 256
        AV_omega0      <- 6        # Morlet; Torrence & Compo (1998)

        # Candidate periodicities reported in the manuscript.
        # Each is tested inside a +/- 10% window in PERIOD, which is the
        # honest way of testing "a cycle near P" without claiming that the
        # exact value of P was known in advance.
        AV_candidates <- c(3.6, 5.9, 8, 12, 15.4, 24, 38, 60, 80, 110)
        AV_cand_tol   <- 0.10

        # Pre-registered octave bands (as in step 7)
        AV_bands <- data.frame(
            Band       = c("2-4", "4-8", "8-16", "16-32", "32-64", "64-128"),
            Period_min = c(2, 4, 8, 16, 32, 64),
            Period_max = c(4, 8, 16, 32, 64, 128),
            stringsAsFactors = FALSE
        )


        AV_read_input <- function() {

          xls <- data_file

          if (file.exists(xls) && requireNamespace("readxl", quietly = TRUE)) {

            list(
                diag = as.data.frame(readxl::read_excel(xls, sheet = "Diagnostic_data")),
                c40  = as.data.frame(readxl::read_excel(xls, sheet = "Constant40"))
            )

          } else {

            list(
                diag = read.csv(in_data("AV_diagnostic_data.csv"), stringsAsFactors = FALSE),
                c40  = read.csv(in_data("AV_constant40.csv"),      stringsAsFactors = FALSE)
            )
          }
        }

        AV_in <- AV_read_input()

        AV_diag <- AV_in$diag
        AV_c40  <- AV_in$c40

        names(AV_diag)[names(AV_diag) == "samp.depth"] <- "samp_depth"

        stopifnot(all(c("Year", "samp_depth", "Rbar101", "std050", "vsc101") %in%
            names(AV_diag)))


        # Analysis windows
        #
        # FULL : longest continuous stretch with n >= 10, i.e. the published
        #        analysis window 1283 BCE - 156 CE (1439 years).
        # CORE : longest continuous stretch with n >= 40, 1209-40 BCE
        #        (1170 years). This is the window in which sample depth can be
        #        held exactly constant, and is therefore the window in which
        #        Reviewer 1's variance-artefact concern can be tested directly.

        AV_longest_run <- function(flag, years) {
          r <- rle(flag)
          i <- which(r$values)
          if (!length(i)) {
            stop("no years meet the criterion")
          }
          j <- i[which.max(r$lengths[i])]
          end <- cumsum(r$lengths)[j]
          start <- end - r$lengths[j] + 1
          list(idx = start:end, start = years[start], end = years[end],
              n = r$lengths[j])
        }

        AV_full <- AV_longest_run(AV_diag$samp_depth >= 10, AV_diag$Year)
        AV_core <- AV_longest_run(AV_diag$samp_depth >= 40, AV_diag$Year)

        cat("\nFULL window:", AV_full$start, "to", AV_full$end,
            "(", AV_full$n, "years )\n")
        cat("CORE window:", AV_core$start, "to", AV_core$end,
            "(", AV_core$n, "years )\n")

        AV_series <- list(
            std050_full = list(
                x = AV_diag$std050[AV_full$idx],
                year = AV_diag$Year[AV_full$idx],
                label = "std050 (1283 BCE-156 CE)"),
            std050_core = list(
                x = AV_diag$std050[AV_core$idx],
                year = AV_diag$Year[AV_core$idx],
                label = "std050 (1209-40 BCE)"),
            vsc101_core = list(
                x = AV_diag$vsc101[AV_core$idx],
                year = AV_diag$Year[AV_core$idx],
                label = "vsc101 (1209-40 BCE)"),
            constant40_core = list(
                x = AV_c40$constant40_median,
                year = AV_c40$Year,
                label = "constant40 (1209-40 BCE)")
        )

        stopifnot(identical(AV_series$std050_core$year,
            as.numeric(AV_series$constant40_core$year)) ||
            all(AV_series$std050_core$year ==
                AV_series$constant40_core$year))


        # Morlet wavelet transform (Torrence & Compo 1998)
        #
        # Implemented directly so that (a) the cone of influence is explicit,
        # (b) any surrogate series can be pushed through exactly the same
        # transform, and (c) power is normalised by the series variance, so
        # that power is comparable between chronologies with different
        # variance (relevant for Reviewer 1, point 1).

        AV_scales <- function(n, dt = 1) {
          s0 <- AV_period_min * dt / (4 * pi / (AV_omega0 + sqrt(2 + AV_omega0^2)))
          smax <- AV_period_max * dt / (4 * pi / (AV_omega0 + sqrt(2 + AV_omega0^2)))
          J <- floor(log2(smax / s0) / AV_dj)
          s0 * 2^((0:J) * AV_dj)
        }

        AV_cwt <- function(x, dt = 1, scales = NULL) {

          n <- length(x)
          x <- x - mean(x)
          sdx <- sd(x)

          if (is.null(scales)) {
            scales <- AV_scales(n, dt)
          }

          # zero padding to the next power of two (reduces wrap-around)
          base2 <- 2^ceiling(log2(n))
          xpad <- c(x, rep(0, base2 - n))
          N <- length(xpad)

          k <- (1:floor(N / 2)) * (2 * pi) / (N * dt)
          k <- c(0, k, -k[(floor((N - 1) / 2)):1])

          f <- fft(xpad)

          W <- matrix(complex(real = 0, imaginary = 0), nrow = length(scales), ncol = n)

          for (i in seq_along(scales)) {

            s <- scales[i]
            expnt <- -((s * k - AV_omega0)^2) / 2 * (k > 0)
            norm <- sqrt(2 * pi * s / dt) * (pi^(-0.25))
            daughter <- norm * exp(expnt) * (k > 0)
            W[i, ] <- fft(f * daughter, inverse = TRUE)[1:n] / N
          }

          period <- scales * (4 * pi) / (AV_omega0 + sqrt(2 + AV_omega0^2))

          # Variance-normalised AND scale-rectified power.
          #
          # The raw wavelet power spectrum |W|^2 is biased towards long periods:
          # at a given scale the wavelet spans more of the series, so power
          # accumulates with scale even for white noise. Liu et al. (2007) show
          # that dividing by the scale removes this bias, and that the rectified
          # spectrum is the one whose peaks correspond to the periodicities that
          # actually carry variance. WaverideR, used for the published figures,
          # rectifies; this implementation now does the same, so that the two
          # agree. Rectification is a per-scale constant, so it cancels in the
          # pointwise tests and affects the band tests only through the relative
          # weighting of periods within a band.
          power <- (Mod(W)^2) / (sdx^2 * scales)

          # cone of influence (e-folding time of the Morlet wavelet)
          coi_fac <- (4 * pi) / (AV_omega0 + sqrt(2 + AV_omega0^2)) * (1 / sqrt(2))
          edge <- pmin(seq_len(n) - 1, n - seq_len(n))
          edge[edge == 0] <- 1e-5
          coi <- coi_fac * sqrt(2) * dt * edge

          valid <- outer(period, coi, "<=")   # TRUE = outside the COI

          list(power = power, period = period, scales = scales,
              coi = coi, valid = valid, n = n)
        }

        # Global wavelet spectrum: time-average of power, using ONLY cells
        # outside the cone of influence.
        AV_gws <- function(w) {
          out <- numeric(length(w$period))
          for (i in seq_along(w$period)) {
            v <- w$valid[i, ]
            if (!any(v)) {
              out[i] <- NA_real_
            } else {
              out[i] <- mean(w$power[i, v])
            }
          }
          out
        }

        # Per-scale quantile of the power, using ONLY cells outside the cone
        # of influence (NA where a scale has no such cell).
        AV_valid_row_quantile <- function(w, prob) {
          out <- numeric(length(w$period))
          for (i in seq_along(w$period)) {
            v <- w$valid[i, ]
            if (!any(v)) {
              out[i] <- NA_real_
            } else {
              out[i] <- quantile(w$power[i, v], prob, names = FALSE)
            }
          }
          out
        }

        # Row-wise quantiles of a matrix: one row per matrix row, one column
        # per probability (columns named "90%", "95%", ... as by quantile()).
        AV_row_quantiles <- function(M, probs) {
          out <- matrix(NA_real_, nrow = nrow(M), ncol = length(probs))
          for (i in seq_len(nrow(M))) {
            out[i, ] <- quantile(M[i, ], probs = probs, na.rm = TRUE, names = FALSE)
          }
          dimnames(out) <- list(NULL, names(quantile(numeric(0), probs = probs)))
          out
        }

        # Column-wise maximum of a matrix
        AV_column_max <- function(M, na.rm = FALSE) {
          out <- numeric(ncol(M))
          for (j in seq_len(ncol(M))) {
            out[j] <- max(M[, j], na.rm = na.rm)
          }
          out
        }

        # Column-wise all() of a logical matrix
        AV_column_all <- function(M) {
          out <- logical(ncol(M))
          for (j in seq_len(ncol(M))) {
            out[j] <- all(M[, j])
          }
          out
        }


        # Background (null) models
        #
        # All four generators return a series with the same length, mean and
        # variance as the observed chronology, so that only the SHAPE of the
        # background spectrum differs between them.

        AV_ar1_phi <- function(x) {
          as.numeric(arima(x, order = c(1, 0, 0), include.mean = TRUE,
              method = "ML")$coef["ar1"])
        }

        # "Signal-corrected" AR(1)
        #
        # Fitting AR(1) to a series that contains the very cycles that are
        # being tested inflates phi, because a periodic component contributes
        # positive lag-1 autocovariance. The inflated phi then raises the
        # red-noise background at exactly the frequencies of interest, so the
        # test is biased against the signal it is meant to detect.
        #
        # Here the candidate bands are notched out of the Fourier spectrum
        # before phi is estimated. The remaining series is what the manuscript
        # treats as "background", and its phi is the honest background phi.

        AV_notch <- function(x, periods, tol = AV_cand_tol) {

          n <- length(x)
          mu <- mean(x)
          f <- fft(x - mu)
          freq <- (0:(n - 1)) / n

          drop <- rep(FALSE, n)
          for (p in periods) {
            f0 <- 1 / p
            lo <- 1 / (p * (1 + tol))
            hi <- 1 / (p * (1 - tol))
            drop <- drop | (freq >= lo & freq <= hi) |
                (freq >= 1 - hi & freq <= 1 - lo)
          }

          f[drop] <- 0
          Re(fft(f, inverse = TRUE) / n) + mu
        }

        AV_sim_white <- function(n, par) rnorm(n, par$mean, par$sd)

        AV_sim_ar1 <- function(n, par) {
          as.numeric(arima.sim(list(ar = par$phi), n = n,
              sd = par$sd * sqrt(1 - par$phi^2))) + par$mean
        }

        # Power-law (1/f^beta) surrogate with prescribed mean and sd
        AV_sim_powerlaw <- function(n, par) {
          nn <- 2^ceiling(log2(n))
          freq <- (1:(nn / 2)) / nn
          amp <- freq^(-par$beta / 2)
          ph <- runif(nn / 2, 0, 2 * pi)
          spec <- complex(modulus = amp, argument = ph)
          full <- c(0, spec, Conj(rev(spec))[-1])
          y <- Re(fft(full, inverse = TRUE))[1:n]
          par$mean + par$sd * (y - mean(y)) / sd(y)
        }

        # beta estimated by log-log regression of the MTM spectrum,
        # with the candidate bands excluded from the fit
        AV_fit_beta <- function(sp, periods = AV_candidates, tol = AV_cand_tol) {
          keep <- sp$Period >= AV_period_min & sp$Period <= AV_period_max &
              is.finite(sp$Power) & sp$Power > 0
          for (p in periods) {
            keep <- keep & !(sp$Period >= p * (1 - tol) & sp$Period <= p * (1 + tol))
          }
          fit <- lm(log(sp$Power[keep]) ~ log(sp$Frequency[keep]))
          -as.numeric(coef(fit)[2])
        }

        AV_null_models <- function(x) {

          phi_full <- AV_ar1_phi(x)
          x_notch  <- AV_notch(x, AV_candidates)
          phi_corr <- AV_ar1_phi(x_notch)

          list(
              N0_white = list(sim = AV_sim_white,
                  par = list(mean = mean(x), sd = sd(x)),
                  label = "white noise"),
              N1_ar1_full = list(sim = AV_sim_ar1,
                  par = list(mean = mean(x), sd = sd(x), phi = phi_full),
                  label = sprintf("AR(1), phi = %.3f (observed series)",
                      phi_full)),
              N2_ar1_corr = list(sim = AV_sim_ar1,
                  par = list(mean = mean(x), sd = sd(x), phi = phi_corr),
                  label = sprintf("AR(1), phi = %.3f (signal-corrected)",
                      phi_corr))
          )
        }


        # MTM spectrum

        AV_mtm <- function(x, nfft) {

          fit <- multitaper::spec.mtm(
              timeSeries = as.numeric(x), nw = AV_mtm_nw, k = AV_mtm_k,
              nFFT = nfft, returnZeroFreq = FALSE, Ftest = FALSE, jackknife = FALSE,
              adaptiveWeighting = TRUE, plot = FALSE, deltat = 1, dtUnits = "year")

          data.frame(Frequency = as.numeric(fit$freq),
              Period = 1 / as.numeric(fit$freq),
              Power = as.numeric(fit$spec))
        }


        # Significance tests
        #
        # Three levels, reported side by side:
        #
        #   pointwise    P(background power at THIS period >= observed)
        #   max-in-band  P(max background power ANYWHERE in the band >= observed
        #                max) - controls for having searched the band
        #   area-wise    P(a contiguous patch of pointwise-significant cells of
        #                THIS size occurs in the background) - controls for the
        #                fact that red noise produces isolated significant
        #                patches by construction (Maraun & Kurths 2004)

        AV_band_index <- function(period, pmin, pmax, last = FALSE) {
          if (last) {
            which(period >= pmin & period <= pmax)
          } else {
            which(period >= pmin & period < pmax)
          }
        }

        AV_spectral_test <- function(obs, period, null, n, transform,
            nsim = AV_nsim_spectral, quiet = FALSE) {

          # obs       : observed spectrum (vector, aligned with `period`)
          # transform : function(series) -> spectrum vector on the same grid
          sims <- matrix(NA_real_, nrow = length(period), ncol = nsim)

          for (s in seq_len(nsim)) {
            sims[, s] <- transform(null$sim(n, null$par))
            if (!quiet && s %% 250 == 0) {
              cat("   surrogate", s, "of", nsim, "\n")
            }
          }

          # pointwise p and envelopes
          p_point <- rowMeans(sims >= obs, na.rm = TRUE)
          ci <- AV_row_quantiles(sims, probs = c(0.90, 0.95, 0.99))

          # max-in-band p
          nb <- nrow(AV_bands)
          band_rows <- vector("list", nb)
          for (b in seq_len(nb)) {
            band_rows[[b]] <- AV_band_index(period, AV_bands$Period_min[b],
                AV_bands$Period_max[b],
                last = b == nb)
          }

          band_list <- vector("list", nb)
          for (b in seq_len(nb)) {
            rr <- band_rows[[b]]
            if (!length(rr)) {
              next
            }
            obs_max <- max(obs[rr], na.rm = TRUE)
            sim_max <- AV_column_max(sims[rr, , drop = FALSE], na.rm = TRUE)
            band_list[[b]] <- data.frame(Band = AV_bands$Band[b],
                Peak_period = period[rr][which.max(obs[rr])],
                Peak_power = obs_max,
                p_max_in_band = mean(sim_max >= obs_max),
                stringsAsFactors = FALSE)
          }
          band_tab <- do.call(rbind, band_list)

          # candidate-period test: max inside a +/-10% period window
          cand_list <- vector("list", length(AV_candidates))
          for (j_cand in seq_along(AV_candidates)) {
            p <- AV_candidates[j_cand]
            rr <- which(period >= p * (1 - AV_cand_tol) &
                period <= p * (1 + AV_cand_tol))
            if (!length(rr)) {
              next
            }
            obs_max <- max(obs[rr], na.rm = TRUE)
            sim_max <- AV_column_max(sims[rr, , drop = FALSE], na.rm = TRUE)
            cand_list[[j_cand]] <- data.frame(Candidate = p,
                Peak_period = period[rr][which.max(obs[rr])],
                Peak_power = obs_max,
                p_pointwise = min(p_point[rr], na.rm = TRUE),
                p_window = mean(sim_max >= obs_max),
                stringsAsFactors = FALSE)
          }
          cand_tab <- do.call(rbind, cand_list)

          list(p_point = p_point, ci = ci, band = band_tab, cand = cand_tab)
        }


        # Area-wise (cluster) test on the local wavelet field

        AV_clusters <- function(mask) {

          # 4-connected labelling of a logical matrix, returns patch sizes
          lab <- matrix(0L, nrow(mask), ncol(mask))
          cur <- 0L
          sizes <- integer(0)

          for (j in seq_len(ncol(mask))) {
            for (i in seq_len(nrow(mask))) {
              if (!mask[i, j] || lab[i, j] != 0L) {
                next
              }
              cur <- cur + 1L
              stack <- matrix(c(i, j), ncol = 2)
              size <- 0L
              while (nrow(stack)) {
                pt <- stack[nrow(stack), ]
                stack <- stack[-nrow(stack), , drop = FALSE]
                a <- pt[1]
                b <- pt[2]
                if (a < 1 || b < 1 || a > nrow(mask) || b > ncol(mask)) {
                  next
                }
                if (!mask[a, b] || lab[a, b] != 0L) {
                  next
                }
                lab[a, b] <- cur
                size <- size + 1L
                stack <- rbind(stack, c(a - 1, b), c(a + 1, b),
                    c(a, b - 1), c(a, b + 1))
              }
              sizes <- c(sizes, size)
            }
          }
          list(label = lab, sizes = sizes)
        }


        # Band-wise persistence (area-wise) test
        #
        # Rationale (Maraun & Kurths 2004): under any red background, isolated
        # patches of pointwise-significant wavelet power occur by construction.
        # What cannot be produced easily by red noise is a LONG, CONTINUOUS
        # stretch of elevated band power. The statistic used here is therefore
        # the longest continuous run of years in which the band-mean wavelet
        # power exceeds the 95th percentile of the background, evaluated only
        # in the years in which the whole band lies outside the cone of
        # influence. The p-value is the fraction of surrogate fields whose
        # longest run equals or exceeds the observed one.

        AV_band_mean_series <- function(power, period, valid, b) {

          rr <- AV_band_index(period, AV_bands$Period_min[b], AV_bands$Period_max[b],
              last = b == nrow(AV_bands))
          if (!length(rr)) {
            return(NULL)
          }

          full_valid <- AV_column_all(valid[rr, , drop = FALSE])
          if (!any(full_valid)) {
            return(NULL)
          }

          list(rows = rr,
              cols = which(full_valid),
              value = colMeans(power[rr, full_valid, drop = FALSE]))
        }

        AV_longest_above <- function(v, crit) {
          r <- rle(v > crit)
          i <- which(r$values)
          if (!length(i)) {
            0L
          } else {
            max(r$lengths[i])
          }
        }


        # Main analysis for one series

        AV_analyse <- function(entry, key, do_local = TRUE) {

          x <- as.numeric(entry$x)
          n <- length(x)

          cat("\n\n============================================================\n")
          cat("SERIES:", entry$label, "\n")
          cat("============================================================\n")
          cat("n =", n, "  mean =", round(mean(x), 4),
              "  sd =", round(sd(x), 4), "\n")

          nfft <- 2^ceiling(log2(n))

          # Observed spectra
          mtm_obs <- AV_mtm(x, nfft)
          mtm_obs <- mtm_obs[mtm_obs$Period >= AV_period_min &
              mtm_obs$Period <= AV_period_max, ]

          w <- AV_cwt(x)
          gws_obs <- AV_gws(w)

          # Background models
          nulls <- AV_null_models(x)

          beta <- AV_fit_beta(mtm_obs)
          nulls$N3_powerlaw <- list(
              sim = AV_sim_powerlaw,
              par = list(mean = mean(x), sd = sd(x), beta = beta),
              label = sprintf("power law, beta = %.3f", beta))

          cat("\nBackground models:\n")
          for (nm in names(nulls)) {
            cat("   ", nm, ":", nulls[[nm]]$label, "\n")
          }

          # Global tests under every background
          mtm_transform <- function(y) {
            sp <- AV_mtm(y, nfft)
            sp$Power[sp$Period >= AV_period_min & sp$Period <= AV_period_max]
          }
          gws_transform <- function(y) AV_gws(AV_cwt(y))

          tests <- list()

          for (nm in names(nulls)) {
            cat("\n   MTM vs", nm, "\n")
            r_mtm <- AV_spectral_test(mtm_obs$Power, mtm_obs$Period, nulls[[nm]], n,
                mtm_transform, quiet = TRUE)
            cat("   GWS vs", nm, "\n")
            r_gws <- AV_spectral_test(gws_obs, w$period, nulls[[nm]], n,
                gws_transform, quiet = TRUE)
            tests[[nm]] <- list(mtm = r_mtm, gws = r_gws)
          }

          # Local wavelet field
          local <- NULL

          if (do_local) {

            local <- list()

            for (nm in names(nulls)) {

              cat("\n   local wavelet field vs", nm, "\n")

              nsim <- AV_nsim_local
              nper <- length(w$period)

              # scale-wise null distribution of local power
              sim_scale <- matrix(NA_real_, nrow = nper, ncol = nsim)

              # band-wise longest-run statistic under the null
              nb <- nrow(AV_bands)
              sim_run <- matrix(NA_integer_, nrow = nb, ncol = nsim)
              band_def <- vector("list", nb)
              for (b in seq_len(nb)) {
                # AV_band_mean_series() can return NULL; keep the slot
                band_def[b] <- list(AV_band_mean_series(w$power, w$period,
                    w$valid, b))
              }

              sim_band_vals <- vector("list", nb)
              for (b in seq_len(nb)) {
                sim_band_vals[[b]] <- numeric(0)
              }

              sim_store <- vector("list", nsim)

              for (s in seq_len(nsim)) {
                ws <- AV_cwt(nulls[[nm]]$sim(n, nulls[[nm]]$par))
                # per-scale 95th percentile of valid cells
                sim_scale[, s] <- AV_valid_row_quantile(ws, 0.95)

                bv <- vector("list", nb)
                for (b in seq_len(nb)) {
                  bd <- band_def[[b]]
                  if (is.null(bd)) {
                    next
                  }
                  bv[[b]] <- colMeans(ws$power[bd$rows, bd$cols, drop = FALSE])
                }
                sim_store[[s]] <- bv
                if (s %% 200 == 0) {
                  cat("      surrogate", s, "of", nsim, "\n")
                }
              }

              # critical level per scale = mean over surrogates of the
              # scale-wise 95th percentile
              crit_scale <- rowMeans(sim_scale, na.rm = TRUE)

              # band critical levels from the pooled surrogate band series
              band_list <- vector("list", nb)
              for (b in seq_len(nb)) {
                bd <- band_def[[b]]
                if (is.null(bd)) {
                  next
                }
                pooled_parts <- vector("list", nsim)
                for (s in seq_len(nsim)) {
                  pooled_parts[s] <- list(sim_store[[s]][[b]])
                }
                pooled <- unlist(pooled_parts)
                crit <- quantile(pooled, 0.95, names = FALSE)
                obs_run <- AV_longest_above(bd$value, crit)
                sim_runs <- numeric(nsim)
                for (s in seq_len(nsim)) {
                  sim_runs[s] <- AV_longest_above(sim_store[[s]][[b]], crit)
                }
                band_list[[b]] <- data.frame(
                    Band = AV_bands$Band[b],
                    N_years_outside_COI = length(bd$cols),
                    Crit95_band_power = crit,
                    Percent_years_above = 100 * mean(bd$value > crit),
                    Longest_run_years = obs_run,
                    Null_median_run = median(sim_runs),
                    Null_p95_run = quantile(sim_runs, 0.95, names = FALSE),
                    p_areawise = mean(sim_runs >= obs_run),
                    stringsAsFactors = FALSE)
              }
              band_tab <- do.call(rbind, band_list)

              local[[nm]] <- list(crit_scale = crit_scale, band = band_tab)
            }
          }

          list(series = entry, key = key, n = n, mtm = mtm_obs, gws = gws_obs,
              wavelet = w, nulls = nulls, tests = tests, local = local)
        }


        # Run every chronology and window

        AV_results <- list()

        AV_results$std050_full <-
            AV_analyse(AV_series$std050_full, "std050_full", do_local = TRUE)

        AV_results$std050_core <-
            AV_analyse(AV_series$std050_core, "std050_core", do_local = TRUE)

        AV_results$vsc101_core <-
            AV_analyse(AV_series$vsc101_core, "vsc101_core", do_local = FALSE)

        AV_results$constant40_core <-
            AV_analyse(AV_series$constant40_core, "constant40_core", do_local = FALSE)


        # AR(1)-prewhitened sensitivity run
        #
        # Reviewer 1 asks explicitly whether the results survive prewhitening.
        # Prewhitening removes variance at ALL frequencies in proportion to the
        # AR(1) transfer function, so the absolute power of the low-frequency
        # components necessarily drops; the question is whether the PEAKS
        # survive relative to the (now white) background.

        AV_prewhiten <- function(x) {
          phi <- AV_ar1_phi(x)
          r <- as.numeric(x[-1] - phi * x[-length(x)])
          list(x = r, phi = phi)
        }

        pw <- AV_prewhiten(AV_series$std050_full$x)

        AV_series$std050_full_pw <- list(
            x = pw$x,
            year = AV_series$std050_full$year[-1],
            label = sprintf("std050 prewhitened (phi = %.3f)", pw$phi))

        AV_results$std050_full_pw <-
            AV_analyse(AV_series$std050_full_pw, "std050_full_pw", do_local = FALSE)


        # Sample-depth diagnostics

        AV_running_var <- function(x, k = 101) {
          n <- length(x)
          out <- rep(NA_real_, n)
          h <- (k - 1) / 2
          for (i in seq_len(n)) {
            lo <- i - h
            hi <- i + h
            if (lo < 1 || hi > n) {
              next
            }
            out[i] <- var(x[lo:hi])
          }
          out
        }

        AV_depth <- data.frame(
            Year = AV_diag$Year[AV_full$idx],
            samp_depth = AV_diag$samp_depth[AV_full$idx],
            Rbar101 = AV_diag$Rbar101[AV_full$idx],
            Var_std050 = AV_running_var(AV_diag$std050[AV_full$idx]),
            Var_vsc101 = AV_running_var(AV_diag$vsc101[AV_full$idx])
        )

        AV_depth_cor <- data.frame(
            Comparison = c("samp_depth vs Var(std050)",
                "samp_depth vs Var(vsc101)",
                "samp_depth vs Rbar101",
                "Var(std050) vs Var(vsc101)"),
            r = c(cor(AV_depth$samp_depth, AV_depth$Var_std050, use = "complete.obs"),
                cor(AV_depth$samp_depth, AV_depth$Var_vsc101, use = "complete.obs"),
                cor(AV_depth$samp_depth, AV_depth$Rbar101,    use = "complete.obs"),
                cor(AV_depth$Var_std050, AV_depth$Var_vsc101, use = "complete.obs")),
            stringsAsFactors = FALSE)

        # EPS from running Rbar and sample depth
        AV_eps <- with(AV_depth,
            samp_depth * Rbar101 / (1 + (samp_depth - 1) * Rbar101))
        AV_depth$EPS <- AV_eps


        # Do the three chronologies give the same global spectrum?

        AV_gws_compare <- data.frame(
            Period = AV_results$std050_core$wavelet$period,
            std050 = AV_results$std050_core$gws,
            vsc101 = AV_results$vsc101_core$gws,
            constant40 = AV_results$constant40_core$gws
        )

        AV_gws_cor <- data.frame(
            Comparison = c("std050 vs vsc101", "std050 vs constant40",
                "vsc101 vs constant40"),
            r = c(cor(AV_gws_compare$std050, AV_gws_compare$vsc101),
                cor(AV_gws_compare$std050, AV_gws_compare$constant40),
                cor(AV_gws_compare$vsc101, AV_gws_compare$constant40)),
            stringsAsFactors = FALSE)


        # Collect the tables

        AV_collect <- function(what) {
          keys <- names(AV_results)
          per_key <- vector("list", length(keys))
          for (i_key in seq_along(keys)) {
            k <- keys[i_key]
            r <- AV_results[[k]]
            null_names <- names(r$tests)
            per_null <- vector("list", length(null_names))
            for (i_null in seq_along(null_names)) {
              nm <- null_names[i_null]
              methods <- c("mtm", "gws")
              per_meth <- vector("list", length(methods))
              for (i_meth in seq_along(methods)) {
                meth <- methods[i_meth]
                z <- r$tests[[nm]][[meth]][[what]]
                if (is.null(z)) {
                  next
                }
                z$Series <- r$series$label
                z$Key <- k
                z$Null <- nm
                z$Null_label <- r$nulls[[nm]]$label
                z$Method <- toupper(meth)
                per_meth[[i_meth]] <- z
              }
              per_null[i_null] <- list(do.call(rbind, per_meth))
            }
            per_key[i_key] <- list(do.call(rbind, per_null))
          }
          do.call(rbind, per_key)
        }

        AV_candidate_table <- AV_collect("cand")
        AV_band_table <- AV_collect("band")

        AV_local_table <- local({
          keys <- names(AV_results)
          per_key <- vector("list", length(keys))
          for (i_key in seq_along(keys)) {
            k <- keys[i_key]
            r <- AV_results[[k]]
            if (is.null(r$local)) {
              next
            }
            null_names <- names(r$local)
            per_null <- vector("list", length(null_names))
            for (i_null in seq_along(null_names)) {
              nm <- null_names[i_null]
              z <- r$local[[nm]]$band
              if (is.null(z)) {
                next
              }
              z$Series <- r$series$label
              z$Key <- k
              z$Null <- nm
              z$Null_label <- r$nulls[[nm]]$label
              per_null[[i_null]] <- z
            }
            per_key[i_key] <- list(do.call(rbind, per_null))
          }
          do.call(rbind, per_key)
        })

        AV_phi_table <- local({
          keys <- names(AV_results)
          rows <- vector("list", length(keys))
          for (i_key in seq_along(keys)) {
            k <- keys[i_key]
            r <- AV_results[[k]]
            rows[[i_key]] <- data.frame(Key = k, Series = r$series$label, N = r$n,
                phi_observed = r$nulls$N1_ar1_full$par$phi,
                phi_signal_corrected = r$nulls$N2_ar1_corr$par$phi,
                beta_powerlaw = r$nulls$N3_powerlaw$par$beta,
                stringsAsFactors = FALSE)
          }
          do.call(rbind, rows)
        })


        list(candidate = AV_candidate_table, band = AV_band_table, local = AV_local_table,
            phi = AV_phi_table, depth_cor = AV_depth_cor, gws_cor = AV_gws_cor,
            gws_compare = AV_gws_compare, depth = AV_depth)
      })))

  # Table S1: background parameters. Tables S2, S3, S7 and S8 are subsets of the candidate
  # table (S2: Method GWS and S3: Method MTM for std050_full; S7: std050_full_pw;
  # S8: the three core-window chronologies under N1_ar1_full)
  write.csv(suite$phi, in_out("tables", "Table_S1_background_parameters.csv"), row.names = FALSE)
  write.csv(suite$candidate, in_out("tables", "Tables_S2_S3_S7_S8_candidate_window_tests.csv"),
      row.names = FALSE)
  write.csv(suite$band, in_out("tables", "Table_S4_octave_band_tests.csv"), row.names = FALSE)
  write.csv(suite$local, in_out("tables", "Table_S6_duration_test_omega6.csv"), row.names = FALSE)
  write.csv(suite$depth_cor, in_out("tables", "S1_sample_depth_correlations.csv"), row.names = FALSE)
  write.csv(suite$gws_cor, in_out("tables", "S5_global_spectrum_correlations.csv"), row.names = FALSE)
  write.csv(suite$depth, in_out("tables", "S1_sample_depth_diagnostics.csv"), row.names = FALSE)

  cat("\nTable S1: background parameters\n")
  print(suite$phi, row.names = FALSE)
  cat("\nTable S6: duration-based persistence test at omega0 = 6, std050 full window\n")
  print(subset(suite$local, Key == "std050_full",
      select = c(Band, Null, Longest_run_years, Null_p95_run, p_areawise)), row.names = FALSE)
}


# S6: maximum-in-band tests of the detrending variants (optional) ####
# Needs data/AV_objects_07.rds, the five detrending variants built from the raw ring widths
# by J. Bazelmans' step scripts 1-7. The section is skipped when that file is absent.
# Output: output/tables/S6_variants_max_in_band_tests.csv, S6_variants_candidate_window_tests.csv,
#         S6_variants_background_parameters.csv
{
  variants_file <- in_data("AV_objects_07.rds")
  if (!file.exists(variants_file)) {
    cat("\nS6 detrending variants: ", variants_file, " not found, section skipped\n", sep = "")
  } else {
    variants <- run_or_load(paste0("S6_variants_n", n_surrogates_variants),
        with_seed(20092026, local({
          AV_nsim_band <- n_surrogates_variants

          AV_mtm_nw <- 4
          AV_mtm_k  <- 7

          AV_period_min <- 2
          AV_period_max <- 128

          AV_omega0 <- 6             # Morlet, as used by WaveletComp

          # Candidate periodicities reported in the manuscript; each is tested
          # inside a +/-10% window in PERIOD.
          AV_candidates <- c(3.6, 5.9, 8, 12, 15.4, 24, 38, 60, 80, 110)
          AV_cand_tol   <- 0.10


          AV07 <- readRDS(variants_file)

          AV_VSC_core <- AV07$AV_VSC_core
          AV_bands    <- AV07$AV_bands

          chronology_names <- c("AVchron", "Spline67", "Spline200",
              "Spline100", "ModNegExp")

          stopifnot(all(chronology_names %in% names(AV_VSC_core)))

          cat("\nYears:", nrow(AV_VSC_core),
              " (", min(AV_VSC_core$Year), "to", max(AV_VSC_core$Year), ")\n")


          # Morlet transform with an explicit cone of influence
          #
          # Implemented here rather than taken from WaveletComp so that the same
          # transform can be applied to observed and surrogate series and so that
          # the COI mask is available as a matrix.

          AV_cwt <- function(x, dt = 1, dj = 1 / 20) {

            n <- length(x)
            x <- x - mean(x)
            sdx <- sd(x)

            fac <- (4 * pi) / (AV_omega0 + sqrt(2 + AV_omega0^2))
            s0 <- AV_period_min * dt / fac
            smax <- AV_period_max * dt / fac
            scales <- s0 * 2^((0:floor(log2(smax / s0) / dj)) * dj)

            N <- 2^ceiling(log2(n))
            xpad <- c(x, rep(0, N - n))

            k <- (1:floor(N / 2)) * (2 * pi) / (N * dt)
            k <- c(0, k, -k[(floor((N - 1) / 2)):1])
            f <- fft(xpad)

            power <- matrix(NA_real_, nrow = length(scales), ncol = n)

            for (i in seq_along(scales)) {
              s <- scales[i]
              daughter <- sqrt(2 * pi * s / dt) * pi^(-0.25) *
                  exp(-((s * k - AV_omega0)^2) / 2 * (k > 0)) * (k > 0)
              w <- fft(f * daughter, inverse = TRUE)[1:n] / N
              # scale-rectified power (Liu et al., 2007); see step 9
              power[i, ] <- (Mod(w)^2) / (sdx^2 * s)
            }

            period <- scales * fac
            edge <- pmin(seq_len(n) - 1, n - seq_len(n))
            edge[edge == 0] <- 1e-5
            coi <- fac / sqrt(2) * sqrt(2) * dt * edge

            list(power = power, period = period, coi = coi,
                valid = outer(period, coi, "<="), n = n)
          }

          # Global wavelet spectrum over non-COI cells only
          AV_gws <- function(w) {
            out <- numeric(length(w$period))
            for (i in seq_along(w$period)) {
              v <- w$valid[i, ]
              if (!any(v)) {
                out[i] <- NA_real_
              } else {
                out[i] <- mean(w$power[i, v])
              }
            }
            out
          }

          # Row-wise quantile of a matrix (one value per row, NA values dropped)
          AV_row_quantile <- function(M, prob) {
            out <- numeric(nrow(M))
            for (i in seq_len(nrow(M))) {
              out[i] <- quantile(M[i, ], probs = prob, na.rm = TRUE, names = FALSE)
            }
            out
          }

          # Column-wise maximum of a matrix
          AV_column_max <- function(M, na.rm = FALSE) {
            out <- numeric(ncol(M))
            for (j in seq_len(ncol(M))) {
              out[j] <- max(M[, j], na.rm = na.rm)
            }
            out
          }


          # MTM spectrum

          AV_mtm <- function(x, nfft) {
            fit <- multitaper::spec.mtm(
                timeSeries = as.numeric(x), nw = AV_mtm_nw, k = AV_mtm_k,
                nFFT = nfft, returnZeroFreq = FALSE, Ftest = FALSE, jackknife = FALSE,
                adaptiveWeighting = TRUE, plot = FALSE, deltat = 1, dtUnits = "year")
            data.frame(Period = 1 / as.numeric(fit$freq),
                Frequency = as.numeric(fit$freq),
                Power = as.numeric(fit$spec))
          }


          # Background models

          AV_ar1_phi <- function(x) {
            as.numeric(arima(x, order = c(1, 0, 0), include.mean = TRUE,
                method = "ML")$coef["ar1"])
          }

          # Remove the candidate bands from the spectrum before estimating phi.
          AV_notch <- function(x, periods, tol = AV_cand_tol) {
            n <- length(x)
            mu <- mean(x)
            f <- fft(x - mu)
            freq <- (0:(n - 1)) / n
            drop <- rep(FALSE, n)
            for (p in periods) {
              lo <- 1 / (p * (1 + tol))
              hi <- 1 / (p * (1 - tol))
              drop <- drop | (freq >= lo & freq <= hi) |
                  (freq >= 1 - hi & freq <= 1 - lo)
            }
            f[drop] <- 0
            Re(fft(f, inverse = TRUE) / n) + mu
          }

          AV_sim_white <- function(n, par) rnorm(n, par$mean, par$sd)

          AV_sim_ar1 <- function(n, par) {
            as.numeric(arima.sim(list(ar = par$phi), n = n,
                sd = par$sd * sqrt(1 - par$phi^2))) + par$mean
          }

          AV_sim_powerlaw <- function(n, par) {
            nn <- 2^ceiling(log2(n))
            freq <- (1:(nn / 2)) / nn
            spec <- complex(modulus = freq^(-par$beta / 2),
                argument = runif(nn / 2, 0, 2 * pi))
            y <- Re(fft(c(0, spec, Conj(rev(spec))[-1]), inverse = TRUE))[1:n]
            par$mean + par$sd * (y - mean(y)) / sd(y)
          }

          AV_fit_beta <- function(sp) {
            keep <- sp$Period >= AV_period_min & sp$Period <= AV_period_max &
                is.finite(sp$Power) & sp$Power > 0
            for (p in AV_candidates) {
              keep <- keep & !(sp$Period >= p * (1 - AV_cand_tol) &
                  sp$Period <= p * (1 + AV_cand_tol))
            }
            -as.numeric(coef(lm(log(sp$Power[keep]) ~ log(sp$Frequency[keep])))[2])
          }

          AV_backgrounds <- function(x, mtm_obs) {

            phi_obs  <- AV_ar1_phi(x)
            phi_corr <- AV_ar1_phi(AV_notch(x, AV_candidates))
            beta     <- AV_fit_beta(mtm_obs)

            list(
                N0_white = list(sim = AV_sim_white,
                    par = list(mean = mean(x), sd = sd(x)),
                    label = "white noise"),
                N1_ar1_observed = list(sim = AV_sim_ar1,
                    par = list(mean = mean(x), sd = sd(x),
                        phi = phi_obs),
                    label = sprintf("AR(1) phi = %.3f (observed)",
                        phi_obs)),
                N2_ar1_corrected = list(sim = AV_sim_ar1,
                    par = list(mean = mean(x), sd = sd(x),
                        phi = phi_corr),
                    label = sprintf("AR(1) phi = %.3f (signal-corrected)",
                        phi_corr)),
                N3_powerlaw = list(sim = AV_sim_powerlaw,
                    par = list(mean = mean(x), sd = sd(x), beta = beta),
                    label = sprintf("power law beta = %.3f", beta))
            )
          }


          # The maximum-in-band test
          #
          # Statistic: the largest power anywhere inside the band (or inside a
          # +/-10% window around a candidate period).
          #
          # p = proportion of surrogates whose band maximum is at least as large.
          # This is a family-wise test over the band: it pays the price of having
          # searched the band for its best frequency, which pointwise testing
          # does not.

          AV_max_in_band <- function(obs, period, background, n, transform,
              nsim = AV_nsim_band, label = "") {

            sims <- matrix(NA_real_, nrow = length(period), ncol = nsim)

            for (s in seq_len(nsim)) {
              sims[, s] <- transform(background$sim(n, background$par))
              if (s %% 500 == 0) {
                cat("      ", label, s, "/", nsim, "\n")
              }
            }

            p_point <- rowMeans(sims >= obs, na.rm = TRUE)

            nb <- nrow(AV_bands)
            band_list <- vector("list", nb)
            for (b in seq_len(nb)) {
              last <- b == nb
              rr <- if (last) {
                which(period >= AV_bands$Period_min[b] &
                    period <= AV_bands$Period_max[b])
              } else {
                which(period >= AV_bands$Period_min[b] &
                    period < AV_bands$Period_max[b])
              }
              if (!length(rr)) {
                next
              }
              obs_max <- max(obs[rr], na.rm = TRUE)
              sim_max <- AV_column_max(sims[rr, , drop = FALSE], na.rm = TRUE)
              band_list[[b]] <- data.frame(Band = AV_bands$Band[b],
                  Peak_period = period[rr][which.max(obs[rr])],
                  Peak_power = obs_max,
                  p_pointwise_at_peak = p_point[rr][which.max(obs[rr])],
                  p_max_in_band = mean(sim_max >= obs_max),
                  stringsAsFactors = FALSE)
            }
            band_tab <- do.call(rbind, band_list)

            cand_list <- vector("list", length(AV_candidates))
            for (j_cand in seq_along(AV_candidates)) {
              p <- AV_candidates[j_cand]
              rr <- which(period >= p * (1 - AV_cand_tol) &
                  period <= p * (1 + AV_cand_tol))
              if (!length(rr)) {
                next
              }
              obs_max <- max(obs[rr], na.rm = TRUE)
              sim_max <- AV_column_max(sims[rr, , drop = FALSE], na.rm = TRUE)
              cand_list[[j_cand]] <- data.frame(Candidate = p,
                  Peak_period = period[rr][which.max(obs[rr])],
                  Peak_power = obs_max,
                  p_pointwise = min(p_point[rr], na.rm = TRUE),
                  p_max_in_window = mean(sim_max >= obs_max),
                  stringsAsFactors = FALSE)
            }
            cand_tab <- do.call(rbind, cand_list)

            list(band = band_tab, cand = cand_tab, p_point = p_point,
                ci95 = AV_row_quantile(sims, 0.95))
          }


          # Run for every chronology and background, MTM and GWS

          band_out <- list()
          cand_out <- list()
          bg_out   <- list()
          k <- 0

          for (nm in chronology_names) {

            x <- as.numeric(AV_VSC_core[[nm]])
            n <- length(x)
            nfft <- 2^ceiling(log2(n))

            cat("\n\n============================================================\n")
            cat("MAX-IN-BAND TEST:", nm, "\n")
            cat("============================================================\n")

            mtm_obs <- AV_mtm(x, nfft)
            mtm_obs <- mtm_obs[mtm_obs$Period >= AV_period_min &
                mtm_obs$Period <= AV_period_max, ]

            w <- AV_cwt(x)
            gws_obs <- AV_gws(w)

            backgrounds <- AV_backgrounds(x, mtm_obs)

            bg_out[[nm]] <- data.frame(
                Chronology = nm,
                phi_observed = backgrounds$N1_ar1_observed$par$phi,
                phi_signal_corrected = backgrounds$N2_ar1_corrected$par$phi,
                beta_powerlaw = backgrounds$N3_powerlaw$par$beta,
                stringsAsFactors = FALSE)

            mtm_transform <- function(y) {
              sp <- AV_mtm(y, nfft)
              sp$Power[sp$Period >= AV_period_min & sp$Period <= AV_period_max]
            }
            gws_transform <- function(y) AV_gws(AV_cwt(y))

            for (bn in names(backgrounds)) {

              cat("\n  background:", backgrounds[[bn]]$label, "\n")

              for (meth in c("MTM", "GWS")) {

                cat("   ", meth, "\n")

                res <- if (meth == "MTM") {
                  AV_max_in_band(mtm_obs$Power, mtm_obs$Period, backgrounds[[bn]],
                      n, mtm_transform, label = paste(nm, bn, meth))
                } else {
                  AV_max_in_band(gws_obs, w$period, backgrounds[[bn]],
                      n, gws_transform, label = paste(nm, bn, meth))
                }

                k <- k + 1
                res$band$Chronology <- nm
                res$band$Background <- bn
                res$band$Background_label <- backgrounds[[bn]]$label
                res$band$Method <- meth
                res$cand$Chronology <- nm
                res$cand$Background <- bn
                res$cand$Background_label <- backgrounds[[bn]]$label
                res$cand$Method <- meth

                band_out[[k]] <- res$band
                cand_out[[k]] <- res$cand
              }
            }
          }

          AV_band_tests <- do.call(rbind, band_out)
          AV_cand_tests <- do.call(rbind, cand_out)
          AV_backgrounds_used <- do.call(rbind, bg_out)
          rownames(AV_band_tests) <- NULL
          rownames(AV_cand_tests) <- NULL
          rownames(AV_backgrounds_used) <- NULL


          # Collect the results

          list(band = AV_band_tests, cand = AV_cand_tests, backgrounds = AV_backgrounds_used)
        })))
    write.csv(variants$band, in_out("tables", "S6_variants_max_in_band_tests.csv"), row.names = FALSE)
    write.csv(variants$cand, in_out("tables", "S6_variants_candidate_window_tests.csv"),
        row.names = FALSE)
    write.csv(variants$backgrounds, in_out("tables", "S6_variants_background_parameters.csv"),
        row.names = FALSE)
  }
}


# Session information ####
# Output: output/tables/sessionInfo_Weesp.txt
{
  capture.output(sessionInfo(),
      file = in_out("tables", "sessionInfo_Weesp.txt"))

  cat("\nComplete. Figures are in output/figures_pdf and output/figures_png, tables in output/tables\n")

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
