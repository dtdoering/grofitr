#'

makePlots <- function(platename,
                      trunctime = 1000,
                      out_dest = "/Users/dtdoering/1_Research/Lab/DATA/phenotyping/plate\ reader/Output/") {
  # makeplots() accepts a single 4-character plate barcode as a string and generates a PDF
  # containing a growth curve for each well from that plate annotated with the
  # lag, growth rate, and saturation point called by GroFit.
  #
  # Example: makeplots("O7ED")
  cat(noquote(paste('Plotting plate ', platename, '...', sep = "")))
  pdf(
    paste(out_dest,
          Sys.time() %>% format("%Y-%m-%d-%H%M"),
          "_",
          platename,
          ".pdf",
          sep = ""
    )
  )
  for (well in seq_along(GroFitResults[[platename]])) {
    times <- GroFitResults[[platename]][[well]]$GroFitResults$gcFit$gcFittedModels[[1]]$raw.time
    od <- GroFitResults[[platename]][[well]]$GroFitResults$gcFit$gcFittedModels[[1]]$raw.data
    temp_df <- data.frame(cbind(times, od))
    colnames(temp_df) = c("Time", "Absorbance")

    Aobs <- summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$A.model
    A.upCI <-  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$ci95.A.model.up
    A.loCI <-  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$ci95.A.model.lo

    muobs <-  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$mu.model
    mu.upCI <-  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$ci95.mu.model.up
    mu.loCI <-  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$ci95.mu.model.lo

    lambdaobs <-  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$lambda.model
    lambda.upCI <-  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$ci95.lambda.model.up
    lambda.loCI <-  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$ci95.lambda.model.lo

    if (is.null(lambdaobs)) {
      lambdaobs <- 0
    }

    # ggplot statements =======================================================
    # growth data points ------------------------------------------------------
    curve <- ggplot(temp_df, aes(x = Time, y = Absorbance))
    curve <- curve + geom_point(pch = 19)

    # saturation point with confidence intervals ------------------------------
    if (!is.null(Aobs) & !is.na(Aobs) & !is.nan(A.upCI) & !is.nan(A.loCI)) {
      curve <- curve +
      geom_hline(
        yintercept = Aobs,
        color = "blue"
      )
    }
    if (!is.nan(A.upCI) & !is.na(A.upCI)) {
      curve <- curve +
      geom_hline(
        yintercept = A.upCI,
        col = "cyan",
        lty = 2
      )
    }
    if (!is.nan(A.loCI) & !is.na(A.loCI)) {
      curve <- curve +
      geom_hline(
        yintercept = A.loCI,
        col = "cyan",
        lty = 2
      )
    }

    # lag time with confidence intervals --------------------------------------
    if (!is.null(lambdaobs) & !is.na(lambdaobs) & !is.nan(lambda.upCI) & !is.nan(lambda.loCI) & 0 < lambdaobs & lambdaobs < trunctime) {
      curve <- curve +
      geom_vline(
        xintercept = lambdaobs,
        color = "green4"
      )
    }
    if (!is.nan(lambda.upCI) & !is.na(lambda.upCI) & 0 < lambda.upCI & lambda.upCI < trunctime) {
      curve <- curve +
      geom_vline(
        xintercept = lambda.upCI,
        col = "green",
        lty = 2
      )
    }
    if (!is.nan(lambda.loCI) & !is.na(lambda.loCI) & 0 < lambda.loCI & lambda.loCI < trunctime) {
      curve <- curve +
      geom_vline(
        xintercept = lambda.loCI,
        col = "green",
        lty = 2
      )
    }

    # max growth rate with confidence intervals -------------------------------
    if (!is.null(muobs) & !is.na(muobs) & !is.nan(mu.upCI) & !is.nan(mu.loCI) & 0 < lambdaobs & lambdaobs < trunctime) {
      curve <- curve +
      geom_abline(
        intercept = -(muobs*lambdaobs),
        slope = muobs,
        col = "red"
      )
    }
    if (!is.nan(mu.upCI) & !is.na(mu.upCI)) {
      curve <- curve +
      geom_abline(
        intercept = -(muobs*lambdaobs),
        slope = mu.upCI,
        col = "orange",
        lty = 2
      )
    }
    if (!is.nan(mu.loCI) & !is.na(mu.loCI) & 0 < lambdaobs & lambdaobs < trunctime) {
      curve <- curve +
      geom_abline(
        intercept = -(muobs*lambdaobs),
        slope = mu.loCI,
        col = "orange",
        lty = 2
      )
    }

    # format plots - axes, max/min, gridlines, etc. ---------------------------
    curve <- curve +
    ylim(0, 2) +
    scale_x_continuous(breaks = seq(0,150,10)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(color = "black")
    ) +
    xlab("Time (h)") +
    ylab("Absorbance") +

    ggtitle(paste(GroFitResults[[platename]][[well]]$Plate, ": ",
      GroFitResults[[platename]][[well]][[colnames(get(PlateNames[[1]]))[1]]],
      " LICM(",
      GroFitResults[[platename]][[well]][[colnames(get(PlateNames[[1]]))[2]]],
      ",",
      GroFitResults[[platename]][[well]][[colnames(get(PlateNames[[1]]))[3]]],
      ")",
      sep = "")
    )
    print(curve)
  }
  dev.off()
  cat(noquote('Done.'), '\n')
}
