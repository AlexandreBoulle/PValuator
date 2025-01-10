####################
#    PVALUATOR     #
# Alexandre BOULLE #
####################

# Determination of the optimal p-value
# Using the euclidean distance to identify the furthest point

#' @importFrom graphics abline lines points title
#' @importFrom stats loess predict var
#' @import caret

#' @param df a dataframe containing p-values associated to a number (occurrence, percentage, ...)
#' @param cv a character ("yes" or "no") to use or not the cross validation method in order to detect the best value for span parameter
#' @param sp an integer or a float to use as the span parameter in the loess() function
#' @param fraction.var keep lines if variance of p-value occurrences is greater than a fraction of total variance
euclidist_treatment <- function(df, cv, sp, fraction.var){

  colnames(df) <- c("PValue", "Number")
  df$PValue <- as.numeric(gsub(",", ".", df$PValue))
  df <- df[order(df$PValue), ]
  df$PValue <- -log10(df$PValue)

  # Remove the lines associated with 0 when n lines contain the value 0
  ind.0 <- which(df$Number %in% 0)
  if (length(ind.0) > 0){
    df <- df[-ind.0, ]
  }
  # Find the maximum value to stop the representation of the curve at this value
  ind.max.df <- which.max(df$Number)
  if (ind.max.df < dim(df)[1]){
    df <- df[-seq((ind.max.df + 1), dim(df)[1], 1), ]
  }

  # Remove the (n - 1) lines associated with 0 when n lines contain the value 0
  # ind.0 <- which(df$Number %in% 0)
  # Find the maximum value to stop the representation of the curve at this value
  # ind.max.df <- which.max(df$Number)
  # if (length(ind.0) > 1 && ind.max.df < dim(df)[1]){
  #   df <- df[-c(ind.0[1 : (length(ind.0) - 1)], seq((ind.max.df + 1), dim(df)[1], 1)), ]
  # }
  # if (length(ind.0) > 1 && ind.max.df == dim(df)[1]){
  #   df <- df[-c(ind.0[1 : (length(ind.0) - 1)]), ]
  # }
  # if (length(ind.0) <= 1 && ind.max.df < dim(df)[1]){
  #   df <- df[-seq((ind.max.df + 1), dim(df)[1], 1), ]
  # }
  # if (length(ind.0) <= 1 && ind.max.df == dim(df)[1]){
  #   df <- df
  # }

  df <- remove.small.var(df, value = fraction.var)
  # Recover values from dataframe
  x <- df$PValue
  y <- df$Number

  # Smoothing of the curve representing the number as a function of the p-value
  if (cv == "yes"){
    min.span <- cv.loess(x, y)
    smooth.loess <- loess(y ~ x, span = min.span)
  }
  if (cv == "no"){
    smooth.loess <- loess(y ~ x, span = sp)
  }
  # lis <- loess(y ~ x, span = 0.5)

  # Study the curve between the first value and the value associated to the highest y (ordinate)
  ind.max.fitted <- which.max(smooth.loess$fitted)
  x.after.max <- x[1:ind.max.fitted]
  y.after.max <- y[1:ind.max.fitted]
  fitted <- smooth.loess$fitted[1:ind.max.fitted]
  # Directing coefficient of the line connecting the extremities of the curve
  a <- (fitted[1] - fitted[length(fitted)]) / (x.after.max[1] - x.after.max[length(x.after.max)])
  # Ordinate at the origin of the line
  b <- fitted[1] - a * x.after.max[1]

  # Directing coefficient of the perpendicular line
  a.perp <- -1 / a
  # Ordinate at the origin of the perpendicular line for each x
  x.smooth <- seq(x.after.max[length(x.after.max)], x.after.max[1], 0.001)
  y.smooth <- predict(smooth.loess, x.smooth)
  b.perp <- y.smooth - a.perp * x.smooth

  # Coordinates of the intersection points between the perpendicular lines and the segment
  x.perp <- (b.perp - b) / (a - a.perp)
  y.perp <- a * x.perp + b

  # Vector containing the distances between each point on the curve and the segment
  # dist.list <- sqrt( ( x.perp - x.smooth )**2 + ( y.perp - y.smooth )**2 )
  dist.list <- dist.euclid(x.perp, x.smooth, y.perp, y.smooth)

  # Value of the maximum distance between the curve and the line
  ind.max <- which.max(dist.list)
  # Value : -log10(p-value)
  pval <- 10**-x.smooth[ind.max]

  # Return results
  return(list(x, y, smooth.loess, a, b, x.smooth, y.smooth, a.perp, b.perp, x.perp, y.perp, ind.max, pval))
}

#' @title Graphic for p-value detection (Euclidean distance)
#' @description Plot graphic to show determination of p-value according to the euclidean distance method
#' @param df a character containing p-values associated to a number (occurrence, percentage, ...)
#' @param cv a string ("yes" or "no") to use or not the cross validation method in order to detect the best value for span parameter
#' @param sp an integer or a float to use as the span parameter in the loess() function
#' @param fraction.var keep lines if variance of p-value occurrences is greater than a fraction of total variance
#' @param main.title a title for the graph
#' @param x.title a title for abscissa axis
#' @param y.title a title for ordinate axis
#' @export
euclidist_plot <- function(df, cv = "no", sp = 0.5, fraction.var = 5e-3, main.title = "Detect the furthest point of the line", x.title = "-log10(p-value)", y.title = "Number / Percentage"){
  # Recovering results in order to show the plot
  results <- euclidist_treatment(df, cv, sp, fraction.var)
  x <- results[[1]]
  y <- results[[2]]
  smooth.loess <- results[[3]]
  a <- results[[4]]
  b <- results[[5]]
  x.smooth <- results[[6]]
  y.smooth <- results[[7]]
  a.perp <- results[[8]]
  b.perp <- results[[9]]
  x.perp <- results[[10]]
  y.perp <- results[[11]]
  ind.max <- results[[12]]
  pval <- results[[13]]

  # The curve representing the number as a function of the p-value (asp = numeric, giving the aspect ratio y/x)
  plot(x, y, xlab = x.title, ylab = y.title)#, asp = 1)
  title(main = main.title)
  # Smoothing of the curve
  lines(smooth.loess$fitted ~ smooth.loess$x, col = "blue3", lwd = 3)
  lines(x, a * x + b, col = "green")
  # Line perpendicular to the segment and passing through the furthest point of the segment
  lines(x.smooth, (a.perp * x.smooth + b.perp[ind.max]), col = "blue")
  # Intersection between the p-value curve and the line perpendicular to the segment
  points(x.smooth[ind.max], y.smooth[ind.max], pch = 19, col = 'red')
  # Intersection between the perpendicular line and the segment
  points(x.perp[ind.max], y.perp[ind.max], pch = 19, col = 'red')
  # Value of the p-value corresponding to the maximum distance between the curve and the segment
  abline(v = x.smooth[ind.max], col = "red")
}

#' @title p-value detection (Euclidean distance)
#' @description Determine p-value thanks to the distance method
#' @param df a character containing p-values associated to a number (occurrence, percentage, ...)
#' @param cv a string ("yes" or "no") to use or not the cross validation method in order to detect the best value for span parameter
#' @param sp an integer or a float to use as the span parameter in the loess() function
#' @param fraction.var keep lines if variance of p-value occurrences is greater than a fraction of total variance
#' @export
euclidist_pvalue <- function(df, cv = "no", sp = 0.5, fraction.var = 5e-3){
  results <- euclidist_treatment(df, cv, sp, fraction.var)
  pvalue <- results[[13]]
  pvalue
}

