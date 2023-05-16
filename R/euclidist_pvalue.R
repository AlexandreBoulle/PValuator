#' Determination of the optimal p-value
#' Using the euclidean distance to identify the furthest point
#' @importFrom graphics abline lines points title
#' @importFrom stats approx lm loess predict smooth.spline


#' @param df as a dataframe containing p-values associated to a number (occurrence, percentage, ...)
euclidist_treatment <- function(df){

  df$PValue <- as.numeric(gsub(",", ".", df$PValue))
  df$PValue <- -log10(df$PValue)

  # Remove the (n - 1) lines associated with 0 when n lines contain the value 0
  ind.0 <- which(df$Number %in% 0)
  # Find the maximum value to stop the representation of the curve at this value
  ind.max.df <- which.max(df$Number)
  if (length(ind.0) > 1 && ind.max.df < dim(df)[1]){
    df <- df[-c(ind.0[1 : (length(ind.0) - 1)], seq((ind.max.df + 1), dim(df)[1], 1)), ]
  }
  if (length(ind.0) > 1 && ind.max.df == dim(df)[1]){
    df <- df[-c(ind.0[1 : (length(ind.0) - 1)]), ]
  }
  if (length(ind.0) <= 1 && ind.max.df < dim(df)[1]){
    df <- df[-seq((ind.max.df + 1), dim(df)[1], 1), ]
  }
  if (length(ind.0) <= 1 && ind.max.df == dim(df)[1]){
    df <- df
  }

  # Recover values from dataframe
  x <- df$PValue
  y <- df$Number
  # Smoothing of the curve representing the number as a function of the p-value
  lis <- loess(y ~ x, span = 0.5)
  # Study the curve between the first value and the value associated to the highest y (ordinate)
  ind.max.fitted <- which.max(lis$fitted)
  x.after.max <- x[1:ind.max.fitted]
  y.after.max <- y[1:ind.max.fitted]
  fitted <- lis$fitted[1:ind.max.fitted]
  # Directing coefficient of the line connecting the extremities of the curve
  a <- (fitted[1] - fitted[length(fitted)]) / (x.after.max[1] - x.after.max[length(x.after.max)])
  # Ordinate at the origin of the line
  b <- fitted[1] - a * x.after.max[1]

  # Directing coefficient of the perpendicular line
  a.perp <- -1 / a
  # Ordinate at the origin of the perpendicular line for each x
  x.smooth <- seq(x.after.max[length(x.after.max)], x.after.max[1], 0.001)
  y.smooth <- predict(lis, x.smooth)
  b.perp <- y.smooth - a.perp * x.smooth

  # Coordinates of the intersection points between the perpendicular lines and the segment
  x.perp <- (b.perp - b) / (a - a.perp)
  y.perp <- a * x.perp + b

  # Vector containing the distances between each point on the curve and the segment
  dist.list <- sqrt( ( x.perp - x.smooth )**2 + ( y.perp - y.smooth )**2 )

  # Value of the maximum distance between the curve and the line
  ind.max <- which.max(dist.list)
  # Value : -log10(p-value)
  pval <- 10**-x.smooth[ind.max]

  # Return results
  return(list(x, y, lis, a, b, x.smooth, y.smooth, a.perp, b.perp, x.perp, y.perp, ind.max, pval))
}

#' Plot graphic to show determination of p-value according to the euclidean distance method
#' @param df as a dataframe containing p-values associated to a number (occurrence, percentage, ...)
#' @export
euclidist_plot <- function(df){
  # Recovering results in order to show the plot
  results <- euclidist_treatment(df)
  x <- results[[1]]
  y <- results[[2]]
  lis <- results[[3]]
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

  # The curve representing the number as a function of the p-value
  plot(x, y, xlab = "-log10(p-value)", ylab = "Percentage of different regions")
  title(main = "Detect the furthest point of the line")
  # Smoothing of the curve
  lines(lis$fitted ~ lis$x, col = "blue3", lwd = 3)
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

#' Show p-value thanks to the distance method
#' @param df as a dataframe containing p-values associated to a number (occurrence, percentage, ...)
#' @export
euclidist_pvalue <- function(df){
  results <- euclidist_treatment(df)
  pvalue <- results[[13]]
  pvalue
}

