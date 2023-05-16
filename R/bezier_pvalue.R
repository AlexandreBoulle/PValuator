#' Determination of the optimal p-value
#' Using the Bezier method to identify the curvature of a graph
#' @importFrom graphics abline lines points title
#' @importFrom stats approx lm loess predict smooth.spline

#' @param x.tangent as a numeric point (abscissa) to determine a tangent
#' @param smooth as the smoothing of the curve
tangent <- function(x.tangent, smooth){
  pred.deriv0 <- predict(smooth, x = x.tangent, deriv = 0)
  pred.deriv1 <- predict(smooth, x = x.tangent, deriv = 1)
  yint <- pred.deriv0$y - (pred.deriv1$y * x.tangent)
  xint <- - yint / pred.deriv1$y
  return(list(pred.deriv0, pred.deriv1, xint, yint))
}

#' @param x1 as the abscissa of point 1
#' @param x2 as the abscissa of point 2
#' @param y1 as the ordinate of point 1
#' @param y2 as the ordinate of point 2
dist.euclid <- function(x1, x2, y1, y2){
  dist <- sqrt((x1 - x2)**2 + (y1 - y2)**2)
  return(dist)
}

#' @param df as a dataframe containing p-values associated to a number (occurrence, percentage, ...)
bezier_treatment <- function(df){

  df$PValue <- as.numeric(gsub(",", ".", df$PValue))

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

  x <- -log10(df$PValue)
  y <- df$Number
  spl <- smooth.spline(y ~ x, spar = 0.5)

  # Check whether smoothing changes the profile
  # If it does, then the maximum value is changed
  ind.diff <- which(diff(spl$y) > 0.5)
  if (length(ind.diff) > 1){
    ind.diff <- ind.diff[ind.diff < 0.5 * length(x)]
  }
  if (length(ind.diff) > 0){
    ind.begin <- length(x) - max(ind.diff)
  }
  if (length(ind.diff) == 0){
    ind.begin <- length(x)
  }
  # Update usable values
  x <- x[1:ind.begin]
  y <- y[1:ind.begin]
  spl <- smooth.spline(y ~ x, spar = 0.5, cv = TRUE)

  # Add points thanks to an approximation
  approximation <- approx(x, y, n = 1000, yleft = min(y), yright = max(y))
  x <- approximation$x
  y <- approximation$y

  # Find points to determine tangents to the curve
  ind.min <- length(spl$y)
  ind.max <- which.max(spl$y)
  x.tangent1 <- spl$x[ind.min]
  x.tangent2 <- spl$x[ind.max]

  # Find points along each tangent
  predictions1 <- tangent(x.tangent1, spl)
  predictions2 <- tangent(x.tangent2, spl)
  y.tangent1 <- predictions1[[4]] + predictions1[[2]]$y * x
  y.tangent2 <- predictions2[[4]] + predictions2[[2]]$y * x

  # Find the intersection of tangents
  # Tangent 1
  xs1 <- c(spl$x[1], spl$x[length(spl$x)])
  ys1 <- c(y.tangent1[1], y.tangent1[length(y.tangent1)])
  # Fit a linear equation that predicts ys using xs
  tangent1 <- lm(ys1 ~ xs1)
  b.t1 <- as.numeric(tangent1$coefficients[1])
  a.t1 <- as.numeric(tangent1$coefficients[2])
  # Tangent 2
  xs2 <- c(spl$x[1], spl$x[length(spl$x)])
  ys2 <- c(y.tangent2[1], y.tangent2[length(y.tangent2)])
  # Fit a linear equation that predicts ys using xs
  tangent2 <- lm(ys2 ~ xs2)
  b.t2 <- as.numeric(tangent2$coefficients[1])
  a.t2 <- as.numeric(tangent2$coefficients[2])
  # Abscissa and ordinate of intersection point
  x.inter <- (b.t2 - b.t1) / (a.t1 - a.t2)
  y.inter <- a.t1 * x.inter + b.t1

  # Middle of tangents
  # Tangent 1
  x.middle1 <- ( x.inter + x.tangent1 ) / 2
  y.middle1 <- ( y.inter + (a.t1 * x.tangent1 + b.t1) ) / 2
  # Tangent 2
  x.middle2 <- ( x.inter + x.tangent2 ) / 2
  y.middle2 <- ( y.inter + (a.t2 * x.tangent2 + b.t2) ) / 2

  # Middle of line binding the middle of each tangent
  x.middle <- ( x.middle1 + x.middle2 ) / 2
  y.middle <- ( y.middle1 + y.middle2 ) / 2

  # Line binding the middle of each tangent
  a.middles <- (y.middle1 - y.middle2) / (x.middle1 - x.middle2)
  b.middles <- y.middle1 - a.middles * x.middle1

  # Line connecting intersection of tangents and the last middle point
  a.final <- (y.middle - y.inter) / (x.middle - x.inter)
  b.final <- y.middle - a.final * x.middle

  # Find intersection between the last line and the curve
  ind.x.interval1 <- which(spl$x > x.inter)
  x.interval.few.points <- spl$x[ind.x.interval1]
  # Increase the precision of detection thanks to the using of more points
  x.interval <- seq(x.interval.few.points[1], x.interval.few.points[length(x.interval.few.points)], 10**-3)
  dist <- dist.euclid(x.interval, x.interval, (a.final * x.interval + b.final), predict(spl, x.interval)$y)
  ind.pval <- which.min(dist)
  pval <- x.interval[ind.pval]

  # Return results
  return(list(x, y, spl, predictions1[[1]], predictions2[[1]], y.tangent1, y.tangent2, x.inter, y.inter,
              x.middle1, y.middle1, x.middle2, y.middle2, a.middles, b.middles, x.middle, y.middle,
              a.final, b.final, pval))
}

#' Plot graphic to show determination of p-value according to the drawing of Bezier curve
#' @param df as a dataframe containing p-values associated to a number (occurrence, percentage, ...)
#' @export
bezier_plot <- function(df){
  # Recovering results in order to show the plot
  results <- bezier_treatment(df)
  x <- results[[1]]
  y <- results[[2]]
  spl <- results[[3]]
  predictions1 <- results[[4]]
  predictions2 <- results[[5]]
  y.tangent1 <- results[[6]]
  y.tangent2 <- results[[7]]
  x.inter <- results[[8]]
  y.inter <- results[[9]]
  x.middle1 <- results[[10]]
  y.middle1 <- results[[11]]
  x.middle2 <- results[[12]]
  y.middle2 <- results[[13]]
  a.middles <- results[[14]]
  b.middles <- results[[15]]
  x.middle <- results[[16]]
  y.middle <- results[[17]]
  a.final <- results[[18]]
  b.final <- results[[19]]
  pval <- results[[20]]

  # The curve representing the number as a function of the p-value
  plot(x, y, main = "Automatic detection of p-value : Bezier curve", xlab = "-log10(p-value)", ylab = "Percentage of different regions")
  # Smoothing of the curve
  lines(spl, col = "blue3", lwd = 3)
  # Point from which the tangent 1 is constructed
  points(predictions1, col = 2, pch = 19)
  # Line of tangent 1
  lines(x, y.tangent1, col = 3)
  # Point from which the tangent 2 is constructed
  points(predictions2, col = 2, pch = 19)
  # Line of tangent 2
  lines(x, y.tangent2, col = 3)
  # Point of intersection between the two tangents
  points(x.inter, y.inter, col = 2, pch = 19)
  # Point in the middle of the tangent 1
  points(x.middle1, y.middle1, col = 3, pch = 19)
  # Point in the middle of the tangent 2
  points(x.middle2, y.middle2, col = 3, pch = 19)
  # Line passing through the two previous middles
  lines(x, (a.middles * x + b.middles), col = 3)
  # Point located in the middle of the segment passing through the two previous middles
  points(x.middle, y.middle, col = 3, pch = 19)
  # Line connecting the last midpoint to the point of intersection of the tangents
  lines(x, (a.final * x + b.final), col = 3)
  # p-value point
  abline(v = pval, col = "red")
  points(pval, predict(spl, pval)$y, col = "red", pch = 19)
}

#' Show p-value thanks to the geometric method detection (Bezier curve)
#' @param df as a dataframe containing p-values associated to a number (occurrence, percentage, ...)
#' @export
bezier_pvalue <- function(df){
  results <- bezier_treatment(df)
  pvalue <- 10**-results[[20]]
  pvalue
}
