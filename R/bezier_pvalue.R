####################
#    PVALUATOR     #
# Alexandre BOULLE #
####################

# Determination of the optimal p-value
# Using the Bezier method to identify the curvature of a graph

#' @importFrom graphics abline lines points title
#' @importFrom stats loess predict var
#' @import caret

#' @param df as a dataframe containing p-values associated to a number (occurrence, percentage, ...)
#' @param value as the fraction of total variance
remove.small.var <- function(df, value){
  # Remove values of the distribution when they are too close
  variance <- c()
  for (i in 1:(dim(df)[1] - 1)){
    variance <- c(variance, var(df$Number[1:(i+1)]))
  }
  # Keep positions of variances
  ind.var <- which(variance > (variance[length(variance)] * value))
  # Keep lines of interest : -1 to recover the previous value and +1 to recover the last value used to calculate the variance
  df <- df[(ind.var[1] - 1):(ind.var[length(ind.var)] + 1), ]
  return(df)
}

#' @param x as the p-values
#' @param y as the occurrence / percentage associated to each p-value
cv.loess <- function(x, y){
  df <- data.frame(PValue = x, Number = y)
  # Define k-fold cross validation method
  ctrl <- trainControl(method = "cv", number = 5)
  grid <- expand.grid(span = seq(0.1, 1.5, len = 30), degree = 1)
  # Perform cross-validation using smoothing spans ranging from 0.1 to 1.5
  model <- train(Number ~ PValue, data = df, method = "gamLoess", tuneGrid = grid, trControl = ctrl)
  # Print results of k-fold cross-validation
  ind.min.RMSE <- which.min(model$results$RMSE)
  min.span <- model$results$span[ind.min.RMSE]
  return(min.span)
}

#' @param x1 as the abscissa of point 1
#' @param x2 as the abscissa of point 2
#' @param y1 as the ordinate of point 1
#' @param y2 as the ordinate of point 2
dist.euclid <- function(x1, x2, y1, y2){
  dist <- sqrt((x1 - x2)**2 + (y1 - y2)**2)
  return(dist)
}

#' @param df a dataframe containing p-values associated to a number (occurrence, percentage, ...)
#' @param cv a character ("yes" or "no") to use or not the cross validation method in order to detect the best value for span parameter
#' @param sp an integer or a float to use as the span parameter in the loess() function
#' @param fraction.var keep lines if variance of p-value occurrences is greater than a fraction of total variance
bezier_treatment <- function(df, cv, sp, fraction.var){

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

  if (cv == "yes"){
    min.span <- cv.loess(x, y)
    smooth.loess <- loess(y ~ x, span = min.span)
  }
  if (cv == "no"){
    smooth.loess <- loess(y ~ x, span = sp)
  }

  # If smoothing changes the profile then the maximum value is changed
  # Check whether all values are positive
  # Negative values indicate a profile changing
  # Tolerance threshold : [-0.5 - 0]
  ind.diff <- which(diff(smooth.loess$fitted) < -0.5)
  if (length(ind.diff) <= 1){
    ind.begin <- length(x)
  }
  if (length(ind.diff) > 1){
    # Keep the indices if they are smallest than the middle position
    # The goal is to avoid removing a large part of data
    ind.diff <- ind.diff[ind.diff > 0.5 * length(x)]
    if (length(ind.diff) > 0){
      # ind.begin <- length(x) - max(ind.diff)
      ind.begin <- min(ind.diff)
    }
    if (length(ind.diff) == 0){
      ind.begin <- length(x)
    }
  }

  # Update usable values
  x <- x[1:ind.begin]
  y <- y[1:ind.begin]

  if (cv == "yes"){
    min.span <- cv.loess(x, y)
    smooth.loess <- loess(y ~ x, span = min.span)
  }
  if (cv == "no"){
    smooth.loess <- loess(y ~ x, span = sp)
  }
  # smooth.loess <- loess(y ~ x, span = 0.5)

  # Find points to determine tangents to the curve
  # ind.min <- length(spl$y)
  # ind.max <- which.max(spl$y)
  # x.tangent1 <- spl$x[ind.min]
  # x.tangent2 <- spl$x[ind.max]

  # Find the intersection of tangents
  # Tangent 1 equation
  a.t1 <- (smooth.loess$fitted[1] - smooth.loess$fitted[3]) / (smooth.loess$x[1] - smooth.loess$x[3])
  b.t1 <- smooth.loess$fitted[2] - a.t1 * smooth.loess$x[2]
  # Tangent 2 equation
  ind.t2 <- length(smooth.loess$x)
  a.t2 <- (smooth.loess$fitted[ind.t2] - smooth.loess$fitted[ind.t2 - 2]) / (smooth.loess$x[ind.t2] - smooth.loess$x[ind.t2 - 2])
  b.t2 <- smooth.loess$fitted[ind.t2 - 1] - a.t2 * smooth.loess$x[ind.t2 - 1]
  # Abscissa and ordinate of intersection point
  x.inter <- (b.t2 - b.t1) / (a.t1 - a.t2)
  y.inter <- a.t1 * x.inter + b.t1

  # Middle of tangents
  # Tangent 1
  x.tangent1 <- smooth.loess$x[2]
  # x.tangent1 <- smooth.loess$x[1] - 0.002
  x.middle1 <- ( x.inter + x.tangent1 ) / 2
  y.middle1 <- ( y.inter + (a.t1 * x.tangent1 + b.t1) ) / 2
  # Tangent 2
  x.tangent2 <- smooth.loess$x[ind.t2 - 1]
  # x.tangent2 <- smooth.loess$x[ind.t2] + 0.001
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
  # Increase the precision of detection thanks to the using of more points
  x.interval <- seq(x[length(x)], x[1], 0.001)
  y.interval <- predict(smooth.loess, x.interval)
  dist <- dist.euclid(x.interval, x.interval, (a.final * x.interval + b.final), y.interval)
  ind.pval <- which.min(dist)
  pval <- x.interval[ind.pval]

  # Return results
  y.tangent1 <- a.t1 * x.interval + b.t1
  y.tangent2 <- a.t2 * x.interval + b.t2
  return(list(x, y, smooth.loess, list(smooth.loess$x[2], smooth.loess$fitted[2]), list(smooth.loess$x[ind.t2 - 1], smooth.loess$fitted[ind.t2 - 1]),
              x.interval, y.tangent1, y.tangent2, x.inter, y.inter,
              x.middle1, y.middle1, x.middle2, y.middle2, a.middles, b.middles, x.middle, y.middle,
              a.final, b.final, pval))
}

#' @title Graphic for p-value detection (Bezier curve)
#' @description Plot graphic to show determination of p-value according to the geometric method (Bezier curve)
#' @param df a dataframe containing p-values associated to a number (occurrence, percentage, ...)
#' @param cv a character ("yes" or "no") to use or not the cross validation method in order to detect the best value for span parameter
#' @param sp an integer or a float to use as the span parameter in the loess() function
#' @param fraction.var keep lines if variance of p-value occurrences is greater than a fraction of total variance
#' @param main.title a title for the graph
#' @param x.title a title for abscissa axis
#' @param y.title a title for ordinate axis
#' @export
bezier_plot <- function(df, cv = "no", sp = 0.5, fraction.var = 5e-3, main.title = "Automatic detection of p-value : Bezier curve", x.title = "-log10(p-value)", y.title = "Percentage of different regions"){
  # Recovering results in order to show the plot
  results <- bezier_treatment(df, cv, sp, fraction.var)
  x <- results[[1]]
  y <- results[[2]]
  smooth.loess <- results[[3]]
  predictions1 <- results[[4]]
  predictions2 <- results[[5]]
  x.interval <- results[[6]]
  y.tangent1 <- results[[7]]
  y.tangent2 <- results[[8]]
  x.inter <- results[[9]]
  y.inter <- results[[10]]
  x.middle1 <- results[[11]]
  y.middle1 <- results[[12]]
  x.middle2 <- results[[13]]
  y.middle2 <- results[[14]]
  a.middles <- results[[15]]
  b.middles <- results[[16]]
  x.middle <- results[[17]]
  y.middle <- results[[18]]
  a.final <- results[[19]]
  b.final <- results[[20]]
  pval <- results[[21]]

  # The curve representing the number as a function of the p-value
  plot(x, y, main = main.title, xlab = x.title, ylab = y.title)
  # Smoothing of the curve
  # lines(spl, col = "blue3", lwd = 3)
  lines(smooth.loess$fitted ~ smooth.loess$x, col = "blue3", lwd = 3)
  # Point from which the tangent 1 is constructed
  points(predictions1[1], predictions1[2], col = 2, pch = 19)
  # Line of tangent 1
  lines(x.interval, y.tangent1, col = 3)
  # Point from which the tangent 2 is constructed
  points(predictions2[1], predictions2[2], col = 2, pch = 19)
  # Line of tangent 2
  lines(x.interval, y.tangent2, col = 3)
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
  # points(pval, predict(spl, pval)$y, col = "red", pch = 19)
  points(pval, predict(smooth.loess, pval), col = "red", pch = 19)
}

#' @title p-value detection (Bezier curve)
#' @description Determine p-value thanks to the geometric method (Bezier curve)
#' @param df a dataframe containing p-values associated to a number (occurrence, percentage, ...)
#' @param cv a character ("yes" or "no") to use or not the cross validation method in order to detect the best value for span parameter
#' @param sp an integer or a float to use as the span parameter in the loess() function
#' @param fraction.var keep lines if variance of p-value occurrences is greater than a fraction of total variance
#' @export
bezier_pvalue <- function(df, cv = "no", sp = 0.5, fraction.var = 5e-3){
  results <- bezier_treatment(df, cv, sp, fraction.var)
  pvalue <- 10**-results[[21]]
  pvalue
}

