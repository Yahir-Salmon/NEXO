rmse <- function(y,p) {stopifnot(length(y)==length(p),all(is.finite(y)),all(is.finite(p)));sqrt(mean((y-p)^2))}
mae <- function(y,p) {stopifnot(length(y)==length(p),all(is.finite(y)),all(is.finite(p)));mean(abs(y-p))}
