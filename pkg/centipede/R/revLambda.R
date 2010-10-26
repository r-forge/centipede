revLambda <-
function(lambda){
  L <- length(lambda)/2;
  NewLambda <- c(rev(lambda[(L+1):(2*L)]),rev(lambda[1:L]));
}

