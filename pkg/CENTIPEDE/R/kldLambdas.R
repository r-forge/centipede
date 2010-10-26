kldLambdas <-
function(lambda1,lambda2){
  sum(lambda1 * (log(lambda1) - log(lambda2)))
}

