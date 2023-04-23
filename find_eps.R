find_eps <- function(r_star, p_grid = seq(0.001, 1, 0.001), 
                     q_grid = seq(0.001, 1, 0.001)){
  min_eps <- Inf
  for(p in p_grid){
    for(q in q_grid){
      if(q == 1){ eps <- log((1-p)/(1/r_star - p)) }
      else{ eps <- log(2*p*(1-q)/(sqrt((1-p)^2 + 4*p*(1-q)*(1/r_star - p*q)) - (1-p)))}
      if(min_eps > eps){ min_eps <- eps }
    }
  }
  return(min_eps)
}
