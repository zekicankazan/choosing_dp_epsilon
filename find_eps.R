find_eps <- function(r_star, p_grid = seq(1e-3, 1, 1e-3), 
                     q_grid = seq(1e-3, 1, 1e-3)){
  # Function to search over a grid to find the minimal epsilon satisfying
  # the risk profile r_star
  min_eps <- Inf
  for(p in p_grid){
    for(q in q_grid){
      if(q == 1){ 
        if(p > 1/r_star(p,1)){ next }
        eps <- log((1-p)/(1/r_star(p,1) - p)) 
      }
      else{ 
        if(4*p*(1-q)*(1/r_star(p,q) - p*q) < 0){ next }
        eps <- log(2*p*(1-q)/(sqrt((1-p)^2 + 4*p*(1-q)*(1/r_star(p,q) - p*q)) - (1-p)))
      }
      if(min_eps > eps){ min_eps <- eps; min_p <- p; min_q <- q }
    }
  }
  return(list(eps = min_eps, p = min_p, q = min_q))
}
