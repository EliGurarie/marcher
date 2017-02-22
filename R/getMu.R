#' Obtain mean vector for a range shift process
#' 
#'  Obtain a mean vector for a movement with one (\code{getMu}) or more (\code{getMu.multi}) range shifts.  This function is mainly used within the likelihood of range shift processes, but is also useful for simulating processes.  
#'  
#' @param T vector of times
#' @param p.m mean parameters. A named vector with elements t1, dt, x1, y1, x2, y2, for a single-shift process.  For multiple (n) shifts, the paramaters are numbered: (x1, x2 ... xn), (y1, y2 ... yn), (t1 .. t[n-1]), (dt1 ... dt[n-1]) 
#' @aliases getMu.multi
#' @seealso \code{\link{sim.mouf}}
#' @example ./examples/getMu.example.r

getMu <- function(T, p.m)
{
    M.x <- ifelse(T <=  p.m["t1"], p.m["x1"], 
                  ifelse(T > p.m["t1"]+p.m["dt"], p.m["x2"], 
                         p.m["x1"] + (p.m["x2"] - p.m["x1"])/ p.m["dt"] * (T - p.m["t1"])))
    
    M.y <- ifelse(T <=  p.m["t1"], p.m["y1"], 
                  ifelse(T > p.m["t1"]+p.m["dt"], p.m["y2"], 
                         p.m["y1"] + (p.m["y2"] - p.m["y1"])/ p.m["dt"] * (T - p.m["t1"])))
    return(cbind(x = M.x, y = M.y))
}


getMu.multi <- function(T, p.m){
  
  xs <- p.m[grep("x", names(p.m))] %>% as.vector
  ys <- p.m[grep("y", names(p.m))] %>% as.vector
  dt <- p.m[grep("dt", names(p.m))] %>% as.vector
  
  n.shifts <- length(xs) - 1
  t.starts <- p.m[paste0("t",1:n.shifts)] %>% as.vector
  t.ends <- t.starts + dt
  
  breaks <- c(min(T) - 1,rowMeans(cbind(t.starts[-1], t.ends[-n.shifts])),max(T))
  
  Mu <- matrix(nrow = length(T), ncol=2, dimnames = list(NULL, c("x","y")))
  for(i in 1:n.shifts){
    which <- which(T > breaks[i] & T <= breaks[i+1])
    
    myp.m = c(t1 = t.starts[i], dt = dt[i],
              x1 = xs[i], y1 = ys[i],
              x2 = xs[i+1], y2 = ys[i+1])
    
    Mu[which,] <- getMu(T[which], myp.m)
  }
  return(Mu)
}



getMu.return <- function(T, p.m){
  x1 <- p.m['x1']
  x2 <- p.m['x2']
  x3 <- p.m['x1']
  y1 <- p.m['y1']
  y2 <- p.m['y2']
  y3 <- p.m['y1']
  t1 <- p.m['t1']
  dt1 <- p.m['dt1']
  t2 <- p.m['t2']
  dt2 <- p.m['dt2']
  
  M.x <- ifelse(T <= t1, x1, 
                ifelse(T < t1+dt1,  x1 + (x2 - x1)/dt1 * (T - t1), 
                       ifelse(T < t2, x2, 
                              ifelse(T < t2 + dt2, x2 + (x3 - x2)/dt2 * (T - t2), 
                                     x3))))
  M.y <- ifelse(T <= t1, y1, 
                ifelse(T < t1+dt1,  y1 + (y2 - y1)/dt1 * (T - t1 ), 
                       ifelse(T < t2, y2, 
                              ifelse(T < t2 + dt2, y2 + (y3 - y2)/dt2 * (T - t2), 
                                     y3))))
  return(cbind(x = M.x, y = M.y))
}