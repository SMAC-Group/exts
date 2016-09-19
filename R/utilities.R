#' Cast Simulation Matrix to Data Frame
#'
#' Creates a \code{data.frame} out of a simulation matrix
#' for use with \code{ggplot2}.
#' @param m     A \code{matrix}
#' @param wide  A \code{boolean} indicating if the simulated data is group
#' by row (\code{TRUE}) or column (\code{FALSE}).
#' @return A \code{data.frame} with three variables:
#' \itemize{
#' \item{Round}{Iteration of the Simulation}
#' \item{Draw}{Draw during the iteration of the simulation}
#' \item{Value}{Value of the statistic at round and draw}
#' }
#' @export
#' @examples
#' # Set Seed
#' set.seed(5812)
#'
#' # Generate data
#' m = matrix(rnorm(10), 2,5)
#'
#' # Organize data.frame by row
#' cast_simdf(m)
#'
#' # Organize by column
#' cast_simdf(m, wide = FALSE)
cast_simdf = function(m, wide = TRUE){

  if(!is.matrix(m)){stop("`m` must be a `data.frame`.")}

  n = nrow(m)
  p = ncol(m)

  if(wide){
    Round = seq_len(n); Draw = seq_len(p)
  } else{
    Draw = seq_len(n); Round = seq_len(p)
  }

  # Hacked from as.data.frame(as.table())
  o = data.frame(expand.grid(Round = Round, Draw = Draw), Values = c(m))

  structure(o, class = c("simdf", "data.frame"))
}


#' Plot Simulation Trials
#'
#' Constructs a line graph containing different simulations
#' @param x,object An \code{\link{simdf}} object.
#' @export
#' @rdname plot.simdf
#' @examples
#' # Set Seed
#' set.seed(5812)
#'
#' # Generate data
#' m = matrix(rnorm(10), 2,5)
#'
#' # Organize data.frame by row
#' sim = cast_simdf(m)
#'
#' # Graph Sim
#' plot(sim)
plot.simdf = function(x, ...){
  autoplot.simdf(object = x, ...)
}


#' @export
#' @rdname plot.simdf
autoplot.simdf = function(object, type = "line", ...){
  ggplot(object, aes(x = Draw, y = Values, group = factor(Round),
                     color = factor(Round))) +
    geom_line(size = 1) +
    theme_bw() + labs(
      x = "Draw",
      y = "Values",
      color = "Round")
}
