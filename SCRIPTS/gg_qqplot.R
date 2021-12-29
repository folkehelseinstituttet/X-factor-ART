#' Create a quantile-quantile plot with ggplot2.
#'
#' Assumptions:
#'   - Expected P values are uniformly distributed.
#'   - Confidence intervals assume independence between tests.
#'     We expect deviations past the confidence intervals if the tests are
#'     not independent.
#'     For example, in a genome-wide association study, the genotype at any
#'     position is correlated to nearby positions. Tests of nearby genotypes
#'     will result in similar test statistics.
#'
#' Adapted by Julia Romanowska to include data grouping
#'
#' @param ps Vector of p-values; if named, the names give groups.
#' @param ci Size of the confidence interval, 95% by default.
#' @param eq.axes Logical: should the axes show equal range? Default: TRUE.
#' @return A ggplot2 plot.
#' @examples
#' library(ggplot2)
#' gg_qqplot(runif(1e2)) + theme_grey(base_size = 24)
gg_qqplot <- function(ps, ci = 0.95, eq.axes = TRUE) {
  n_groups <- 1
  if(!is.null(names(ps))){
    n_groups <- length(unique(names(ps)))
  }
  df <- create_df_qq(ps, ci)
  
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  max.lims <- range( c( df$observed, df$expected ) ) + 0.5
  
  plot.out <- ggplot(df) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2) +
    geom_line(aes(expected, clower), linetype = 2) +
    xlab(log10Pe) +
    ylab(log10Po)
  if( n_groups > 1 ){
    plot.out <- plot.out +
      geom_point(aes(expected, observed, color = group), shape = 16, size = 1.5)
  } else {
    plot.out <- plot.out +
      geom_point(aes(expected, observed), shape = 16, size = 1.5)
  }
  if( eq.axes ){
    plot.out <- plot.out +
      ylim( 0, round( max.lims[2] ) ) +
      xlim( 0, round( max.lims[2] ) )
  }
  # plot( plot.out )
  return( invisible( plot.out ) )
}

create_df_qq <- function(ps, ci = 0.95){
  n_groups <- 1
  if(!is.null(names(ps))){
    n_groups <- length(unique(names(ps)))
    message("Detected ", n_groups, " groups in the data.")
  }
  df_list <- purrr::map(seq_len(n_groups), function(group_no){
    cur_name <- unique(names(ps))[group_no]
    cur_ps <- ps[which(names(ps) == cur_name)]
    n  <- length(cur_ps)
    ps.sorted <- sort(cur_ps, index.return = TRUE)
    df_out <- tibble(
      observed = -log10(ps.sorted$x),
      expected = -log10(ppoints(n)),
      clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
      cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1)),
      group = cur_name
    )
  })
  df <- dplyr::bind_rows(df_list)
  return(df)
}

get_x_coords <- function( n.points, n.out ){
  if( n.out > n.points ){
    stop( "Too large number of points expected!" )
  }
  if( n.out <= 0 ){
    stop( "Wrong number n.out given!" )
  }
  x.coords <- -log10( ppoints( n.points ) )
  return( x.coords[ 1:n.out ] )
}
