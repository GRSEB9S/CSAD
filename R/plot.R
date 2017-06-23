#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @export
plot.admm <- function(x, y, ...) {
  tolerance <- x$tolerance

  theme_set(theme_bw())
  g_primal <- ggplot(subset(tolerance, type == "primal" & value != 0), aes(iter, value)) +
    geom_line(aes(color = tolerance)) +
    scale_y_log10() +
    xlab("iteration") +
    ylab(expression(paste("||", Theta - Z, "||"[F]))) +
    theme(legend.title=element_blank())

  g_dual <- ggplot(subset(tolerance, type == "dual" & value != 0), aes(iter, value)) +
    geom_line(aes(color = tolerance)) +
    scale_y_log10() +
    xlab("iteration") +
    ylab(expression(paste("||",rho(Z - Z[old]), "||"[F]))) +
    theme(legend.title=element_blank())

  plot_grid(g_primal, g_dual, nrow = 2)
}
