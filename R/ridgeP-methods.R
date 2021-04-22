#' @export
format.ridgeP <- function(x, ..., n = 6, m = 6) {
  N <- NROW(x)
  printed_rows <- min(n, N)
  printed_cols <- min(m, N)
  missing_rows <- N - printed_rows
  missing_cols <- N - printed_cols
  header <- sprintf("A %d x %d ridge precision matrix estimate with lambda = %f",
                    N, N, attr(x, "lambda"))

  body <- utils::capture.output(print(x[seq_len(printed_rows),
                                        seq_len(printed_cols)]))#, ...)

  # Defermine if more space when missing cols
  if (missing_cols > 0 & max(nchar(body)) < options()$width) {
    body <- paste(body, "\u2026")
  }

  # Create footer
  if (missing_rows > 0) {
    footer <- sprintf("\u2026 %d more rows", missing_rows)
    if (missing_cols > 0) {
      footer <- paste(footer, sprintf("and %d more columns", missing_cols))
    }
  } else {
    footer <- NULL
    if (missing_cols > 0) {
      footer <- sprintf("with %d more columns", missing_cols)
    }
  }

  fmt <- c(header, body, footer)
  return(fmt)
}

#' @export
print.ridgeP <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  return(invisible(x))
}

#' @export
summary.ridgeP <- function(object, ...) {
  cat(format(object)[1], "\n")
}

#' @export
as.matrix.ridgeP <- function(x, ...) {
  class(x) <- class(x)[-1]
  x
}

#' @export
plot.ridgeP <- function(x, ...) {
  edgeHeat(x, ...)
}
