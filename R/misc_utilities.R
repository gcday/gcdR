# Get the number of threads provided by the current plan
#
# @return The number of threads (workers) for the current future plan, or 1 if no workers detected
#
#' @importFrom future plan
#' @export
PlanThreads <- function() {
  nthreads <- eval(expr = formals(fun = plan())$workers)
  return(nthreads %||% 1)
}