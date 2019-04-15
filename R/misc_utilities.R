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

# Check the length of components of a list
#
# @param values A list whose components should be checked
# @param cutoff A minimum value to check for
#
# @return a vector of logicals
#
LengthCheck <- function(values, cutoff = 0) {
  return(vapply(
    X = values,
    FUN = function(x) {
      return(length(x = x) > cutoff)
    },
    FUN.VALUE = logical(1)
  ))
}