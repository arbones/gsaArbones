
# Generic GSAfisher function


GSAfisher <- function(p_values) {
    UseMethod("GSAfisher")
}

###############################################################################################
# Default method for GSAfisher

#' Gene-Set Analysis Using Fisher's Method
#' @aliases GSAfisher-default
#' @description The default method for GSAfisher, handling a single vector of p-values.
#'
#' @param p_values A vector of p-values
#' This function performs Fisher's combined probability test on a vector or list of p-values.
#' @return Returns a list containing Fisher's combined statistic and the combined p-value. It also returns a vector of all combined p-values.
#' @examples
#' p_values_single <- runif(10, 0, 0.25)
#' result_single <- GSAfisher(p_values_single)
#' print(result_single)
#' @export
GSAfisher <- function(p_values) {
    # Validate input
    if (length(p_values) == 0) {
        stop("No p-values provided")
    }
    if (any(p_values < 0 | p_values > 1)) {
        stop("P-values must be between 0 and 1")
    }

    # Calculate Fisher's Combined Statistic
    fisher_statistic <- -2 * sum(log(p_values))

    # Number of tests
    k <- length(p_values)

    # Calculate the p-value for the Fisher's statistic
    p_value_combined <- stats::pchisq(fisher_statistic, df = 2 * k, lower.tail = FALSE)

    # Create result as a list
    result <- list("Fisher_statistic" = fisher_statistic, "P_value" = p_value_combined, "p_values"=p_values)

    # Assign class for the result
    class(result) <- "GSAfisherResult"

    return(result)
}

###############################################################################################
# Print method for GSAfisherResult objects

#' Print GSAfisher Results
#'
#' Prints the GSAfisher analysis results directly to the console.
#'
#' @param x An object of class `GSAfisherResult`.
#' @param ... other arguments
#' @return returns the input object.
#' @export
print.GSAfisherResult <- function(x,...) {
    cat("Fisher's Combined Test Results:\n")
    cat("Fisher's Combined Statistic:", x$Fisher_statistic, "\n")
    cat("Combined P-value:", format.pval(x$P_value, digits = 3, eps=0.001), "\n")
}

###############################################################################################
# Summary method for GSAfisherResult objects

#' Summary of GSAfisher Analysis
#'
#' Generates a concise summary of the GSAfisher analysis results.
#'
#' @param object An object of class `GSAfisherResult`.
#' @param ... other arguments
#' @return Invisible returns a GSAfisherSummary object with key statistics from the Fisher's combined test.
#' @export
summary.GSAfisherResult <- function(object,...) {
    summary_result <- list(
        n.p_values=length(object$p_values),
        n.significant=sum(object$p_values<0.05),
        Fisher_statistic = object$Fisher_statistic,
        P_value = object$P_value,
        Significant = object$P_value < 0.05  # Example: a simple significance test
    )

    # Assign a class to the summary result
    class(summary_result) <- "GSAfisherSummary"

    return(summary_result)
}
###############################################################################################
# Print method for GSAfisherSummary objects

#' Print GSAfisher Results
#'
#' Prints the GSAfisher analysis results directly to the console.
#'
#' @param x An object of class `GSAfisherSummary`.
#' @param ... other arguments
#' @return Summary.
#' @export
print.GSAfisherSummary <- function(x,...) {
    cat("Number of p-values:",x$n.p_values,"\n")
    cat("Number of significant (< 0.05) p-values:",x$n.significant,"\n")
    cat("*********************************************\n")
    cat("Summary of Fisher's Combined Test Results:\n")
    cat("Fisher's Combined Statistic:", x$Fisher_statistic, "\n")
    cat("Combined P-value:", format.pval(x$P_value, digits = 3, eps=0.001), "\n")
    cat("Significant (P < 0.05):", ifelse(x$Significant, "Yes", "No"), "\n")
}


##########################################################################################
# Method to compute Fisher combined test for multiple data vectors ######################
##########################################################################################

###### p-values must be in a list
#' Gene-Set Analysis Using Fisher's Method for a list of vectors
#'
#' @aliases GSAfisher-multiple
#' @description The default method for GSAfisher, handling a single vector of p-values.
#'
#' This function performs Fisher's combined probability test on a list of p-values.
#' It can handle a list of p-values.
#' @param p_values_list  A list of p-values
#' @return Returns a vector with the p-value of all Fisher's combined statistic tests.
#' @examples
#' p_values_matrix <- matrix(replicate(500,runif(100, 0, 1)),byrow = TRUE,nrow = 500)
#' p_values_list <- split(p_values_matrix,row(p_values_matrix))
#' result_multiple <- GSAfisher.multiple(p_values_list)
#' plot(result_multiple)
#' @export
GSAfisher.multiple <- function(p_values_list) {
    # Check if input is a list
    if (!is.list(p_values_list)) {
        stop("Input must be a list of numeric vectors.")
    }

    # Initialize a vector to store combined p-values
    combined_p_values <- numeric(length(p_values_list))

    # Loop over each set of p-values in the list
    for (i in seq_along(p_values_list)) {
        # Check if each element in the list is a numeric vector
        if (!is.numeric(p_values_list[[i]])) {
            stop("Each element in the list must be a numeric vector.")
        }

        # Apply GSAfisher to each vector of p-values
        combined_p_values[i] <- GSAfisher(p_values_list[[i]])$P_value
    }

    class(combined_p_values) <- "GSAfisherMultipleResult"
    return(combined_p_values)
}


## Implement a method to visualize the obtained results.

#' Plot of GSAfisher.multiple Analysis
#'
#' Generates a concise summary of the GSAfisher analysis results.
#'
#' @param x An object of class `GSAfisherMultipleResult`.
#' @param ...  Other plot arguments
#' @return An histogram of all combined p-values
#' @export
plot.GSAfisherMultipleResult <- function(x, ...) {
    # Check if input p_values is numeric
    if (!is.numeric(x)) {
        stop("Input must be numeric.")
    }

    graphics::hist(x,
         main="Histogram of all p-values for the GSAfisher")

}





###############################################################################################
