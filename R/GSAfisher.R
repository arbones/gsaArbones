
# Generic GSAfisher function
GSAfisher <- function(p_values) {
    UseMethod("GSAfisher")
}


# Default method for GSAfisher
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
    p_value_combined <- pchisq(fisher_statistic, df = 2 * k, lower.tail = FALSE)

    # Create result as a list
    result <- list("Fisher_statistic" = fisher_statistic, "P_value" = p_value_combined, "p_values"=p_values)

    # Assign class for the result
    class(result) <- "GSAfisherResult"

    return(result)
}


# Print method for GSAfisherResult objects
print.GSAfisherResult <- function(object) {
    cat("Fisher's Combined Test Results:\n")
    cat("Fisher's Combined Statistic:", object$Fisher_statistic, "\n")
    cat("Combined P-value:", format.pval(object$P_value, digits = 3, eps=0.001), "\n")
}

# Summary method for GSAfisherResult objects
summary.GSAfisherResult <- function(object) {
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

# Print method for GSAfisherSummary objects
print.GSAfisherSummary <- function(object) {
    cat("Number of p-values:",object$n.p_values,"\n")
    cat("Number of significant (< 0.05) p-values:",object$n.significant,"\n")
    cat("*********************************************\n")
    cat("Summary of Fisher's Combined Test Results:\n")
    cat("Fisher's Combined Statistic:", object$Fisher_statistic, "\n")
    cat("Combined P-value:", format.pval(object$P_value, digits = 3, eps=0.001), "\n")
    cat("Significant (P < 0.05):", ifelse(object$Significant, "Yes", "No"), "\n")
}
