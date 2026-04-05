# Internal utility functions for agristat
#
# This file contains shared S3 methods and helper utilities used across
# all six modules. Full implementations will be added in Phase 2 onwards.
#
# Planned internal helpers:
#   - .validate_data_frame() : Validate that an object is a non-empty data frame
#   - .check_column_exists() : Verify required columns are present
#   - .format_pvalue()       : Format p-values for publication tables
#   - .stars_pvalue()        : Convert p-values to significance stars
#   - .is_whole_number()     : Check if a numeric value is a whole number

# S3 generic stubs (to be registered via NAMESPACE in later phases)
# print.agristat_result    <- function(x, ...) NextMethod()
# summary.agristat_result  <- function(object, ...) NextMethod()
# plot.agristat_result     <- function(x, ...) NextMethod()
