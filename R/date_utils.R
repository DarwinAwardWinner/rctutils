# Convert datetimes that are actually dates to dates. Also works on data frames,
# where it runs on each. A "datetime that's actually a date" is defined as one
# where the time is midnight (i.e. hours, minutes, and seconds are all 0).

#' Convert all-midnight POSIXt to Date
#'
#' Dates in databases are often stored as date-time objects with all the time
#' values set to midnight (i.e. 00:00:00.0).
#'
#' @param x A POSIXt vector.
#' @param force If TRUE, then forcibly convert `x` to a Date even if its times
#'   are not all midnight. This also speeds up the conversion by skipping the
#'   check, so it mayu be useful if you already know that all times are
#'   midnight.
#' @param eps_seconds "Fudge factor" for detecting midnight. As long as every
#'   date-time is within this many seconds of midnight, it is considered
#'   equivalent to midnight.
#'
#' @returns Either a Date vector if the vector was converted, or the original
#'   vector if not.
#'
#' @export
#'
#' @examples
#' library(lubridate)
#'
#' actually_a_date <- ymd_hms("2000-01-01 00:00:00")
#' print(actually_a_date)
#' class(actually_a_date)
#' actually_a_date <- datetime_to_date(actually_a_date)
#' print(actually_a_date)
#' class(actually_a_date)
#'
#' actually_a_datetime <- ymd_hms("2000-01-01 01:27:13")
#' print(actually_a_datetime)
#' class(actually_a_datetime)
#' actually_a_datetime <- datetime_to_date(actually_a_datetime)
#' print(actually_a_datetime)
#' class(actually_a_datetime)
datetime_to_date <- function(x, force = FALSE, eps_seconds = 1) {
    req_ns("lubridate")
    if (lubridate::is.POSIXt(x)) {
        if (force) {
            # If forcing, we don't care if it's midnight
            return(as.Date(x))
        } else {
            # We're not going to try to round 23:59:59 to the next day's
            # midnight, because that would change the date. So technically this
            # is "seconds after midnight".
            secs_from_midnight <-
                lubridate::hour(x) * 3600 +
                lubridate::minute(x) * 60 +
                lubridate::second(x)
            if (max(secs_from_midnight, na.rm = TRUE) <= eps_seconds) {
                return(as.Date(x))
            } else {
                # If it's not midnight, keep as date time
                return(x)
            }
        }
    } else {
        return(x)
    }
}

# Apply the above function to every datetime column in a data frame

#' Convert all-midnight POSIXt columns to Date columns in a data frame
#'
#' @param df A data frame, potentially containing POSIXt columns.
#' @param ... Additional arguments to [datetime_to_date()]. All other columns
#'   are left alone.
#'
#' @returns The original data frame, with all POSIXt columns mutated by
#'   [datetime_to_date()].
#' @export
#'
#' @examples
#' library(lubridate)
#'
#' mydf <- data.frame(
#'     chr_col = "a",
#'     num_col = 1,
#'     fac_col = factor(letters)[1],
#'     actually_a_date = ymd_hms("2000-01-01 00:00:00"),
#'     actually_a_datetime = ymd_hms("2000-01-01 01:27:13")
#' )
#' sapply(mydf, \(x) class(x)[1])
#' sapply(datetimes_to_dates(mydf), \(x) class(x)[1])
datetimes_to_dates <- function(df, ...) {
    req_ns("lubridate")
    mutate(df, across(where(lubridate::is.POSIXt), \(x) datetime_to_date(x, ...)))
}
