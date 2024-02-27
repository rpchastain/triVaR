#' Magnitude, Duration, and Maximum
#'
#' @param x A numeric vector. The vector should be standardized so that values
#' above the threshold are positive and values below are negative. NA values
#' should be stripped from the vector.
#' @param below Boolean indicating if peak above or below threshold should be
#' calculated. If `TRUE` will calculate MDM values
#' for events above the threshold. Otherwise will calculate values below
#' threshold. Default is `TRUE`.
#' @param duration Minimum duration required for the MDM of and event to be
#' calculated. Default is `1`.
#'
#' @return Returns a `data.frame` with trivariate event summaries for all events
#' found in the vector. See `summarize_event` for further details.
#' @export
#'
#' @examples
#' observations <- c(5, 4, 3, 4, 2, 1, 1, 2, 4)
#' threshold <- 2
#'
#' centered_obs <- observations - threshold
#' mdm(centered_obs)
#' mdm(centered_obs, duration = 3)
mdm <- function(x, below = FALSE, duration = 1) {

  # Check to makes sure the user provided a vector that is at least as long as
  # the minimum duration

  ## TODO: Might want to make this a return
  if (length(x) < duration) {
    error_msg <- sprintf(
      'Insufficient data. Please provide a vector with length greater
          than %s.', duration
    )
    stop(strwrap(error_msg))
  }

  # Check to make sure that the threshold is crossed at least one time
  if (all(sign(x) < 0) | all(sign(x) > 0) | all(sign(x) == 0)) {
    return('No variability in data.')
  }

  if (below) {
    x <- -x
  }

  # Note where each change in sign occurs, indicating a transition across
  # the threshold
  x[x == 0] <- -1
  change_ind <- sign(x * c(x[2:length(x)], x[length(x)]))

  # Locate each possible start of an event
  start_index <- setdiff(which(change_ind == -1 & x <= 0), length(x)) + 1
  end_index <- setdiff(which(change_ind == -1 & x > 0), 1)

  # If there is only one event, and the event starts at the first observation or
  # ends at the last observation then both the start and end indices will be
  # integer(0). This will cause R to issue a warning. This check accounts for
  # this scenario and instead appropriately assigns the start and/or end index.
  if (length(start_index) == 0) {
    start_index <- 1
  }

  if (length(end_index) == 0) {
    end_index <- length(x)
  }

  # If the first event is a transition below the threshold then it can be
  # assumed the first event in the sequence was above the threshold and should
  # be included as a start event
  if (min(end_index) < min(start_index)) {
    start_index <- c(1, start_index)
  }

  # Similarly if the last event is a transition above the threshold then it
  # can be assumed the last event in the sequence was below the threshold and
  # should be included as an end event
  if (max(start_index) > max(end_index)) {
    end_index <- c(end_index, length(x))
  }

  bookend_index <- rbind.data.frame(start = start_index, end = end_index)

  # Filter out short duration events so that only events that meet minimum
  # duration requirements are counted
  if (!is.null(duration)) {
    long_events <- (bookend_index[2, ] - bookend_index[1, ]) >= (duration - 1)
    if (!any(long_events)) {
      warning('No events longer than minimum duration.')
    }
    bookend_index <- bookend_index[, long_events, drop = FALSE]
  }

  event_summary <- do.call(
    rbind.data.frame,
    lapply(bookend_index, summarize_event, x = x)
  )
  rownames(event_summary) <- c()

  event_summary
}

#' Summarize Peak-Over-Threshold event
#'
#' @param x Numeric vector of observations above/below a designated threshold.
#' @param ind Index indicating the start and end of single event.
#'
#' @return Returns a data frame summarizing the Trivaraite eveng including the
#' `sum`, `maximum`, `length`, and `start` and `end` indices for the events.
#' @export
#'
#' @examples
#' events <- c(2, 4, 5, -2, -2)
#' event_index <- c(1, 3)
#' summarize_event(events, event_index)
summarize_event <- function(x, ind) {
  event <- x[ind[1]:ind[2]]

  data.frame(
    start = ind[1],
    end = ind[2],
    magnitude = sum(event),
    duration = length(event),
    max = max(event)
  )
}
