# Messages are classed conditions (deconflate_message), so that callers can
# muffle them by class.

test_that("the notes of cm_read_inputs() are a classed message with the problems", {
  # An unrecognised column is a note.
  read <- function() cm_read_inputs(
    diseases = data.frame(disease = c("d1", "d2"), value = c(0.10, 0.15), colour = "red"),
    associations = data.frame(disease = "d1", associated_with = "d2", value = 2),
    impacts = data.frame(disease = c("d1", "d2"), value = c(2.5, 5)))
  expect_message(m <- read(), class = "deconflate_input_notes")
  msg <- tryCatch(read(), message = function(e) e)
  expect_s3_class(msg, c("deconflate_input_notes", "deconflate_message", "message", "condition"), exact = TRUE)
  expect_equal(msg$problems, attr(m, "problems"))
  expect_match(conditionMessage(msg), "\n$")
  expect_silent(suppressMessages(read()))
  # Muffled by class, nothing else is shown.
  expect_silent(withCallingHandlers(read(), deconflate_input_notes = function(e) invokeRestart("muffleMessage")))
})

test_that("cm_template() reports the files it wrote as a classed message", {
  dir <- file.path(tempdir(), "deconflate-test-messages")
  on.exit(unlink(dir, recursive = TRUE))
  expect_message(cm_template(dir), "Wrote 4 files", class = "deconflate_message")
})

test_that("cm_inform() ends the message with a newline", {
  expect_output(withCallingHandlers(cm_inform("a"), message = function(e) {
    cat(conditionMessage(e))
    invokeRestart("muffleMessage")
  }), "^a$")
  msg <- tryCatch(cm_inform("a", class = "x", n = 1), message = function(e) e)
  expect_equal(conditionMessage(msg), "a\n")
  expect_equal(msg$n, 1)
  expect_s3_class(msg, c("x", "deconflate_message"))
})
