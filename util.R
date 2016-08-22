# util.R
#
# Assorted utilities.

# Creates a directory at given path, optionally logging the action.
# If it already exists, the directory is not created.
dir.create.safe <- function(path, log=T) {
  if (!dir.exists(path)) {
    if (log) {
      loginfo("Creating directory: %s", path)
    }

    dir.create(path)
  }
}

# Confidence interval stat summary
ci.fun <- function(d) {
  test <- t.test(d)
  ci <- test$conf.int
  m <- as.double(test$estimate)

  data.frame(ymin=ci[1], ymax=ci[2], y=m)
}
