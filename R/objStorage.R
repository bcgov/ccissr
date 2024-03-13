#' Download a file from Object Storage
#'
#' @param prefix passed to [`aws.s3::s3sync()`]. Must be a *path*
#' to a remote folder. Use trailing slash (`~/`).
#' @param bucket passed to [`aws.s3::s3sync()`].
#' @param path passed to [`aws.s3::s3sync()`]. Must be an *existing path*
#' to a local folder.
#' 
#' @details
#'   It is assumed that th object storage connection details and
#'   credentials have already been set up as system environment
#'   variables (notably "AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY" 
#'   and "AWS_S3_ENDPOINT").
#'   Note that all objects inside `path` will be synced to `path`
#'
#' @return NULL
#' 
#' @export
dwnldFromObjSto <- function(prefix, bucket, path) {
  if (requireNamespace("aws.s3")) {
    invisible({
      aws.s3::s3sync(path = dirname(path), 
                     bucket = bucket, prefix = prefix, region = "", 
                     direction = "download", create = TRUE)
    })
  } else {
    stop("Please install 'aws.s3' package to use this function")
  }
}