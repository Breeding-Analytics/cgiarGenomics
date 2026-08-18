#' Find the bundled Beagle JAR file
#'
#' Searches for a Beagle JAR in the package's inst/extdata directory.
#' Returns the path to the first .jar file found that matches "beagle" in the name.
#'
#' @return Path to the Beagle JAR file, or NULL if not found
#' @export
#'
#' @examples
#' find_beagle_jar()
find_beagle_jar <- function() {
  extdata_dir <- system.file("extdata", package = "cgiarGenomics")
  if (extdata_dir == "") {
    return(NULL)
  }
  jars <- list.files(extdata_dir, pattern = "beagle.*\\.jar$", 
                     full.names = TRUE, ignore.case = TRUE)
  if (length(jars) == 0) {
    return(NULL)

  }
  # Return the first match (most recent if multiple exist)
  jars <- sort(jars, decreasing = TRUE)
  return(jars[1])
}


#' Find a working Java executable
#'
#' Attempts to locate Java on the system PATH by running `java -version`.
#' Returns the path to the java executable if found and functional.
#'
#' @return Path to the java executable (typically "java"), or NULL if not found
#' @export
#'
#' @examples
#' find_java()
find_java <- function() {
  # First try JAVA_HOME
  java_home <- Sys.getenv("JAVA_HOME", unset = "")
  if (nzchar(java_home)) {
    java_exe <- file.path(java_home, "bin", "java")
    if (.Platform$OS.type == "windows") {
      java_exe <- paste0(java_exe, ".exe")
    }
    if (file.exists(java_exe)) {
      # Verify it actually runs
      status <- tryCatch(
        system2(java_exe, "-version", stdout = FALSE, stderr = FALSE),
        error = function(e) 1L
      )
      if (status == 0) return(java_exe)
    }
  }
  
  # Try java on PATH
  status <- tryCatch(
    system2("java", "-version", stdout = FALSE, stderr = FALSE),
    error = function(e) 1L
  )
  if (status == 0) return("java")
  
  return(NULL)
}


#' Get the Java version string
#'
#' @param java_exe Path to java executable (default: "java")
#' @return Character string with the Java version, or NULL if unavailable
#' @export
get_java_version <- function(java_exe = "java") {
  tryCatch({
    out <- system2(java_exe, "-version", stdout = TRUE, stderr = TRUE)
    # Java version is typically on the first line
    version_line <- grep("version", out, value = TRUE, ignore.case = TRUE)[1]
    return(version_line)
  }, error = function(e) {
    return(NULL)
  })
}


#' Get default number of threads for Beagle
#'
#' Returns half the available cores, capped at 4.
#' Safe default that leaves the machine responsive.
#'
#' @return Integer number of threads
#' @export
get_default_beagle_threads <- function() {
  n_cores <- tryCatch(
    parallel::detectCores(),
    error = function(e) 2L
  )
  if (is.na(n_cores) || is.null(n_cores)) n_cores <- 2L
  min(4L, max(1L, n_cores %/% 2L))
}


#' Check Beagle requirements
#'
#' Verifies that Java, the Beagle JAR and the vcfR package are all available.
#' Returns a list with status information useful for UI feedback.
#'
#' @return A list with elements:
#'   - java_ok: logical, whether Java was found
#'   - java_exe: path to java executable (or NULL)
#'   - java_version: version string (or NULL)
#'   - beagle_ok: logical, whether Beagle JAR was found
#'   - beagle_path: path to Beagle JAR (or NULL)
#'   - vcfr_ok: logical, whether the vcfR package is installed
#'   - ready: logical, TRUE if all requirements are met
#'   - message: human-readable status message
#' @export
check_beagle_requirements <- function() {
  java_exe <- find_java()
  java_ok <- !is.null(java_exe)
  java_version <- if (java_ok) get_java_version(java_exe) else NULL

  beagle_path <- find_beagle_jar()
  beagle_ok <- !is.null(beagle_path)

  # Beagle round-trips through VCF, so vcfR is required to read results back.
  vcfr_ok <- requireNamespace("vcfR", quietly = TRUE)

  ready <- java_ok && beagle_ok && vcfr_ok

  problems <- character(0)
  if (!java_ok) {
    problems <- c(problems, "Java not found. Install Java 8+ from https://adoptium.net/temurin/releases/")
  }
  if (!beagle_ok) {
    problems <- c(problems, sprintf(
      "Beagle JAR not found. Place it in: %s",
      system.file("extdata", package = "cgiarGenomics")
    ))
  }
  if (!vcfr_ok) {
    problems <- c(problems, 'R package vcfR not installed. Run: install.packages("vcfR")')
  }

  msg <- if (ready) {
    sprintf(
      "Beagle ready. Java: %s | JAR: %s",
      if (is.null(java_version)) "found" else java_version,
      basename(beagle_path)
    )
  } else {
    paste(problems, collapse = " | ")
  }

  list(
    java_ok = java_ok,
    java_exe = java_exe,
    java_version = java_version,
    beagle_ok = beagle_ok,
    beagle_path = beagle_path,
    vcfr_ok = vcfr_ok,
    ready = ready,
    message = msg
  )
}
