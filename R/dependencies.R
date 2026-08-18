# Optional dependency management -------------------------------------------
#
# cgiarGenomics keeps its hard dependency list to CRAN-only packages so that a
# plain `install.packages()` always succeeds. Features that need heavier or
# non-CRAN packages load them on demand and, if they are missing, fail with a
# copy-pasteable install command instead of a cryptic namespace error.

# Registry of optional packages: what they unlock and how to install them.
.optional_deps <- list(
  vcfR = list(
    feature = "reading and writing VCF files",
    functions = "read_vcf(), impute_beagle()",
    source = "CRAN",
    install = 'install.packages("vcfR")'
  ),
  randomForest = list(
    feature = "Random Forest imputation",
    functions = 'impute_gl(method = "random_forest")',
    source = "CRAN",
    install = 'install.packages("randomForest")'
  ),
  dartR.base = list(
    feature = "reading DArTSeq SNP files",
    functions = "read_DArTSeq_SNP()",
    source = "CRAN, but requires two Bioconductor packages first",
    install = paste0(
      'install.packages("BiocManager")\n',
      'BiocManager::install(c("SNPRelate", "snpStats"))\n',
      'install.packages("dartR.base")'
    )
  )
)


#' Require an optional package, with an actionable error message
#'
#' Internal helper. Checks that an optional dependency is installed and aborts
#' with installation instructions if it is not.
#'
#' @param pkg Name of the package to check.
#' @return Invisibly TRUE if available; aborts otherwise.
#' @keywords internal
require_optional <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(invisible(TRUE))
  }

  info <- .optional_deps[[pkg]]

  if (is.null(info)) {
    cli::cli_abort(c(
      "Package {.pkg {pkg}} is required but not installed.",
      "i" = 'Install it with: {.code install.packages("{pkg}")}'
    ))
  }

  cli::cli_abort(c(
    "Package {.pkg {pkg}} is required for {info$feature}.",
    "*" = "Needed by: {info$functions}",
    "*" = "Available from: {info$source}",
    "i" = "Install with:",
    " " = "{info$install}",
    "i" = "Or install every optional dependency at once with {.code cgiarGenomics::install_optional_deps()}."
  ))
}


#' Report the status of cgiarGenomics optional dependencies
#'
#' Shows which optional packages are installed and which features are currently
#' unavailable. Useful as a first troubleshooting step.
#'
#' @return Invisibly, a data.frame with one row per optional dependency.
#' @export
#'
#' @examples
#' \dontrun{
#' check_dependencies()
#' }
check_dependencies <- function() {
  status <- data.frame(
    package = names(.optional_deps),
    feature = vapply(.optional_deps, `[[`, character(1), "feature"),
    installed = vapply(
      names(.optional_deps),
      function(p) requireNamespace(p, quietly = TRUE),
      logical(1)
    ),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  cli::cli_h3("cgiarGenomics optional dependencies")
  for (i in seq_len(nrow(status))) {
    if (status$installed[i]) {
      cli::cli_alert_success("{.pkg {status$package[i]}} - {status$feature[i]}")
    } else {
      cli::cli_alert_danger("{.pkg {status$package[i]}} - {status$feature[i]} {.emph (not installed)}")
    }
  }

  if (any(!status$installed)) {
    cli::cli_alert_info(
      "Run {.code cgiarGenomics::install_optional_deps()} to install what is missing."
    )
  }

  invisible(status)
}


#' Install cgiarGenomics optional dependencies
#'
#' Installs the optional packages that unlock VCF reading, Random Forest
#' imputation and DArTSeq reading. The Bioconductor packages required by
#' \pkg{dartR.base} (\pkg{SNPRelate} and \pkg{snpStats}) are handled through
#' \pkg{BiocManager}.
#'
#' @param which Character vector of optional packages to install. Defaults to
#'   all of them.
#' @param ask Logical. If TRUE (default in interactive sessions), asks for
#'   confirmation before installing.
#'
#' @return Invisibly, the character vector of packages that were installed.
#' @export
#'
#' @examples
#' \dontrun{
#' # Everything
#' install_optional_deps()
#'
#' # Only what you need
#' install_optional_deps("vcfR")
#' }
install_optional_deps <- function(which = names(.optional_deps),
                                  ask = interactive()) {
  which <- match.arg(which, choices = names(.optional_deps), several.ok = TRUE)

  missing <- which[!vapply(which, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) == 0) {
    cli::cli_alert_success("All requested optional dependencies are already installed.")
    return(invisible(character(0)))
  }

  cli::cli_alert_info("Missing optional packages: {.pkg {missing}}")

  if (ask) {
    answer <- readline("Install them now? [y/N] ")
    if (!tolower(trimws(answer)) %in% c("y", "yes")) {
      cli::cli_alert_info("Aborted. Nothing was installed.")
      return(invisible(character(0)))
    }
  }

  # dartR.base needs two Bioconductor packages before it can be installed.
  if ("dartR.base" %in% missing) {
    bioc_deps <- c("SNPRelate", "snpStats")
    bioc_missing <- bioc_deps[
      !vapply(bioc_deps, requireNamespace, logical(1), quietly = TRUE)
    ]

    if (length(bioc_missing) > 0) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        cli::cli_alert_info("Installing {.pkg BiocManager} (needed for Bioconductor packages)...")
        utils::install.packages("BiocManager")
      }
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        cli::cli_abort(c(
          "Could not install {.pkg BiocManager}.",
          "i" = "Install it manually with {.code install.packages(\"BiocManager\")} and try again."
        ))
      }
      cli::cli_alert_info("Installing Bioconductor dependencies: {.pkg {bioc_missing}}")
      BiocManager <- asNamespace("BiocManager")
      BiocManager$install(bioc_missing, update = FALSE, ask = FALSE)
    }
  }

  cli::cli_alert_info("Installing from CRAN: {.pkg {missing}}")
  utils::install.packages(missing)

  still_missing <- missing[
    !vapply(missing, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(still_missing) > 0) {
    cli::cli_alert_warning(
      "These packages could not be installed: {.pkg {still_missing}}"
    )
    cli::cli_alert_info("Manual instructions:")
    for (p in still_missing) {
      cli::cli_text("{.strong {p}}:")
      cli::cli_verbatim(.optional_deps[[p]]$install)
    }
  } else {
    cli::cli_alert_success("All requested optional dependencies installed.")
  }

  invisible(setdiff(missing, still_missing))
}
