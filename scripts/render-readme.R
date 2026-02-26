quarto <- Sys.which("quarto")
if (identical(quarto, "")) {
  stop(
    "Quarto is not installed or not on PATH. Install Quarto from https://quarto.org/ and retry.",
    call. = FALSE
  )
}

status <- system2(quarto, c("render", "README.qmd"))
if (!identical(status, 0L)) {
  stop("Quarto render failed.", call. = FALSE)
}

message("Rendered README.md and README.html from README.qmd")
