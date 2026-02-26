quarto <- Sys.which("quarto")
if (identical(quarto, "")) {
  stop(
    "Quarto is not installed or not on PATH. Install Quarto from https://quarto.org/ and retry.",
    call. = FALSE
  )
}

targets <- c("README.qmd", "technical-notes.qmd")
for (target in targets) {
  status <- system2(quarto, c("render", target))
  if (!identical(status, 0L)) {
    stop(sprintf("Quarto render failed for %s.", target), call. = FALSE)
  }
}

message("Rendered README.md, README.html, and technical-notes.html")
