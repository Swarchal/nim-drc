# Package

version = "0.1.0"
author = "Scott Warchal"
description = "Dose-response curve fitting with no dependencies"
license = "MIT"
srcDir = "src"

# Dependencies

requires "nim >= 2.0.0"

task figures, "Render the README figures to docs/images (needs plotter and vl-convert)":
  # plotter is not a dependency of drc, so point --path at a checkout of it.
  let plotterPath = getEnv("PLOTTER_PATH", "../nim-plotter/src")
  exec "nim c -r --hints:off --path:src --path:" & plotterPath & " examples/figures.nim"
