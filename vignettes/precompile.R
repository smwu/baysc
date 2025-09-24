# Pre-compiled vignettes that have a long runtime
# Must manually move image files from baysc/ to baysc/vignettes/ after knit
knitr::knit("vignettes/baysc.Rmd.orig", output = "vignettes/baysc.Rmd")

