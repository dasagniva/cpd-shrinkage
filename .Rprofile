# Preload local shared libraries globally so packages link correctly
libs <- c(
  "/home/shiva/.local/usr/lib/x86_64-linux-gnu/libgslcblas.so.0",
  "/home/shiva/.local/usr/lib/x86_64-linux-gnu/libgsl.so.27",
  "/home/shiva/.local/usr/lib/x86_64-linux-gnu/libsuitesparseconfig.so.7",
  "/home/shiva/.local/usr/lib/x86_64-linux-gnu/libamd.so.3",
  "/home/shiva/.local/usr/lib/x86_64-linux-gnu/libcolamd.so.3",
  "/home/shiva/.local/usr/lib/x86_64-linux-gnu/libglpk.so.40"
)
for (l in libs) {
  if (file.exists(l)) dyn.load(l, local = FALSE)
}

options(repos = c(
  stan = "https://stan-dev.r-universe.dev",
  CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"
))

source("renv/activate.R")
