files <- Sys.glob(paste0("*", SHLIB_EXT))
xfiles <- Sys.glob(file.path("..","inst","lib","*"))
dest <- file.path(R_PACKAGE_DIR, paste0('libs', R_ARCH))
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

file.copy(files, dest, overwrite = TRUE)
file.copy(xfiles, dest, overwrite = TRUE)
file.copy(file.path("..", "NLOPT_LICENSE"), R_PACKAGE_DIR)
if(file.exists("symbols.rds"))
    file.copy("symbols.rds", dest, overwrite = TRUE)
