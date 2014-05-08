files <- Sys.glob(paste0("*", SHLIB_EXT))
xfiles <- Sys.glob("../inst/lib/*")
dest <- file.path(R_PACKAGE_DIR, paste0('libs', R_ARCH))

print(sprintf("selected dest = %s", dest))

dir.create(dest, recursive = TRUE, showWarnings = TRUE)

print(list.files(dest))

file.copy(files, dest, overwrite = TRUE)
file.copy(xfiles, dest, overwrite = TRUE)
if(file.exists("symbols.rds"))
    file.copy("symbols.rds", dest, overwrite = TRUE)

print(sprintf("contents of %s:", dest)
print(list.files(dest))
