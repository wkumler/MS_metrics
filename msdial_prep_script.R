
msdial_dir <- "msdial"

if(dir.exists(msdial_dir)){
  unlink(msdial_dir, recursive = TRUE, force = TRUE)
} else {
  dir.create(msdial_dir)
}
file.copy(from = "mzMLs", to = msdial_dir, recursive = TRUE)
file.copy("msdialparams.txt", paste0(msdial_dir, "/msdialparams.txt"))

msdial_cmd <- paste(
  r"("C:\Program Files\MSDIAL ver.4.9.221218 Windowsx64\MsdialConsoleApp.exe" lcmsdda)",
  paste0("-i .\\", msdial_dir, "\\"),
  paste0("-o .\\", msdial_dir, "\\"),
  paste0("-m .\\", msdial_dir, "\\msdialparams.txt"),
  "-p"
)

system(msdial_cmd)
