
msdial_cmd <- paste(
  r"("C:\Program Files\MSDIAL ver.4.9.221218 Windowsx64\MsdialConsoleApp.exe" lcmsdda)",
  "-i .\\msdial\\",
  "-o .\\msdial\\",
  "-m .\\msdial\\msdialparams.txt",
  "-p"
)

system(msdial_cmd)
