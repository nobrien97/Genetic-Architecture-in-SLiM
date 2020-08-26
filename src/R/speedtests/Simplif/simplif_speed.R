

# Worst case scenario test run for 8T Null


TMPDIR <- Sys.getenv('TMPDIR')
USER <- Sys.getenv('USER')
PBSIND <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))




if (PBSIND == 1) {
  system.time(system("/home/$USER/SLiM/slim -s 123 -d locisigma=10.0 -d pleiorate=1.0 -d delmu=1.0 -d rwide=1.241e-4 -d pleio_cov=0.5 /home/$USER/SLiM/Scripts/speedtest/Simplif/null8T100L_nowrite.slim"))
} else {
  system.time(system("/home/$USER/SLiM/slim -s 123 -d locisigma=10.0 -d pleiorate=1.0 -d delmu=1.0 -d rwide=1.241e-4 -d pleio_cov=0.5 /home/$USER/SLiM/Scripts/speedtest/Simplif/stabsel_recom_8T100L_nowrite.slim"))
}