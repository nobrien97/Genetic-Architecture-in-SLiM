
# Worst case scenario test run for 8T Selection


TMPDIR <- Sys.getenv('TMPDIR')
USER <- Sys.getenv('USER')
PBSIND <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))




if (PBSIND == 1) {
	system.time(system("/home/$USER/SLiM/slim -s 123 -d Ne=100 -d nloci=500 -d locisigma=10.0 -d pleiorate=1.0 -d delmu=1.0 -d rwide=1.241e-4 -d pleio_cov=0.5 -d delchr=10 -d wsd=10.0 /home/$USER/SLiM/Scripts/memorytest/stabsel_recom_8T100L_nowrite.slim"))
} else {
	system.time(system("/home/$USER/SLiM/slim -s 123 -d Ne=100 -d nloci=500 -d locisigma=10.0 -d pleiorate=1.0 -d delmu=1.0 -d rwide=1.241e-4 -d pleio_cov=-0.5 -d delchr=10 -d wsd=10.0 /home/$USER/SLiM/Scripts/memorytest/stabsel_recom_8T100L_nowrite_mem.slim"))
}