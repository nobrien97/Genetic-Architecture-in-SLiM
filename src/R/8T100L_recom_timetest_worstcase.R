
# Worst case scenario test run for 8T Selection

system("/home/$USER/SLiM/slim -s 123 -d Ne=10000 -d QTL_cov=0.5 -d nloci=500 -d locisigma=10.0 -d pleiorate=1.0 -d delmu=1.0 -d rwide=1.346e-5 -d pleio_cov=0.5 -d delchr=10 /home/$USER/SLiM/Scripts/null_recom_8T100L_nowrite.slim")
