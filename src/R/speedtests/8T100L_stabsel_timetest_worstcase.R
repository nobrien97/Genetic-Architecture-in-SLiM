
# Worst case scenario test run for 8T Selection

system.time(system("/home/$USER/SLiM/slim -s 123 -d rregion=0.5 -d Ne=10000 -d QTL_cov=0.5 -d pleiorate=1.0 -d delmu=0.0 -d rwide=1.346e-5 /home/$USER/SLiM/Scripts/stabsel_recom_8T100L_nowrite.slim"))
