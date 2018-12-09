mkdir out
mkdir ./out/mean_std
qsub -t 1-22 step1_mean_std.sh
