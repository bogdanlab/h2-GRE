for chr_i in $(seq 1 22)
do
    if [ $chr_i -le 6 ]
    then
        qsub -l h_data=60G,h_rt=4:00:00,highp step3_merge.sh $chr_i
    elif [ $chr_i -ge 7 -a $chr_i -le 12 ]
    then
        qsub -l h_data=40G,h_rt=3:00:00,highp step3_merge.sh $chr_i
    elif [ $chr_i -ge 13 -a $chr_i -le 20 ]
    then
        qsub -l h_data=24G,h_rt=2:00:00,highp step3_merge.sh $chr_i
    else
        qsub -l h_data=6G,h_rt=1:00:00,highp step3_merge.sh $chr_i
    fi
done

