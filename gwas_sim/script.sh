#####################################

    export PATH=~/project-pasaniuc/software/anaconda3/bin:$PATH
    export PYTHONNOUSERSITE=True
    
#$ -S /bin/bash
#$ -cwd
#$ -N zsc
#$ -e log/error.$JOB_ID
#$ -o log/output.$JOB_ID
#$ -l h_data=6G,h_rt=01:00:00
#$ -t 1-3
#####################################
echo "------------------------------------------------------------------------"
echo "Job started on" `date`
echo "------------------------------------------------------------------------"
python gwas_sim.py --step zsc --root_dir ./out/model_15 --part_i $((SGE_TASK_ID - 1))
echo "------------------------------------------------------------------------"
echo "Job ended on" `date`
echo "------------------------------------------------------------------------"