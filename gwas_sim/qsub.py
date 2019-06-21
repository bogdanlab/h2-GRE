import os
import subprocess
from gwas_sim import SumherModel
from os.path import join
import yaml
import inspect

TEMPLATE = """
#####################################
{premable}
#$ -S /bin/bash
#$ -cwd
#$ -N {name}
#$ -e {errfile}
#$ -o {logfile}
#$ -l h_data={memory},h_rt={duration}
#$ -t 1-{jobnum}
#####################################
echo "------------------------------------------------------------------------"
echo "Job started on" `date`
echo "------------------------------------------------------------------------"
{script}
echo "------------------------------------------------------------------------"
echo "Job ended on" `date`
echo "------------------------------------------------------------------------"
"""


def submit_script(script, premable, jobnum, memory, duration,
                name="job", logfile="log/output.$JOB_ID", 
                errfile="log/error.$JOB_ID", cleanup=True):
    """
    content of script
    """
    script_path = 'script.sh'
    
    script_content = TEMPLATE.format(script=script, premable=premable, 
                                     jobnum=jobnum, name=name, 
                                     logfile=logfile, errfile=errfile,memory=memory, duration=duration)
    script_content = inspect.cleandoc(script_content)

    open(script_path, 'w').write(script_content)
    
    try:
        subprocess.call('qsub ' + script_path, shell=True)
    finally:
        if cleanup:
            pass
            # os.remove(script_path)

def check_step(root_dir, config, step):
    assert(step in ['mean_std', 'beta', 'phe_g', 'phe', 'zsc'])
    print(join(root_dir, step))
    if not os.path.exists(join(root_dir, step)):
        return False
    else:
        files = os.listdir(join(root_dir, config['{}_dir'.format(step)]))
        if step != 'phe':
            return len(files) == config['num_part']
        else:
            return len(files) == 1


def submit_step(root_dir, config, step):
    """
    submit the job to the cluster, this script is used for UCLA Hoffman2 cluster
    the parameters are chosen to compute the UK Biobank typed data
    """

    premable = """
    export PATH=~/project-pasaniuc/software/anaconda3/bin:$PATH
    export PYTHONNOUSERSITE=True
    """

    assert(step in ['mean_std', 'beta', 'phe_g', 'phe', 'zsc'])
    script ='python gwas_sim.py --step {step} --root_dir {root_dir}'.format(step=step, root_dir=root_dir)
     
    if step in ['mean_std', 'phe_g', 'zsc']:
        script += ' --part_i $((SGE_TASK_ID - 1))'
        job_num = config['num_part']
    else:
        job_num = 1
    print(script)
    
    if step == 'mean_std':
        memory = '4G'
        duration='00:20:00'
    elif step == 'beta':
        memory = '4G'
        duration='00:05:00'
    elif step == 'phe_g' or step == 'zsc':
        memory = '6G'
        duration='01:00:00'
    elif step == 'phe':
        memory = '6G'
        duration='00:10:00'
    else:
        assert(False)
    
    submit_script(script=script, 
                    premable=premable,
                    jobnum=job_num,
                    name=step,
                    memory=memory,
                    duration=duration)


def submit(root_dir):
    """
    submit jobs for one model
    will check the status of each step and submit the corresponding step
    """
    # initialize the model
    model = SumherModel(root_dir)
    print(model.config)
    
    # mean_std
    for step in ['mean_std', 'beta', 'phe_g', 'phe', 'zsc']:
        if check_step(root_dir, model.config, step) == False:
            print('submitting step {}'.format(step))
            submit_step(root_dir, model.config, step)
            return
    print('all steps are finished!')
    
if __name__ == '__main__':
    for i in range(16): 
        submit('./out/model_{}'.format(i))
    
