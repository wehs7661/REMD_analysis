#!/usr/bin/env python
import os
import sys

def main():
    s = 'This Python script prepares all files needed to run a replica exchange, including' \
        ' .mdp files and .tpr files, given the .gro file(s) and a template .mdp file.' \
        ' The .gro file could be only one common .gro file for all the replicas, or one for' \
        ' each replica. (Check REMD_conf_prep.py for methods for preparing one .gro for each' \
        ' replica.) In addition, a job submission script can also be generated.' 

    print(s)

    
    # First request user input
    platform = input('Are you using this script locally or in HPC? (local/HPC): ')
    jobname = input("Please input the name of the job: ")
    prefix = input("Please input the prefix of the file names: ")
    N = input("Please input the number of replicas: ")
    nodes = input("Please input the number of requested nodes: ")
    simutime = input("Please input the simulation time (in hour): ")
    HPC = input("Will the job be submitted to Summit or Bridges? ")
    conf = input("Are the configurations for different states the same? (yes/no) ") 
    jobscript = input("Do you want to write a job submission script? (yes/no) ")
    sub = input("Do you want to submit the job? (yes/no) ")
    # REMD = input('Are you preparaing HREMD or TREMD? (HREMD/TREMD)')
    

    # Identify the current working platform
    if platform == 'HPC':
        gmx_cmd = 'gmx_mpi'
        if HPC == 'Summit' or HPC == 'summit':
            os.system('ml gromacs/2018.3')
        elif HPC == 'Bridges' or HPC == 'bridges':
            os.system('module load gromacs/2018_cpu')
            os.system('module load mpi/intel_mpi')
        else:
            print('Error: Wrong HPC specified. Available options are Summit and Bridges.')
            sys.exit()
    elif platform == 'local':
        gmx_cmd = 'gmx'
    else:
        print('Error: Wrong platform specified. Available options are local and HPC.')

    
    # Prepare all the input files
    print('Preparaing the files needed for the Hamiltonian replica exchange ...')
    for i in range(int(N)):
        if conf == 'yes':
            os.system('mkdir state_%s' % i)
            os.system('cp *gro *top *ndx state_%s' % i)
        elif conf == 'no':
            os.system('cp *top *ndx *state_%s' % i)
        else:
            print('Bad input for the question: Are the configurations for different states the same? (yes/no')
            sys.exit()

        os.system('cp *template.mdp state_%s/%s.mdp' % (i, prefix))
        os.system('sed -i -e "s/init-lambda-state        = X/init-lambda-state        = %s/g" state_%s/%s.mdp' %(i, i, prefix))
        os.chdir('state_%s' % i)
        os.system('%s grompp -f %s.mdp -c %s.gro -p %s.top -n %s.ndx -o %s.tpr -maxwarn 4' % (gmx_cmd, prefix, prefix, prefix, prefix, prefix))
        os.chdir('../')
    

    # Write the job submission script and submit it if needed
    if jobscript == 'yes':
        print('Writing job submission script ...')

        job_dir = ''
        for i in range(int(N)):
            job_dir += 'state_%s ' %i
        
        job_dict = {'jobname': jobname, 'N': N, 'N_cube': int(N) ** 3, 'nodes': nodes, 
                    'simutime': simutime, 'prefix': prefix, 'job_dir': job_dir, 
                    'task-per-nodes': int(int(N) / int(nodes)), 
                    'cpus-per-task-s': int(24 / (int(N) / int(nodes))), 
                    'cpus-per-tasks-b': int(28 / (int(N) / int(nodes)))}
        
        if HPC == 'Summit' or HPC == 'summit':
            text = """#!/bin/sh
#SBATCH --job-name {jobname}
#SBATCH --qos normal
#SBATCH --nodes={nodes}
#SBATCH --tasks-per-node={tasks-per-nodes}
#SBATCH --cpus-per-task={cpus-per-tasks-s}
#SBATCH --time {simutime}:00:00
#SBATCH --partition shas

module load gromacs/2018.3

# export OMP_NUM_THREADS=24

mpirun -np {N} gmx_mpi mdrun -deffnm {prefix} -dhdl {prefix}_dhdl.xvg -replex 100 -nex {N_cube} -multidir {job_dir}
            """.format(**job_dict)
            outfile = open('HREMD_run.sh', 'w')
            outfile.write(text)

            if sub == 'yes':
                os.system('sbatch HREMD_run.sh')
                print('Job submitted!')
            elif sub == 'no':
                pass
            else: 
                print('Error: Bad input for the question: Do you want to submit the job? (yes/no).')
        
        if HPC == 'Bridges' or HPC == 'bridges':
            text = """#!/bin/sh
#SBATCH -N {nodes}
#SBATCH -p RM
#SBATCH -t {simutime}:00:00
#SBATCH --ntasks-per-node={task-per-nodes}

module load gromacs/2018_gpu
module load mpi/intel_mpi

# export OMP_NUM_THREADS=24
mpirun -np {N} gmx_mpi mdrun -deffnm {prefix} -dhdl {prefix}_dhdl.xvg -replex 100 -nex {N_cube} -multidir {job_dir}
            """.format(**job_dict)
            outfile = open('HREMD_run.sh', 'w')
            outfile.write(text)

            if sub == 'yes':
                os.system('sbatch HREMD_run.sh')
                print('Job submitted!')
            elif sub == 'no':
                pass
            else: 
                print('Error: Bad input for the question: Do you want to submit the job? (yes/no).')


    # Check if all the files were successfully generated
    # mdout.mdp, PLCpep7.mdp, PLCpep.ndx, PLCpep7.top, PLCpep7.tpr
