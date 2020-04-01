#!/bin/sh
echo This shell script prepares all the files needed to run a Hamiltonian replica exchange, including .mdp files and .tpr files, given the common .gro file and a template .mdp file. A job submission script will also be generated. 

set -e   # exit upon error
read -p "Please input the name of the job: " jobname
read -p "Please input the prefix of the file names: " name
read -p "Please input the number of replicas: " N
read -p "Please input the number of requested nodes: " nodes
read -p "Please input the simulation time (in hour): " simutime
read -p "Will the job be submitted to Summit or Bridges? " HPC
read -p "Are the configurations for different states the same? (yes/no) " conf 
read -p "Do you want to write a job submission script? (yes/no) " jobscript
read -p "Do you want to submit the job? (yes/no) " sub

if [ ${HPC} == "Summit" ] || [ ${HPC} == "summit" ]
then
    ml gromacs/2018.3
elif [ ${HPC} == "Bridges" ] || [ ${HPC} == "bridges" ]
then 
    module load gromacs/2018_cpu
    module load mpi/intel_mpi
else
    echo "Wrong selected HPC. Available options are Summit and Bridges."
fi


echo Preparing the files needed for the Hamiltonian replica exchange ...

for (( i=0; i<$N; i=i+1 ))
do 
    if [ ${conf} == "yes" ]
    then
        mkdir state_${i}
        cp *gro *top *ndx state_${i}
    elif [ ${conf} == "no" ]
    then
        cp *top *ndx state_${i}
    else
        echo "Bad input for the last question: Are the configurations for different states the same? (yes/no). "
    fi
    cp *template.mdp state_${i}/${name}_restr4_${i}.mdp
    sed -i -e "s/init-lambda-state        = X/init-lambda-state        = ${i}/g" state_${i}/${name}_restr4_${i}.mdp
    cd state_${i}
    gmx_mpi grompp -f ${name}_restr4_${i}.mdp -c ${name}.gro -p ${name}.top -n ${name}.ndx -o ${name}.tpr -maxwarn 4
    cd ../
done

if [ ${jobscript} == "yes" ]
then
    echo Writing job submission script ...

    dir=""
    for (( i=0; i<$N; i=i+1 ))
    do 
        dir+="state_${i} "
    done

    if [ ${HPC} == "Summit" ] || [ ${HPC} == "summit" ]
    then
        touch HREX_run.sh
        d="#!/bin/sh
    #SBATCH --job-name ${jobname}
    #SBATCH --qos normal
    #SBATCH --nodes=$nodes
    #SBATCH --tasks-per-node=$(( ${N} / ${nodes} ))
    #SBATCH --cpus-per-task=$(( 24 / (( ${N} / ${nodes} )) ))
    #SBATCH --time ${simutime}:00:00
    #SBATCH --partition shas

    module load gromacs/2018.3

    # export OMP_NUM_THREADS=24

    mpirun -np ${N} gmx_mpi mdrun -deffnm ${name} -dhdl ${name}_dhdl.xvg -replex 100 -nex $(( ${N} * ${N} * ${N})) -multidir ${dir}" 

        echo "${d}" >> HREX_run.sh 
        
        if [ ${sub} == "yes" ]
        then
            sbatch HREX_run.sh
            echo "Job submitted!"
        elif [ ${sub} == "no" ]
        then
            true
        else
            echo "Bad input for the question: Do you want to submit the job? (yes/no). "
        fi

        echo Complete!

    elif [ ${HPC} == "Bridges" ] || [ ${HPC} == "bridges" ]
    then 
        touch HREX_run.sh
        d="#!/bin/sh
    #SBATCH -N $nodes
    #SBATCH -p RM
    #SBATCH -t $simutime:00:00
    #SBATCH --ntasks-per-node=$(( ${N} / ${nodes} ))

    module load gromacs/2018_gpu
    module load mpi/intel_mpi

    # export OMP_NUM_THREADS=24
    mpirun -np ${N} gmx_mpi mdrun -deffnm ${name} -dhdl ${name}_dhdl.xvg -replex 100 -nex $(( ${N} * ${N} * ${N})) -multidir ${dir}"

        echo "${d}" >> HREX_run.sh 

        if [ ${sub} == "yes" ]
        then
            sbatch HREX_run.sh
            echo "Job submitted"
        elif [ ${sub} == "no" ]
        then
            true
        else
            echo "Bad input for the question: Do you want to submit the job? (yes/no). "
        fi

        echo Complete!
    else
        echo "Wrong selected HPC. Available options are Summit and Bridges."
    fi
else
    echo Complete!
fi
