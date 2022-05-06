# load in packages
import sys
import os

# Handle command line arguments
## load cell

## create output directories
out_path = '/home/ckelley/SDinSlice/batch_files/'
try:
    os.mkdir(out_path)
except:
    pass

# main code for generating sbatch files
file_name = sys.argv[-1] + '.sbatch'
file = open(out_path+file_name,'w')

file.write('#!/bin/bash\n')
job_name = '#SBATCH --job-name=' + sys.argv[-1] + '\n'
file.write(job_name)
file.write('#SBATCH -t 24:00:00\n')
file.write('#SBATCH -p shared\n')
file.write('#SBATCH --nodes=1\n')
file.write('#SBATCH --ntasks-per-node=48\n')
file.write('#SBATCH --account csd403')
log_line = '#SBATCH -o /home/ckelley/SDinSlice/logs/' + sys.argv[-1] + '.log\n'
file.write(log_line)
err_line = '#SBATCH -e /home/ckelley/SDinSlice/errs/' + sys.argv[-1] + '.err\n'
file.write(err_line)
file.write('#SBATCH --mail-user=craig.kelley@downstate.edu\n')
file.write('#SBATCH --mail-type=end\n')
file.write('source /home/ckelley/.bashrc\n')
file.write('cd /home/ckelley/SDinSlice/\n')
file.write('module purge\n')
file.write('module load slurm\n')
file.write('module load cpu\n')
file.write('module load gcc/10.2.0\n')
file.write('module load openmpi/4.0.4\n')
run_line = 'mpiexec -n 48 nrniv -python -mpi SpatialModel.py cfgs/' + sys.argv[-1] + '.json\n'
file.write(run_line)
file.close()
