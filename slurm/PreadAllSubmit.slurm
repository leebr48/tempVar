#!/bin/bash
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -J PreadAllSubmit
#SBATCH -A hammond
#SBATCH -t 0:05:00
#SBATCH -n 1
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=bfl3m8@mail.missouri.edu
#SBATCH --mail-type=ALL
#SBATCH -o PreadAllSubmit.out
#SBATCH -e PreadAllSubmit.error

# Activate the virtual enironment. 
source /home/bfl3m8/virtenv/bin/activate

# Run the Pread.py program.
python3 /group/hammond/blee/Project2/Pread.py -f /group/hammond/blee/Project2/comp_runs/run_*comp0.98/log_* -n comp0.98 -w 10
python3 /group/hammond/blee/Project2/Pread.py -f /group/hammond/blee/Project2/comp_runs/run_*comp0.99/log_* -n comp0.99 -w 10

# Deactivate the virtual environment. 
deactivate 
