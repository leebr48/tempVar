#!/bin/bash
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -J PreadAllSubmit
#SBATCH -A hammond
#SBATCH -t 0:01:00
#SBATCH -n 1
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=bfl3m8@mail.missouri.edu
#SBATCH --mail-type=ALL
#SBATCH -o PreadAllSubmit.out
#SBATCH -e PreadAllSubmit.error

# Activate the virtual enironment. 
source /home/bfl3m8/virtenv/bin/activate

# Run the Pread.py program.
python3 /group/hammond/blee/Project2/Pread.py -f /group/hammond/blee/Project2/comp_runs/run_*/log_* -n comp -w 10

# Deactivate the virtual environment. 
deactivate 
