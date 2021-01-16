#!/bin/bash
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -J PreadAllSubmit
#SBATCH -A hammond
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=bfl3m8@mail.missouri.edu
#SBATCH --mail-type=ALL
#SBATCH -o PreadAllSubmit.out
#SBATCH -e PreadAllSubmit.error

# Activate the virtual enironment. 
source /group/hammond/blee/Project2/virtenv/bin/activate

# Run the admin.py program.
python3 /group/hammond/blee/Project2/Pread.py -f /group/hammond/blee/Project2/runs/run_*/log_* -n all

# Deactivate the virtual environment. 
deactivate 
