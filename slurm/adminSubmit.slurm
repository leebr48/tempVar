#!/bin/bash
#SBATCH -p Lewis
#SBATCH -N 1
#SBATCH -J adminSubmit
#SBATCH -A hammond
#SBATCH -t 0:05:00
#SBATCH -n 1
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=bfl3m8@mail.missouri.edu
#SBATCH --mail-type=FAIL
#SBATCH -o adminSubmit.out
#SBATCH -e adminSubmit.error

# Activate the virtual enironment. 
source ~/virtenv/bin/activate

# Run the admin.py program.
python3 /group/hammond/blee/Project2/admin.py -s

# Deactivate the virtual environment. 
deactivate 
