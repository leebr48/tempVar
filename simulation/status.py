import argparse 
import re
import os 
import shutil 

runPath = '/group/hammond/blee/Project2/runs/'
discrim = '00000'

#-------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-d','--dirname',type = str, nargs = '*', help = 'Input directory names to read from.',required = True)
parser.add_argument('-e','--error', action = 'store_true', default = False, help = 'Toggle on to look for files with errors.', required = False)
parser.add_argument('-g','--good', action = 'store_true', default = False, help = 'Toggle on to look for files that seem to have run correctly.', required = False)
parser.add_argument('--clear', action = 'store_true', default = False, help = 'Optional. Toggle on to remove the pesky restart files generated when LAMMPS performs a minimization.', required = False)
parser.add_argument('--kill', action = 'store_true', default = False, help = 'Optional. Toggle on to delete ALL restart files for runs which the program detects as being complete.', required = False)
parser.add_argument('--delete', action = 'store_true', default = False, help = 'Optional. Toggle on to remove DIRECTORIES for runs that do not have a log_ file.', required = False)
args = parser.parse_args()
#------------------------------------------------
# Set read flag.
error_match = re.compile('^ERROR:.*').match
end_match = re.compile('^quit.*').match
time_match = re.compile('^Time Limit Reached.*').match
# Scan log files. 
for dirname in args.dirname:
    if dirname[-1] != '/':
        actualName = dirname.split('/')[-1]
    if dirname[-1] == '/':
        actualName = dirname.split('/')[-2]
    splitName = actualName.partition('_')
    fileID = splitName[2] 
    goodFlag = True
    if dirname[-1] != '/':
        logname = dirname+'/log_'+fileID
    if dirname[-1] == '/':
        logname = dirname+'log_'+fileID
    fullPath = runPath+'run_'+fileID
    if args.clear:
        clearcmd = "find "+fullPath+" -name 'restart*' -type f ! -name 'restart*{:}' -delete".format(discrim)
        os.system(clearcmd)
    if args.error or args.good or args.kill or args.delete:
        try:
            for line in open(logname):
                if error_match(line):
                    goodFlag = False
                    if args.error:
                        print('Error found in '+fileID)
                if time_match(line):
                    goodFlag = False
                    if args.error:
                        print('Time ran out in '+fileID)
                last = line
            if not end_match(last):
                goodFlag = False
                if args.error:
                    print('Premature end found in '+fileID)
            if args.good and goodFlag == True:
                print(fileID+' looks good')
            if args.kill and goodFlag == True:
               killcmd = "find "+fullPath+" -name 'restart*' -type f -delete"
               os.system(killcmd)
        except FileNotFoundError:
            if args.error:
                print('No logfile found in '+fileID)
            if args.delete:
                print('Deleting '+dirname)
                shutil.rmtree(dirname)
