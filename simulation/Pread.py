import argparse
import re
import pickle as pkl
import numpy as np
import pandas as pd
import os 

#NOTE: The TS and logInterval variables may need to be updated manually!!!!!!
TS = 0.0005 #Timestep in picoseconds.
logInterval = 100 #Interval between each written line in the log file. 
#-------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-f','--filename', type = str, nargs = '*', help = 'Input logfile names to read from.', required = True)
parser.add_argument('-n','--name', type = str, nargs = 1, help = 'Input name to assign the pickled dataframe output.', required = True)
parser.add_argument('-w','--window', type = int, nargs = 1, help = 'Size (in ps) of window over which to average pressure data. Sefta used 100 ps in her work, but 10 ps is probably more reasonable.', required = True)
parser.add_argument('-a','--append', type = str, nargs = 1, help = 'Optional. Append data from the filenames flag to the end of the pickle specified in this flag, and store the new dataframe with the name specified by the name flag.', required = False) 
args = parser.parse_args()
#------------------------------------------------
# Pulls bubble pressure data from log file(s) and processes it.
N = int(np.ceil(args.window[0]/TS/logInterval))
if args.append:
    df = pd.read_pickle(args.append[0])
gate_match = re.compile('.*runLoop.*').match
start_match = re.compile('Inserted.*He so far.*').match
list_match = re.compile('Step Temp Press TotEng v_countHe v_bubPress.*').match
stop_match = re.compile('Loop time of.*').match
burst_match = re.compile('^Apparent Relief Method: Bursting.*').match
punch_match = re.compile('^Apparent Relief Method: Loop Punching.*').match
bottom_match = re.compile('^Apparent Relief Method: Punching Through Bottom.*').match
error_match = re.compile('^ERROR.*').match
nv_match = re.compile('^nv:.*').match
time_match = re.compile('^Time Limit Reached.*').match
quit_match = re.compile('^quit.*').match
n = 1
for filename in args.filename:
    useful = True 
    # File Identification
    actualName = filename.split('/')[-1]
    splitName = actualName.partition('_')
    fileID = splitName[2]
    splitID = fileID.split('_')
    runNum = int(splitID[0])
    almostLig = [i for i in splitID[5] if i.isdigit() or i=='.']
    lig = ''.join(map(str,almostLig))
    almostTemp = [int(i) for i in splitID[1] if i.isdigit()]
    temp = int(''.join(map(str,almostTemp)))
    surfOrient = splitID[2]
    bubShape = splitID[4]
    almostRs = [int(i) for i in splitID[3] if i.isdigit()]
    Rs = int(''.join(map(str,almostRs)))
    # Read data
    PressVec = []
    StepVec = []
    HeVec = []
    lines = []
    gate_flag = False
    read_flag = False
    list_flag = False
    with open(filename, 'rb') as f:
        f.seek(-2, os.SEEK_END)
        while f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
        last = f.readline().decode()
    if not quit_match(last): 
        print('Premature end found in '+fileID) 
        useful = False
        continue
    for line in open(filename):
        if nv_match(line):
            nvLine = line
        if time_match(line):
            reliefType = 'timeOut'
            print('Time ran out in ' + fileID)
            break 
        if error_match(line):
            reliefType = 'ERROR'
            print('Error found in ' + fileID)
            break 
        if burst_match(line):
            reliefType = 'burst'
            break 
        if punch_match(line):
            reliefType = 'loopPunch'
            break 
        if bottom_match(line):
            reliefType = 'bottomPunch'
            break
        if gate_flag == True:
            if stop_match(line) or line == '\n':
                read_flag = False
                list_flag = False
                continue 
            elif start_match(line):
               read_flag = True 
            elif read_flag == True:
                if list_flag == True:
                    lines.append(str(line))	
                elif list_match(line):
                    list_flag = True		
        if gate_match(line):
            gate_flag = True
    almostnv = [int(i) for i in nvLine if i.isdigit()]
    nv = int(''.join(map(str,almostnv)))
    lines2 = []
    for item in lines:
        item = item.strip()
        item = item.split()
        lines2.append(item)
    for bundle in lines2:
        pitem = bundle[-1]
        sitem = bundle[0]
        hitem = bundle[-2]
        PressVec.append(float(pitem))
        StepVec.append(float(sitem))
        HeVec.append(float(hitem))
    # Convert steps to ps (rough way for now)
    TimeVec = [i*TS for i in StepVec]
    rawP = np.vstack((TimeVec,PressVec,HeVec))
    PressVec = []
    StepVec = []
    HeVec = []
    lines = []
    lines2 = []
    TimeVec = []
    #Average raw pressure data. 
    avgTime = np.convolve(rawP[0],np.ones((N,))/N,mode='valid')
    avgPress = np.convolve(rawP[1],np.ones((N,))/N,mode='valid')
    avgHe = np.convolve(rawP[2],np.ones((N,))/N,mode='valid') 
    avgP = np.vstack((avgTime,avgPress,avgHe))
    # Find maximum (average) pressure during the run. 
    maxPress = np.amax(avgP[1])
    maxPTime = avgP[0][np.argmax(avgP[1])]
    maxPHe = avgP[2][np.argmax(avgP[1])] #This line gives maxPHe as a float. Should be within one atom of the actual number of heliums at relief time, provided the window is 10 ps.  
    #maxPHe = max(rawP[2]) #This should work okay, but be careful.  
    dfPart = pd.DataFrame([(int(runNum),float(temp),surfOrient,float(Rs),bubShape,float(lig),float(maxPTime),float(maxPress),float(maxPHe),int(nv),reliefType)],index=[fileID],columns=('runNum','temp(K)','surfOrient','Rs(ao)','bubShape','lig(ao)','maxPTime(ps)','maxPress(bar)','maxPHe','nv','reliefType'))
    if args.append:
        df = df.append(dfPart)
    else:
        if n == 1:
            df = dfPart
        else:
            df = df.append(dfPart)
    print('Processed data for '+fileID)
    n = n + 1
if useful == True:
    df.to_pickle(args.name[0]+'.df.pkl')
