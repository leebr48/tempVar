import re
import random
import math
import os
import sys 
import subprocess
import argparse 

# Variables to alter.
numRun = list(range(1,2,1)) #Number of runs to execute for each variable set. Enumerate these!!! 
t = [1500] #temperatures in Kelvin
orients = ['001'] #Surface orientations as a LIST. Choices: '001', '011', '111', and '310'.
Rs = [3] #Spherical bubble radii. MULTIPLY by ao to get radii in Angstroms. Focusing on smaller bubbles (below continuum limit) is probably best. 
l = [12] #Ligament thicknesses. MULTIPLY by ao to get thickness in Angstroms. Go down to 20 ao in the end!  
singleLayer = False #Type 'True' to append an extra run to the 'l' list above, creating runs whose ligament is one atom thick. (Note that this thickness is surface-orientation dependent.) Type 'False' to use only the vaules for 'l' given above. 
bubShape = ['sph'] #Specifies spherical bubble shape. ADD other bubble shapes to the list and CHANGE the switch case below. Hamond asked you to ignore anything but spherical bubbles for now.  
machine = 'Lewis'

#----------------------------------------------------------------------------------------

#Infile information.
if machine == 'Lewis':
    pots = '/group/hammond/blee/Project2/inputs' #MODIFY for the computer you are running on. 
    LatParamFile = '/group/hammond/blee/Project2/inputs/LatticeParams.txt' #MODIFY for the computer you are running on. 

# Mechanical components of the simulation. 
runLength = 10000 #Longest allowable run length in picoseconds. 
k = 8.617333262145e-5 #Boltzmann's constant in eV/K
m1 = 183.841 #mass of tungsten in amu
m2 = 4.0026032497 #mass of helium in amu
ts = 0.0005 #timestep for units metal, should be short enough for accuracy
stepsToEq = 1000000 #Timesteps to let the system run with the thermostat on without inserting any helium atoms. Note that this number is adjusted by multipliers in the code.    
betweenInsertionsPS = 5 #picoseconds between insertions
betweenInsertionsTS = int(betweenInsertionsPS/ts) #timesteps between insertions
betweenDumpsTS = int(10*betweenInsertionsTS) #timesteps between dumps
betweenDumpsTSFast = int(math.floor(betweenDumpsTS/100))
betweenRestartsTS = int(10*betweenInsertionsTS) #timesteps between restarts
betweenThermoTS = int(betweenInsertionsTS/100) #Timesteps between writing thermo data
fullRunTS = int(runLength/betweenInsertionsPS*betweenInsertionsTS) #Total length of run in timesteps
numberOfInsertions = int(fullRunTS/betweenInsertionsTS) #Total number of insertions 'scheduled' for the run
finalRunSeqTS = 2*betweenInsertionsTS #Number of timesteps to  run for determining if bursting or punching occured. 
statDamp = 100 #Damping parameter for nvt thermostat in timesteps.
volRat = 0.014 #Ratio of bubble volume to box volume. 0.014 was used by Faiza. 
WChunkMinIter = 100000 #Number of iterations to allow for the solid tungsten chunk to minimize before the spherical region is deleted. 
WChunkMinEvals = WChunkMinIter*10
HeInitPlaceIter = 1000 #Number of iterations to allow for an inserted He atom to move during the initial bubble creation. 
HeInitPlaceEvals = HeInitPlaceIter*10
bubMaxIter = 1000000 #Number of iterations to allow for the full bubble to minimize before the actual runs begin.
bubMaxEvals = bubMaxIter*10
HeMinIter = 100 #Number of minimizer iterations to allow for the helium atoms during the main run. 
HeMinEvals = HeMinIter*10 
botLigRat = 2 #Minimum ratio of the depth the of the tungsten block beneath the bubble to the ligament above the bubble. Hammond spitballed this as 'at least 2'.      
sideDepthRat = 2 #Minimum ratio of the side width of the tungsten block to the depth of the bubble. 2 is the 'magic number'. 
punchMultip = 0.3 #ao. How far to allow tungsten atoms to move from their normal positions before the program catches it as a 'loop punch'. 
punchMultipAdj = 0.5 #Multiplier to adjust the punchMultip variable for the read_restart fileType. 
#countWThresh = 40 #Number of tungsten atoms that must be within the loop punch 'catch' region for the behavior to be characterized as loop punching. 
HeLatParMultip = 1.45 #Factor used to make the HeLatpar larger than LatPar. 
xymod = 1 #Fraction of lattice parameter that the bubble can be randomly moved in the x and y direction.
verymin = 10 #Minimum box side length in lattice spacings. 
#maxHe = 5 #Factor that multiplies the number of originally inserted helium atoms to determine the most helium atoms that can be inserted before the run terminates. 
xrepeat = {'001':'1','011':'1','111':'3','310':'5'} #x axis repeat modulus for various surface orientations. 
yrepeat = {'001':'1', '011':'1', '111':'1', '310':'1'} #y axis repeat modulus for various surface orientations. 
zrepeat = {'001':'1', '011':'1', '111':'1', '310':'1'} #z axis repeat modulus for various surface orientations. 
xSpacing = {'001':'1', '011':'sqrt(2)', '111':'2/3*sqrt(6)', '310':'sqrt(8/5)'} #Multipliers for ao to get the latttice units right. 
ySpacing = {'001':'1', '011':'1', '111':'sqrt(2)', '310':'1'}
zSpacing = {'001':'1', '011':'sqrt(2)', '111':'sqrt(3)', '310':'sqrt(8/5)'}
fileTypes = ['in','read_restart1']
discrim = ''.join([i for i in str(betweenRestartsTS) if i == '0'])
runTime = {'short':'12:00:00','med':'36:00:00','long':'48:00:00'} #Lewis run lengths for small, medium, and large bubbles.
minThickness = {'001':0.4,'011':0.5,'111':0.75,'310':0.75} #Ligament thickness to specify for ligament to be one atom thick. 
maxJobs = 2000 #Maximum number of jobs for any single user on Lewis. 

#----------------------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument('-s','--submit', action = 'store_true', help = 'Optional. If used, this flag will cause the program to submit one job for each infile in the set to LEWIS and submit a dependent read_restart program as well.', required = False)
parser.add_argument('-c','--compress', type = float, help = 'Optional. If used, applies the given multiplicative factor to scale both horizontal dimensions of the box (>1 is expansion, <1 is compression).', required = False, default = 0) #If using, values of 0.99 or 0.98 will probably suffice for our purposes.
args = parser.parse_args()

#----------------------------------------------------------------------------------------

if args.compress != 0:
    compFac = args.compress
else:
    compFac = False

orientX001 = [1,0,0]; orientY001 = [0,1,0]; orientZ001 = [0,0,1]
orientX011 = [0,-1,1]; orientY011 = [1,0,0]; orientZ011 = [0,1,1]
orientX111 = [-1,-1,2]; orientY111 = [1,-1,0]; orientZ111 = [1,1,1]
orientX310 = [-1,3,0]; orientY310 = [0,0,1]; orientZ310 = [3,1,0]

orientsX = []
orientsY = []
orientsZ = []

if '001' in orients:
    orientsX.append(orientX001)
    orientsY.append(orientY001)
    orientsZ.append(orientZ001)
if '011' in orients:
    orientsX.append(orientX011)
    orientsY.append(orientY011)
    orientsZ.append(orientZ011)
if '111' in orients:
    orientsX.append(orientX111)
    orientsY.append(orientY111)
    orientsZ.append(orientZ111)
if '310' in orients:
    orientsX.append(orientX310)
    orientsY.append(orientY310)
    orientsZ.append(orientZ310)

def threshAtomsSphere(r):
    return int(math.ceil(1.335*math.pi*r**2)) #CHECK this empirical approximation for 9ao bubbles!!!

for runNum in numRun: 
    for temp in t:
        for oX in orientsX:
            index = orientsX.index(oX)
            oY = orientsY[index]
            oZ = orientsZ[index]
            crystOrient = str(oZ[0])+str(oZ[1])+str(oZ[2])
            if singleLayer == True:
                appendor = []
                newl = []
                appendor.append(minThickness[crystOrient])
                newl = appendor+l
            else:
                newl = l
            for sphR in Rs:
                for shp in bubShape:
                    for lig in newl:
                        if lig <= 4:
                            if lig <= 2:
                                PeqAdj = 0.6
                            else: 
                                PeqAdj = 0.9
                        else: 
                            PeqAdj = 1
                        fileline = []
                        for line in open(LatParamFile):
                            if re.match('^'+str(temp)+'.*', line):
                                fileline.append(str(line))
                        fileline = fileline[0].split()
                        ao = float(fileline[1]) #Temperature-specific lattice parameter in Angstroms.    
                        fileID = str(runNum) + '_' + str(temp)+'K_'+crystOrient+'_Rs'+str(sphR)+'ao_'+shp+'_lig'+str(lig)+'ao' #Unique name based on variables.
                        if compFac:
                            fileID = fileID + '_comp' + str(compFac)
                        # Put file in appropriate directory. 
                        parent = os.getcwd()
                        subdir = parent + '/run_' + fileID 
                        try:
                            ################################################################                        
                            currentJobs = int(os.popen('squeue -u $USER | wc -l').read())-1
                            if (currentJobs + len(fileTypes)) > maxJobs:
                                print('Your queue is now full.\n')
                                quit()
                            ################################################################
                            os.mkdir(subdir)
                            os.chdir(subdir)
                            for fileType in fileTypes: 
                                # Write the LAMMPS input file.
                                fileTitle = fileType+'.'+fileID
                                with open(fileTitle,'w+') as f: 
                                    if fileType == 'in':
                                        f.write('# Initialization\n')
                                        f.write('units metal \n')
                                        f.write('dimension 3 \n')
                                        f.write('timestep {:}\n'.format(ts))
                                        f.write('neigh_modify every 1 delay 0 check yes\n')
                                        f.write('boundary p p f \n\n')
                                    if fileType != 'in':
                                        f.write('print "RESTART BEGINS HERE."\n')
                                        f.write('read_restart restart_{:}.*{:}\n\n'.format(fileID,discrim))
                                        f.write('#Specify things that the restart file may not contain\n\n')
                                    f.write('# Initial Setup\n')
                                    f.write('variable latPar equal {:}\n'.format(ao))
                                    f.write('variable tempSet equal {:}\n'.format(temp))
                                    f.write('variable sphereRadius equal {:}\n'.format(sphR))
                                    f.write('variable ligament equal {:}\n'.format(lig))
                                    f.write('variable lig equal ${ligament}*${latPar}\n')
                                    f.write('variable sphR equal ${sphereRadius}*${latPar}\n')
                                    f.write('variable volRat equal {:}\n'.format(volRat))
                                    f.write('variable botLigRat equal {:}\n'.format(botLigRat))
                                    f.write('variable sideDepthRat equal {:}\n'.format(sideDepthRat))
                                    f.write('variable sIdeal equal ${sphR}*(4*PI/3/${volRat})^(1/3)\n')
                                    if fileType == 'in':
                                        f.write('lattice bcc ${{latPar}} origin 0 0 0 orient x {:} {:} {:} orient y {:} {:} {:} orient z {:} {:} {:}\n'.format(oX[0],oX[1],oX[2],oY[0],oY[1],oY[2],oZ[0],oZ[1],oZ[2]))
                                    f.write('variable xSpacing equal {:}*${{latPar}}\n'.format(xSpacing[crystOrient]))
                                    f.write('variable ySpacing equal {:}*${{latPar}}\n'.format(ySpacing[crystOrient]))
                                    f.write('variable zSpacing equal {:}*${{latPar}}\n'.format(zSpacing[crystOrient]))
                                    f.write('variable veryminx equal {:}*${{xSpacing}}\n'.format(verymin))
                                    f.write('variable veryminy equal {:}*${{ySpacing}}\n'.format(verymin))
                                    f.write('variable veryminz equal {:}*${{zSpacing}}\n'.format(verymin))
                                    f.write('variable minSide equal ${sideDepthRat}*(${lig}+${sphR})\n')
                                    f.write('variable xSideinit equal ceil(${sIdeal}/${xSpacing})*${xSpacing}\n')
                                    f.write('variable ySideinit equal ceil(${sIdeal}/${ySpacing})*${ySpacing}\n')
                                    f.write("""if "${xSideinit} >= ${veryminx}" then &\n""")
                                    f.write("""\t"variable xSide equal ${xSideinit}" &\n""")
                                    f.write("""else &\n""")
                                    f.write("""\t"variable xSide equal ${veryminx}"\n""")
                                    f.write("""if "${ySideinit} >= ${veryminy}" then &\n""")
                                    f.write("""\t"variable ySide equal ${ySideinit}" &\n""")
                                    f.write("""else &\n""")
                                    f.write("""\t"variable ySide equal ${veryminy}"\n""")
                                    f.write("""if "${xSide} < ${ySide}" then &\n""")
                                    f.write("""\t"variable sideComp equal ${xSide}" &\n""")
                                    f.write("""else &\n""")
                                    f.write("""\t"variable sideComp equal ${ySide}"\n""")
                                    f.write("""if "${sideComp} >= ${minSide}" then &\n""")
                                    f.write("""\t"variable almostBlockSideX equal ${xSide}" &\n""")
                                    f.write("""\t"variable almostBlockSideY equal ${ySide}" &\n""")
                                    f.write("""else &\n""")
                                    f.write("""\t"variable almostBlockSideX equal ceil(${minSide}/${xSpacing})*${xSpacing}" &\n""")
                                    f.write("""\t"variable almostBlockSideY equal ceil(${minSide}/${ySpacing})*${ySpacing}"\n""")
                                    f.write('variable xUnits equal ${almostBlockSideX}/${xSpacing}\n')
                                    f.write('variable xremainder equal ${{xUnits}}%{:}\n'.format(xrepeat[crystOrient]))
                                    f.write("""if "${xremainder} == 0" then &\n""")
                                    f.write("""\t"variable xmultiplier equal ${xUnits}" &\n""")
                                    f.write('else &\n')
                                    f.write("""\t"variable xmultiplier equal ${{xUnits}}+{:}-${{xremainder}}"\n""".format(xrepeat[crystOrient]))
                                    f.write('variable blockSideX equal ${xmultiplier}*${xSpacing}\n') 
                                    f.write('variable yUnits equal ${almostBlockSideY}/${ySpacing}\n')
                                    f.write('variable yremainder equal ${{yUnits}}%{:}\n'.format(yrepeat[crystOrient]))
                                    f.write("""if "${yremainder} == 0" then &\n""")
                                    f.write("""\t"variable ymultiplier equal ${yUnits}" &\n""")
                                    f.write("""else &\n""")
                                    f.write("""\t"variable ymultiplier equal ${{yUnits}}+{:}-${{yremainder}}"\n""".format(yrepeat[crystOrient]))
                                    f.write('variable blockSideY equal ${ymultiplier}*${ySpacing}\n')                                 
                                    f.write('variable bottom equal ${sIdeal}-${lig}-2*${sphR}\n')
                                    f.write('variable minBottom equal ${botLigRat}*${lig}\n')
                                    f.write('if "${bottom} >= ${minBottom}" then &\n')
                                    f.write('\t"variable zSide equal ceil(${sIdeal}/${zSpacing})*${zSpacing}" &\n')
                                    f.write('else &\n')
                                    f.write('\t"variable sZadjusted equal 2*${sphR}+${lig}+${botLigRat}*${lig}" &\n')
                                    f.write('\t"variable zSide equal ceil(${sZadjusted}/${zSpacing})*${zSpacing}"\n')
                                    f.write("""if "${zSide} >= ${veryminz}" then &\n""")
                                    f.write("""\t"variable almostBlockSideZ equal ${zSide}" &\n""")
                                    f.write("""else &\n""")
                                    f.write("""\t"variable almostBlockSideZ equal ${veryminz}"\n""")
                                    f.write('variable zUnits equal ${almostBlockSideZ}/${zSpacing}\n')
                                    f.write('variable zremainder equal ${{zUnits}}%{:}\n'.format(zrepeat[crystOrient]))
                                    f.write("""if "${zremainder} == 0" then &\n""")
                                    f.write("""\t"variable zmultiplier equal ${zUnits}" &\n""")
                                    f.write("""else &\n""")
                                    f.write("""\t"variable zmultiplier equal ${{zUnits}}+{:}-${{zremainder}}"\n""".format(zrepeat[crystOrient]))
                                    f.write('variable blockSideZ equal ${zmultiplier}*${zSpacing}+0.1\n')
                                    f.write('variable boxLow equal -1/4*${blockSideZ}\n')
                                    f.write('variable boxHigh equal 1.25*${blockSideZ}\n')
                                    f.write('region whole block 0 ${blockSideX} 0 ${blockSideY} ${boxLow} ${boxHigh} units box\n')
                                    if fileType == 'in':
                                        f.write('create_box 2 whole\n')
                                    f.write('region Wcube block INF INF INF INF -0.1 ${blockSideZ} units box\n')
                                    if fileType == 'in':
                                        f.write('create_atoms 1 region Wcube\n')
                                        f.write('group beforeDel type 1\n')
                                        f.write('variable countBefore equal count(beforeDel)\n')
                                        f.write('variable countBeforeUse equal ${countBefore}\n\n')
                                    if fileType != 'in':
                                        f.write('\n')
             
                                    f.write('# Further Setup\n')
                                    f.write('pair_style hybrid eam/fs table linear 10000\n')
                                    f.write('pair_coeff * * eam/fs {:}/W_Juslin2010_AT_mod.eam.fs W NULL\n'.format(pots))
                                    f.write('pair_coeff 1 2 table {:}/W-He-Juslin2013.table WHe\n'.format(pots))
                                    f.write('pair_coeff 2 2 table {:}/He-Beck1968_modified.table HeHe\n'.format(pots))
                                    f.write('mass 1 {:}\n'.format(m1))
                                    f.write('mass 2 {:}\n'.format(m2))
                                    f.write('newton on\n')
                                    f.write('fix bal all balance 20 1.02 shift z 30 1.02\n') #ASK if okay, could slow simulation a lot.
                                    if fileType == 'in':
                                        f.write('minimize 0 0 {:} {:}\n'.format(WChunkMinIter,WChunkMinEvals))
                                        f.write('velocity all create ${{tempSet}} {:} dist gaussian\n'.format(random.randint(1,99999)))
                                    f.write('fix thermostat all nvt temp ${{tempSet}} ${{tempSet}} $({:}*dt)\n'.format(statDamp))
                                    if fileType == 'in':
                                        f.write('run {:}\n'.format(int(math.ceil(1.5*stepsToEq))))
                                        Xpos = random.uniform(xymod*ao/10,xymod*ao)
                                        Xneg = random.uniform(-1*xymod*ao,-1*xymod*ao/10)
                                        Ypos = random.uniform(xymod*ao/10,xymod*ao)
                                        Yneg = random.uniform(-1*xymod*ao,-1*xymod*ao/10)
                                        Xshift = random.choice((Xpos,Xneg))
                                        Yshift = random.choice((Ypos,Yneg))
                                    if Xshift > 0:
                                        f.write('variable SphCentX equal ${{blockSideX}}/2+{:}\n'.format(Xshift))
                                    else:
                                        f.write('variable SphCentX equal ${{blockSideX}}/2{:}\n'.format(Xshift))
                                    if Yshift > 0:    
                                        f.write('variable SphCentY equal ${{blockSideY}}/2+{:}\n'.format(Yshift))
                                    else: 
                                        f.write('variable SphCentY equal ${{blockSideY}}/2{:}\n'.format(Yshift))
                                    f.write('variable tippity equal bound(all,zmax)\n')
                                    f.write('variable SphCentZ equal ${tippity}-${sphR}-${lig}\n')
                                    f.write('variable pressR equal ${sphR}*0.8\n')
                                    f.write('variable pressVol equal 4/3*PI*(${pressR}^3)\n')
                                    if shp == 'sph':
                                        f.write('region bubble sphere ${SphCentX} ${SphCentY} ${SphCentZ} ${sphR} units box\n')
                                        f.write('region pressReg sphere ${SphCentX} ${SphCentY} ${SphCentZ} ${pressR} units box\n') #FIXME your region might get messed up after rescaling, check!
                                    #elif shp == 'ell':
                                        #Fix this section only if necessary. 
                                        #Rh = sphR*3**(1/3) #Ellipsoid horizontal radii in lattice units.
                                        #Rv = sphR/3 #Ellipsoid vertical radius in lattice units. 
                                        #f.write('region bubble ellipsoid {:} {:} {:} {:} {:} {:} units lattice\n'.format(EllCentX,EllCenty,EllCentZ,Rh,Rh,Rv))
                                    f.write('group pressRegAtoms dynamic all region pressReg every 1\n')
                                    if fileType == 'in':
                                        f.write('delete_atoms region bubble\n')
                                        f.write('group afterDel type 1\n')
                                        f.write('variable countAfter equal count(afterDel)\n')
                                        f.write('variable countAfterUse equal ${countAfter}\n')
                                        f.write('variable deleted equal ${countBeforeUse}-${countAfterUse}\n') 
                                        f.write('print "nv: ${deleted}"\n')
                                    f.write('compute stresses pressRegAtoms stress/atom NULL\n') #Finds per-atom stress tensor in relevant region. 
                                    f.write('compute diagSums pressRegAtoms reduce sum c_stresses[1] c_stresses[2] c_stresses[3]\n') #Sums xx, yy, and zz components (respectively) for all atoms in system. Outputs 1x3 vector.
                                    f.write('variable bubPress equal -(c_diagSums[1]+c_diagSums[2]+c_diagSums[3])/(3*${pressVol})\n') #Computes pressure within bubble. '3' is present since system is 3D. 
                                    f.write('compute Wstress afterDel stress/atom NULL\n\n')
                                    
                                    # Use equation of state (Eq 25) given in "Theoretical Model of Helium Buble Growth and Density in Plasma-Facing Metals" by Hammond, Maroudas, and Wirth. This equation is a simplification of a more complicated one, but should still give solid results for our purposes here.
                                    
                                    if fileType == 'in':
                                        f.write('# Fill Bubble Initially\n')
                                        #Values that you either have or MAYBE SHOULD calculate. CHECK.
                                        f.write('variable ao equal {:}E-10\n'.format(ao)) #m (lattice parameter).
                                        f.write('variable gamma equal 2.67\n') #J/m^2 (surface tension). REPLACE with your own number (maybe). Also adjust for temperature (probably not).
                                        f.write('variable omega equal (v_ao^3)/2\n') #m^3 (atomic volume). This formula is given in the paper above and is explained well at https://en.wikipedia.org/wiki/Atomic_packing_factor#Body-centered_cubic.
                                        f.write('variable G equal 156.1E9\n') #Pa (shear modulus). REPLACE with your own number (maybe).
                                        f.write('variable b equal v_ao*(3^(1/2))/2\n') #m (magnitude of dislocation loop Burgers vector). The formula used was ao/2*|<111>|. I THINK this is correct in general since <111> is the closest-packed direction. ASK.
                                        f.write('variable k equal 1.38064852E-23\n') #J/K (Boltzmann constant)
                                        f.write('variable Fs equal (3*v_omega/4/PI)^(1/3)\n') #'geometric parameter for spheres' that makes the formulas more graceful
                                        #Values from the cited paper fit at 933 K. Will introduce some error, so BE READY TO MAKE ADJUSTMENTS. Underbars should be read as negative signs for A's.
                                        f.write('variable A_11 equal -5.5991E-27\n') #m^3*K^(1/2)*Pa^(1/3)
                                        f.write('variable A01 equal 1.7400E-26\n') #m^3*Pa^(1/3)
                                        f.write('variable A21 equal 4.9833E-30\n') #m^3/K*Pa^(1/3)
                                        f.write('variable A02 equal -4.4658E-24\n') #m^3*Pa^(2/3)
                                        f.write('variable A22 equal -8.7825E-27\n') #m^3/K*Pa^(2/3)
                                        f.write('variable A03 equal 1.7595E-22\n') #J 
                                        f.write('variable A23 equal 1.7608E-23\n') #J/K
                                        f.write('variable A_13 equal -3.2615E-21\n') #J*K^(1/2)
                                        f.write('variable A_23 equal 3.1524E-20\n') #J*K
                                        #The function evaluation begins here.
                                        f.write('variable f1 equal v_A_11*(${tempSet}^(-1/2))+v_A01+v_A21*${tempSet}\n')
                                        f.write('variable f2 equal v_A02+v_A22*${tempSet}\n')
                                        f.write('variable f3 equal v_A_23/${tempSet}+v_A_13*(${tempSet}^(-1/2))+v_A03+v_A23*${tempSet}\n')
                                        f.write('variable par equal 4*v_gamma^2+2*v_gamma*v_G*v_b\n')
                                        f.write('variable t1 equal v_f1*v_Fs^(1/3)*v_par^(-1/6)*v_deleted^(1/9)\n')
                                        f.write('variable t2 equal v_f2*v_Fs^(2/3)*v_par^(-1/3)*v_deleted^(2/9)\n')
                                        f.write('variable t3 equal v_f3*v_Fs*v_par^(-1/2)*v_deleted^(1/3)\n')
                                        f.write('variable brac equal v_t1+v_t2+v_t3\n')
                                        f.write('variable nHe equal floor({:}*v_deleted*v_omega/v_brac)\n'.format(PeqAdj))
                                        f.write('print "Helium atoms to be inserted: ${nHe}"\n')
                                        #f.write('variable maxHe equal {:}*${{nHe}}\n'.format(maxHe))
                                        f.write('variable HeLatPar equal ${{latPar}}*{:}\n'.format(HeLatParMultip))
                                        f.write('lattice fcc ${HeLatPar}\n')
                                        f.write('create_atoms 2 region bubble\n')
                                        f.write('group He type 2\n')
                                        f.write('variable HeCount equal count(He)\n')
                                        f.write('variable HeCountUse equal ${HeCount}\n')
                                        f.write('print "Created ${HeCountUse} fcc helium atoms of the target ${nHe} helium atoms."\n')
                                        f.write('group He clear\n')
                                        f.write('minimize 0 0 {:} {:}\n'.format(bubMaxIter,bubMaxEvals)) #Might not need this (or can make it much shorter) if the minimization inside the loop works well enough.
                                        f.write('variable HeDispMax equal 0.95/sqrt(3)*${sphR}\n')
                                        f.write('variable HeDispMin equal -1*${HeDispMax}\n')
                                        f.write('variable loopVar loop ${HeCountUse} ${nHe}\n')
                                        f.write('variable rando1 equal random(${{HeDispMin}},${{HeDispMax}},{:})\n'.format(random.randint(1,99999)))
                                        f.write('variable rando2 equal random(${{HeDispMin}},${{HeDispMax}},{:})\n'.format(random.randint(1,99999)))
                                        f.write('variable rando3 equal random(${{HeDispMin}},${{HeDispMax}},{:})\n'.format(random.randint(1,99999)))
                                        f.write('group W type 1\n')
                                        f.write('fix chilly W setforce 0 0 0\n')
                                        f.write('label HeFill\n')
                                        f.write('variable posx equal ${SphCentX}+v_rando1\n')
                                        f.write('variable posy equal ${SphCentY}+v_rando2\n')
                                        f.write('variable posz equal ${SphCentZ}+v_rando3\n')
                                        f.write('create_atoms 2 single ${posx} ${posy} ${posz} units box\n')
                                        f.write('minimize 0 0 {:} {:}\n'.format(HeInitPlaceIter,HeInitPlaceEvals))
                                        f.write('next loopVar\n')
                                        f.write('print "Inserted ${loopVar} of ${nHe} He atoms up to this point."\n')
                                        f.write('if "${loopVar} < ${nHe}" then "jump SELF HeFill" else "jump SELF continue"\n')
                                        f.write('label continue\n') 
                                        f.write('unfix chilly\n')
                                        f.write('minimize 0 0 {:} {:}\n'.format(bubMaxIter,bubMaxEvals)) #Might not need this (or can make it much shorter) if the minimization inside the loop works well enough.
                                        f.write('group W clear\n\n')
                                    
                                    f.write('# Final System Prep\n')
                                    f.write('group W type 1\n')
                                    f.write('group He type 2\n')
                                    f.write('variable countHe equal count(He)\n')
                                    f.write('fix linearP W momentum 1 linear 1 1 1\n')
                                    if fileType == 'in':
                                        if compFac:
                                            f.write('fix rescale all deform 1 x scale {:} y scale {:} remap x\n'.format(compFac,compFac))
                                            f.write('run {:}\n'.format(int(math.ceil(2*stepsToEq))))
                                            f.write('unfix rescale\n')
                                        f.write('run {:}\n'.format(int(math.ceil(2*stepsToEq))))
                                    f.write('variable topW equal bound(all,zmax)\n')
                                    if fileType == 'in':    
                                        f.write('variable punchyMin equal ${{topW}}+{:}*${{latPar}}\n'.format(punchMultip))
                                    if fileType != 'in': 
                                        f.write('variable punchyMin equal ${{topW}}+{:}*{:}*${{latPar}}\n'.format(punchMultipAdj,punchMultip)) #NOTE the adjusted 'catch' threshold. This may need to be adjusted. 
                                    f.write('region reliefReg block 0 ${blockSideX} 0 ${blockSideY} ${punchyMin} ${boxHigh} units box\n')
                                    f.write('group punchyGroup dynamic W region reliefReg every 1\n')
                                    f.write('variable bottomW equal bound(all,zmin)\n')
                                    if fileType == 'in':    
                                        f.write('variable punchyBottom equal ${{bottomW}}-{:}*${{latPar}}\n'.format(punchMultip))
                                    if fileType != 'in': 
                                        f.write('variable punchyBottom equal ${{bottomW}}-{:}*{:}*${{latPar}}\n'.format(punchMultipAdj,punchMultip)) #NOTE the adjusted 'catch' threshold. This may need to be adjusted. 
                                    f.write('region bottomReg block 0 ${blockSideX} 0 ${blockSideY} ${boxLow} ${punchyBottom} units box\n')
                                    f.write('group bottomGroup dynamic W region bottomReg every 1\n\n')

                                    f.write('# Outputs\n')
                                    f.write('thermo_style custom step temp press etotal v_countHe v_bubPress\n') 
                                    f.write('thermo {:}\n'.format(betweenThermoTS))
                                    f.write('thermo_modify lost ignore flush yes\n')
                                    f.write('dump atomDump all custom {:} dump_{:}.* id type xs ys zs c_Wstress[1] c_Wstress[2] c_Wstress[3] c_Wstress[4] c_Wstress[5] c_Wstress[6]\n'.format(betweenDumpsTS,fileID))  
                                    f.write('dump_modify atomDump pad 8\n')
                                    f.write('restart {:} restart_{:}\n'.format(betweenRestartsTS,fileID))
                                    if fileType == 'in':
                                        f.write('log log_{:}\n\n'.format(fileID))   
                                    if fileType != 'in':
                                        f.write('log log_{:} append\n\n'.format(fileID))    

                                    f.write('# Actual Run\n')
                                    if fileType == 'in':
                                        f.write('print "nv: ${deleted}"\n')
                                        f.write('reset_timestep 0\n')
                                    f.write('variable loopVar2 loop 1 {:}\n'.format(numberOfInsertions)) 
                                    f.write('label runLoop\n') #This begins the insert/run loop.
                                    f.write('undump atomDump\n')
                                    f.write('group almostAll region whole\n')
                                    f.write('variable timestep equal step\n')
                                    f.write('variable current_ts equal ${timestep}\n')
                                    f.write('fix freeze almostAll setforce 0 0 0\n')
                                    f.write('create_atoms 2 single ${SphCentX} ${SphCentY} ${SphCentZ} units box\n') #FIXME changing the box might make atoms insert outside of bubble, be careful!
                                    f.write('minimize 0 0 {:} {:}\n'.format(HeMinIter,HeMinEvals))
                                    f.write('unfix freeze\n')
                                    f.write('print "Inserted ${loopVar2} He so far during the run phase."\n')
                                    f.write('variable countAll equal count(all)\n')
                                    f.write('variable countAllBefore equal ${countAll}\n')
                                    f.write('variable countAll delete\n')
                                    f.write('group He clear\n')
                                    f.write('variable countHe delete\n')
                                    f.write('group He type 2\n')
                                    f.write('variable countHe equal count(He)\n')
                                    f.write('next loopVar2\n')
                                    f.write('group almostAll clear\n')
                                    f.write('reset_timestep ${current_ts}\n') 
                                    f.write('dump atomDump all custom {:} dump_{:}.* id type xs ys zs c_Wstress[1] c_Wstress[2] c_Wstress[3] c_Wstress[4] c_Wstress[5] c_Wstress[6]\n'.format(betweenDumpsTS,fileID))  
                                    f.write('dump_modify atomDump pad 8\n')
                                    f.write('run {:}\n'.format(betweenInsertionsTS))
                                    f.write('variable countAll equal count(all)\n')
                                    f.write('variable countAllAfter equal ${countAll}\n')
                                    f.write('variable countAll delete\n')
                                    f.write("""if "${countAllAfter} < ${countAllBefore}" then &\n""")
                                    f.write("""\t"print 'Apparent Relief Method: Bursting.'" &\n""")
                                    f.write("""\t"undump atomDump" &\n""")
                                    f.write("""\t"dump atomDump all custom {:} dump_{:}.* id type xs ys zs c_Wstress[1] c_Wstress[2] c_Wstress[3] c_Wstress[4] c_Wstress[5] c_Wstress[6]" &\n""".format(betweenDumpsTSFast,fileID))  
                                    f.write("""\t"dump_modify atomDump pad 8" &\n""")
                                    f.write("""\t"run {:}" &\n""".format(betweenInsertionsTS))
                                    f.write("""\t"quit 0"\n""")
                                    f.write('variable countallBefore delete\n')
                                    f.write('variable countallAfter delete\n')
                                    f.write('variable countWBottom equal count(bottomGroup)\n')
                                    f.write('print "Number of tungsten atoms in bottom region: ${countWBottom}"\n')
                                    f.write("""if "${{countWBottom}} > {:}" then &\n""".format(threshAtomsSphere(sphR)))
                                    f.write("""\t"print 'Apparent Relief Method: Punching Through Bottom.'" &\n""")
                                    f.write("""\t"undump atomDump" &\n""")
                                    f.write("""\t"dump atomDump all custom {:} dump_{:}.* id type xs ys zs c_Wstress[1] c_Wstress[2] c_Wstress[3] c_Wstress[4] c_Wstress[5] c_Wstress[6]" &\n""".format(betweenDumpsTSFast,fileID))  
                                    f.write("""\t"dump_modify atomDump pad 8" &\n""")
                                    f.write("""\t"run {:}" &\n""".format(betweenInsertionsTS))
                                    f.write("""\t"quit 0"\n""")
                                    f.write('variable countW equal count(punchyGroup)\n')
                                    f.write('print "Number of tungsten atoms in loop punch detection region: ${countW}"\n')
                                    f.write("""if "${{countW}} > {:}" then "jump SELF lastLoop"\n""".format(threshAtomsSphere(sphR)))
                                    #f.write("""if "${countHe} >= ${maxHe}" then &\n""")
                                    #f.write("""\t"print 'Many heliums inserted - loop surely floating around.'" &\n""")
                                    #f.write("""\t"print 'Apparent Relief Method: Punching Through Bottom.'" &\n""")
                                    #f.write("""\t"quit 0"\n""")
                                    f.write("""if "${{timestep}} < {:}" then &\n""".format(fullRunTS))
                                    f.write("""\t"jump SELF runLoop" &\n""")
                                    f.write("""\telse &\n""")
                                    f.write("""\t"print 'Time Limit Reached.'" &\n""")
                                    f.write("""\t"quit 0"\n""")
                                    f.write('label lastLoop\n')
                                    f.write('print "Entering final run sequence to determine if loop punching or bursting occured."\n')
                                    f.write('variable countAll equal count(all)\n')
                                    f.write('variable countAllBefore equal ${countAll}\n')
                                    f.write('variable countAll delete\n')
                                    f.write('undump atomDump\n')
                                    f.write('dump atomDump all custom {:} dump_{:}.* id type xs ys zs c_Wstress[1] c_Wstress[2] c_Wstress[3] c_Wstress[4] c_Wstress[5] c_Wstress[6]\n'.format(betweenDumpsTSFast,fileID))  
                                    f.write('dump_modify atomDump pad 8\n')
                                    f.write('run {:}\n'.format(betweenInsertionsTS))
                                    f.write('variable countAll equal count(all)\n')
                                    f.write('variable countAllAfter equal ${countAll}\n')
                                    f.write('variable countAll delete\n')
                                    f.write("""if "${countAllAfter} < ${countAllBefore}" then &\n""")
                                    f.write("""\t"print 'Apparent Relief Method: Bursting.'" &\n""")
                                    f.write("""\t"quit 0"\n""")
                                    f.write('print "Apparent Relief Method: Loop Punching."\n')
                                    f.write('quit 0') 

                                    # Create slurm files. 
                                    jobname = fileType + '_' + fileID
                                    filename = jobname + '.slurm'
                                    with open(filename, 'w') as f:
                                        f.write('#!/bin/bash\n')
                                        f.write('#SBATCH -p Lewis\n')
                                        f.write('#SBATCH -N 1\n')
                                        f.write('#SBATCH -J ' + jobname + '\n')
                                        f.write('#SBATCH -A hammond\n')
                                        if sphR <= 1 and lig <= 20:
                                            f.write('#SBATCH -t {:}\n'.format(runTime['short']))
                                        elif sphR <= 3 and lig <= 10:
                                            f.write('#SBATCH -t {:}\n'.format(runTime['short']))
                                        elif sphR <= 3 and lig <= 20:
                                            f.write('#SBATCH -t {:}\n'.format(runTime['med']))
                                        else:
                                            f.write('#SBATCH -t {:}\n'.format(runTime['long']))
                                        f.write('#SBATCH -n 24\n')
                                        f.write('#SBATCH --mem-per-cpu=8G\n')
                                        f.write('#SBATCH --mail-user=bfl3m8@mail.missouri.edu\n')
                                        f.write('#SBATCH --mail-type=ALL\n')
                                        f.write('#SBATCH -o ' + jobname + '.out\n')
                                        f.write('#SBATCH -e ' + jobname + '.error\n')
                                        if fileType != 'in':
                                            f.write('#SBATCH --kill-on-invalid-dep=yes\n\n')
                                        else:
                                            f.write('\n')

                                        f.write('unset MODULEPATH\n') #FIXME: will these bits still work on Lewis?
                                        f.write('source /cluster/spack-2020/opt/spack/linux-centos7-x86_64/gcc-6.3.0/lmod-8.1.5-w5bf3f6lvfzxk6vmeaagqkjof47zl3yn/lmod/lmod/init/bash\n')
                                        f.write('module use /cluster/spack-2020/share/spack/lmod/linux-centos7-x86_64/Core\n')		
                                        f.write('module load intel/19.0.1-x86_64\n')
                                        f.write('module load openmpi\n\n')

                                        f.write('#Remove unnecessary files:\n')
                                        f.write("find -name 'restart*' -type f ! -name 'restart*{:}' -delete\n\n".format(discrim))
                                        
                                        f.write('#Run the program:\n')
                                        f.write('srun /home/bfl3m8/bin/lmp_mpi -in {:}\n\n'.format(fileTitle))
                                        
                                        f.write('#Remove unnecessary files:\n')
                                        f.write("find -name 'restart*' -type f ! -name 'restart*{:}' -delete".format(discrim))
                                    # Submit jobs to LEWIS if instructed to do so. 
                                    if args.submit:
                                        if fileType == 'in':
                                            command = 'sbatch ' + filename
                                            p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
                                            (output, err) = p.communicate()
                                            p.wait()
                                            text = output.decode('ascii')
                                            job_id = text.split(' ')
                                            job_id = job_id[-1].strip()
                                        if fileType != 'in':
                                            command = 'sbatch -d afternotok:' + job_id + ' ' + filename 
                                            p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
                                            (output, err) = p.communicate()
                                            p.wait()
                                            text = output.decode('ascii')
                                            job_id = text.split(' ')
                                            job_id = job_id[-1].strip()
                            # Move back into parent directory
                            os.chdir(parent)
                        except FileExistsError:
                            continue
print('All requested runs have been created.\n')
