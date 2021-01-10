import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import argparse 
import pickle as pkl
from itertools import cycle

minbulklig = 10 #ao
meancol = 'k'
datacol = 'grey'
linecol = 'k'
opacity = 1 #for grey data points
tickstep = 20
markerArea = 8 #Squared pixels

plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['mathtext.it'] = 'serif:italic'
plt.rcParams['mathtext.bf'] = 'serif:bold'
plt.rcParams['lines.linewidth'] = 0.8
plt.rcParams.update({'font.size':10}) 

ao = {'500':3.17,'1000':3.183,'1200':3.189,'1500':3.199,'2000':3.214,'2500':3.229} #'temperature(K)':latPar(Angs)

#-------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-f','--filename', type = str, nargs = 1, help = 'Input dataframe pickle to read from.', required = True)
parser.add_argument('-r','--radius', type = float, nargs = '*', help = 'Optional filter operation. Include only datapoints for bubbles whose radii match the given radius/radii in the analysis. Give the radius in units of lattice parameters (e.g. for a r = 3 ao bubble, type "3").', required = False)
parser.add_argument('-t','--temp', type = float, nargs = '*', help = 'Optional filter operation. Include only datapoints for systems whose temperatures match the given temperature(s) in the analysis. Give the temperature in Kelvin, but do not type "K".', required = False)
parser.add_argument('-d','--display', action = 'store_true', default = False, help = 'Display the produced graph on the screen.', required = False)
parser.add_argument('-i','--image', type = str, nargs = '*', help = 'Save the produced graph with the given name. INCLUDE a file extension (png or pdf). If you used the -c flag, you must supply three names. Otherwise, supply one name.', required = False) 
parser.add_argument('--dualcol', action = 'store_true', default = False, help = 'Optional. Toggle on to split the figure legend into two columns.', required = False)
parser.add_argument('--figheight', type = float, nargs = 1, default = 3, help = 'Optional. Specify the figure height in inches. Defaults to 10. Note that exact height may vary due to screen DPI.', required = False)
parser.add_argument('--figwidth', type = float, nargs = 1, default = 7, help = 'Optional. Specify the figure width in inches. Defaults to 18. Note that exact width may vary due to screen DPI.', required = False)
parser.add_argument('--color', action = 'store_true', default = False, help = 'Optional. Make the plots in color rather than black and white.', required = False)
args = parser.parse_args()
#------------------------------------------------

if isinstance(args.figheight,int):
    figheight = args.figheight
else:
    figheight = args.figheight[0]
if isinstance(args.figwidth,int):
    figwidth = args.figwidth
else:
    figwidth = args.figwidth[0]

if args.image:
    assert len(args.image)==1, 'Single images must be given a single name.'

def magic(rad,temp):
    A_P = 99.33750903
    B_P = 0.664284758
    C_P = 0.000240928
    Pb = A_P*rad**(-B_P)*np.exp(-C_P*temp)
    return Pb

def Varshni(temp):
    'This is equation 4 from https://doi.org/10.1103%2FPhysRevB.2.3952. For us, C44 = G.'
    C440 = 16.31577E11 #dyn/cm^2
    s44 = 0.568995E11 #dyn/cm^2
    t44 = 346 #K
    G = C440-s44/(np.exp(t44/temp)-1) #dyn/cm^2
    conversion_factor = 1E-10
    return G*conversion_factor #GPa

def cont(rad,temp):
    'This is the "Trinkhaus equation", Pb=(2*(gamma)+G*b)/R. Use constant value for gamma from https://doi.org/10.1038/s41598-020-58581-8.'
    gamma = 2.67 #J/m^2
    G = Varshni(temp)*10**9 #Pa
    latpar = ao[str(int(temp))]/10 #nm
    b = latpar*np.sqrt(3)/2/10**9 #m
    R = rad*latpar/10**9 #m
    conversion_factor = 10**9
    return (2*gamma+G*b)/R/conversion_factor #GPa

def boring(x):
    T = tuple(x)
    if len(T) > 1:
        return T
    else:
        return T[0]

df = pd.read_pickle(args.filename[0])

a = pd.DataFrame(columns = df.columns) #Create a temporary blank dataframe. 
for index,series in df.iterrows():
    truth = {'radius':False,'temp':False,'lig':False}
    if args.radius:
        if series['Rs(ao)'] in args.radius:
            truth['radius'] = True
    else:
        truth['radius'] = True
    if args.temp:
        if series['temp(K)'] in args.temp:
            truth['temp'] = True
    else:
        truth['temp'] = True
    if series['lig(ao)'] >= minbulklig:
        truth['lig'] = True
    if all(value == True for value in truth.values()):
        a = a.append(series)
    truth = None

a['reliefType'] = a['reliefType'].map({'ERROR':'ERROR', 'loopPunch':'loopPunch', 'bottomPunch':'bottomPunch', 'burst':'burst', 'timeOut':'bottomPunch'})

df = a
a = None

aggIn = 'maxPress(bar)' 
ylab = r'Relief Pressure (GPa)' 

moreGrouping = ['temp(K)','Rs(ao)'] 
dfMoreGrouped = df.groupby(by=moreGrouping)

raddf = dfMoreGrouped[aggIn].agg(boring) #The mean is just a dummy function

temps = raddf.unstack().index
radii = raddf.unstack().columns
pressures = raddf.unstack().values

colors = cm.rainbow(np.linspace(0,1,len(temps)))

symlist = cycle(['d','s','^','o','v','*','2','X'])

################################################ Plotting starts here. 

radrange = np.linspace(min(radii),max(radii),1000)
maxP = []
plt.figure(figsize=(figwidth,figheight))
for tempind in range(len(temps)):
    symplot = next(symlist)
    if args.color:   
        col = colors[tempind]
        linecol = col
        datacol = col
        meancol = col
    temp = temps[tempind]
    label = str(int(temp))+' K'
    for radind in range(len(radii)):
        rad = radii[radind]
        pseries = pressures[tempind][radind]
        if not isinstance(pseries,tuple):
            continue
        pseries = [i/10000 for i in pseries] #converts bar to GPa
        avgP = np.mean(pseries)
        maxP.append(max(pseries))
        radseries = [rad]*len(pseries)
        plt.scatter(radseries,pseries,color=datacol,facecolors='none',marker=symplot,s=markerArea,alpha=opacity,zorder=1)
        plt.scatter(rad,avgP,color=meancol,marker=symplot,s=markerArea,label=label if radind==0 else '',zorder=2)
    mypredy = magic(radrange,temp)  
    contpredy = cont(radrange,temp)
    maxP.append(max(contpredy))
    plt.plot(radrange,contpredy,color=linecol,linestyle='--',zorder=3)
    plt.plot(radrange,mypredy,color=linecol,linestyle='-',zorder=4)
plt.ylabel(ylab)
plt.xlim([min(radii)-0.1,max(radii)+0.1])
vmaxP = max(maxP)
ymax = round(vmaxP+5,-1)
plt.yticks(np.arange(0,ymax+1,tickstep))
plt.ylim([0,ymax])
plt.xlabel(r'Bubble Radius $(a_0)$')
if args.dualcol:
    plt.legend(loc='upper right',ncol=2,borderpad=0.2,labelspacing=0.4,handletextpad=0.1,borderaxespad=0.2,columnspacing=0.1)
else:
    plt.legend(loc='upper right',borderpad=0.2,labelspacing=0.4,handletextpad=0.1,borderaxespad=0.2,columnspacing=0.1)
if args.image:
    plt.savefig(args.image[0], bbox_inches = 'tight')
if args.display:
    plt.show()
