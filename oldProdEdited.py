#! /usr/bin/python
# This version of the script contains edits from Dr. Hammond. 
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import argparse 
from scipy import stats
from scipy.optimize import curve_fit
import pickle as pkl
from itertools import cycle
import math

ao = {'500':3.17,'1000':3.183,'1200':3.189,'1500':3.199,'2000':3.214,'2500':3.229} #'temperature(K)':latPar(Angs)

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ('Times, Times New Roman, Nimbus Roman No9 L, STIX')
#plt.rcParams['mathtext.fontset'] = 'custom'
#plt.rcParams['mathtext.rm'] = 'Times'
#plt.rcParams['mathtext.it'] = 'Times:italic'
#plt.rcParams['mathtext.bf'] = 'Times:bold'
#plt.rcParams['font.cursive'] = 'Times'
plt.rcParams['axes.unicode_minus']=False
plt.rcParams['lines.linewidth'] = 0.8
plt.rcParams.update({'font.size':10}) 

plt.rc('text',usetex=True)

markerArea = 8 #Squared pixels

#-------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-f','--filename', type = str, nargs = 1, help = 'Input dataframe pickle to read from.', required = True)
parser.add_argument('-s','--shape', type = str, nargs = '*', help = 'Optional filter operation. Include only datapoints for bubbles whose shapes matches the given shape(s) in the analysis. For spheres, type "sph".', required = False)
parser.add_argument('-o','--orient', type = str, nargs = '*', help = 'Optional filter operation. Include only datapoints for surfaces whose surface orientations match the given orientation(s) in the analysis. If you have one or more zeros in the surface orientation, they should lead (e.g. type "011" intead of "110").', required = False)
parser.add_argument('-r','--radius', type = float, nargs = '*', help = 'Optional filter operation. Include only datapoints for bubbles whose radii match the given radius/radii in the analysis. Give the radius in units of lattice parameters (e.g. for a r = 3 ao bubble, type "3").', required = False)
parser.add_argument('-t','--temp', type = float, nargs = '*', help = 'Optional filter operation. Include only datapoints for systems whose temperatures match the given temperature(s) in the analysis. Give the temperature in Kelvin, but do not type "K".', required = False)
parser.add_argument('-l','--lig', type = float, nargs = '*', help = 'Optional filter operation. Include only datapoints for bubbles whose ligaments match the given ligament(s) in the analysis. Give the ligament thickness in units of lattice parameters (e.g. for a l = 4 ao bubble, type "4").', required = False)
parser.add_argument('-d','--display', action = 'store_true', default = False, help = 'Display the produced graph on the screen.', required = False)
parser.add_argument('-i','--image', type = str, nargs = '*', help = 'Save the produced graph with the given name. INCLUDE a file extension (png or pdf). If you used the -c flag, you must supply three names. Otherwise, supply one name.', required = False) 
parser.add_argument('-z','--zlabel', type = str, nargs = '*', help = "Optional. Specify what labels to include for each line in the legend, and the order in which they appear. 's' -> shape, 'o' -> surface orientation, 'r' -> bubble radius, 't' -> temperature", required = False) 
parser.add_argument('--dualcol', action = 'store_true', default = False, help = 'Optional. Toggle on to split the figure legend into two columns.', required = False)
parser.add_argument('--figheight', type = float, nargs = 1, default = 10, help = 'Optional. Specify the figure height in inches. Defaults to 10. Note that exact height may vary due to screen DPI.', required = False)
parser.add_argument('--figwidth', type = float, nargs = 1, default = 18, help = 'Optional. Specify the figure width in inches. Defaults to 18. Note that exact width may vary due to screen DPI.', required = False)
parser.add_argument('-p','--printparam', action = 'store_true', default = False, help = 'Display the fitted model parameters.', required = False)
args = parser.parse_args()
#------------------------------------------------

al = 0.5*(1-95/100) #This script doesn't actually make error bars anymore, I just don't want to retype all the statistical parts of the code. 

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

#This is the fitting function.
def f(l, K, Y0, v):
    Q = ((K/Y0)**v-1).astype(np.float128)
    return (K/(1+Q*np.exp(-(l-l0)))**(1/v)).astype(np.float128)
    
df = pd.read_pickle(args.filename[0])

a = pd.DataFrame(columns = df.columns) #Create a temporary blank dataframe. 
for index,series in df.iterrows():
    truth = {'shape':False,'orient':False,'radius':False,'temp':False,'lig':False}
    if args.shape:
        if series['bubShape'] in args.shape:
           truth['shape'] = True 
    else: 
        truth['shape'] = True
    if args.orient:
        if series['surfOrient'] in args.orient:
            truth['orient'] = True
    else:
        truth['orient'] = True
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
    if args.lig:
        if series['lig(ao)'] in args.lig:
            truth['lig'] = True
    else:
        truth['lig'] = True
    if all(value == True for value in truth.values()):
        a = a.append(series)
    truth = None

a['reliefType'] = a['reliefType'].map({'ERROR':'ERROR', 'loopPunch':'loopPunch', 'bottomPunch':'bottomPunch', 'burst':'burst', 'timeOut':'bottomPunch'})

df = a
a = None

df['dens(nHe/nv)'] = pd.to_numeric(df['maxPHe']/df['nv'])
df['rad'] = df['Rs(ao)'] 

grouping = ['bubShape','surfOrient','Rs(ao)','temp(K)','lig(ao)']  
dfGrouped = df.groupby(by=grouping)

moreGrouping = ['bubShape','surfOrient','Rs(ao)','temp(K)','lig(ao)','runNum'] 
dfMoreGrouped = df.groupby(by=moreGrouping)

choice1 = ['maxPress(bar)','dens(nHe/nv)']
choice2 = [r'Relief Pressure (GPa)',r'Relief Density ($n_{\scriptsize\textrm{He}}$/$n_{\scriptsize\textrm{V}}$)']

for pick in [0,1]: 
    aggIn = choice1[pick]
    ylab = choice2[pick]
    statdf = dfGrouped[aggIn].agg(['mean','count','std']) #Can change aggIn to custom (eg calculated) quantity if desired. 
    rawdf = dfMoreGrouped[aggIn].agg('mean') #The mean is just a dummy function because of how this is grouped - sum should work as well. 
    raddf = dfMoreGrouped['rad'].agg('mean') #The mean is just a dummy function

    def boring(x):
        T = tuple(x)
        if len(T) > 1:
            return T
        else:
            return T[0]

    burstdf = dfMoreGrouped['reliefType'].agg(boring)

    ciHi = []
    ciLow = []
    for i in statdf.index:
        m,n,s = statdf.loc[i]
        t = stats.t.ppf(1-al,n-1) #T distribution output. 
        ciHi.append(t*s/np.sqrt(n))
        ciLow.append(t*s/np.sqrt(n)) #NOTE: The error is symmetric when calculated in this way. You could theoretically just have one of these columns, but two is more general and (may???) come in usefule later.

    statdf['ciLow'] = ciLow
    statdf['ciHi'] = ciHi
    statdf['reliefType'] = dfGrouped['reliefType'].agg(boring)

    outVec = []
    annotVec = []
    for i in statdf['reliefType']:
        a = i
        marked = False
        if 'ERROR' in a or 'timeOut' in a:
            out = 'ERROR'
            annotfrac = ''
        else:
            if all(x == a[0] for x in a):
                out = a[0]
                annotfrac = ''
                marked = True
            else:
                if not marked:
                    out = 'Mix'
                    bursts = a.count('burst')
                    tops = a.count('loopPunch')
                    bottoms = a.count('bottomPunch')
                    annotfrac = str(bursts)+'Bst/'+str(tops)+'TP/'+str(bottoms)+'BP'
        outVec.append(out)
        annotVec.append(annotfrac)
    statdf['reliefMethod'] = outVec
    statdf['annotation'] = annotVec

    plotmark = []
    for i in outVec:
        if i == 'ERROR':
            mark = 'x'#Mark errors with an X on the plot. 
        if i == 'loopPunch':
            mark = '^'#Mark loop punching with an upwards triangle on the plot. 
        if i == 'bottomPunch':
            mark = 'v' #Mark bottom punching with an downwards triangle on the plot. 
        if i == 'burst': 
            mark = 'o' #Mark bursting with a circle on the plot.
        if i == 'Mix':
            mark = '*' #Mark mixed pressure relief behavior with a star on the plot. 
        plotmark.append(mark)
        mark = None
    statdf['figMark'] = plotmark

    plotmark2 = []
    for i in burstdf:
        if i == 'ERROR':
            mark = 'x'#Mark errors with an X on the plot. 
        if i == 'loopPunch':
            mark = '^'#Mark loop punching with an updwards triangle on the plot. 
        if i == 'bottomPunch':
            mark = 'v'#Mark loop punching with a downwards triangle on the plot. 
        if i == 'burst': 
            mark = 'o' #Mark bursting with a circle on the plot.
        plotmark2.append(mark)
        mark = None
    markdf = pd.Series(plotmark2,burstdf.index)

    meanplot = statdf['mean'].unstack().values
    ciLowplot = statdf['ciLow'].unstack().values
    ciHiplot = statdf['ciHi'].unstack().values
    markplot = statdf['figMark'].unstack().values
    annotplot = statdf['annotation'].unstack().values
    rawplot = rawdf.unstack().unstack().values
    symplot = markdf.unstack().unstack().values
    radplot = raddf.unstack().unstack().values

    new = []
    for i in markplot:
        a = ['None' if str(x)=='nan' else x for x in i]
        new.append(a)
    markplot = new
    new = None

    new = []
    for i in symplot:
        inner = []
        for j in i:
            if str(j) != 'nan':
                inner.append(j)
        new.append(inner)
    symplot = new
    new = None
            
    lab = statdf['mean'].unstack().transpose().columns.to_list()
    ind = statdf['mean'].unstack().transpose().index.to_list()
    indplot = np.array([ind,]*len(lab)) 

    rawlab = rawdf.unstack().unstack().transpose().columns.to_list()
    rawind = [i[1] for i in rawdf.unstack().unstack().transpose().index.to_list()]
    rawindplot = np.array([rawind,]*len(rawlab))

    symlist = cycle(['d','s','^','o','v','*','2','X'])

    colors = cm.rainbow(np.linspace(0,1,len(meanplot)))
    if aggIn == 'maxPress(bar)':
        fig, axs = plt.subplots(2)
        axind = 0
    if aggIn == 'dens(nHe/nv)':
        axind = 1
    ################################################ Plotting starts here. 
    col = 'k'
    pointcol = 'k'
    for i in range(len(meanplot)):
        symplot = next(symlist)
        if not args.zlabel:
            label = str(lab[i][0])+','+grouping[1]+'='+str(lab[i][1])+','+grouping[2]+'='+str(lab[i][2])+','+grouping[3]+'='+str(lab[i][3])
        else:
            label = ''
            for tag in args.zlabel:
                if tag == 's':
                    if lab[i][0] == 'sph':
                        label += 'Sphere'
                    elif lab[i][0] == 'ell':
                        label += 'Ellipse'
                if tag == 'o':
                    broken = list(str(lab[i][1]))
                    inted = [int(i) for i in broken]
                    shuffled = sorted(inted, reverse=True)
                    restringed = [str(i) for i in shuffled]
                    rejoined = ''.join(restringed)
                    label += 'W('+rejoined+')'
                if tag == 'r':
                    label += r'$R={}$'+str(int(lab[i][2]))+r' $a_0^{}$'
                if tag == 't':
                    label += r'$T={}$'+str(int(lab[i][3]))+' K'
                if tag != args.zlabel[-1]:
                    label += ', '
        stacked = np.vstack((rawindplot[i],rawplot[i],radplot[i])).T
        cleaned = []
        for row in stacked:
            flag = True
            for item in row:
                if str(item) == 'nan':
                    flag = False
            if flag == True:
                if np.size(cleaned) == 0:
                    cleaned = row
                else: 
                    cleaned = np.vstack((cleaned,row))
        x = []
        y = []
        r = []
        if np.ndim(cleaned) > 1:
            for item in cleaned:
                x = np.append(x,item[0])
                y = np.append(y,item[1])
                r = np.append(r,item[2])
        else:
            x = cleaned[0]
            y = cleaned[1]
            r = cleaned[2]
        if aggIn == 'maxPress(bar)':
            y = y/10000 #converts bar to GPa
        assert all(x==r[0] for x in r), 'Radius vector has multiple values - something is very wrong.'
        R = r[0]
        l0 = min(x)
        #popt, pcov = curve_fit(f,x,y,p0=[y[-1],y[0],1],bounds=(0,np.inf))          
        popt, pcov = curve_fit(f,x,y,p0=[y[-1],y[0],1],bounds=(0,np.inf))         
        if args.printparam:
            print('The model parameters for '+label+' are: ' + str(popt))
        intermed = []
        for item in x:
            if item not in intermed:
                intermed = np.append(intermed,item)
        xclean = np.array(sorted(intermed))
        xclean = np.linspace(min(xclean),max(xclean),1000)
        predy = f(xclean,*popt)
        axs[axind].scatter(x,y,color=pointcol,marker=symplot,label=label,s=markerArea)
        axs[axind].plot(xclean,predy,color=col,linestyle='-')
        axs[axind].set_ylabel(ylab)
        axs[axind].xaxis.set_major_locator(MaxNLocator(integer=True,nbins=int(np.ceil(max(indplot[i])))))
        axs[axind].set_xlim([0,np.ceil(max(indplot[i]))+0.1])
        # FIXME
        ybound = axs[axind].get_ybound()
        ylim = 5*math.ceil(ybound[1]/5)
        axs[axind].set_ylim([0,ylim])
        # END FIXME
fig.set_figheight(figheight)
fig.set_figwidth(figwidth)
axs[1].set_xlabel(r'Ligament Thickness $(a_0^{})$')
if args.dualcol:
    axs[1].legend(loc='lower right',ncol=2)
else:
    axs[1].legend(loc='lower right')
plt.subplots_adjust(hspace=0.2)
if args.image:
    plt.savefig(args.image[0], bbox_inches = 'tight')
if args.display:
    plt.show()
