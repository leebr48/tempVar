import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import argparse 
import pickle as pkl
from itertools import cycle

datacol = 'grey'
meancol = 'k'
linecol = 'k'
opacity = 1
Ptickstep = 5
rhotickstep = 1
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
parser.add_argument('-y','--yaxis', type = str, nargs = '*', help = 'Optional. Specify "dens" to create a density plot and "press" to create a pressure plot. You can also type them in order to set which plot goes on top of the other. If you do not use this option, you will get a combo plot with pressure on top and density on bottom.', required = False)
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
 
if args.yaxis:
    assert len(args.yaxis)<=2, 'You have specified too many plots with -y.'

def magic(yax,lig,rad,temp,surforient):
    if yax == 'maxPress(bar)': 
        A_P = 99.33750903
        B_P = 0.664284758
        C_P = 0.000240928
        if surforient == '001':
            D_P = 117.294327
            E_P = 0.845432613
            F_P = 0.000446182
            G_P = 0.149804059
            H_P = 0.599528882
            I_P = 0.276902518
        if surforient == '011':
            D_P = 130.7352995
            E_P = 0.648058316
            F_P = 0.000508496
            G_P = 0.644859196
            H_P = 0.65904144
            I_P = 0.061593012
        Pb = A_P*rad**(-B_P)*np.exp(-C_P*temp)
        P0 = D_P*np.exp(-E_P*rad-F_P*temp)
        a_P = G_P*rad**(H_P)*temp**(I_P)
        return Pb+(P0-Pb)*np.exp((l0-lig)/a_P)
    if yax == 'dens(nHe/nv)':
        A_rho = 6.363621915
        B_rho = 0.379481286
        C_rho = 0.000166676
        if surforient == '001':
            D_rho = 2.682163389
            E_rho = 0.127923571
            F_rho = 0.000532157
            G_rho = 1.004882167
            H_rho = 0.322211731
            I_rho = 0.000208461
        if surforient == '011':
            D_rho = 3.56716549
            E_rho = 0.1249488
            F_rho = 0.000519787
            G_rho = 1.418142836
            H_rho = 0.324761283
            I_rho = -1.58791E-05
        rhob = A_rho*rad**(-B_rho)*np.exp(-C_rho*temp)
        rho0 = D_rho*np.exp(-E_rho*rad-F_rho*temp)
        a_rho = G_rho*rad**(H_rho)*np.exp(I_rho*temp)
        return rhob+(rho0-rhob)*np.exp((l0-lig)/a_rho)

def boring(x):
    T = tuple(x)
    if len(T) > 1:
        return T
    else:
        return T[0]

df = pd.read_pickle(args.filename[0])

a = pd.DataFrame(columns = df.columns) #Create a temporary blank dataframe. 
for index,series in df.iterrows():
    truth = {'orient':False,'radius':False,'temp':False,'lig':False}
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
choice2 = [r'Relief Pressure (GPa)',r'Relief Density $\mathrm{\left(n_{He}/n_{V}\right)}$']

choices = []
if not args.yaxis:
    choices = [0,1]
else:
    for item in args.yaxis:
        if item == 'press':
            choices.append(0)
        if item == 'dens':
            choices.append(1)
    
axind = 0
if len(choices) == 2:
    fig, axs = plt.subplots(2)
    plt.subplots_adjust(hspace=0.2)
else:
    fig, axs = plt.subplots(1)
    axs = [axs]

for pick in choices: 
    aggIn = choice1[pick]
    ylab = choice2[pick]
    statdf = dfGrouped[aggIn].agg(['mean','count','std']) #Can change aggIn to custom (eg calculated) quantity if desired. 
    rawdf = dfMoreGrouped[aggIn].agg('mean') #The mean is just a dummy function because of how this is grouped - sum should work as well. 
    raddf = dfMoreGrouped['rad'].agg('mean') #The mean is just a dummy function

    meanplot = statdf['mean'].unstack().values
    rawplot = rawdf.unstack().unstack().values
    radplot = raddf.unstack().unstack().values
            
    lab = statdf['mean'].unstack().transpose().columns.to_list()
    ind = statdf['mean'].unstack().transpose().index.to_list()
    indplot = np.array([ind,]*len(lab)) 

    rawlab = rawdf.unstack().unstack().transpose().columns.to_list()
    rawind = [i[1] for i in rawdf.unstack().unstack().transpose().index.to_list()]
    rawindplot = np.array([rawind,]*len(rawlab))

    symlist = cycle(['d','s','^','o','v','*','2','X'])

    colors = cm.rainbow(np.linspace(0,1,len(meanplot)))
    
    ################################################ Plotting starts here. 

    ydatmaxes = []
    for i in range(len(meanplot)):
        symplot = next(symlist)
        if args.color:   
            col = colors[i]
            linecol = col
            datacol = col
            meancol = col
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
                    label += r'$\mathrm{R=}$'+str(int(lab[i][2]))+r' $a_0$'
                if tag == 't':
                    label += r'$\mathrm{T=}$'+str(int(lab[i][3]))+' K'
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
        intermed = []
        for item in x:
            if item not in intermed:
                intermed = np.append(intermed,item)
        xclean = np.array(sorted(intermed))
        yavs = []
        for item in xclean:
            indices = np.where(x==item)
            matchingys = []
            for index in indices:
                matchingys.append(y[index])
            yavs.append(np.mean(matchingys)) 
        xcleaner = np.linspace(min(xclean),max(xclean),1000)
        predy = magic(aggIn,xcleaner,R,float(lab[i][3]),str(lab[i][1]))  
        axs[axind].scatter(x,y,color=datacol,facecolors='none',marker=symplot,s=markerArea,alpha=opacity)
        axs[axind].scatter(xclean,yavs,color=meancol,marker=symplot,label=label,s=markerArea)
        axs[axind].plot(xcleaner,predy,color=linecol,linestyle='-')
        axs[axind].set_ylabel(ylab)
        axs[axind].xaxis.set_major_locator(MaxNLocator(integer=True,nbins=int(np.ceil(max(indplot[i])))))
        axs[axind].set_xlim([0,np.ceil(max(indplot[i]))+0.1])
        ydatmaxes.append(max(y))
    ydatmax = max(ydatmaxes)
    if ydatmax > 10:
        ymax = round(ydatmax+5,-1)
        tickstep = Ptickstep
    else:
        ymax = round(ydatmax+0.5,0)
        tickstep = rhotickstep
    axs[axind].set_yticks(np.arange(0,ymax+1,tickstep))
    axs[axind].set_ylim([0,ymax])
    axind +=1

fig.set_figheight(figheight)
fig.set_figwidth(figwidth)
axs[-1].set_xlabel(r'Ligament Thickness $(a_0)$')
if args.dualcol:
    axs[-1].legend(loc='lower right',ncol=2)
else:
    axs[-1].legend(loc='lower right')
plt.subplots_adjust(hspace=0.2)
if args.image:
    plt.savefig(args.image[0], bbox_inches = 'tight')
if args.display:
    plt.show()
