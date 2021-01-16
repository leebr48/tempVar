import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import argparse 
from scipy import stats
from scipy.optimize import curve_fit
import pickle as pkl

ao = {'500':3.17,'1000':3.183,'1200':3.189,'1500':3.199,'2000':3.214,'2500':3.229} #'temperature(K)':latPar(Angs)
#-------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-f','--filename', type = str, nargs = 1, help = 'Input dataframe pickle to read from.', required = True)
parser.add_argument('-y','--yaxis', type = str, nargs = 1, help = "Choose the variable plotted on the y axis. 'press' -> bubble pressure, 'dens' -> bubble density.", required = True)
parser.add_argument('-s','--shape', type = str, nargs = '*', help = 'Optional filter operation. Include only datapoints for bubbles whose shapes matches the given shape(s) in the analysis. For spheres, type "sph".', required = False)
parser.add_argument('-o','--orient', type = str, nargs = '*', help = 'Optional filter operation. Include only datapoints for surfaces whose surface orientations match the given orientation(s) in the analysis. If you have one or more zeros in the surface orientation, they should lead (e.g. type "011" intead of "110").', required = False)
parser.add_argument('-r','--radius', type = float, nargs = '*', help = 'Optional filter operation. Include only datapoints for bubbles whose radii match the given radius/radii in the analysis. Give the radius in units of lattice parameters (e.g. for a r = 3 ao bubble, type "3").', required = False)
parser.add_argument('-t','--temp', type = float, nargs = '*', help = 'Optional filter operation. Include only datapoints for systems whose temperatures match the given temperature(s) in the analysis. Give the temperature in Kelvin, but do not type "K".', required = False)
parser.add_argument('-l','--lig', type = float, nargs = '*', help = 'Optional filter operation. Include only datapoints for bubbles whose ligaments match the given ligament(s) in the analysis. Give the ligament thickness in units of lattice parameters (e.g. for a l = 4 ao bubble, type "4").', required = False)
parser.add_argument('-m','--model', action = 'store_true', default = False, help = 'Enable fitting curves to simulation data.', required = False)
parser.add_argument('-d','--display', action = 'store_true', default = False, help = 'Display the produced graph on the screen.', required = False)
parser.add_argument('-i','--image', type = str, nargs = '*', help = 'Save the produced graph with the given name. INCLUDE a file extension (png or pdf). If you used the -c flag, you must supply three names. Otherwise, supply one name.', required = False) 
parser.add_argument('-p','--printparam', action = 'store_true', default = False, help = 'Display the fitted model parameters.', required = False)
parser.add_argument('-c','--chunk', action = 'store_true', default = False, help = 'Divide output into separate graphs rather than one subplot.', required = False) 
parser.add_argument('-a','--annotate', action = 'store_true', default = False, help = 'Display burst/topPunch/bottomPunch annotations on the means plot.', required = False)
parser.add_argument('-b','--bound', action = 'store_true', default = False, help = 'Display horizontal lines where continuum mechanics roughly predicts the bulk loop punching pressure should be.', required = False)
parser.add_argument('-z','--zlabel', type = str, nargs = '*', help = "Optional. Specify what labels to include for each line in the legend, and the order in which they appear. 's' -> shape, 'o' -> surface orientation, 'r' -> bubble radius, 't' -> temperature", required = False) 
parser.add_argument('-g','--goodness', action = 'store_true', default = False, help = 'Optional. Toggle to display R^2 values.', required = False)
parser.add_argument('-x','--extra', action = 'store_true', default = False, help = 'Optional. Toggle to display model output on the mean/error-bar graph. NOTE that the -m option MUST be used for this option to work.', required = False)
parser.add_argument('-p2','--printparamfig', action = 'store_true', default = False, help = 'Optional. Print the solver parameters on the plot itself.', required = False)
parser.add_argument('-w','--wow', action = 'store_true', default = False, help = 'Optional. Toggle on to use the newer model that implicitely specifies Pb rather than explicitly.', required = False) 
parser.add_argument('--CI', type = float, nargs = 1, default = 95, help = 'Optional. Specify the confidence interval for the graphs as an implicit percent (e.g. for a 90%% confidence interval, type 90). Defaults to 95.', required = False)
parser.add_argument('--figheight', type = float, nargs = 1, default = 10, help = 'Optional. Specify the figure height in inches. Defaults to 10. Note that exact height may vary due to screen DPI.', required = False)
parser.add_argument('--figwidth', type = float, nargs = 1, default = 18, help = 'Optional. Specify the figure width in inches. Defaults to 18. Note that exact width may vary due to screen DPI.', required = False)
args = parser.parse_args()
#------------------------------------------------

if isinstance(args.CI,int):
    CI = args.CI
else:
    CI = args.CI[0]
if isinstance(args.figheight,int):
    figheight = args.figheight
else:
    figheight = args.figheight[0]
if isinstance(args.figwidth,int):
    figwidth = args.figwidth
else:
    figwidth = args.figwidth[0]

al = 0.5*(1-CI/100) #SciPy defines alpha a bit differently than your statistics book (1-tail vs 2-tail). This form should work fine for SciPy. 

if args.image:
    if not args.chunk:
        assert len(args.image)==1, 'Single images must be given a single name.'
    if args.chunk and args.model:
        assert len(args.image)==3, 'When three separate graphs are created, they must be given three separate names.'
    if args.chunk and not args.model: 
        assert len(args.image)==1, 'Single images must be given a single name.'

#This is the fitting function.
if any((args.model,args.printparam,args.extra)) and args.wow:
    def f(l, A, B): 
        return A*np.exp(-B*R)*(1-np.exp(-B*l)) 
elif any((args.model,args.printparam,args.extra)):
    def f(x, a, b):
        return a*(1-np.exp(-b*x))  
    
if args.bound:
    #This function finds the shear modulus based on the MTS model. 
    #See https://doi.org/10.1103%2FPhysRevB.2.3952
    #Some modifications are built in to make the results more in line with MS simulations - kinda janky, but the results should be 'close enough'. 
    C440 = 16.31577e11 #dyn/cm^2
    s44 = 0.568995e11 #dyn/cm^2
    t44 = 346 #K
    fudge = 1.03731889282655 #This (janky) empirical factor makes the MTS model line up with MS predictions.
    def G(T):
        Cij = C440-s44/(np.exp(t44/T)-1) #dyn/cm^2
        #return Cij/10*fudge #Pa
        return Cij/10 #Pa
    for key in ao:
        ao[key] *= 10**(-10)
    #This is just the burgers vector at different temperatures.
    def b(T):
        return ao[str(T)]*np.sqrt(3)/2
    #This bit is from Wolfer. Input bubble radius in ao and temperature in Kelvin.  
    gamma = 2.67 #J/m^2. Approximated as constant - should make less of an impact than if G were constant. 
    def Pcont(r,T):
        latpar = ao[str(int(T))]
        rad = r*latpar #m
        shear = G(np.float(T))
        burg = b(int(T))
        Pbulk = (2*gamma+shear*burg)/rad #Pa
        return Pbulk/100000 #bar

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

if args.yaxis[0] == 'press':
    aggIn = 'maxPress(bar)'
    ylab = 'Relief Pressure (bar)'
if args.yaxis[0] == 'dens':
    aggIn = 'dens(nHe/nv)'
    ylab = r'Bubble Relief Density $\mathrm{(n_{He}/n_{V})}$'

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

if args.display or args.image or args.model or args.printparam:
    colors = cm.rainbow(np.linspace(0,1,len(meanplot)))
    if not args.chunk:
        if args.model: 
            fig1, (ax1,ax2,ax3) = plt.subplots(3,1)
        else:
            fig1,ax1 = plt.subplots()
    else:
        if args.model:
            fig1, ax1 = plt.subplots()
            fig2, ax2 = plt.subplots()
            fig3, ax3 = plt.subplots()
        else: 
            fig1, ax1 = plt.subplots()

    for i in range(len(meanplot)):
        col = colors[i] 
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
                    label += '('+rejoined+')'
                if tag == 'r':
                    label += r'$\mathrm{R=}$'+str(int(lab[i][2]))+r' $\mathrm{a_o}$'
                if tag == 't':
                    label += r'$\mathrm{T=}$'+str(int(lab[i][3]))+' K'
                if tag != args.zlabel[-1]:
                    label += ', '
        if any((args.model,args.printparam,args.extra)):
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
            assert all(x==r[0] for x in r), 'Radius vector has multiple values - something is very wrong.'
            R = r[0]
            popt, pcov = curve_fit(f,x,y,p0=[y[-1],1],bounds=([0,0],[np.inf,np.inf]))   
            if args.printparam and args.wow:
                print('The model parameters (A,B) for '+label+' are: ' + str(popt))
            elif args.printparam:
                print('The model parameters (Pr,C) for '+label+' are: ' + str(popt))
            if args.model or args.extra:
                residuals = y - f(x,*popt)
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((y-np.mean(y))**2)
                r_squared = 1 - (ss_res/ss_tot) 
                std_res = residuals/np.std(residuals) 
                xclean = []
                for item in x:
                    if item not in xclean:
                        xclean = np.append(xclean,item)
                predy = f(xclean,*popt)
                predydouble = f(x,*popt)
                if args.model:
                    for xin,yin,symin in zip(x,y,symplot[i]):
                        ax2.scatter(xin,yin,color=col,marker=symin)
                    if not args.chunk:
                        label2 = r'$\mathrm{R^2}=$'+str(round(r_squared,4))
                    else:
                        if args.printparamfig and args.wow:
                            newlab = label + r', $\mathrm{A=}$'+str(round(popt[0]))+' bar, '+r'$\mathrm{B=}$'+str(round(popt[1],4))+r' $\mathrm{a_{o}^{-1}}$'
                        elif args.printparamfig:
                            newlab = label + r', $\mathrm{P_{b}=}$'+str(round(popt[0]))+' bar, '+r'$\mathrm{C=}$'+str(round(popt[1],4))+r' $\mathrm{a_{o}^{-1}}$'
                        else:
                            newlab = label
                        if args.goodness:
                            label2 = newlab +r', $\mathrm{R^2}=$'+str(round(r_squared,4))
                        else:
                            label2 = newlab
                    ax2.plot(xclean,predy,color=col,linestyle='--',label=label2)
                    for xin,yin,symin in zip(predydouble,std_res,symplot[i]):
                        ax3.scatter(xin,yin,color=col,marker=symin)
                if args.extra:
                    ax1.plot(xclean,predy,color=col,linestyle='--') 
        if all((args.extra,args.printparamfig,args.wow)):
            label = label + r', $\mathrm{A=}$'+str(round(popt[0]))+' bar, '+r'$\mathrm{B=}$'+str(round(popt[1],4))+r' $\mathrm{a_{o}^{-1}}$'
        elif args.extra and args.printparamfig:
            label = label + r', $\mathrm{P_{b}=}$'+str(round(popt[0]))+' bar, '+r'$\mathrm{C=}$'+str(round(popt[1],4))+r' $\mathrm{a_{o}^{-1}}$'
        if args.extra and args.goodness:
            label = label + r', $\mathrm{R^2}=$'+str(round(r_squared,4))
        ax1.errorbar(indplot[i],meanplot[i],yerr=[ciLowplot[i],ciHiplot[i]],color=col,label=label)
        if args.bound:
            ax1.axhline(y=Pcont(lab[i][2],lab[i][3]), color=col, linestyle='--')
        ax1.xaxis.set_major_locator(MaxNLocator(integer=True,nbins=int(np.ceil(max(indplot[i])))))
        ax1.set_xlim([0,np.ceil(max(indplot[i]))+0.1])
        if args.annotate:
            for j in range(len(meanplot[i])):
                ax1.annotate(annotplot[i][j],(indplot[i][j],meanplot[i][j]),xytext=(indplot[i][j]+0.03,meanplot[i][j]-7500))
        for x,y,m in zip(indplot[i],meanplot[i],markplot[i]):
            ax1.scatter(x,y,marker=m,color=col) 
    ax1.set_xlabel(r'Ligament Thickness $\mathrm{(a_o)}$')
    ax1.set_ylabel('Mean '+ylab)
    if args.chunk or not args.model:
        ax1.legend(loc='lower right')
    else:
        ax1.legend(loc='best')
    fig1.set_figheight(figheight)
    fig1.set_figwidth(figwidth)
    if args.image and args.chunk:
        fig1.savefig(args.image[0], bbox_inches = 'tight')
    if args.chunk and args.model:
        fig2.set_figheight(figheight)
        fig2.set_figwidth(figwidth)
        fig3.set_figheight(figheight)
        fig3.set_figwidth(figwidth)
    if args.model:
        ax2.xaxis.set_major_locator(MaxNLocator(integer=True,nbins=int(np.ceil(max(indplot[i])))))
        ax2.set_xlim([0,np.ceil(max(indplot[i]))+0.1])
        ax2.set_xlabel(r'Ligament Thickness $\mathrm{(a_o)}$')
        ax2.set_ylabel('Predicted and Actual '+ylab)
        if args.chunk or not args.model:
            ax2.legend(loc='lower right')
        else:
            ax2.legend(loc='best')
        if args.chunk and args.image:
            fig2.savefig(args.image[1], bbox_inches = 'tight')
        ax3.set_xlabel('Predicted '+ylab)
        ax3.set_ylabel('Standardized Residuals')
        ax3.axhline(y=0, color='k', linestyle='-')
        if args.chunk and args.image:
            fig3.savefig(args.image[2], bbox_inches = 'tight')
        if not args.chunk:
            plt.subplots_adjust(hspace=0.3)
    if args.image and not args.chunk:
        plt.savefig(args.image[0], bbox_inches = 'tight')
    if args.display:
        plt.show()
