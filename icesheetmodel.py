import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import tkinter as tk
import math
from tkinter import *

# Universal constants ===================================================================================================================

dt = 1                                                              # Timestep (years)
tStepsNum = 40000                                                   # Number of timesteps
R = 6371000                                                         # Radius of the Earth (meters)
m = 3                                                               # Constant Nye's Law
A = 2                                                               # Diffusivity constant (m^(-3/2)/yr)
dICE = 830; dBED = dICE*3                                           # Density of glacier ice and bedrock (kg/m^3)
dRatio = dICE/dBED                                                  # Ratio of glacier ice and bedrock density
tscale = 1000                                                       # Time scale (years) --> 1000
alpha = -2e-6; beta = 1.2e-3                                        # ELA slope parameters

# Setting grid ===========================================================================================================================

numCells_x = 25                                                     # x grid cells for latitude
xEdge = np.linspace(0,5000*1000,numCells_x+1)                       # x grid cell edges (m)
                                                                    # Continental boundaries: 
                                                                    # 0 km = polar ocean shore --> 5000 km = gulf coast shore
x = (xEdge[0:numCells_x]+xEdge[1:numCells_x+1])/2                   # Grid centers (meters)
dx = xEdge[1] - xEdge[0]                                            # Incremental change in x (meters)

# Initial Conditions ======================================================================================================================

H = np.zeros_like(x)                                                # Thickness output (meters)
HUD = np.zeros_like(x)                                              # Thickness above undisturbed bedrock (meters)
h = np.zeros_like(x)                                                # Bedrock depression (meters)
M = np.zeros_like(xEdge)                                            # Vertically integrated mass flux (m^2/yr) ...I think
D = np.zeros_like(xEdge)                                            # Diffusion (m^2/yr)
conv1 = np.zeros_like(x)                                            # First part of the convergence calculation (if needed)
conv2 = np.zeros_like(x)                                            # Second part of the convergence calculation
B = np.zeros_like(x)                                                # Mass balance (m/yr)
ELA = np.zeros_like(x)                                              # Equilibrium Line Altitude (m)


# Functions
def diffCalc(H2,H1,HUD2,HUD1,A,m,dx):
    # Diffusion calculation
    #   HUD2 = HUD[j], HUD1 = HUD[j-1]
    #   H2 = H[j], H1 = H[j-1]
    Hint = (H2 + H1)/2                                              # Interpolated thickness (m)
    D = A*Hint**(m+1)*(((HUD2-HUD1)/(dx))**2)**((m-1)/2)            # Diffusion equation
    return D
def fluxCalc(HUD2,HUD1,D,dx):
    # Integrated mass flux calculation
    #   HSL2 = HSL[j], HSL1 = HSL[j-1]
    M = D*((HUD2-HUD1)/dx)                                          # Vertically integrated mass flux
    return M
def balanceCalc(x,HUD,alpha,beta,P):
    # Mass balance
    #   la = latRad[j]
    #   H = H[j]
    B = alpha*(x - P) + beta*HUD                                    # Mass balance (m/yr)
    return B
def convCalc(M2,M1,dx):
    # Convergence calculation
    #   M2 = M[j+1], M1 = M[j]
    #   cosLatE2 = cosLatE[j+1], cosLatE1 = cosLatE[j]
    #   cosLa = cosLat[j]
    conv2 = (M2 - M1)/dx                                            # Convergence
    return conv2
def ELAcalc(x,P):
    # Equilibrium Line Altitude calculation
    #   la = latRad[j]
    ELA = -alpha/beta*(x - P)                                       # ELA
    return ELA

# Run

def hystRun_noGIA(tStepsNum):
    dt = 1
    print('Grid Information ================================================')
    print('  {0} grid cells'.format(numCells_x))
    print('Time Information ================================================')
    print('  Integration Time: {0} years'.format(tStepsNum))
    print('  Time Step = {0} years'.format(dt))
    print('Glen\'s Flow Law ================================================')
    print('  A = ',A)
    print('  m = ',m)
    print('GIA Information =================================================')
    print('  No bedrock adjustment')
    print('ELA Paramaters ==================================================')
    print('  alpha = ',alpha)
    print('  beta = ',beta)
    print('  ELA Slope (-alpha/beta) = {0:.5}'.format(-alpha/beta))
    print('  Climate Point = varying')
    print('')
    print('Changing climate...')
    print('Crushing snow into glacier ice...')
    P_points = [-500000*np.cos((2*np.pi*t)/20000) for t in range(tStepsNum)]

    H_tseries = []
    glaciation = []
    deglaciation = []

    for (n,p) in zip(range(tStepsNum),P_points):
        #print(p)
        for j in range(1,numCells_x):
            
            D[j] = diffCalc(H[j],H[j-1],HUD[j],HUD[j-1],A,m,dx)
        
            M[j] = fluxCalc(HUD[j],HUD[j-1],D[j],dx)
        
        for j in range(numCells_x):

            B[j] = balanceCalc(x[j],HUD[j],alpha,beta,p)

            ELA[j] = ELAcalc(x[j],p)
            if ELA[j] < 0:
                ELA[j] = np.nan

            conv2[j] = convCalc(M[j+1],M[j],dx)

            H[j] = (conv2[j]+B[j])*dt + H[j]
            if H[j] < 0:
                H[j] = 0
            HUD[j] = H[j]
        H_tseries.append(max(H))

        H[numCells_x-1] = 0
        H[0] = 0
        h[0] = 0

        if H_tseries[n] > H_tseries[n-1]:
            glaciation.append(0)
        elif H_tseries[n] < H_tseries[n-1]:
            deglaciation.append(0)
    
    print('')
    print('Complete!')
    print('')
    print('Glaciation length = {0} years'.format(len(glaciation)))
    print('Deglaciation length = {0} years'.format(len(deglaciation)))

    # Plotting

    fig1,axs = plt.subplots(2,sharex=True)
    fig1.set_figheight(8)
    fig1.set_figwidth(8)
    axs[0].plot(range(tStepsNum),P_points)
    axs[0].set_title('Climate Point Position')
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1,1))
    axs[0].yaxis.set_major_formatter(formatter)
    axs[0].set_ylabel('Climate Point\n(m south from\npolar ocean)')
    axs[1].plot(range(tStepsNum),H_tseries)
    axs[1].set_title('Ice Sheet Timeseries')
    axs[1].set_xlabel('Integration Time (yr)')
    axs[1].set_ylabel('Ice Sheet\nThickness (meters)')
    axs[1].text(0,1500,'No bedrock adjustment\nGlaciation: {0} years\nDeglaciation: {1} years'.format(len(glaciation),len(deglaciation)))

    return fig1

def hystRun(tStepsNum,dRatio):
    dt = 1
    print('Grid Information ================================================')
    print('  {0} grid cells'.format(numCells_x))
    print('Time Information ================================================')
    print('  Integration Time: {0} years'.format(tStepsNum))
    print('  Time Step = {0} years'.format(dt))
    print('Glen\'s Flow Law ================================================')
    print('  A = ',A)
    print('  m = ',m)
    print('GIA Information =================================================')
    print('  Ice-to-bedrock density ratio = {0} kg/m^3 / {1} kg/m^3 = {2:.2}'
    .format(dICE,dBED,dRatio))
    print('  Bedrock Response Time = {0} ka'.format(tscale/1000))
    print('ELA Paramaters ==================================================')
    print('  alpha = ',alpha)
    print('  beta = ',beta)
    print('  ELA Slope (-alpha/beta) = {0:.5}'.format(-alpha/beta))
    print('  Climate Point = varying')
    print('')
    print('Changing climate...')
    print('Crushing snow into glacier ice...')
    P_points = [-500000*np.cos((2*np.pi*t)/20000) for t in range(tStepsNum)]

    H_tseries = []
    glaciation = []
    deglaciation = []

    for (n,p) in zip(range(tStepsNum),P_points):
        #print(p)
        for j in range(1,numCells_x):
            
            D[j] = diffCalc(H[j],H[j-1],HUD[j],HUD[j-1],A,m,dx)
        
            M[j] = fluxCalc(HUD[j],HUD[j-1],D[j],dx)
        
        for j in range(numCells_x):

            B[j] = balanceCalc(x[j],HUD[j],alpha,beta,p)

            ELA[j] = ELAcalc(x[j],p)
            if ELA[j] < 0:
                ELA[j] = np.nan

            conv2[j] = convCalc(M[j+1],M[j],dx)

            H[j] = (conv2[j]+B[j])*dt + H[j]
            if H[j] < 0:
                H[j] = 0
            
            HUD[j] = H[j] + h[j]
            h[j] = -((dRatio*H[j] + h[j] - 0)/tscale)*dt + h[j]
            
        H_tseries.append(max(H))

        H[numCells_x-1] = 0
        H[0] = 0
        h[0] = 0

        if H_tseries[n] > H_tseries[n-1]:
            glaciation.append(0)
        elif H_tseries[n] < H_tseries[n-1]:
            deglaciation.append(0)

    print('')
    print('Complete!')
    print('')
    print('Glaciation length = {0} years'.format(len(glaciation)))
    print('Deglaciation length = {0} years'.format(len(deglaciation)))

    # Plotting

    fig1,axs = plt.subplots(2,sharex=True)
    fig1.set_figheight(8)
    fig1.set_figwidth(8)
    axs[0].plot(range(tStepsNum),P_points)
    axs[0].set_title('Climate Point Position')
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1,1))
    axs[0].yaxis.set_major_formatter(formatter)
    axs[0].set_ylabel('Climate Point\n(m south from\npolar ocean)')
    axs[1].plot(range(tStepsNum),H_tseries)
    axs[1].set_title('Ice Sheet Timeseries')
    axs[1].set_xlabel('Integration Time (yr)')
    axs[1].set_ylabel('Ice Sheet Thickness (meters)')
    axs[1].text(0,1500,'With bedrock adjustment\nGlaciation: {0} years\nDeglaciation: {1} years'.format(len(glaciation),len(deglaciation)))

    return fig1

def fixedRun(x,H,ELA,h,HUD,P):
    print('Grid Information ================================================')
    print('  {0} grid cells'.format(numCells_x))
    print('Time Information ================================================')
    print('  Integration Time: {0} years'.format(tStepsNum))
    print('  Time Step = {0} years'.format(dt))
    print('Glen\'s Flow Law ================================================')
    print('  A = ',A)
    print('  m = ',m)
    print('GIA Information =================================================')
    print('  Ice-to-bedrock density ratio = {0} kg/m^3 / {1} kg/m^3 = {2:.2}'
      .format(dICE,dBED,dRatio))
    print('  Bedrock Response Time = {0} ka'.format(tscale/1000))
    print('ELA Paramaters ==================================================')
    print('  alpha = ',alpha)
    print('  beta = ',beta)
    print('  ELA Slope (-alpha/beta) = {0:.5}'.format(-alpha/beta))
    print('  Climate Point = {0} km south from the polar ocean'.format(P/1000))
    print('')
    print('Crushing snow into glacier ice...')
    print('')
    h_tseries = []
    HUD_tseries = []
    for n in range(tStepsNum):
        # Fluxes
        for j in range(1,numCells_x):
            
            D[j] = diffCalc(H[j],H[j-1],HUD[j],HUD[j-1],A,m,dx)
        
            M[j] = fluxCalc(HUD[j],HUD[j-1],D[j],dx)
        
        for j in range(numCells_x):

            B[j] = balanceCalc(x[j],HUD[j],alpha,beta,P)

            ELA[j] = ELAcalc(x[j],P)
            if ELA[j] < 0:
                ELA[j] = np.nan

            conv2[j] = convCalc(M[j+1],M[j],dx)

            H[j] = (conv2[j]+B[j])*dt + H[j]
            if H[j] < 0:
                H[j] = 0
            
            HUD[j] = H[j] + h[j]
            h[j] = -((dRatio*H[j] + h[j] - 0)/tscale)*dt + h[j]
            
        
        h_tseries.append(min(h))
        HUD_tseries.append(max(HUD))

        H[numCells_x-1] = 0
        H[0] = 0
        h[0] = 0
    print('Complete!')

    # Plotting
    fig1 = plt.figure()
    plt.title('North American Ice Sheet\n$x_p$ = {0} km south'.format(P/1000))
    plt.xlabel('Distance (meters south from polar ocean shoreline)')
    plt.ylabel('Altitude (meters)')
    plt.plot(x,ELA,'--',color='black')
    plt.plot(x,h,color='green')
    plt.plot(x,HUD,color='blue')
    plt.legend(['Equilibrium Line Altitude (B = 0)','Bedrock','Ice Surface'])

    fig2 = plt.figure()
    plt.title('Bedrock Depression\n$x_p$ = {0} km south'.format(P/1000))
    plt.xlabel('Time (years)')
    plt.ylabel('Altitude (meters)')
    plt.plot(range(tStepsNum),h_tseries)
    
    fig3 = plt.figure()
    plt.title('Ice Surface Above Undisturbed Bed\n$x_p$ = {0} km south'.format(P/1000))
    plt.xlabel('Time (years)')
    plt.ylabel('Altitude (meters)')
    plt.plot(range(tStepsNum),HUD_tseries)
    
    fig4 = plt.figure()
    plt.title('Diffusivity\n$x_p$ = {0} km south'.format(P/1000))
    plt.xlabel('Distance (meters south from polar ocean shoreline)')
    plt.ylabel('Diffusivity (m$^2$/s)')
    plt.plot(xEdge,D)
    
    fig5 = plt.figure()
    plt.title('Ice Flux  $x_p$ = {0} km south'.format(P/1000))
    plt.xlabel('Distance (meters south from polar ocean shoreline)')
    plt.ylabel('Flux (m$^2$/s)')
    plt.plot(xEdge,M)
    
    fig6 = plt.figure()
    plt.title('Mass Balance  $x_p$ = {0} km south'.format(P/1000))
    plt.xlabel('Distance (meters south from polar ocean shoreline)')
    plt.ylabel('Mass Balance (m/yr)')
    plt.plot(x,B)

    return fig1, fig2, fig3, fig4, fig5, fig6

def onecycle():
    tStepsNum = 21000
    hystRun(tStepsNum,dRatio)
    plt.show()

def noGIA():
    tStepsNum = 21000
    hystRun_noGIA(tStepsNum)
    plt.show()

def multicycle():
    tStepsNum = 120000
    hystRun(tStepsNum,dRatio)
    plt.show()

def varyOption():
    winVary = tk.Tk()
    winVary.title('Varying Climate Point')
    varyInstruct = tk.Label(master = winVary, text = 'Choose Desired Integration Time')
    winVary.rowconfigure(0, minsize = 800, weight = 1)
    winVary.columnconfigure(1, minsize = 800, weight = 1)
    btn_OneCycle = tk.Button(master = winVary, text = '20 ka (one cycle)', command = onecycle)
    btn_noGIA = tk.Button(master = winVary, text = '20 ka (without bed adjustment)',command = noGIA)
    btn_MultiCycle = tk.Button(master = winVary, text = '120 ka (multiple cycles)', command = multicycle)
    
    varyInstruct.pack()
    btn_OneCycle.pack(pady=5)
    btn_noGIA.pack(pady=5)
    btn_MultiCycle.pack(pady=5)

def fixedOption():

    def fixedChoice():
        P = ent_P.get()
        fixedRun(x,H,ELA,h,HUD,float(P))
        plt.show()

    winFixed = tk.Tk()
    winFixed.title('Fixed Climate Point')
    fixInstruct = tk.Label(master = winFixed, text = 'Insert Desired Climate Point (m South from polar ocean)')
    winFixed.rowconfigure(0, minsize = 800, weight = 1)
    winFixed.columnconfigure(1, minsize = 800, weight = 1)
    P_entry = tk.Frame(master = winFixed)
    ent_P = tk.Entry(master = P_entry, width = 8)
    btn_FixedRun = tk.Button(master = winFixed, text = 'Run Model', command = fixedChoice)

    fixInstruct.pack(pady=10)
    P_entry.pack(pady=10)
    btn_FixedRun.pack(pady=3)
    ent_P.pack()

def moreInfo():
    winMore = tk.Tk()
    winMore.title('More Information')
    moretxt = 'Oerlemans, J. "Some basic experiments with a vertically-integrated ice sheet model." Tellus 33.1 (1981): 1-11'
    moreVar = tk.Message(winMore, text = moretxt)
    moreVar.pack()


win = tk.Tk()
win.title('1-dimensional Ice Sheet Model (North America)')
greeting = tk.Label(text = 'Introduction')

win.rowconfigure(0, minsize = 800, weight = 1)
win.columnconfigure(1, minsize = 800, weight = 1)

welcometxt = 'This 1-dimensional ice sheet model simulates the Laurentide ice sheet, ' \
    'which once covered much of North America during the last glacial. ' \
    'For simplicity, this model assumes a flat bedrock at sea level before isostatic adjustment. ' \
    'At the coast of the polar ocean, it is assumed all mass is lost through calving. \n \n' \
    'This model is based on Oerlemans (1981). \n \n' \
    'INSTRUCTIONS: \n' \
    '   (1) Choose whether you want the climate point to stay the same or vary. \n' \
    '   (2) Choose your desired climate point (if required). \n   (3) Click \"Run Model\" and wait for a graph to show up. \n' \
    '   (4) End the program to reset (I\'ll add a reset button at some point.) \n \n'\
    '(Important Note) The shoreline of the polar ocean is at 0 km South (75 degrees North), and negative numbers are points in the polar ocean. '\
    'The shoreline for the Gulf of Mexico is at 5000 km South (30 degrees North), and numbers greater than 5000 km are in the Gulf'\
    ' of Mexico.'
welcomeVar = tk.Message(win, text = welcometxt)
btn_fixed = tk.Button(master = win, text = 'Fixed', command = fixedOption)
btn_vary = tk.Button(master = win, text = 'Vary', command = varyOption)
btn_more = tk.Button(master= win, text = 'More Information', command = moreInfo)

greeting.pack()
welcomeVar.pack()
btn_fixed.pack(pady=5)
btn_vary.pack(pady=5)
btn_more.pack(pady=8)
win.mainloop()