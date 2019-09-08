import numpy as np
import matplotlib.pyplot as plt
import Library as Lib 
import Outgassing as Out
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import math

params = {'text.usetex':True,'font.size':12,'font.family':'serif'}
plt.rcParams.update(params) 
colors = ['#1f78b4', '#e66101', '#33a02c', '#984ea3', 'black']
Linestyles = ['-', '--', ':', ';']

class System(): 
    def __init__(self, Setup, Material, Solute, Version):
        self.Name = Setup
        self.Material = Material
        self.Solute = Solute
        self.Version = Version
        
        self.Diffusion = Lib.Material.get(Material).get(Solute).get('Diffusion Constant')
        self.Solubility = Lib.Material.get(Material).get(Solute).get('Solubility')
        self.ActivationEnergy = Lib.Material.get(Material).get(Solute).get('Activation Energy')
        self.Abundance =  Lib.Gas.get(Solute).get('Abundance in Air')
        self.MolarMass =  Lib.Gas.get(Solute).get('Molar Mass')

        self.XeMass = Lib.System.get(Setup).get('Xenon Mass')
        self.Volume = Lib.System.get(Setup).get(Material).get(Version).get('Volume')
        self.Area = Lib.System.get(Setup).get(Material).get(Version).get('Area')
        self.Thickness = Lib.System.get(Setup).get(Material).get(Version).get('Thickness')

    def Print(self):
        Attributes = vars(self)
        print('\n '.join("%s: %s" % item for item in Attributes.items()))

def GetTimeStamps(Points, Spacing, TimeScale): 
    X = []
    for ii, x in enumerate(Points): 
        if ii == len(Points)-1: break 
        X.append(np.linspace(Points[ii], Points[ii+1], int((Points[ii+1] - Points[ii]) / Spacing + 1)))
    return np.array(X)

def PlotImpuritiesVsTime(Data, XRange=0, YRange=0, XTicks=0)):
    fig = plt.figure(figsize=(10,7))
    ax = fig.gca()
    plt.xlabel('Time [hours]', fontsize=16)
    plt.ylabel(r'Total Number of Impurities', fontsize=16)
    plt.yscale('log')
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    ax.grid(b=True, which='major', color='grey', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')

    for jj, data in enumerate(Data):
        for ii,(X,Y) in enumerate(zip(data.Time, data.Impurities)):
            plt.plot(X, Y, label=data.Labels[ii], color=colors[jj], linewidth=2.0, linestyle=Linestyles[ii])

    if XTicks == 0:
        XTicks = int(np.max([np.max(x.Time) for x in Data ])/10)
    else:
        XTicks = XTicks 
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.xaxis.set_major_locator(MultipleLocator(XTicks))

    YMaxExp = np.max([np.max(x.FlowRate) for x in Data ])
    YMaxExp = math.ceil(math.log10(YMaxExp))
    if XRange == 0:
        plt.xlim(np.min(Data[0].Time[0]), np.max(Data[0].Time[-1]))
    else:
        plt.xlim(XRange[0], XRange[1])
    if YRange == 0:
        plt.ylim(ymin=1E10, ymax = 1*10**YMaxExp)
    else:
        plt.ylim(YRange[0], YRange[1])
    ax.legend(loc='lower right', fontsize=10)
    fig.tight_layout()

def PlotFlowRateVsTime(Data, XRange=0, YRange=0, XTicks=0):
    fig = plt.figure(figsize=(10,7))
    ax = fig.gca()
    plt.xlabel('Time [Hours]', fontsize=16)
    plt.ylabel(r'Outgassing Rate [mBar$\,\cdot\,$Liter/s]', fontsize=16)
    plt.yscale('log')
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    ax.grid(b=True, which='major', color='grey', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')

    for jj,data in enumerate(Data):
        for ii,(X,Y) in enumerate(zip(data.Time, data.FlowRate)):
            plt.plot(X, Y, label=data.Labels[ii], color=colors[jj], linewidth=2.0, linestyle=Linestyles[ii])

    if XTicks == 0:
        XTicks = int(np.max([np.max(x.Time) for x in Data ])/10)
    else:
        XTicks = XTicks 
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.xaxis.set_major_locator(MultipleLocator(XTicks))

    YMaxExp = np.max([np.max(x.FlowRate) for x in Data ])
    YMaxExp = math.ceil(math.log10(YMaxExp))
    if XRange == 0:
        plt.xlim(np.min(Data[0].Time[0]), np.max(Data[0].Time[-1]))
    else:
        plt.xlim(XRange[0], XRange[1])
    if YRange == 0:
        plt.ylim(ymin=1E-22, ymax = 1*10**YMaxExp)
    else:
        plt.ylim(YRange[0], YRange[1])

    ax.legend(loc='lower right', fontsize=10)
    fig.tight_layout()

def DoModelling(Systems, Labels, Temperature, Time):
    for ii, System in enumerate(Systems): 
        # This will output all important parameters that have been imported into each model
        # System.Print()

        # Define the different temperatures for which to calculate outgassing
        System.Temp = Temperature

        # Calculate the diffusion constants for the above defined temperatures using the Arrhenius equation
        System.DiffConstants = Out.GetDiffTemp(System, Temperatures=System.Temp)    

        # Get the initial number of impurities from model parameters. Can define '#', 'ppm,ppb,ppt' or 'Mass' for units 
        System.InitialImpurities = Out.GetInitialImpurities(System, '#')

        # Forward above defined labels, that will be later used for the plot legend.
        System.Labels = Labels[ii]

        System.Time = Time

        # Calculate the number of impurities left in the sample vs time by using the solution to the diffusion equation.
        System.Impurities = Out.GetImpuritiesVsTime(Data=System, TimeScale='Hours')

        # Calculate the outgassing rate vs time using Fick's 1st law 
        System.FlowRate = Out.GetFlowRateVsTime(Data=System, Units='mBar Liter', TimeScale='Hours')

def GetLabels(Systems, Temperature): 
    Labels = []
    for ii, System in enumerate(Systems): 
        label = [] 
        for jj, temp in enumerate(Temperature):
            first = '%s, Thickness: %s cm, Temperature: %d K' % (System.Version, System.Thickness, temp)
            label.append(first)
        Labels.append(label)
    return Labels 


if __name__ == '__main__':

    # Define parameters for each model. Details about which options are available are in Library.py
    # Units for all parameters are defined in Library.py.
    S1 = System(Setup='YLXPS', Material='Teflon', Solute='Oxygen', Version='EXO-Teflon')
    S2 = System(Setup='YLXPS', Material='Teflon', Solute='Oxygen', Version='Stock-Teflon')
    S3 = System(Setup='YLXPS', Material='Teflon', Solute='Oxygen', Version='Columbia-Teflon')

    # Bundle above defined models together to do calculations for all of them 
    Systems = [S1,S2,S3]
    

    # Define the different temperatures for which to calculate outgassing
    Temperature = [295, 200]

    Labels = [['EXO-200 Teflon'], ['Stock Room Teflon'], ['Columbia Setup Teflon']]

    # Sets time regions for the above defined temperature values
    # Points goes from 0 to first value, first value to second value and so on. Ex: [0,10,15] gives two regions, 0-10 and 10-15
    # TimeScale defines in which units to plot the x-Axis later.
    Time = GetTimeStamps(Points=[0,100, 200], Spacing=1.0, TimeScale='Hours')

    Labels = GetLabels(Systems, Temperature)

    # Executing all steps to get impuritie numbers and outgassing rate as a function of time 
    DoModelling(Systems, Labels, Temperature, Time)

    # Plotting the impurities and outgassing rates for all above defined models and saving the plots. 
    PlotImpuritiesVsTime(Systems)
    # plt.savefig('impurities.pdf')
    PlotFlowRateVsTime(Systems)
    # plt.savefig('outgassing_rate.pdf')
    plt.show()
