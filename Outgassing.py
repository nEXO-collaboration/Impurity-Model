import numpy as np 

IdealGasMolarVolume = 22.4 # liter/mol
Avogadro = 6.022E23 # particles/mol
Boltzmann = 8.6173303E-5 # eV/K

def DoTimeConversion(TimeScale): 
    if TimeScale == 'Seconds': return 1
    if TimeScale == 'Minutes': return 60.0
    if TimeScale == 'Hours': return 3600.0
    if TimeScale == 'Days': return 3600.0*24.0
    if TimeScale == 'Weeks': return 3600.0*24.0*7.0

def GetImpuritiesVsTime(Data, TimeScale='Seconds'): 
    Impurities = []
    InitialImpurities = Data.InitialImpurities
    for ii, (Time, DiffConstant) in enumerate(zip(Data.Time, Data.DiffConstants)): 
        DiffConstant = DiffConstant * DoTimeConversion(TimeScale)
        Y = SolveDiffusionEquation(Time, DiffConstant, Data.Thickness, InitialImpurities)
        Impurities.append(Y)
        InitialImpurities = Y[-1]
    for ii in range(1,len(Impurities)):
        Impurities[ii] = Impurities[ii] * (Impurities[ii-1][-1]/Impurities[ii][0])
    return np.array(Impurities)

def SolveDiffusionEquation(Time, Diff, Thickness, Conc): 
    value = 0.0 
    for N in range(0,100,1): 
        factor1 = 1.0/((2.0*N+1.0)**2)
        factor2 = np.exp(-(np.pi*(2.0*N+1.0)/Thickness)**2*Diff*Time)
        comp = factor1*factor2 * (Conc*8.0*Thickness)/(np.pi**2*2.0)
        value = value + comp
    return value 

def GetFlowRateVsTime(Data, Units='#', TimeScale='Seconds'): 
    FlowRate = []
    InitialConcentration = [x[0]/(Data.Volume*1000) for x in Data.Impurities]
    for ii, (Time, DiffConstant) in enumerate(zip(Data.Time, Data.DiffConstants)): 
        DiffConstant = DiffConstant * DoTimeConversion(TimeScale)
        Y = GetFlowRate(Time, DiffConstant, Data.Thickness, InitialConcentration[ii], Data.Area)
        FlowRate.append(Y)
    if Units == '#': pass 
    if Units == 'mBar Liter': 
        for ii,x in enumerate(FlowRate):
            FlowRate[ii] = FlowRate[ii] / (Boltzmann * Data.Temp[ii]) / Avogadro
    return np.array(FlowRate)

def GetFlowRate(Time, Diff, Thickness, Conc, Area): 
    value = 0.0 
    for N in range(0,2,1): 
        factor = np.exp(-(np.pi*(2.0*N+1.0)/Thickness)**2*Diff*Time)
        value = value + factor*(4.0*Conc*Diff)/Thickness
    value = value * Area 
    return value

def GetPartsConversion(Units):
    if Units == 'ppm': return 1E6
    if Units == 'ppb': return 1E9
    if Units == 'ppt': return 1E12

def GetInitialImpurities(Data, Units): 
    ImpurityVolume = Data.Volume * Data.Solubility * Data.Abundance
    ImpurityMass = ImpurityVolume/IdealGasMolarVolume * Data.MolarMass
    if Units == 'Mass':
        pass
    if 'pp' in Units: 
        ImpurityMass = ImpurityMass / Data.XeMass * GetPartsConversion(Units)
    if Units == '#': 
        ImpurityMass = ImpurityMass / Data.MolarMass * Avogadro
    return ImpurityMass

def GetDiffTemp(Data, Temperatures): 
    DiffTemp = []
    if len(Temperatures) == 1: 
        Temp = Temperatures[0]
        DiffTemp = Data.Diffusion * np.exp(Data.ActivationEnergy/Boltzmann * ((1.0/294.15) - (1.0/Temp)))
    else:
        for Temp in Temperatures: 
            Diff = Data.Diffusion * np.exp(Data.ActivationEnergy/Boltzmann * ((1.0/294.15) - (1.0/Temp)))
            DiffTemp.append(Diff)
        DiffTemp = np.asarray(DiffTemp)
    return DiffTemp