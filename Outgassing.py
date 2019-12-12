import numpy as np 

IdealGasMolarVolume = 22.4 # liter/mol
Avogadro = 6.022E23 # particles/mol
Boltzmann = 8.6173303E-5 # eV/K
BoltzmannL = 83.14 # mBar*Liter/(mol*K)

def DoTimeConversion(TimeScale): 
    if TimeScale == 'Seconds': return 1
    if TimeScale == 'Minutes': return 60.0
    if TimeScale == 'Hours': return 3600.0
    if TimeScale == 'Days': return 3600.0*24.0
    if TimeScale == 'Weeks': return 3600.0*24.0*7.0

def GetImpuritiesVsTime(Data, TimeScale='Seconds'): 
    Impurities = []
    InitialImpurities = Data.InitialImpurities
    Time = [np.array([item for sublist in Data.Time for item in sublist])]
    DiffConstants = [np.array([item for sublist in Data.DiffConstants for item in sublist])]
    for ii, (Time, DiffConstant) in enumerate(zip(Time, DiffConstants)):   
        
        if ii==0:
            pass 
        else:
            Time = Time - np.max(Data.Time[ii-1])
        DiffConstant = DiffConstant * DoTimeConversion(TimeScale)
        Y = SolveDiffusionEquation(Time, DiffConstant, Data.Thickness, InitialImpurities)
        print('='*20)
        print(ii)
        print(InitialImpurities)
        print('time', Time)
        print('diff', DiffConstant)
        print(Y)
        YIndex = np.where(Y<Y[0]/Data.Constraints[ii])
        print(YIndex)
        print(Y[0]/Data.Constraints[ii])
        if np.shape(YIndex)[1] > 0: 
            YIndex = YIndex[0][0]
            Y[YIndex:] = Y[0]/Data.Constraints[ii]
            Data.ConstraintIndex.append(YIndex)
        else: 
            pass 
        Impurities.append(Y)
        InitialImpurities = Y[-1]
    
    for ii in range(1,len(Impurities)):
        Impurities[ii] = Impurities[ii] * (Impurities[ii-1][-1]/Impurities[ii][0])
    return np.array(Impurities)

def SolveDiffusionEquation(Time, Diff, Thickness, Conc): 
    value = 0.0 
    for N in range(0,1000,1): 
        factor1 = 1.0/((2.0*N+1.0)**2)
        factor2 = np.exp(-(np.pi*(2.0*N+1.0)/Thickness)**2*Diff*Time)
        comp = factor1*factor2 * (Conc*8.0*Thickness)/(np.pi**2*2.0)
        value = value + comp
    return value 

def GetFlowRateVsTime(Data, Units='#', TimeScale='Seconds'): 
    FlowRate = []
    InitialConcentration = [x[0]/(Data.Volume*1E3) for x in Data.Impurities]
    for ii, (Time, DiffConstant) in enumerate(zip(Data.Time, Data.DiffConstants)):
        if ii==0:
            pass 
        else:
            Time = Time - np.max(Data.Time[ii-1])

        DiffConstant = DiffConstant * DoTimeConversion(TimeScale)
        Y = GetFlowRate(Time, DiffConstant, Data.Thickness, InitialConcentration[ii], Data.Area)
        Y = Y/DoTimeConversion(TimeScale)
    
        if ii < len(Data.ConstraintIndex): 
            
            frac = 2.0
            YIndex = Data.ConstraintIndex[ii]
            # Y[:YIndex] += Y[YIndex]
            YFrac = int(YIndex//frac)
            XIndex = np.linspace(0, len(Y[YIndex//3:]), len(Y[YFrac:]))
            Y[YFrac:] = Y[YFrac]*np.exp(-0.7E-2*XIndex) + Y[YIndex]
            # Y[:YIndex] -= Y[YIndex]
            # Y[YFrac:] = 0
        FlowRate.append(Y)
    if Units == '#': pass 
    if Units == 'mBar Liter': 
        
        for ii,x in enumerate(FlowRate):
            FlowRate[ii] = FlowRate[ii] * (BoltzmannL * Data.Temp[ii]) / Avogadro
    return np.asarray(FlowRate)

def GetFlowRate(Time, Diff, Thickness, Conc, Area): 
    value = 0.0
    for N in range(0,3,1): 
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
    for Temp in Temperatures: 
        Diff = Data.Diffusion * np.exp(Data.ActivationEnergy/Boltzmann * ((1.0/293.15) - (1.0/Temp)))
        DiffTemp.append(Diff)
    return np.asarray(DiffTemp)