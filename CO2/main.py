from matplotlib import pyplot as plt
import numpy as np

c_const = 3.0e10    #---in centimeters per second---#
T_ref = 296.0       #---in Kelvins---#
P_ref = 1.0         #---in atms---#
R_const = 8.31e7    #---in CGS system---#
Mass = 44           #---in gramm per mole---#

T_custom = 300
P_custom = 0.007

LeftLimit = 4500    #---in cm^-1---#
RightLimit = 6000   #---in cm^-1---#
IntensityLimit = 1.0e-26
GridStep = 0.1
GridSize = int((RightLimit - LeftLimit) / GridStep)

#-----------------------------#
R_const_SI = 8.31
Mass_SI = 0.044
Gamma = 7.0 / 5.0

T1 = 220
T2 = 180
T3 = 130

P1 = 700            #---in Pa---#
C1 = P1 * ((Mass_SI * P1) / (R_const_SI * T1)) ** (-Gamma)

def PressureFunc(TemperatureP):
    return ((R_const_SI * TemperatureP / Mass_SI) ** Gamma / C1) ** (1.0 / (Gamma - 1))

def PaToAtm(PressureP):
    return PressureP / 101325.0

P1 = PaToAtm(PressureFunc(T1))
P2 = PaToAtm(PressureFunc(T2))
P3 = PaToAtm(PressureFunc(T3))

print("Pressure (atm):", P1, P2, P3)
#-----------------------------#

FrequiencyArray, IntensityArray, GammaSelfArray = [], [], []

for line in open('TheFuckingNewest.out', 'r'):
    values = [float(s) for s in line.split()]

    if (values[3] > LeftLimit and values[3] < RightLimit and values[4] > IntensityLimit):
        plt.plot((values[3], values[3]), (0.0, values[4]), 'k-', linewidth=0.5)

        FrequiencyArray.append(values[3])
        IntensityArray.append(values[4])
        GammaSelfArray.append(values[5])

plt.xlabel('Wavelength, cm^-1')
plt.ylabel('Intensity, cm^-1/(molec*cm^-2)')
plt.xlim([LeftLimit, RightLimit])
plt.ylim([0.0, 1.3e-20])
plt.show()
#-----------------------------#

FrequiencyArray = np.asarray(FrequiencyArray)
IntensityArray = np.asarray(IntensityArray)
GammaSelfArray = np.asarray(GammaSelfArray)

print("Number of frequencies:", np.size(FrequiencyArray))

#---Profiles---#
def LorentzHWHM(TemperatureP, PressureP, GammaSelfP):
    return GammaSelfP * PressureP * np.sqrt(T_ref / TemperatureP)

def DopplerHWHM(TemperatureP, FrequencyP):
    return FrequencyP / c_const * np.sqrt(2 * R_const * TemperatureP * np.log(2) / Mass)

def VoigtHWHM(LorenzP, DopplerP):
    TempL = LorenzP + np.sqrt(LorenzP ** 2 + 4.0 * DopplerP ** 2)
    return  0.5 * TempL + 0.05 * (1 - 2.0 * LorenzP / TempL)
#---Profiles---#

def Absorption(MiddleFrequencyP, IntensityP, GammaLorentzP, GammaVoigtP, OutputArray):
    for i in range(0, GridSize):
        XiL = GammaLorentzP / GammaVoigtP
        EtaL = (LeftLimit + GridStep * i - MiddleFrequencyP) / GammaVoigtP

        CoeffitientL =  (IntensityP * np.sqrt(np.log(2))) / (np.sqrt(np.pi) * GammaVoigtP) *    \
                        (1.0 - XiL) * np.exp(-EtaL ** 2 * np.log(2)) +                          \
                        (IntensityP * XiL) / (np.pi * GammaVoigtP * (1.0 + EtaL ** 2))

        OutputArray[i] += CoeffitientL

    return

TemperatureArray = [T1, T2, T3]
PressureArray = [P1, P2, P3]

for i in range(0, 3):
    SumCoeff = np.zeros(GridSize)

    for LineI in range(0, np.size(FrequiencyArray)):
        LorentzValueL = LorentzHWHM(TemperatureArray[i], PressureArray[i], GammaSelfArray[LineI])
        DoppleValueL = DopplerHWHM(TemperatureArray[i], FrequiencyArray[LineI])
        VoigtValueL = VoigtHWHM(LorentzValueL, DoppleValueL)

        Absorption(FrequiencyArray[LineI], IntensityArray[LineI], LorentzValueL, VoigtValueL, SumCoeff)

        if (LineI % 100 == 0):
            print(LineI, "out of", np.size(FrequiencyArray))

    plt.title("Temperature:" + str(TemperatureArray[i]))

    plt.xlabel('Wavelength, cm^-1')
    plt.ylabel('Intensity, cm^-1/(molec*cm^-2)')

    plt.xlim([LeftLimit, RightLimit])
    plt.ylim([0.0, 1.3e-20])

    X_values = np.arange(LeftLimit, RightLimit, GridStep)

    plt.plot(X_values, SumCoeff)
    plt.show()