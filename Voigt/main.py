from matplotlib import pyplot as plt
import numpy as np

c_const = 3.0e10    #---in centimeters per second---#
T_ref = 296.0       #---in Kelvins---#
P_ref = 1.0         #---in atms---#
R_const = 8.31e7    #---in CGS system---#
Mass = 44           #---in gramm per mole---#

LeftLimit = 4500    #---in cm^-1---#
RightLimit = 6000   #---in cm^-1---#
IntensityLimit = 1.0e-24
GridStep = 0.005
GridSize = int((RightLimit - LeftLimit) / GridStep)

#-----------------------------#
R_const_SI = 8.31
Mass_SI = 0.044
Gamma = 7.0 / 5.0
g = 3.75

T1 = 220
T2 = 220 - 2.5 * 15
T3 = 220 - 2.5 * 20

P1 = 700            #---in Pa---#
C1 = P1 * ((Mass_SI * P1) / (R_const_SI * T1)) ** (-Gamma)

def PressureFunc(TemperatureP):
    return ((R_const_SI * TemperatureP / Mass_SI) ** Gamma / C1) ** (1.0 / (Gamma - 1))

def PaToAtm(PressureP):
    return PressureP / 101325.0

P1_Pa = PressureFunc(T1)
P2_Pa = PressureFunc(T2)
P3_Pa = PressureFunc(T3) * np.exp(-30000 / ((R_const_SI * T3) / (Mass_SI * g)))

def ConcentrationCGS(TemperatreP, PressureP):
    return PressureP / (1.23e-23 * TemperatreP) * 1.0e-6

print("%10.3e %10.3e %10.3e" % (ConcentrationCGS(T1, P1_Pa), ConcentrationCGS(T2, P2_Pa), ConcentrationCGS(T3, P3_Pa)))

P1 = PaToAtm(P1_Pa)
P2 = PaToAtm(P2_Pa)
P3 = PaToAtm(P3_Pa)

print("Pressure (atm):", P1, P2, P3)
# print("Pressure (atm):", P1, P2)
#-----------------------------#

FrequiencyArray, IntensityArray, GammaSelfArray = [], [], []

for line in open('Data.out', 'r'):
    values = [float(s) for s in line.split()]

    if (values[3] > LeftLimit and values[3] < RightLimit and values[4] > IntensityLimit):
        plt.plot((values[3], values[3]), (0.0, values[4]), 'k-', linewidth=0.5)

        FrequiencyArray.append(values[3])
        IntensityArray.append(values[4])
        GammaSelfArray.append(values[5])

plt.xlabel('Wavelength, cm^-1')
plt.ylabel('Intensity, cm^-1/(molec*cm^-2)')
plt.xlim([LeftLimit, RightLimit])
plt.ylim([0.0, 1.4e-20])
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

        if (LineI % 10 == 0):
            print(LineI, "out of", np.size(FrequiencyArray))

    SumCoeff *= ConcentrationCGS(TemperatureArray[i], PressureArray[i] * 101325.0)

    plt.title("Temperature:" + str(TemperatureArray[i]))

    plt.xlabel('Wavelength, cm^-1')
    plt.ylabel('Absorbtion, cm^-1')

    plt.xlim([LeftLimit, RightLimit])
    plt.ylim([0.0, 5.0e-3])

    X_values = np.arange(LeftLimit, RightLimit, GridStep)

    with open(str(TemperatureArray[i]) + '.txt', 'w') as f:
        for j in range(0, np.size(X_values)):
            f.write("%s\t" "%s\n" % (X_values[j], SumCoeff[j]))

    plt.plot(X_values, SumCoeff)
    plt.show()