# vSolar - External voltage of solar PV in Volts

# temperaturePV - temperature of solar PV in Kelvin



# ------------------------------------------------------------------------------------------------ #
def ElectricalPower(vExt):
# Convenience function to help find the maximum power point of a single junction solar photovoltaic
	
	global gTemperaturePV
	global gBandGapEnergy
	global gESpectrum
	global gEnergies
	
	photoCurrent = PhotoCurrent(gESpectrum, gEnergies, gBandGapEnergy, vExt, gTemperaturePV)
	power = photoCurrent * vExt
	return -1.0 * power
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def FindSolarPhotovoltaicPeakPowerPoint(energySpectrum, energies, bandGapEnergy, vExtStart, \
vExtBounds, temperaturePV, plotPowerCurve=False):
	
	from scipy.optimize import minimize
	from matplotlib.pyplot import figure, show, plot
	from numpy import arange

	global gTemperaturePV
	global gBandGapEnergy
	global gESpectrum
	global gEnergies
	
	if plotPowerCurve == True:
		vExtArray = arange(vExtBounds[0], vExtBounds[1], 0.01)
		powerArray = []

		for vExt in vExtArray:

			photoCurrent = PhotoCurrent(energySpectrum, energies, bandGapEnergy, vExt, temperaturePV)
			powerArray.append(photoCurrent * vExt)
	
		figure()
		plot(vExtArray, powerArray)
		show()
	
	# Reinitialize global variables
	gTemperaturePV = temperaturePV
	gBandGapEnergy = bandGapEnergy
	gESpectrum = energySpectrum
	gEnergies = energies

	
	# Find the best voltage for peak power production form the solar photovoltaic
	sol = minimize(ElectricalPower, vExtStart, bounds=[vExtBounds])


	vExtMax = sol.x[0]
	
	return vExtMax
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def SolarToElectricalEfficiency(energySpectrum, energies, bandGapEnergy, vExtIn, temperaturePV):
# 	Calculates the efficiency of an ideal photovoltaic in converting solar energy into something else	
	import pdb
	
	jPhoto = PhotoCurrent(energySpectrum, energies, bandGapEnergy, vExtIn, temperaturePV)
	
	power = IntegratePowerDensitySpectrum(energySpectrum, energies)
	
	solarToElectricalEfficiency = (jPhoto * abs(vExtIn))/power
	
	return solarToElectricalEfficiency
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def SolarAndCellCurrents(energySpectrum, energies, bandGapEnergy, vSolar, temperaturePV, \
vRedox, vBias):
	
	# Convenience function to calculate solar and electrochemical cell currents for
	# Solar to Hydrogen to Fuel Efficiency functions
	
	solarPower = IntegratePowerDensitySpectrum(energySpectrum, energies)
	
	# Calculate the current in the solar cell
	# From Physics of Solar Cells by Nelson
	photoCurrent = PhotoCurrent(energySpectrum, energies, bandGapEnergy, vSolar, temperaturePV)
	
	# Calculate the current in the electrochemical cell with the transformer equation
	# See Farshid #1 page 34.
	currentCell = (photoCurrent * vSolar) / (vRedox + vBias)
	
	return solarPower, currentCell
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ConvertWavelengthSpectrumToEnergySpectrum(wavelengthSpectrum, wavelengths):
	# Converts a wavelength spectrum into an energy spectrum
	
	import pdb
	from electrosynthesis.physicalconstants import speedOfLight as c, planckConstant as h
	from numpy import array
	
	energies = h*c/array(wavelengths)
	energies = energies[::-1]
	energySpectrum = array(wavelengthSpectrum)*array(wavelengths)**2/(h*c)
	energySpectrum = energySpectrum[::-1]
	
	return [energies, energySpectrum]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ImportSolarSpectrum(fileName):
# 	Imports the NREL solar spectrum. This code will likely not work for other files
	import pdb
	
	dataHandle = open(fileName, 'r')
	data = dataHandle.readlines()
	directSpectrum = []
	globalTilt = []
	etr = []
	wavelengths = []

	i = 2
	while i < len(data):
		line = data[i]
		lineData = line.strip().split(',')
		try:
			wavelengths.append(float(lineData[0])*1E-9)
			globalTilt.append(float(lineData[2])*1E9)
			etr.append(float(lineData[1])*1E9)
			directSpectrum.append(float(lineData[3])*1E9)
		except:
			pdb.set_trace()
		
		i += 1
		
	dataHandle.close()
	
	[energies, directESpectrum] = \
	ConvertWavelengthSpectrumToEnergySpectrum(directSpectrum, wavelengths)
	
	[energies, globalTiltESpectrum] = \
	ConvertWavelengthSpectrumToEnergySpectrum(globalTilt, wavelengths)
	
	[energies, etrESpectrum] = \
	ConvertWavelengthSpectrumToEnergySpectrum(etr, wavelengths)
	
	return [wavelengths, directSpectrum, globalTilt, etr, energies, directESpectrum, \
	globalTiltESpectrum, etrESpectrum]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def IntegratePowerDensitySpectrum(spectrum, wavelengths):
# This function will integrate a power density spectrum with respect to wavelength and return the
# total integrated power. Currently only integrates over the entire spectrum. If we want to change 
# limits of integration we will need to get clever. 

	from scipy.integrate import trapz
	
	totalPower = trapz(spectrum, x=wavelengths)
	
	return totalPower
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def CutoffWavelength(bandGapEnergy):
# Calculates the wavelength corresponding to a particular bandgap energy.
	
	from electrosynthesis.physicalconstants import speedOfLight, planckConstant
	import pdb
	
	energy = bandGapEnergy
	cutoff = (planckConstant*speedOfLight)/(energy)
		
	return cutoff
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def DarkCurrentDensity(energy, temperaturePV, chemicalPotential):
# Used to calculate the bias current to subtract from with the photovoltaic out of equilibrium 'Dark Current'
	
	import pdb
	from electrosynthesis.physicalconstants import speedOfLight as c, planckConstant as h, \
	elementaryCharge as e, boltzmannConstant as kB
	from numpy import pi, exp
	
	part1 = e*(2*pi*energy**2)/(h**3 * c**2)
	part2 = (1/(exp((energy - chemicalPotential)/(kB * temperaturePV)) - 1))
	part3 = 1/(exp(energy/(kB * temperaturePV)) - 1)
	
	darkCurrentDensity = e*(2*pi*energy**2)/(h**3 * c**2)*\
	(1/(exp((energy - chemicalPotential)/(kB * temperaturePV)) - 1) - \
	1/(exp(energy/(kB * temperaturePV)) - 1))
	
	return darkCurrentDensity
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def PhotoCurrent(intensities, energies, bandGapEnergy, voltageC, temperaturePV):
# Calculates the current output given the energy spectrum, the photon energies, bandgap energy, 
# the output potential voltageC and temperature of the cell temperaturePV, 
	
	import pdb
	from electrosynthesis.physicalconstants import speedOfLight as c, planckConstant as h, \
	elementaryCharge as e
	from scipy.integrate import trapz
	from numpy import array, float
	
	cutoffWavelength = CutoffWavelength(bandGapEnergy)
	cutoffEnergy = h * c/cutoffWavelength
	chemicalPotential = e * voltageC
	length = len(intensities)
	intensity2D = array([energies, intensities], float)
	intensity2D = intensity2D.transpose()
	intensity2D = sorted(intensity2D, key=lambda x:x[0])
	truncatedIntensity2D = []
	
	i = 0
	while i < len(intensity2D):	
		if intensity2D[i][0] > cutoffEnergy:
			truncatedIntensity2D.append(intensity2D[i])
		i += 1
	truncatedIntensity2D = array(truncatedIntensity2D, float)
	
	try:
		truncatedEnergies = truncatedIntensity2D.transpose()[0]
	except:
		pdb.set_trace()
	
	truncatedIntensities = truncatedIntensity2D.transpose()[1]
	darkCurrentD = DarkCurrentDensity(truncatedEnergies, temperaturePV, chemicalPotential)
 	
	try:
		derivative = e * truncatedIntensities/truncatedEnergies \
		- DarkCurrentDensity(truncatedEnergies, temperaturePV, chemicalPotential)
	except:
		pdb.set_trace()
		
	photoCurrent = trapz(derivative, x=truncatedEnergies)

	return photoCurrent
# ------------------------------------------------------------------------------------------------ #

