# ------------------------------------------------------------------------------------------------ #
# scaleup.py
# 
# Buz Barstow
# Created: 2019-06-14
#
# Functions for calculating the effect of scale up on efficiency
#
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Calculate_Film_Thickness_for_Diffusion(pressureHydrogen, \
diffusionHydrogen, solubilityHydrogen, electronRatePerCell, cellDensity):
# Calculating the required thickness of biofilm assuming H2 is required and limiting

	import pdb
	from rewiredcarbon.physicalconstants import elementaryCharge as e, avogadroNumber as N_A
	from math import sqrt
	from numpy import float
		
	hydFilmThickness = sqrt ((2 *pressureHydrogen * diffusionHydrogen * N_A)/ \
	(solubilityHydrogen * cellDensity * electronRatePerCell))
	
	#pdb.set_trace()
	return hydFilmThickness
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Calculate_Precise_Conductive_Matrix_Biofilm_Thickness(approximateBiofilmThickness, \
biofilmLayerHeight):
	
	from math import floor
	import pdb
	
	try:
		numberCellLayers = floor(approximateBiofilmThickness / biofilmLayerHeight)
	except:
		pdb.set_trace()
	
	preciseBiofilmThicnkess = numberCellLayers * biofilmLayerHeight
	
	return preciseBiofilmThicnkess, numberCellLayers
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Overpotential(energySpectrum, energies, bandGapEnergy, vExtIn, temperaturePV, temperatureCell, \
alpha, RedoxElectron, vRedox, pressureHydrogen, \
diffusionHydrogen, solubilityHydrogen, \
electronsPerFuel, carbonsPerFuel, co2FixEnzymeRate, co2FixEnzymePerCell, cellDensity, \
SelfExchangeCurrentDesnity):
# 	Calculates additional bias due to overpotential at H2 electrode, calculated by the Tafel 
# 	Equation:

	import pdb, sys
	from rewiredcarbon.physicalconstants import elementaryCharge as e, faradayConstant as F, \
	idealGasConstant as R
	from math import log
	
	solarToElectricalEfficiency= SolarToElectricalEfficiency(energySpectrum, energies, \
	bandGapEnergy, vExtIn, temperaturePV)
	
	electronRatePerCell= ElectronRatePerCell(electronsPerFuel, carbonsPerFuel, co2FixEnzymeRate, \
	co2FixEnzymePerCell)
	
	hydFilmThickness = HydFilmThickness(pressureHydrogen, diffusionHydrogen, solubilityHydrogen, \
	electronsPerFuel, carbonsPerFuel, co2FixEnzymeRate, co2FixEnzymePerCell, cellDensity)
	
	
	CellCurrent = ((solarToElectricalEfficiency*1000)/vRedox)/ e
	cellCultureVolume = (CellCurrent/(cellDensity * electronRatePerCell))
	
	supplyFaceArea = cellCultureVolume / hydFilmThickness

	currentDensityHydrogen = ((solarToElectricalEfficiency*1000)/vRedox)/ supplyFaceArea
	
	overpotential = abs(((R * temperatureCell)/(alpha * F* RedoxElectron)) \
	* log (currentDensityHydrogen/SelfExchangeCurrentDesnity))
	
	#pdb.set_trace()
	return overpotential
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def WireOverpotential(layerArea, packingFraction, diameterCell, lengthCell, diameterWire, \
resistanceFactor, energySpectrum, energies, bandGapEnergy, vExtIn, temperaturePV, \
electronsPerFuel, carbonsPerFuel, co2FixEnzymeRate, co2FixEnzymePerCell, pressureHydrogen, \
diffusionHydrogen, solubilityHydrogen, cellDensity, vRedox):
# 	Calculated the additional bias due to nanowire
	
	import pdb
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	from math import pi
	
	surfaceAreaCell= diameterCell * lengthCell
	wireArea = pi * (diameterWire/2)**2
	numberCellPerLayers = (layerArea/ surfaceAreaCell)* packingFraction
	
	solarToElectricalEfficiency= SolarToElectricalEfficiency(energySpectrum, energies, \
	bandGapEnergy, vExtIn, temperaturePV)
	
	electronRatePerCell= ElectronRatePerCell(electronsPerFuel, carbonsPerFuel, co2FixEnzymeRate, \
	co2FixEnzymePerCell)
	
	hydFilmThickness = HydFilmThickness(pressureHydrogen, diffusionHydrogen, solubilityHydrogen, \
	electronsPerFuel, carbonsPerFuel, co2FixEnzymeRate, co2FixEnzymePerCell, cellDensity)
	
	CellCurrent = ((solarToElectricalEfficiency*1000)/vRedox)/ e
	cellCultureVolume = (CellCurrent/(cellDensity * electronRatePerCell))
	supplyFaceArea = cellCultureVolume/ hydFilmThickness
	
	numberLayers = (solarToElectricalEfficiency*1000)/ (numberCellPerLayers * \
	electronRatePerCell * e)
	currentDensity= electronRatePerCell * e
	lengthWire= (diameterCell * numberLayers)/ supplyFaceArea
	wireOverpotential = (currentDensity * resistanceFactor * lengthWire)/ wireArea
	
	
	#pdb.set_trace()

	
	return wireOverpotential, numberCellPerLayers
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
####################################################################################################
# Functions for stirring power estimations
####################################################################################################
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ImportPowerNumberCurve(powerCurveFileName):

	import csv
	import pdb
	from numpy import float

	fHandle = open(powerCurveFileName, 'r')

	poolColumnToHeaderIndexDict = {}	
	datareader = csv.reader(fHandle)

	headers = next(datareader, None)

	column = {}
	for h in headers:
		column[h] = []

	for row in datareader:
		for h, v in zip(headers, row):
			column[h].append(float(v))

	fHandle.close()
	
	return column
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ImportAndExtraplotePowerNumberCurve(powerNumberCurveFileName, turbulentRegimeStart=1e4, \
showPowerNumberCurvePlot=False):

	from matplotlib.pyplot import show, loglog, xlabel, ylabel, legend, figure, title
	from numpy import arange, mean
	from scipy.interpolate import interp1d

	powerNumberCurveDataColumns = ImportPowerNumberCurve(powerNumberCurveFileName)
	reynoldsNumber = powerNumberCurveDataColumns['x']
	powerNumber = powerNumberCurveDataColumns['y']


	reynoldsNumberSorted = sorted(reynoldsNumber)
	lowestReynoldsNumber = reynoldsNumberSorted[0]
	lastReynoldsNumber = reynoldsNumberSorted[-1]

	reynoldsNumberInterp = arange(lowestReynoldsNumber,lastReynoldsNumber,1)
	powerNumberInterp = interp1d(reynoldsNumber, powerNumber)
	powerNumberInterpArray = powerNumberInterp(reynoldsNumberInterp)


	# Calculate the average of the last part of the power curve to enable extrapolation
	turbulentRegimeReynoldsNumber = arange(turbulentRegimeStart, lastReynoldsNumber, 1)
	turbulentPowerNumberInterp = powerNumberInterp(turbulentRegimeReynoldsNumber)
	turbulentPowerNumberMean = mean(turbulentPowerNumberInterp)

	reynoldsNumberExtrapolate = arange(lowestReynoldsNumber,lastReynoldsNumber*10,1)
	powerNumberExtrapolate = interp1d(reynoldsNumber, powerNumber, fill_value='extrapolate') 
	powerNumberExtrapolateArray = powerNumberExtrapolate(reynoldsNumberExtrapolate)

	if showPowerNumberCurvePlot == True:
		figure()
		loglog(reynoldsNumberExtrapolate, powerNumberExtrapolate(reynoldsNumberExtrapolate), \
		label='Interpolated and Extrapolated Function')
		loglog(reynoldsNumberInterp, powerNumberInterpArray, label='Original Data')
		xlabel("Reynolds Number")
		ylabel("Power Number")
		title("Extrapolation of Power Number")
		legend()

		show()
	
	return powerNumberExtrapolate
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Calculate_Power_From_Power_Curve_and_Rotational_Velocity(powerCurveFunction, \
rotationalVelocity, agitatorDiameter, fluidDensity=997, fluidViscosity=9.1e-4):

	agitatorReynoldsNumber = Calculate_Agitator_Reynolds_Number(rotationalVelocity, \
	agitatorDiameter, fluidDensity=fluidDensity, fluidViscosity=fluidViscosity)
	
	powerNumber = powerCurveFunction(agitatorReynoldsNumber)
	
	power = powerNumber * fluidDensity * (rotationalVelocity**3) * (agitatorDiameter**5)

	return agitatorReynoldsNumber, power
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Calculate_Power_From_Power_Curve_and_Reynolds_Number(powerCurveFunction, \
agitatorReynoldsNumber, agitatorDiameter, fluidDensity=997, fluidViscosity=9.1e-4):
	
	import sys
	import pdb
	from numpy import isclose
	import math
	import numpy
	
	agitatorRotationalVelocity = \
	Calculate_Agitator_Rotational_Velocity_from_Reynolds_Number(agitatorReynoldsNumber, \
	agitatorDiameter, fluidDensity=fluidDensity, fluidViscosity=fluidViscosity)
	
	# Cross check Reynolds number
	crossCheckAgitatorReynoldsNumber = Calculate_Agitator_Reynolds_Number(\
	agitatorRotationalVelocity, agitatorDiameter, fluidDensity=fluidDensity, \
	fluidViscosity=fluidViscosity)
	
	if not isclose(crossCheckAgitatorReynoldsNumber, agitatorReynoldsNumber):
		print("Reynolds number cross check is incorrect.")
		pdb.set_trace()
		sys.exit()
	
	powerNumber = powerCurveFunction(agitatorReynoldsNumber)
	
	try:
		power = powerNumber * fluidDensity * (agitatorRotationalVelocity**3) * (agitatorDiameter**5)
	except:
		pdb.set_trace()

	return agitatorRotationalVelocity, power
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Calculate_Agitator_Reynolds_Number(rotationalVelocity, agitatorDiameter, fluidDensity=997, \
fluidViscosity=9.1e-4):

	# rotationalVelocity = rotational velocity of agitator in revolutions per second
	# agitatorDiameter = diameter of agitator in meters
	# fluidDensity = density of the fluid that the agitator is agitating in kg m^-3. Default is 
	# density of water. 
	# fluidViscosity = viscosity of the fluid that the agitator is agitating in Pascal seconds.
	# Default is the viscosity of water at 25 ˚C. 
	
	agitatorReynoldsNumber = ((agitatorDiameter**2) * rotationalVelocity * fluidDensity) \
	/ fluidViscosity
	
	
	return agitatorReynoldsNumber
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Calculate_Agitator_Rotational_Velocity_from_Reynolds_Number(agitatorReynoldsNumber, \
agitatorDiameter, fluidDensity=997, fluidViscosity=9.1e-4):

	# agitatorDiameter = diameter of agitator in meters
	# fluidDensity = density of the fluid that the agitator is agitating in kg m^-3. Default is 
	# density of water. 
	# fluidViscosity = viscosity of the fluid that the agitator is agitating in Pascal seconds.
	# agitatorRotationalVelocity = returns agitator rotational velocity in rotations per second
	
	agitatorRotationalVelocity = (agitatorReynoldsNumber * fluidViscosity) \
	/ ((agitatorDiameter**2) * fluidDensity)

	return agitatorRotationalVelocity
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Calculate_Agitator_Diameter_from_Volume(agitatorVolume, tankDiameterToHeightRatio=1, \
agitatorDiameterToTankDiameterRatio=0.5):
	
	# See Geankoplis 3rd edition page 145
	# agitatorDiameterToTankDiameterRatio = ratio of agitator and tank diameters, 
	# default is agitator diameter = 0.5 * tank diameter
	
	from numpy import pi
	import math
	import numpy
	import pdb
	
	
	tankDiameter = (numpy.complex(4 * agitatorVolume / pi))**(1/3)
	
	agitatorDiameter = agitatorDiameterToTankDiameterRatio * numpy.complex(tankDiameter)
	
	# try:
# 		if math.isnan(agitatorDiameter) or numpy.isnan(agitatorDiameter):
# 			pdb.set_trace()
# 	except:
# 		pdb.set_trace()

	
	return agitatorDiameter
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Calculate_Tank_Cross_Sectional_Area_from_Volume(agitatorVolume, tankDiameterToHeightRatio=1):
	
	from numpy import pi
	import math
	import numpy
	import pdb
	
	tankDiameter = (4 * tankDiameterToHeightRatio * agitatorVolume / pi)**(1/3)
	tankCrossSectionalArea = pi * (tankDiameter/2)**2
	
	return tankCrossSectionalArea
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Calculate_Agitator_Power_from_Reynolds_Number_and_Tank_Volume(powerCurveFunction, \
agitatorReynoldsNumber, agitatorVolume, fluidDensity=997, fluidViscosity=9.1e-4, \
tankDiameterToHeightRatio=1, agitatorDiameterToTankDiameterRatio=0.5):

	import pdb
	import math
	import numpy

	agitatorDiameter = Calculate_Agitator_Diameter_from_Volume(agitatorVolume, \
	tankDiameterToHeightRatio=tankDiameterToHeightRatio, \
	agitatorDiameterToTankDiameterRatio=agitatorDiameterToTankDiameterRatio)
	
	# try:
# 		if math.isnan(agitatorDiameter) or numpy.isnan(agitatorDiameter):
# 			pdb.set_trace()
# 	except:
# 		pdb.set_trace()
	
	agitatorRotationalVelocity, power = Calculate_Power_From_Power_Curve_and_Reynolds_Number(\
	powerCurveFunction, agitatorReynoldsNumber, agitatorDiameter, fluidDensity=fluidDensity, \
	fluidViscosity=fluidViscosity)
	
# 	pdb.set_trace()
	
	# if math.isnan(agitatorDiameter) or numpy.isnan(agitatorDiameter):
# 		pdb.set_trace()

	return agitatorRotationalVelocity, power, agitatorDiameter
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Calculate_Agitator_Power_from_Hydrogen_Current(hydrogenCurrent, cellCultureVolume, \
pressureHydrogen=5066.25, temperature=293, bulkToPeakHydrogenConcentrationRatio=0.5, \
Aconstant=0.026, Bconstant=0.4, Cconstant=0.5, tankDiameterToHeightRatio=1):

	# See Buz Medium #2 pages 241-243 and Buz Medium #3 pages 1-2
	# Based on Richard Grenville's notes
	
	# pressureHydrogen - Partial pressure of hydrogen in Pascals in gas bubbles. Default is 5%.
	# bulkToPeakHydrogenConcentrationRatio - Ratio of bulk to peak concentration of hydrogen 
	# in liquid
	
	from rewiredcarbon.physicalconstants import avogadroNumber as NA
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	from rewiredcarbon.physicalconstants import henrysLawHydrogen
	from rewiredcarbon.physicalconstants import boltzmannConstant as kB
	from numpy import log
	import pdb
	
	
	# Peak concentration of Hydrogen in moles per cubic meter
	peakConcentrationHydrogen = pressureHydrogen / henrysLawHydrogen
	bulkConcentrationHydrogen = peakConcentrationHydrogen * bulkToPeakHydrogenConcentrationRatio
	deltaConcentrationHydrogen = peakConcentrationHydrogen - bulkConcentrationHydrogen
	
	# Transfer rate of H2 in moles per second
	massTransferRateHydrogen = hydrogenCurrent / ( 2 * e * NA )
	
	# Transfer of H2 in terms of cubic meters per second
	volumeTransferRateHydrogen = massTransferRateHydrogen * NA * kB * temperature / pressureHydrogen
	
	# Calculate the tank cross sectional area in square meters
	tankCrossSectionalArea = \
	Calculate_Tank_Cross_Sectional_Area_from_Volume(cellCultureVolume, \
	tankDiameterToHeightRatio=tankDiameterToHeightRatio)
	
	# Calculate the superficial gas velocity, the volume of gas passing through each square
	# meter of the top of the liquid. 
	superficialGasVelocity = volumeTransferRateHydrogen / tankCrossSectionalArea
	
	# Calculate the volumetric mass transfer coefficient
	kLaHydrogen = massTransferRateHydrogen \
	/ ( cellCultureVolume * (peakConcentrationHydrogen - bulkConcentrationHydrogen) )
	
	# Calculate the stirring power density in watts per cubic meter of liquid
	stirPowerDensity = \
	(kLaHydrogen / ( Aconstant * superficialGasVelocity**Cconstant ))**(1/Bconstant)
		
	stirPower = stirPowerDensity * cellCultureVolume
	
# 	pdb.set_trace()
	
	outputDict = {}
	outputDict['stirPower'] = stirPower
	outputDict['stirPowerDensity'] = stirPowerDensity
	outputDict['superficialGasVelocity'] = superficialGasVelocity
	outputDict['kLaHydrogen'] = kLaHydrogen
	outputDict['tankCrossSectionalArea'] = tankCrossSectionalArea
	outputDict['massTransferRateHydrogen'] = massTransferRateHydrogen
	outputDict['peakConcentrationHydrogen'] = peakConcentrationHydrogen
	outputDict['bulkConcentrationHydrogen'] = bulkConcentrationHydrogen
	outputDict['volumeTransferRateHydrogen'] = volumeTransferRateHydrogen
	
	
	return outputDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Calculate_Gas_Mixing_Parameters(hydrogenCurrent, cellCultureVolume, stirPower, \
pressureHydrogen=5066.25, temperature=293, bulkToPeakHydrogenConcentrationRatio=0.5, \
Aconstant=0.026, Bconstant=0.4, Cconstant=0.5, tankDiameterToHeightRatio=1):

	# See Buz Medium #2 pages 241-243 and Buz Medium #3 pages 1-2
	# Based on Richard Grenville's notes
	
	# pressureHydrogen - Partial pressure of hydrogen in Pascals in gas bubbles. Default is 5%.
	# bulkToPeakHydrogenConcentrationRatio - Ratio of bulk to peak concentration of hydrogen 
	# in liquid
	
	from rewiredcarbon.physicalconstants import avogadroNumber as NA
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	from rewiredcarbon.physicalconstants import henrysLawHydrogen
	from rewiredcarbon.physicalconstants import boltzmannConstant as kB
	import pdb
	
	
	# Peak concentration of Hydrogen in moles per cubic meter
	peakConcentrationHydrogen = pressureHydrogen / henrysLawHydrogen
	bulkConcentrationHydrogen = peakConcentrationHydrogen * bulkToPeakHydrogenConcentrationRatio
	deltaConcentrationHydrogen = peakConcentrationHydrogen - bulkConcentrationHydrogen
	
	# Transfer rate of H2 in moles per second
	massTransferRateHydrogen = hydrogenCurrent / ( 2 * e * NA )
	
	# Transfer of H2 in terms of cubic meters per second
	volumeTransferRateHydrogen = massTransferRateHydrogen * NA * kB * temperature / pressureHydrogen
	
	# Calculate the tank cross sectional area in square meters
	tankCrossSectionalArea = \
	Calculate_Tank_Cross_Sectional_Area_from_Volume(cellCultureVolume, \
	tankDiameterToHeightRatio=tankDiameterToHeightRatio)
	
	# Calculate the superficial gas velocity (usg), the volume of gas passing through each square
	# meter of the top of the liquid. 
	superficialGasVelocity = volumeTransferRateHydrogen / tankCrossSectionalArea
	
	# Calculate the stirring power density in watts per cubic meter of liquid
	stirPowerDensity = stirPower / cellCultureVolume
	
	# Calculate the volumetric mass transfer coefficient using the geometry independent formula
	# (concentration dependent)
	kLaHydrogenGeoInd = massTransferRateHydrogen \
	/ ( cellCultureVolume * (peakConcentrationHydrogen - bulkConcentrationHydrogen) )
	
	# Calculate the volumetric mass transfer coefficient (kLa) using the geometry dependent formula
	kLaHydrogenGeoDep = Aconstant * (stirPowerDensity**Bconstant) * (superficialGasVelocity**Cconstant) 
	
	
	kLaHydrogenGeoIndToDepRatio = kLaHydrogenGeoInd / kLaHydrogenGeoDep
	
	outputDict = {}
	outputDict['stirPowerDensity'] = stirPowerDensity
	outputDict['superficialGasVelocity'] = superficialGasVelocity
	outputDict['kLaHydrogen'] = kLaHydrogenGeoInd
	outputDict['tankCrossSectionalArea'] = tankCrossSectionalArea
	outputDict['massTransferRateHydrogen'] = massTransferRateHydrogen
	outputDict['peakConcentrationHydrogen'] = peakConcentrationHydrogen
	outputDict['bulkConcentrationHydrogen'] = bulkConcentrationHydrogen
	outputDict['volumeTransferRateHydrogen'] = volumeTransferRateHydrogen
	outputDict['kLaHydrogenGeoInd'] = kLaHydrogenGeoInd
	outputDict['kLaHydrogenGeoDep'] = kLaHydrogenGeoDep
	outputDict['kLaHydrogenGeoIndToDepRatio'] = kLaHydrogenGeoIndToDepRatio
	
	return outputDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Calculate_Density_and_Area_for_Target_Efficiency_to_Peak_Efficiency_Ratio(efficiencyDict, \
targetEfficiencyRatio, peakEfficiency, solarPVArea):

	# See Buz Medium #3 page 69
	
	# This function calculates the cell density, and tank area at a target efficiency (a fraction
	# of the peak efficiency that you could achieve if you didn't have to stir). 
	# Also calculates how big this tank cross section is relative to the area solar PV of a 
	# solar PV supplying it.

	from scipy.interpolate import interp1d
	from numpy import float, array
	import pdb
	
	effTotalElectricalToFuel = array(efficiencyDict['effTotalElectricalToFuel'])
	
	try:
		effTotalToPeakEffRatio = effTotalElectricalToFuel / float(peakEfficiency)
	except:
		pdb.set_trace()
		
		
	cellDensity = efficiencyDict['cellDensity']

	cellDensityToEffInterpFunction = interp1d(effTotalToPeakEffRatio, cellDensity)
	
	cellDensityAtTargetEfficiencyRatio = \
	float(cellDensityToEffInterpFunction(targetEfficiencyRatio))
	
	tankCrossSectionalArea = efficiencyDict['tankCrossSectionalArea']
	
	tankAreaInterpFunction = interp1d(cellDensity, tankCrossSectionalArea)
	
	tankAreaAtTargetEfficiencyRatio = \
	tankAreaInterpFunction(cellDensityAtTargetEfficiencyRatio)
	
	tankAreaRelativeToSolarPVArea = tankAreaAtTargetEfficiencyRatio / solarPVArea
	
	outputDict = {}
	outputDict['tankAreaRelativeToSolarPVArea'] = tankAreaRelativeToSolarPVArea
	outputDict['tankCrossSectionalAreaAtTargetEfficiencyRatio'] = tankAreaAtTargetEfficiencyRatio
	outputDict['cellDensityAtTargetEfficiencyRatio'] = cellDensityAtTargetEfficiencyRatio
	
	return outputDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Calculate_Density_and_Area_for_Target_to_Peak_Efficiency_Ratio_for_Array_of_Input_Powers(\
scenarioDict, efficienciesDict, efficienciesDictKeysForCalc, targetEfficiencyRatio, \
peakElectricalToFuelEfficiency, solarConstant=1000):

	# Solar constant in W m^-2

	# This function repeatedly calls 
	# Calculate_Density_and_Area_for_Target_Efficiency_to_Peak_Efficiency_Ratio
	# and calculates the system footprint needed to achieve a target efficiency for a range 
	# of input powers.
	
	import pdb
	
	inputPowerArray = []
	totalElectricalPowerArray = []
	solarPVAreaArray = []
	tankAreaRelativeToSolarPVAreaArray = []
	tankCrossSectionalAreaArray = []
	cellDensityAtTargetEfficiencyRatioArray = []
	
	collectedOutputDict = {}


	for key in efficienciesDictKeysForCalc:
	
		solarPower = float(scenarioDict[key]['totalInputPower'])
		totalElectricalPower = float(scenarioDict[key]['totalElectricalPower'])
		solarPVArea = solarPower / solarConstant
	
		efficiencyDict = efficienciesDict[key]
	
		outputDict = \
		Calculate_Density_and_Area_for_Target_Efficiency_to_Peak_Efficiency_Ratio(efficiencyDict, \
		targetEfficiencyRatio, peakElectricalToFuelEfficiency, solarPVArea)
	
		tankAreaRelativeToSolarPVArea = outputDict['tankAreaRelativeToSolarPVArea']
		tankCrossSectionalArea = outputDict['tankCrossSectionalAreaAtTargetEfficiencyRatio']
		cellDensityAtTargetEfficiencyRatio = outputDict['cellDensityAtTargetEfficiencyRatio']
	
		inputPowerArray.append(solarPower)
		solarPVAreaArray.append(solarPVArea)
		totalElectricalPowerArray.append(totalElectricalPower)
		tankAreaRelativeToSolarPVAreaArray.append(tankAreaRelativeToSolarPVArea)
		tankCrossSectionalAreaArray.append(tankCrossSectionalArea)
		cellDensityAtTargetEfficiencyRatioArray.append(cellDensityAtTargetEfficiencyRatio)
		
	
	collectedOutputDict['inputPowerArray'] = inputPowerArray
	collectedOutputDict['solarPVAreaArray'] = solarPVAreaArray
	collectedOutputDict['totalElectricalPowerArray'] = totalElectricalPowerArray
	collectedOutputDict['tankAreaRelativeToSolarPVAreaArray'] = tankAreaRelativeToSolarPVAreaArray
	collectedOutputDict['tankCrossSectionalAreaArray'] = tankCrossSectionalAreaArray
	collectedOutputDict['cellDensityAtTargetEfficiencyRatioArray'] = \
	cellDensityAtTargetEfficiencyRatioArray
		
	
	return collectedOutputDict
# ------------------------------------------------------------------------------------------------ #
