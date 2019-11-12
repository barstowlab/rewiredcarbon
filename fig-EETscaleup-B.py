# ------------------------------------------------------------------------------------------------ #
def Calculate_BiolfimThickness_and_Area_for_Target_Efficiency_to_Peak_Efficiency_Ratio(\
efficiencyDict, targetEfficiencyRatio, peakEfficiency, solarPVArea):

	# See Buz Medium #3 page 85, 88-89
	
	# This function calculates the biofilm thickness and area at a target efficiency (a fraction
	# of the peak efficiency that you could achieve if you didn't have the biofilm). 
	# Also calculates how big this biofilm area is relative to the area of a 
	# solar PV supplying it.

	from scipy.interpolate import interp1d
	from numpy import interp
	from numpy import float, array, log10
	import pdb
	
	effTotalElectricalToFuel = array(efficiencyDict['effTotalElectricalToFuel'])
	
	try:
		effTotalToPeakEffRatio = effTotalElectricalToFuel / float(peakEfficiency)
	except:
		pdb.set_trace()
		
		
	preciseBiofilmThickness = efficiencyDict['preciseBiofilmThickness']
	
	logPreciseBiofilmThickness = log10(preciseBiofilmThickness)

	logBiofilmThicknessToEffInterpFunction = interp1d(effTotalToPeakEffRatio, logPreciseBiofilmThickness)
	
	logBiofilmThicknessAtTargetEfficiencyRatio = \
	float(logBiofilmThicknessToEffInterpFunction(targetEfficiencyRatio))
	
	biofilmThicknessAtTargetEfficiencyRatio = 10**(logBiofilmThicknessAtTargetEfficiencyRatio)
	
	areaBiofilm = efficiencyDict['areaBiofilm']
	
	areaBiofilmInterpFunction = interp1d(preciseBiofilmThickness, areaBiofilm)
	
	areaBiofilmAtTargetEfficiencyRatio = \
	float(areaBiofilmInterpFunction(biofilmThicknessAtTargetEfficiencyRatio))
	
	areaBiofilmRelativeToSolarPVArea = areaBiofilmAtTargetEfficiencyRatio / solarPVArea
	
	outputDict = {}
	outputDict['areaBiofilmRelativeToSolarPVArea'] = areaBiofilmRelativeToSolarPVArea
	outputDict['areaBiofilmAtTargetEfficiencyRatio'] = areaBiofilmAtTargetEfficiencyRatio
	outputDict['biofilmThicknessAtTargetEfficiencyRatio'] = biofilmThicknessAtTargetEfficiencyRatio
	outputDict['areaBiofilmInterpFunction'] = areaBiofilmInterpFunction
	outputDict['logBiofilmThicknessToEffInterpFunction'] = logBiofilmThicknessToEffInterpFunction
	
	
	
	return outputDict
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Calculate_BiolfimThickness_and_Area_for_Target_to_Peak_Eff_for_Array_of_Input_Resitivities(\
scenarioDict, efficienciesDict, efficienciesDictKeysForCalc, targetEfficiencyRatio, \
peakElectricalToFuelEfficiency, solarConstant=1000):

	# Solar constant in W m^-2

	# This function repeatedly calls 
	# Calculate_BiolfimThickness_and_Area_for_Target_Efficiency_to_Peak_Efficiency_Ratio
	# and calculates the system footprint (biofilm area) and biofilm thickness needed to achieve 
	# a target efficiency for a range of input biofilm resistitivies.
	
	import pdb
	
	inputResistivityArray = []
	conductivityArray = []
	solarPVAreaArray = []
	inputPowerArray = []
	areaBiofilmRelativeToSolarPVAreaArray = []
	areaBiofilmAtTargetEfficiencyRatioArray = []
	biofilmThicknessAtTargetEfficiencyRatioArray = []
	totalElectricalPowerArray = []
	
	
	collectedOutputDict = {}


	for key in efficienciesDictKeysForCalc:
		
		efficiencyDict = efficienciesDict[key]
		
		solarPower = float(scenarioDict[key]['totalInputPower'])
		totalElectricalPower = float(scenarioDict[key]['totalElectricalPower'])
		solarPVArea = solarPower / solarConstant
		
		# Convert the input resistivity into ohm meters from ohm centimeters
		resistivity = float(scenarioDict[key]['resistivityBiofilm'])*1e-2
		# Convert to Siemens per meter
		conductivity = 1/resistivity
	
		outputDict = \
		Calculate_BiolfimThickness_and_Area_for_Target_Efficiency_to_Peak_Efficiency_Ratio(\
		efficiencyDict, targetEfficiencyRatio, peakElectricalToFuelEfficiency, solarPVArea)
		
		areaBiofilmRelativeToSolarPVArea = outputDict['areaBiofilmRelativeToSolarPVArea']
		areaBiofilmAtTargetEfficiencyRatio = outputDict['areaBiofilmAtTargetEfficiencyRatio']
		biofilmThicknessAtTargetEfficiencyRatio = \
		outputDict['biofilmThicknessAtTargetEfficiencyRatio']
		
		inputPowerArray.append(solarPower)
		solarPVAreaArray.append(solarPVArea)
		totalElectricalPowerArray.append(totalElectricalPower)
		inputResistivityArray.append(resistivity)
		conductivityArray.append(conductivity)
		
		areaBiofilmRelativeToSolarPVAreaArray.append(areaBiofilmRelativeToSolarPVArea)
		areaBiofilmAtTargetEfficiencyRatioArray.append(areaBiofilmAtTargetEfficiencyRatio)
		biofilmThicknessAtTargetEfficiencyRatioArray.append(biofilmThicknessAtTargetEfficiencyRatio)
		
	
	collectedOutputDict['inputPowerArray'] = inputPowerArray
	collectedOutputDict['solarPVAreaArray'] = solarPVAreaArray
	collectedOutputDict['totalElectricalPowerArray'] = totalElectricalPowerArray
	collectedOutputDict['inputResistivityArray'] = inputResistivityArray
	collectedOutputDict['conductivityArray'] = conductivityArray
	
	collectedOutputDict['areaBiofilmRelativeToSolarPVAreaArray'] = \
	areaBiofilmRelativeToSolarPVAreaArray
		
	collectedOutputDict['areaBiofilmAtTargetEfficiencyRatioArray'] = \
	areaBiofilmAtTargetEfficiencyRatioArray
	
	collectedOutputDict['biofilmThicknessAtTargetEfficiencyRatioArray'] = \
	biofilmThicknessAtTargetEfficiencyRatioArray
	
	
	return collectedOutputDict
# ------------------------------------------------------------------------------------------------ #


from rewiredcarbon.scenario import ImportScenarioTable, CalculateScenarioEfficiencies, \
Plot_Efficiency_Scattergraphs, Plot_Efficiency_Scattergraphs_2,  Export_Efficiency_Scattergraphs

from rewiredcarbon.vectorOutput import generateOutputMatrixWithHeaders, writeOutputMatrix
from rewiredcarbon.utils import ensure_dir

from copy import deepcopy

scenarioTableFileName = 'input/fig-EETscaleup-B.csv'

outputFileDirname = 'output/fig-EETscaleup/fig-EETscaleup-B/'
outputFilename = outputFileDirname + '/' + 'fig-EETscaleup-B.csv'
ensure_dir(outputFilename)


# Solar constant in W m^-2
solarConstant = 1000

# Calculate the scenario efficiencies
scenarioDict = ImportScenarioTable(scenarioTableFileName)
efficienciesDict = CalculateScenarioEfficiencies(scenarioDict, mode='scattergraph')

# Calculate out the peak efficiency
peakElectricalToFuelEfficiency = max(efficienciesDict['Reference']['effAvailElectricalToFuel'])

# Remove the reference scenario from consideration
efficienciesDictKeysForCalc = list(efficienciesDict.keys())
efficienciesDictKeysForCalc.remove('Reference')


# Define target efficiencies
targetEfficiencyRatios = [0.5, 0.75, 0.95]
collectedOutputDict = {}

# Calculate the system footprint and biofilm thickness data for different target efficiency ratios
for targetEfficiencyRatio in targetEfficiencyRatios:

	outputDict = \
	Calculate_BiolfimThickness_and_Area_for_Target_to_Peak_Eff_for_Array_of_Input_Resitivities(\
	scenarioDict, efficienciesDict, efficienciesDictKeysForCalc, targetEfficiencyRatio, \
	peakElectricalToFuelEfficiency, solarConstant=solarConstant)

	collectedOutputDict[str(targetEfficiencyRatio)] = deepcopy(outputDict)


# ------------------------------------------------------------------------------------------------ #
# Plot out the system footprint and biofilm thickness data for different target efficiency ratios



keys = collectedOutputDict.keys()
vectorList = []
headerList = []

# ------------------------------------------------------------------------------------------------ #
# Plot out conductivity versus biofilm area
figure()
for key in keys:
	outputDict = collectedOutputDict[key]
	# Convert to S cm^-1 from S m^-1 
	conductivityArray = array(outputDict['conductivityArray']) / 100 
	areaBiofilmAtTargetEfficiencyRatioArray = outputDict['areaBiofilmAtTargetEfficiencyRatioArray']
	areaBiofilmRelativeToSolarPVAreaArray = outputDict['areaBiofilmRelativeToSolarPVAreaArray']
		
	loglog(conductivityArray, areaBiofilmAtTargetEfficiencyRatioArray, label=key)

	headerList.append('conductivityArray_' + key)
	headerList.append('areaBiofilmAtTargetEfficiencyRatio_' + key)
	headerList.append('areaBiofilmRelativeToSolarPVArea_' + key)
	vectorList.append(conductivityArray)
	vectorList.append(areaBiofilmAtTargetEfficiencyRatioArray)
	vectorList.append(areaBiofilmRelativeToSolarPVAreaArray)
	
xlabel('Biofilm Conductivity (S cm^{-1})')
ylabel('Biofilm Area (m^2)')
grid()
legend()
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Plot out resistivity versus biofilm area
figure()
for key in keys:
	outputDict = collectedOutputDict[key]
	# Convert from Ω m to Ω cm
	inputResistivityArray = array(outputDict['inputResistivityArray']) * 100 
	areaBiofilmAtTargetEfficiencyRatioArray = outputDict['areaBiofilmAtTargetEfficiencyRatioArray']
	areaBiofilmRelativeToSolarPVAreaArray = outputDict['areaBiofilmRelativeToSolarPVAreaArray']
		
	loglog(inputResistivityArray, areaBiofilmAtTargetEfficiencyRatioArray, label=key)
	
	headerList.append('inputResistivityArray_' + key)
	vectorList.append(inputResistivityArray)
	

xlabel('Biofilm Resistivity (Ohm cm)')
ylabel('Biofilm Area (m^2)')
grid()
legend()
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Plot out conductivity versus biofilm thickness
figure()
for key in keys:
	outputDict = collectedOutputDict[key]
	# Convert to S cm^-1 from S m^-1 
	conductivityArray = array(outputDict['conductivityArray']) / 100
		
	biofilmThicknessAtTargetEfficiencyRatioArray \
	= outputDict['biofilmThicknessAtTargetEfficiencyRatioArray']
	
	loglog(conductivityArray, biofilmThicknessAtTargetEfficiencyRatioArray, label=key)
	
	headerList.append('biofilmThicknessAtTargetEfficiencyRatioArray_' + key)
	vectorList.append(biofilmThicknessAtTargetEfficiencyRatioArray)
	
xlabel('Biofilm Conductivity (S cm^{-1})')
ylabel('Biofilm Thickness (m)')
grid()
legend()
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Plot out resistivity versus biofilm thickness
figure()
for key in keys:
	outputDict = collectedOutputDict[key]
	# Convert from Ω m to Ω cm
	inputResistivityArray = array(outputDict['inputResistivityArray']) * 100
	areaBiofilmAtTargetEfficiencyRatioArray = outputDict['areaBiofilmAtTargetEfficiencyRatioArray']
	areaBiofilmRelativeToSolarPVAreaArray = outputDict['areaBiofilmRelativeToSolarPVAreaArray']
	
	biofilmThicknessAtTargetEfficiencyRatioArray \
	= outputDict['biofilmThicknessAtTargetEfficiencyRatioArray']
	
	loglog(inputResistivityArray, biofilmThicknessAtTargetEfficiencyRatioArray, label=key)
	
	
xlabel('Biofilm Resistivity (Ohm cm)')
ylabel('Biofilm Thickness (m)')
grid()
legend()
# ------------------------------------------------------------------------------------------------ #


oMatrix = generateOutputMatrixWithHeaders(vectorList, headerList, delimeter=',')
writeOutputMatrix(outputFilename, oMatrix)
show()
