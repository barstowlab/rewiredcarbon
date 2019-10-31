from rewiredcarbon.utils import ensure_dir
from rewiredcarbon.vectorOutput import generateOutputMatrixWithHeaders, writeOutputMatrix
from rewiredcarbon.scaleup import \
Calculate_Density_and_Area_for_Target_to_Peak_Efficiency_Ratio_for_Array_of_Input_Powers
from numpy import max
from copy import deepcopy
from rewiredcarbon.scenario import ImportScenarioTable, CalculateScenarioEfficiencies, \
Plot_Efficiency_Scattergraphs, Export_Efficiency_Scattergraphs


# Solar constant in W m^-2
solarConstant = 1000

outputFileDirname = 'output/sfig-h2scaleup-B/'

scenarioTableFileName = 'input/Scenario Tables/sfig-h2scaleup-B.csv'

scenarioDict = ImportScenarioTable(scenarioTableFileName)

efficienciesDict = CalculateScenarioEfficiencies(scenarioDict, mode='scattergraph')


# Calculate out the peak efficiency
peakElectricalToFuelEfficiency = max(efficienciesDict['Reference']['effAvailElectricalToFuel'])

# Remove the reference scenario from consideration
efficienciesDictKeysForCalc = list(efficienciesDict.keys())
efficienciesDictKeysForCalc.remove('Reference')



targetEfficiencyRatios = [0.57,0.855]
collectedOutputDict = {}

# Calculate the system footprint data for different target efficiency ratios
for targetEfficiencyRatio in targetEfficiencyRatios:

	outputDict = \
	Calculate_Density_and_Area_for_Target_to_Peak_Efficiency_Ratio_for_Array_of_Input_Powers(\
	scenarioDict, efficienciesDict, efficienciesDictKeysForCalc, targetEfficiencyRatio, \
	peakElectricalToFuelEfficiency, solarConstant=solarConstant)

	collectedOutputDict[str(targetEfficiencyRatio)] = deepcopy(outputDict)


# Plot out the system footprint calculations and export them to a CSV file
keys = collectedOutputDict.keys()
vectorList_tankAreaRelativeToSolarPVArea = []
headerList_tankAreaRelativeToSolarPVArea = []
outputFilename_tankAreaRelativeToSolarPVArea \
= outputFileDirname + '/' + 'fig-h2scaleup-E-area_relative_to_PV_area.csv'
ensure_dir(outputFilename_tankAreaRelativeToSolarPVArea)

figure()
for key in keys:
	outputDict = collectedOutputDict[key]
	
	totalElectricalPowerArray = outputDict['totalElectricalPowerArray']
	tankAreaRelativeToSolarPVAreaArray = outputDict['tankAreaRelativeToSolarPVAreaArray']
	
	loglog(totalElectricalPowerArray, tankAreaRelativeToSolarPVAreaArray, label=key)
	
	headerList_tankAreaRelativeToSolarPVArea.append('totalElectricalPower_' + key)
	headerList_tankAreaRelativeToSolarPVArea.append('tankAreaRelativeToSolarPVArea_' \
	+ key)
	
	vectorList_tankAreaRelativeToSolarPVArea.append(totalElectricalPowerArray)
	vectorList_tankAreaRelativeToSolarPVArea.append(tankAreaRelativeToSolarPVAreaArray)

	
xlabel('Total Electrical Power (W)')
ylabel('Tank Area Relative To Solar PV Area')
grid()
legend()

oMatrix_tankAreaRelativeToSolarPVArea \
= generateOutputMatrixWithHeaders(vectorList_tankAreaRelativeToSolarPVArea, \
headerList_tankAreaRelativeToSolarPVArea, delimeter=',')

writeOutputMatrix(outputFilename_tankAreaRelativeToSolarPVArea, \
oMatrix_tankAreaRelativeToSolarPVArea)



# Plot out the cell density at target efficient calculations and export them to a CSV file
keys = collectedOutputDict.keys()
vectorList_cellDensityAtTargetEfficiencyRatio = []
headerList_cellDensityAtTargetEfficiencyRatio = []
outputFilename_cellDensityAtTargetEfficiencyRatio \
= outputFileDirname + '/' + 'fig-h2scaleup-E-cell_density.csv'
ensure_dir(outputFilename_cellDensityAtTargetEfficiencyRatio)

figure()
for key in keys:
	outputDict = collectedOutputDict[key]
	totalElectricalPowerArray = outputDict['totalElectricalPowerArray']
	cellDensityAtTargetEfficiencyRatioArray = outputDict['cellDensityAtTargetEfficiencyRatioArray']
	
	loglog(totalElectricalPowerArray, cellDensityAtTargetEfficiencyRatioArray, label=key)
	
	headerList_cellDensityAtTargetEfficiencyRatio.append('totalElectricalPower_' + key)
	headerList_cellDensityAtTargetEfficiencyRatio.append('cellDensityAtTargetEfficiencyRatioArray' \
	+ key)
	
	vectorList_cellDensityAtTargetEfficiencyRatio.append(totalElectricalPowerArray)
	vectorList_cellDensityAtTargetEfficiencyRatio.append(cellDensityAtTargetEfficiencyRatioArray)

xlabel('Total Electrical Power (W)')
ylabel('Cell Density At Target Efficiency Ratio (cells m^-3)')
grid()
legend()


oMatrix_cellDensityAtTargetEfficiencyRatio \
= generateOutputMatrixWithHeaders(vectorList_cellDensityAtTargetEfficiencyRatio, \
headerList_cellDensityAtTargetEfficiencyRatio, delimeter=',')

writeOutputMatrix(outputFilename_cellDensityAtTargetEfficiencyRatio, \
oMatrix_cellDensityAtTargetEfficiencyRatio)


show()


