from rewiredcarbon.scenario import ImportScenarioTable, CalculateScenarioEfficiencies, \
ProcessSingleValueScenario, Plot_Efficiency_Scattergraphs_with_More_than_1_Multi_Value_Variable, \
Combine_Efficiency_Array_with_More_than_1_Multi_Value_Variable, \
Export_Efficiency_Scattergraphs_with_More_than_1_Multi_Value_Variable

from rewiredcarbon.utils import ensure_dir
from copy import deepcopy


scenarioTableFileName = 'input/sfig-vMtr-quinone-A-and-C.csv'
outputFileDirname = 'output/sfig-vMtr-quinone/'
outputFilePrefix = 'sfig-vMtr-quinone-C'

ensure_dir(outputFileDirname)

scenarioDict = ImportScenarioTable(scenarioTableFileName)

multiValueVariable1 = 'vMtr'
multiValueVariable2 = 'voltageCellTwoCathode'
multiValueVariable3 = 'vQuinone'

multiValueVariableDict = {'vMtr':arange(0.2, -0.32, -0.001), \
'voltageCellTwoCathode':arange(0.2, -0.32, -0.001), \
'vQuinone':arange(0.2, -0.32, -0.001)}


scenarioDataArray = []
efficiencyDictArray = []


scenarioKeys = list(scenarioDict.keys())
scenarioData = scenarioDict[scenarioKeys[0]]

multiValueVariableDictKeys = list(multiValueVariableDict.keys())
firstMultiValueVariableArray = multiValueVariableDict[multiValueVariableDictKeys[0]]

i = 0
while i < len(firstMultiValueVariableArray):
	
	scenarioDataSingleValue = deepcopy(scenarioData)
	
	j = 0
	while j < len(multiValueVariableDictKeys):
		multiValueVariableDictKey = multiValueVariableDictKeys[j]
		valueOfMultiValueVariable = multiValueVariableDict[multiValueVariableDictKey][i]
		scenarioDataSingleValue[multiValueVariableDictKey] = valueOfMultiValueVariable
	
		j += 1
		
	scenarioDataArray.append(scenarioDataSingleValue)
	efficiencyDictArray.append(ProcessSingleValueScenario(scenarioDataSingleValue))
	
	i += 1
	
efficienciesDict = {}
 
efficienciesDict[scenarioKeys[0]] = \
Combine_Efficiency_Array_with_More_than_1_Multi_Value_Variable(efficiencyDictArray, \
multiValueVariableDict)


efficienciesDict[scenarioKeys[0]]['independentVariableScale'] \
= scenarioDict[scenarioKeys[0]]['independentVariableScale']

efficienciesDict[scenarioKeys[0]]['dependentVariableScale'] \
= scenarioDict[scenarioKeys[0]]['dependentVariableScale']




Plot_Efficiency_Scattergraphs_with_More_than_1_Multi_Value_Variable(efficienciesDict, \
'effTotalElectricalToFuel')

Export_Efficiency_Scattergraphs_with_More_than_1_Multi_Value_Variable(outputFileDirname, \
outputFilePrefix, efficienciesDict, 'effTotalElectricalToFuel', \
keysToPlot=None, addKeyToHeader=True)
