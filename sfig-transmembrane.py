from rewiredcarbon.scenario import ImportScenarioTable, CalculateScenarioEfficiencies, \
Plot_Efficiency_Scattergraphs, Export_Efficiency_Scattergraphs

from rewiredcarbon.utils import ensure_dir

import os


scenarioTableFileName = 'input/Scenario Tables/Scenario Table Transmembrane Voltage.csv'

scenarioDict = ImportScenarioTable(scenarioTableFileName)

efficienciesDict = CalculateScenarioEfficiencies(scenarioDict, mode='scattergraph')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'effTotalElectricalToFuel')

inputFileBasename = os.path.basename(scenarioTableFileName)
outputFileDirname = os.path.dirname(scenarioTableFileName)
[outputFilePrefix, inputFileExtension] = os.path.splitext(inputFileBasename)


outputFileDirname = \
'output/Supplementary Figure 1 {sfig-transmembrane} - Efficiency vs Transmembrane Voltage'

ensure_dir(outputFileDirname + '/')

outputFilePrefix = 'transmembrane'

Export_Efficiency_Scattergraphs(outputFileDirname, outputFilePrefix, efficienciesDict, \
'effTotalElectricalToFuel', keysToPlot=None)