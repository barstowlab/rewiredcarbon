#!/opt/local/bin/python3.7

# ------------------------------------------------------------------------------------------------ #
# sfig-transmembrane.py
# Calculates rewired carbon fixation efficiency as a function of transmembrane voltage
# 
# Farshid Salimijazi and Buz Barstow
# Last updated by Buz Barstow on 2019-11-13
# ------------------------------------------------------------------------------------------------ #


from rewiredcarbon.scenario import ImportScenarioTable, CalculateScenarioEfficiencies, \
Plot_Efficiency_Scattergraphs, Export_Efficiency_Scattergraphs

from rewiredcarbon.utils import ensure_dir

import os


scenarioTableFileName = 'input/sfig-transmembrane.csv'
outputFileDirname = 'output/sfig-transmembrane/'

ensure_dir(outputFileDirname)

outputFilePrefix = 'transmembrane'


scenarioDict = ImportScenarioTable(scenarioTableFileName)

efficienciesDict = CalculateScenarioEfficiencies(scenarioDict, mode='scattergraph')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'effTotalElectricalToFuel')

Export_Efficiency_Scattergraphs(outputFileDirname, outputFilePrefix, efficienciesDict, \
'effTotalElectricalToFuel', keysToPlot=None)