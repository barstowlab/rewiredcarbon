#!/opt/local/bin/python3.7

# ------------------------------------------------------------------------------------------------ #
# fig-EETscaleup-A.py
# Calculates efficiency of EET based rewired carbon fixation
# 
# Farshid Salimijazi and Buz Barstow
# Last updated by Buz Barstow on 2019-10-28
# ------------------------------------------------------------------------------------------------ #


from rewiredcarbon.scenario import ImportScenarioTable, CalculateScenarioEfficiencies, \
Plot_Efficiency_Scattergraphs, Export_Efficiency_Scattergraphs

from rewiredcarbon.utils import ensure_dir

scenarioTableFileName = 'input/fig-EETscaleup-A.csv'

outputFileDirname = 'output/fig-EETscaleup/fig-EETscaleup-A/'
ensure_dir(outputFileDirname)

outputFilePrefix = 'VariableResistivity'


scenarioDict = ImportScenarioTable(scenarioTableFileName)

efficienciesDict = CalculateScenarioEfficiencies(scenarioDict, mode='scattergraph')



Plot_Efficiency_Scattergraphs(efficienciesDict, 'effTotalElectricalToFuel')



Export_Efficiency_Scattergraphs(outputFileDirname, outputFilePrefix, efficienciesDict, \
'effTotalElectricalToFuel', keysToPlot=None)