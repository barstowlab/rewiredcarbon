#!/opt/local/bin/python3.7

# ------------------------------------------------------------------------------------------------ #
# sfig-vMtr-quinone-B.py
# Calculates how changing the quinone pool redox potential affects efficiency of electricity to 
# butanol conversion
# 
# Farshid Salimijazi and Buz Barstow
# Last updated by Buz Barstow on 2019-11-13
# ------------------------------------------------------------------------------------------------ #


from rewiredcarbon.scenario import ImportScenarioTable, CalculateScenarioEfficiencies, \
Plot_Efficiency_Scattergraphs, Export_Efficiency_Scattergraphs

scenarioTableFileName = 'input/sfig-vMtr-quinone-B.csv'
outputFileDirname = 'output/sfig-vMtr-quinone/'
outputFilePrefix = 'sfig-vMtr-quinone-B'

ensure_dir(outputFileDirname)


scenarioDict = ImportScenarioTable(scenarioTableFileName)

efficienciesDict = CalculateScenarioEfficiencies(scenarioDict, mode='scattergraph')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'effTotalElectricalToFuel')

Export_Efficiency_Scattergraphs(outputFileDirname, 'sfig-vMtr-quinone-B', efficienciesDict, \
'effTotalElectricalToFuel', keysToPlot=None)

