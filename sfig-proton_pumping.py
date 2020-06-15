#!/opt/local/bin/python3.7

# ------------------------------------------------------------------------------------------------ #
# sfig-proton_pumping.py
# Calculates rewired carbon fixation efficiency as a function of allowed number of proton
# pumping events
# 
# Farshid Salimijazi and Buz Barstow
# Last updated by Buz Barstow on 2019-11-13
# Last updated by Farshid Salimijazi on 2020-04-17
# Last updated by Farshid Salimijazi on 2020-05-31
# ------------------------------------------------------------------------------------------------ #


from rewiredcarbon.scenario import ImportScenarioTable, CalculateScenarioEfficiencies, \
Plot_Efficiency_Scattergraphs, Export_Efficiency_Scattergraphs, Plot_Efficiency_Scattergraphs_2

from rewiredcarbon.utils import ensure_dir



scenarioTableFileName = 'input/sfig-proton_pumping.csv'
outputFileDirname = 'output/sfig-proton_pumping/'

ensure_dir(outputFileDirname)


scenarioDict = ImportScenarioTable(scenarioTableFileName)

efficienciesDict = CalculateScenarioEfficiencies(scenarioDict, mode='scattergraph')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'effTotalElectricalToFuel')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'protonsPumpedOut')

# Plot_Efficiency_Scattergraphs_2(efficienciesDict, 'protonsPumpedOut','effTotalElectricalToFuel')

outputFilePrefix = 'eff'
Export_Efficiency_Scattergraphs(outputFileDirname, outputFilePrefix, efficienciesDict, \
'effTotalElectricalToFuel', keysToPlot=None)

outputFilePrefix = 'protons'
Export_Efficiency_Scattergraphs(outputFileDirname, outputFilePrefix, efficienciesDict, \
'protonsPumpedOut', keysToPlot=None)