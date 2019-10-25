# ------------------------------------------------------------------------------------------------ #
# fig-co2fixation.py
# Calculates area and thickness of thin film in hydrogen-transport by diffusion scenario
# 
# Farshid Salimijazi and Buz Barstow
# Last updated by Buz Barstow on 2019-10-25
# ------------------------------------------------------------------------------------------------ #


from rewiredcarbon.scenario import ImportScenarioTable, CalculateScenarioEfficiencies, \
Plot_Efficiency_Bargraph, Generate_EfficienciesDict_Keys_Sorted_by_Efficiency, \
Export_Efficiency_Bargraph

from rewiredcarbon.utils import ensure_dir

import os



scenarioTableFileName = 'input/fig-co2fixation.csv'
outputFilename = 'output/fig-co2fixation/fig-co2fixation.csv'
ensure_dir(outputFilename)


scenarioDict = ImportScenarioTable(scenarioTableFileName)

efficienciesDict = CalculateScenarioEfficiencies(scenarioDict)


# keysArray = \
# Generate_EfficienciesDict_Keys_Sorted_by_Efficiency(efficienciesDict, 'effTotalElectricalToFuel')


keysArray = list(efficienciesDict.keys())

Plot_Efficiency_Bargraph(efficienciesDict, 'effTotalElectricalToFuel', \
'effTotalElectricalToFuel_lowerError', 'effTotalElectricalToFuel_upperError', keysToPlot=keysArray)



Export_Efficiency_Bargraph(outputFilename, efficienciesDict, scenarioDict, \
'effTotalElectricalToFuel', 'effTotalElectricalToFuel_lowerError', \
'effTotalElectricalToFuel_upperError', keysToPlot=keysArray)