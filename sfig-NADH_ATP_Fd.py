#!/opt/local/bin/python3.7

# ------------------------------------------------------------------------------------------------ #
# sfig-NADH_ATP_Fd.py
# Calculates rewired carbon fixation efficiency as a function of NADH, ATP and Fd numbers required
# for fuel synthesis
# 
# Farshid Salimijazi and Buz Barstow
# Last updated by Buz Barstow on 2019-11-13
# ------------------------------------------------------------------------------------------------ #


from rewiredcarbon.scenario import ImportScenarioTable, CalculateScenarioEfficiencies, \
Plot_Efficiency_Scattergraphs, Export_Efficiency_Scattergraphs

from rewiredcarbon.utils import ensure_dir



scenarioTableFileName = 'input/sfig-NADH_ATP_Fd.csv'
outputFileDirname = 'output/sfig-NADH_ATP_Fd/'
outputFilePrefix = 'sfig-NADH_ATP_Fd'

ensure_dir(outputFileDirname)

scenarioDict = ImportScenarioTable(scenarioTableFileName)

efficienciesDict = CalculateScenarioEfficiencies(scenarioDict, mode='scattergraph')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'effTotalElectricalToFuel')

Export_Efficiency_Scattergraphs(outputFileDirname, outputFilePrefix, efficienciesDict, \
'effTotalElectricalToFuel', keysToPlot=None, addKeyToHeader=True)