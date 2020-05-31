#!/opt/local/bin/python3.7

# ------------------------------------------------------------------------------------------------ #
# sfig-proton_pumping_atp.py
# Calculates number of protons that need to be pumped in to regenerate ATP
#
# Farshid Salimijazi and Buz Barstow
# Last updated by Buz Barstow on 2020-05-31
# ------------------------------------------------------------------------------------------------ #


from rewiredcarbon.scenario import ImportScenarioTable, CalculateScenarioEfficiencies, \
Plot_Efficiency_Scattergraphs, Export_Efficiency_Scattergraphs

from rewiredcarbon.utils import ensure_dir



scenarioTableFileName = 'input/sfig-proton_pumping_atp.csv'
outputFileDirname = 'output/sfig-proton_pumping_atp/'

ensure_dir(outputFileDirname)

outputFilePrefix = 'proton_pumping_atp'


scenarioDict = ImportScenarioTable(scenarioTableFileName)

efficienciesDict = CalculateScenarioEfficiencies(scenarioDict, mode='scattergraph')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'numberOfProtonsPumpedInForATP')

Export_Efficiency_Scattergraphs(outputFileDirname, outputFilePrefix, efficienciesDict, \
'numberOfProtonsPumpedInForATP', keysToPlot=None)