#!/opt/local/bin/python3.7

# ------------------------------------------------------------------------------------------------ #
# fig-h2agitation-B-and-C.py
# Calculates stirring power and tank size for hydrogen-transport by agitation
# 
# Farshid Salimijazi and Buz Barstow
# Last updated by Buz Barstow on 2019-10-26
# ------------------------------------------------------------------------------------------------ #


from rewiredcarbon.scenario import ImportScenarioTable, CalculateScenarioEfficiencies, \
Plot_Efficiency_Scattergraphs, Export_Efficiency_Scattergraphs

from rewiredcarbon.utils import ensure_dir

outputFileDirname = 'output/fig-h2agitation/fig-h2agitation-B-and-C/'
ensure_dir(outputFileDirname)


scenarioTableFileName = 'input/fig-h2agitation-B-and-C.csv'

scenarioDict = ImportScenarioTable(scenarioTableFileName)

efficienciesDict = CalculateScenarioEfficiencies(scenarioDict, mode='scattergraph')



Plot_Efficiency_Scattergraphs(efficienciesDict, 'kLaHydrogen', \
overridePlotScalesInScenarioFile=True, yScale='Logarithmic', xScale='Logarithmic')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'kLaHydrogenGeoInd', \
overridePlotScalesInScenarioFile=True, yScale='Logarithmic', xScale='Logarithmic')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'kLaHydrogenGeoDep', \
overridePlotScalesInScenarioFile=True, yScale='Logarithmic', xScale='Logarithmic')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'kLaHydrogenGeoIndToDepRatio', \
overridePlotScalesInScenarioFile=True, yScale='Linear', xScale='Logarithmic')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'stirPowerDensity', \
overridePlotScalesInScenarioFile=True, yScale='Linear', xScale='Logarithmic')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'cellCultureVolume', \
overridePlotScalesInScenarioFile=True, yScale='Linear', xScale='Logarithmic')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'tankCrossSectionalArea', \
overridePlotScalesInScenarioFile=True, yScale='Linear', xScale='Logarithmic')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'hydrogenCurrent', \
overridePlotScalesInScenarioFile=True, yScale='Linear', xScale='Logarithmic')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'initialStirPowerGuess', \
overridePlotScalesInScenarioFile=True, yScale='Linear', xScale='Logarithmic')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'stirPowerBestGuess', \
overridePlotScalesInScenarioFile=True, yScale='Linear', xScale='Logarithmic')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'stirPower', \
overridePlotScalesInScenarioFile=True, yScale='Linear', xScale='Logarithmic')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'effTotalElectricalToFuel', \
overridePlotScalesInScenarioFile=True, yScale='Linear', xScale='Logarithmic')



Export_Efficiency_Scattergraphs(outputFileDirname, 'effTotalElectricalToFuel', efficienciesDict, \
'effTotalElectricalToFuel', keysToPlot=None, addKeyToHeader=True)

Export_Efficiency_Scattergraphs(outputFileDirname, 'stirPower', efficienciesDict, \
'stirPower', keysToPlot=None, addKeyToHeader=True)

Export_Efficiency_Scattergraphs(outputFileDirname, 'cellCultureVolume', efficienciesDict, \
'cellCultureVolume', keysToPlot=None, addKeyToHeader=True)

Export_Efficiency_Scattergraphs(outputFileDirname, 'kLaHydrogen', efficienciesDict, \
'kLaHydrogen', keysToPlot=None, addKeyToHeader=True)

Export_Efficiency_Scattergraphs(outputFileDirname, 'hydrogenCurrent', efficienciesDict, \
'hydrogenCurrent', keysToPlot=None, addKeyToHeader=True)

Export_Efficiency_Scattergraphs(outputFileDirname, 'tankCrossSectionalArea', efficienciesDict, \
'tankCrossSectionalArea', keysToPlot=None, addKeyToHeader=True)