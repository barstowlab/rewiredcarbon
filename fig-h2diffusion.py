#!/opt/local/bin/python3.7

# ------------------------------------------------------------------------------------------------ #
# fig-h2scaleup-C.py
# Calculates area and thickness of thin film in hydrogen-transport by diffusion scenario
# 
# Farshid Salimijazi and Buz Barstow
# Last updated by Buz Barstow on 2019-10-25
# ------------------------------------------------------------------------------------------------ #


from rewiredcarbon.scenario import ImportScenarioTable, CalculateScenarioEfficiencies, \
Plot_Efficiency_Scattergraphs, Export_Efficiency_Scattergraphs

from rewiredcarbon.utils import ensure_dir

outputFileDirname = 'output/fig-h2diffusion/'
ensure_dir(outputFileDirname)


scenarioTableFileName = 'input/fig-h2diffusion.csv'

scenarioDict = ImportScenarioTable(scenarioTableFileName)

efficienciesDict = CalculateScenarioEfficiencies(scenarioDict, mode='scattergraph')



Plot_Efficiency_Scattergraphs(efficienciesDict, 'hydFilmThickness', \
overridePlotScalesInScenarioFile=True, yScale='Logarithmic', xScale='Logarithmic')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'hydFilmArea', \
overridePlotScalesInScenarioFile=True, yScale='Logarithmic', xScale='Logarithmic')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'effTotalElectricalToFuel', \
overridePlotScalesInScenarioFile=True, yScale='Logarithmic', xScale='Logarithmic')



Export_Efficiency_Scattergraphs(outputFileDirname, 'hydFilmThickness', efficienciesDict, \
'hydFilmThickness', keysToPlot=None, addKeyToHeader=True)


Export_Efficiency_Scattergraphs(outputFileDirname, 'hydFilmArea', efficienciesDict, \
'hydFilmArea', keysToPlot=None, addKeyToHeader=True)

Export_Efficiency_Scattergraphs(outputFileDirname, 'effTotalElectricalToFuel', efficienciesDict, \
'effTotalElectricalToFuel', keysToPlot=None, addKeyToHeader=True)
