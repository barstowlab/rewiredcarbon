from rewiredcarbon.scenario import ImportScenarioTable, CalculateScenarioEfficiencies, \
Plot_Efficiency_Scattergraphs, Export_Efficiency_Scattergraphs

scenarioTableFileName = 'input/Scenario Tables/sfig-vMtr-quinone-B.csv'
outputFileDirname = 'output/sfig-vMtr-quinone'
outputFilePrefix = 'sfig-vMtr-quinone-B'



scenarioDict = ImportScenarioTable(scenarioTableFileName)

efficienciesDict = CalculateScenarioEfficiencies(scenarioDict, mode='scattergraph')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'effTotalElectricalToFuel')

Export_Efficiency_Scattergraphs(outputFileDirname, 'sfig-vMtr-quinone-B', efficienciesDict, \
'effTotalElectricalToFuel', keysToPlot=None)
