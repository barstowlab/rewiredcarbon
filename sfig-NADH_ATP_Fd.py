from rewiredcarbon.scenario import ImportScenarioTable, CalculateScenarioEfficiencies, \
Plot_Efficiency_Scattergraphs, Export_Efficiency_Scattergraphs

from rewiredcarbon.utils import ensure_dir


scenarioTableFileName = 'input/Scenario Tables/sfig-NADH_ATP_Fd.csv'
outputFileDirname = 'output/sfig-NADH_ATP_Fd'
outputFilePrefix = 'sfig-NADH_ATP_Fd'

ensure_dir(outputFileDirname + '/')

scenarioDict = ImportScenarioTable(scenarioTableFileName)

efficienciesDict = CalculateScenarioEfficiencies(scenarioDict, mode='scattergraph')

Plot_Efficiency_Scattergraphs(efficienciesDict, 'effTotalElectricalToFuel')

Export_Efficiency_Scattergraphs(outputFileDirname, outputFilePrefix, efficienciesDict, \
'effTotalElectricalToFuel', keysToPlot=None, addKeyToHeader=True)