#!/opt/local/bin/python3.7

# ------------------------------------------------------------------------------------------------ #
# fig-growth_maintenance_energy.py
# Calculates minimum growth and maintenance energy requirements for electromicrobial production


# Farshid Salimijazi and Buz Barstow
# Last updated by Buz Barstow on 2020-07-28
# ------------------------------------------------------------------------------------------------ #


from rewiredcarbon.scenario import ImportScenarioTable, CalculateScenarioEfficiencies
from rewiredcarbon.physicalconstants import avogadroNumber
from rewiredcarbon.utils import ensure_dir

from numpy import average


scenarioTableFileName = 'input/snote-maintenance.csv'
outputFilename = 'output/fig-co2fixation/snote-maintenance.csv'
ensure_dir(outputFilename)


scenarioDict = ImportScenarioTable(scenarioTableFileName)

efficienciesDict = CalculateScenarioEfficiencies(scenarioDict)

key = list(efficienciesDict.keys())[0]

scenarioResults = efficienciesDict[key]

totalCells = efficienciesDict[key]['totalCells']
totalEnergyContentOfFuel = efficienciesDict[key]['totalEnergyContentOfFuel']
energyOfFuelPerCell = totalEnergyContentOfFuel / totalCells
effTotalElectricalToFuel = efficienciesDict[key]['effTotalElectricalToFuel']
inputEnergyPerCell = energyOfFuelPerCell/effTotalElectricalToFuel
totalPowerInput = inputEnergyPerCell*totalCells

averageDryCellWeight = average([148e-15,211e-15,280e-15]) # E. coli dry cell weight from Bionumbers
totalDryCellMass = averageDryCellWeight*totalCells

# Using data from Kliphuis2011a. ATP consumption per cell per second
deltaGADPATPMol = average([57000,28000]) # Bionumbers 101989 and 100775
deltaGADPATP = deltaGADPATPMol/avogadroNumber

# 2.85 mmol atp per gram dry cell weight per hour, from Kliphuis2011a
atpForMaintenanceCellRate = (2.85e-3*avogadroNumber/3600)*averageDryCellWeight 
atpForMaintenancePowerPerCell = atpForMaintenanceCellRate*deltaGADPATP
atpForMaintenance_to_totalPower_perCell = atpForMaintenancePowerPerCell / inputEnergyPerCell

# ATP, and ATP energy needed to make biomass
atpForGramDryCellMass = 109e-3 # 109 mmol per gram dry cell weight, from Kliphuis2011a
atpEnergyForGramCellMass = atpForGramDryCellMass * deltaGADPATPMol # now in joules per gram
atpEnergyForCellMass = atpEnergyForGramCellMass*totalDryCellMass

# Calculate the time needed for input energy to make biomass to become negligible
timeForInputEnergyToEqualGrowthEnergy = atpEnergyForCellMass/totalPowerInput

print("Total chemical power output per cell (Watts): " + str(energyOfFuelPerCell))
print("Total electrical power available per cell (Watts): " + str(inputEnergyPerCell))
print("Total dry cell mass (grams): " + str(totalDryCellMass))
print("Maintenance power for single cell (Watts):" + str(atpForMaintenancePowerPerCell))
print("Ratio of maintenance to total power: " + str(atpForMaintenance_to_totalPower_perCell))
print("Total energy of ATP needed to produce total dry cell mass (Joules): " + str(atpEnergyForCellMass))


