# ------------------------------------------------------------------------------------------------ #
# efficiency.py
# Buz Barstow
# Created: 2017-10-26
# Last Update: 2018-06-27 by Farshid Salimijazi
# Last Update: 2018-07-22 by Buz Barstow
# Last Update: 2018-07-26 by Farshid Salimijazi
# Update: 2019-06-01 by Buz Barstow
# Update: 2019-06-15 by Buz Barstow
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# For symbol definitions see the file SymbolTable.pdf, included with this repository.
# Note that LaTeX macro names, are the same as symbols used here, with some minor mismatches. 
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Variables used in this file

# Input parameters
# vRedoxCellOne - Potential difference across first electrochemical cell in Volts
# vBiasCellOne - Additional external bias voltage across cell one in Volts
# vRedoxCellTwo - Potential difference across second electrochemical cell (bio-cell) in Volts
# vBiasCellTwo - Additional external bias voltage across cell 2 in Volts

# vCellTwo - Potential across second electrochemical cell in Volts

# efficiencyCurrentToFirstCell - Efficiency of conversion of current to first cell
# efficiencyCurrentToSecondCell - Efficiency of conversion of current to second cell
# efficiencyCarbonToSecondCell - Efficiency of conversion of fixed carbon from first to second cell
# Assumed to be 1.0.

# carbonsPerPrimaryFix - Number of carbon atoms per primary fixation


# electronsPerPrimaryFix - Number of electrons needed per carbon to reduce CO2 to 
# primary fixation product in first cell

# electronsNADHforFuel - Number of electrons needed to make NADH to reduce CO2 to fuel

# energyPerFuelMolecule - Energy extracted when fuel molecule is combusted in standard conditions
# in Joules

# carbonsPerFuel - Number of carbon atoms per fuel

# efficiencyElectronsToHyd - Efficiency of conversion of current to H2. Assumed to be 1.0. 

# electricalPower - Scaling power number for calculation in Watts. Efficiency should not change
# if this is changed unless you get close to a very small number (electronsPerFuelInH) number of 
# electrons. This is around 24 to 40 electrons. 

# vRedox - Redox potential difference of electrochemical cell in Volts

# vBias - Additional external bias voltage needed to make H2 in Volts

# vMembrane - Potential across inner membrane in microbe in Volts

# vAcceptor - Potential of terminal electron acceptor in Volts vs. Standard Hydrogen Electrode

# vMtr - Potential of Mtr EET complex in Volts vs. Standard Hydrogen Electrode

# vQuinone - Potential of quinone molecule in Volts vs. Standard Hydrogen Electrode

# NADHforFuel - Number of NADH to reduce CO2 to fuel

# FdForFuel - Number of Ferredoxin to reduce CO2 to fuel

# ATPforFuel - Number of ATP to make fuel from CO2
# energyPerFuelMolecule - Energy extracted when fuel molecule is combusted in standard conditions
# in Joules

# efficiencyElectronsToMtr - Efficiency of conversion of current to reduced Mtr. Assumed to be 1.0.

# efficiencyElectronsMtrToQuinone - Efficiency of conversion of current to reduced Mtr. 
# Assumed to be 1.0.

# numberOfProtonsPumpedInForATP - Number of protons needed to be pumped through ATP synthase to 
# catalyze ADP + Pi -> ATP. When set to -1, we calculate based on the membrane potential and ...

# deltaGADPATP - The free energy difference between ATP and ADP

# cellDensity - Density of planktonic cells in electrochemical cell (usually used for H2-mediated
# cases. Unit is cells per m^3. 
# ------------------------------------------------------------------------------------------------ #


 
 # ------------------------------------------------------------------------------------------------ #
def Electron_Rate_Per_Cell(electronsPerFuel, carbonsFixedPerFuel, co2FixEnzymeRate, co2FixEnzymePerCell):
# Calculates electron demand of a biotic fixation system per microbe cell
# Note, the variable carbonsFixedPerFuel, this is the number of carbons you need to fix to make 
# the fuel molecule, not the number that eventually wind up in it. In some cases, carbon will be
# lost in the production of the fuel (propanol is a good example of this). 
# Also, even though we say co2FixEnzymeRate and co2FixEnzymePerCell, this could all equally apply to any 
# carbon fixing enzyme. 

	fuelRate = (co2FixEnzymeRate * co2FixEnzymePerCell)/carbonsFixedPerFuel
	electronRatePerCell = electronsPerFuel * fuelRate

	return electronRatePerCell
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Energy_Per_Fuel_Molecule(specificEnergyFuelIn, molecularWeightFuelIn):
# Calculates the combustion energy stored in each fuel molecule
# specificEnergyFuelIn - Energy density of fuel in Megajoules per Kg
# molecularWeightFuelIn - Molecular weight of compound in grams per mole


	from rewiredcarbon.physicalconstants import avogadroNumber as N_a
	
	energyPerMol = (specificEnergyFuelIn * 1000) * molecularWeightFuelIn
	energyPerFuelMolecule = energyPerMol/N_a

	return energyPerFuelMolecule
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def ProtonPumping_and_Electron_Requirements_for_EET_to_Fuel(vMembrane, vAcceptor, vMtr, vQuinone, \
deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, numberOfProtonsPumpedInForATP=-1):

	from math import floor, ceil
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	from rewiredcarbon.electrochemicalconstants import vNADH, vFd
	import pdb
	
	# -------------------------------------------------------------------------------------------- #
	# ATP
	# Calculate number of protons needed to be pumped in to generate 1 ATP, and number of electrons
	# needed to be sent downhill to generate proton gradient
	# See Farshid #1 pages 29 and 30

	if numberOfProtonsPumpedInForATP == -1:
		# If the number of protons pumped into make a single ATP is set to -1, we calculate it
		# from the membrane potential
		numberOfProtonsPumpedInForATP = ceil(abs(deltaGADPATP/(e*vMembrane)))
		# print("Number of protons pumped in for ATP: " + str(numberOfProtonsPumpedInForATP))
# 		print("Transmembrane voltage: " + str(vMembrane))
		# Otherwise, you can override this by setting numberOfProtonsPumpedInForATP to a positive
		# value
	
	numberOfProtonsPumpedOutPerElectronDownInEET = floor(abs((vQuinone - vAcceptor) / vMembrane))
	
	try:
		electronsDownATPforFuelInEET = (ATPforFuel * numberOfProtonsPumpedInForATP)\
		/ numberOfProtonsPumpedOutPerElectronDownInEET
	except:
		pdb.set_trace()
	# -------------------------------------------------------------------------------------------- #

	# -------------------------------------------------------------------------------------------- #
	# Reductants: NADH and Ferredoxin
	# See Farshid #1 pages 29 and 30
	# Calculate proton pumping needed to send electrons uphill to NADH and Ferredoxin, and electrons
	# that need to be sent downhill to generate the necessary proton gradient. 
	
	numberOfProtonsPumpedInPerElectronForNADH = ceil(abs((vNADH - vQuinone)/vMembrane))
	electronsNADHforFuel = 2*NADHforFuel
	
	numberOfProtonsPumpedInPerElectronForFd = ceil(abs((vFd - vQuinone)/vMembrane))
	numberOfProtonsPumpedInPerFd = 2*numberOfProtonsPumpedInPerElectronForFd
	electronsFdforFuel = 2*FdForFuel

	electronsDownNADHforFuel = (electronsNADHforFuel * numberOfProtonsPumpedInPerElectronForNADH) \
	/numberOfProtonsPumpedOutPerElectronDownInEET
	
	electronsDownFdForFuel = (electronsFdforFuel * numberOfProtonsPumpedInPerElectronForFd) \
	/numberOfProtonsPumpedOutPerElectronDownInEET
	# -------------------------------------------------------------------------------------------- #

	# -------------------------------------------------------------------------------------------- #
	# Calculate the total number of electrons we need to make one fuel molecule through EET.
	# Defined in Farshid #1 page 30.
	electronsPerFuelInEET = electronsNADHforFuel + electronsFdforFuel \
	+ electronsDownATPforFuelInEET + electronsDownNADHforFuel + electronsDownFdForFuel
	
	# Calculate the number of protons that get pumped out and back in to make the ATP and reductants
	numberOfProtonsPumpedOut = \
	(electronsDownATPforFuelInEET * numberOfProtonsPumpedOutPerElectronDownInEET) \
	+ (electronsDownNADHforFuel * numberOfProtonsPumpedOutPerElectronDownInEET) \
	+ (electronsDownFdForFuel * numberOfProtonsPumpedOutPerElectronDownInEET)
	# -------------------------------------------------------------------------------------------- #


	return numberOfProtonsPumpedOut, electronsPerFuelInEET
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ProtonPumping_and_Electron_Requirements_for_Hydrogen_to_Fuel(vMembrane, vAcceptor, \
deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, numberOfProtonsPumpedInForATP=-1):
	
	from rewiredcarbon.electrochemicalconstants import vH2
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	from math import floor, ceil
	
	# Calculate the number of protons that you get pumped outside the membrane to store energy by 
	# sending one electron downhill.
	# See Farshid #32
	numberOfProtonsPumpedOutPerElectronDownInH = floor(abs((vH2 - vAcceptor) / vMembrane))
	
	# Calculate the number of electrons the bug has to send downhill to the acceptor
	# in order to make all of the ATP it needs for one fuel molecule.

	# Number of ATPs per fuel molecule * number of protons that you pump inside the membrane to 
	# make one ATP / number of protons pumped out by sending one electron downhill to acceptor
	# See Farshid #1 page 30.
	
	if numberOfProtonsPumpedInForATP == -1:
		# If the number of protons pumped into make a single ATP is set to -1, we calculate it
		# from the membrane potential
		numberOfProtonsPumpedInForATP = ceil(abs(deltaGADPATP/(e * vMembrane)))
		# Otherwise, you can override this by setting numberOfProtonsPumpedInForATP to a positive
		# value
	
	electronsPerFuelInH = 2*NADHforFuel + 2*FdForFuel \
	+ ATPforFuel*(numberOfProtonsPumpedInForATP/numberOfProtonsPumpedOutPerElectronDownInH)
	
	electronsDownATPforFuelInH = (ATPforFuel * numberOfProtonsPumpedInForATP)\
	/ numberOfProtonsPumpedOutPerElectronDownInH
	
	numberOfProtonsPumpedOut = \
	electronsDownATPforFuelInH * numberOfProtonsPumpedOutPerElectronDownInH

	
	return numberOfProtonsPumpedOut, electronsPerFuelInH
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Stir_Power_Convenience_Function(powerNumberExtrapolate, agitatorReynoldsNumber, \
cellCultureVolume, fluidDensity, fluidViscosity, tankDiameterToHeightRatio, \
agitatorDiameterToTankDiameterRatio):

	# Convenience function that calls Calculate_Agitator_Power_from_Reynolds_Number_and_Tank_Volume
	# but returns only one variable rather than three.

	from rewiredcarbon.scaleup import Calculate_Agitator_Power_from_Reynolds_Number_and_Tank_Volume
	import pdb
	
	agitatorRotationalVelocity, power, diameter = \
	Calculate_Agitator_Power_from_Reynolds_Number_and_Tank_Volume(powerNumberExtrapolate, \
	agitatorReynoldsNumber, cellCultureVolume, fluidDensity=fluidDensity, \
	fluidViscosity=fluidViscosity, tankDiameterToHeightRatio=tankDiameterToHeightRatio, \
	agitatorDiameterToTankDiameterRatio=agitatorDiameterToTankDiameterRatio)
	
# 	pdb.set_trace()
	
	return power
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def H2_CellCultureVolume_StirringPower_Coupled_Equations(X, totalElectricalPower, vRedox, vBias, \
electronRatePerCell, cellDensity, powerNumberExtrapolate, agitatorReynoldsNumber, fluidDensity, \
fluidViscosity, tankDiameterToHeightRatio, agitatorDiameterToTankDiameterRatio):

	# Function for solving set of coupled equations for scale up of electrosynthesis reactor

	from rewiredcarbon.physicalconstants import elementaryCharge as e
	import pdb
	from numpy import complex
	
	availableElectricalPowerRe, availableElectricalPowerIm, stirPowerRe, stirPowerIm, \
	hydrogenCurrentRe, hydrogenCurrentIm, totalCellsRe, totalCellsIm, \
	cellCultureVolumeRe, cellCultureVolumeIm = X
	
	
	f = [\
	availableElectricalPowerRe - totalElectricalPower + stirPowerRe, \
	hydrogenCurrentRe - (availableElectricalPowerRe / (vRedox + vBias)), \
	totalCellsRe - (hydrogenCurrentRe / (e * electronRatePerCell)), \
	cellCultureVolumeRe - (totalCellsRe / cellDensity), \
	stirPowerRe - Stir_Power_Convenience_Function(powerNumberExtrapolate, agitatorReynoldsNumber, \
	cellCultureVolumeRe, fluidDensity, fluidViscosity, tankDiameterToHeightRatio, \
	agitatorDiameterToTankDiameterRatio), \
	availableElectricalPowerIm, \
	stirPowerIm, \
	hydrogenCurrentIm, \
	totalCellsIm, \
	cellCultureVolumeIm \
	]


	return f
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def Tafel_OverPotential(currentDensity, selfExchangeCurrentDensity, temperature, alpha):
	from rewiredcarbon.physicalconstants import boltzmannConstant as kB
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	from numpy import log
	
	overpotential = ((kB * temperature) / (e * alpha)) \
	* log(currentDensity/selfExchangeCurrentDensity)

	return overpotential
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def H2_FilmArea_Overpotential_Coupled_Equations(X, availableElectricalPower, vRedox, vBias, \
electronRatePerCell, cellDensity, cellFilmThickness, efficiencyElectronsToHyd, cathodeTafelAlpha, \
cathodeSelfExchangeCurrentDensity, cathodeTemperature):
 
    # Function for solving set of coupled equations for scale up of electrosynthesis reactor
    from rewiredcarbon.physicalconstants import elementaryCharge as e
    import pdb
     
    currentCell, cathodeOverpotential, cathodeArea, totalCells, cellCultureVolume  = X
     
    f = [\
    currentCell - (availableElectricalPower / (vRedox + vBias + cathodeOverpotential)), \
    totalCells - ((currentCell * efficiencyElectronsToHyd) / (e * electronRatePerCell)), \
    cellCultureVolume - (cathodeArea * cellFilmThickness), \
    cellCultureVolume - (totalCells / cellDensity ), \
    cathodeOverpotential \
    - Tafel_OverPotential(currentCell/cathodeArea, cathodeSelfExchangeCurrentDensity, \
    cathodeTemperature, cathodeTafelAlpha)
    ]
 
    return f
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def EET_Conductive_Matrix_Area_Overpotential_Coupled_Equations(X, availableElectricalPower, \
vRedox, vBias, electronRatePerCell, biofilmThickness, resistivityBiofilm, efficiencyElectronsToMtr, \
efficiencyElectronsMtrToQuinone, biofilmLayerFillFactor, layerHeight, crossSectionalAreaCell):

	from rewiredcarbon.physicalconstants import elementaryCharge as e
	import pdb

	currentCell, biofilmOverpotential, biofilmArea, totalCells, biofilmVolume  = X
	
	f = [\
	currentCell - (availableElectricalPower / (vRedox + vBias + biofilmOverpotential)), \
	\
	biofilmOverpotential - ((resistivityBiofilm * biofilmThickness *  currentCell) / biofilmArea), \
	\
	totalCells \
	- ((currentCell * efficiencyElectronsMtrToQuinone * efficiencyElectronsToMtr) \
	/ (e * electronRatePerCell) ), \
	\
	biofilmVolume - (biofilmArea * biofilmThickness), \
	\
	biofilmVolume - ( (layerHeight * totalCells * crossSectionalAreaCell) / biofilmLayerFillFactor )
	]
	

	return f
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def EET_Direct_Contact_Area_Overpotential_Coupled_Equations(X, availableElectricalPower, \
vRedox, vBias, electronRatePerCell, contactResistance, contactsPerCell, efficiencyElectronsToMtr, \
efficiencyElectronsMtrToQuinone, monolayerFillFactor, crossSectionalAreaCell):

	from rewiredcarbon.physicalconstants import elementaryCharge as e
	import pdb

	currentCell, contactResistanceOverpotential, areaMonolayer, totalCells  = X
	
	f = [\
	currentCell - (availableElectricalPower / (vRedox + vBias + contactResistanceOverpotential)), \
	contactResistanceOverpotential \
	- ((currentCell * contactResistance ) / (totalCells * contactsPerCell)), \
	totalCells \
	- ((currentCell * efficiencyElectronsMtrToQuinone * efficiencyElectronsToMtr) \
	/ (e * electronRatePerCell) ), \
	areaMonolayer - ( (totalCells * crossSectionalAreaCell) / monolayerFillFactor )
	]

	return f
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Efficiencies_Current_TotalEnergyContentOfFuel_Hydrogen_BioCO2(availableElectricalPower, \
totalElectricalPower, totalInputPower, vRedox, vBias, efficiencyElectronsToHyd, \
energyPerFuelMolecule, electronsPerFuelInH):

	from rewiredcarbon.physicalconstants import elementaryCharge as e
	
	# Calculate the currents in the electrochemical cell. Defined in Buz Medium #1 page 22 and 
	# Farshid #1 page 32.
	currentCell = availableElectricalPower / (vRedox + vBias)
	hydrogenCurrent = efficiencyElectronsToHyd * currentCell

	
	# Calculate the total rate of fuel production (fuel molecules per second)
	# See Farshid #1 page 32 and Buz Medium #1 pages 22-23. 
	totalFuelProductionRateInHyd = hydrogenCurrent / (e * electronsPerFuelInH)
	totalEnergyContentOfFuel = totalFuelProductionRateInHyd * energyPerFuelMolecule
	
	
	# Calculate the EET to fuel efficiencies. Defined in Farshid #1 page 37.
	
	available_Electrical_to_Hyd_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / availableElectricalPower

	total_Electrical_to_Hyd_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / totalElectricalPower

	total_Input_Power_to_Hyd_to_Fuel_Efficiency = totalEnergyContentOfFuel / totalInputPower

	outputDict = {}
	outputDict['effAvailElectricalToFuel'] = available_Electrical_to_Hyd_to_Fuel_Efficiency
	outputDict['effTotalElectricalToFuel'] = total_Electrical_to_Hyd_to_Fuel_Efficiency
	outputDict['effTotalPowerToFuel'] = total_Input_Power_to_Hyd_to_Fuel_Efficiency
	outputDict['totalEnergyContentOfFuel'] = totalEnergyContentOfFuel
	outputDict['hydrogenCurrent'] = hydrogenCurrent
	outputDict['currentCell'] = currentCell

	return outputDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Efficiencies_Current_TotalEnergyContentOfFuel_EET_BioCO2(availableElectricalPower, \
totalElectricalPower, totalInputPower, vRedox, vBias, efficiencyElectronsToMtr, \
efficiencyElectronsMtrToQuinone, energyPerFuelMolecule, electronsPerFuelInEET):

	from rewiredcarbon.physicalconstants import elementaryCharge as e
	
	# Calculate the current in the electrochemical cell.
	currentCell = availableElectricalPower / (vRedox + vBias)

	# This the electrical current (think electrons) carried from the electrochemical cell to the
	# Mtr EET complex, and then to quinone pool. Measured in Amps. 
	# Defined in Farshid #1 page 28. 
	mtrCurrent = efficiencyElectronsToMtr * currentCell
	quinoneCurrent = efficiencyElectronsMtrToQuinone * mtrCurrent

	# Calculate the total rate of fuel production (fuel molecules per second)
	# See Farshid #1 page 30. 
	totalFuelProductionRateInEET = quinoneCurrent / (e * electronsPerFuelInEET)
	
	# Calculate the energy content of the fuel made in one second (in Joules)
	totalEnergyContentOfFuel = totalFuelProductionRateInEET * energyPerFuelMolecule


	# Calculate the EET to fuel efficiencies 
	# Defined in Farshid #1 page 37.

	available_Electrical_to_EET_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / availableElectricalPower

	total_Electrical_to_EET_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / totalElectricalPower
	
	total_Input_Power_to_EET_to_Fuel_Efficiency = totalEnergyContentOfFuel / totalInputPower


	outputDict = {}
	outputDict['effAvailElectricalToFuel'] = available_Electrical_to_EET_to_Fuel_Efficiency
	outputDict['effTotalElectricalToFuel'] = total_Electrical_to_EET_to_Fuel_Efficiency
	outputDict['effTotalPowerToFuel'] = total_Input_Power_to_EET_to_Fuel_Efficiency
	outputDict['totalEnergyContentOfFuel'] = totalEnergyContentOfFuel
	outputDict['mtrCurrent'] = mtrCurrent
	outputDict['quinoneCurrent'] = quinoneCurrent
	outputDict['currentCell'] = currentCell

	return outputDict
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def Efficiencies_Current_TotalEnergyContentOfFuel_Hydrogen_ElectrochemicalCO2(\
availableElectricalPower, totalElectricalPower, totalInputPower, vCellOne, vCellTwo, \
energyPerFuelMolecule, carbonsPerFuel, electronsToAddInHyd, electronsPerPrimaryFix, \
carbonsPerPrimaryFix, efficiencyElectronsToHyd, efficiencyCurrentToFirstCell, \
efficiencyCurrentToSecondCell, efficiencyCarbonToSecondCell):

	from rewiredcarbon.physicalconstants import elementaryCharge as e
	
	# Calculate the currents through first and second electrochemical cells.
    # Defined in Farshid #1 pages 36 and 37.
 	
	if carbonsPerPrimaryFix > 0:
		ICellOneDivICellTwo = \
		(carbonsPerFuel * electronsPerPrimaryFix * efficiencyCurrentToSecondCell)\
		/(electronsToAddInHyd * efficiencyCurrentToFirstCell * carbonsPerPrimaryFix \
		* efficiencyCarbonToSecondCell)
	elif carbonsPerPrimaryFix == 0:
		ICellOneDivICellTwo = 0
	
	ICellTwo = abs(totalElectricalPower)/((ICellOneDivICellTwo * vCellOne) + vCellTwo)
	ICellOne = ICellTwo * ICellOneDivICellTwo
	
	hydrogenCurrent = ICellTwo * efficiencyElectronsToHyd
	
	# Calculate the total rate of fuel production (fuel molecules per second)
	# See Farshid #1 page 30. 
	totalFuelProductionRateInHyd = hydrogenCurrent / (e * electronsToAddInHyd)
	
	# Calculate the energy content of the fuel made in one second (in Joules)
	totalEnergyContentOfFuel = totalFuelProductionRateInHyd * energyPerFuelMolecule
	
	
	# Calculate the EET to fuel efficiencies. Defined in Farshid #1 page 37.
	available_Electrical_to_Hyd_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / availableElectricalPower
	
	total_Electrical_to_Hyd_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / totalElectricalPower
	
	total_Input_Power_to_Hyd_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / totalInputPower
	
	
	outputDict = {}
	outputDict['effAvailElectricalToFuel'] = available_Electrical_to_Hyd_to_Fuel_Efficiency
	outputDict['effTotalElectricalToFuel'] = total_Electrical_to_Hyd_to_Fuel_Efficiency
	outputDict['effTotalPowerToFuel'] = total_Input_Power_to_Hyd_to_Fuel_Efficiency
	outputDict['totalEnergyContentOfFuel'] = totalEnergyContentOfFuel
	outputDict['ICellOne'] = ICellOne
	outputDict['ICellTwo'] = ICellTwo
	outputDict['hydrogenCurrent'] = hydrogenCurrent
	outputDict['currentCell'] = ICellTwo

	return outputDict
# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
def Efficiencies_Current_TotalEnergyContentOfFuel_EET_ElectrochemicalCO2(availableElectricalPower, \
totalElectricalPower, totalInputPower, vCellOne, vCellTwo, energyPerFuelMolecule, carbonsPerFuel, \
electronsToAddInEET, \
electronsPerPrimaryFix, carbonsPerPrimaryFix, efficiencyElectronsToMtr, \
efficiencyElectronsMtrToQuinone, efficiencyCurrentToFirstCell, efficiencyCurrentToSecondCell, \
efficiencyCarbonToSecondCell):

	from rewiredcarbon.physicalconstants import elementaryCharge as e
	
	if 	carbonsPerPrimaryFix > 0:
		ICellOneDivICellTwo = \
		(carbonsPerFuel * electronsPerPrimaryFix * efficiencyCurrentToSecondCell)\
		/(electronsToAddInEET * efficiencyCurrentToFirstCell * carbonsPerPrimaryFix \
		* efficiencyCarbonToSecondCell)
	elif carbonsPerPrimaryFix == 0:
		ICellOneDivICellTwo = 0
	
	ICellTwo = abs(availableElectricalPower)/((ICellOneDivICellTwo * vCellOne) + vCellTwo)
	ICellOne = ICellTwo * ICellOneDivICellTwo
	
	mtrCurrent = ICellTwo * efficiencyElectronsToMtr
	quinoneCurrent = mtrCurrent * efficiencyElectronsMtrToQuinone

	# Calculate the total rate of fuel production (fuel molecules per second)
	# See Farshid #1 page 30. 
	totalFuelProductionRateInEET = quinoneCurrent / (e * electronsToAddInEET)
	
	# Calculate the energy content of the fuel made in one second (in Joules)
	totalEnergyContentOfFuel = totalFuelProductionRateInEET * energyPerFuelMolecule
	
	
	# Calculate the EET to fuel efficiencies. Defined in Farshid #1 page 37.
	available_Electrical_to_EET_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / availableElectricalPower
	
	total_Electrical_to_EET_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / totalElectricalPower
	
	total_Input_Power_to_EET_to_Fuel_Efficiency = \
	totalEnergyContentOfFuel / totalInputPower
	
	
	outputDict = {}
	outputDict['effAvailElectricalToFuel'] = available_Electrical_to_EET_to_Fuel_Efficiency
	outputDict['effTotalElectricalToFuel'] = total_Electrical_to_EET_to_Fuel_Efficiency
	outputDict['effTotalPowerToFuel'] = total_Input_Power_to_EET_to_Fuel_Efficiency
	outputDict['totalEnergyContentOfFuel'] = totalEnergyContentOfFuel
	outputDict['ICellOne'] = ICellOne
	outputDict['ICellTwo'] = ICellTwo
	outputDict['mtrCurrent'] = mtrCurrent
	outputDict['quinoneCurrent'] = quinoneCurrent
	outputDict['currentCell'] = ICellTwo

	return outputDict
# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
def Efficiency_Hydrogen_BioCO2_to_Fuel_No_ScaleUp(vRedox, vBias, vMembrane, vAcceptor, \
NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, efficiencyElectronsToHyd=1.0, 
numberOfProtonsPumpedInForATP=-1, stirPower=0.0, totalElectricalPower=330.0, \
totalInputPower=1000.0):
# Calculate electrical to hydrogen to fuel efficiency 
# This is a really simple version where we don't do any scale up calculations beyond considering
# what the stir power does to efficiency

	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	import sys
	from copy import deepcopy

	# First, calculate proton pumping and electron per fuel molecule requirements. This is 
	# independent of any scale up considerations. 

	numberOfProtonsPumpedOut, electronsPerFuelInH = \
	ProtonPumping_and_Electron_Requirements_for_Hydrogen_to_Fuel(vMembrane, vAcceptor, \
	deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP)

	availableElectricalPower = totalElectricalPower - stirPower

	outputDict = \
	Efficiencies_Current_TotalEnergyContentOfFuel_Hydrogen_BioCO2(availableElectricalPower, \
	totalElectricalPower, totalInputPower, vRedox, vBias, efficiencyElectronsToHyd, \
	energyPerFuelMolecule, electronsPerFuelInH)
	
	efficiencyDict = deepcopy(outputDict)
	
	efficiencyDict['protonsPumpedOut'] = numberOfProtonsPumpedOut
	efficiencyDict['electronsPerFuel'] = electronsPerFuelInH
	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Efficiency_Hydrogen_BioCO2_to_Fuel_Simple_ScaleUp(vRedox, vBias, vMembrane, vAcceptor, \
NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, carbonsFixedPerFuel, \
efficiencyElectronsToHyd=1.0, numberOfProtonsPumpedInForATP=-1, stirPower=0.0, \
totalElectricalPower=330.0, totalInputPower=1000.0, co2FixEnzymeRate=12, co2FixEnzymePerCell=1e6):
# Calculate electrical to hydrogen to fuel efficiency
# This is still a really simple version where we don't do any scale up calculations beyond 
# considering what the stir power does to efficiency and calculating the number of cells needed to 
# process the current

	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	import pdb
	from copy import deepcopy
	

	# First, calculate proton pumping and electron per fuel molecule requirements. This is 
	# independent of any scale up considerations. 

	numberOfProtonsPumpedOut, electronsPerFuelInH = \
	ProtonPumping_and_Electron_Requirements_for_Hydrogen_to_Fuel(vMembrane, vAcceptor, \
	deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP)

	availableElectricalPower = totalElectricalPower - stirPower

	outputDict = \
	Efficiencies_Current_TotalEnergyContentOfFuel_Hydrogen_BioCO2(availableElectricalPower, \
	totalElectricalPower, totalInputPower, vRedox, vBias, efficiencyElectronsToHyd, \
	energyPerFuelMolecule, electronsPerFuelInH)
	
	hydrogenCurrent = outputDict['hydrogenCurrent']
	
	electronRatePerCell = Electron_Rate_Per_Cell(electronsPerFuelInH, carbonsFixedPerFuel, \
	co2FixEnzymeRate, co2FixEnzymePerCell)
	
	totalCells = hydrogenCurrent / (e * electronRatePerCell)
	
	efficiencyDict = deepcopy(outputDict)
	
	efficiencyDict['protonsPumpedOut'] = numberOfProtonsPumpedOut
	efficiencyDict['electronsPerFuel'] = electronsPerFuelInH
	efficiencyDict['stirPower'] = stirPower
	efficiencyDict['totalCells'] = totalCells
	efficiencyDict['electronRatePerCell'] = electronRatePerCell
	efficiencyDict['hydrogenCurrent'] = hydrogenCurrent
	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Efficiency_Hydrogen_BioCO2_to_Fuel_ScaleUp_Stirring_Fixed_CellDensity(vRedox, vBias, vMembrane,\
vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, carbonsFixedPerFuel, \
powerNumberFunction, agitatorReynoldsNumber, \
efficiencyElectronsToHyd=1.0, numberOfProtonsPumpedInForATP=-1, totalElectricalPower=330.0, \
totalInputPower=1000.0, co2FixEnzymeRate=12, co2FixEnzymePerCell=1e6, cellDensity=1e15, \
fluidDensity=997, fluidViscosity=9.1e-4, tankDiameterToHeightRatio=1, \
agitatorDiameterToTankDiameterRatio=0.5):
# Calculate electrical to hydrogen to fuel efficiency
# This is still a more complicated version where we calculate the number of cells needed to 
# process a given current, and self-consistently solve for the power needed to stir that volume of 
# cells using the scipy optimization toolkit. 

	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	import pdb
	from scipy.optimize import root
	from copy import deepcopy
	from numpy import complex
	
	# First, calculate proton pumping and electron per fuel molecule requirements. This is 
	# independent of any scale up considerations. 

	numberOfProtonsPumpedOut, electronsPerFuelInH = \
	ProtonPumping_and_Electron_Requirements_for_Hydrogen_to_Fuel(vMembrane, vAcceptor, \
	deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP)
	
	electronRatePerCell = Electron_Rate_Per_Cell(electronsPerFuelInH, carbonsFixedPerFuel, \
	co2FixEnzymeRate, co2FixEnzymePerCell)

	# First, make a guess about the solution
	availableElectricalPower0 = totalElectricalPower
	
	outputDict0 = \
	Efficiencies_Current_TotalEnergyContentOfFuel_Hydrogen_BioCO2(availableElectricalPower0, \
	totalElectricalPower, totalInputPower, vRedox, vBias, efficiencyElectronsToHyd, \
	energyPerFuelMolecule, electronsPerFuelInH)
	
	hydrogenCurrent0 = outputDict0['hydrogenCurrent']
	
	totalCells0 = hydrogenCurrent0 / (e * electronRatePerCell)
	cellCultureVolume0 = totalCells0 / cellDensity
	
	stirPower0 = Stir_Power_Convenience_Function(powerNumberFunction, agitatorReynoldsNumber, \
	cellCultureVolume0, fluidDensity, fluidViscosity, tankDiameterToHeightRatio, \
	agitatorDiameterToTankDiameterRatio)
	
	initialSolution = [\
	complex(availableElectricalPower0).real, complex(availableElectricalPower0).imag, \
	complex(stirPower0).real, complex(stirPower0).imag, \
	complex(hydrogenCurrent0).real, complex(hydrogenCurrent0).imag, \
	complex(totalCells0).real, complex(totalCells0).imag, \
	complex(cellCultureVolume0).real, complex(cellCultureVolume0).imag ]
	
	# Now, do a self-consistent solution using the starting solution
	
	solution = root(H2_CellCultureVolume_StirringPower_Coupled_Equations, initialSolution, \
	args=(totalElectricalPower, vRedox, vBias, electronRatePerCell, cellDensity, \
	powerNumberFunction, agitatorReynoldsNumber, fluidDensity, fluidViscosity, \
	tankDiameterToHeightRatio, agitatorDiameterToTankDiameterRatio))
	
	availableElectricalPowerRe, availableElectricalPowerIm, stirPowerRe, stirPowerIm, \
	hydrogenCurrentRe, hydrogenCurrentIm, totalCellsRe, totalCellsIm, \
	cellCultureVolumeRe, cellCultureVolumeIm = solution.x
	
	outputDict = \
	Efficiencies_Current_TotalEnergyContentOfFuel_Hydrogen_BioCO2(availableElectricalPowerRe, \
	totalElectricalPower, totalInputPower, vRedox, vBias, efficiencyElectronsToHyd, \
	energyPerFuelMolecule, electronsPerFuelInH)

	efficiencyDict = deepcopy(outputDict)
	efficiencyDict['protonsPumpedOut'] = numberOfProtonsPumpedOut
	efficiencyDict['electronsPerFuel'] = electronsPerFuelInH
	efficiencyDict['stirPower'] = stirPowerRe
	efficiencyDict['totalCells'] = totalCellsRe
	efficiencyDict['electronRatePerCell'] = electronRatePerCell
	efficiencyDict['cellCultureVolume'] = cellCultureVolumeRe
	efficiencyDict['cellDensity'] = cellDensity
	
	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Guess_Initial_Stir_Power(vRedox, vBias, vMembrane, vAcceptor, NADHforFuel, FdForFuel, \
ATPforFuel, energyPerFuelMolecule, carbonsFixedPerFuel, pressureHydrogen, temperature, \
bulkToPeakHydrogenConcentrationRatio, tankDiameterToHeightRatio, Aconstant, Bconstant, \
Cconstant, efficiencyElectronsToHyd, totalElectricalPower, totalInputPower, \
numberOfProtonsPumpedInForATP, co2FixEnzymeRate, co2FixEnzymePerCell, cellDensity, \
stirPowerGuessIncrement=1):

	# This function guesses the stir power needed to stir in the H2 produced by an electrochemical
	# cell operating at a defined voltage. 
	
	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	from rewiredcarbon.scaleup import Calculate_Agitator_Power_from_Hydrogen_Current
	import pdb
	from scipy.optimize import root
	from copy import deepcopy
	from numpy import complex, mean, std, min, max, log, exp, array
	
	# First, calculate proton pumping and electron per fuel molecule requirements. This is 
	# independent of any scale up considerations. 

	numberOfProtonsPumpedOut, electronsPerFuelInH = \
	ProtonPumping_and_Electron_Requirements_for_Hydrogen_to_Fuel(vMembrane, vAcceptor, \
	deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP)
	
	electronRatePerCell = Electron_Rate_Per_Cell(electronsPerFuelInH, carbonsFixedPerFuel, \
	co2FixEnzymeRate, co2FixEnzymePerCell)

	# First, make a guess about the solution
	
	stirPowerGuess = 0
	stirPowerGuessArray = []
	stirPowerCalculationArray = []
	
	while stirPowerGuess < totalElectricalPower:
		
		availableElectricalPower0 = totalElectricalPower - stirPowerGuess
		
		stirPowerGuessArray.append(stirPowerGuess)
	
		outputDict0 = \
		Efficiencies_Current_TotalEnergyContentOfFuel_Hydrogen_BioCO2(availableElectricalPower0, \
		totalElectricalPower, totalInputPower, vRedox, vBias, efficiencyElectronsToHyd, \
		energyPerFuelMolecule, electronsPerFuelInH)
	
		hydrogenCurrent0 = outputDict0['hydrogenCurrent']
	
		totalCells0 = hydrogenCurrent0 / (e * electronRatePerCell)
		cellCultureVolume0 = totalCells0 / cellDensity
	
		agitatorPowerOutputDict0 = \
		Calculate_Agitator_Power_from_Hydrogen_Current(hydrogenCurrent0, cellCultureVolume0, \
		pressureHydrogen=pressureHydrogen, temperature=temperature, \
		bulkToPeakHydrogenConcentrationRatio=bulkToPeakHydrogenConcentrationRatio, \
		Aconstant=Aconstant, Bconstant=Bconstant, Cconstant=Cconstant, \
		tankDiameterToHeightRatio=tankDiameterToHeightRatio)
	
		stirPower0 = agitatorPowerOutputDict0['stirPower']
		
		stirPowerCalculationArray.append(stirPower0)
		
		stirPowerGuess += stirPowerGuessIncrement
	
	stirPowerGuessCalcDiffArray = array(stirPowerGuessArray) - array(stirPowerCalculationArray)
	stirPowerGuessCalcAbsDiffArray = \
	list(abs(array(stirPowerGuessArray) - array(stirPowerCalculationArray)))
	
	stirPowerBestGuessIndex = \
	stirPowerGuessCalcAbsDiffArray.index(min(stirPowerGuessCalcAbsDiffArray))
	
	stirPowerBestGuess = stirPowerGuessArray[stirPowerBestGuessIndex]
	
	stirPowerGuessDict = {}
	stirPowerGuessDict['stirPowerGuessArray'] = stirPowerGuessArray
	stirPowerGuessDict['stirPowerCalculationArray'] = stirPowerCalculationArray
	stirPowerGuessDict['stirPowerGuessCalcDiffArray'] = stirPowerGuessCalcDiffArray
	stirPowerGuessDict['stirPowerGuessCalcDiffArray'] = stirPowerGuessCalcDiffArray
	stirPowerGuessDict['stirPowerBestGuess'] = stirPowerBestGuess
		
	return stirPowerGuessDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Stir_Power_Guess_to_Calc_Difference(stirPowerGuess, vRedox, vBias, vMembrane, vAcceptor, \
NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, carbonsFixedPerFuel, pressureHydrogen, \
temperature, bulkToPeakHydrogenConcentrationRatio, tankDiameterToHeightRatio, Aconstant, \
Bconstant, Cconstant, efficiencyElectronsToHyd, totalElectricalPower, totalInputPower, \
numberOfProtonsPumpedInForATP, co2FixEnzymeRate, co2FixEnzymePerCell, cellDensity):

	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	from rewiredcarbon.scaleup import Calculate_Agitator_Power_from_Hydrogen_Current
	
	# First, calculate proton pumping and electron per fuel molecule requirements. This is 
	# independent of any scale up considerations. 

	numberOfProtonsPumpedOut, electronsPerFuelInH = \
	ProtonPumping_and_Electron_Requirements_for_Hydrogen_to_Fuel(vMembrane, vAcceptor, \
	deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP)
	
	electronRatePerCell = Electron_Rate_Per_Cell(electronsPerFuelInH, carbonsFixedPerFuel, \
	co2FixEnzymeRate, co2FixEnzymePerCell)

	# Make a guess about the solution
	
	availableElectricalPower0 = totalElectricalPower - stirPowerGuess
	
	outputDict0 = \
	Efficiencies_Current_TotalEnergyContentOfFuel_Hydrogen_BioCO2(availableElectricalPower0, \
	totalElectricalPower, totalInputPower, vRedox, vBias, efficiencyElectronsToHyd, \
	energyPerFuelMolecule, electronsPerFuelInH)

	hydrogenCurrent0 = outputDict0['hydrogenCurrent']

	totalCells0 = hydrogenCurrent0 / (e * electronRatePerCell)
	cellCultureVolume0 = totalCells0 / cellDensity


	agitatorPowerOutputDict = \
	Calculate_Agitator_Power_from_Hydrogen_Current(hydrogenCurrent0, cellCultureVolume0, \
	pressureHydrogen=pressureHydrogen, temperature=temperature, \
	bulkToPeakHydrogenConcentrationRatio=bulkToPeakHydrogenConcentrationRatio, \
	Aconstant=Aconstant, Bconstant=Bconstant, Cconstant=Cconstant, \
	tankDiameterToHeightRatio=tankDiameterToHeightRatio)

	stirPowerCalc = agitatorPowerOutputDict['stirPower']
		
	stirPowerGuess_to_Calc_Difference = stirPowerGuess - stirPowerCalc
			
	return stirPowerGuess_to_Calc_Difference
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def Efficiency_Hydrogen_BioCO2_to_Fuel_ScaleUp_Stirring_Fixed_CellDensity_Grenville(vRedox, vBias,\
vMembrane, vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, \
carbonsFixedPerFuel, \
pressureHydrogen=5066.25, temperature=293, bulkToPeakHydrogenConcentrationRatio=0.5, \
tankDiameterToHeightRatio=1, Aconstant=0.026, Bconstant=0.4, Cconstant=0.5, \
efficiencyElectronsToHyd=1.0, totalElectricalPower=330.0, totalInputPower=1000.0, \
numberOfProtonsPumpedInForATP=-1, co2FixEnzymeRate=12, co2FixEnzymePerCell=1e6, cellDensity=1e15):

	# Calculates the efficiency of an H2-mediated system that uses stirring for H2 transport.
	# Calculates the solution to a coupled set of equations that relate stir power, remaining power
	# for H2 production and cell culture volume (determined by cell density, and total number of 
	# cells required for H2 produced by available power).
	
	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	from rewiredcarbon.scaleup import Calculate_Gas_Mixing_Parameters
	from scipy.optimize import minimize
	from copy import deepcopy
	import pdb
	
	# Calculate internal cellular parameters that are independent of available power considerations
	numberOfProtonsPumpedOut, electronsPerFuelInH = \
	ProtonPumping_and_Electron_Requirements_for_Hydrogen_to_Fuel(vMembrane, vAcceptor, \
	deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP)
	
	electronRatePerCell = Electron_Rate_Per_Cell(electronsPerFuelInH, carbonsFixedPerFuel, \
	co2FixEnzymeRate, co2FixEnzymePerCell)

	# Make an initial guess for stirring power by making a guess for stir power, calculating
	# how much H2 is produced with the remaining power, and calculating how much power 
	# is needed to stir that H2 into solution. When the difference is smallest, this is our
	# best first guess. 
	
	initialStirPowerGuessDict = Guess_Initial_Stir_Power(vRedox, vBias, vMembrane, vAcceptor, \
	NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, carbonsFixedPerFuel, \
	pressureHydrogen, temperature, bulkToPeakHydrogenConcentrationRatio, tankDiameterToHeightRatio,\
	Aconstant, Bconstant, Cconstant, efficiencyElectronsToHyd, totalElectricalPower, \
	totalInputPower, numberOfProtonsPumpedInForATP, co2FixEnzymeRate, co2FixEnzymePerCell, \
	cellDensity)
	
	initialStirPowerGuess = initialStirPowerGuessDict['stirPowerBestGuess']
	
	# Take the best guess and refine it with the optimization algorithm. 
	
	args0 = (vRedox, vBias, vMembrane, vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, \
	energyPerFuelMolecule, carbonsFixedPerFuel, pressureHydrogen, temperature, \
	bulkToPeakHydrogenConcentrationRatio, tankDiameterToHeightRatio, Aconstant, Bconstant, \
	Cconstant, efficiencyElectronsToHyd, totalElectricalPower, totalInputPower, \
	numberOfProtonsPumpedInForATP, co2FixEnzymeRate, co2FixEnzymePerCell, cellDensity)
	
	constraints = (\
	{'type': 'eq', "fun": lambda x: \
	Stir_Power_Guess_to_Calc_Difference(x, \
	vRedox, vBias, vMembrane, vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, \
	energyPerFuelMolecule, carbonsFixedPerFuel, pressureHydrogen, temperature, \
	bulkToPeakHydrogenConcentrationRatio, tankDiameterToHeightRatio, Aconstant, Bconstant, \
	Cconstant, efficiencyElectronsToHyd, totalElectricalPower, totalInputPower, \
	numberOfProtonsPumpedInForATP, co2FixEnzymeRate, co2FixEnzymePerCell, cellDensity) })
		
	solution = minimize(Stir_Power_Guess_to_Calc_Difference, initialStirPowerGuess, args=args0, \
	constraints=constraints)
	
	stirPowerBestGuess = solution.x[0]
		
	availableElectricalPower = totalElectricalPower - stirPowerBestGuess
	
	# Use the refined estimated of the available electrical power to calculate electricity to H2 
	# to biofuel efficiency
	
	outputDict = \
	Efficiencies_Current_TotalEnergyContentOfFuel_Hydrogen_BioCO2(availableElectricalPower, \
	totalElectricalPower, totalInputPower, vRedox, vBias, efficiencyElectronsToHyd, \
	energyPerFuelMolecule, electronsPerFuelInH)
	
	hydrogenCurrent = outputDict['hydrogenCurrent']
	
	totalCells = hydrogenCurrent / (e * electronRatePerCell)
	cellCultureVolume = totalCells / cellDensity
	
	# Use the volume and current to back-calculate gas mixing parameters
	
	gasMixingDict = Calculate_Gas_Mixing_Parameters(hydrogenCurrent, cellCultureVolume, \
	stirPowerBestGuess, \
	pressureHydrogen=pressureHydrogen, temperature=temperature, \
	bulkToPeakHydrogenConcentrationRatio=bulkToPeakHydrogenConcentrationRatio, \
	Aconstant=Aconstant, Bconstant=Bconstant, Cconstant=Cconstant, \
	tankDiameterToHeightRatio=tankDiameterToHeightRatio)
	
		
	efficiencyDict = deepcopy(outputDict)
	efficiencyDict['stirPowerBestGuess'] = stirPowerBestGuess
	efficiencyDict['stirPower'] = stirPowerBestGuess
	efficiencyDict['initialStirPowerGuess'] = initialStirPowerGuess
	efficiencyDict['availableElectricalPower'] = availableElectricalPower	
	efficiencyDict['protonsPumpedOut'] = numberOfProtonsPumpedOut
	efficiencyDict['electronsPerFuel'] = electronsPerFuelInH
	efficiencyDict['electronRatePerCell'] = electronRatePerCell
	efficiencyDict['cellDensity'] = cellDensity
	efficiencyDict['totalCells'] = totalCells
	efficiencyDict['cellCultureVolume'] = cellCultureVolume
	
	efficiencyDict['kLaHydrogen'] = gasMixingDict['kLaHydrogen']
	efficiencyDict['kLaHydrogenGeoInd'] = gasMixingDict['kLaHydrogenGeoInd']
	efficiencyDict['kLaHydrogenGeoDep'] = gasMixingDict['kLaHydrogenGeoDep']
	efficiencyDict['kLaHydrogenGeoIndToDepRatio'] = gasMixingDict['kLaHydrogenGeoIndToDepRatio']
	
	efficiencyDict['superficialGasVelocity'] = gasMixingDict['superficialGasVelocity']
	efficiencyDict['stirPowerDensity'] = gasMixingDict['stirPowerDensity']
	
	efficiencyDict['volumeTransferRateHydrogen'] = gasMixingDict['volumeTransferRateHydrogen']
	efficiencyDict['tankCrossSectionalArea'] = gasMixingDict['tankCrossSectionalArea']
	efficiencyDict['massTransferRateHydrogen'] = gasMixingDict['massTransferRateHydrogen']
	efficiencyDict['peakConcentrationHydrogen'] = gasMixingDict['peakConcentrationHydrogen']
	efficiencyDict['bulkConcentrationHydrogen'] = gasMixingDict['bulkConcentrationHydrogen']
	
# 	pdb.set_trace()
	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Efficiency_Hydrogen_BioCO2_to_Fuel_ScaleUp_Diffusion_Fixed_CellDensity(vRedox, vBias, \
vMembrane, vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, \
carbonsFixedPerFuel, \
efficiencyElectronsToHyd=1.0, numberOfProtonsPumpedInForATP=-1, totalElectricalPower=330.0, \
totalInputPower=1000.0, co2FixEnzymeRate=12, co2FixEnzymePerCell=1e6, cellDensity=1e15, \
cathodeSelfExchangeCurrentDensity=10, cathodeTafelAlpha=1, cathodeTemperature=293, \
pressureHydrogen=5000, solubilityHydrogen=1.282e5, diffusionHydrogen=4.5e-9):
# Calculate electrical to hydrogen to fuel efficiency
# In this case, we're going to spread the cells out over a very thin film with a large
# area, and calculate the current density for an electrode with the same area, and from 
# this calculate the overpotential from the Tafel equation. This will be self-consistently solved 
# using the scipy optimization toolkit.  
# Default value for solubilityHydrogen is in mol Pa^-1 m^-3
# Default value for diffusionHydrogen is in m^2 s^-1
# Default value for pressureHydrogen is in Pa. Default value is 5% H2 at atmospheric pressure 
# (0.1 MPa)

	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	from rewiredcarbon.scaleup import Calculate_Film_Thickness_for_Diffusion
	import pdb
	import sys
	from scipy.optimize import root
	from copy import deepcopy
	
	# First, calculate proton pumping and electron per fuel molecule requirements. This is 
	# independent of any scale up considerations. 

	numberOfProtonsPumpedOut, electronsPerFuelInH = \
	ProtonPumping_and_Electron_Requirements_for_Hydrogen_to_Fuel(vMembrane, vAcceptor, \
	deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP)
	
	electronRatePerCell = Electron_Rate_Per_Cell(electronsPerFuelInH, carbonsFixedPerFuel, \
	co2FixEnzymeRate, co2FixEnzymePerCell)
	
	# First, make a guess about the solution	
	outputDict0 = \
	Efficiencies_Current_TotalEnergyContentOfFuel_Hydrogen_BioCO2(totalElectricalPower, \
	totalElectricalPower, totalInputPower, vRedox, vBias, efficiencyElectronsToHyd, \
	energyPerFuelMolecule, electronsPerFuelInH)
	
	totalCells0 = outputDict0['hydrogenCurrent'] / (e * electronRatePerCell)
	cellCultureVolume0 = totalCells0 / cellDensity
	
	currentCell0 = outputDict0['hydrogenCurrent'] / efficiencyElectronsToHyd
				
	cellFilmThickness = Calculate_Film_Thickness_for_Diffusion(pressureHydrogen, \
	diffusionHydrogen, solubilityHydrogen, electronRatePerCell, cellDensity)
	
	cathodeArea0 = cellCultureVolume0 / cellFilmThickness
	
	cathodeOverpotential0 = Tafel_OverPotential(currentCell0/cathodeArea0, \
	cathodeSelfExchangeCurrentDensity, cathodeTemperature, cathodeTafelAlpha)
	
	initialSolution = [currentCell0, cathodeOverpotential0, cathodeArea0, totalCells0, \
	cellCultureVolume0]
	
# 	pdb.set_trace()
	
	solution = root(H2_FilmArea_Overpotential_Coupled_Equations, initialSolution, \
	args=(totalElectricalPower, vRedox, vBias, electronRatePerCell, cellDensity, \
	cellFilmThickness, efficiencyElectronsToHyd, cathodeTafelAlpha, \
	cathodeSelfExchangeCurrentDensity, cathodeTemperature))

	currentCell, cathodeOverpotential, cathodeArea, totalCells, cellCultureVolume = solution.x
		
	vBiasTotal = vBias + cathodeOverpotential
		
	outputDict = \
	Efficiencies_Current_TotalEnergyContentOfFuel_Hydrogen_BioCO2(totalElectricalPower, \
	totalElectricalPower, totalInputPower, vRedox, vBiasTotal, efficiencyElectronsToHyd, \
	energyPerFuelMolecule, electronsPerFuelInH)
	
	efficiencyDict = deepcopy(outputDict)
	efficiencyDict['protonsPumpedOut'] = numberOfProtonsPumpedOut
	efficiencyDict['electronsPerFuel'] = electronsPerFuelInH
	efficiencyDict['totalCells'] = totalCells
	efficiencyDict['electronRatePerCell'] = electronRatePerCell
	efficiencyDict['cellCultureVolume'] = cellCultureVolume
	efficiencyDict['cellDensity'] = cellDensity
	efficiencyDict['cathodeOverpotential'] = cathodeOverpotential
	efficiencyDict['cellFilmThickness'] = cellFilmThickness
	efficiencyDict['cathodeArea'] = cathodeArea
	efficiencyDict['cellCurrentDensity'] = efficiencyDict['currentCell'] / cathodeArea
	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Efficiency_Hydrogen_BioCO2_to_Fuel_ScaleUp_Diffusion_Fixed_CellDensity_Fixed_Overpotential(\
vRedox, vBias, vMembrane, vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, \
carbonsFixedPerFuel, \
efficiencyElectronsToHyd=1.0, numberOfProtonsPumpedInForATP=-1, totalElectricalPower=330.0, \
totalInputPower=1000.0, co2FixEnzymeRate=12, co2FixEnzymePerCell=1e6, cellDensity=1e15, \
pressureHydrogen=5000, solubilityHydrogen=1.282e5, diffusionHydrogen=4.5e-9):
# Calculate electrical to hydrogen to fuel efficiency
# In this case, we're going to spread the cells out over a very thin film with a large
# area. 
# There should be no variation in efficiency with cell density, but a change in area and thickness
# of the cell film 
# Default value for solubilityHydrogen is in mol Pa^-1 m^-3
# Default value for diffusionHydrogen is in m^2 s^-1
# Default value for pressureHydrogen is in Pa. Default value is 5% H2 at atmospheric pressure 
# (0.1 MPa)

	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	from rewiredcarbon.scaleup import Calculate_Film_Thickness_for_Diffusion
	import pdb
	import sys
	from scipy.optimize import root
	from copy import deepcopy
	
	# First, calculate proton pumping and electron per fuel molecule requirements. This is 
	# independent of any scale up considerations. 

	numberOfProtonsPumpedOut, electronsPerFuelInH = \
	ProtonPumping_and_Electron_Requirements_for_Hydrogen_to_Fuel(vMembrane, vAcceptor, \
	deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP)
	
	electronRatePerCell = Electron_Rate_Per_Cell(electronsPerFuelInH, carbonsFixedPerFuel, \
	co2FixEnzymeRate, co2FixEnzymePerCell)
	
	
	outputDict = \
	Efficiencies_Current_TotalEnergyContentOfFuel_Hydrogen_BioCO2(totalElectricalPower, \
	totalElectricalPower, totalInputPower, vRedox, vBias, efficiencyElectronsToHyd, \
	energyPerFuelMolecule, electronsPerFuelInH)
	
	totalCells = outputDict['hydrogenCurrent'] / (e * electronRatePerCell)
	cellCultureVolume = totalCells / cellDensity
	
	
	hydFilmThickness = Calculate_Film_Thickness_for_Diffusion(pressureHydrogen, \
	diffusionHydrogen, solubilityHydrogen, electronRatePerCell, cellDensity)
	
	hydFilmArea = cellCultureVolume / hydFilmThickness
		
	efficiencyDict = deepcopy(outputDict)
	efficiencyDict['protonsPumpedOut'] = numberOfProtonsPumpedOut
	efficiencyDict['electronsPerFuel'] = electronsPerFuelInH
	efficiencyDict['totalCells'] = totalCells
	efficiencyDict['electronRatePerCell'] = electronRatePerCell
	efficiencyDict['cellCultureVolume'] = cellCultureVolume
	efficiencyDict['hydFilmThickness'] = hydFilmThickness
	efficiencyDict['hydFilmArea'] = hydFilmArea
	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def Efficiency_EET_BioCO2_to_Fuel_No_ScaleUp(vRedox, vBias, vMembrane, vAcceptor, vMtr, vQuinone, \
NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, efficiencyElectronsToMtr=1.0, \
efficiencyElectronsMtrToQuinone=1.0, numberOfProtonsPumpedInForATP=-1, totalElectricalPower=330.0, \
totalInputPower=1000.0):

	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	import pdb
	from copy import deepcopy
	
	if vMtr > vQuinone:
		print('Cannot transfer electron down due to cytochrome Mtr being more positive than' \
		+ ' menaquinone')
		sys.exit()

	availableElectricalPower = totalElectricalPower
	
	numberOfProtonsPumpedOut, electronsPerFuelInEET = \
	ProtonPumping_and_Electron_Requirements_for_EET_to_Fuel(vMembrane, vAcceptor, vMtr, vQuinone, \
	deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP)
	
	outputDict = \
	Efficiencies_Current_TotalEnergyContentOfFuel_EET_BioCO2(availableElectricalPower, \
	totalElectricalPower, totalInputPower, vRedox, vBias, efficiencyElectronsToMtr, \
	efficiencyElectronsMtrToQuinone, energyPerFuelMolecule, electronsPerFuelInEET)
	
	efficiencyDict = deepcopy(outputDict)

	efficiencyDict['protonsPumpedOut'] = numberOfProtonsPumpedOut
	efficiencyDict['electronsPerFuel'] = electronsPerFuelInEET
	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Efficiency_EET_BioCO2_to_Fuel_Simple_ScaleUp(vRedox, vBias, vMembrane, vAcceptor, vMtr, \
vQuinone, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, carbonsFixedPerFuel, \
efficiencyElectronsToMtr=1.0, \
efficiencyElectronsMtrToQuinone=1.0, numberOfProtonsPumpedInForATP=-1, totalElectricalPower=330.0, \
totalInputPower=1000.0, co2FixEnzymeRate=12, co2FixEnzymePerCell=1e6):

	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	import pdb
	from copy import deepcopy
	
	if vMtr > vQuinone:
		print('Cannot transfer electron down due to cytochrome Mtr being more positive than' \
		+ ' menaquinone')
		sys.exit()

	availableElectricalPower = totalElectricalPower
	
	numberOfProtonsPumpedOut, electronsPerFuelInEET = \
	ProtonPumping_and_Electron_Requirements_for_EET_to_Fuel(vMembrane, vAcceptor, vMtr, vQuinone, \
	deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP)
	
	outputDict = \
	Efficiencies_Current_TotalEnergyContentOfFuel_EET_BioCO2(availableElectricalPower, \
	totalElectricalPower, totalInputPower, vRedox, vBias, efficiencyElectronsToMtr, \
	efficiencyElectronsMtrToQuinone, energyPerFuelMolecule, electronsPerFuelInEET)

	quinoneCurrent = outputDict['quinoneCurrent']
	
	electronRatePerCell = Electron_Rate_Per_Cell(electronsPerFuelInEET, carbonsFixedPerFuel, \
	co2FixEnzymeRate, co2FixEnzymePerCell)
	
	totalCells = quinoneCurrent / (e * electronRatePerCell)
	
	efficiencyDict = deepcopy(outputDict)
	efficiencyDict['protonsPumpedOut'] = numberOfProtonsPumpedOut
	efficiencyDict['electronsPerFuel'] = electronsPerFuelInEET
	efficiencyDict['electronRatePerCell'] = electronRatePerCell
	efficiencyDict['totalCells'] = totalCells
	efficiencyDict['quinoneCurrent'] = quinoneCurrent

	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
def Efficiency_EET_BioCO2_to_Fuel_Conductive_Matrix_Fixed_Height(vRedox, vBias, vMembrane, \
vAcceptor, vMtr, vQuinone, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, \
carbonsFixedPerFuel, approximateBiofilmThickness, resistivityBiofilm, cellLength, cellDiameter, \
biofilmLayerHeight, \
biofilmLayerFillFactor=0.25, efficiencyElectronsToMtr=1.0, efficiencyElectronsMtrToQuinone=1.0, \
numberOfProtonsPumpedInForATP=-1, totalElectricalPower=330.0, totalInputPower=1000.0, \
co2FixEnzymeRate=12, co2FixEnzymePerCell=1e6):

	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	from rewiredcarbon.scaleup import Calculate_Precise_Conductive_Matrix_Biofilm_Thickness
	import pdb
	from scipy.optimize import root
	from numpy import pi
	from copy import deepcopy
	
	if vMtr > vQuinone:
		print('Cannot transfer electron down due to cytochrome Mtr being more positive than' \
		+ ' menaquinone')
		sys.exit()

	availableElectricalPower = totalElectricalPower
	
	numberOfProtonsPumpedOut, electronsPerFuelInEET = \
	ProtonPumping_and_Electron_Requirements_for_EET_to_Fuel(vMembrane, vAcceptor, vMtr, vQuinone, \
	deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP)
	
	electronRatePerCell = Electron_Rate_Per_Cell(electronsPerFuelInEET, carbonsFixedPerFuel, \
	co2FixEnzymeRate, co2FixEnzymePerCell)
	
	# Calculate the cell cross sectional area. Assume a sphero-spherical model and calculate the 
	# projected cross section. 
	crossSectionalAreaCell = (cellLength - cellDiameter)*cellDiameter + pi*((cellDiameter/2)**2)
	
	# First, make a guess about the solution
	outputDict0 = \
	Efficiencies_Current_TotalEnergyContentOfFuel_EET_BioCO2(availableElectricalPower, \
	totalElectricalPower, totalInputPower, vRedox, vBias, efficiencyElectronsToMtr, \
	efficiencyElectronsMtrToQuinone, energyPerFuelMolecule, electronsPerFuelInEET)

	totalCells0 = outputDict0['quinoneCurrent'] / (e * electronRatePerCell)
	currentCell0 = outputDict0['currentCell']
	
	preciseBiofilmThickness, numberCellLayers = \
	Calculate_Precise_Conductive_Matrix_Biofilm_Thickness(approximateBiofilmThickness, \
	biofilmLayerHeight)
	
	volumeBiofilm0 = biofilmLayerHeight * totalCells0 * crossSectionalAreaCell \
	/ biofilmLayerFillFactor
	
	areaBiofilm0 = volumeBiofilm0 / preciseBiofilmThickness
	
	biofilmOverpotential0 = \
	(currentCell0 * resistivityBiofilm * preciseBiofilmThickness) / areaBiofilm0
	
	initialSolution = [currentCell0, biofilmOverpotential0, areaBiofilm0, totalCells0, \
	volumeBiofilm0]
	
	solution = root(EET_Conductive_Matrix_Area_Overpotential_Coupled_Equations, initialSolution, \
	args=(availableElectricalPower, vRedox, vBias, electronRatePerCell, preciseBiofilmThickness, \
	resistivityBiofilm, efficiencyElectronsToMtr, efficiencyElectronsMtrToQuinone, \
	biofilmLayerFillFactor, biofilmLayerHeight, crossSectionalAreaCell))
	
	currentCell, biofilmOverpotential, areaBiofilm, totalCells, volumeBiofilm = solution.x
	
	vBiasTotal = vBias + biofilmOverpotential
	
	outputDict = \
	Efficiencies_Current_TotalEnergyContentOfFuel_EET_BioCO2(availableElectricalPower, \
	totalElectricalPower, totalInputPower, vRedox, vBiasTotal, efficiencyElectronsToMtr, \
	efficiencyElectronsMtrToQuinone, energyPerFuelMolecule, electronsPerFuelInEET)
	
	efficiencyDict = deepcopy(outputDict)
	efficiencyDict['protonsPumpedOut'] = numberOfProtonsPumpedOut
	efficiencyDict['electronsPerFuel'] = electronsPerFuelInEET
	efficiencyDict['totalCells'] = totalCells
	efficiencyDict['electronRatePerCell'] = electronRatePerCell
	efficiencyDict['volumeBiofilm'] = volumeBiofilm
	efficiencyDict['preciseBiofilmThickness'] = preciseBiofilmThickness
	efficiencyDict['biofilmOverpotential'] = biofilmOverpotential
	efficiencyDict['areaBiofilm'] = areaBiofilm
	efficiencyDict['cellCurrentDensity'] = efficiencyDict['currentCell'] / areaBiofilm
	
# 	pdb.set_trace()
	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Efficiency_EET_BioCO2_to_Fuel_Direct_Contact(vRedox, vBias, vMembrane, \
vAcceptor, vMtr, vQuinone, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, \
carbonsFixedPerFuel, contactResistance, contactsPerCell, cellLength, cellDiameter, \
monoLayerFillFactor=0.25, efficiencyElectronsToMtr=1.0, efficiencyElectronsMtrToQuinone=1.0, \
numberOfProtonsPumpedInForATP=-1, totalElectricalPower=330.0, totalInputPower=1000.0, \
co2FixEnzymeRate=12, co2FixEnzymePerCell=1e6):

	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	from rewiredcarbon.physicalconstants import elementaryCharge as e
	from scipy.optimize import root
	import pdb
	from copy import deepcopy
	from numpy import pi
	
	if vMtr > vQuinone:
		print('Cannot transfer electron down due to cytochrome Mtr being more positive than' \
		+ ' menaquinone')
		sys.exit()

	availableElectricalPower = totalElectricalPower
	
	numberOfProtonsPumpedOut, electronsPerFuelInEET = \
	ProtonPumping_and_Electron_Requirements_for_EET_to_Fuel(vMembrane, vAcceptor, vMtr, vQuinone, \
	deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP)
	
	
	electronRatePerCell = Electron_Rate_Per_Cell(electronsPerFuelInEET, carbonsFixedPerFuel, \
	co2FixEnzymeRate, co2FixEnzymePerCell)
	
	# Calculate the cell cross sectional area. Assume a sphero-spherical model and calculate the 
	# projected cross section. 
	crossSectionalAreaCell = (cellLength - cellDiameter)*cellDiameter + pi*((cellDiameter/2)**2)
	
	# First, make a guess about the solution
	outputDict0 = \
	Efficiencies_Current_TotalEnergyContentOfFuel_EET_BioCO2(availableElectricalPower, \
	totalElectricalPower, totalInputPower, vRedox, vBias, efficiencyElectronsToMtr, \
	efficiencyElectronsMtrToQuinone, energyPerFuelMolecule, electronsPerFuelInEET)

	totalCells0 = outputDict0['quinoneCurrent'] / (e * electronRatePerCell)
	currentCell0 = outputDict0['currentCell']
	
	areaMonolayer0 = totalCells0 / monoLayerFillFactor
	
	contactResistance0 = contactsPerCell * totalCells0 * contactResistance
	
	contactResistanceOverpotential0 = contactResistance0 * currentCell0
	
	initialSolution = [currentCell0, contactResistanceOverpotential0, areaMonolayer0, totalCells0]
	
	solution = root(EET_Direct_Contact_Area_Overpotential_Coupled_Equations, initialSolution, \
	args=(availableElectricalPower, vRedox, vBias, electronRatePerCell, contactResistance, \
	contactsPerCell, efficiencyElectronsToMtr, efficiencyElectronsMtrToQuinone, \
	monoLayerFillFactor, crossSectionalAreaCell))
	
	currentCell, contactResistanceOverpotential, areaMonolayer, totalCells = solution.x
	
	outputDict = \
	Efficiencies_Current_TotalEnergyContentOfFuel_EET_BioCO2(availableElectricalPower, \
	totalElectricalPower, totalInputPower, vRedox, vBias, efficiencyElectronsToMtr, \
	efficiencyElectronsMtrToQuinone, energyPerFuelMolecule, electronsPerFuelInEET)
	
	efficiencyDict = deepcopy(outputDict)

	efficiencyDict['protonsPumpedOut'] = numberOfProtonsPumpedOut
	efficiencyDict['electronsPerFuel'] = electronsPerFuelInEET
	efficiencyDict['totalCells'] = totalCells
	efficiencyDict['electronRatePerCell'] = electronRatePerCell
	efficiencyDict['contactResistanceOverpotential'] = contactResistanceOverpotential
	efficiencyDict['areaMonolayer'] = areaMonolayer
	efficiencyDict['cellCurrentDensity'] = efficiencyDict['currentCell'] / areaMonolayer
	
	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Efficiency_Hydrogen_ElectrochemCO2_to_Fuel_No_ScaleUp(vRedoxCellOne, vBiasCellOne, \
vRedoxCellTwo, vBiasCellTwo, vMembrane, vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, \
energyPerFuelMolecule, carbonsPerFuel, electronsPerPrimaryFix, carbonsPerPrimaryFix, \
efficiencyCurrentToFirstCell=1, efficiencyCurrentToSecondCell=1, efficiencyCarbonToSecondCell=1, \
numberOfProtonsPumpedInForATP=-1, totalElectricalPower=330.0, totalInputPower=1000.0, \
stirPower=0.0, efficiencyElectronsToHyd=1.0):
 
    from rewiredcarbon.electrochemicalconstants import deltaGADPATP
    from copy import deepcopy
    import pdb
     
    if efficiencyCurrentToFirstCell == 0 or efficiencyCurrentToSecondCell == 0 or \
    efficiencyCarbonToSecondCell == 0:
        return 0
  
    availableElectricalPower = totalElectricalPower - stirPower
    vCellOne = vRedoxCellOne + vBiasCellOne
    vCellTwo = vRedoxCellTwo + vBiasCellTwo
 
    numberOfProtonsPumpedOut, electronsToAddInHyd = \
    ProtonPumping_and_Electron_Requirements_for_Hydrogen_to_Fuel(vMembrane, vAcceptor, \
    deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
    numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP)
         
    outputDict = \
    Efficiencies_Current_TotalEnergyContentOfFuel_Hydrogen_ElectrochemicalCO2(\
	availableElectricalPower, totalElectricalPower, totalInputPower, vCellOne, vCellTwo, \
	energyPerFuelMolecule, carbonsPerFuel, electronsToAddInHyd, electronsPerPrimaryFix, \
	carbonsPerPrimaryFix, efficiencyElectronsToHyd, efficiencyCurrentToFirstCell, \
	efficiencyCurrentToSecondCell, efficiencyCarbonToSecondCell)
 
    efficiencyDict = deepcopy(outputDict)
    efficiencyDict['protonsPumpedOut'] = numberOfProtonsPumpedOut
    efficiencyDict['electronsPerFuel'] = electronsToAddInHyd
     
    return efficiencyDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Efficiency_Hydrogen_ElectrochemCO2_to_Fuel_Simple_ScaleUp(vRedoxCellOne, vBiasCellOne, \
vRedoxCellTwo, vBiasCellTwo, vMembrane, vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, \
energyPerFuelMolecule, carbonsPerFuel, electronsPerPrimaryFix, carbonsPerPrimaryFix, \
efficiencyCurrentToFirstCell=1, efficiencyCurrentToSecondCell=1, efficiencyCarbonToSecondCell=1, \
numberOfProtonsPumpedInForATP=-1, totalElectricalPower=330.0, totalInputPower=1000.0, \
stirPower=0.0, efficiencyElectronsToHyd=1.0):

	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	import pdb
 
	if efficiencyCurrentToFirstCell == 0 or efficiencyCurrentToSecondCell == 0 or \
	efficiencyCarbonToSecondCell == 0:
		return 0

	availableElectricalPower = totalElectricalPower - stirPower
	vCellOne = vRedoxCellOne + vBiasCellOne
	vCellTwo = vRedoxCellTwo + vBiasCellTwo

	numberOfProtonsPumpedOut, electronsToAddInHyd = \
	ProtonPumping_and_Electron_Requirements_for_Hydrogen_to_Fuel(vMembrane, vAcceptor, \
	deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP)

	outputDict = \
	Efficiencies_Current_TotalEnergyContentOfFuel_Hydrogen_ElectrochemicalCO2(\
	availableElectricalPower, totalElectricalPower, totalInputPower, vCellOne, vCellTwo, \
	energyPerFuelMolecule, carbonsPerFuel, electronsToAddInHyd, electronsPerPrimaryFix, \
	carbonsPerPrimaryFix, efficiencyElectronsToHyd, efficiencyCurrentToFirstCell, \
	efficiencyCurrentToSecondCell, efficiencyCarbonToSecondCell)

	electronRatePerCell = Electron_Rate_Per_Cell(electronsToAddInHyd, primaryFixPerFuel, \
	primaryFixEnzymeRate, primaryFixEnzymesPerCell)

	hydrogenCurrent = outputDict['hydrogenCurrent']
 
	totalCells = hydrogenCurrent / (e * electronRatePerCell)

	efficiencyDict['protonsPumpedOut'] = numberOfProtonsPumpedOut
	efficiencyDict['electronsPerFuel'] = electronsToAddInHyd
	efficiencyDict['electronRatePerCell'] = electronRatePerCell
	efficiencyDict['totalCells'] = totalCells

	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Efficiency_EET_ElectrochemCO2_to_Fuel_No_ScaleUp(vRedoxCellOne, vBiasCellOne, vRedoxCellTwo, \
vBiasCellTwo, vMembrane, vAcceptor, vMtr, vQuinone, NADHforFuel, FdForFuel, ATPforFuel, \
energyPerFuelMolecule, carbonsPerFuel, electronsPerPrimaryFix, carbonsPerPrimaryFix, \
efficiencyCurrentToFirstCell=1, efficiencyCurrentToSecondCell=1, efficiencyCarbonToSecondCell=1, \
numberOfProtonsPumpedInForATP=-1, totalElectricalPower=330.0, totalInputPower=1000.0, \
stirPower=0.0, efficiencyElectronsToMtr=1.0, efficiencyElectronsMtrToQuinone=1.0):

	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	import pdb
	from copy import deepcopy
		
	if vMtr > vQuinone:
		print('Cannot transfer electron down due to cytochrome Mtr being more positive than' \
		+ ' menaquinone')
		sys.exit()
	
	if efficiencyCurrentToFirstCell == 0 or efficiencyCurrentToSecondCell == 0 or \
	efficiencyCarbonToSecondCell == 0:
 		return 0

	availableElectricalPower = totalElectricalPower - stirPower
	vCellOne = vRedoxCellOne + vBiasCellOne
	vCellTwo = vRedoxCellTwo + vBiasCellTwo

	numberOfProtonsPumpedOut, electronsToAddInEET = \
	ProtonPumping_and_Electron_Requirements_for_EET_to_Fuel(vMembrane, vAcceptor, vMtr, vQuinone, \
	deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP)
 	
	outputDict = \
	Efficiencies_Current_TotalEnergyContentOfFuel_EET_ElectrochemicalCO2(availableElectricalPower, \
	totalElectricalPower, totalInputPower, vCellOne, vCellTwo, energyPerFuelMolecule, \
	carbonsPerFuel, electronsToAddInEET, electronsPerPrimaryFix, carbonsPerPrimaryFix, \
	efficiencyElectronsToMtr, efficiencyElectronsMtrToQuinone, efficiencyCurrentToFirstCell, \
	efficiencyCurrentToSecondCell, efficiencyCarbonToSecondCell)

	efficiencyDict = deepcopy(outputDict)	
	efficiencyDict['protonsPumpedOut'] = numberOfProtonsPumpedOut
	efficiencyDict['electronsPerFuel'] = electronsToAddInEET

	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Efficiency_EET_ElectrochemCO2_to_Fuel_Simple_ScaleUp(vRedoxCellOne, vBiasCellOne, vRedoxCellTwo, \
vBiasCellTwo, vMembrane, vAcceptor, vMtr, vQuinone, NADHforFuel, FdForFuel, ATPforFuel, \
energyPerFuelMolecule, carbonsPerFuel, electronsPerPrimaryFix, primaryFixPerFuel, \
efficiencyCurrentToFirstCell=1, efficiencyCurrentToSecondCell=1, \
efficiencyCarbonToSecondCell=1, numberOfProtonsPumpedInForATP=-1, totalElectricalPower=330.0, \
totalInputPower=1000.0, stirPower=0.0, primaryFixEnzymeRate=12, primaryFixEnzymesPerCell=1e6):

	from rewiredcarbon.electrochemicalconstants import deltaGADPATP
	import pdb
		
	if vMtr > vQuinone:
		print('Cannot transfer electron down due to cytochrome Mtr being more positive than' \
		+ ' menaquinone')
		sys.exit()
	
	if efficiencyCurrentToFirstCell == 0 or efficiencyCurrentToSecondCell == 0 or \
	efficiencyCarbonToSecondCell == 0:
 		return 0

	availableElectricalPower = totalElectricalPower - stirPower
	vCellOne = vRedoxCellOne + vBiasCellOne
	vCellTwo = vRedoxCellTwo + vBiasCellTwo

	numberOfProtonsPumpedOut, electronsToAddInEET = \
	ProtonPumping_and_Electron_Requirements_for_EET_to_Fuel(vMembrane, vAcceptor, vMtr, vQuinone, \
	deltaGADPATP, NADHforFuel, FdForFuel, ATPforFuel, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP)
 	
	outputDict = \
	Efficiencies_Current_TotalEnergyContentOfFuel_EET_ElectrochemicalCO2(availableElectricalPower, \
	totalElectricalPower, totalInputPower, vRedox, vBias, energyPerFuelMolecule, carbonsPerFuel, \
	electronsToAddInEET, electronsPerPrimaryFix, carbonsPerPrimaryFix, efficiencyElectronsToMtr, \
	efficiencyElectronsMtrToQuinone, efficiencyCurrentToSecondCell, efficiencyCarbonToSecondCell)
	
	electronRatePerCell = Electron_Rate_Per_Cell(electronsToAddInEET, primaryFixPerFuel, \
	primaryFixEnzymeRate, primaryFixEnzymesPerCell)
	
	quinoneCurrent = outputDict['quinoneCurrent']
			
	totalCells = quinoneCurrent / (e * electronRatePerCell)
		
	efficiencyDict['protonsPumpedOut'] = numberOfProtonsPumpedOut
	efficiencyDict['electronsPerFuel'] = electronsToAddInEET
	efficiencyDict['electronRatePerCell'] = electronRatePerCell
	efficiencyDict['totalCells'] = totalCells
	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Process_Hydrogen_with_BioCO2_Scenario(scenarioData):
	
	import pdb
	from rewiredcarbon.scaleup import ImportAndExtraplotePowerNumberCurve
	import sys
	
	# These are variables that every scale up scenario needs
	NADHforFuel = float(scenarioData['NADHforFuel'])
	FdForFuel = float(scenarioData['FdForFuel'])
	ATPforFuel = float(scenarioData['ATPforFuel'])	
	vAcceptor = float(scenarioData['vAcceptor'])
	energyPerFuelMolecule = float(scenarioData['energyPerFuelMolecule'])
	voltageCellTwoCathode = float(scenarioData['voltageCellTwoCathode'])
	voltageCellTwoAnode = float(scenarioData['voltageCellTwoAnode'])
	voltageCellTwoCathodeBias = float(scenarioData['voltageCellTwoCathodeBias'])
	voltageCellTwoAnodeBias = float(scenarioData['voltageCellTwoAnodeBias'])
	vMembrane = float(scenarioData['voltageMembrane'])/1000.
	totalElectricalPower = float(scenarioData['totalElectricalPower'])
	totalInputPower = float(scenarioData['totalInputPower'])
	
	try:
		stirPower = float(scenarioData['stirPower'])
	except:
		stirPower = 0
			
	vRedox = abs(float(voltageCellTwoCathode) - float(voltageCellTwoAnode))
	vBias = abs(float(voltageCellTwoCathodeBias) + float(voltageCellTwoAnodeBias))
	
	scaleUpMode = scenarioData['scaleUpMode']
	
	if scaleUpMode == 'None':
		efficiencyDict = \
		Efficiency_Hydrogen_BioCO2_to_Fuel_No_ScaleUp(vRedox, vBias, vMembrane, \
		vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, \
		efficiencyElectronsToHyd=1.0, numberOfProtonsPumpedInForATP=-1, 
		stirPower=stirPower, totalElectricalPower=totalElectricalPower, \
		totalInputPower=totalInputPower)
	
	elif scaleUpMode == 'Simple':
		carbonsFixedPerFuel = float(scenarioData['carbonsFixedPerFuel']) 
		co2FixEnzymeRate = float(scenarioData['co2FixEnzymeRate'])
		co2FixEnzymePerCell = float(scenarioData['co2FixEnzymePerCell'])
	
		efficiencyDict = \
		Efficiency_Hydrogen_BioCO2_to_Fuel_Simple_ScaleUp(vRedox, vBias, vMembrane, vAcceptor, \
		NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, carbonsFixedPerFuel, \
		efficiencyElectronsToHyd=1.0, numberOfProtonsPumpedInForATP=-1, stirPower=stirPower, \
		totalElectricalPower=totalElectricalPower, totalInputPower=totalInputPower, \
		co2FixEnzymeRate=co2FixEnzymeRate, co2FixEnzymePerCell=co2FixEnzymePerCell)
	
	elif scaleUpMode == 'Stirring with Fixed Density':
		carbonsFixedPerFuel = float(scenarioData['carbonsFixedPerFuel']) 
		co2FixEnzymeRate = float(scenarioData['co2FixEnzymeRate'])
		co2FixEnzymePerCell = float(scenarioData['co2FixEnzymePerCell'])
		cellDensity = float(scenarioData['cellDensity'])
		agitatorReynoldsNumber = float(scenarioData['agitatorReynoldsNumber'])

		powerNumberFunction = \
		ImportAndExtraplotePowerNumberCurve(scenarioData['powerNumberCurveFileName'])

		efficiencyDict = \
		Efficiency_Hydrogen_BioCO2_to_Fuel_ScaleUp_Stirring_Fixed_CellDensity(vRedox, vBias, \
		vMembrane, vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, \
		carbonsFixedPerFuel, powerNumberFunction, agitatorReynoldsNumber, \
		efficiencyElectronsToHyd=1.0, numberOfProtonsPumpedInForATP=-1, \
		totalElectricalPower=totalElectricalPower, totalInputPower=totalInputPower, \
		co2FixEnzymeRate=co2FixEnzymeRate, co2FixEnzymePerCell=co2FixEnzymePerCell, \
		cellDensity=cellDensity, \
		fluidDensity=997, fluidViscosity=9.1e-4, tankDiameterToHeightRatio=1, \
		agitatorDiameterToTankDiameterRatio=0.5)
		
	elif scaleUpMode == 'Stirring with Fixed Density Grenville':
		carbonsFixedPerFuel = float(scenarioData['carbonsFixedPerFuel']) 
		co2FixEnzymeRate = float(scenarioData['co2FixEnzymeRate'])
		co2FixEnzymePerCell = float(scenarioData['co2FixEnzymePerCell'])
		cellDensity = float(scenarioData['cellDensity'])
		pressureHydrogen = float(scenarioData['pressureHydrogen'])
		temperature = float(scenarioData['cathodeTemperature'])
		bulkToPeakHydrogenConcentrationRatio = float(scenarioData['bulkToPeakHydrogenConcentrationRatio'])
		Aconstant = float(scenarioData['Aconstant'])
		Bconstant = float(scenarioData['Bconstant'])
		Cconstant = float(scenarioData['Cconstant'])
		
		efficiencyDict = \
		Efficiency_Hydrogen_BioCO2_to_Fuel_ScaleUp_Stirring_Fixed_CellDensity_Grenville(vRedox, \
		vBias, vMembrane, vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, \
		carbonsFixedPerFuel, \
		pressureHydrogen=pressureHydrogen, temperature=temperature, \
		bulkToPeakHydrogenConcentrationRatio=bulkToPeakHydrogenConcentrationRatio, \
		tankDiameterToHeightRatio=1, Aconstant=Aconstant, Bconstant=Bconstant, \
		Cconstant=Cconstant, \
		efficiencyElectronsToHyd=1.0, totalElectricalPower=totalElectricalPower, \
		totalInputPower=totalInputPower, \
		numberOfProtonsPumpedInForATP=-1, co2FixEnzymeRate=co2FixEnzymeRate, \
		co2FixEnzymePerCell=co2FixEnzymePerCell, \
		cellDensity=cellDensity)
		
	elif scaleUpMode == 'Initial Stir Power Guess':
		carbonsFixedPerFuel = float(scenarioData['carbonsFixedPerFuel']) 
		co2FixEnzymeRate = float(scenarioData['co2FixEnzymeRate'])
		co2FixEnzymePerCell = float(scenarioData['co2FixEnzymePerCell'])
		cellDensity = float(scenarioData['cellDensity'])
		pressureHydrogen = float(scenarioData['pressureHydrogen'])
		temperature = float(scenarioData['cathodeTemperature'])
		bulkToPeakHydrogenConcentrationRatio = float(scenarioData['bulkToPeakHydrogenConcentrationRatio'])
		Aconstant = float(scenarioData['Aconstant'])
		Bconstant = float(scenarioData['Bconstant'])
		Cconstant = float(scenarioData['Cconstant'])
		
		efficiencyDict = \
		Guess_Initial_Stir_Power(vRedox, \
		vBias, vMembrane, vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, \
		carbonsFixedPerFuel, \
		pressureHydrogen=pressureHydrogen, temperature=temperature, \
		bulkToPeakHydrogenConcentrationRatio=bulkToPeakHydrogenConcentrationRatio, \
		tankDiameterToHeightRatio=1, Aconstant=Aconstant, Bconstant=Bconstant, \
		Cconstant=Cconstant, \
		efficiencyElectronsToHyd=1.0, totalElectricalPower=totalElectricalPower, \
		totalInputPower=totalInputPower, \
		numberOfProtonsPumpedInForATP=-1, co2FixEnzymeRate=co2FixEnzymeRate, \
		co2FixEnzymePerCell=co2FixEnzymePerCell, \
		cellDensity=cellDensity)
	
	
	elif scaleUpMode == 'Diffusion with Fixed Density':
		carbonsFixedPerFuel = float(scenarioData['carbonsFixedPerFuel']) 
		co2FixEnzymeRate = float(scenarioData['co2FixEnzymeRate'])
		co2FixEnzymePerCell = float(scenarioData['co2FixEnzymePerCell'])
		cellDensity = float(scenarioData['cellDensity'])
		cathodeSelfExchangeCurrentDensity = float(scenarioData['cathodeSelfExchangeCurrentDensity'])
		cathodeTafelAlpha = float(scenarioData['cathodeTafelAlpha'])
		cathodeTemperature = float(scenarioData['cathodeTemperature'])
		pressureHydrogen = float(scenarioData['pressureHydrogen'])
		
		
		efficiencyDict = \
		Efficiency_Hydrogen_BioCO2_to_Fuel_ScaleUp_Diffusion_Fixed_CellDensity(vRedox, vBias, \
		vMembrane, vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, \
		carbonsFixedPerFuel, \
		efficiencyElectronsToHyd=1.0, numberOfProtonsPumpedInForATP=-1, \
		totalElectricalPower=totalElectricalPower, totalInputPower=totalInputPower, \
		co2FixEnzymeRate=co2FixEnzymeRate, co2FixEnzymePerCell=co2FixEnzymePerCell, \
		cellDensity=cellDensity, \
		cathodeSelfExchangeCurrentDensity=cathodeSelfExchangeCurrentDensity, \
		cathodeTafelAlpha=cathodeTafelAlpha, cathodeTemperature=cathodeTemperature, \
		pressureHydrogen=pressureHydrogen, solubilityHydrogen=1.282e5, diffusionHydrogen=4.5e-9)
		
		
	elif scaleUpMode == 'Diffusion with Fixed Density and Fixed Electrode Overpotentials':
		carbonsFixedPerFuel = float(scenarioData['carbonsFixedPerFuel']) 
		co2FixEnzymeRate = float(scenarioData['co2FixEnzymeRate'])
		co2FixEnzymePerCell = float(scenarioData['co2FixEnzymePerCell'])
		cellDensity = float(scenarioData['cellDensity'])
		pressureHydrogen = float(scenarioData['pressureHydrogen'])
		
		
		efficiencyDict = \
		Efficiency_Hydrogen_BioCO2_to_Fuel_ScaleUp_Diffusion_Fixed_CellDensity_Fixed_Overpotential(\
		vRedox, vBias, vMembrane, vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, \
		energyPerFuelMolecule, carbonsFixedPerFuel, \
		efficiencyElectronsToHyd=1.0, numberOfProtonsPumpedInForATP=-1, \
		totalElectricalPower=totalElectricalPower, totalInputPower=totalInputPower, \
		co2FixEnzymeRate=co2FixEnzymeRate, co2FixEnzymePerCell=co2FixEnzymePerCell, \
		cellDensity=cellDensity, pressureHydrogen=pressureHydrogen, solubilityHydrogen=1.282e5, \
		diffusionHydrogen=4.5e-9)
		
	
		
		
	else:
		print("Can't do scale up modes other than None; Simple; Stirring with Fixed Density;" \
		+ " or Diffusion with Fixed Density right now")
		sys.exit()
	
	
		
	return efficiencyDict 
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Process_EET_with_BioCO2_Scenario(scenarioData):
	
	import pdb
	
	NADHforFuel = float(scenarioData['NADHforFuel'])
	FdForFuel = float(scenarioData['FdForFuel'])
	
	ATPforFuel = float(scenarioData['ATPforFuel'])	
	vAcceptor = float(scenarioData['vAcceptor'])
	energyPerFuelMolecule = float(scenarioData['energyPerFuelMolecule'])
	vMtr = float(scenarioData['vMtr'])
	vQuinone = float(scenarioData['vQuinone'])
	
	voltageCellTwoCathode = float(scenarioData['voltageCellTwoCathode'])
	voltageCellTwoAnode = float(scenarioData['voltageCellTwoAnode'])
	voltageCellTwoCathodeBias = float(scenarioData['voltageCellTwoCathodeBias'])
	voltageCellTwoAnodeBias = float(scenarioData['voltageCellTwoAnodeBias'])
	vMembrane = float(scenarioData['voltageMembrane'])/1000.
	
	totalElectricalPower = float(scenarioData['totalElectricalPower'])
	totalInputPower = float(scenarioData['totalInputPower'])
	
	vRedox = abs(float(voltageCellTwoCathode) - float(voltageCellTwoAnode))
	vBias = abs(float(voltageCellTwoCathodeBias) + float(voltageCellTwoAnodeBias))
	
	scaleUpMode = scenarioData['scaleUpMode']
	
	if scaleUpMode.lower() == 'none':
		efficiencyDict = \
		Efficiency_EET_BioCO2_to_Fuel_No_ScaleUp(vRedox, vBias, vMembrane, \
		vAcceptor, vMtr, vQuinone, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, \
		efficiencyElectronsToMtr=1.0, efficiencyElectronsMtrToQuinone=1.0, \
		numberOfProtonsPumpedInForATP=-1, totalElectricalPower=totalElectricalPower, \
		totalInputPower=totalInputPower)
		
	elif scaleUpMode.lower() == 'simple':
		carbonsFixedPerFuel = float(scenarioData['carbonsFixedPerFuel']) 
		co2FixEnzymeRate = float(scenarioData['co2FixEnzymeRate'])
		co2FixEnzymePerCell = float(scenarioData['co2FixEnzymePerCell'])
		
		efficiencyDict = \
		Efficiency_EET_BioCO2_to_Fuel_Simple_ScaleUp(vRedox, vBias, vMembrane, vAcceptor, vMtr, \
		vQuinone, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, carbonsFixedPerFuel, \
		efficiencyElectronsToMtr=1.0, efficiencyElectronsMtrToQuinone=1.0, \
		numberOfProtonsPumpedInForATP=-1, totalElectricalPower=totalElectricalPower, \
		totalInputPower=totalInputPower, co2FixEnzymeRate=co2FixEnzymeRate, \
		co2FixEnzymePerCell=co2FixEnzymePerCell)
			
	
	elif scaleUpMode.lower() == 'conductive matrix with fixed height':
		carbonsFixedPerFuel = float(scenarioData['carbonsFixedPerFuel']) 
		co2FixEnzymeRate = float(scenarioData['co2FixEnzymeRate'])
		co2FixEnzymePerCell = float(scenarioData['co2FixEnzymePerCell'])
		approximateBiofilmThickness = float(scenarioData['approximateBiofilmThickness']) * 1e-6
		resistivityBiofilm = float(scenarioData['resistivityBiofilm']) / 100
		cellLength = float(scenarioData['cellLength']) * 1e-6
		cellDiameter = float(scenarioData['cellDiameter']) * 1e-6
		biofilmLayerHeight = float(scenarioData['biofilmLayerHeight']) * 1e-6
		biofilmLayerFillFactor = float(scenarioData['biofilmLayerFillFactor'])
		
		efficiencyDict = \
		Efficiency_EET_BioCO2_to_Fuel_Conductive_Matrix_Fixed_Height(vRedox, vBias, vMembrane, \
		vAcceptor, vMtr, vQuinone, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, \
		carbonsFixedPerFuel, approximateBiofilmThickness, resistivityBiofilm, \
		cellLength, cellDiameter, biofilmLayerHeight, \
		biofilmLayerFillFactor=biofilmLayerFillFactor, efficiencyElectronsToMtr=1.0, \
		efficiencyElectronsMtrToQuinone=1.0, numberOfProtonsPumpedInForATP=-1, \
		totalElectricalPower=totalElectricalPower, totalInputPower=totalInputPower, \
		co2FixEnzymeRate=co2FixEnzymeRate, co2FixEnzymePerCell=co2FixEnzymePerCell)
	
	elif scaleUpMode.lower() == 'direct contact':
		carbonsFixedPerFuel = float(scenarioData['carbonsFixedPerFuel']) 
		co2FixEnzymeRate = float(scenarioData['co2FixEnzymeRate'])
		co2FixEnzymePerCell = float(scenarioData['co2FixEnzymePerCell'])
		cellLength = float(scenarioData['cellLength']) * 1e-6
		cellDiameter = float(scenarioData['cellDiameter']) * 1e-6
		contactResistance = float(scenarioData['contactResistance'])
		contactsPerCell = float(scenarioData['contactsPerCell'])
		monoLayerFillFactor = float(scenarioData['monoLayerFillFactor'])
	
		efficiencyDict = \
		Efficiency_EET_BioCO2_to_Fuel_Direct_Contact(vRedox, vBias, vMembrane, \
		vAcceptor, vMtr, vQuinone, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, \
		carbonsFixedPerFuel, contactResistance, contactsPerCell, cellLength, cellDiameter, \
		monoLayerFillFactor=monoLayerFillFactor, efficiencyElectronsToMtr=1.0, \
		efficiencyElectronsMtrToQuinone=1.0, numberOfProtonsPumpedInForATP=-1, \
		totalElectricalPower=330.0, totalInputPower=1000.0, co2FixEnzymeRate=12, \
		co2FixEnzymePerCell=1e6)
	
	else:
		print("Can't do scale up modes other than None; Simple; Direct Contact;" \
		+ " or Conductive Matrix with Fixed Height right now")
		sys.exit()
	
	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Process_Hydrogen_with_ElectrochemicalCO2_Scenario(scenarioData):
	
	
	NADHforFuel = float(scenarioData['NADHforFuel'])
	FdForFuel = float(scenarioData['FdForFuel'])
	
	ATPforFuel = float(scenarioData['ATPforFuel'])	
	vAcceptor = float(scenarioData['vAcceptor'])
	energyPerFuelMolecule = float(scenarioData['energyPerFuelMolecule'])
	carbonsPerFuel = float(scenarioData['carbonsPerFuel'])
	
	electronsPerPrimaryFix = float(scenarioData['electronsPerPrimaryFix'])
	carbonsPerPrimaryFix = float(scenarioData['carbonsPerPrimaryFix'])
	
	voltageCellTwoCathode = float(scenarioData['voltageCellTwoCathode'])
	voltageCellTwoAnode = float(scenarioData['voltageCellTwoAnode'])
	voltageCellTwoCathodeBias = float(scenarioData['voltageCellTwoCathodeBias'])
	voltageCellTwoAnodeBias = float(scenarioData['voltageCellTwoAnodeBias'])
	
	voltageCellOneCathode = float(scenarioData['voltageCellOneCathode'])
	voltageCellOneAnode = float(scenarioData['voltageCellOneAnode'])
	voltageCellOneCathodeBias = float(scenarioData['voltageCellOneCathodeBias'])
	voltageCellOneAnodeBias = float(scenarioData['voltageCellOneAnodeBias'])
	
	vMembrane = float(scenarioData['voltageMembrane'])/1000.
	
	stirPower = float(scenarioData['stirPower'])
	totalElectricalPower = float(scenarioData['totalElectricalPower'])
	totalInputPower = float(scenarioData['totalInputPower'])

	vRedoxCellTwo = abs(float(voltageCellTwoCathode) - float(voltageCellTwoAnode))
	vBiasCellTwo = abs(float(voltageCellTwoCathodeBias) + float(voltageCellTwoAnodeBias))
	
	vRedoxCellOne = abs(float(voltageCellOneCathode) - float(voltageCellOneAnode))
	vBiasCellOne = abs(float(voltageCellOneCathodeBias) + float(voltageCellOneAnodeBias))
	
	efficiencyCurrentToFirstCell = float(scenarioData['efficiencyCurrentToFirstCell'])

	
	scaleUpMode = scenarioData['scaleUpMode']
	
	if scaleUpMode.lower() == 'simple':
		primaryFixPerFuel = float(scenarioData['primaryFixPerFuel']) 
		primaryFixEnzymeRate = float(scenarioData['primaryFixEnzymeRate'])
		primaryFixEnzymesPerCell = float(scenarioData['primaryFixEnzymesPerCell'])
	
		efficiencyDict = \
		Efficiency_Hydrogen_ElectrochemCO2_to_Fuel_Simple_ScaleUp(vRedoxCellOne, vBiasCellOne, \
		vRedoxCellTwo, vBiasCellTwo, vMembrane, vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, \
		energyPerFuelMolecule, carbonsPerFuel, electronsPerPrimaryFix, carbonsPerPrimaryFix, \
		efficiencyCurrentToFirstCell=efficiencyCurrentToFirstCell, efficiencyCurrentToSecondCell=1, \
		efficiencyCarbonToSecondCell=1, numberOfProtonsPumpedInForATP=-1, \
		totalElectricalPower=totalElectricalPower, totalInputPower=totalInputPower, \
		stirPower=stirPower, efficiencyElectronsToHyd=1.0)
		
		
	elif scaleUpMode.lower() == 'none':
		efficiencyDict = \
		Efficiency_Hydrogen_ElectrochemCO2_to_Fuel_No_ScaleUp(vRedoxCellOne, vBiasCellOne, \
		vRedoxCellTwo, vBiasCellTwo, vMembrane, vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, \
		energyPerFuelMolecule, carbonsPerFuel, electronsPerPrimaryFix, carbonsPerPrimaryFix, \
		efficiencyCurrentToFirstCell=efficiencyCurrentToFirstCell, efficiencyCurrentToSecondCell=1,\
		efficiencyCarbonToSecondCell=1, numberOfProtonsPumpedInForATP=-1, \
		totalElectricalPower=totalElectricalPower, totalInputPower=totalInputPower, \
		stirPower=stirPower, efficiencyElectronsToHyd=1.0)
	
	else:
		print('We can\'t do electrochemical CO2 fixation scale up modes other than ' \
		+ 'None or Simple right now')
		sys.exit

	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Process_EET_with_ElectrochemicalCO2_Scenario(scenarioData):
	
	import pdb
	import sys
	
	NADHforFuel = float(scenarioData['NADHforFuel'])
	FdForFuel = float(scenarioData['FdForFuel'])
	ATPforFuel = float(scenarioData['ATPforFuel'])	
	
	vAcceptor = float(scenarioData['vAcceptor'])
	energyPerFuelMolecule = float(scenarioData['energyPerFuelMolecule'])
	carbonsPerFuel = float(scenarioData['carbonsPerFuel'])
	vMtr = float(scenarioData['vMtr'])
	vQuinone = float(scenarioData['vQuinone'])
	
	electronsPerPrimaryFix = float(scenarioData['electronsPerPrimaryFix'])
	carbonsPerPrimaryFix = float(scenarioData['carbonsPerPrimaryFix'])
	
	voltageCellTwoCathode = float(scenarioData['voltageCellTwoCathode'])
	voltageCellTwoAnode = float(scenarioData['voltageCellTwoAnode'])
	voltageCellTwoCathodeBias = float(scenarioData['voltageCellTwoCathodeBias'])
	voltageCellTwoAnodeBias = float(scenarioData['voltageCellTwoAnodeBias'])
	
	voltageCellOneCathode = float(scenarioData['voltageCellOneCathode'])
	voltageCellOneAnode = float(scenarioData['voltageCellOneAnode'])
	voltageCellOneCathodeBias = float(scenarioData['voltageCellOneCathodeBias'])
	voltageCellOneAnodeBias = float(scenarioData['voltageCellOneAnodeBias'])
	
	
	vMembrane = float(scenarioData['voltageMembrane'])/1000.
	
	stirPower = float(scenarioData['stirPower'])
	totalElectricalPower = float(scenarioData['totalElectricalPower'])
	totalInputPower = float(scenarioData['totalInputPower'])
	
	vRedoxCellTwo = abs(float(voltageCellTwoCathode) - float(voltageCellTwoAnode))
	vBiasCellTwo = abs(float(voltageCellTwoCathodeBias) + float(voltageCellTwoAnodeBias))
	
	vRedoxCellOne = abs(float(voltageCellOneCathode) - float(voltageCellOneAnode))
	vBiasCellOne = abs(float(voltageCellOneCathodeBias) + float(voltageCellOneAnodeBias))
	
	efficiencyCurrentToFirstCell = float(scenarioData['efficiencyCurrentToFirstCell'])

	scaleUpMode = scenarioData['scaleUpMode']
		
	if scaleUpMode.lower() == 'simple':	
		primaryFixPerFuel = float(scenarioData['primaryFixPerFuel']) 
		primaryFixEnzymeRate = float(scenarioData['primaryFixEnzymeRate'])
		primaryFixEnzymesPerCell = float(scenarioData['primaryFixEnzymesPerCell'])
		
		
		efficiencyDict = \
		Efficiency_EET_ElectrochemCO2_to_Fuel_Simple_ScaleUp(vRedoxCellOne, vBiasCellOne, \
		vRedoxCellTwo, vBiasCellTwo, vMembrane, vAcceptor, vMtr, vQuinone, NADHforFuel, FdForFuel, \
		ATPforFuel, energyPerFuelMolecule, carbonsPerFuel, electronsPerPrimaryFix, \
		primaryFixPerFuel, \
		efficiencyCurrentToFirstCell=efficiencyCurrentToFirstCell, efficiencyCurrentToSecondCell=1, \
		efficiencyCarbonToSecondCell=1, numberOfProtonsPumpedInForATP=-1, \
		totalElectricalPower=totalElectricalPower, \
		totalInputPower=totalInputPower, stirPower=stirPower, \
		primaryFixEnzymeRate=primaryFixEnzymeRate, \
		primaryFixEnzymesPerCell=primaryFixEnzymesPerCell)
	
	elif scaleUpMode.lower() == 'none':
		primaryFixPerFuel = -1 
		primaryFixEnzymeRate = -1
		primaryFixEnzymesPerCell = -1
		
		efficiencyDict = \
		Efficiency_EET_ElectrochemCO2_to_Fuel_No_ScaleUp(vRedoxCellOne, vBiasCellOne, \
		vRedoxCellTwo, vBiasCellTwo, vMembrane, vAcceptor, vMtr, vQuinone, NADHforFuel, FdForFuel, \
		ATPforFuel, energyPerFuelMolecule, carbonsPerFuel, electronsPerPrimaryFix, \
		carbonsPerPrimaryFix, \
		efficiencyCurrentToFirstCell=efficiencyCurrentToFirstCell, efficiencyCurrentToSecondCell=1,\
		efficiencyCarbonToSecondCell=1, numberOfProtonsPumpedInForATP=-1, \
		totalElectricalPower=totalElectricalPower, totalInputPower=totalInputPower, \
		stirPower=stirPower)
	
	else:
		print('We can\'t do electrochemical CO2 fixation modes other than None or Simple right now')
		sys.exit
		
	
	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Efficiency_Hydrogen_BioCO2_to_Fuel_ScaleUp_Stirring_Fixed_CellDensity_with_Packed_ScenarioData(\
scenarioData, outputKeyToReport):
	
	from rewiredcarbon.efficiency import \
	Efficiency_Hydrogen_BioCO2_to_Fuel_ScaleUp_Stirring_Fixed_CellDensity
	
	from rewiredcarbon.scaleup import ImportAndExtraplotePowerNumberCurve
	
	import pdb

	# Variables that are set to defaults and aren't included in the scenarioData
	efficiencyElectronsToHyd = 1
	numberOfProtonsPumpedInForATP = -1
	
	# Figure out if we need to vary the dependent variable on a logarithmic or linear scale
	dependentVariableKey = scenarioData['dependentVariableKey']
	dependentVariableScale = scenarioData['dependentVariableScale']
	
	if dependentVariableScale.lower() == 'logarithmic':
		scenarioData[dependentVariableKey] = 10**(float(scenarioData[dependentVariableKey]))
	
	# Unpack the scenarioData
	voltageCellTwoCathode = float(scenarioData['voltageCellTwoCathode'])
	voltageCellTwoCathodeBias = float(scenarioData['voltageCellTwoCathodeBias'])
	voltageCellTwoAnode = float(scenarioData['voltageCellTwoAnode'])
	voltageCellTwoAnodeBias = float(scenarioData['voltageCellTwoAnodeBias'])
	vMembrane = float(scenarioData['voltageMembrane'])/1000
	vAcceptor = float(scenarioData['vAcceptor'])
	NADHforFuel = float(scenarioData['NADHforFuel'])
	FdForFuel = float(scenarioData['FdForFuel'])
	ATPforFuel = float(scenarioData['ATPforFuel'])
	energyPerFuelMolecule = float(scenarioData['energyPerFuelMolecule'])
	carbonsFixedPerFuel = float(scenarioData['carbonsFixedPerFuel'])
	totalElectricalPower = float(scenarioData['totalElectricalPower'])
	totalInputPower = float(scenarioData['totalInputPower'])
	co2FixEnzymeRate = float(scenarioData['co2FixEnzymeRate'])
	co2FixEnzymePerCell = float(scenarioData['co2FixEnzymePerCell'])
	cellDensity = float(scenarioData['cellDensity'])
	agitatorReynoldsNumber = float(scenarioData['agitatorReynoldsNumber'])

	simulatePowerCurve = bool(scenarioData['simulatePowerCurve'])

	if simulatePowerCurve == False:
		powerNumberFunction = \
		ImportAndExtraplotePowerNumberCurve(scenarioData['powerNumberCurveFileName'])
	else:
		powerNumberFunction = lambda x:float(scenarioData['powerNumberConstant'])

	vRedox = abs(float(voltageCellTwoCathode) - float(voltageCellTwoAnode))
	vBias = abs(float(voltageCellTwoCathodeBias) + float(voltageCellTwoAnodeBias))
	
	efficiencyDict = Efficiency_Hydrogen_BioCO2_to_Fuel_ScaleUp_Stirring_Fixed_CellDensity(vRedox, \
	vBias, vMembrane,\
	vAcceptor, NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, carbonsFixedPerFuel, \
	powerNumberFunction, agitatorReynoldsNumber, \
	efficiencyElectronsToHyd=efficiencyElectronsToHyd, \
	numberOfProtonsPumpedInForATP=numberOfProtonsPumpedInForATP, \
	totalElectricalPower=totalElectricalPower, \
	totalInputPower=totalInputPower, co2FixEnzymeRate=co2FixEnzymeRate, \
	co2FixEnzymePerCell=co2FixEnzymePerCell, cellDensity=cellDensity, \
	fluidDensity=997, fluidViscosity=9.1e-4, tankDiameterToHeightRatio=1, \
	agitatorDiameterToTankDiameterRatio=0.5)

	efficiency = efficiencyDict[outputKeyToReport]
	
	return efficiency
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Efficiency_Hydrogen_BioCO2_to_Fuel_ScaleUp_Stirring_Fixed_CellDensity_to_Target_Difference(X, \
dependentVariableKey, efficiencyTarget, scenarioData):

	import pdb

	scenarioData[dependentVariableKey] = X
	
	efficiency = \
	Efficiency_Hydrogen_BioCO2_to_Fuel_ScaleUp_Stirring_Fixed_CellDensity_with_Packed_ScenarioData(\
	scenarioData, 'effTotalElectricalToFuel')
	
	try:
		diff = efficiencyTarget - efficiency
	except:
		pdb.set_trace()

	return diff
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Efficiency_Hydrogen_BioCO2_to_Fuel_Simple_ScaleUp_with_Packed_ScenarioData(scenarioData, \
outputKeyToReport):

	from rewiredcarbon.efficiency import Efficiency_Hydrogen_BioCO2_to_Fuel_Simple_ScaleUp
	import pdb
	
	# Unpack the scenarioData
	voltageCellTwoCathode = float(scenarioData['voltageCellTwoCathode'])
	voltageCellTwoCathodeBias = float(scenarioData['voltageCellTwoCathodeBias'])
	voltageCellTwoAnode = float(scenarioData['voltageCellTwoAnode'])
	voltageCellTwoAnodeBias = float(scenarioData['voltageCellTwoAnodeBias'])
	vMembrane = float(scenarioData['voltageMembrane'])/1000
	vAcceptor = float(scenarioData['vAcceptor'])
	NADHforFuel = float(scenarioData['NADHforFuel'])
	FdForFuel = float(scenarioData['FdForFuel'])
	ATPforFuel = float(scenarioData['ATPforFuel'])
	energyPerFuelMolecule = float(scenarioData['energyPerFuelMolecule'])
	carbonsFixedPerFuel = float(scenarioData['carbonsFixedPerFuel'])
	totalElectricalPower = float(scenarioData['totalElectricalPower'])
	totalInputPower = float(scenarioData['totalInputPower'])
	co2FixEnzymeRate = float(scenarioData['co2FixEnzymeRate'])
	co2FixEnzymePerCell = float(scenarioData['co2FixEnzymePerCell'])
	
	try:
		stirPower = float(scenarioData['stirPower'])
	except:
		stirPower = 0

	vRedox = abs(float(voltageCellTwoCathode) - float(voltageCellTwoAnode))
	vBias = abs(float(voltageCellTwoCathodeBias) + float(voltageCellTwoAnodeBias))

# 	pdb.set_trace()


	efficiencyDict = \
	Efficiency_Hydrogen_BioCO2_to_Fuel_Simple_ScaleUp(vRedox, vBias, vMembrane, vAcceptor, \
	NADHforFuel, FdForFuel, ATPforFuel, energyPerFuelMolecule, carbonsFixedPerFuel, \
	efficiencyElectronsToHyd=1.0, numberOfProtonsPumpedInForATP=-1, stirPower=stirPower, \
	totalElectricalPower=totalElectricalPower, totalInputPower=totalInputPower, \
	co2FixEnzymeRate=co2FixEnzymeRate, co2FixEnzymePerCell=co2FixEnzymePerCell)

	efficiency = efficiencyDict[outputKeyToReport]
	
# 	pdb.set_trace()
	
	return efficiency
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Efficiency_Hydrogen_BioCO2_to_Fuel_Simple_ScaleUp_to_Target_Difference(X, dependentVariableKey,\
efficiencyTarget, scenarioData):

	import pdb

	scenarioData[dependentVariableKey] = X
	
	efficiency = \
	Efficiency_Hydrogen_BioCO2_to_Fuel_Simple_ScaleUp_with_Packed_ScenarioData(scenarioData, \
	'effTotalElectricalToFuel')
		
	try:
		diff = efficiencyTarget - efficiency
	except:
		pdb.set_trace()

	return diff
# ------------------------------------------------------------------------------------------------ #
