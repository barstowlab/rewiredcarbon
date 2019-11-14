# ------------------------------------------------------------------------------------------------ #
def ImportScenarioTable(scenarioTableFileName):

	import csv
	import ast
	import pdb

	fHandle = open(scenarioTableFileName, 'r')

	poolColumnToHeaderIndexDict = {}
	i = 0
	datareader = csv.reader(fHandle)

	headers = next(datareader, None)
	headerNicknames = next(datareader, None)

	column = {}
	for h in headerNicknames:
		column[h] = []

	for row in datareader:
		for h, v in zip(headerNicknames, row):
			column[h].append(v)

	fHandle.close()

	# Compile everything together into a single line
	scenarioDict = {}
	columnKeys = list(column.keys())

	i = 0
	while i < len(column['Nickname']):
		scenarioDict[column['Nickname'][i]] = {}
	
		j = 0
		scenarioData = {}
		
		while j < len(columnKeys):
			scenarioData[columnKeys[j]] = column[columnKeys[j]][i]
			j += 1
	
		scenarioDict[column['Nickname'][i]] = scenarioData
		i += 1

	scenarioNames = scenarioDict.keys()

	for scenarioName in scenarioNames:
		scenarioData = scenarioDict[scenarioName]
		scenarioDataKeys = scenarioData.keys()
	
		for key in scenarioDataKeys:
			entry = scenarioData[key]
		
			if len(entry) > 0 and entry[0] == '[':
				entry_evaluated = ast.literal_eval(entry)
			else:
				entry_evaluated = entry
			
			if type(entry_evaluated) == list:
				scenarioData[key] = entry_evaluated

	return scenarioDict
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def CalculateScenarioEfficiencies(scenarioDict, mode='bargraph', \
generateZeroErrorsForBargraphMode=False):

	import pdb

	scenarioNames = scenarioDict.keys()
	
	efficienciesDict = {}

	for scenarioName in scenarioNames:
		scenarioData = scenarioDict[scenarioName]
		
		try:
			efficiencyDict = ProcessScenario(scenarioData, mode=mode, \
			generateZeroErrorsForBargraphMode=generateZeroErrorsForBargraphMode)
		except:
			pdb.set_trace()
		
		efficienciesDict[scenarioName] = efficiencyDict
		
	return efficienciesDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ProcessScenario(scenarioData, mode='bargraph', generateZeroErrorsForBargraphMode=False):
	
	import pdb
	import sys
	
	# Figure out if any of the input variables have multiple variables
	scenarioKeys = scenarioData.keys()
	keysWithMultiValueVariable = []
	
	for key in scenarioKeys:
		if type(scenarioData[key]) is list:
			keysWithMultiValueVariable.append(key)
	
	
	if len(keysWithMultiValueVariable) == 1:
		multiValueVariableKey = keysWithMultiValueVariable[0]
		multiValueVariable = scenarioData[multiValueVariableKey]
		
		if mode == 'bargraph':	
			if len(multiValueVariable) != 3:
				print("Currently only dealing with 3 value multi-value variables in bargraph mode.")
			else:
				efficiencyDict = \
				Process_Scenario_and_Calculate_Errors_for_BargraphMode(scenarioData, \
				multiValueVariableKey=multiValueVariableKey, \
				generateZeroErrors=generateZeroErrorsForBargraphMode)
		
		elif mode == 'scattergraph':
			
# 			pdb.set_trace()
			
			efficiencyDict = \
			Process_Scenario_for_ScattergraphMode(scenarioData, multiValueVariableKey)
				
		else:
			print("We only do bargraph and scattergraph modes right now.")
	
	elif len(keysWithMultiValueVariable) == 0:
		if mode == 'bargraph':
			efficiencyDict = \
			Process_Scenario_and_Calculate_Errors_for_BargraphMode(scenarioData, \
			multiValueVariableKey=None, \
			generateZeroErrors=generateZeroErrorsForBargraphMode)
		else:
			print("Currently only dealing with single-value variable scenarios in bargraph mode.")
			pdb.set_trace()
		
	else:
		print("Can't yet deal with more than one multi-value variable per scenario.")
		sys.exit()
	

	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Process_Scenario_and_Calculate_Errors_for_BargraphMode(scenarioData, \
multiValueVariableKey=None, generateZeroErrors=False):
	
	from copy import deepcopy
	
	# If we are dealing with a scenario with no multi-value variables
	if multiValueVariableKey == None and generateZeroErrors == False:
		efficiencyDict = ProcessSingleValueScenario(scenarioData)
		returnValue = efficiencyDict
	
	elif multiValueVariableKey == None and generateZeroErrors == True:
		efficiencyDict = ProcessSingleValueScenario(scenarioData)
		midValue = deepcopy(efficiencyDict)
		lowerValue = deepcopy(efficiencyDict)
		upperValue = deepcopy(efficiencyDict)
		
		efficiencyErrorUpperDict = Calculate_EfficiencyDict_Error(midValue, upperValue, \
		'upperError')
		efficiencyErrorLowerDict = Calculate_EfficiencyDict_Error(midValue, lowerValue, \
		'lowerError')
	
		efficiencyDictWithErrors = \
		Combine_Efficiency_Errors(midValue, efficiencyErrorLowerDict, efficiencyErrorUpperDict)
		
		returnValue = efficiencyDictWithErrors
		
	# On the other hand, if we are dealing with a scenario with a multi-value variable
	elif multiValueVariableKey != None:		
		multiValueVariable = scenarioData[multiValueVariableKey]
		i = 0
		scenarioDataArray = []
		efficiencyDictArray = []
	
		while i < len(multiValueVariable):
			# Get the instance of the multi-value variable
			valueOfMultiValueVariable = multiValueVariable[i]
			scenarioDataSingleValue = deepcopy(scenarioData)
			scenarioDataSingleValue[multiValueVariableKey] = valueOfMultiValueVariable
			scenarioDataArray.append(scenarioDataSingleValue)
			efficiencyDictArray.append(ProcessSingleValueScenario(scenarioDataSingleValue))
			i += 1
	
		efficiencyDictArray = sorted(efficiencyDictArray, \
		key = lambda i: i['effTotalElectricalToFuel'])
	
		lowerValue = efficiencyDictArray[0]
		midValue = efficiencyDictArray[1]
		upperValue = efficiencyDictArray[2]

		efficiencyErrorUpperDict = Calculate_EfficiencyDict_Error(midValue, upperValue, \
		'upperError')
		efficiencyErrorLowerDict = Calculate_EfficiencyDict_Error(midValue, lowerValue, \
		'lowerError')
	
		efficiencyDictWithErrors = \
		Combine_Efficiency_Errors(midValue, efficiencyErrorLowerDict, efficiencyErrorUpperDict)
		
		returnValue = efficiencyDictWithErrors
	
	
	return returnValue
# ------------------------------------------------------------------------------------------------ #
			

# ------------------------------------------------------------------------------------------------ #
def Process_Scenario_for_ScattergraphMode(scenarioData, multiValueVariableKey):
	
	from copy import deepcopy
	from numpy import arange, float, logspace
	import pdb
	
	multiValueVariable = scenarioData[multiValueVariableKey]
	
	try:
		independentVariableScale = scenarioData['independentVariableScale']
	except:
		independentVariableScale = 'Linear'
		
	try:
		dependentVariableScale = scenarioData['dependentVariableScale']
	except:
		dependentVariableScale = 'Linear'
		

	# Generate an array		
	if independentVariableScale.lower() == 'logarithmic':
		multiValueVariableArray = logspace(float(multiValueVariable[0]), \
		float(multiValueVariable[1]), num=float(multiValueVariable[2]))
	elif independentVariableScale.lower() == 'linear':
		multiValueVariableArray = arange(float(multiValueVariable[0]), \
		float(multiValueVariable[1]), float(multiValueVariable[2]))
	else:
		print('We can only do Linear and Logarithmic independent variable scales right now.')
		sys.exit()

# 	rpdb.set_trace()


	scenarioDataArray = []
	efficiencyDictArray = []
	
	scenarioDataMaster = deepcopy(scenarioData)
	
	i = 0 
	while i < len(multiValueVariableArray):
		valueOfMultiValueVariable = multiValueVariableArray[i]
		scenarioDataSingleValue = deepcopy(scenarioData)
		scenarioDataSingleValue[multiValueVariableKey] = valueOfMultiValueVariable
		scenarioDataArray.append(scenarioDataSingleValue)
		efficiencyDictArray.append(ProcessSingleValueScenario(scenarioDataSingleValue))
		i += 1
# 		
	consolidatedEfficiencyDict = Combine_Efficiency_Array(efficiencyDictArray, \
	multiValueVariableArray, multiValueVariableKey)
	
	consolidatedEfficiencyDict['independentVariableScale'] = independentVariableScale
	consolidatedEfficiencyDict['dependentVariableScale'] = dependentVariableScale
	
	return consolidatedEfficiencyDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Combine_Efficiency_Array(efficiencyDictArray, multiValueVariableArray, multiValueVariableKey):
	
	consolidatedEfficiencyDict = {}
	
	consolidatedEfficiencyDict[multiValueVariableKey] = multiValueVariableArray
	consolidatedEfficiencyDict['multiValueVariableKey'] = multiValueVariableKey
	
	
	# Break out the efficiencyDictArray
	efficiencyDictArrayFirstElement = efficiencyDictArray[0]
	efficiencyDictArrayFirstElementKeys = list(efficiencyDictArrayFirstElement.keys())
	
	for key in efficiencyDictArrayFirstElementKeys:
		consolidatedEfficiencyDict[key] = []
	
	i = 0
	while i < len(efficiencyDictArray):
		keys = list(efficiencyDictArray[i].keys())
		
		for key in keys:
			consolidatedEfficiencyDict[key].append(efficiencyDictArray[i][key])
	
		i += 1
	
	return consolidatedEfficiencyDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ProcessSingleValueScenario(scenarioData):
	
	from rewiredcarbon.efficiency import Process_Hydrogen_with_ElectrochemicalCO2_Scenario, \
	Process_EET_with_ElectrochemicalCO2_Scenario, Process_Hydrogen_with_BioCO2_Scenario, \
	Process_EET_with_BioCO2_Scenario
	
	co2Method = scenarioData['CO2Method']
	mediator = scenarioData['Mediator']
	

	if co2Method == 'Electrochemical':
		if mediator == 'H2':
			efficiencyDict = \
			Process_Hydrogen_with_ElectrochemicalCO2_Scenario(scenarioData)			
		elif mediator == 'EET':
			efficiencyDict = \
			Process_EET_with_ElectrochemicalCO2_Scenario(scenarioData)
		else:
			print("We can't do mediators other than H2 or EET, at least not yet.")
			efficiencyDict = {}

	elif co2Method == 'Enzymatic':
		if mediator == 'H2':
			efficiencyDict = Process_Hydrogen_with_BioCO2_Scenario(scenarioData)
		elif mediator == 'EET':			
			efficiencyDict = Process_EET_with_BioCO2_Scenario(scenarioData)
		else:
			print("We can't do mediators other than H2 or EET, at least not yet.")
			efficiencyDict = {}

	else:
		print("We can't do CO2 fixation methods that aren't Electrochemical or Enzymatic, "\
		+ "at least not yet.")
		efficiencyDict = {}

	return efficiencyDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Calculate_EfficiencyDict_Error(refValueDict, extremumValueDict, errorLabel, separator='_'):
	
	import pdb
	
	efficiencyErrorDict = {}
	refValueDictKeys = refValueDict.keys()
	
	for key in refValueDictKeys:
		errorKey = key + separator + errorLabel
		try:
			efficiencyErrorDict[errorKey] = abs(refValueDict[key] - extremumValueDict[key])
		except:
			pdb.set_trace()


	return efficiencyErrorDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Combine_Efficiency_Errors(midValue, lowerValue, upperValue):
	
	from copy import deepcopy
	import pdb
	
	lowerValueKeys = lowerValue.keys()
	upperValueKeys = upperValue.keys()
	
	
	efficiencyDictWithErrors = deepcopy(midValue)
	
	for key in lowerValueKeys:
		efficiencyDictWithErrors[key] = lowerValue[key]
	
	for key in upperValue:
		efficiencyDictWithErrors[key] = upperValue[key]
	
# 	pdb.set_trace()
	
	return efficiencyDictWithErrors
# ------------------------------------------------------------------------------------------------ #








# ------------------------------------------------------------------------------------------------ #
def Plot_Efficiency_Bargraph(efficienciesDict, efficiencyKeyToPlot, efficiencyErrorLowerKeyToPlot, \
efficiencyErrorUpperKeyToPlot, keysToPlot=None):
	
	from matplotlib.pyplot import bar, figure, show, xticks, subplots_adjust
	import pdb

	if keysToPlot == None:
		scenarioKeys = list(efficienciesDict.keys())
	else:
		scenarioKeys = keysToPlot

	efficiencyValuesArray = []
	efficiencyValuesUpperErrorArray = []
	efficiencyValuesLowerErrorArray = []

	for key in scenarioKeys:
		efficiencyValuesArray.append(efficienciesDict[key][efficiencyKeyToPlot])
	
		lowerError = efficienciesDict[key][efficiencyErrorLowerKeyToPlot]
		upperError = efficienciesDict[key][efficiencyErrorUpperKeyToPlot]
		efficiencyValuesLowerErrorArray.append(lowerError)
		efficiencyValuesUpperErrorArray.append(upperError)

	efficiencyValuesErrorArray = [efficiencyValuesLowerErrorArray, efficiencyValuesUpperErrorArray]

	figure()
	bar(scenarioKeys, efficiencyValuesArray, yerr=efficiencyValuesErrorArray, capsize=5)
	xticks(scenarioKeys, scenarioKeys, rotation='vertical')
	subplots_adjust(bottom=0.4)
	show()

	
	return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Plot_Efficiency_Scattergraphs(efficienciesDict, efficiencyKeyToPlot, keysToPlot=None, \
overridePlotScalesInScenarioFile=False, xScale='Linear', yScale='Linear'):
	
	from matplotlib.pyplot import plot, figure, show, xlabel, ylabel, title, loglog, semilogx, \
	semilogy
	import pdb

	if keysToPlot == None:
		scenarioKeys = list(efficienciesDict.keys())
	else:
		scenarioKeys = keysToPlot
		
	for key in scenarioKeys:		
		
		multiValueVariableKey = efficienciesDict[key]['multiValueVariableKey']
		
		if overridePlotScalesInScenarioFile == False:
			dependentVariableScale = efficienciesDict[key]['dependentVariableScale']
			independentVariableScale = efficienciesDict[key]['independentVariableScale']
		else:
			dependentVariableScale = yScale
			independentVariableScale = xScale

		
		xAxis = efficienciesDict[key][multiValueVariableKey]
		yAxis = efficienciesDict[key][efficiencyKeyToPlot]
		
		figure()
		
		if independentVariableScale.lower() == 'logarithmic' \
		and dependentVariableScale.lower() == 'logarithmic':
			loglog(xAxis, yAxis)
		elif independentVariableScale.lower() == 'logarithmic' \
		and dependentVariableScale.lower() == 'linear':
			semilogx(xAxis, yAxis)
		elif independentVariableScale.lower() == 'linear' \
		and dependentVariableScale.lower() == 'logarithmic':
			semilogy(xAxis, yAxis)
		elif independentVariableScale.lower() == 'linear' \
		and dependentVariableScale.lower() == 'linear':
			plot(xAxis, yAxis)
		
		title(key)
		xlabel(multiValueVariableKey)
		ylabel(efficiencyKeyToPlot)
		
		show()
	
	return
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Plot_Efficiency_Scattergraphs_2(efficienciesDict, xAxisKey, efficiencyKeyToPlot, \
keysToPlot=None, overridePlotScalesInScenarioFile=False, xScale='Linear', yScale='Linear'):
	
	# Like Plot_Efficiency_Scattergraphs, but allows the user to choose the x variable they want
	# to plot as well. 
	
	from matplotlib.pyplot import plot, figure, show, xlabel, ylabel, title, loglog, semilogx, \
	semilogy
	import pdb

	if keysToPlot == None:
		scenarioKeys = list(efficienciesDict.keys())
	else:
		scenarioKeys = keysToPlot
		
	for key in scenarioKeys:		
		
		if overridePlotScalesInScenarioFile == False:
			dependentVariableScale = efficienciesDict[key]['dependentVariableScale']
			independentVariableScale = efficienciesDict[key]['independentVariableScale']
		else:
			dependentVariableScale = yScale
			independentVariableScale = xScale

		
		xAxis = efficienciesDict[key][xAxisKey]
		yAxis = efficienciesDict[key][efficiencyKeyToPlot]
		
		figure()
		
		if independentVariableScale.lower() == 'logarithmic' \
		and dependentVariableScale.lower() == 'logarithmic':
			loglog(xAxis, yAxis)
		elif independentVariableScale.lower() == 'logarithmic' \
		and dependentVariableScale.lower() == 'linear':
			semilogx(xAxis, yAxis)
		elif independentVariableScale.lower() == 'linear' \
		and dependentVariableScale.lower() == 'logarithmic':
			semilogy(xAxis, yAxis)
		elif independentVariableScale.lower() == 'linear' \
		and dependentVariableScale.lower() == 'linear':
			plot(xAxis, yAxis)
		
		title(key)
		xlabel(xAxisKey)
		ylabel(efficiencyKeyToPlot)
		
		show()
	
	return
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Generate_EfficienciesDict_Keys_Sorted_by_Efficiency(efficienciesDict, efficiencyKeyToSortBy):
	
	scenarioKeys = list(efficienciesDict.keys())
	efficiencyValuesArray = []
	
	for key in scenarioKeys:
		efficiencyValuesArray.append([key, efficienciesDict[key][efficiencyKeyToSortBy]])
		
	efficiencyValuesArray = sorted(efficiencyValuesArray, key = lambda i: i[1])
	
	keysArray = []
	i = 0
	while i < len(efficiencyValuesArray):
		keysArray.append(efficiencyValuesArray[i][0])
		i += 1
	
	
	return keysArray
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Export_Efficiency_Bargraph(filename, efficienciesDict, scenarioDict, efficiencyKeyToPlot, \
efficiencyErrorLowerKeyToPlot, efficiencyErrorUpperKeyToPlot, keysToPlot=None):
	
	from rewiredcarbon.vectorOutput import generateOutputMatrixWithHeaders, writeOutputMatrix
	import pdb

	if keysToPlot == None:
		scenarioKeys = list(efficienciesDict.keys())
	else:
		scenarioKeys = keysToPlot

	efficiencyValuesArray = []
	efficiencyValuesUpperErrorArray = []
	efficiencyValuesLowerErrorArray = []
	labelsArray = []

	for key in scenarioKeys:
		efficiencyValuesArray.append(efficienciesDict[key][efficiencyKeyToPlot])
	
		lowerError = efficienciesDict[key][efficiencyErrorLowerKeyToPlot]
		upperError = efficienciesDict[key][efficiencyErrorUpperKeyToPlot]
		efficiencyValuesLowerErrorArray.append(lowerError)
		efficiencyValuesUpperErrorArray.append(upperError)
		
		labelsArray.append('"' + scenarioDict[key]['BargraphLabel'] + '"')
		
	headers = ['scenario', 'label', efficiencyKeyToPlot, efficiencyErrorLowerKeyToPlot, \
	efficiencyErrorUpperKeyToPlot]

	vectorList = [scenarioKeys, labelsArray, efficiencyValuesArray, \
	efficiencyValuesLowerErrorArray, efficiencyValuesUpperErrorArray]

	oMatrix = generateOutputMatrixWithHeaders(vectorList, headers, delimeter=',')
	
	writeOutputMatrix(filename, oMatrix)

	return
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Export_Efficiency_Scattergraphs(outputDirname, outputFilePrefix, efficienciesDict, \
efficiencyKeyToPlot, keysToPlot=None, addKeyToHeader=False):
	
	from rewiredcarbon.vectorOutput import generateOutputMatrixWithHeaders, writeOutputMatrix
	import pdb

	if keysToPlot == None:
		scenarioKeys = list(efficienciesDict.keys())
	else:
		scenarioKeys = keysToPlot

	for key in scenarioKeys:		
		
		multiValueVariableKey = efficienciesDict[key]['multiValueVariableKey']
		
		xAxis = efficienciesDict[key][multiValueVariableKey]
		yAxis = efficienciesDict[key][efficiencyKeyToPlot]
		
		if addKeyToHeader == False:
			headers = [multiValueVariableKey, efficiencyKeyToPlot]
		elif addKeyToHeader == True:
			headers = [multiValueVariableKey + '_' + key, efficiencyKeyToPlot + '_' + key]

		vectorList = [xAxis, yAxis]

		oMatrix = generateOutputMatrixWithHeaders(vectorList, headers, delimeter=',')
		
		
		if outputDirname == '':
			outputFilename = outputFilePrefix + '_' + key + '.csv'
		else:
			outputFilename = outputDirname + '/' + outputFilePrefix + '_' + key + '.csv'
	
		writeOutputMatrix(outputFilename, oMatrix)


	return
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def Export_Efficiency_Scattergraphs_with_More_than_1_Multi_Value_Variable(outputDirname, \
outputFilePrefix, efficienciesDict, efficiencyKeyToPlot, keysToPlot=None, addKeyToHeader=False):
	
	from rewiredcarbon.vectorOutput import generateOutputMatrixWithHeaders, writeOutputMatrix
	import pdb

	if keysToPlot == None:
		scenarioKeys = list(efficienciesDict.keys())
	else:
		scenarioKeys = keysToPlot
		
	for key in scenarioKeys:		
		
		multiValueVariableKeys = efficienciesDict[key]['multiValueVariableKeys']
		
		for multiValueVariableKey in multiValueVariableKeys:
		
			xAxis = efficienciesDict[key][multiValueVariableKey]
			yAxis = efficienciesDict[key][efficiencyKeyToPlot]
			
			if addKeyToHeader == False:
				headers = [multiValueVariableKey, efficiencyKeyToPlot]
			elif addKeyToHeader == True:
				headers = [multiValueVariableKey + '_' + key, efficiencyKeyToPlot + '_' + key]

			
			vectorList = [xAxis, yAxis]
			
			oMatrix = generateOutputMatrixWithHeaders(vectorList, headers, delimeter=',')
		
				
			if outputDirname == '':
				outputFilename = outputFilePrefix + '_' + key + '.csv'
			else:
				outputFilename = outputDirname + '/' + outputFilePrefix + '_' \
				+ multiValueVariableKey + '_' + key + '.csv'
	
			writeOutputMatrix(outputFilename, oMatrix)
		
	
	return
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def Plot_Efficiency_Scattergraphs_with_More_than_1_Multi_Value_Variable(efficienciesDict, \
efficiencyKeyToPlot, keysToPlot=None, \
overridePlotScalesInScenarioFile=False, xScale='Linear', yScale='Linear'):
	
	from matplotlib.pyplot import plot, figure, show, xlabel, ylabel, title, loglog, semilogx, \
	semilogy
	import pdb

	if keysToPlot == None:
		scenarioKeys = list(efficienciesDict.keys())
	else:
		scenarioKeys = keysToPlot
		
	for key in scenarioKeys:		
		
		multiValueVariableKeys = efficienciesDict[key]['multiValueVariableKeys']
		
		if overridePlotScalesInScenarioFile == False:
			dependentVariableScale = efficienciesDict[key]['dependentVariableScale']
			independentVariableScale = efficienciesDict[key]['independentVariableScale']
		else:
			dependentVariableScale = yScale
			independentVariableScale = xScale

		
		for multiValueVariableKey in multiValueVariableKeys:
		
			xAxis = efficienciesDict[key][multiValueVariableKey]
			yAxis = efficienciesDict[key][efficiencyKeyToPlot]
		
			figure()
		
			if independentVariableScale.lower() == 'logarithmic' \
			and dependentVariableScale.lower() == 'logarithmic':
				loglog(xAxis, yAxis)
			elif independentVariableScale.lower() == 'logarithmic' \
			and dependentVariableScale.lower() == 'linear':
				semilogx(xAxis, yAxis)
			elif independentVariableScale.lower() == 'linear' \
			and dependentVariableScale.lower() == 'logarithmic':
				semilogy(xAxis, yAxis)
			elif independentVariableScale.lower() == 'linear' \
			and dependentVariableScale.lower() == 'linear':
				plot(xAxis, yAxis)
		
			title(key)
			xlabel(multiValueVariableKey)
			ylabel(efficiencyKeyToPlot)
		
			show()
	
	return
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Combine_Efficiency_Array_with_More_than_1_Multi_Value_Variable(efficiencyDictArray, \
multiValueVariableDict):
	
	import pdb
	
	consolidatedEfficiencyDict = {}
	
	multiValueVariableKeys = list(multiValueVariableDict.keys())
	
	consolidatedEfficiencyDict['multiValueVariableKeys'] = multiValueVariableKeys
	
	i = 0
	while i < len(multiValueVariableKeys):
		
		multiValueVariableKey = multiValueVariableKeys[i]
		multiValueVariableArray = multiValueVariableDict[multiValueVariableKey]
		
		consolidatedEfficiencyDict[multiValueVariableKey] = multiValueVariableArray
		
		i += 1
	
	
	# Break out the efficiencyDictArray
	efficiencyDictArrayFirstElement = efficiencyDictArray[0]
	efficiencyDictArrayFirstElementKeys = list(efficiencyDictArrayFirstElement.keys())
	
	for key in efficiencyDictArrayFirstElementKeys:
		consolidatedEfficiencyDict[key] = []
	
	i = 0
	while i < len(efficiencyDictArray):
		keys = list(efficiencyDictArray[i].keys())
		
		for key in keys:
			consolidatedEfficiencyDict[key].append(efficiencyDictArray[i][key])
	
		i += 1
	
	return consolidatedEfficiencyDict
# ------------------------------------------------------------------------------------------------ #

