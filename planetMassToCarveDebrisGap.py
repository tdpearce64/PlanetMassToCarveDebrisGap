
'''Python program to calculate the gap carved by a planet embedded in a
massive debris disc, by Tim D. Pearce. The model is described in 
Friebe, Pearce & Lohne 2022, and produces a figure akin to their Fig. 5.
The main difference to conventional gap width estimates (which use e.g. 
the Wisdom overlap criterion) is that, if the debris disc has mass, then
an embedded planet would migrate as it scattered debris. This means that 
a low-mass, strongly migrating planet could carve as wide a gap as a 
high-mass, barely migrating planet. Given the parameters of the star, the
disc edges, the central location of a gap, and the (pre-interaction) 
debris disc mass, the program plots gap width as a function of the mass 
of an embedded planet. Additionally, if a target gap width is provided by
the user (optional), then the program returns the parameters of the 
possible planets that could carve that gap. 

IMPORTANT: The model assumes that the interaction has finished, i.e. the 
planet has stopped migrating. If instead migration is still ongoing, then
the disc gap would widen in time and these results could be invalid (see 
Sect. 3.3 in Friebe, Pearce & Lohne 2022 for ways to check and account 
for this).

To use the program, simply change the values in the 'User Inputs' section
just below. You should not have to change anything outside of that 
section. The default settings are for HD107146.

Feel free to use this code, and if the results go into a publication,
then please cite Friebe, Pearce & Lohne 2022. Finally, let me know if you
find any bugs or have any requests!'''

############################### Libraries ###############################
import sys
import numpy as np
import math
from scipy.optimize import minimize
import matplotlib.pyplot as plt

############################## User Inputs ##############################
'''Parameters of the system. You should not need to change any part of 
the code other than that in this section. The units of each parameter are
given at the end of thier name, e.g. mStar_mSun is the star mass in Solar 
masses.'''

# Star mass
mStar_mSun = 1.0

# Disc values (NOTE: mDisc_mEarth is the mass of the pre-interaction 
# disc, before material has been removed by planet scattering)
mDisc_mEarth = 50.
discInnerEdge_au = 36.6		
discOuterEdge_au = 145.6
gapCentralRadius_au = 75.5

# Initial disc surface density index: SD ~ r^{-sdIndex}. Note: strictly
# cannot be exactly 1 or 2. If 1 or 2 asked for, then the code will 
# change it to e.g. 1.000001 instead
sdIndex = 0.25

# Number of Hill radii that planet clears either side of it (k in 
# Friebe, Pearce & Lohne 2022)
hillRadiiClearedByPlt = 3

# Planet masses to consider
minPltMass_mJup = 1e-3
maxPltMass_mJup = 1e2
numberOfPltMasses = 300

# Target gap width to find planet parameters for (optional: enter None if
# not required)
targetGapWidth_au = 38.6

############################### Constants ###############################
# Common functions of the surface density index
if sdIndex == 1: sdIndex = 1.000001
if sdIndex == 2: sdIndex = 2.000001
twoMinusY = 2.-sdIndex
YMinusOne = sdIndex-1.

# Mass conversions
mSun_kg = 1.988435e30
mJup_kg = 1.898e27
mEarth_kg = 5.97e24

mStar_mJup = mStar_mSun * mSun_kg / mJup_kg
mDisc_mJup = mDisc_mEarth * mEarth_kg / mJup_kg

# Log planet masses
log10MinPltMass_mJup = np.log10(minPltMass_mJup)
log10MaxPltMass_mJup = np.log10(maxPltMass_mJup)

# The maximum possible gap width (gap must not go outside disc edges)
maxGapWidth_au = 2.*min(gapCentralRadius_au-discInnerEdge_au, discOuterEdge_au-gapCentralRadius_au)

# Edges of the target gap
if targetGapWidth_au is not None:
	targetGapOuterEdge_au = gapCentralRadius_au + targetGapWidth_au/2.
	targetGapInnerEdge_au = gapCentralRadius_au - targetGapWidth_au/2.

############################ Maths functions ############################
def RoundNumberToDesiredSigFigs(num, sigFigs=None):
	'''Returns a float of a given number to the specified significant 
	figures'''

	# Default number of decimal digits (in case precision unspecified)
	defaultSigFigs = 2

	# Get sig figs if not defined
	if sigFigs is None:
		sigFigs = defaultSigFigs
			
	# Catch case if number is zero
	if num == 0:
		exponent = 0
	
	# Otherwise number is non-zero
	else:
		exponent = GetBase10OrderOfNumber(num)
	
	# Get the coefficient
	coefficient = round(num / float(10**exponent), sigFigs-1)

	roundedNumber = coefficient * 10**exponent
	
	# Get the decimal places to round to (to avoid computer rounding errors)
	decimalPlacesToRoundTo = sigFigs-exponent
	
	roundedNumber = round(roundedNumber, decimalPlacesToRoundTo)
	
	return roundedNumber
	
#------------------------------------------------------------------------
def GetBase10OrderOfNumber(number):
	'''Return the order of the positive number in base 10, e.g. inputting
	73 returns 1 (since 10^1 < 73 < 10^2).'''
	
	if number <= 0: return np.nan
	
	return int(math.floor(np.log10(abs(number))))

############################ Print functions ############################
def PrintUserInputs():
	'''Print the user inputs'''
	
	print('User inputs:')
	print('     Star mass: %s MSun' % GetValueString(mStar_mSun))
	print('     Pre-interaction disc mass: %s MEarth' % GetValueString(mDisc_mEarth))
	print('     Disc inner edge: %s au' % GetValueString(discInnerEdge_au))
	print('     Disc outer edge: %s au' % GetValueString(discOuterEdge_au))
	print('     Gap central radius: %s au' % GetValueString(gapCentralRadius_au))	
	print('     Disc surface density: ~r^(%s)' % GetValueString(sdIndex))
	print('     Hill radii cleared by planet (k): %s' % GetValueString(hillRadiiClearedByPlt))

	if targetGapWidth_au is not None:
		print('     Target gap width: %s au' % GetValueString(targetGapWidth_au))
	else:
		print('     No target gap width entered')

	PrintEmptyLine()

	print('Considering %s planet masses between %s and %s MJup' % (numberOfPltMasses, GetValueString(minPltMass_mJup), GetValueString(maxPltMass_mJup)))
			
	PrintEmptyLine()
	
#------------------------------------------------------------------------	
def CheckUserInputsOK():
	'''Check the user inputs are OK. All values should be positive, and 
	all parameters except sdIndex non-zero'''
	
	areUserInputsOK = True
	reasonsInputsAreBad = []
	
	# Some values must be positive and non-zero
	for value in [mStar_mSun, mDisc_mEarth, discInnerEdge_au, discOuterEdge_au, gapCentralRadius_au, hillRadiiClearedByPlt, minPltMass_mJup, maxPltMass_mJup, numberOfPltMasses]:
		if math.isnan(value) or value <= 0:
			areUserInputsOK = False
			reasonsInputsAreBad.append('Some parameters must be non-zero, positive, and not nan.')
			break

	# Some values can be negative, but must not be nan
	for value in [sdIndex]:
		if math.isnan(value):
			areUserInputsOK = False
			reasonsInputsAreBad.append('A value is nan that cannot be.')
			break	

	# Some values can be None, but otherwise must be positive and non-zero
	for value in [targetGapWidth_au]:
		if value is not None:
			if math.isnan(value) or value <= 0:
				areUserInputsOK = False
				reasonsInputsAreBad.append('Some parameters can be None, but but must otherwise be non-zero, positive, and not nan.')
				break	

	# Disc apo and peri must be defined correctly
	if discOuterEdge_au < discInnerEdge_au:
		reasonsInputsAreBad.append('Disc inner edge apocentre is greater than outer edge pericentre.')
		areUserInputsOK = False

	# Gap central radius must be within disc
	if gapCentralRadius_au <= discInnerEdge_au or gapCentralRadius_au >= discOuterEdge_au:
		reasonsInputsAreBad.append('Gap central radius is not within disc.')
		areUserInputsOK = False		

	# Gap edges must be within disc
	if targetGapWidth_au is not None:
		if gapCentralRadius_au-targetGapWidth_au/2. <= discInnerEdge_au or gapCentralRadius_au+targetGapWidth_au/2. >= discOuterEdge_au:
			reasonsInputsAreBad.append('At least one edge of the target gap is outside the disc.')
			areUserInputsOK = False		

	# Must search real range of planet masses
	if minPltMass_mJup >= maxPltMass_mJup:
		reasonsInputsAreBad.append('Minimum planet mass to consider must be smaller than maximum to consider.')
		areUserInputsOK = False		
				
	# Print errors if the inputs are bad
	if areUserInputsOK == False:
		print('***ERROR*** Problem(s) with user inputs:')
		for reasonInputsAreBad in reasonsInputsAreBad:
			print('     -%s' % reasonInputsAreBad)
		PrintEmptyLine()	

	# If no error, there may still be warnings
	if areUserInputsOK == True:
	
		areUserInputsWarningFree = True
		
		# Should search a reasonable number of planet masses
		orderOfMinPltMass = GetBase10OrderOfNumber(minPltMass_mJup)
		orderOfMaxPltMass = GetBase10OrderOfNumber(maxPltMass_mJup)

		reasonableNumberOfPltMassesToTry = (orderOfMaxPltMass-orderOfMinPltMass)*3
		if numberOfPltMasses < reasonableNumberOfPltMassesToTry:
			reasonsInputsAreBad.append('The number of tested planet masses may be too small. Consider increasing this.')
			areUserInputsWarningFree = False
			
		# Print warnings if the inputs are bad
		if areUserInputsWarningFree == False:
			print('***WARNING*** Possible problem(s) with user inputs:')
			for reasonInputsAreBad in reasonsInputsAreBad:
				print('     -%s' % reasonInputsAreBad)
			PrintEmptyLine()			
					
	return areUserInputsOK

#------------------------------------------------------------------------
def PrintEmptyLine():
	'''Print an empty line (done this way to enable compatibility across 
	Python versions)'''
	
	if sys.version_info[0] < 3: print
	else: print()

#------------------------------------------------------------------------
def GetValueString(value):
	'''Get a string neatly showing the value'''
	
	# Significant figures to round number to
	sigFigsToRoundValueTo = 3
	
	# The orders of the largest and smallest values that should be 
	# written in non-SI notation
	minOrderForNonSINotation = -3
	maxOrderForNonSINotation = 3
		
	# If the value is non-zero and not nan, proceed. Some of the 
	# functions below would otherwise fail
	if math.isnan(value) == False and value > 0:

		# Get the order of the value
		orderOfValue = GetBase10OrderOfNumber(value)

		# Round the value
		valueRounded = RoundNumberToDesiredSigFigs(value, sigFigsToRoundValueTo)
		orderOfRoundedValue = GetBase10OrderOfNumber(valueRounded)

		# If the rounded value has gone up an order, will round value 
		# to an extra significant figure (e.g. 0.099 +/- 0.03 -> 0.10 +/- 0.03)
		if orderOfRoundedValue > orderOfValue:
			sigFigsToRoundValueTo += 1
		
		# If the value is very small or large, divide it by its order, 
		# and later quote order in string. Use rounding function to 
		# remove rounding errors, and note that uncertainties are 
		# always quoted to 1 sig fig
		if orderOfRoundedValue > maxOrderForNonSINotation or orderOfRoundedValue < minOrderForNonSINotation:
			wasPowerAdjustmentDone = True
			
			#valueRounded /= 10**orderOfRoundedValue
			valueRounded = RoundNumberToDesiredSigFigs(valueRounded/10**orderOfRoundedValue, sigFigsToRoundValueTo)

			orderOfRoundedValueAfterPowerAdjust = GetBase10OrderOfNumber(valueRounded)
			
		else:
			wasPowerAdjustmentDone = False
			orderOfRoundedValueAfterPowerAdjust = orderOfRoundedValue

		# If all significant figures of the value are to the left of 
		# the decimal point, the value is an integer
		numberOfValueFiguresLeftOfDecimalPoint = max(orderOfRoundedValueAfterPowerAdjust + 1, 0)

		if numberOfValueFiguresLeftOfDecimalPoint - sigFigsToRoundValueTo >= 0:
			valueRounded = int(valueRounded)

		# Convert the value to a string. If there are value figures to 
		# the right of the decimal point, and the final one(s) should 
		# be zero, append zeros to the string 
		valueRoundedString = str(valueRounded)
		
		numberOfValueFiguresNeededRightOfDecimalPoint = max(sigFigsToRoundValueTo - (orderOfRoundedValueAfterPowerAdjust+1), 0)
		
		if numberOfValueFiguresNeededRightOfDecimalPoint > 0:

			while True:
				indexOfPointInString = valueRoundedString.index('.')
				numberOfValueFiguresRightOfDecimalPointInString = len(valueRoundedString[indexOfPointInString+1:])
									
				if numberOfValueFiguresRightOfDecimalPointInString == numberOfValueFiguresNeededRightOfDecimalPoint:
					break
				
				valueRoundedString += '0'				
		
		# If a power adjustment was done
		if wasPowerAdjustmentDone:
			valueString = '%s * 10^%s' % (valueRoundedString, orderOfRoundedValue)

		# Otherwise no power adjusment was done
		else:
			valueString = '%s' % (valueRoundedString)

	# Otherwise the value is zero or nan
	else:
		valueString = '%s' % (value)

	return valueString

#------------------------------------------------------------------------
def PrintGapCarvingPlanets(pltMass1_mJup, startLocOfPlt1_au, endLocOfPlt1_au, pltMass2_mJup, startLocOfPlt2_au, endLocOfPlt2_au):
	'''Print the parameters of the possible gap-carving planets'''
	
	print('RESULT: Two different planets could carve the gap, if the disc has initial mass %s MEarth:' % GetValueString(mDisc_mEarth))
	PrintEmptyLine()
	
	print('     First planet:')
	print('          Mass: %s MJup' % GetValueString(pltMass1_mJup))
	print('          Pre-migration semimajor axis: %s au' % GetValueString(startLocOfPlt1_au))	
	print('          Post-migration semimajor axis: %s au' % GetValueString(endLocOfPlt1_au))
	PrintEmptyLine()
		
	print('     Second planet:')
	print('          Mass: %s MJup' % GetValueString(pltMass2_mJup))
	print('          Pre-migration semimajor axis: %s au' % GetValueString(startLocOfPlt2_au))	
	print('          Post-migration semimajor axis: %s au' % GetValueString(endLocOfPlt2_au))
	PrintEmptyLine()	
	
########################## Dynamics functions ###########################
def GetGapWidthsAtPlanetMasses_au(pltMasses_mJup):
	'''Solve Friebe2022 eq. 6 for the gap width as a function of planet 
	mass. This is the exact equation, and is solved numerically'''
	
	print('Getting exact gap width at each planet mass...')
	
	gapWidths_au = []
	
	for pltMass_mJup in pltMasses_mJup:
		gapWidth_au = GetGapWidthAtPlanetMass_au(pltMass_mJup)
		gapWidths_au.append(gapWidth_au)
	
	print('Complete.')
	PrintEmptyLine()
	
	return gapWidths_au
		
#------------------------------------------------------------------------
def GetGapWidthAtPlanetMass_au(pltMass_mJup):
	''''''

	h = GetHillRadius_PlanetSemimajorAxis(pltMass_mJup)

	# Get the inital guess. Can't be exactly the no-migration gap width,
	# because this gives an infinite value in the solver
	initialGapWidth_au = 2*GetGapWidthIfNoMigrationAtPlanetMass_au(pltMass_mJup, h=h)

	# Numerically solve for the gap width
	minimisationFunction = minimize(FunctionToMinimise, initialGapWidth_au, args=(pltMass_mJup,h), bounds=[(1e-6,maxGapWidth_au)])
	gapWidth_au = minimisationFunction.x[0]
	
	# If the fit is bad, return nan
	gapInnerEdge_au, gapOuterEdge_au, Gamma_aplt0_perAu = GetGapParsForGivenWidth(gapWidth_au)
	planetMassToCarveGap_mJup = GetPlanetMassToCarveGap_mJup(gapWidth_au, h, gapInnerEdge_au, gapOuterEdge_au, Gamma_aplt0_perAu)
	if abs(planetMassToCarveGap_mJup/pltMass_mJup - 1) > 1e-2: return np.nan
	
	return gapWidth_au

#------------------------------------------------------------------------
def FunctionToMinimise(gapWidth_au, pltMass_mJup, h):
	
	gapInnerEdge_au, gapOuterEdge_au, Gamma_aplt0_perAu = GetGapParsForGivenWidth(gapWidth_au)
	
	planetMassToCarveGap_mJup = GetPlanetMassToCarveGap_mJup(gapWidth_au, h, gapInnerEdge_au, gapOuterEdge_au, Gamma_aplt0_perAu)
	
	largerPlanetMass_mJup = max(pltMass_mJup, planetMassToCarveGap_mJup)
	smallerPlanetMass_mJup = min(pltMass_mJup, planetMassToCarveGap_mJup)
	functionToMinimise = abs(largerPlanetMass_mJup/smallerPlanetMass_mJup - 1)
	
	return functionToMinimise

#------------------------------------------------------------------------
def GetHillRadius_PlanetSemimajorAxis(pltMass_mJup):

	hillRadius_aPlt = (pltMass_mJup/(3.*mStar_mJup))**(1./3.)
	
	return hillRadius_aPlt
	
#------------------------------------------------------------------------
def GetGapWidthIfNoMigrationAtPlanetMass_au(pltMass_mJup, h=None):
	''''''

	if h is None:
		h = GetHillRadius_PlanetSemimajorAxis(pltMass_mJup)
	
	gapWidth_au = 2.*hillRadiiClearedByPlt*gapCentralRadius_au*h

	if gapWidth_au > maxGapWidth_au: return np.nan
	
	return gapWidth_au
	
#------------------------------------------------------------------------
def GetGapParsForGivenWidth(gapWidth_au):

	# Get the inner and outer edges of the gap
	gapInnerEdge_au = gapCentralRadius_au - gapWidth_au/2.
	gapOuterEdge_au = gapCentralRadius_au + gapWidth_au/2.

	# Catch numeric error resultin from approximate gap width equation 	
	if gapInnerEdge_au < 0 or gapOuterEdge_au < 0:
		return gapInnerEdge_au, gapOuterEdge_au, np.nan
	
	# Get Gamma divided by the initial planet semimajor axis
	Gamma_aplt0_perAu = twoMinusY/YMinusOne * (gapInnerEdge_au**-YMinusOne-gapOuterEdge_au**-YMinusOne)/(discOuterEdge_au**twoMinusY-discInnerEdge_au**twoMinusY)

	return gapInnerEdge_au, gapOuterEdge_au, Gamma_aplt0_perAu

#------------------------------------------------------------------------
def GetPlanetMassToCarveGap_mJup(gapWidth_au, h, gapInnerEdge_au, gapOuterEdge_au, Gamma_aplt0_perAu):
	
	planetMassToCarveGap_mJup = mDisc_mJup*Gamma_aplt0_perAu*gapInnerEdge_au*gapOuterEdge_au/(gapWidth_au-2*hillRadiiClearedByPlt*gapCentralRadius_au*h)
	
	return planetMassToCarveGap_mJup
	
#------------------------------------------------------------------------
def GetApproxGapWidths_au(pltMasses_mJup):
	'''Solve Friebe2022 eq. 7 for the gap width as a function of planet
	mass. This is an approximate equation that does not need to be solved 
	numerically'''

	print('Getting approximate gap width at each planet mass...')
		
	approxGapWidths_au = []
		
	for pltMass_mJup in pltMasses_mJup:
		approxGapWidth_au = GetApproxGapWidth_au(pltMass_mJup)
		approxGapWidths_au.append(approxGapWidth_au)

	print('Complete.')
	PrintEmptyLine()
	
	return approxGapWidths_au	
	
#------------------------------------------------------------------------
def GetApproxGapWidth_au(pltMass_mJup, h=None):
	'''Approximate form of the gap width that doesn't need to be 
	evaluated numerically. Valid if gap width << 2*gap radius'''
	
	if h is None:
		h = GetHillRadius_PlanetSemimajorAxis(pltMass_mJup)
	
	approxGapWidth_au = 2.*gapCentralRadius_au*hillRadiiClearedByPlt*h / (1.-mDisc_mJup/pltMass_mJup*twoMinusY*gapCentralRadius_au**twoMinusY/(discOuterEdge_au**twoMinusY - discInnerEdge_au**twoMinusY))

	# Approx gap width must be positive
	if approxGapWidth_au <0: return np.nan
	
	# Only valid if gap width does not take gap outside disc
	gapInnerEdge_au, gapOuterEdge_au, Gamma_aplt0_perAu = GetGapParsForGivenWidth(approxGapWidth_au)

	if gapInnerEdge_au < discInnerEdge_au or gapOuterEdge_au > discOuterEdge_au or approxGapWidth_au <0:
		return np.nan
	
	return approxGapWidth_au	

#------------------------------------------------------------------------
def InsertSpecialPointsIntoPlanetMassList(pltMasses_mJup, gapWidths_au, approxGapWidths_au):

	# Make dictionaries
	gapWidthsByPlanetMass_au = {}
	gapWidthsByPlanetMassApprox_au = {}
	
	for i in range(len(pltMasses_mJup)):
		pltMass_mJup = pltMasses_mJup[i]
		gapWidth_au = gapWidths_au[i]
		approxGapWidth_au = approxGapWidths_au[i]
	
		gapWidthsByPlanetMass_au[pltMass_mJup] = gapWidth_au
		gapWidthsByPlanetMassApprox_au[pltMass_mJup] = approxGapWidth_au
	
	# Get min mass that doesn't result in nan widths
	minPltMassNotNanWidthExact_mJup = GetSmallestPlanetMassGivingNonNanWidth_mJup(pltMasses_mJup, gapWidths_au)
	minPltMassNotNanWidthApprox_mJup = GetSmallestPlanetMassGivingNonNanWidth_mJup(pltMasses_mJup, approxGapWidths_au)
	
	maxWidth_au = 2*min(gapCentralRadius_au-discInnerEdge_au, discOuterEdge_au-gapCentralRadius_au)
	
	# Exact equation	
	RPlusWOver2_au = gapCentralRadius_au + maxWidth_au/2.
	RMinusWOver2_au = gapCentralRadius_au - maxWidth_au/2.

	A = 2.*gapCentralRadius_au*hillRadiiClearedByPlt/(3*mStar_mJup)**(1./3.)
	B_exact = mDisc_mJup*twoMinusY/YMinusOne*(-RMinusWOver2_au*RPlusWOver2_au**twoMinusY + RPlusWOver2_au*RMinusWOver2_au**twoMinusY)/(discOuterEdge_au**twoMinusY-discInnerEdge_au**twoMinusY)
	
	minFunExact = lambda mPlt_mJup: abs(A*mPlt_mJup**(1./3.)/(maxWidth_au - B_exact/mPlt_mJup)-1.)

	if minPltMass_mJup != minPltMassNotNanWidthExact_mJup:
		mPltAtMaxWidthExact_mJup = minimize(minFunExact, minPltMassNotNanWidthExact_mJup, bounds=[(minPltMass_mJup,minPltMassNotNanWidthExact_mJup)]).x[0]
	else:
		mPltAtMaxWidthExact_mJup = minPltMassNotNanWidthExact_mJup
	gapWidthsByPlanetMass_au[mPltAtMaxWidthExact_mJup] = maxWidth_au
		
	# Approx equation
	B_approx = mDisc_mJup*twoMinusY*gapCentralRadius_au**twoMinusY/(discOuterEdge_au**twoMinusY-discInnerEdge_au**twoMinusY)
	
	minFunApprox = lambda mPlt_mJup: abs(A*mPlt_mJup**(1./3.)/(maxWidth_au*(1.-B_approx/mPlt_mJup))-1.)
	mPltAtMaxWidthApprox_mJup = minimize(minFunApprox, minPltMassNotNanWidthApprox_mJup, bounds=[(minPltMass_mJup,minPltMassNotNanWidthApprox_mJup)]).x[0]
	gapWidthsByPlanetMassApprox_au[mPltAtMaxWidthApprox_mJup] = maxWidth_au
	
	# Convert to lists to plot
	pltMasses_mJup = list(pltMasses_mJup)
	pltMasses_mJup.append(mPltAtMaxWidthExact_mJup)
	pltMasses_mJup.append(mPltAtMaxWidthApprox_mJup)
	pltMasses_mJup.sort()
	
	gapWidths_au, approxGapWidths_au = [], []
	
	for pltMass_mJup in pltMasses_mJup:
		if pltMass_mJup in gapWidthsByPlanetMass_au:
			gapWidths_au.append(gapWidthsByPlanetMass_au[pltMass_mJup])
		else: gapWidths_au.append(GetGapWidthAtPlanetMass_au(pltMass_mJup))
		
		if pltMass_mJup in gapWidthsByPlanetMassApprox_au:
			approxGapWidths_au.append(gapWidthsByPlanetMassApprox_au[pltMass_mJup])
		else: approxGapWidths_au.append(GetApproxGapWidth_au(pltMass_mJup))
		
	return pltMasses_mJup, gapWidths_au, approxGapWidths_au
	
#------------------------------------------------------------------------
def GetSmallestPlanetMassGivingNonNanWidth_mJup(pltMasses_mJup, gapWidths_au):
	'''Get the planet carving the maximum gap width in the low mass limit'''

	minPltMassNotNanWidth_mJup = 1e99
	for i in range(len(pltMasses_mJup)):
		pltMass_mJup = pltMasses_mJup[i]
		gapWidth_au = gapWidths_au[i]
		if math.isnan(gapWidth_au) == False and pltMass_mJup < minPltMassNotNanWidth_mJup:
			minPltMassNotNanWidth_mJup = pltMass_mJup
	
	return minPltMassNotNanWidth_mJup
		
#------------------------------------------------------------------------
def GetGapWidthsIfNoMigrationAtPlanetMasses_au(pltMasses_mJup):
	''''''
	
	gapWidthsIfNoMigration_au = []
	
	for pltMass_mJup in pltMasses_mJup:
		gapWidthIfNoMigration_au = GetGapWidthIfNoMigrationAtPlanetMass_au(pltMass_mJup)
		gapWidthsIfNoMigration_au.append(gapWidthIfNoMigration_au)
		
	return gapWidthsIfNoMigration_au
	
#------------------------------------------------------------------------
def GetPlanetMassesToCarveGap_mJup(pltMasses_mJup, gapWidths_au):

	# Get the minimum gap width and the corresponding index
	minGapWidth_au = 1e99
	indexAtMinGapWidth = np.nan
	
	for i in range(len(gapWidths_au)):
		gapWidth_au = gapWidths_au[i]

		if gapWidth_au < minGapWidth_au:
			minGapWidth_au = gapWidth_au
			indexAtMinGapWidth = i
	
	# Continue only if the minimum gap width is less that the target gap 
	# width
	if minGapWidth_au > targetGapWidth_au: return np.nan, np.nan

	# Get the planet masses where the gap width is most similar to the 
	# target width 
	pltMass1_mJup = GetPlanetMassBetweenIndiciesWhereGapWidthClosestToTarget(0, indexAtMinGapWidth, pltMasses_mJup, gapWidths_au)
	pltMass2_mJup = GetPlanetMassBetweenIndiciesWhereGapWidthClosestToTarget(indexAtMinGapWidth, len(gapWidths_au), pltMasses_mJup, gapWidths_au)

	return pltMass1_mJup, pltMass2_mJup
	
#------------------------------------------------------------------------
def GetPlanetMassBetweenIndiciesWhereGapWidthClosestToTarget(startIndex, endIndex, pltMasses_mJup, gapWidths_au):

	smallestDifferenceBetweenGapWidthAndTarget_au = 1e99
	
	for i in range(startIndex, endIndex):
		gapWidth_au = gapWidths_au[i]
		differenceBetweenGapWidthAndTarget_au = abs(gapWidth_au - targetGapWidth_au)
		
		if differenceBetweenGapWidthAndTarget_au < smallestDifferenceBetweenGapWidthAndTarget_au:
			smallestDifferenceBetweenGapWidthAndTarget_au = differenceBetweenGapWidthAndTarget_au
			indexOfBestPltMass = i
	
	return pltMasses_mJup[indexOfBestPltMass]
	
#------------------------------------------------------------------------	
def GetInitialAndFinalLocationsOfSolutionPlanet(mPlt_mJup):
	'''Having found the mass of a planet required to carve the gap, find
	its initial and final semimajor axes'''
	
	h = GetHillRadius_PlanetSemimajorAxis(mPlt_mJup)

	startLocOfPlt_au = targetGapOuterEdge_au / (1. + hillRadiiClearedByPlt*h)

	endLocOfPlt_au = targetGapInnerEdge_au / (1. - hillRadiiClearedByPlt*h)

	return startLocOfPlt_au, endLocOfPlt_au

############################# Plot functions ############################
def MakePlot(pltMasses_mJup, gapWidths_au, gapWidthsIfNoMigration_au, approxGapWidths_au, pltMass1_mJup, pltMass2_mJup):

	print('Making plot...')
	
	# Title
	plt.suptitle(r'Pre-interaction disc mass: %s M$_\oplus$' % GetValueString(mDisc_mEarth))
	
	# Make the plot
	plt.loglog(pltMasses_mJup, gapWidths_au, label = 'Exact (Eq. 6)')
	plt.loglog(pltMasses_mJup, approxGapWidths_au, '--', label = 'Approximation (Eq. 7)')
	plt.loglog(pltMasses_mJup, gapWidthsIfNoMigration_au, 'k:', label = 'No migration')
										
	plt.xlabel(r'Planet mass / $M_{Jup}$')
	plt.ylabel(r'Gap width / au')

	# Plot the target gap width
	if targetGapWidth_au is not None:
		plt.loglog([pltMasses_mJup[0], pltMasses_mJup[-1]], [targetGapWidth_au, targetGapWidth_au], '-.', label = 'Target gap width')
	
		# Plot the planet masses that yield the target gap width
		plt.loglog(pltMass1_mJup, targetGapWidth_au, 'ko', label='Solution for target gap')
		plt.loglog(pltMass2_mJup, targetGapWidth_au, 'ko')
	
	plt.legend()
	
	# Axis limits
	plt.xlim(minPltMass_mJup, maxPltMass_mJup)

	print('Complete.')
	PrintEmptyLine()
	
	plt.show()
	
################################ Program ################################
PrintEmptyLine()

# Print the user inputs
PrintUserInputs()
	
# Check user inputs fine
areUserInputsOK = CheckUserInputsOK()

# Continue if the user inputs are OK
if areUserInputsOK:

	# Get the planet masses to consider
	pltMasses_mJup = np.logspace(log10MinPltMass_mJup, log10MaxPltMass_mJup, numberOfPltMasses)

	# Get the gap widths at the planet masses
	gapWidths_au = GetGapWidthsAtPlanetMasses_au(pltMasses_mJup)

	# Get the approximate gap width, assuming gap width << 2*gap radius
	approxGapWidths_au = GetApproxGapWidths_au(pltMasses_mJup)
	
	# Add some special points to the plot
	pltMasses_mJup, gapWidths_au, approxGapWidths_au = InsertSpecialPointsIntoPlanetMassList(pltMasses_mJup, gapWidths_au, approxGapWidths_au)
	
	# Get the gap widths if there is no migration
	gapWidthsIfNoMigration_au = GetGapWidthsIfNoMigrationAtPlanetMasses_au(pltMasses_mJup)

	# Get the solutions where the desired gap width is created
	if targetGapWidth_au is not None:
		pltMass1_mJup, pltMass2_mJup = GetPlanetMassesToCarveGap_mJup(pltMasses_mJup, gapWidths_au)

		if math.isnan(pltMass1_mJup) or math.isnan(pltMass2_mJup):
			print('RESULT: No planet between %s and %s MJup can carve a %s au gap centred on %s au, if the disc has initial mass %s mEarth.' % (GetValueString(minPltMass_mJup),GetValueString(maxPltMass_mJup),GetValueString(targetGapWidth_au), GetValueString(gapCentralRadius_au), GetValueString(mDisc_mEarth)))
			PrintEmptyLine()
					
		else:
			startLocOfPlt1_au, endLocOfPlt1_au = GetInitialAndFinalLocationsOfSolutionPlanet(pltMass1_mJup)
			startLocOfPlt2_au, endLocOfPlt2_au = GetInitialAndFinalLocationsOfSolutionPlanet(pltMass2_mJup)

			PrintGapCarvingPlanets(pltMass1_mJup, startLocOfPlt1_au, endLocOfPlt1_au, pltMass2_mJup, startLocOfPlt2_au, endLocOfPlt2_au)
	
	# Otherwise no target gap width has been specified						
	else:
		pltMass1_mJup, pltMass2_mJup = None, None

	# Plot the graph
	MakePlot(pltMasses_mJup, gapWidths_au, gapWidthsIfNoMigration_au, approxGapWidths_au, pltMass1_mJup, pltMass2_mJup)


