#!/usr/bin/env python3
import sys, os, argparse, copy
import configHelper, spectrumClasses
import matplotlib.pyplot, numpy, scipy
import generalUtils

class rvdata:
	def __init__(self, filename = 'rvs'):
		self.filename = filename
		self.rvs = []

	def setFilename(self, filename):
		self.filename = filename

	def sort(self):
		self.rvs = sorted(self.rvs, key=lambda k: k['HJD'])

	def addData(self, data):
		self.rvs.append(data)

	def saveDataTSV(self):
		outfile = open(self.filename + ".tsv", 'wt')
		print("Writing file: %s"%self.filename + ".tsv")
		for r in self.rvs:
			if r['good']: outfile.write("%f\t%f\t%f\n"%(r['HJD'], r['vel'], r['velError']))
		outfile.close()

	def saveData(self):
		outfile = open(self.filename + ".dat", 'wt')
		print("Writing file: %s"%self.filename + ".dat")
		outfile.write("# HJD, wavelength, wavelengthError, velocity, velocityError, redChiSq, goodFit\n")
		for r in self.rvs:
			outfile.write("%f, %f, %f, %f, %f, %f, "%(r['HJD'], r['wavelength'], r['wavelengthError'], r['vel'], r['velError'], r['redChiSq']))
			if r['good']: outfile.write("1")
			else: outfile.write("0")
			outfile.write("\n")
		outfile.close()


defaultConfiguration = {
	'width': 1.0,
	'blueRestWavelength': 8183.2556,
	'separation': 11.5349,
	'lowerWavelength': 8150,
	'upperWavelength': 8240,
	'fitLower': 8175,
	'fitUpper': 8205,
	'plotBG': 'black',
	'sigma': 4,
	'plot' : { 'background' : 'black', 'foreground' : 'white'}
	}

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads a spectrum JSON file and fits a double gaussian to a double line as defined in the config file.')
	parser.add_argument('inputFiles', type=str, nargs='+', help='JSON files containing the spectra')
	parser.add_argument('--list', action='store_true', help='Specify this option if the input file is actually a list of input files.')
	parser.add_argument('--show', action='store_true', help='Showed stored configuration and exit.')
	parser.add_argument('--clear', action='store_true', help='Clear stored configuration and exit.')
	parser.add_argument('--set', type=str, nargs='*', help='Set a parameter')
	arg = parser.parse_args()
	print(arg)

  	# Set up the matplotlib environment
	generalUtils.setMatplotlibDefaults()

	config = configHelper.config(debug = False)
	if not config._loaded:
		config.setProperties(defaultConfiguration)
		config.save()

	if arg.clear:
		config.clear()
		sys.exit()

	if arg.show:
		print(config)
		sys.exit()

	if arg.set is not None:
		print(arg.set)
		if len(arg.set)!=2:
			print("Please specify a parameter and a value.")
			sys.exit()
		config.set(arg.set[0], float(arg.set[1]))
		config.save()
		sys.exit()

	filenames = []
	if arg.list:
		# Load the list of files.
		if len(arg.inputFiles)>1:
			print("You can only give me one list of filenames.")
			sys.exit()
		filename = arg.inputFiles[0]
		fileList = open(filename, 'r')
		for line in fileList:
			filenames.append(str(line.strip()))
	else:
		filenames = arg.inputFiles

	spectra = []
	for fileIndex, f in enumerate(filenames):
		spectrum = spectrumClasses.spectrumObject()
		spectrum.loadFromJSON(f)
		print("%d: \t%s, contains object: %s."%(fileIndex+1, f, spectrum.objectName))
		spectra.append(spectrum)

	numSpectra = len(spectra)

	# Trim out the important region of the spectrum
	print("Discarding all info outside of the range %f to %f Angstroms."%(config.lowerWavelength, config.upperWavelength))
	for s in spectra:
		s.trimWavelengthRange(config.lowerWavelength, config.upperWavelength)

	rvInfo = rvdata()
	plotWidth = 8
	plotHeight = plotWidth/1.62
	spectrumOverview = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
	spectrumOverview.canvas.set_window_title('Overview')
	plotWidth = 6
	plotHeight = plotWidth/1.62
	fitZoom = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
	fitZoom.canvas.set_window_title('Fit')
	velocity = 0
	velocityError = 0
	wavelength = 0
	wavelengthError = 0
	reducedChiSq = 0
	for s in spectra:
		wavelengthGuess = config.blueRestWavelength
		repeatFit = True
		while repeatFit:
			wavelengths = s.getWavelengths()
			fluxes = s.getFlux()
			fluxErrors = s.getErrors()
			matplotlib.pyplot.figure(spectrumOverview.number)
			matplotlib.pyplot.title(s.HJD)
			matplotlib.pyplot.step(wavelengths, fluxes, where='mid', color='blue')
			matplotlib.pyplot.step(wavelengths, fluxErrors, color='red', where='mid')
			axes = matplotlib.pyplot.gca()
			(yLower, yUpper) = axes.get_ylim()
			matplotlib.pyplot.plot([config.fitLower, config.fitLower], [yLower, yUpper], linestyle=':', color='grey')
			matplotlib.pyplot.plot([config.fitUpper, config.fitUpper], [yLower, yUpper], linestyle=':', color='grey')
			matplotlib.pyplot.show(block = False)

			# Get the mean value of the continuum outside of the fitting region
			continuumSpectrum = copy.deepcopy(s)
			continuumSpectrum.snipWavelengthRange(config.fitLower, config.fitUpper)
			continuumMean = numpy.mean(continuumSpectrum.flux)
			matplotlib.pyplot.plot([numpy.min(wavelengths), numpy.max(wavelengths)], [continuumMean, continuumMean], linestyle=':', color='green')

			# Extract just the small region of the spectrum required for the fit
			fitSpectrum = copy.deepcopy(s)
			fitSpectrum.trimWavelengthRange(config.fitLower, config.fitUpper)
			wavelengths = fitSpectrum.getWavelengths()
			fluxes = fitSpectrum.getFlux()
			fluxErrors = fitSpectrum.getErrors()
			depthGuess = -1.0 * ((continuumMean - numpy.min(fluxes)) * 0.8)

			# Plot the spectrum
			matplotlib.pyplot.figure(fitZoom.number)
			matplotlib.pyplot.step(wavelengths, fluxes, where='mid', color='blue')

			# Fit the doubleGaussian
			a0 = continuumMean        			# Constant
			a1 = 0								# Slope
			a2 = depthGuess						# Depth
			a3 = wavelengthGuess				# Position of blue line
			bounds = ( [0, -0.2, -a0, config.fitLower], [2*a0, 0.2, 0, config.fitUpper] )
			# print("Bounds:", bounds)
			separation = config.separation
			width = config.width
			def doubleGaussian(x, a0, a1, a2, a3):
				global width, separation
				w = width
				s = separation
				y = a0 + a1 * (x-a3) + a2 * numpy.exp(-.5 * ((x-a3)/w)**2) + a2 * numpy.exp(-.5 * (((x-(a3+s))/w)**2) )
				return y

			goodValues = numpy.full(len(wavelengths), True, dtype=bool)
			rejectedPoints = 0
			newRejectedPoints = True
			iteration = 0
			while newRejectedPoints:
				# Plot the starting value
				xValues = numpy.arange(numpy.min(wavelengths), numpy.max(wavelengths), 0.1)
				yValues = doubleGaussian(xValues, a0, a1, a2, a3)
				matplotlib.pyplot.plot(xValues, yValues, color='green', linestyle=':')

				# Plot the fit and the spectrum
				matplotlib.pyplot.figure(fitZoom.number)
				matplotlib.pyplot.step(wavelengths, fluxes, where='mid', color='blue')

				(yLower, yUpper) = matplotlib.pyplot.gca().get_ylim()
				matplotlib.pyplot.step(wavelengths, [fe+yLower for fe in fluxErrors], color='red', where='mid')
				matplotlib.pyplot.show(block = False)

				guess = [a0, a1, a2, a3]
				x = numpy.array(wavelengths)
				y = numpy.array(fluxes)
				ye = numpy.array(fluxErrors)
				try:
					results, covariance = scipy.optimize.curve_fit(doubleGaussian, x[goodValues], y[goodValues], guess, ye[goodValues], absolute_sigma = True, bounds=bounds)
				except ValueError:
					print("Fit failed. Try to tweak")
					newRejectedPoints = False
					continue
				errors = numpy.sqrt(numpy.diag(covariance))
				(a0, a1, a2, a3) = results
				wavelength = a3
				wavelengthError = errors[3]
				print("Centroid blueward wavelength %f [%f] A"%(wavelength, wavelengthError))
				velocity = (wavelength - config.blueRestWavelength)/config.blueRestWavelength * 3E5
				velocityError = 3E5 / config.blueRestWavelength * wavelengthError
				print("Velocity %f [%f] km/s"%(velocity, velocityError))

				# Plot the fitted curve
				xValues = numpy.arange(numpy.min(wavelengths), numpy.max(wavelengths), 0.1)
				yValues = doubleGaussian(xValues, a0, a1, a2, a3)
				matplotlib.pyplot.plot(xValues, yValues, color='green', linestyle='-')

				# Calculate the chiSquared
				chiSq = 0
				sigma = 0
				for w, flux, fluxError in zip(x[goodValues], y[goodValues], ye[goodValues]):
					fittedFlux = doubleGaussian(w, a0, a1, a2, a3)
					chiSq+= ((flux - fittedFlux)/fluxError)**2
				print("Chi squared", chiSq)
				reducedChiSq = chiSq / (numpy.sum(goodValues) - 4)
				print("Reduced Chi squared", reducedChiSq)

				# Calculate the residuals

				residuals = numpy.array([abs(f - doubleGaussian(w, a0, a1, a2, a3)) for w, f, fe in zip(wavelengths, fluxes, fluxErrors)])
				scaledResiduals = [r/(fe*numpy.sqrt(reducedChiSq)) for r, fe in zip(residuals, fluxErrors)]
				for w, r, rs, fe in zip(wavelengths, residuals, scaledResiduals, fluxErrors):
					print(w, r, rs, fe)
				worstResidual = numpy.max(numpy.array(scaledResiduals)[goodValues])
				index = numpy.argmax(numpy.array(scaledResiduals)[goodValues])
				print("Worst residual: ", wavelengths[index], worstResidual, scaledResiduals[index])
				if worstResidual < config.sigma:
					newRejectedPoints = False
				else:
					goodValues[index] = False

				# goodValues = scaledResiduals < numpy.full(len(wavelengths), config.sigma)
				goodResiduals = residuals[goodValues]
				matplotlib.pyplot.scatter(x[goodValues], [r + yLower for r in goodResiduals], marker='+', color='green')
				badResiduals = residuals[numpy.logical_not(goodValues)]
				matplotlib.pyplot.scatter(x[numpy.logical_not(goodValues)], [r + yLower for r in badResiduals], marker='x', color='red')

				iteration+= 1
				print("Iteration: %d    %d rejected points lying outside %2.1f sigma."%(iteration, numpy.sum(numpy.logical_not(goodValues)), config.sigma))
				if iteration > 10:
					print("The fit is not converging, rejecting fit.")
					break
				# generalUtils.query_yes_no("Continue?")

			# Pause for input
			print("Happy with the fit? ([Y], [n] or retry with a tweak to the [l]eft or [r]ight. ")
			choice = input().lower()

			if choice=='n':
				goodFit = False
				repeatFit = False
			elif choice=='l':
				wavelengthGuess-= 2
				repeatFit = True
			elif choice=='r':
				wavelengthGuess+= 2
				repeatFit = True
			else:
				goodFit = True
				repeatFit = False


			# Clear the old plots
			matplotlib.pyplot.figure(spectrumOverview.number)
			matplotlib.pyplot.clf()
			matplotlib.pyplot.figure(fitZoom.number)
			matplotlib.pyplot.clf()

		saveRV = { 'HJD': s.HJD,
					'vel': velocity, 'velError': velocityError,
					'wavelength': wavelength, 'wavelengthError': wavelengthError,
					'redChiSq': reducedChiSq,
					'good': goodFit }
		rvInfo.addData(saveRV)
		rvInfo.setFilename(s.objectName)
		rvInfo.sort()
		rvInfo.saveDataTSV()
		rvInfo.saveData()
