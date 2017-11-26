import os, sys
import json, numpy

class config:
	def __init__(self, filename = None, debug = False):
		self._debug = debug
		self._homeDir = os.path.expanduser('~')
		self._appName = sys.argv[0].split('/')[-1]
		self._appName = self._appName.split('.')[0]
		self._loaded = False
		if filename is None:
			if self._debug: print("No filename specified, will make my own")
			configFilename = self._homeDir + "/.config/" + self._appName + "/" + self._appName + ".conf"
			self._configFilename = configFilename
		else:
			self._configFilename = filename

		if self._debug: print("The configuration information will be stored in: %s"%self._configFilename)

		if os.path.exists(self._configFilename):
			self.load()
			self._loaded = True
		else:
			print("No existing config found.")
			self._loaded = False

	def __str__(self):
		retString = ""
		for key in self.__dict__.keys():
			if key[0]=='_': continue         # Don't return properties that start with an underscore
			retString+= str(key) + ": " + str(getattr(self, key)) + "\n"
		return retString


	def createConfigFolder(self, appName):
		configPath = self._homeDir + "/.config"
		if not os.path.exists(configPath):
			os.mkdir(configPath)
			if self._debug: print("DEBUG: Creating directory %s"%configPath)

		fullConfigPath = configPath + "/" + appName
		if not os.path.exists(fullConfigPath):
			os.mkdir(fullConfigPath)
			if self._debug: print("DEBUG: Creating directory %s"%fullConfigPath)

	def setProperties(self, properties):
		for key in properties:
			if self._debug: print("Adding:", key, properties[key])
			setattr(self, key, properties[key])

	def set(self, key, value):
		setattr(self, key, value)

	def load(self):
		filename = self._configFilename
		inputfile = open(filename, "r")
		jsonObject = json.load(inputfile)
		for key in jsonObject.keys():
			keyString = str(key)
			value = jsonObject[key]
			# print(type(value))
			# if type(value) is unicode:
			#  	value = str(value)
			if type(value) is list:
				value = numpy.array(value)
			if self._debug: print("Loading", key, value)
			setattr(self, key, value)
		inputfile.close()
		return True

	def save(self):
		if self._debug: print("Writing config info to %s"%(self._configFilename))
		filename = self._configFilename
		object = {}
		for key in self.__dict__.keys():
			if key[0]=='_': continue         # Don't write properties that start with an underscore
			data = getattr(self, key)
			if type(data)==numpy.float32:
				data = float(data)
			if type(data)==numpy.ndarray:
				data = numpy.array(data).tolist()
			if type(data)==list:
				data = numpy.array(data).tolist()
			object[key] = data

		outputfile = open(filename, 'w')
		json.dump(object, outputfile, indent = 4)
		outputfile.write("\n")
		outputfile.close()

	def __str__(self):
		retString = ""
		for key in self.__dict__.keys():
			if key[0]=='_': continue         # Don't return properties that start with an underscore
			retString+= str(key) + ": " + str(getattr(self, key)) + "\n"
		return retString

	def clear(self):
		if os.path.exists(self._configFilename):
			os.remove(self._configFilename)
