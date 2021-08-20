// Copyright (C) 2016 Vicini Fabio
//
// This file is part_fracture of the dissertation of the author.
//
// This is a free program: you can redistribute it and/or modify
// it under the terms of the author.
//
// This program is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.
//
// Modified by Vicini Fabio 2016
//
// First added:  2016-05-13

#include "MacroDefinitions.hpp"

#if USE_MPI == 1
#include <mpi.h>
#endif

#include <sys/stat.h>
#include <sys/types.h>
#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>
#include <sstream>
#include <unistd.h>

#include "ConfigFile.hpp"
#include "Output.hpp"


using namespace std;
using namespace MainApplication;

namespace Configuration
{
	// ***************************************************************************
	string ConfigFile::defaultPropertyValue = "<default>";
	InputParameters ConfigFile::inputParam;
	CGParameters ConfigFile::cGParam;
	GeometricParameters ConfigFile::geomParam;
	DiscretizationParameters ConfigFile::discParam;
	ParallelParameters ConfigFile::paralParam;
	GraphicParameters ConfigFile::graphicParam;
	OutputParameters ConfigFile::outputParam;
	DebugParameters ConfigFile::debugParam;
	FileNameParameters ConfigFile::filesParam;
	LogParameters ConfigFile::logParam;
  ModelParameters ConfigFile::modelParam;
	// **************************************************************************
	bool ConfigFile::ConfigurationFileExists(const string& nameFile, const string& nameFolder, const string& description)
	{
		#if USE_MPI == 1
		int rank = MPI::COMM_WORLD.Get_rank();
		#elif USE_MPI == 0
		int rank = 0;
		#endif // USE_MPI

		bool fileExists = false;

		if(rank == 0)
		{
			ostringstream fileName;
			fileName<< nameFolder<< nameFile;

			fileExists = Output::FileExists(fileName.str());
			if(!fileExists)
			{
				ostringstream errorMessage;
				errorMessage<< "Configuration file '"<< nameFile<< "' created in folder '"<< nameFolder<< "'";
				Output::PrintWarningMessage(errorMessage.str(), true);

				ofstream createdFile(fileName.str().c_str());
				if(!description.empty())
				{
					vector<string> descriptions = ConfigFile::StringSplit(description, '\n');
					for(unsigned int i=0;i<descriptions.size();i++)
						createdFile<< "# "<< descriptions[i]<< endl;
				}
				createdFile<< "#####################################################"<< endl;
			}

			ifstream inFile;
			inFile.open(fileName.str().c_str());
			if(!inFile.is_open())
			{
				Output::PrintErrorMessage("Configuration file '" + nameFile + "' cannot be open", true);
				fileExists = false;
			}
			else
			{
				inFile.close();
				fileExists = true;
			}
		}

		#if USE_MPI == 1
		MPI::COMM_WORLD.Bcast(&fileExists, 1, MPI::BOOL, 0);
		#endif // USE_MPI

		return fileExists;
	}
	// *************************************************************************
	vector<string> ConfigFile::StringSplit(const string& stringToSplit, const char character)
	{
		stringstream splitter(stringToSplit);
		string tempString;
		vector<string> strings;

		while(getline(splitter, tempString, character))
			strings.push_back(tempString);

		return strings;
	}
	// *************************************************************************
	bool ConfigFile::FindStringProperty(const string& nameProperty, const string defaultValue, string& propertyValue, int argc, char** argv)
	{
		bool propertyFound = false;
		propertyValue = defaultValue;

    for(int i=0; i<argc && !propertyFound;i++)
		{
			vector<string> strings = ConfigFile::StringSplit(argv[i], '=');
			if(strings.size() < 2)
				continue;

			if(strings[0] == nameProperty)
			{
				propertyValue = defaultValue.empty() ? strings[1] : (strings[1].empty() ? defaultValue : strings[1]);
				propertyFound = true;
			}
		}

		return propertyFound;
	}
	// *************************************************************************
	bool ConfigFile::ReadStringProperty(const string& nameProperty, const string defaultValue, string& propertyValue, const string& nameFile, const string& nameFolder, const string& description)
	{
		#if USE_MPI == 1
		int rank = MPI::COMM_WORLD.Get_rank();
		#elif USE_MPI == 0
		int rank = 0;
		#endif // USE_MPI

		bool valueFound = false;

		if(rank == 0)
		{
			ostringstream fileName;
			fileName<< nameFolder<< nameFile;

			ifstream inFile;
			inFile.open(fileName.str().c_str());

			if(!inFile.is_open())
			{
        Output::PrintErrorMessage("Cannot open file '%s' in folder '%s'", false, nameFile.c_str(), nameFolder.c_str());
        exit(-1);
			}

			// Find Property in File
			bool propertyFound = false;
			istringstream convert;
			string line, fileNameProperty;
			string fileValue;
			while (!inFile.eof() && !propertyFound)
			{
				getline(inFile, line);
				// Skip Comment Line
				if(line[0] == '#')
					continue;

				convert.clear();
				convert.str(line);

				if(!(convert >> fileNameProperty))
					continue;
				if(nameProperty != fileNameProperty)
					continue;
				propertyFound = true;

				if(!(convert >> fileValue))
					continue;

				valueFound = true;
			}

			inFile.close();

			if(propertyFound && valueFound)
			{
				if(!defaultValue.empty() && fileValue.empty())
				{
					ostringstream errorMessage;
					errorMessage<< "Property '"<< nameProperty<< "' is empty in file '"<< fileName.str()<< "'.\nSet default value: "<< defaultValue;
					Output::PrintWarningMessage(errorMessage.str(), true);

					propertyValue = defaultValue;
				}
				else
					propertyValue = fileValue;
			}
			else if(propertyFound)
			{
				ostringstream errorMessage;
				errorMessage<< "Wrong format of property '"<< nameProperty<< "' in file '"<< fileName.str()<< "'.\nSet default value: "<< defaultValue;
				Output::PrintWarningMessage(errorMessage.str(), true);

				propertyValue = defaultValue;
			}
			else
			{
				// Write Default Value in File
				ofstream outFile(fileName.str().c_str(), std::ios_base::app);
				if(!description.empty())
				{
					vector<string> descriptions = ConfigFile::StringSplit(description, '\n');
					for(unsigned int i=0;i<descriptions.size();i++)
						outFile<< "# "<< descriptions[i]<< endl;
				}
				outFile<< nameProperty<< " "<< defaultValue<< "\n";
				outFile<< "#####################################################"<< endl;

				outFile.close();
				propertyValue = defaultValue;
			}
		}

    int stringLenght = propertyValue.size();
    #if USE_MPI == 1
		MPI::COMM_WORLD.Bcast(&stringLenght, 1, MPI::INT, 0);
		#endif // USE_MPI

		char* valueCStr = new char[stringLenght + 1];
		if(rank==0)
			strcpy(valueCStr, propertyValue.c_str());

		#if USE_MPI == 1
		MPI::COMM_WORLD.Bcast(valueCStr, stringLenght + 1, MPI::CHAR, 0);
		#endif // USE_MPI
		propertyValue.assign(valueCStr);

		#if USE_MPI == 1
		MPI::COMM_WORLD.Bcast(&valueFound, 1, MPI::BOOL, 0);
		#endif // USE_MPI

		delete[] valueCStr;

		return valueFound;
	}
	// *************************************************************************
	bool ConfigFile::FindIntProperty(const string& nameProperty, const int defaultValue, int& propertyValue, int argc, char** argv)
	{
		bool propertyFound = false;
		propertyValue = defaultValue;

    for(int i=0; i<argc && !propertyFound;i++)
		{
			istringstream convert;
			vector<string> strings = ConfigFile::StringSplit(argv[i], '=');
			if(strings.size() < 2)
				continue;

			if(strings[0] == nameProperty)
			{
				convert.str(strings[1]);
				if(!(convert >> propertyValue))
				{
					propertyValue = defaultValue;
					continue;
				}
				propertyFound = true;
			}
		}

		return propertyFound;
	}
	// *************************************************************************
	bool ConfigFile::ReadIntProperty(const string& nameProperty, const int defaultValue, int& propertyValue, const string& nameFile, const string& nameFolder, const string& description)
	{
		#if USE_MPI == 1
		int rank = MPI::COMM_WORLD.Get_rank();
		#elif USE_MPI == 0
		int rank = 0;
		#endif // USE_MPI

		bool valueFound = false;

		if(rank == 0)
		{
			ostringstream fileName;
			fileName<< nameFolder<< nameFile;

			ifstream inFile;
			inFile.open(fileName.str().c_str());

			// Find Property in File
			bool propertyFound = false;
			istringstream convert;
			string line, fileNameProperty;
			int fileValue;
			while (!inFile.eof() && !propertyFound)
			{
				getline(inFile, line);
				// Skip Comment Line
				if(line[0] == '#')
					continue;

				convert.clear();
				convert.str(line);
				if(!(convert >> fileNameProperty))
					continue;
				if(nameProperty != fileNameProperty)
					continue;
				propertyFound = true;

				if(!(convert >> fileValue))
					continue;

				valueFound = true;
			}

			inFile.close();

			if(propertyFound && valueFound)
				propertyValue = fileValue;
			else if(propertyFound)
			{
				ostringstream errorMessage;
				errorMessage<< "Wrong format of property '"<< nameProperty<< "' in file '"<< fileName.str()<< "'.\nSet default value: "<< defaultValue;
				Output::PrintWarningMessage(errorMessage.str(), true);

				propertyValue = defaultValue;
			}
			else
			{
				// Write Default Value in File
				ofstream outFile(fileName.str().c_str(), std::ios_base::app);
				if(!description.empty())
				{
					vector<string> descriptions = ConfigFile::StringSplit(description, '\n');
					for(unsigned int i=0;i<descriptions.size();i++)
						outFile<< "# "<< descriptions[i]<< endl;
				}
				outFile<< nameProperty<< " "<< defaultValue<< endl;
				outFile<< "#####################################################"<< endl;

				outFile.close();
				propertyValue = defaultValue;
			}
		}

		#if USE_MPI == 1
		MPI::COMM_WORLD.Bcast(&propertyValue, 1, MPI::INT, 0);
		MPI::COMM_WORLD.Bcast(&valueFound, 1, MPI::BOOL, 0);
		#endif // USE_MPI

		return valueFound;
	}
	// *************************************************************************
	bool ConfigFile::FindIntProperty(const string& nameProperty, const vector<int> defaultValues, vector<int>& propertyValues, int argc, char** argv)
	{
		bool propertyFound = false;
		propertyValues.resize(defaultValues.size());
		memcpy(&propertyValues[0], &defaultValues[0], defaultValues.size() * sizeof(int));

    for(int i=0; i<argc && !propertyFound;i++)
		{
			vector<string> strings = ConfigFile::StringSplit(argv[i], '=');
			if(strings.size() < 2)
				continue;

			if(strings[0] == nameProperty)
			{
        vector<string> values = ConfigFile::StringSplit(&strings[1][1], ',');
				int tempValue;
				propertyValues.clear();
				propertyValues.reserve(0);
				for(unsigned int j=0; j< values.size(); j++)
				{
					istringstream convert;
					convert.str(values[j]);
					if(!(convert >> tempValue))
						continue;

					propertyValues.push_back(tempValue);
				}

				if(propertyValues.size() > 0)
					propertyFound = true;
				else
				{
					propertyValues.resize(defaultValues.size());
					memcpy(&propertyValues[0], &defaultValues[0], defaultValues.size() * sizeof(int));
				}
			}
		}

		return propertyFound;
	}
	// *************************************************************************
	bool ConfigFile::ReadIntProperty(const string& nameProperty, const vector<int> defaultValues, vector<int>& propertyValues, const string& nameFile, const string& nameFolder, const string& description)
	{
		#if USE_MPI == 1
		int rank = MPI::COMM_WORLD.Get_rank();
		#elif USE_MPI == 0
		int rank = 0;
		#endif // USE_MPI

		bool valueFound = false;

		if(rank == 0)
		{
			ostringstream fileName;
			fileName<< nameFolder<< nameFile;

			ifstream inFile;
			inFile.open(fileName.str().c_str());

			// Find Property in File
			bool propertyFound = false;
			istringstream convert;
			string line, fileNameProperty;
			int fileValue;
			while (!inFile.eof() && !propertyFound)
			{
				getline(inFile, line);
				// Skip Comment Line
				if(line[0] == '#')
					continue;

				convert.clear();
				convert.str(line);
				if(!(convert >> fileNameProperty))
					continue;
				if(nameProperty != fileNameProperty)
					continue;
				propertyFound = true;

				propertyValues.resize(0);
				while(convert >> fileValue)
					propertyValues.push_back(fileValue);

				valueFound = true;
			}

			inFile.close();

			if(propertyFound && valueFound)
			{
				if(defaultValues.size() > 0 && propertyValues.size() == 0)
				{
					ostringstream errorMessage;
					errorMessage<< "Property '"<< nameProperty<< "' empty in file '"<< fileName.str()<< "'.\nSet default values: "<< defaultValues;
					Output::PrintWarningMessage(errorMessage.str(), true);

					propertyValues.resize(defaultValues.size());
					memcpy(&propertyValues[0], &defaultValues[0], defaultValues.size() * sizeof(int));
				}
			}
			else if(propertyFound)
			{
				ostringstream errorMessage;
				errorMessage<< "Wrong format of property '"<< nameProperty<< "' in file '"<< fileName.str()<< "'.\nSet default values: "<< defaultValues;
				Output::PrintWarningMessage(errorMessage.str(), true);

				propertyValues.resize(defaultValues.size());
				memcpy(&propertyValues[0], &defaultValues[0], defaultValues.size() * sizeof(int));
			}
			else
			{
				// Write Default Value in File
				ofstream outFile(fileName.str().c_str(), std::ios_base::app);
				if(!description.empty())
				{
					vector<string> descriptions = ConfigFile::StringSplit(description, '\n');
					for(unsigned int i=0;i<descriptions.size();i++)
						outFile<< "# "<< descriptions[i]<< endl;
				}
				outFile<< nameProperty;
				for(unsigned int i=0; i<defaultValues.size(); i++)
					outFile<< " "<< defaultValues[i];
				outFile<< endl;
				outFile<< "#####################################################"<< endl;

				outFile.close();
				propertyValues.resize(defaultValues.size());
				memcpy(&propertyValues[0], &defaultValues[0], defaultValues.size() * sizeof(int));
			}
		}

		int propertyValueSize = propertyValues.size();
		#if USE_MPI == 1
		MPI::COMM_WORLD.Bcast(&propertyValueSize, 1, MPI::INT, 0);
		#endif // USE_MPI
		if(propertyValueSize > 0)
		{
			propertyValues.resize(propertyValueSize);
			#if USE_MPI == 1
			MPI::COMM_WORLD.Bcast(&propertyValues[0], propertyValueSize, MPI::INT, 0);
			#endif // USE_MPI
		}

		#if USE_MPI == 1
		MPI::COMM_WORLD.Bcast(&valueFound, 1, MPI::BOOL, 0);
		#endif // USE_MPI

		return valueFound;
	}
	// *************************************************************************
	bool ConfigFile::FindDoubleProperty(const string& nameProperty, const double defaultValue, double& propertyValue, int argc, char** argv)
	{
		bool propertyFound = false;
		propertyValue = defaultValue;

    for(int i=0; i<argc && !propertyFound;i++)
		{
			istringstream convert;
			vector<string> strings = ConfigFile::StringSplit(argv[i], '=');
			if(strings.size() < 2)
				continue;

			if(strings[0] == nameProperty)
			{
				convert.str(strings[1]);
				if(!(convert >> propertyValue))
				{
					propertyValue = defaultValue;
					continue;
				}
				propertyFound = true;
			}
		}

		return propertyFound;
	}
	// *************************************************************************
	bool ConfigFile::ReadDoubleProperty(const string& nameProperty, const double defaultValue, double& propertyValue, const string& nameFile, const string& nameFolder, const string& description)
	{
		#if USE_MPI == 1
		int rank = MPI::COMM_WORLD.Get_rank();
		#elif USE_MPI == 0
		int rank = 0;
		#endif // USE_MPI

		bool valueFound = false;

		if(rank == 0)
		{
			ostringstream fileName;
			fileName<< nameFolder<< nameFile;

			ifstream inFile;
			inFile.open(fileName.str().c_str());

			// Find Property in File
			bool propertyFound = false;
			istringstream convert;
			string line, fileNameProperty;
			double fileValue;
			while (!inFile.eof() && !propertyFound)
			{
				getline(inFile, line);
				// Skip Comment Line
				if(line[0] == '#')
					continue;

				convert.clear();
				convert.str(line);
				if(!(convert >> fileNameProperty))
					continue;
				if(nameProperty != fileNameProperty)
					continue;
				propertyFound = true;

				if(!(convert >> fileValue))
					continue;

				valueFound = true;
			}

			inFile.close();

			if(propertyFound && valueFound)
				propertyValue = fileValue;
			else if(propertyFound)
			{
				ostringstream errorMessage;
				errorMessage<< "Wrong format of property '"<< nameProperty<< "' in file '"<< fileName.str()<< "'.\nSet default value: "<< defaultValue;
				Output::PrintWarningMessage(errorMessage.str(), true);

				propertyValue = defaultValue;
			}
			else
			{
				// Write Default Value in File
				ofstream outFile(fileName.str().c_str(), std::ios_base::app);
				if(!description.empty())
				{
					vector<string> descriptions = ConfigFile::StringSplit(description, '\n');
					for(unsigned int i=0;i<descriptions.size();i++)
						outFile<< "# "<< descriptions[i]<< endl;
				}
				outFile<< nameProperty<< " "<< defaultValue<< "\n";
				outFile<< "#####################################################"<< endl;

				outFile.close();
				propertyValue = defaultValue;
			}
		}

		#if USE_MPI == 1
		MPI::COMM_WORLD.Bcast(&propertyValue, 1, MPI::DOUBLE, 0);
		MPI::COMM_WORLD.Bcast(&valueFound, 1, MPI::BOOL, 0);
		#endif // USE_MPI

		return valueFound;
	}
	// *************************************************************************
	bool ConfigFile::FindDoubleProperty(const string& nameProperty, const vector<double> defaultValues, vector<double>& propertyValues, int argc, char** argv)
	{
		bool propertyFound = false;
		propertyValues.resize(defaultValues.size());
		memcpy(&propertyValues[0], &defaultValues[0], defaultValues.size() * sizeof(double));

    for(int i=0; i<argc && !propertyFound;i++)
		{
			vector<string> strings = ConfigFile::StringSplit(argv[i], '=');
			if(strings.size() < 2)
				continue;

			if(strings[0] == nameProperty)
			{
        vector<string> values = ConfigFile::StringSplit(&strings[1][1], ',');

				double tempValue;
				propertyValues.clear();
				propertyValues.reserve(values.size());
				for(unsigned int j=0; j< values.size(); j++)
				{
					istringstream convert;
					convert.str(values[j]);
					if(!(convert >> tempValue))
						continue;

					propertyValues.push_back(tempValue);
				}

				if(propertyValues.size() > 0)
					propertyFound = true;
				else
				{
					propertyValues.resize(defaultValues.size());
					memcpy(&propertyValues[0], &defaultValues[0], defaultValues.size() * sizeof(double));
				}
			}
		}

		return propertyFound;
	}
	// *************************************************************************
	bool ConfigFile::ReadDoubleProperty(const string& nameProperty, const vector<double> defaultValues, vector<double>& propertyValues, const string& nameFile, const string& nameFolder, const string& description)
	{
		#if USE_MPI == 1
		int rank = MPI::COMM_WORLD.Get_rank();
		#elif USE_MPI == 0
		int rank = 0;
		#endif // USE_MPI

		bool valueFound = false;

		if(rank == 0)
		{
			ostringstream fileName;
			fileName<< nameFolder<< nameFile;

			ifstream inFile;
			inFile.open(fileName.str().c_str());

			// Find Property in File
			bool propertyFound = false;
			istringstream convert;
			string line, fileNameProperty;
			double fileValue;
			while (!inFile.eof() && !propertyFound)
			{
				getline(inFile, line);
				// Skip Comment Line
				if(line[0] == '#')
					continue;

				convert.clear();
				convert.str(line);
				if(!(convert >> fileNameProperty))
					continue;
				if(nameProperty != fileNameProperty)
					continue;
				propertyFound = true;

				propertyValues.resize(0);
				while(convert >> fileValue)
					propertyValues.push_back(fileValue);

				valueFound = true;
			}

			inFile.close();

			if(propertyFound && valueFound)
			{
				if(defaultValues.size() > 0 && propertyValues.size() == 0)
				{
					ostringstream errorMessage;
					errorMessage<< "Property '"<< nameProperty<< "' empty in file '"<< fileName.str()<< "'.\nSet default values: "<< defaultValues;
					Output::PrintWarningMessage(errorMessage.str(), true);

					propertyValues = defaultValues;
				}
			}
			else if(propertyFound)
			{
				ostringstream errorMessage;
				errorMessage<< "Wrong format of property '"<< nameProperty<< "' in file '"<< fileName.str()<< "'.\nSet default values: "<< defaultValues;
				Output::PrintWarningMessage(errorMessage.str(), true);

				propertyValues = defaultValues;
			}
			else
			{
				// Write Default Value in File
				ofstream outFile(fileName.str().c_str(), std::ios_base::app);
				if(!description.empty())
				{
					vector<string> descriptions = ConfigFile::StringSplit(description, '\n');
					for(unsigned int i=0;i<descriptions.size();i++)
						outFile<< "# "<< descriptions[i]<< endl;
				}
				outFile<< nameProperty;
				for(unsigned int i=0; i<defaultValues.size(); i++)
					outFile<< " "<< defaultValues[i];
				outFile<< endl;
				outFile<< "#####################################################"<< endl;

				outFile.close();
				propertyValues = defaultValues;
			}
		}

		int propertyValueSize = propertyValues.size();
		#if USE_MPI == 1
		MPI::COMM_WORLD.Bcast(&propertyValueSize, 1, MPI::INT, 0);
		#endif // USE_MPI
		if(propertyValueSize > 0)
		{
			propertyValues.resize(propertyValueSize);
			#if USE_MPI == 1
			MPI::COMM_WORLD.Bcast(&propertyValues[0], propertyValueSize, MPI::DOUBLE, 0);
			#endif // USE_MPI
		}

		#if USE_MPI == 1
		MPI::COMM_WORLD.Bcast(&valueFound, 1, MPI::BOOL, 0);
		#endif // USE_MPI

		return valueFound;
	}
	// *************************************************************************
	bool ConfigFile::FindBoolProperty(const string& nameProperty, const bool defaultValue, bool& propertyValue, int argc, char** argv)
	{
		bool propertyFound = false;
		propertyValue = defaultValue;

    for(int i=0; i<argc && !propertyFound;i++)
		{
			istringstream convert;
			vector<string> strings = ConfigFile::StringSplit(argv[i], '=');
			if(strings.size() < 2)
				continue;

			if(strings[0] == nameProperty)
			{
				convert.str(strings[1]);
				if(!(convert >> propertyValue))
				{
					propertyValue = defaultValue;
					continue;
				}
				propertyFound = true;
			}
		}

		return propertyFound;
	}
	// *************************************************************************
	bool ConfigFile::ReadBoolProperty(const string& nameProperty, const bool defaultValue, bool& propertyValue, const string& nameFile, const string& nameFolder, const string& description)
	{
		#if USE_MPI == 1
		int rank = MPI::COMM_WORLD.Get_rank();
		#elif USE_MPI == 0
		int rank = 0;
		#endif // USE_MPI

		bool valueFound = false;

		if(rank == 0)
		{
			ostringstream fileName;
			fileName<< nameFolder<< nameFile;

			ifstream inFile;
			inFile.open(fileName.str().c_str());

			// Find Property in File
			bool propertyFound = false;
			istringstream convert;
			string line, fileNameProperty;
			bool fileValue;
			while (!inFile.eof() && !propertyFound)
			{
				getline(inFile, line);
				// Skip Comment Line
				if(line[0] == '#')
					continue;

				convert.clear();
				convert.str(line);
				if(!(convert >> fileNameProperty))
					continue;
				if(nameProperty != fileNameProperty)
					continue;
				propertyFound = true;

				if(!(convert >> fileValue))
					continue;

				valueFound = true;
			}

			inFile.close();

			if(propertyFound && valueFound)
				propertyValue = fileValue;
			else if(propertyFound)
			{
				ostringstream errorMessage;
				errorMessage<< "Wrong format of property '"<< nameProperty<< "' in file '"<< fileName.str()<< "'.\nSet default value: "<< defaultValue;
				Output::PrintWarningMessage(errorMessage.str(), true);

				propertyValue = defaultValue;
			}
			else
			{
				// Write Default Value in File
				ofstream outFile(fileName.str().c_str(), std::ios_base::app);
				if(!description.empty())
				{
					vector<string> descriptions = ConfigFile::StringSplit(description, '\n');
					for(unsigned int i=0;i<descriptions.size();i++)
						outFile<< "# "<< descriptions[i]<< endl;
				}
				outFile<< nameProperty<< " "<< defaultValue<< endl;
				outFile<< "#####################################################"<< endl;

				outFile.close();
				propertyValue = defaultValue;
			}
		}

		#if USE_MPI == 1
		MPI::COMM_WORLD.Bcast(&propertyValue, 1, MPI::BOOL, 0);
		MPI::COMM_WORLD.Bcast(&valueFound, 1, MPI::BOOL, 0);
		#endif // USE_MPI

		return valueFound;
	}
	// *************************************************************************
	bool ConfigFile::CheckIntProperty(const string& nameFile, const string& nameProperty, const int defaultValue, int& propertyValue, const int& minValue, const int& maxValue)
	{
		bool isPropertyCorrect = true;

		if(minValue > -INT_MAX && propertyValue < minValue)
		{
			ostringstream errorMessage;
			errorMessage<< "Property '"<< nameProperty<< "' not correct in file '"<< nameFile<< "'.\nSet default value: "<< defaultValue;
			Output::PrintWarningMessage(errorMessage.str(), true);
			isPropertyCorrect = false;
			propertyValue = defaultValue;
		}

		if(maxValue < INT_MAX && propertyValue > maxValue)
		{
			ostringstream errorMessage;
			errorMessage<< "Property '"<< nameProperty<< "' not correct in file '"<< nameFile<< "'.\nSet default value: "<< defaultValue;
			Output::PrintWarningMessage(errorMessage.str(), true);
			isPropertyCorrect = false;
			propertyValue = defaultValue;
		}

		return isPropertyCorrect;
	}
	// *************************************************************************
	bool ConfigFile::CheckIntProperty(const string& nameFile, const string& nameProperty, const vector<int> defaultValues, vector<int>& propertyValues, const int& minValue, const int& maxValue)
	{
		bool isPropertyCorrect = true;

		if(defaultValues.size() > 0 && propertyValues.empty())
		{
			ostringstream errorMessage;
			errorMessage<< "Property '"<< nameProperty<< "' empty in file '"<< nameFile<< "'.\nSet default value: "<< defaultValues;
			Output::PrintWarningMessage(errorMessage.str(), true);
			isPropertyCorrect = false;
			propertyValues.resize(defaultValues.size());
			memcpy(&propertyValues[0], &defaultValues[0], defaultValues.size() * sizeof(int));
		}

		for(unsigned int i=0;i<propertyValues.size() && isPropertyCorrect; i++)
		{
			if(minValue > -INT_MAX && propertyValues[i] < minValue)
			{
				ostringstream errorMessage;
				errorMessage<< "Property '"<< nameProperty<< "' not correct in file '"<< nameFile<< "'.\nSet default value: "<< defaultValues;
				Output::PrintWarningMessage(errorMessage.str(), true);
				isPropertyCorrect = false;
				propertyValues.resize(defaultValues.size());
				memcpy(&propertyValues[0], &defaultValues[0], defaultValues.size() * sizeof(int));
			}

			if(maxValue < INT_MAX && propertyValues[i] > maxValue)
			{
				ostringstream errorMessage;
				errorMessage<< "Property '"<< nameProperty<< "' not correct in file '"<< nameFile<< "'.\nSet default value: "<< defaultValues;
				Output::PrintWarningMessage(errorMessage.str(), true);
				isPropertyCorrect = false;
				propertyValues.resize(defaultValues.size());
				memcpy(&propertyValues[0], &defaultValues[0], defaultValues.size() * sizeof(int));
			}
		}

		return isPropertyCorrect;
	}
	// *************************************************************************
	bool ConfigFile::CheckDoubleProperty(const string& nameFile, const string& nameProperty, const double defaultValue, double& propertyValue, const double& minValue, const double& maxValue)
	{
		bool isPropertyCorrect = true;

		if(minValue > -DBL_MAX && propertyValue < minValue)
		{
			ostringstream errorMessage;
			errorMessage<< "Property '"<< nameProperty<< "' not correct in file '"<< nameFile<< "'.\nSet default value: "<< defaultValue;
			Output::PrintWarningMessage(errorMessage.str(), true);
			isPropertyCorrect = false;
			propertyValue = defaultValue;
		}

		if(maxValue < DBL_MAX && propertyValue > maxValue)
		{
			ostringstream errorMessage;
			errorMessage<< "Property '"<< nameProperty<< "' not correct in file '"<< nameFile<< "'.\nSet default value: "<< defaultValue;
			Output::PrintWarningMessage(errorMessage.str(), true);
			isPropertyCorrect = false;
			propertyValue = defaultValue;
		}

		return isPropertyCorrect;
	}
	// *************************************************************************
	bool ConfigFile::CheckDoubleProperty(const string& nameFile, const string& nameProperty, const vector<double> defaultValues, vector<double>& propertyValues, const double& minValue, const double& maxValue)
	{
		bool isPropertyCorrect = true;

		if(defaultValues.size() > 0 && propertyValues.empty())
		{
			ostringstream errorMessage;
			errorMessage<< "Property '"<< nameProperty<< "' empty in file '"<< nameFile<< "'.\nSet default value: "<< defaultValues;
			Output::PrintWarningMessage(errorMessage.str(), true);
			isPropertyCorrect = false;
			propertyValues.resize(defaultValues.size());
			memcpy(&propertyValues[0], &defaultValues[0], defaultValues.size() * sizeof(double));
		}

		for(unsigned int i=0;i<propertyValues.size() && isPropertyCorrect; i++)
		{
			if(minValue > -DBL_MAX && propertyValues[i] < minValue)
			{
				ostringstream errorMessage;
				errorMessage<< "Property '"<< nameProperty<< "' not correct in file '"<< nameFile<< "'.\nSet default value: "<< defaultValues;
				Output::PrintWarningMessage(errorMessage.str(), true);
				isPropertyCorrect = false;
				propertyValues.resize(defaultValues.size());
				memcpy(&propertyValues[0], &defaultValues[0], defaultValues.size() * sizeof(double));
			}

			if(maxValue < DBL_MAX && propertyValues[i] > maxValue)
			{
				ostringstream errorMessage;
				errorMessage<< "Property '"<< nameProperty<< "' not correct in file '"<< nameFile<< "'.\nSet default value: "<< defaultValues;
				Output::PrintWarningMessage(errorMessage.str(), true);
				isPropertyCorrect = false;
				propertyValues.resize(defaultValues.size());
				memcpy(&propertyValues[0], &defaultValues[0], defaultValues.size() * sizeof(double));
			}
		}

		return isPropertyCorrect;
	}
	// *************************************************************************
	void ConfigFile::InitializeConfigurations(const string& configurationFolder, int argc, char** argv)
	{
		string nameFolder = configurationFolder + "Configurations/";
		Output::CreateFolder(nameFolder);

		inputParam.ReadConfiguration("Input.ini", nameFolder, argc, argv);
    cGParam.ReadConfiguration("ConjugateGradient.ini", nameFolder, argc, argv);
    geomParam.ReadConfiguration("Geometry.ini", nameFolder, argc, argv);
		discParam.ReadConfiguration("Discretization.ini", nameFolder, argc, argv);
		paralParam.ReadConfiguration("Parallel.ini", nameFolder, argc, argv);
		graphicParam.ReadConfiguration("Graphic.ini", nameFolder, argc, argv);
		outputParam.ReadConfiguration("Output.ini", nameFolder, argc, argv);
		debugParam.ReadConfiguration("Debug.ini", nameFolder, argc, argv);
		filesParam.ReadConfiguration("FileNames.ini", nameFolder, argc, argv);
		logParam.ReadConfiguration("Log.ini", nameFolder, argc, argv);
		modelParam.ReadConfiguration("Model.ini", nameFolder, argc, argv);
	}
	// *************************************************************************
	InputParameters::InputParameters()
	{
		// Set Properties Default Values
		inputFolder = "file";
		inputData = "Frac6";
		runId = "Run1";
		computeConnectivity = true;
		discretizeDFN = 2;
		assembleSystem = true;
		computeSolution = true;
		activateProfiler = false;
		computeBoundaryConditions = 0;
	}
	// *************************************************************************
	void InputParameters::ReadConfiguration(const string& nameFile, const string& nameFolder, int argc, char** argv)
	{
		if(!ConfigFile::ConfigurationFileExists(nameFile, nameFolder, "Input Parameters"))
			return;

		ConfigFile::ReadStringProperty("InputFolder", inputFolder, inputFolder, nameFile, nameFolder, "Folder where all inputs are taken");
    ConfigFile::FindStringProperty("InputFolder", inputFolder, inputFolder, argc, argv);

    ConfigFile::ReadStringProperty("InputData", inputData, inputData, nameFile, nameFolder, "Input data folder name");
    ConfigFile::FindStringProperty("InputData", inputData, inputData, argc, argv);

    ConfigFile::ReadStringProperty("RunId", runId, runId, nameFile, nameFolder, "Id of the program run");
    ConfigFile::FindStringProperty("RunId", runId, runId, argc, argv);

    ConfigFile::ReadBoolProperty("ComputeConnectivity", computeConnectivity, computeConnectivity, nameFile, nameFolder, "Compute DFN connectivity");
    ConfigFile::FindBoolProperty("ComputeConnectivity", computeConnectivity, computeConnectivity, argc, argv);

    ConfigFile::ReadIntProperty("DiscretizeDFN", discretizeDFN, discretizeDFN, nameFile, nameFolder, "Discretize the DFN\n- 0 No Discretization\n- 1 Only Fracture Discretization\n- 2 Both Fracture and Trace Discretization\n Active only if ComputeConnectivity is active");
		ConfigFile::FindIntProperty("DiscretizeDFN", discretizeDFN, discretizeDFN, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "DiscretizeDFN", 2, discretizeDFN, 0, 2);

		ConfigFile::ReadBoolProperty("AssembleSystem", assembleSystem, assembleSystem, nameFile, nameFolder, "Assemble the system\n Active only if DiscretizeDFN is active");
    ConfigFile::FindBoolProperty("AssembleSystem", assembleSystem, assembleSystem, argc, argv);

		ConfigFile::ReadBoolProperty("ComputeSolution", computeSolution, computeSolution, nameFile, nameFolder, "Compute the solution of the system\n Active only if AssembleSystem is active");
    ConfigFile::FindBoolProperty("ComputeSolution", computeSolution, computeSolution, argc, argv);

    ConfigFile::ReadBoolProperty("ActivateProfiler", activateProfiler, activateProfiler, nameFile, nameFolder, "Compute time/memory/CPU for program profiling");
    ConfigFile::FindBoolProperty("ActivateProfiler", activateProfiler, activateProfiler, argc, argv);

    ConfigFile::ReadIntProperty("ComputeBoundaryConditions", computeBoundaryConditions, computeBoundaryConditions, nameFile, nameFolder, "Used to compute the boundary conditions of DFN using tess fractures in .fab file\n 0 - No computation\n 1 - Computation and remove tess fracture\n 2 - Computation and mantain tess fracture");
		ConfigFile::FindIntProperty("ComputeBoundaryConditions", computeBoundaryConditions, computeBoundaryConditions, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "ComputeBoundaryConditions", 0, computeBoundaryConditions, 0);
	}
	// *************************************************************************
	CGParameters::CGParameters()
	{
		// Set Properties Default Values
		maxIterations.push_back(1000);
		convergence.push_back(1e-12);
		alpha = 1.0;
		penaltyFactorDirichlet = 1.0;
		domainSize = 1.0;
		diffusionTermAmplificationType = 2;
		diffusionTermAmplificationConstant = 0.0;
		traceLengthScaling = false;
    preconditioner = 0;
    restart = false;
    refine = 0;
	}
	// *************************************************************************
	void CGParameters::ReadConfiguration(const string& nameFile, const string& nameFolder, int argc, char** argv)
	{
		if(!ConfigFile::ConfigurationFileExists(nameFile, nameFolder, "Conjugate Gradient Parameters"))
			return;

		vector<int> defaultMaxIterations(1, 1000);
		ConfigFile::ReadIntProperty("MaxIterations", maxIterations, maxIterations, nameFile, nameFolder, "Max iterations");
		ConfigFile::FindIntProperty("MaxIterations", maxIterations, maxIterations, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "MaxIterations", defaultMaxIterations, maxIterations, 1);

		vector<double> defaultConvergence(1, 1e-12);
		ConfigFile::ReadDoubleProperty("Convergence", convergence, convergence, nameFile, nameFolder, "Convergence");
		ConfigFile::FindDoubleProperty("Convergence", convergence, convergence, argc, argv);
		ConfigFile::CheckDoubleProperty(nameFile, "Convergence", defaultConvergence, convergence, 1e-16, 1e2);

    ConfigFile::ReadDoubleProperty("Alpha", alpha, alpha, nameFile, nameFolder, "Alpha");
    ConfigFile::FindDoubleProperty("Alpha", alpha, alpha, argc, argv);
    ConfigFile::CheckDoubleProperty(nameFile, "Alpha", 1.0, alpha, 0);

    ConfigFile::ReadDoubleProperty("PenaltyFactorDirichlet", penaltyFactorDirichlet, penaltyFactorDirichlet, nameFile, nameFolder, "Penalty Factor used in the flow computation for boundary traces");
    ConfigFile::FindDoubleProperty("PenaltyFactorDirichlet", penaltyFactorDirichlet, penaltyFactorDirichlet, argc, argv);
    ConfigFile::CheckDoubleProperty(nameFile, "PenaltyFactorDirichlet", 1.0, penaltyFactorDirichlet, 0);

    ConfigFile::ReadDoubleProperty("DomainSize", domainSize, domainSize, nameFile, nameFolder, "Qualitative size of the domain of the DFN (Used for the DiffusionTermAmplification)");
    ConfigFile::FindDoubleProperty("DomainSize", domainSize, domainSize, argc, argv);
    ConfigFile::CheckDoubleProperty(nameFile, "DomainSize", 1.0, domainSize, 0.0);

    ConfigFile::ReadIntProperty("DiffusionTermAmplificationType", diffusionTermAmplificationType, diffusionTermAmplificationType, nameFile, nameFolder, "Amplification type for the diffusion term\n- 0 No amplification\n- 1 Amplification with DomainSize/minK\n- 2 Amplification with DomainSize/meanK\n- 3 Amplification with DiffusionTermAmplificationConstant");
		ConfigFile::FindIntProperty("DiffusionTermAmplificationType", diffusionTermAmplificationType, diffusionTermAmplificationType, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "DiffusionTermAmplificationType", 2, diffusionTermAmplificationType, 0, 3);

    ConfigFile::ReadDoubleProperty("DiffusionTermAmplificationConstant", diffusionTermAmplificationConstant, diffusionTermAmplificationConstant, nameFile, nameFolder, "Constant amplification for diffusion term.\nUsed only if DiffusionTermAmplificationType is 3");
    ConfigFile::FindDoubleProperty("DiffusionTermAmplificationConstant", diffusionTermAmplificationConstant, diffusionTermAmplificationConstant, argc, argv);
    ConfigFile::CheckDoubleProperty(nameFile, "DiffusionTermAmplificationConstant", 0.0, diffusionTermAmplificationConstant, 0.0);

		ConfigFile::ReadBoolProperty("TraceLengthScaling", traceLengthScaling, traceLengthScaling, nameFile, nameFolder, "Activate the scaling factor by traces length");
		ConfigFile::FindBoolProperty("TraceLengthScaling", traceLengthScaling, traceLengthScaling, argc, argv);

		ConfigFile::ReadIntProperty("Precond", preconditioner, preconditioner, nameFile, nameFolder, "Activate the preconditioner.\n- 0 No Preconditioner\n- 1 Diagonal Precoditioner\n- 2 Block Diagonal Preconditioner");
		ConfigFile::FindIntProperty("Precond", preconditioner, preconditioner, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "Precond", 0, preconditioner, 0, 2);

		ConfigFile::ReadBoolProperty("Restart", restart, restart, nameFile, nameFolder, "Restart the CG from a previous state [Unused by now]");
		ConfigFile::FindBoolProperty("Restart", restart, restart, argc, argv);

		ConfigFile::ReadIntProperty("Refine", refine, refine, nameFile, nameFolder, "Used to refine the solution every Refine iterations");
		ConfigFile::FindIntProperty("Refine", refine, refine, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "Refine", 0, refine, 0);
	}
	// *************************************************************************
	GeometricParameters::GeometricParameters()
	{
		// Set Properties Default Values
		numEdgesRoundedFractures = 0;
		fractureIntersectionTolerance = 1e-12;
		fractureRotationTolerance = 1e-5;
		fractureCellTolerance = 1e-12;
		traceFractureEdgeTolerance = 1e-8;
		boundaryTraceFractureEdgeTolerance = 1e-8;
		pointCorrespondenceTolerance = 1e-08;
		lineParallelismTolerance = 1e-7;
		smoothingTolerance = 1e-12;
	}
	// *************************************************************************
	void GeometricParameters::ReadConfiguration(const string& nameFile, const string& nameFolder, int argc, char** argv)
	{
		if(!ConfigFile::ConfigurationFileExists(nameFile, nameFolder, "Geometry Parameters"))
			return;

		ConfigFile::ReadIntProperty("NumEdgesRoundedFractures", numEdgesRoundedFractures, numEdgesRoundedFractures, nameFile, nameFolder, "Number of edges of the rounded fractures.\nIf 0 then no rounded algorithm is applyed");
		ConfigFile::FindIntProperty("NumEdgesRoundedFractures", numEdgesRoundedFractures, numEdgesRoundedFractures, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "NumEdgesRoundedFractures", 0, numEdgesRoundedFractures, 0);

		ConfigFile::ReadDoubleProperty("FractureIntersectionTolerance", fractureIntersectionTolerance, fractureIntersectionTolerance, nameFile, nameFolder, "Tolerance used for trace research");
    ConfigFile::FindDoubleProperty("FractureIntersectionTolerance", fractureIntersectionTolerance, fractureIntersectionTolerance, argc, argv);
    ConfigFile::CheckDoubleProperty(nameFile, "FractureIntersectionTolerance", 1e-12, fractureIntersectionTolerance, 1e-16, 1e2);

		ConfigFile::ReadDoubleProperty("FractureRotationTolerance", fractureRotationTolerance, fractureRotationTolerance, nameFile, nameFolder, "Tolerance used for 2D points on fracture rotation");
    ConfigFile::FindDoubleProperty("FractureRotationTolerance", fractureRotationTolerance, fractureRotationTolerance, argc, argv);
    ConfigFile::CheckDoubleProperty(nameFile, "FractureRotationTolerance", 1e-5, fractureRotationTolerance, 1e-16, 1e2);

		ConfigFile::ReadDoubleProperty("FractureCellTolerance", fractureCellTolerance, fractureCellTolerance, nameFile, nameFolder, "Tolerance used to search points inside a fracture cell");
    ConfigFile::FindDoubleProperty("FractureCellTolerance", fractureCellTolerance, fractureCellTolerance, argc, argv);
    ConfigFile::CheckDoubleProperty(nameFile, "FractureCellTolerance", 1e-12, fractureCellTolerance, 1e-16, 1e2);

    ConfigFile::ReadDoubleProperty("TraceFractureEdgeTolerance", traceFractureEdgeTolerance, traceFractureEdgeTolerance, nameFile, nameFolder, "Tolerance used to search trace intersection with fracture edge");
    ConfigFile::FindDoubleProperty("TraceFractureEdgeTolerance", traceFractureEdgeTolerance, traceFractureEdgeTolerance, argc, argv);
    ConfigFile::CheckDoubleProperty(nameFile, "TraceFractureEdgeTolerance", 1e-8, traceFractureEdgeTolerance, 1e-16, 1e2);

    ConfigFile::ReadDoubleProperty("BoundaryTraceFractureEdgeTolerance", boundaryTraceFractureEdgeTolerance, boundaryTraceFractureEdgeTolerance, nameFile, nameFolder, "Tolerance used to search boundary trace intersection with fracture edge");
    ConfigFile::FindDoubleProperty("BoundaryTraceFractureEdgeTolerance", boundaryTraceFractureEdgeTolerance, boundaryTraceFractureEdgeTolerance, argc, argv);
    ConfigFile::CheckDoubleProperty(nameFile, "BoundaryTraceFractureEdgeTolerance", 1e-8, boundaryTraceFractureEdgeTolerance, 1e-16, 1e2);

		ConfigFile::ReadDoubleProperty("PointCorrespondenceTolerance", pointCorrespondenceTolerance, pointCorrespondenceTolerance, nameFile, nameFolder, "Tolerance used for points identification");
    ConfigFile::FindDoubleProperty("PointCorrespondenceTolerance", pointCorrespondenceTolerance, pointCorrespondenceTolerance, argc, argv);
    ConfigFile::CheckDoubleProperty(nameFile, "pointCorrespondenceTolerance", 1e-8, pointCorrespondenceTolerance, 1e-16, 1e2);

		ConfigFile::ReadDoubleProperty("LineParallelismTolerance", lineParallelismTolerance, lineParallelismTolerance, nameFile, nameFolder, "Tolerance used for straight lines parallelism");
    ConfigFile::FindDoubleProperty("LineParallelismTolerance", lineParallelismTolerance, lineParallelismTolerance, argc, argv);
    ConfigFile::CheckDoubleProperty(nameFile, "LineParallelismTolerance", 1e-7, lineParallelismTolerance, 1e-16, 1e2);

    ConfigFile::ReadDoubleProperty("SmoothingTolerance", smoothingTolerance, smoothingTolerance, nameFile, nameFolder, "Tolerance used for the smoothing process");
    ConfigFile::FindDoubleProperty("SmoothingTolerance", smoothingTolerance, smoothingTolerance, argc, argv);
    ConfigFile::CheckDoubleProperty(nameFile, "SmoothingTolerance", 1e-14, smoothingTolerance, 1e-16, 1e2);
	}
	// *************************************************************************
	DiscretizationParameters::DiscretizationParameters()
	{
		// Set Properties Default Values
		fractureMethodType = 0;
		traceMethodType = 4;
		boundaryTraceMethodType = 0;
		removeTraceEndPoints = 0;
		fracturePolynomialDegree = 1;
		fractureQuadratureOrder = 5;
		fractureMeshParameter.resize(2);
		fractureMeshParameter[0] = 1;
		fractureMeshParameter[1] = 100;
		fractureMeshMinNumberElement = 0;
		fractureMeshMaxNumberElement = 0;
		traceContinuosElements = false;
		tracePolynomialDegree = 1;
		traceQuadratureOrder = 3;
		traceMeshStep = 0.5;
		traceMeshAmplificationConstant = 1.0;
		boundaryTraceMeshStep = 0.5;
		boundaryTracePolynomialDegree = 0;
		boundaryTraceQuadratureOrder = 3;
		smoothFractureMesh = false;
		smoothFractureMeshType = 1;
		fractureMeshRefinementThreshold = 0.0;
		fractureMeshRefinementMaxStep = 6;
		fractureMeshRefinementThresholdType = 0;
		fractureMeshRefinementAreaDivision = 3;
	}
	// *************************************************************************
	void DiscretizationParameters::ReadConfiguration(const string& nameFile, const string& nameFolder, int argc, char** argv)
	{
		if(!ConfigFile::ConfigurationFileExists(nameFile, nameFolder, "Discretization Parameters"))
			return;

		ConfigFile::ReadIntProperty("FractureMethodType", fractureMethodType, fractureMethodType, nameFile, nameFolder, "Method type used for the fracture discretization\n- 0 FemClassic : Classic FEM Method");
		ConfigFile::FindIntProperty("FractureMethodType", fractureMethodType, fractureMethodType, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "FractureMethodType", 0, fractureMethodType, 0, 0);

		ConfigFile::ReadIntProperty("TraceMethodType", traceMethodType, traceMethodType, nameFile, nameFolder, "Method type used for the trace discretization\n- 0 InducedBase : Trace discretization induced by fracture triangulation\n- 1 FemClassic : Trace discretization equispaced\n- 2 FemClassicDep : Trace discretization equispaced fracture dependent (MIN)\n- 3 FemClassicDep : Trace discretization equispaced fracture dependent (MAX)\n- 4 FemClassicDep : Trace discretization equispaced fracture dependent (MEAN)");
		ConfigFile::FindIntProperty("TraceMethodType", traceMethodType, traceMethodType, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "TraceMethodType", 4, traceMethodType, 0, 4);

		ConfigFile::ReadIntProperty("BoundaryTraceMethodType", boundaryTraceMethodType, boundaryTraceMethodType, nameFile, nameFolder, "Method type used for the boundary trace discretization\n- 0 InducedBase : Boundary Trace discretization induced by fracture triangulation\n- 1 FemClassic : Boundary Trace discretization equispaced ");
		ConfigFile::FindIntProperty("BoundaryTraceMethodType", boundaryTraceMethodType, boundaryTraceMethodType, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "BoundaryTraceMethodType", 0, boundaryTraceMethodType, 0, 1);

		ConfigFile::ReadIntProperty("RemoveTraceEndPoints", removeTraceEndPoints, removeTraceEndPoints, nameFile, nameFolder, "Method type used for the trace discretization\n- 0 Not remove end points\n- 1 Remove end points\n- 2 Remove end points and crossing trace points");
		ConfigFile::FindIntProperty("RemoveTraceEndPoints", removeTraceEndPoints, removeTraceEndPoints, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "RemoveTraceEndPoints", 0, removeTraceEndPoints, 0, 2);

		ConfigFile::ReadIntProperty("FracturePolynomialDegree", fracturePolynomialDegree, fracturePolynomialDegree, nameFile, nameFolder, "Number of degrees of freedom on the fracture method");
		ConfigFile::FindIntProperty("FracturePolynomialDegree", fracturePolynomialDegree, fracturePolynomialDegree, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "FracturePolynomialDegree", 1, fracturePolynomialDegree, 1, 3);

		ConfigFile::ReadIntProperty("FractureQuadratureOrder", fractureQuadratureOrder, fractureQuadratureOrder, nameFile, nameFolder, "Order of the quadrature formula on the fracture");
		ConfigFile::FindIntProperty("FractureQuadratureOrder", fractureQuadratureOrder, fractureQuadratureOrder, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "FractureQuadratureOrder", 5, fractureQuadratureOrder, 1, 10);

		vector<double> defaultFractureMeshParameter(2);
		defaultFractureMeshParameter[0] = 1;
		defaultFractureMeshParameter[1] = 100;
		ConfigFile::ReadDoubleProperty("FractureMeshParameter", fractureMeshParameter, fractureMeshParameter, nameFile, nameFolder, "Fracture Mesh Parameter. Its size is 2.\n The first number is the type of the parameter\n  - 0 Area\n  - 1 N. Triangles\n The second is the parameter value");
		ConfigFile::FindDoubleProperty("FractureMeshParameter", fractureMeshParameter, fractureMeshParameter, argc, argv);
		ConfigFile::CheckDoubleProperty(nameFile, "FractureMeshParameter", defaultFractureMeshParameter, fractureMeshParameter, 0);

		ConfigFile::ReadIntProperty("FractureMeshMinNumberElement", fractureMeshMinNumberElement, fractureMeshMinNumberElement, nameFile, nameFolder, "Minimum number of element in the fracture mesh\nnote: inactive if 0");
		ConfigFile::FindIntProperty("FractureMeshMinNumberElement", fractureMeshMinNumberElement, fractureMeshMinNumberElement, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "FractureMeshMinNumberElement", 0, fractureMeshMinNumberElement, 0);

		ConfigFile::ReadIntProperty("FractureMeshMaxNumberElement", fractureMeshMaxNumberElement, fractureMeshMaxNumberElement, nameFile, nameFolder, "Maximum number of element in the fracture mesh\nnote: inactive if 0");
		ConfigFile::FindIntProperty("FractureMeshMaxNumberElement", fractureMeshMaxNumberElement, fractureMeshMaxNumberElement, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "FractureMeshMaxNumberElement", 0, fractureMeshMaxNumberElement, 0);

		ConfigFile::ReadBoolProperty("TraceContinuosElements", traceContinuosElements, traceContinuosElements, nameFile, nameFolder, "Activate trace continous elements");
		ConfigFile::FindBoolProperty("TraceContinuosElements", traceContinuosElements, traceContinuosElements, argc, argv);

		ConfigFile::ReadIntProperty("TracePolynomialDegree", tracePolynomialDegree, tracePolynomialDegree, nameFile, nameFolder, "Number of degrees of freedom on the trace method");
		ConfigFile::FindIntProperty("TracePolynomialDegree", tracePolynomialDegree, tracePolynomialDegree, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "TracePolynomialDegree", 0, tracePolynomialDegree, 0, 3);

		ConfigFile::ReadIntProperty("TraceQuadratureOrder", traceQuadratureOrder, traceQuadratureOrder, nameFile, nameFolder, "Order of the quadrature formula on the trace");
		ConfigFile::FindIntProperty("TraceQuadratureOrder", traceQuadratureOrder, traceQuadratureOrder, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "TraceQuadratureOrder", 3, traceQuadratureOrder, 1, 10);

		ConfigFile::ReadDoubleProperty("TraceMeshStep", traceMeshStep, traceMeshStep, nameFile, nameFolder, "Mesh step of trace mesh.\nUse only if TraceMethodType is higher than 1");
		ConfigFile::FindDoubleProperty("TraceMeshStep", traceMeshStep, traceMeshStep, argc, argv);
		ConfigFile::CheckDoubleProperty(nameFile, "TraceMeshStep", 0.5, traceMeshStep, 1e-15);

		ConfigFile::ReadDoubleProperty("TraceMeshAmplificationConstant", traceMeshAmplificationConstant, traceMeshAmplificationConstant, nameFile, nameFolder, "Mesh amplification constant applied to traceMeshStep. Use only if TraceMethodType is higher than 2");
		ConfigFile::FindDoubleProperty("TraceMeshAmplificationConstant", traceMeshAmplificationConstant, traceMeshAmplificationConstant, argc, argv);
		ConfigFile::CheckDoubleProperty(nameFile, "TraceMeshAmplificationConstant", 1.0, traceMeshAmplificationConstant, 1e-15);

		ConfigFile::ReadIntProperty("BoundaryTracePolynomialDegree", boundaryTracePolynomialDegree, boundaryTracePolynomialDegree, nameFile, nameFolder, "Number of degrees of freedom on the boundary trace method");
		ConfigFile::FindIntProperty("BoundaryTracePolynomialDegree", boundaryTracePolynomialDegree, boundaryTracePolynomialDegree, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "BoundaryTracePolynomialDegree", 0, boundaryTracePolynomialDegree, 0, 3);

		ConfigFile::ReadIntProperty("BoundaryTraceQuadratureOrder", boundaryTraceQuadratureOrder, boundaryTraceQuadratureOrder, nameFile, nameFolder, "Order of the quadrature formula on the boundary trace");
		ConfigFile::FindIntProperty("BoundaryTraceQuadratureOrder", boundaryTraceQuadratureOrder, boundaryTraceQuadratureOrder, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "BoundaryTraceQuadratureOrder", 3, boundaryTraceQuadratureOrder, 1, 10);

		ConfigFile::ReadDoubleProperty("BoundaryTraceMeshStep", boundaryTraceMeshStep, boundaryTraceMeshStep, nameFile, nameFolder, "Mesh step of boundary trace mesh.\nUse only  if the BoundaryTraceMethodType is 1");
		ConfigFile::FindDoubleProperty("BoundaryTraceMeshStep", boundaryTraceMeshStep, boundaryTraceMeshStep, argc, argv);
		ConfigFile::CheckDoubleProperty(nameFile, "BoundaryTraceMeshStep", 0.5, boundaryTraceMeshStep, 1e-15);

		ConfigFile::ReadBoolProperty("SmoothFractureMesh", smoothFractureMesh, smoothFractureMesh, nameFile, nameFolder, "Smooth the fractures discretization");
		ConfigFile::FindBoolProperty("SmoothFractureMesh", smoothFractureMesh, smoothFractureMesh, argc, argv);

		ConfigFile::ReadIntProperty("SmoothFractureMeshType", smoothFractureMeshType, smoothFractureMeshType, nameFile, nameFolder, "Method type used to do the smoothing of fracture mesh\n- 0 Classic search (slower)\n- 1 Search using trace direction (faster)");
		ConfigFile::FindIntProperty("SmoothFractureMeshType", smoothFractureMeshType, smoothFractureMeshType, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "SmoothFractureMeshType", 1, smoothFractureMeshType, 0, 1);

		ConfigFile::ReadDoubleProperty("FractureMeshRefinementThreshold", fractureMeshRefinementThreshold, fractureMeshRefinementThreshold, nameFile, nameFolder, "Threshold that activate the refining process");
		ConfigFile::FindDoubleProperty("FractureMeshRefinementThreshold", fractureMeshRefinementThreshold, fractureMeshRefinementThreshold, argc, argv);
		ConfigFile::CheckDoubleProperty(nameFile, "FractureMeshRefinementThreshold", 0.0, fractureMeshRefinementThreshold, 0.0, 1.0);

		ConfigFile::ReadIntProperty("FractureMeshRefinementMaxStep", fractureMeshRefinementMaxStep, fractureMeshRefinementMaxStep, nameFile, nameFolder, "Max number of steps of refining process");
		ConfigFile::FindIntProperty("FractureMeshRefinementMaxStep", fractureMeshRefinementMaxStep, fractureMeshRefinementMaxStep, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "FractureMeshRefinementMaxStep", 6, fractureMeshRefinementMaxStep, 1);

		ConfigFile::ReadIntProperty("FractureMeshRefinementThresholdType", fractureMeshRefinementThresholdType, fractureMeshRefinementThresholdType, nameFile, nameFolder, "Method type used for refinement threshold computation\n- 0 APosterioriErrorSum_0\n- 1 NormSolution_0\n- 2 NormSolution_s");
		ConfigFile::FindIntProperty("FractureMeshRefinementThresholdType", fractureMeshRefinementThresholdType, fractureMeshRefinementThresholdType, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "FractureMeshRefinementThresholdType", 0, fractureMeshRefinementThresholdType, 0, 2);

		ConfigFile::ReadIntProperty("FractureMeshRefinementAreaDivision", fractureMeshRefinementAreaDivision, fractureMeshRefinementAreaDivision, nameFile, nameFolder, "Number of division of fracture area cells in refinement");
		ConfigFile::FindIntProperty("FractureMeshRefinementAreaDivision", fractureMeshRefinementAreaDivision, fractureMeshRefinementAreaDivision, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "FractureMeshRefinementAreaDivision", 3, fractureMeshRefinementAreaDivision, 1);
	}
	// *************************************************************************
	ParallelParameters::ParallelParameters()
	{
		// Set Properties Default Values
		dfnConnectivityDivisionType = 1;
		masterFracturePartition = 0.0;
		partitionType = 1;
		metisDivisionType = 0;
		cudaStreamNumber = 1;
		cudaVersionSolver = 4;
	}
	// *************************************************************************
	void ParallelParameters::ReadConfiguration(const string& nameFile, const string& nameFolder, int argc, char** argv)
	{
		if(!ConfigFile::ConfigurationFileExists(nameFile, nameFolder, "Parallel Parameters"))
			return;

		ConfigFile::ReadIntProperty("DfnConnectivityDivisionType", dfnConnectivityDivisionType, dfnConnectivityDivisionType, nameFile, nameFolder, "Type of division of fractures among processes on DFN Connectivity computation\n- 0 Simple division, each process same amount of fractures\n- 1 Weight division, fractures are divided in processes depending on theirs index");
		ConfigFile::FindIntProperty("DfnConnectivityDivisionType", dfnConnectivityDivisionType, dfnConnectivityDivisionType, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "DfnConnectivityDivisionType", 1, dfnConnectivityDivisionType, 0, 1);

		ConfigFile::ReadDoubleProperty("MasterFracturePartition", masterFracturePartition, masterFracturePartition, nameFile, nameFolder, "Partition weight on Master, from 0.0 to 1.0\n0.0 means that the master will not participate to the partition\n1.0 means that the master participate to partition as the slaves");
    ConfigFile::FindDoubleProperty("MasterFracturePartition", masterFracturePartition, masterFracturePartition, argc, argv);
    ConfigFile::CheckDoubleProperty(nameFile, "MasterFracturePartition", 0.0, masterFracturePartition, 0.0, 1.0);

		ConfigFile::ReadIntProperty("PartitionType", partitionType, partitionType, nameFile, nameFolder, "Partition type\n 0 - Uniform Partition\n 1 - METIS Partition");
    ConfigFile::FindIntProperty("PartitionType", partitionType, partitionType, argc, argv);
    ConfigFile::CheckIntProperty(nameFile, "PartitionType", 1, partitionType, 0, 1);

		ConfigFile::ReadIntProperty("MetisDivisionType", metisDivisionType, metisDivisionType, nameFile, nameFolder, "METIS division type\n 0 - No weight (METIS DEFAULT)\n 1 - Vertex weigth (METIS_OBJTYPE_CUT)\n 2 - Communication weigth (METIS_OBJTYPE_VOL)");
    ConfigFile::FindIntProperty("MetisDivisionType", metisDivisionType, metisDivisionType, argc, argv);
    ConfigFile::CheckIntProperty(nameFile, "MetisDivisionType", 0, metisDivisionType, 0, 2);

		ConfigFile::ReadIntProperty("CudaStreamNumber", cudaStreamNumber, cudaStreamNumber, nameFile, nameFolder, "Number of the Streams in GPU");
		ConfigFile::FindIntProperty("CudaStreamNumber", cudaStreamNumber, cudaStreamNumber, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "CudaStreamNumber", 1, cudaStreamNumber, 1, 4);

		ConfigFile::ReadIntProperty("CudaVersionSolver", cudaVersionSolver, cudaVersionSolver, nameFile, nameFolder, "Version of the Triangular Solver in GPU");
		ConfigFile::FindIntProperty("CudaVersionSolver", cudaVersionSolver, cudaVersionSolver, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "CudaVersionSolver", 1, cudaVersionSolver, 1, 5);
	}
	// *************************************************************************
	GraphicParameters::GraphicParameters()
	{
		// Set Properties Default Values
		plotGlobalDFN = false;
		plotLocalDFN = false;
		plotFractures = false;
		plotTraces = false;
		plotBoundaryTraces = false;
		scaleFactor = 5;
		rotationGraphic.resize(3);
		rotationGraphic[0] = 50.0; rotationGraphic[1] = 60.0; rotationGraphic[2] = 0.0;
	}
	// *************************************************************************
	void GraphicParameters::ReadConfiguration(const string& nameFile, const string& nameFolder, int argc, char** argv)
	{
		if(!ConfigFile::ConfigurationFileExists(nameFile, nameFolder, "Graphic Parameters\nUsed only if ENABLE_GRAPHIC option is 1"))
			return;

		ConfigFile::ReadBoolProperty("PlotGlobalDFN", plotGlobalDFN, plotGlobalDFN, nameFile, nameFolder, "Plot Global DFN");
		ConfigFile::FindBoolProperty("PlotGlobalDFN", plotGlobalDFN, plotGlobalDFN, argc, argv);

		ConfigFile::ReadBoolProperty("PlotLocalDFN", plotLocalDFN, plotLocalDFN, nameFile, nameFolder, "Plot Local DFN");
		ConfigFile::FindBoolProperty("PlotLocalDFN", plotLocalDFN, plotLocalDFN, argc, argv);

		ConfigFile::ReadBoolProperty("PlotFractures", plotFractures, plotFractures, nameFile, nameFolder, "Plot Fractures");
		ConfigFile::FindBoolProperty("PlotFractures", plotFractures, plotFractures, argc, argv);

		vector<int> defaultFracturesToPlot;
		ConfigFile::ReadIntProperty("FracturesToPlot", fracturesToPlot, fracturesToPlot, nameFile, nameFolder, "Filter on fractures to plot.\nActive only if PlotFractures is 1\nIf empty then all fractures are plot.");
		ConfigFile::FindIntProperty("FracturesToPlot", fracturesToPlot, fracturesToPlot, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "FracturesToPlot", defaultFracturesToPlot, fracturesToPlot, 1);

		ConfigFile::ReadBoolProperty("PlotTraces", plotTraces, plotTraces, nameFile, nameFolder, "Plot Traces");
		ConfigFile::FindBoolProperty("PlotTraces", plotTraces, plotTraces, argc, argv);

		vector<int> defaultTracesToPlot;
		ConfigFile::ReadIntProperty("TracesToPlot", tracesToPlot, tracesToPlot, nameFile, nameFolder, "Filter on traces to plot.\nActive only if PlotTraces is 1\nIf empty then all traces are plot.");
		ConfigFile::FindIntProperty("TracesToPlot", tracesToPlot, tracesToPlot, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "TracesToPlot", defaultTracesToPlot, tracesToPlot, 1);

		ConfigFile::ReadBoolProperty("PlotBoundaryTraces", plotBoundaryTraces, plotBoundaryTraces, nameFile, nameFolder, "Plot Boundary Traces");
		ConfigFile::FindBoolProperty("PlotBoundaryTraces", plotBoundaryTraces, plotBoundaryTraces, argc, argv);

		vector<int> defaultBoundaryTracesToPlot;
		ConfigFile::ReadIntProperty("BoundaryTracesToPlot", boundaryTracesToPlot, boundaryTracesToPlot, nameFile, nameFolder, "Filter on boundary traces to plot.\nActive only if PlotBoundaryTraces is 1\nIf empty then all the boundary traces are plot.");
		ConfigFile::FindIntProperty("BoundaryTracesToPlot", boundaryTracesToPlot, boundaryTracesToPlot, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "BoundaryTracesToPlot", defaultBoundaryTracesToPlot, boundaryTracesToPlot, 1);

		ConfigFile::ReadDoubleProperty("ScaleFactor", scaleFactor, scaleFactor, nameFile, nameFolder, "Scale factor of the output jpg graphics");
		ConfigFile::FindDoubleProperty("ScaleFactor", scaleFactor, scaleFactor, argc, argv);
		ConfigFile::CheckDoubleProperty(nameFile, "ScaleFactor", 5, scaleFactor, 0.1);

		vector<double> defaultRotationGrafic(3);
		defaultRotationGrafic[0] = 50.0; defaultRotationGrafic[1] = 60.0; defaultRotationGrafic[2] = 0.0;
		ConfigFile::ReadDoubleProperty("RotationGraphics", rotationGraphic, rotationGraphic, nameFile, nameFolder, "Angles for graphics rotation");
		ConfigFile::FindDoubleProperty("RotationGraphics", rotationGraphic, rotationGraphic, argc, argv);
		ConfigFile::CheckDoubleProperty(nameFile, "RotationGraphics", defaultRotationGrafic, rotationGraphic);
	}
	// *************************************************************************
	OutputParameters::OutputParameters()
	{
		// Set Properties Default Values
		printSolution = 1;
		exportFormat = 0;
		exportAdjacencyMatrix = 0;
		exportConnectedNetwork = false;
		exportMesh = false;
		exportModel = false;
		exportToParaview = 0;
		exportToParaviewFormat = 1;
		exportToParaviewEachIteration = 1;
		postProcessSolutionMarker = 0;
		exportCudaSystem = false;
	}
	// *************************************************************************
	void OutputParameters::ReadConfiguration(const string& nameFile, const string& nameFolder, int argc, char** argv)
	{
		if(!ConfigFile::ConfigurationFileExists(nameFile, nameFolder, "Output Parameters"))
			return;

    ConfigFile::ReadIntProperty("PrintSolution", printSolution, printSolution, nameFile, nameFolder, "Prints Solution:\n- 0 no print solution\n- 1 print solution summary\n- 2 print solutions");
		ConfigFile::FindIntProperty("PrintSolution", printSolution, printSolution, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "PrintSolution", 1, printSolution, 0, 2);

    vector<int> defaultFracturesToPrintSolution;
		ConfigFile::ReadIntProperty("FracturesPrintSolution", fracturesPrintSolution, fracturesPrintSolution, nameFile, nameFolder, "Filter on fractures to print solution.\nActive only if PrintSolution is 1\nIf empty then all fractures are plot.");
		ConfigFile::FindIntProperty("FracturesPrintSolution", fracturesPrintSolution, fracturesPrintSolution, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "FracturesPrintSolution", defaultFracturesToPrintSolution, fracturesPrintSolution, 1);

    ConfigFile::ReadIntProperty("ExportFormat", exportFormat, exportFormat, nameFile, nameFolder, "Format of the export:\n- 0 no export\n- 1 general info exported\n- 2 extended info exported");
		ConfigFile::FindIntProperty("ExportFormat", exportFormat, exportFormat, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "ExportFormat", 0, exportFormat, 0, 2);

		ConfigFile::ReadIntProperty("ExportAdjacencyMatrix", exportAdjacencyMatrix, exportAdjacencyMatrix, nameFile, nameFolder, "Export Adjacency Matrix:\n- 0 no export\n- 1 global ids\n- 2 main connected component ids");
		ConfigFile::FindIntProperty("ExportAdjacencyMatrix", exportAdjacencyMatrix, exportAdjacencyMatrix, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "ExportAdjacencyMatrix", 0, exportAdjacencyMatrix, 0, 2);

		ConfigFile::ReadBoolProperty("ExportConnectedNetwork", exportConnectedNetwork, exportConnectedNetwork, nameFile, nameFolder, "Export the connected network into folder ./Export");
		ConfigFile::FindBoolProperty("ExportConnectedNetwork", exportConnectedNetwork, exportConnectedNetwork, argc, argv);

		ConfigFile::ReadBoolProperty("ExportMesh", exportMesh, exportMesh, nameFile, nameFolder, "Export mesh into folder ./Export");
		ConfigFile::FindBoolProperty("ExportMesh", exportMesh, exportMesh, argc, argv);

		ConfigFile::ReadBoolProperty("ExportModel", exportModel, exportModel, nameFile, nameFolder, "Export model into folder ./Export");
		ConfigFile::FindBoolProperty("ExportModel", exportModel, exportModel, argc, argv);

		ConfigFile::ReadIntProperty("ExportToParaview", exportToParaview, exportToParaview, nameFile, nameFolder, "Activate export to Paraview:\n- 0 no export\n- 1 geometry\n- 2 point solution\n- 3 communicated data\n- 4 fracture crossing trace type data\n- 5 a posteriori error\n- 6 exact error");
		ConfigFile::FindIntProperty("ExportToParaview", exportToParaview, exportToParaview, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "ExportToParaview", 0, exportToParaview, 0, 6);

		ConfigFile::ReadIntProperty("ExportToParaviewFormat", exportToParaviewFormat, exportToParaviewFormat, nameFile, nameFolder, "Format of the Paraveiew export:\n- 0 XDMF\n- 1 VTU\n- 2 XDMF + VTU");
		ConfigFile::FindIntProperty("ExportToParaviewFormat", exportToParaviewFormat, exportToParaviewFormat, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "ExportToParaviewFormat", 1, exportToParaviewFormat, 0, 2);

		ConfigFile::ReadIntProperty("ExportToParaviewEachIteration", exportToParaviewEachIteration, exportToParaviewEachIteration, nameFile, nameFolder, "Export to paraview et each iteration\nnote: active only in non stationary mode");
		ConfigFile::FindIntProperty("ExportToParaviewEachIteration", exportToParaviewEachIteration, exportToParaviewEachIteration, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "ExportToParaviewEachIteration", 1, exportToParaviewEachIteration, 1);

		ConfigFile::ReadIntProperty("PostProcessSolutionMarker", postProcessSolutionMarker, postProcessSolutionMarker, nameFile, nameFolder, "Neumann marker used to identify the PostProcessSolution:\nnote: 0 no post process");
		ConfigFile::FindIntProperty("PostProcessSolutionMarker", postProcessSolutionMarker, postProcessSolutionMarker, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "PostProcessSolutionMarker", 0, postProcessSolutionMarker, 0);

		ConfigFile::ReadBoolProperty("ExportCudaSystem", exportCudaSystem, exportCudaSystem, nameFile, nameFolder, "Export Cuda System into folder ./Export");
		ConfigFile::FindBoolProperty("ExportCudaSystem", exportCudaSystem, exportCudaSystem, argc, argv);
	}
	// *************************************************************************
	DebugParameters::DebugParameters()
	{
		// Set Properties Default Values
		startingIndexToPrint = 0;
		maxNumElementToPrint = 100;
		printApplicationInfo = false;
		printProcessInfo = false;
		printFractureInfo = false;
		printTraceInfo = false;
		printBoundTraceInfo = false;
		mantainAlteredFiles = false;
		printEachIteration = 1;
	}
	// *************************************************************************
	void DebugParameters::ReadConfiguration(const string& nameFile, const string& nameFolder, int argc, char** argv)
	{
		if(!ConfigFile::ConfigurationFileExists(nameFile, nameFolder, "Debug Parameters\nUsed only if DEBUG option is 1"))
			return;

		ConfigFile::ReadIntProperty("MaxNumElementToPrint", maxNumElementToPrint, maxNumElementToPrint, nameFile, nameFolder, "Max number of elements to print for vectors, lists ...\n 0 means all the element must be printed");
		ConfigFile::FindIntProperty("MaxNumElementToPrint", maxNumElementToPrint, maxNumElementToPrint, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "MaxNumElementToPrint", 100, maxNumElementToPrint, 0);

		ConfigFile::ReadIntProperty("StartingIndexToPrint", startingIndexToPrint, startingIndexToPrint, nameFile, nameFolder, "Index where to start to print vectors, lists ...");
		ConfigFile::FindIntProperty("StartingIndexToPrint", startingIndexToPrint, startingIndexToPrint, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "StartingIndexToPrint", 0, startingIndexToPrint, 0);

		ConfigFile::ReadBoolProperty("PrintApplicationInfo", printApplicationInfo, printApplicationInfo, nameFile, nameFolder, "Print application general information");
		ConfigFile::FindBoolProperty("PrintApplicationInfo", printApplicationInfo, printApplicationInfo, argc, argv);

		ConfigFile::ReadBoolProperty("PrintProcessInfo", printProcessInfo, printProcessInfo, nameFile, nameFolder, "Print process general information");
		ConfigFile::FindBoolProperty("PrintProcessInfo", printProcessInfo, printProcessInfo, argc, argv);

		vector<int> defaultProcessesToPrint;
		ConfigFile::ReadIntProperty("ProcessesToPrint", processesToPrint, processesToPrint, nameFile, nameFolder, "Filter on processes to print.\nActive only if PrintProcessInfo is 1\nIf empty then all processes are plot.");
		ConfigFile::FindIntProperty("ProcessesToPrint", processesToPrint, processesToPrint, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "ProcessesToPrint", defaultProcessesToPrint, processesToPrint, 0);

		ConfigFile::ReadBoolProperty("PrintFractureInfo", printFractureInfo, printFractureInfo, nameFile, nameFolder, "Print fracture general information");
		ConfigFile::FindBoolProperty("PrintFractureInfo", printFractureInfo, printFractureInfo, argc, argv);

		vector<int> defaultFracturesToPrint;
		ConfigFile::ReadIntProperty("FracturesToPrint", fracturesToPrint, fracturesToPrint, nameFile, nameFolder, "Filter on fractures to print.\nActive only if PrintFractureInfo is 1\nIf empty then all fractures are plot.");
		ConfigFile::FindIntProperty("FracturesToPrint", fracturesToPrint, fracturesToPrint, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "FracturesToPrint", defaultFracturesToPrint, fracturesToPrint, 1);

		ConfigFile::ReadBoolProperty("PrintTraceInfo", printTraceInfo, printTraceInfo, nameFile, nameFolder, "Print trace general information");
		ConfigFile::FindBoolProperty("PrintTraceInfo", printTraceInfo, printTraceInfo, argc, argv);

		vector<int> defaultTracesToPrint;
		ConfigFile::ReadIntProperty("TracesToPrint", tracesToPrint, tracesToPrint, nameFile, nameFolder, "Filter on traces to print.\nActive only if PrintTraceInfo is 1\nIf empty then all traces are plot.");
		ConfigFile::FindIntProperty("TracesToPrint", tracesToPrint, tracesToPrint, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "TracesToPrint", defaultTracesToPrint, tracesToPrint, 1);

		ConfigFile::ReadBoolProperty("PrintBoundTraceInfo", printBoundTraceInfo, printBoundTraceInfo, nameFile, nameFolder, "Print boundary trace general information");
		ConfigFile::FindBoolProperty("PrintBoundTraceInfo", printBoundTraceInfo, printBoundTraceInfo, argc, argv);

		vector<int> defaultBoundTracesToPrint;
		ConfigFile::ReadIntProperty("BoundTracesToPrint", boundTracesToPrint, boundTracesToPrint, nameFile, nameFolder, "Filter on boundary traces to print.\nActive only if PrintBoundTraceInfo is 1\nIf empty then all boundary traces are plot.");
		ConfigFile::FindIntProperty("BoundTracesToPrint", boundTracesToPrint, boundTracesToPrint, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "BoundTracesToPrint", defaultBoundTracesToPrint, boundTracesToPrint, 1);

		ConfigFile::ReadBoolProperty("MantainAlteredFiles", mantainAlteredFiles, mantainAlteredFiles, nameFile, nameFolder, "Mantain DFN and BoundaryConditions altered files");
		ConfigFile::FindBoolProperty("MantainAlteredFiles", mantainAlteredFiles, mantainAlteredFiles, argc, argv);

		ConfigFile::ReadIntProperty("PrintEachIteration", printEachIteration, printEachIteration, nameFile, nameFolder, "Print solver iterations every printEachIteration");
		ConfigFile::FindIntProperty("PrintEachIteration", printEachIteration, printEachIteration, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "PrintEachIteration", 1, printEachIteration, 1);
	}
	// *************************************************************************
	FileNameParameters::FileNameParameters()
	{
		// Set Properties Default Values
		dfnInputFormat = 0;
		dfnFileName = ConfigFile::DefaultPropertyValue();
		outputFolder = "Output";
		outputData = ConfigFile::DefaultPropertyValue();
		discretizationFileName = ConfigFile::DefaultPropertyValue();
		boundaryConditionFileName = ConfigFile::DefaultPropertyValue();
		dirichletFileName = ConfigFile::DefaultPropertyValue();
		neumannFileName = ConfigFile::DefaultPropertyValue();
		forcingTermFileName = ConfigFile::DefaultPropertyValue();
		diffusionTermFileName = ConfigFile::DefaultPropertyValue();
		exactSolutionFileName = ConfigFile::DefaultPropertyValue();
		exactSolutionOnTraceFileName = ConfigFile::DefaultPropertyValue();
		dfnStructFile = "general_info.dat";
		dfnDofFile = "dof.dat";
		dfnFracFile = "F";
		triangleMeshOptions = "-Qzpqnea";
		debugFolder = "Debug";
		graphicFolder = "Plots";
		processFolder = "Process";
		fracturesFolder = "Fractures";
		tracesFolder = "Traces";
		boundTracesFolder = "BoundTraces";
		solutionsFolder = "Solutions";
		errorsFile = "Summary.dat";
		iterationFolder = "Solution_it_";
		residualFolder = "Solution_res_";
		finalSolutionFolder = "Solution";
		timeFolder = "Time";
		cGLogFile = "Gradient.log";
		jFile = "J.dat";
		resFile = "Res.dat";
		exportFolder = "Export";
		exportMeshFile = "Mesh";
		importFolder = "./Import";
	}
	// *************************************************************************
	void FileNameParameters::ReadConfiguration(const string& nameFile, const string& nameFolder, int argc, char** argv)
	{
		if(!ConfigFile::ConfigurationFileExists(nameFile, nameFolder, "FileName Parameters"))
			return;

		ConfigFile::ReadIntProperty("DfnInputFormat", dfnInputFormat, dfnInputFormat, nameFile, nameFolder, "Input format for the DFN file\n- 0 .dat\n- 1 .bin");
		ConfigFile::FindIntProperty("DfnInputFormat", dfnInputFormat, dfnInputFormat, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "DfnInputFormat", 0, dfnInputFormat, 0, 1);

		ConfigFile::ReadStringProperty("DfnFileName", dfnFileName, dfnFileName, nameFile, nameFolder, "Input DFNfile name (without extension)\nDefault value is 'DFN' in InputFolder");
    ConfigFile::FindStringProperty("DfnFileName", dfnFileName, dfnFileName, argc, argv);

		ConfigFile::ReadStringProperty("OutputFolder", outputFolder, outputFolder, nameFile, nameFolder, "Output folder name");
    ConfigFile::FindStringProperty("OutputFolder", outputFolder, outputFolder, argc, argv);

    ConfigFile::ReadStringProperty("OutputData", outputData, outputData, nameFile, nameFolder, " Output folder data name\nDefault value is InputData value");
    ConfigFile::FindStringProperty("OutputData", outputData, outputData, argc, argv);

		ConfigFile::ReadStringProperty("DiscretizationFileName", discretizationFileName, discretizationFileName, nameFile, nameFolder, "Input discretization condition file name\nDefault value is 'Discretization' in InputFolder");
    ConfigFile::FindStringProperty("DiscretizationFileName", discretizationFileName, discretizationFileName, argc, argv);

    ConfigFile::ReadStringProperty("BoundaryConditionFileName", boundaryConditionFileName, boundaryConditionFileName, nameFile, nameFolder, "Input boundary condition file name\nDefault value is 'BoundaryConditions' in InputFolder");
    ConfigFile::FindStringProperty("BoundaryConditionFileName", boundaryConditionFileName, boundaryConditionFileName, argc, argv);

    ConfigFile::ReadStringProperty("DirichletFileName", dirichletFileName, dirichletFileName, nameFile, nameFolder, "Input Dirichlet file name\nDefault value is 'Dirichlet' in InputFolder");
    ConfigFile::FindStringProperty("DirichletFileName", dirichletFileName, dirichletFileName, argc, argv);

    ConfigFile::ReadStringProperty("NeumannFileName", neumannFileName, neumannFileName, nameFile, nameFolder, "Input Neumann file name\nDefault value is 'Neumann' in InputFolder");
    ConfigFile::FindStringProperty("NeumannFileName", neumannFileName, neumannFileName, argc, argv);

    ConfigFile::ReadStringProperty("ForcingTermFileName", forcingTermFileName, forcingTermFileName, nameFile, nameFolder, "Input forcing term file name\nDefault value is 'ForcingTerm' in InputFolder");
    ConfigFile::FindStringProperty("ForcingTermFileName", forcingTermFileName, forcingTermFileName, argc, argv);

    ConfigFile::ReadStringProperty("DiffusionTermFileName", diffusionTermFileName, diffusionTermFileName, nameFile, nameFolder, "Input diffusion term file name\nDefault value is 'DiffusionTerm' in InputFolder");
    ConfigFile::FindStringProperty("DiffusionTermFileName", diffusionTermFileName, diffusionTermFileName, argc, argv);

		ConfigFile::ReadStringProperty("ExactSolutionFileName", exactSolutionFileName, exactSolutionFileName, nameFile, nameFolder, "Input exact solution file name\nDefault value is 'ExactSolution' in InputFolder");
    ConfigFile::FindStringProperty("ExactSolutionFileName", exactSolutionFileName, exactSolutionFileName, argc, argv);

		ConfigFile::ReadStringProperty("ExactSolutionOnTraceFileName", exactSolutionOnTraceFileName, exactSolutionOnTraceFileName, nameFile, nameFolder, "Input exact solution on Trace file name\nDefault value is 'ExactSolutionTrace' in InputFolder");
		ConfigFile::FindStringProperty("ExactSolutionOnTraceFileName", exactSolutionOnTraceFileName, exactSolutionOnTraceFileName, argc, argv);

		ConfigFile::ReadStringProperty("DfnStructFile", dfnStructFile, dfnStructFile, nameFile, nameFolder, "Input File containing the DFN structure.\nUsed in Old program");
    ConfigFile::FindStringProperty("DfnStructFile", dfnStructFile, dfnStructFile, argc, argv);

    ConfigFile::ReadStringProperty("DfnDofFile", dfnDofFile, dfnDofFile, nameFile, nameFolder, "Input File containing the DFN dofs.\nUsed on in Old program");
    ConfigFile::FindStringProperty("DfnDofFile", dfnDofFile, dfnDofFile, argc, argv);

    ConfigFile::ReadStringProperty("DfnFracFile", dfnFracFile, dfnFracFile, nameFile, nameFolder, "Input File containing the fracture informations.\nUsed in Old program");
    ConfigFile::FindStringProperty("DfnFracFile", dfnFracFile, dfnFracFile, argc, argv);

    ConfigFile::ReadStringProperty("TriangleMeshOptions", triangleMeshOptions, triangleMeshOptions, nameFile, nameFolder, "Option string for Triangle mesh generator");
		ConfigFile::FindStringProperty("TriangleMeshOptions", triangleMeshOptions, triangleMeshOptions, argc, argv);

		ConfigFile::ReadStringProperty("DebugFolder", debugFolder, debugFolder, nameFile, nameFolder, "Folder name where debug outputs are put");
    ConfigFile::FindStringProperty("DebugFolder", debugFolder, debugFolder, argc, argv);

    ConfigFile::ReadStringProperty("GraphicFolder", graphicFolder, graphicFolder, nameFile, nameFolder, "Folder name where graphic outputs are put");
    ConfigFile::FindStringProperty("GraphicFolder", graphicFolder, graphicFolder, argc, argv);

    ConfigFile::ReadStringProperty("ProcessFolder", processFolder, processFolder, nameFile, nameFolder, "Folder name where process data are put");
    ConfigFile::FindStringProperty("processFolder", processFolder, processFolder, argc, argv);

		ConfigFile::ReadStringProperty("FracturesFolder", fracturesFolder, fracturesFolder, nameFile, nameFolder, "Folder name where fractures data are put");
    ConfigFile::FindStringProperty("FracturesFolder", fracturesFolder, fracturesFolder, argc, argv);

    ConfigFile::ReadStringProperty("TracesFolder", tracesFolder, tracesFolder, nameFile, nameFolder, "Folder name where traces data are put");
    ConfigFile::FindStringProperty("TracesFolder", tracesFolder, tracesFolder, argc, argv);

    ConfigFile::ReadStringProperty("BoundTracesFolder", boundTracesFolder, boundTracesFolder, nameFile, nameFolder, "Folder name where boundary traces data are put");
    ConfigFile::FindStringProperty("BoundTracesFolder", boundTracesFolder, boundTracesFolder, argc, argv);

    ConfigFile::ReadStringProperty("SolutionsFolder", solutionsFolder, solutionsFolder, nameFile, nameFolder, "Folder name where solutions files of the program are put");
    ConfigFile::FindStringProperty("SolutionsFolder", solutionsFolder, solutionsFolder, argc, argv);

		ConfigFile::ReadStringProperty("ErrorsFile", errorsFile, errorsFile, nameFile, nameFolder, "File name of Solution Summary");
    ConfigFile::FindStringProperty("ErrorsFile", errorsFile, errorsFile, argc, argv);

    ConfigFile::ReadStringProperty("ResidualFolder", residualFolder, residualFolder, nameFile, nameFolder, "Folder name where solution outputs for residuals are put");
    ConfigFile::FindStringProperty("ResidualFolder", residualFolder, residualFolder, argc, argv);

    ConfigFile::ReadStringProperty("FinalSolutionFolder", finalSolutionFolder, finalSolutionFolder, nameFile, nameFolder, "Folder name where last solution outputs are put");
    ConfigFile::FindStringProperty("FinalSolutionFolder", finalSolutionFolder, finalSolutionFolder, argc, argv);

    ConfigFile::ReadStringProperty("TimeFolder", timeFolder, timeFolder, nameFile, nameFolder, "Folder name where time data are put");
    ConfigFile::FindStringProperty("TimeFolder", timeFolder, timeFolder, argc, argv);

    ConfigFile::ReadStringProperty("CGLogFile", cGLogFile, cGLogFile, nameFile, nameFolder, "File name of the CG output");
    ConfigFile::FindStringProperty("CGLogFile", cGLogFile, cGLogFile, argc, argv);

		ConfigFile::ReadStringProperty("JFile", jFile, jFile, nameFile, nameFolder, "File with the values of the functional J in solutions");
    ConfigFile::FindStringProperty("JFile", jFile, jFile, argc, argv);

    ConfigFile::ReadStringProperty("ResFile", resFile, resFile, nameFile, nameFolder, "File with the values of the residuals in solutions");
    ConfigFile::FindStringProperty("ResFile", resFile, resFile, argc, argv);

		ConfigFile::ReadStringProperty("ExportFolder", exportFolder, exportFolder, nameFile, nameFolder, "Folder name where export data are put");
    ConfigFile::FindStringProperty("ExportFolder", exportFolder, exportFolder, argc, argv);

		ConfigFile::ReadStringProperty("ExportMeshFile", exportMeshFile, exportMeshFile, nameFile, nameFolder, "File name for the mesh export");
    ConfigFile::FindStringProperty("ExportMeshFile", exportMeshFile, exportMeshFile, argc, argv);

    ConfigFile::ReadStringProperty("ImportFolder", importFolder, importFolder, nameFile, nameFolder, "Folder name where external input data are taken");
    ConfigFile::FindStringProperty("ImportFolder", importFolder, importFolder, argc, argv);
	}
	// *************************************************************************
	LogParameters::LogParameters()
	{
		logFile = "LogFile.log";
		createAnalysisFile = false;
		analysisFile = "Analysis";
		logMaxFileSize = 2;
		logMaxNumFiles = 5;
		timeFile = "Time";
		memoryFile = "Memory.csv";
		logFolder = "Log";
	}
	// *************************************************************************
	void LogParameters::ReadConfiguration(const string& nameFile, const string& nameFolder, int argc, char** argv)
	{
		if(!ConfigFile::ConfigurationFileExists(nameFile, nameFolder, "Log Parameters"))
			return;

		ConfigFile::ReadStringProperty("LogFile", logFile, logFile, nameFile, nameFolder, "File name of Log");
    ConfigFile::FindStringProperty("LogFile", logFile, logFile, argc, argv);

    ConfigFile::ReadBoolProperty("CreateAnalysisFile", createAnalysisFile, createAnalysisFile, nameFile, nameFolder, "Save Analysis file");
		ConfigFile::FindBoolProperty("CreateAnalysisFile", createAnalysisFile, createAnalysisFile, argc, argv);

    ConfigFile::ReadStringProperty("AnalysisFile", analysisFile, analysisFile, nameFile, nameFolder, "File name of Summary");
    ConfigFile::FindStringProperty("AnalysisFile", analysisFile, analysisFile, argc, argv);

		ConfigFile::ReadIntProperty("LogMaxFileSize", logMaxFileSize, logMaxFileSize, nameFile, nameFolder, "Max file log size (MB)");
		ConfigFile::FindIntProperty("LogMaxFileSize", logMaxFileSize, logMaxFileSize, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "LogMaxFileSize", 2, logMaxFileSize, 1);

		ConfigFile::ReadIntProperty("LogMaxNumFiles", logMaxNumFiles, logMaxNumFiles, nameFile, nameFolder, "Max number of backup log files to maintain");
		ConfigFile::FindIntProperty("LogMaxNumFiles", logMaxNumFiles, logMaxNumFiles, argc, argv);
		ConfigFile::CheckIntProperty(nameFile, "LogMaxNumFiles", 5, logMaxNumFiles, 1);

		ConfigFile::ReadStringProperty("TimeFile", timeFile, timeFile, nameFile, nameFolder, "File where time profiling data are stored");
    ConfigFile::FindStringProperty("TimeFile", timeFile, timeFile, argc, argv);

    ConfigFile::ReadStringProperty("MemoryFile", memoryFile, memoryFile, nameFile, nameFolder, "File where memory profiling data are stored");
    ConfigFile::FindStringProperty("MemoryFile", memoryFile, memoryFile, argc, argv);

    ConfigFile::ReadStringProperty("LogFolder", logFolder, logFolder, nameFile, nameFolder, "Folder name where log files of the program are put");
    ConfigFile::FindStringProperty("LogFolder", logFolder, logFolder, argc, argv);

	}
	// *************************************************************************
	ModelParameters::ModelParameters()
	{
		maxTime = 0.0;
		timeStep = 0.0;
		steadyState = false;
	}
	// *************************************************************************
	void ModelParameters::ReadConfiguration(const string& nameFile, const string& nameFolder, int argc, char** argv)
	{
		if(!ConfigFile::ConfigurationFileExists(nameFile, nameFolder, "Model Parameters"))
			return;

		ConfigFile::ReadDoubleProperty("MaxTime", maxTime, maxTime, nameFile, nameFolder, "Maximum time for non stationary problem\n0 value means stationary problem");
    ConfigFile::FindDoubleProperty("MaxTime", maxTime, maxTime, argc, argv);
    ConfigFile::CheckDoubleProperty(nameFile, "MaxTime", 0.0, maxTime, 0.0);

    ConfigFile::ReadDoubleProperty("TimeStep", timeStep, timeStep, nameFile, nameFolder, "Time step used for the non stationary problem discretization\n0 value means stationary problem");
    ConfigFile::FindDoubleProperty("TimeStep", timeStep, timeStep, argc, argv);
    ConfigFile::CheckDoubleProperty(nameFile, "TimeStep", 0.0, timeStep, 0.0);

    ConfigFile::ReadBoolProperty("SteadyState", steadyState, steadyState, nameFile, nameFolder, "Computes the steady state of the time evolving problem");
		ConfigFile::FindBoolProperty("SteadyState", steadyState, steadyState, argc, argv);
	}
	// *************************************************************************
}
