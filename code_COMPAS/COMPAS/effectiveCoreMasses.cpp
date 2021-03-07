//--------------------- EFFECTIVE CORE MASS READER ----------------------//

//-- Isobel Romero-Shaw                                   23.06.2017

//-- This file contains a function capable of providing the effective 
//-- core mass of a MIST star. (.cpp)

#include "effectiveCoreMasses.h"

bool updateEffectiveCoreMass(	 vector<double> *effCoreMass,
				 vector<double> *effCoreType,
				 std::string file_location,
				 double M,
				 double ZdivH
			     )
{
	//-- Define the name of the file from which to read
	string sign = "some string";
	
	if (ZdivH != 0.0)
	{
        sign = "m";
	}
   	else
	{
        sign = "p";
	}

	string stringM, stringZ;

	ostringstream convertM, convertZ;

	convertM << fixed << setprecision(2) << M;
	convertZ << fixed << setprecision(2) << abs(ZdivH);

	stringM = convertM.str();
	stringZ = convertZ.str();

	string Filename = 
		"effCoreMass_" + stringM + "M_" + sign + stringZ + ".dat";

	//-- Read from the file
	string totalPath = file_location + Filename;
	ifstream mass_in(totalPath.c_str());
	string de = "\t\t\t\t";

	double effCM, effCT;

	if (mass_in.fail())
	{
		cout << "Could not find file " +  totalPath << endl;
		return false;
	}
	
	if (mass_in.is_open())
	{
		string Head1, Head2;

		mass_in >> Head1 >> Head2;

		cout << Head1 << "\t" << Head2 << endl;

		for(int linecount = 1; !mass_in.eof() && mass_in.good(); linecount++)
		{
			mass_in >> effCM >> effCT;
			effCoreMass->push_back(effCM);
			effCoreType->push_back(effCT);
			
		}

	}

 	return true;
}

