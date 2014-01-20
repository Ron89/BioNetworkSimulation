#ifndef __BCNetwork_H__
#define __BCNetwork_H__

#include<iostream>
#include<cstdlib>
#include<vector>
#include"basicDef.h"

using namespace std;

template <typename compType, typename rateType, typename powerType >
class BCNetwork
//this class is designed to be used as parent class for different iteration schemes or
//network processing class creation.
//In the class: 1) the network is stored reaction-wisely;
// 				2) basic reset(forcing each reactant to return to their original value)
// 					and saving(save the current time and amount of each reactant into
// 					file) are defined;
// 				3) basic parameter/reactant amount changing methods are provided.
// 					But as it is still inconvenient to use, the whole class is defined
// 					as public so that user may freely change each value they think
// 					necessary.
{
//data saving related
private:
	bool networkSpecified;
	ofstream resultFile; 	//file pointer for storying data
// create and erase network
	void eraseNetwork();
	void createNetwork(int numComponent, int numRate);

public:
//core module
	double t; 				//current timepoint
	int nComp; 				//number of time dependent variables
	compType * comp;
	compType * compBackup; 	//vector to store variables
	int nRate; 				//number of reaction rates 
	rateType * rate;//vector to store reaction rates 
	powerType * rateMatrix; 	//matrix to store dependence relation
												//(left side of reaction i)
	powerType * updateMatrix; 	//matrix to store update relation
												//(net change after reaction i happens)

	void errorWarning(int errorIdentifier); 	//error message

	int loadParameter(vector<compType> & initComp, vector<rateType> & inputRate,
			powerType * inputRateMatrix, powerType * inputUpdateMatrix);
	//load parameters from a certain model. Note that this module is not supposed to read
	//models directly from file. To load a model from file, a seperate module called
	//modelLoader is defined in file modelLoader.h

//data saving functions
	int fileOpen(string & condition); 	//specificFileName should be specified in condition
	virtual void saveData();
	int fileClose();

//operating methods
	int changeRate(int identifier, rateType specifiedValue);
	int changeComp(int identifier, compType specifiedValue);
	int changeRateMatrix(int rateIdentifier, int compIdentifier, powerType specifiedValue);
	int changeUpdMatrix(int rateIdentifier, int compIdentifier, powerType specifiedValue);

//reset
	virtual void reset()
	// the simplest reset,
	{
		for (int i=0;i<nComp;i++) 	comp[i]=compBackup[i];
		t=0;
	}

	BCNetwork()
	{
		networkSpecified=0;
	}
	BCNetwork(vector<compType> & initComp, vector<rateType> & inputRate,  
			powerType * inputRateMatrix, powerType * inputUpdateMatrix)
	{
		networkSpecified=0;
		loadParameter(initComp, inputRate, inputRateMatrix, inputUpdateMatrix);
	}
	~BCNetwork()
	{
		if (networkSpecified) 	eraseNetwork();
	}
};

template <typename compType, typename rateType, typename powerType >
void BCNetwork<TYPENAME>::eraseNetwork()
{
	if (networkSpecified)
	{
		nComp=0;
		nRate=0;
		delete [] comp;
		delete [] compBackup;
		delete [] rate;
		delete [] rateMatrix;
		delete [] updateMatrix;
		networkSpecified=0;
	}
	else return;
}
template <typename compType, typename rateType, typename powerType >
void BCNetwork<TYPENAME>::createNetwork(int numComponent, int numRate)
{
	if (networkSpecified) 	eraseNetwork();
	nComp=numComponent;
	nRate=numRate;
	comp= new compType [nComp];
	compBackup= new compType [nComp];
	rate= new rateType [nRate];
	rateMatrix= new powerType [nRate*nComp];
	updateMatrix= new powerType [nRate*nComp];
	networkSpecified=1;
}

template <typename compType, typename rateType, typename powerType >
void BCNetwork<TYPENAME>::errorWarning(int errorIdentifier)
{
	switch(errorIdentifier)
	{
		case 11: 	//overflow
			cout<<endl<<"Error: overflow. user assigning a space larger than expected."
				<<endl;
		case 12: 	//pointer point to an unexpected place
			cout<<endl<<"Error: exceed. pointer point to an undefined area."<<endl;
		case 21: 	//unexpected fstream open state warning
			cout<<endl<<"Warning: unpected fstream open state."<<endl;
		case 22:
			cout<<endl<<"Warning: unpected fstream closed state."<<endl;
		
		default:
			cout<<endl<<"Error: undefined error."<<endl;
	}
}

template <typename compType, typename rateType, typename powerType >
int BCNetwork<TYPENAME>::changeRate(int identifier, rateType specifiedValue)
{
	if (identifier>=nRate)
	{
		errorWarning(2);
		return 1;
	}
	rate[identifier]=specifiedValue;	
	return 0;
}

template <typename compType, typename rateType, typename powerType >
int BCNetwork<TYPENAME>::changeComp(int identifier, compType specifiedValue)
{
	if (identifier>=nComp)
	{
		errorWarning(2);
		return 1;
	}
	comp[identifier]=specifiedValue;	
	return 0;
}

template <typename compType, typename rateType, typename powerType >
int BCNetwork<TYPENAME>::changeRateMatrix(int rateIdentifier, int compIdentifier, powerType specifiedValue)
{
	if (rateIdentifier>=nRate||compIdentifier>=nComp)
	{
		errorWarning(12);
		return 1;
	}
	rateMatrix[rateIdentifier*nRate+compIdentifier]=specifiedValue;
}

template <typename compType, typename rateType, typename powerType >
int BCNetwork<TYPENAME>::changeUpdMatrix(int rateIdentifier, int compIdentifier, powerType specifiedValue)
{
	if (rateIdentifier>=nRate||compIdentifier>=nComp)
	{
		errorWarning(12);
		return 1;
	}
	updateMatrix[rateIdentifier*nRate+compIdentifier]=specifiedValue;
}

template <typename compType, typename rateType, typename powerType >
int BCNetwork<TYPENAME>::loadParameter(vector<compType> & initComp, 
		vector<rateType> & inputRate,
		powerType * inputRateMatrix, powerType * inputUpdateMatrix)
{
//	nComp=initComp.size();
//	nRate=inputRate.size();
//
	createNetwork(initComp.size(),inputRate.size());

//	if (nComp>MAXNCOMP||nRate>MAXNRATE) 
//	{
//		errorWarning(11);
//		return 1;
//	}

	for (int i=0;i<nComp;i++) 	compBackup[i]=comp[i]=initComp[i];
	for (int i=0;i<nRate;i++) 	rate[i]=inputRate[i];
	for (int i=0;i<nComp*nRate;i++) 	
	{
		rateMatrix[i]=inputRateMatrix[i];
		updateMatrix[i]=inputUpdateMatrix[i];
	}
	return 0;
}

template <typename compType, typename rateType, typename powerType >
int BCNetwork<TYPENAME>::fileOpen(string & condition)
{
	if(resultFile.is_open())
	{
		errorWarning(21);
		resultFile.close();
		return 1;
	}
	string fileName(RESULTFOLDER);
	fileName+=RESULTNAME;
	fileName+=condition;
	resultFile.open(fileName.c_str());
	return 0;
}

template <typename compType, typename rateType, typename powerType >
int BCNetwork<TYPENAME>::fileClose()
{
	if(!resultFile.is_open())
	{
		errorWarning(22);
		return 1;
	}
	resultFile.close();
	return 0;
}

template <typename compType, typename rateType, typename powerType >
void BCNetwork<TYPENAME>::saveData()
{
	if(!resultFile.is_open())
	{
		errorWarning(22);
		string condition("Emergency");
		fileOpen(condition);
	}
	resultFile<<t<<'\t';
	for (int i=0;i<nComp;i++) 	resultFile<<comp[i]<<'\t';
	resultFile<<endl;
}

#endif 	// __BCNetwork_H__ defined
