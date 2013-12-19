#ifndef __timeDerivative_H__
#define __timeDerivative_H__

#include<iostream>
#include<cstdlib>
#include<vector>
#include"basicDef.h"

using namespace std;

template <typename compType, typename rateType, typename powerType >
class timeDerivative
{
public:
//data saving related
	ofstream resultFile;

//core module
	double t;
	int nComp;
	compType comp[MAXNCOMP], compBackup[MAXNCOMP];
	int nRate;
	rateType rate[MAXNRATE]; 
	powerType rateMatrix[MAXNCOMP*MAXNRATE];
	powerType updateMatrix[MAXNCOMP*MAXNRATE];

	void errorWarning(int errorIdentifier);

	int loadParameter(vector<compType> & initComp, vector<rateType> & inputRate,
			powerType * inputRateMatrix, powerType * inputUpdateMatrix);

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
	inline void reset()
	{
		for (int i=0;i<nComp;i++) 	comp[i]=compBackup[i];
		t=0;
	}
};

template <typename compType, typename rateType, typename powerType >
void timeDerivative<TYPENAME>::errorWarning(int errorIdentifier)
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
int timeDerivative<TYPENAME>::changeRate(int identifier, rateType specifiedValue)
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
int timeDerivative<TYPENAME>::changeComp(int identifier, compType specifiedValue)
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
int timeDerivative<TYPENAME>::changeRateMatrix(int rateIdentifier, int compIdentifier, powerType specifiedValue)
{
	if (rateIdentifier>=nRate||compIdentifier>=nComp)
	{
		errorWarning(12);
		return 1;
	}
	rateMatrix[rateIdentifier*nRate+compIdentifier]=specifiedValue;
}

template <typename compType, typename rateType, typename powerType >
int timeDerivative<TYPENAME>::changeUpdMatrix(int rateIdentifier, int compIdentifier, powerType specifiedValue)
{
	if (rateIdentifier>=nRate||compIdentifier>=nComp)
	{
		errorWarning(12);
		return 1;
	}
	updateMatrix[rateIdentifier*nRate+compIdentifier]=specifiedValue;
}

template <typename compType, typename rateType, typename powerType >
int timeDerivative<TYPENAME>::loadParameter(vector<compType> & initComp, 
		vector<rateType> & inputRate,
		powerType * inputRateMatrix, powerType * inputUpdateMatrix)
{
	nComp=initComp.size();
	nRate=inputRate.size();
	if (nComp>MAXNCOMP||nRate>MAXNRATE) 
	{
		errorWarning(11);
		return 1;
	}

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
int timeDerivative<TYPENAME>::fileOpen(string & condition)
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
int timeDerivative<TYPENAME>::fileClose()
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
void timeDerivative<TYPENAME>::saveData()
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

#endif 	// __timeDerivative_H__ defined
