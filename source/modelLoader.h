#ifndef __MODELLOADER_H_INCLUDED__
#define __MODELLOADER_H_INCLUDED__
//===================

#include<string>
#include<vector>
#include<fstream>
#include"basicDef.h"

using namespace std;

template <typename compType, typename rateType, typename powerType >
class modelLoader
{
public:
	vector<rateType> rate;
	vector<compType> initComp;
	powerType *rateMatrix;
	powerType *updateMatrix;
	string modelName;
	int nComp;
	int nRate;
	

	void loadParameter(string & mName);
	
//used to generate/destroy matrices
	bool matricesState;
	inline void generateMatrices(int matrixSize)
	{
		if (matricesState==1)
		{
			destroyMaTrices();
		}
		rateMatrix=new powerType [matrixSize];
		updateMatrix=new powerType [matrixSize];
		matricesState=1;
	}
	inline void destroyMaTrices()
	{
		if (matricesState==0) 	return;
		delete [] rateMatrix;
		delete [] updateMatrix;
	}

//constructor&destructor
	modelLoader(string & mName)
	{
		matricesState=0;
		loadParameter(mName);
	}
	~modelLoader()
	{
		destroyMaTrices();	
	}
};

template <typename compType, typename rateType, typename powerType >
void modelLoader<TYPENAME>::loadParameter(string & mName)
{
//load the name of model. prepare to load the parameters
	modelName=mName;
	string folderName=modelName+"/";
	ifstream sourceFile;

//read ratefile
	string fileName=folderName+RATEFILE;
	sourceFile.open(fileName.c_str());
	sourceFile>>nRate;
	rate.resize(nRate,0);
	for (int i=0;i<nRate;i++) 	sourceFile>>rate[i];
	sourceFile.close();

//read initial condition
	fileName=folderName+INITCONDFILE;
	sourceFile.open(fileName.c_str());
	sourceFile>>nComp;
	initComp.resize(nComp,0);
	for (int i=0;i<nComp;i++) 	sourceFile>>initComp[i];
	sourceFile.close();

//create matrices
	generateMatrices(nComp*nRate);

//read rate matrix
	fileName=folderName+RATEMATRIXFILE;
	sourceFile.open(fileName.c_str());
	for (int i=0;i<nComp*nRate;i++) 	sourceFile>>rateMatrix[i];
	sourceFile.close();

//read update matrix
	fileName=folderName+UPDATEMATRIXFILE;
	sourceFile.open(fileName.c_str());
	for (int i=0;i<nComp*nRate;i++) 	sourceFile>>updateMatrix[i];
	sourceFile.close();
}

#endif 	//__MODELLOADER_H_INCLUDED__
