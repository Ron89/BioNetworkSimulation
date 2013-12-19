#ifndef __GILLESPIE_H_INCLUDED__
#define __GILLESPIE_H_INCLUDED__
//====================

#include<iostream>
#include<cstdlib>
#include<ctime>
#include<fstream>
#include<string>
#include<algorithm>
#include<cmath>
#include<vector>
// mersenne twister random generator

#include<random>

//used to generate directory
#include<sys/types.h>
#include<sys/stat.h>
#include<unistd.h>

#include"timeDerivative.h"
#include"basicDef.h"
//#include "dSFMT.h"


using namespace std;

class gillespie: public timeDerivative<int,double,int>
{
public:
//constructor
	gillespie(){}

	gillespie(vector<int> & initComp, vector<double> & inputRate,  
			int * inputRateMatrix, int * inputUpdateMatrix, 
			double runTime ,bool sMethod, double saveInterval);

//random generator
//	dsfmt_t dsfmt;
	inline void reseedRandom(int seed)
	{
		srand48(seed);
//		randGen.seed(seed);
//		dsfmt_init_gen_rand(&dsfmt, seed);
	}

	inline double popRandom()
	{
		return drand48();
//		return dsfmt_genrand_close_open(&dsfmt);
//		return randGen.operator();
	}

//algorithm functional parts
	void simulate(); 	//simulation module;

//data saving&stopping related
	bool saveMethod; 				//0-> save every savePointInterval steps;
									//1-> save every saveTimeInterval of time;
	bool noSave; 					//if noSave signal is 1. No data saving is allowed;	
	long long int savePointInterval;
	long long int nOfNewStep;
	double saveTimeInterval;
	double lastSavedTime;
	double stoppingTime;

//reset
	inline void reset()
	{
		for (int i=0;i<nComp;i++) 	comp[i]=compBackup[i];
		t=lastSavedTime=0;
		nOfNewStep=0;
	}
};

gillespie::gillespie(vector<int> & initComp, vector<double> & inputRate,
		int * inputRateMatrix, int * inputUpdateMatrix, 
			double runTime ,bool sMethod=0, double saveInterval=1)
{
	loadParameter(initComp, inputRate, inputRateMatrix, inputUpdateMatrix);

	stoppingTime=runTime;
	saveMethod=sMethod;
	if(saveMethod==0) 	savePointInterval=int(saveInterval);
	else 	saveTimeInterval=saveInterval;
	reset();
	mkdir(RESULTFOLDER, 0755);
	reseedRandom(1);

	noSave=0; 	//allow data saving after running period
}

void gillespie::simulate()
{
	double haltSig=(lastSavedTime=t)+stoppingTime;
	double r1, r2;
	double lambda;
	double * lambdaSig=new double[nRate+1];
	double dt;
	int ll,rl,templ;
	
	do
	{
		r1=popRandom();
		lambda=0;
		lambdaSig[0]=0;
		for (int i=1;i<=nRate;i++)
		{
			lambdaSig[i]=rate[i-1];
			for (int j=0;j<nComp;j++)
				lambdaSig[i]*=pow(double(comp[j]),double(rateMatrix[nComp*(i-1)+j]));
			lambdaSig[i]+=lambdaSig[i-1];
		}
		lambda=lambdaSig[nRate];
		for (int i=1;i<nRate+1;i++) 	lambdaSig[i]/=lambda;
		dt=-log(r1)/lambda;
		
		lambda=0;
		while(lambda==0)
		{
			lambda=1;
			r2=popRandom();
			for (int i=0;i<=nRate;i++) 	lambda*=(r2-lambdaSig[i]);
		}

		ll=0;
		rl=nRate;

		while(rl-ll>1)
		{
			templ=(ll+rl)/2;
			if(r2<lambdaSig[templ]) 	rl=templ;
			else 	ll=templ;
		}
		for (int i=0;i<nComp;i++) 	comp[i]+=updateMatrix[ll*nComp+i];
		t+=dt;
		nOfNewStep+=1;
		
// following is for deciding whether to save the current step.
		if (noSave==0&&saveMethod)		
		{
			if(saveTimeInterval<=t-lastSavedTime)
			{	saveData();
 				lastSavedTime=t;
			}
		}
		else if (noSave==0)
		{
			if(savePointInterval<=nOfNewStep) 
			{
				saveData();
				nOfNewStep=0;
			}
		}
	}	
	while (t<=haltSig);
	delete [] lambdaSig;
}



#endif 	//__GILLESPIE_H_INCLUDED__
