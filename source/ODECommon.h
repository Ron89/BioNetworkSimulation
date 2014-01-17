#ifndef __ODECOMMON_H__
#define __ODECOMMON_H__

#include "BCNetwork.h"

class ODENetwork : public BCNetwork<double, double, int>
{
protected:
	double h0; 	//initial stepsize
public:
	ODENetwork(vector<double> & initComp, vector<double> & inputRate,  
			int * inputRateMatrix, int * inputUpdateMatrix, double 
			initTimeStep) : 
		BCNetwork(initComp, inputRate, inputRateMatrix, inputUpdateMatrix)
	{
		h0=initTimeStep;
	}
	void ODETimeDeri(double * timeDeri, double * component);
};

template<typename modelClassType>
class ODEIVPCommon
{
public:
	double ht;
	int varNumber;
	bool iteratorPrepared; 

	void (modelClassType::*ODEs)(double *, double *);
	int (*Normalizer)(double *);

	ODEIVPCommon()
	{
		iteratorPrepared=0;
	}
	ODEIVPCommon(int sysSize,
		double initTimeStep,
		void (modelClassType::*targetODEs)(double *, double *))
	{
		varNumber=sysSize;
		ht=initTimeStep;
		ODEs=targetODEs;
	}
	ODEIVPCommon(int sysSize,
		double initTimeStep, 
		void (modelClassType::*targetODEs)(double *, double *),
		int (*targetNormalizer)(double *))
	{
		varNumber=sysSize;
		ht=initTimeStep;
		ODEs=targetODEs;
		Normalizer=targetNormalizer;
	}
};

int blankNormalizer(double * a)
{
		return 0;
}

void ODENetwork::ODETimeDeri(double * timeDeri, double * component)
{
	double tempDepend=1;
	for (int i=0;i<nComp;i++)
	{
		timeDeri[i]=0;
		for (int j=0;j<nRate;j++)
		{
			tempDepend=updateMatrix[j*nComp+i]*rate[j];
			for (int k=0;k<nComp;k++) 	tempDepend*=pow(component[k],
					rateMatrix[j*nComp+k]);
			timeDeri[i]+=tempDepend;
		}
	}
}


#endif 	
