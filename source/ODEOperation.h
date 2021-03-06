#ifndef __ODESIMULATE_H__
#define __ODESIMULATE_H__

#include"basicDef.h"
#include"ODECommon.h"
#include "rungeKutta.h"

class ODESimulate : public ODENetwork , public RKmethod<ODESimulate>
{
//friend RKmethod<ODESimulate>;
public:
//runtime & saving related
	double stoppingTime;
	double saveTimeInterval;
	double lastSavedTime;

//simulation
	
	void reset();
	void simulate(string & identifier);

	ODESimulate(vector<double> & initComp, vector<double> & inputRate,  
			int * inputRateMatrix, int * inputUpdateMatrix, 
			double initTimeStep, double maxTimeStep_alias, double runTime, double saveInterval) : 
		ODENetwork(initComp, inputRate, inputRateMatrix, inputUpdateMatrix,
				initTimeStep), RKmethod<ODESimulate>(nComp, h0, &ODESimulate::ODETimeDeri, maxTimeStep_alias)
	{
		stoppingTime=runTime;
		saveTimeInterval=saveInterval;
		lastSavedTime=0;
	}
};

void ODESimulate::reset()
{
	BCNetwork<double,double,int>::reset();
	lastSavedTime=0;
	ht=h0;
}

void ODESimulate::simulate(string & identifier)
{
	fileOpen(identifier);
	while(t<stoppingTime)
	{
		t+=iterate(comp);	
		if (t>=lastSavedTime+saveTimeInterval)
		{
			saveData();
			lastSavedTime=t;
		}
	}
	fileClose();
}

#endif
