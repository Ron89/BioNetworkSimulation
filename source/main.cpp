#include<sstream>
#include<cmath>
#include<string>

#include"coarseGrainedOperation.h"
#include"coarseGrainedCommon.h"

using namespace std;


int main()
{
	double simTime=10.;
	double saveInterval=0.01;
	string modelName="fCycle";
	string saveName="fCycleT10N10Data";
//	coarseGrainedStochastic sim(modelName,simTime,saveInterval);
	coarseGrainedDeterministic sim(modelName,saveInterval/1000,saveInterval/10.,simTime,saveInterval);
//	sim.comp[8]*=.3;
//	sim.comp[9]*=.3;
	sim.simulate(saveName);
	return 0;
}
