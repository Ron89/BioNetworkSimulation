#include<sstream>
#include<cmath>
#include<string>

#include"coarseGrainedOperation.h"

using namespace std;


int main()
{
	double simTime=200;
	double saveInterval=0.1;
	string modelName="p53Stochastic";
	string saveName="saveTesting";
	coarseGrainedStochastic sim(modelName,simTime,saveInterval);
	sim.comp[8]*=.3;
	sim.comp[9]*=.3;
	sim.simulate(saveName);
	return 0;
}
