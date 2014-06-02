#include<sstream>
#include<cmath>
#include<string>

#include"coarseGrainedOperation.h"

using namespace std;


int main()
{
	double simTime=2000;
	double saveInterval=1;
	string modelName="p53Stochastic";
	string saveName="saveTesting";
	coarseGrainedStochastic sim(modelName,simTime,saveInterval);
	sim.simulate(saveName);
	return 0;
}
