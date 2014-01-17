#include<sstream>
#include<cmath>

#include"basicDef.h"
#include"modelLoader.h"
#include "ODEOperation.h"

#include<string>


using namespace std;


int main()
{
	string modelName("futileCycle");
	modelLoader<double, double, int> a(modelName);
	stringstream ss;

	ODESimulate sim(a.initComp, a.rate, a.rateMatrix, a.updateMatrix,
			0.000001, 1000, 0.01);
	string resultName("ODET1000");
	sim.simulate(resultName);
	return 0;
}
