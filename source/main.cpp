#include<sstream>
#include<cmath>

#include"marginalDistriExtractor.h"
#include"modelLoader.h"

#include<string>

using namespace std;


int main()
{
	string modelName("futileCycle");
	modelLoader<int, double, int> a(modelName);
	stringstream ss;

	a.initComp[7]=26;

	marginalDistriExtractor sim(a.initComp, a.rate, a.rateMatrix, a.updateMatrix,
			10, 0.001);
	string resultName("N26T10");
	sim.extract(10000);
	sim.saveDistri(resultName);
	return 0;
}
