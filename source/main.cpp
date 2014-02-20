#include<sstream>
#include<cmath>

#include"marginalDistriExtractor.h"
#include"ODEOperation.h"
#include"modelLoader.h"

#include<string>

using namespace std;


int main()
{
	string modelName("futileCycle");
	modelLoader<double, double, int> a(modelName);
	stringstream ss;

	a.initComp[6]=26;

//	marginalDistriExtractor sim(a.initComp, a.rate, a.rateMatrix, a.updateMatrix,
//			10, 0.001);
// ODE Simulation	
	
	double hInit=1e-6;
	ODESimulate sim(a.initComp, a.rate, a.rateMatrix, a.updateMatrix,hInit,0.2,0.001);
	for (int Nvalue=0;Nvalue<=a.initComp[6]+a.initComp[1];Nvalue+=4)
	{
		sim.compBackup[6]=Nvalue;
		sim.compBackup[1]=a.initComp[6]+a.initComp[1]-sim.compBackup[6];
		for (int Xvalue=0;Xvalue<=a.initComp[0]+a.initComp[3];Xvalue+=200)
		{
			sim.compBackup[0]=Xvalue;
			sim.compBackup[3]=a.initComp[0]+a.initComp[3]-Xvalue;
			sim.reset();
			string resultName("ODENE46T10N");
			ss<<Nvalue;
			resultName+=ss.str();
			ss.str("");
			resultName+="X";
			ss<<Xvalue;
			resultName+=ss.str();
			ss.str("");
			sim.simulate(resultName);
		}
	}
	return 0;
}
