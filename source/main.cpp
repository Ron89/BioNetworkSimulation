#include<sstream>
#include<cmath>
#include"modelLoader.h"
#include "gillespie.h"
#include"histoGen.h"
#include"basicDef.h"

#include<string>


using namespace std;


int main()
{
	string modelName("futileCycle");
	modelLoader<int, double, int> a(modelName);
	stringstream ss;

	double equivE;

//	histogramGenrator sim(a.initComp, a.rate, a.rateMatrix, a.updateMatrix);
//	sim.process();

	gillespie sim(a.initComp, a.rate, a.rateMatrix, a.updateMatrix, 1000, 1, 0.001);
//	a.initComp[6]=13;
	for (int i=24;i<60;i++)
//	for (double i=1.1;i<2.;i+=0.1)
	{
//		sim.rate[9]=1.3;
		sim.reset();
		sim.comp[6]=i-sim.comp[1];
		
//		sim.comp[1]=i;
		sim.rate[6]=sim.rate[7];
		equivE=(125.*i+51.-sqrt(pow(125.*i+51.,2.)-4.*75.*
					50.*i*(i+1.)))/(2.*75.);	
		sim.rate[9]=int(equivE/(1-equivE)*10);
		
		string condition("T1000p1N");
		ss<<i;
		condition+=ss.str();
		ss.str(std::string());
//		condition+=to_string(i);
		sim.fileOpen(condition);
		sim.simulate();
		sim.fileClose();
	}
}
