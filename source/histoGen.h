#ifndef __HISTOGEN_H__
#define __HISTOGEN_H__

#include<iostream>
#include<cstdlib>
#include<vector>
#include<algorithm>
#include "gillespie.h"
#include "basicDef.h"

using namespace std;

class histogramGenrator: public gillespie
{
public:

	void saveData();
	void process();
	vector<long long int> histoGram;

//constructor
	histogramGenrator(vector<int> & initComp, vector<double> & inputRate,
			int * inputRateMatrix, int * inputUpdateMatrix)
	{
		loadParameter(initComp, inputRate, inputRateMatrix, inputUpdateMatrix);

		stoppingTime=10;
		saveMethod=0;
		savePointInterval=1;

		reset();
		mkdir(RESULTFOLDER, 0755);
		reseedRandom(1);

		histoGram.resize(initComp[0]+initComp[3]+1,0);
		noSave=1; 	//allow data saving after running period
	}
};

//histogramGenrator::histogramGenrator(vector<int> & initComp, vector<double> & inputRate,
//			int * inputRateMatrix, int * inputUpdateMatrix)

void histogramGenrator::saveData()
{
	histoGram[comp[0]]++;
}

void histogramGenrator::process()
{
	string histogramFile("HistogramE7ExtendedSample");
	fileOpen(histogramFile);
	resultFile<<"#Condition, vary initial N form 4 to 20 stepping 1"<<endl;
	
	for (int N=4;N<64;N++)
	{
		reset();
		comp[6]=N;
		noSave=1;

		while(*max_element(histoGram.begin(),histoGram.end())<1E7)
		{
			simulate();
			noSave=0;
		}
		
		cout<<histoGram.size()<<endl;
		
		for (unsigned int i=0;i<histoGram.size();i++) 	
		{
			resultFile<<histoGram[i]<<"\t";	
			histoGram[i]=0;
		}
		resultFile<<endl;
	}
	fileClose();
}

#endif 	//__HISTOGEN_H__ defined
