#ifndef __MARGINALDISTRIEXTRACTOR_H__
#define __MARGINALDISTRIEXTRACTOR_H__

#include "gillespie.h"

struct mDistriRecord
{
	bool state;
	long int * distri;
	double * timer;
	long int  index;
	long int numOfRecordPoint;
	int minValue, maxValue;
	double timeInterval;

	mDistriRecord()
	{
		state=0;
	}

	mDistriRecord(int minValueInput, int maxValueInput, long int nOfRecordPointInput, 
			double timeIntervalInput)
	{
		state=0;
		initiateRecorder(minValueInput, maxValueInput, nOfRecordPointInput, timeIntervalInput);
	}

	void initiateRecorder(int minValueInput, int maxValueInput, long int nOfRecordPointInput, 
			double timeIntervalInput)
	{
		minValue=minValueInput;
		maxValue=maxValueInput;
		numOfRecordPoint=nOfRecordPointInput;
		timeInterval=timeIntervalInput;
		timer = new double [numOfRecordPoint];
		for (int i=0;i<numOfRecordPoint;i++) timer[i]=i*timeInterval;
		distri= new long int [(maxValue-minValue+1)*numOfRecordPoint];
		reset();
		state=1;
	}
	void reset()
	{
		for (int i=0;i<(maxValue-minValue+1)*numOfRecordPoint;i++) distri[i]=0;
		index=0;
	}
	void pushPoint(int value)
	{
		distri[(maxValue-minValue+1)*index+value-minValue]++;
		index++;
	}
	void resetIndex()
	{
		index=0;	
	}
	~mDistriRecord()
	{
		if (state)
		{
			delete [] distri;
			delete [] timer;
		}
	}
};

//tempModule, to extract X and E+

class marginalDistriExtractor: public gillespie
{
public:
	mDistriRecord recordX;
	mDistriRecord recordE;

	void saveData();
	marginalDistriExtractor(vector<int> & initComp, vector<double> & inputRate,
			int * inputRateMatrix, int * inputUpdateMatrix, double runTime, 
			double saveInterval): 
		gillespie(initComp,inputRate,inputRateMatrix,inputUpdateMatrix, runTime,
				saveInterval)
	{
		long int numRecordPoint=(long int)(runTime/saveInterval+1);
		recordX.initiateRecorder(0,comp[0]+comp[3],numRecordPoint,saveInterval);
		recordE.initiateRecorder(0,comp[1]+comp[6],numRecordPoint,saveInterval);
	}
	void extract(int simTimes);
	virtual void reset()
	{
		gillespie::reset();
		recordX.index=0;
		recordE.index=0;
	}
	void saveDistri(string condition);
};

void marginalDistriExtractor::saveData()
	{
		recordX.pushPoint(comp[0]);
		recordE.pushPoint(comp[1]);
	}
void marginalDistriExtractor::extract(int simTimes)
{
	for (int i=0;i<simTimes;i++)
	{
		reset();
		simulate();
	}
}

void marginalDistriExtractor::saveDistri(string condition)
{
	string summaryFile;
	string distriFile;

	summaryFile+=condition;
	distriFile+=condition;
	summaryFile+="DistriXSummary";
	distriFile+="DistriX";

	fileOpen(summaryFile);
	resultFile<<recordX.minValue<<'\t'<<recordX.maxValue<<endl;
	resultFile<<recordX.numOfRecordPoint<<'\t'<<recordX.timeInterval;
	fileClose();
	fileOpen(distriFile);
	for (int i=0;i<recordX.numOfRecordPoint;i++)
	{
		for (int j=0;j<recordX.maxValue-recordX.minValue+1;j++)
		{
			resultFile<<recordX.distri[i*(recordX.maxValue-recordX.minValue+1)
				+j]<<'\t';
		}
		resultFile<<endl;
	}

	summaryFile.clear();
	distriFile.clear();
	summaryFile+=condition;
	distriFile+=condition;
	summaryFile+="DistriESummary";
	distriFile+="DistriE";

	fileOpen(summaryFile);
	resultFile<<recordE.minValue<<'\t'<<recordE.maxValue<<endl;
	resultFile<<recordE.numOfRecordPoint<<'\t'<<recordE.timeInterval;
	fileClose();
	fileOpen(distriFile);
	for (int i=0;i<recordE.numOfRecordPoint;i++)
	{
		for (int j=0;j<recordE.maxValue-recordE.minValue+1;j++)
		{
			resultFile<<recordE.distri[i*(recordE.maxValue-recordE.minValue+1)
				+j]<<'\t';
		}
		resultFile<<endl;
	}
}

#endif
