#ifndef __COARSEGRAINEDOPERATION_H_INCLUDED__
#define __COARSEGRAINEDOPERATION_H_INCLUDED__

#include"coarseGrainedCommon.h"
#include"gillespie.h"
#include<fstream>
#include<string>
#include<iostream>

using namespace std;

/* (not implimented yet)
 * here we support 2 types of data storage,
 * 	memSaving = 0: 	data will be stored in memory. This method is designed to be multithread
 * 					friendly;
 * 	memSaving = 1; 	data will be stored in file. As fstream is not designed for multithread,
 * 					use such mothod under multithread scenario might cause unexpected behavior.
 */

template<typename compType>
struct trajectory
{
private:
	bool trajectoryDefined;
	long trajectorySize;

	void assign( trajectory & dummy);
public:
	int nComp;
	long trajectoryPointer;
	double * time;
	compType * comp;
	void reallocate(int nComp_alias, unsigned long trajectorySize)
	{
		if (trajectoryDefined) 	erase();
		time=new double [trajectorySize];
		comp=new compType [trajectorySize*nComp_alias];
		for (int i=0;i<trajectorySize;i++)
		{
			time[i]=0;
			for (int j=0;j<nComp_alias;j++) 	comp[i*nComp_alias+j]=0;
		}
		trajectoryDefined=1;
	}
	void erase()
	{
		if (trajectoryDefined)
		{
			delete [] time;
			delete [] comp;
		}
		trajectoryDefined=0;
	}
	
// constructors
	trajectory()
	{
		trajectoryDefined=0;
	}
	trajectory(int nComp_alias, unsigned long trajectorySize_alias)
	{
		trajectoryDefined=0;
		trajectoryPointer=0;
		nComp=nComp_alias;
		trajectorySize=trajectorySize_alias;
		reallocate(nComp, trajectorySize);
	}
// copy constructor;
	trajectory( trajectory<compType> & dummy)
	{
		trajectoryDefined=0;
		assign(dummy);
	}
	trajectory & operator=( trajectory<compType> & dummy)
	{
		if (this!=&dummy) 	assign(dummy);
		return *this;
	}

// operations
	void append(double time_alias, compType * comp_alias)
	{
		time[trajectoryPointer]=time_alias;
		for (int i=0; i<nComp; i++) 	comp[trajectoryPointer*nComp+i]=comp_alias[i];
		trajectoryPointer++;
	}

// savings
	void save( string & outputFileName, bool append=1)
	{
		ofstream filePointer;
		if (append==1)
			filePointer.open(outputFileName.c_str(), ios::out | ios::app);
		else
			filePointer.open(outputFileName.c_str(), ios::out | ios::trunc);
		for (int i=0;i<trajectoryPointer;i++)
		{
			filePointer<<time[i]<<'\t';
			for (int j=0;j<nComp;j++) 	filePointer<<comp[i*nComp+j]<<'\t';
			filePointer<<endl;	
		}
		filePointer.close();
	}
};

template<typename compType>
void trajectory<compType>::assign( trajectory<compType> & dummy)
{
	trajectorySize=dummy.trajectorySize;
	trajectoryPointer=dummy.trajectoryPointer;
	if (dummy.trajectoryDefined) 	
	{
		reallocate(nComp, trajectorySize);
		for (int i=0;i<trajectorySize;i++)
		{
			time[i]=dummy.time[i];
			for (int j=0;j<nComp;j++) 	comp[i*nComp+j]=dummy.comp[i*nComp+j];
		}
	}
}

class coarseGrainedStochastic: public coarseGrainedModel<int,double> , 
							   public gillespie<coarseGrainedStochastic>
{
public:

// update reset function
	void reset()
	{
		coarseGrainedModel<int,double>::reset();
	}

// update assign function
	void assign(const coarseGrainedStochastic & dummy)
	{
		coarseGrainedModel<int,double>::assign(dummy);
	}

// constructors
// load from file 
	coarseGrainedStochastic(string & modelName, double stoppingTime_alias,
			double saveTimeInterval_alias):
		coarseGrainedModel<int,double>(modelName), gillespie<coarseGrainedStochastic>
												   (nReact,
													&coarseGrainedModel::rateDetermine,
													&coarseGrainedModel::reactantUpdate)
	{
		stoppingTime=stoppingTime_alias;
		saveTimeInterval=saveTimeInterval_alias;
		reset();
	}
	
// load from same type of model
	coarseGrainedStochastic(coarseGrainedStochastic & dummy):
		coarseGrainedModel<int,double>(dummy),
		gillespie<coarseGrainedStochastic>(nReact,
				&coarseGrainedModel::rateDetermine,
				&coarseGrainedModel::reactantUpdate)
	{
	//	assign(dummy);
	}
	~coarseGrainedStochastic()
	{
	}

// function member
	void simulate(string fileName);
};

void coarseGrainedStochastic::simulate(string fileName)
{
	double currentProgress=0;
	double timeStart=time;
	double stage=0;
	trajectory<int> trajStorage(nComp, int((stoppingTime-time)/saveTimeInterval)+1);
	while (time<stoppingTime)
	{
		time+=iterate();
		if ((time-lastSavedTime)>=saveTimeInterval)
		{
			trajStorage.append(time,comp);
			lastSavedTime=time;
		}
		currentProgress=(time-timeStart)/(stoppingTime-timeStart);
		if (currentProgress>stage)
		{
			cout<<"progress: "<<currentProgress<<endl;
			stage+=0.01;
		}
	}
	trajStorage.save(fileName);
}


#endif
