#ifndef __COARSEGRAINEDOPERATION_H_INCLUDED__
#define __COARSEGRAINEDOPERATION_H_INCLUDED__

#include"coarseGrainedCommon.h"
#include"gillespie.h"
#include"rungeKutta.h"
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
	void reallocate(int nComp_alias, long trajectorySize);
	void erase();
	
// constructors
	trajectory();
	trajectory(int nComp_alias, long trajectorySize_alias);
// copy constructor;
	trajectory( trajectory<compType> & dummy);
	trajectory & operator=( trajectory<compType> & dummy);

// destructor
	~trajectory();

// operations
	void append(double time_alias, compType * comp_alias);

// savings
	void save( string & outputFileName, bool append=1);
};

// --------------------------

class coarseGrainedStochastic: public coarseGrainedModel<int,double> , 
							   public gillespie<coarseGrainedStochastic>
{
public:

// update reset function
	void reset();

// update assign function
	void assign(const coarseGrainedStochastic & dummy);

// constructors
	coarseGrainedStochastic():coarseGrainedModel<int,double>(),gillespie<coarseGrainedStochastic>
							  (nReact,&coarseGrainedModel<int,double>::rateDetermine,
													&coarseGrainedModel<int,double>::reactantUpdate)
	{}
// load from file 
	coarseGrainedStochastic(string & modelName, double stoppingTime_alias,
			double saveTimeInterval_alias):
		coarseGrainedModel<int,double>(modelName), gillespie<coarseGrainedStochastic>
												   (nReact,
													&coarseGrainedModel<int,double>::rateDetermine,
													&coarseGrainedModel<int,double>::reactantUpdate)
	{
		stoppingTime=stoppingTime_alias;
		saveTimeInterval=saveTimeInterval_alias;
		reset();
	}
	
// load from same type of model
	coarseGrainedStochastic(coarseGrainedStochastic & dummy):
		coarseGrainedModel<int,double>(dummy),
		gillespie<coarseGrainedStochastic>(nReact,
				&coarseGrainedModel<int,double>::rateDetermine,
				&coarseGrainedModel<int,double>::reactantUpdate)
	{}
	~coarseGrainedStochastic(){}

// assign operator
	coarseGrainedStochastic & operator=(coarseGrainedStochastic & dummy)
	{
		if (this!=&dummy) 	assign(dummy);
		return *this;
	}

// function member
	void simulate(string fileName);
	void simulate(bool (*controller)(double, int *)); 	// This simulation module will hand over
															// saving interval to dataCollector.
};

// --------------------------

class coarseGrainedDeterministic: 	public coarseGrainedModel<double,double>,
									public RKmethod<coarseGrainedDeterministic>
{
private:
	bool initState;
	double * reactionRate;

// allocate and free memory	
	void generateModel();
	void eraseModel();

public:
// update reset function
	void reset();

// update assign function 
	void assign(const coarseGrainedDeterministic & dummy);

	void modelODE(double * timeDeri_alias, double * comp_alias);

// constructors
// load from file 
	coarseGrainedDeterministic(string & modelName, double initTimeStep_alias, 
			double maxTimeStep_alias,
			double stoppingTime_alias,
			double saveTimeInterval_alias):
		coarseGrainedModel<double,double>(modelName), RKmethod<coarseGrainedDeterministic>
												   (nComp, initTimeStep_alias, 
													&coarseGrainedDeterministic::modelODE, 
													maxTimeStep_alias)
	{
		initState=0;
		generateModel();
		stoppingTime=stoppingTime_alias;
		saveTimeInterval=saveTimeInterval_alias;
		reset();
	}
	
// load from same type of model 
	coarseGrainedDeterministic(coarseGrainedDeterministic & dummy):
		coarseGrainedModel<double,double>(dummy),
		RKmethod<coarseGrainedDeterministic>(nComp,
				dummy.ht, &coarseGrainedDeterministic::modelODE, dummy.hMax)
	{
		initState=0;
		assign(dummy);
	}
	~coarseGrainedDeterministic()
	{
		eraseModel();
	}

// function member
	void simulate(string fileName);
};

// -------------------------

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

template<typename compType>
void trajectory<compType>::reallocate(int nComp_alias, long trajectorySize)
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
	trajectoryPointer=0;
}

template<typename compType>
void trajectory<compType>::erase()
{
	if (trajectoryDefined)
	{
		delete [] time;
		delete [] comp;
	}
	trajectoryDefined=0;
}

template<typename compType>
trajectory<compType>::trajectory()
{
	trajectoryDefined=0;
}

template<typename compType>
trajectory<compType>::trajectory(int nComp_alias, long trajectorySize_alias)
{
	trajectoryDefined=0;
	trajectoryPointer=0;
	nComp=nComp_alias;
	trajectorySize=trajectorySize_alias;
	reallocate(nComp, trajectorySize);
}

template<typename compType>
trajectory<compType>::trajectory( trajectory<compType> & dummy)
{
	trajectoryDefined=0;
	assign(dummy);
}

template<typename compType>
trajectory<compType> & trajectory<compType>::operator=(trajectory<compType> & dummy) 
{
	if (this!=&dummy) 	assign(dummy);
	return *this;
}

template<typename compType>
trajectory<compType>::~trajectory()
{
	erase();
}

template<typename compType>
void trajectory<compType>::append(double time_alias, compType * comp_alias)
{
	time[trajectoryPointer]=time_alias;
	for (int i=0; i<nComp; i++) 	comp[trajectoryPointer*nComp+i]=comp_alias[i];
	trajectoryPointer++;
}

template<typename compType>
void trajectory<compType>::save( string & outputFileName, bool append)
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

// --------------------------

void coarseGrainedStochastic::reset()
{
	coarseGrainedModel<int,double>::reset();
}

void coarseGrainedStochastic::assign(const coarseGrainedStochastic & dummy)
{
	coarseGrainedModel<int,double>::assign(dummy);
	gillespie::assign(dummy);
}

void coarseGrainedStochastic::simulate(string fileName)
{
	trajectory<int> trajStorage(nComp, int((stoppingTime-time)/saveTimeInterval)+1);
	while (time<stoppingTime)
	{
		time+=iterate(comp);
		if ((time-lastSavedTime)>=saveTimeInterval)
		{
			trajStorage.append(time,comp);
			lastSavedTime+=saveTimeInterval;
		}
	}
	trajStorage.save(fileName,0);
}

void coarseGrainedStochastic::simulate(bool (*controller)(double, int *))
{
	do time+=iterate(comp);
	while (controller(time, comp));
}

// --------------------------

void coarseGrainedDeterministic::simulate(string fileName)
{
	trajectory<double> trajStorage(nComp, int((stoppingTime-time)/saveTimeInterval)+1);
	while (time<stoppingTime)
	{
		time+=iterate(comp);
		if ((time-lastSavedTime)>=saveTimeInterval)
		{
			trajStorage.append(time,comp);
			lastSavedTime+=saveTimeInterval;
		}
	}
	trajStorage.save(fileName,0);
}

void coarseGrainedDeterministic::generateModel()
{
	if (initState) 	eraseModel();
	reactionRate=new double [nReact];
	initState=1;
}

void coarseGrainedDeterministic::eraseModel()
{
	if (initState) 	delete [] reactionRate;
	initState=0;
}

void coarseGrainedDeterministic::reset()
{
	coarseGrainedModel<double,double>::reset();
	RKmethod<coarseGrainedDeterministic>::reset();
}

void coarseGrainedDeterministic::assign(const coarseGrainedDeterministic & dummy)
{
	coarseGrainedModel<double,double>::assign(dummy);
//		for (int i=0;i<nReact;i++) 	reactionRate[i]=dummy.reactionRate[i];
	RKmethod<coarseGrainedDeterministic>::assign(dummy);
	generateModel();
}

void coarseGrainedDeterministic::modelODE(double * timeDeri_alias, double * comp_alias)
{
	rateDetermine(reactionRate, comp_alias);
	for (int i=0;i<nComp;i++) 	timeDeri_alias[i]=0;
	reactantUpdate(timeDeri_alias, reactionRate);
}

#endif
