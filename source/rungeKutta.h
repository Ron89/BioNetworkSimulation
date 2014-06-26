#ifndef __rungeKutta_H__
#define __rungeKutta_H__

#include <fstream>
#include <string>
#include <algorithm>
#include <iostream>
#include "basicDef.h"
#include "ODECommon.h"

using namespace std;

template<typename modelClassName>
class RKmethod: public ODEIVPCommon<modelClassName>
{
public:
	double hMax; 			//setup a maximum value for timestep
	double initTimeStep_backup;
private:
	double * K1, * K2,* K3, * K4, * K5, * K6, * temp, * tempDeri;
	bool iteratorPrepared; 

	inline void RKAccumulate(double * fState, double * iState, double * deri,
		double step)		// Basic iteration from one state to the next state.
	{
		for (int i=0;i<ODEIVPCommon<modelClassName>::varNumber;i++)	fState[i]=iState[i]+step*deri[i];
	}

// constructor and destructor for the dynamically allocated memory for the algorithm
// Constructor should be called after the network is properly defined
	void constructor();
	void destructor();

public:
// Algorithm core, read current time and value for each variable, return
// the time and
// value of the next time step.
	double iterate(double * var);
	void reset()
	{
		ODEIVPCommon<modelClassName>::ht=initTimeStep_backup;
	}

// assign
	void assign(const RKmethod<modelClassName> dummy)
	{
		ODEIVPCommon<modelClassName>::assign(dummy);
		destructor();
		constructor();
		hMax=dummy.hMax;
		initTimeStep_backup=dummy.initTimeStep_backup;
	}

public:
//constructor and destructor
	RKmethod(double maxTimeStep=0)
	{
		hMax=maxTimeStep;
		iteratorPrepared=0;
	}
	RKmethod(int sysSize,
		double initTimeStep, 
		void (modelClassName::*targetODEs)(double *, double *), double maxTimeStep=0) :
		ODEIVPCommon<modelClassName>::ODEIVPCommon(sysSize, initTimeStep, targetODEs)
	{
		initTimeStep_backup=initTimeStep;
		hMax=maxTimeStep;
		iteratorPrepared=0;
		ODEIVPCommon<modelClassName>::Normalizer=blankNormalizer;
		constructor();
	}
	RKmethod(int sysSize,
		double initTimeStep, 
		void (modelClassName::*targetODEs)(double *, double *),
		void (*targetNormalizer)(double *), double maxTimeStep=0):
		ODEIVPCommon<modelClassName>::ODEIVPCommon(sysSize, initTimeStep, 
				targetODEs,targetNormalizer)
	{
		initTimeStep_backup=initTimeStep;
		hMax=maxTimeStep;
		iteratorPrepared=0;
		constructor();
	}

	~RKmethod()
	{
		destructor();
	}
};

template<typename modelClassName>
void RKmethod<modelClassName>::constructor()
{
	if(iteratorPrepared==0)
	{
		iteratorPrepared=1;
		K1=new double [ODEIVPCommon<modelClassName>::varNumber];
		K2=new double [ODEIVPCommon<modelClassName>::varNumber];
		K3=new double [ODEIVPCommon<modelClassName>::varNumber];
		K4=new double [ODEIVPCommon<modelClassName>::varNumber];
		K5=new double [ODEIVPCommon<modelClassName>::varNumber];
		K6=new double [ODEIVPCommon<modelClassName>::varNumber];
		temp=new double [ODEIVPCommon<modelClassName>::varNumber];
		tempDeri=new double [ODEIVPCommon<modelClassName>::varNumber];
		hMax=0;
	}
}

template<typename modelClassName>
void RKmethod<modelClassName>::destructor()
{
	if(iteratorPrepared==1)
	{
		iteratorPrepared=0;
		delete [] K1;
		delete [] K2;
		delete [] K3;
		delete [] K4;
		delete [] K5;
		delete [] K6;
		delete [] temp;
		delete [] tempDeri;
	}
}

template<typename modelClassName>
double RKmethod<modelClassName>::iterate(double * var)
{
	double maxDelta=0,tempDelta,deltaTime;

	for (int i=0;i<ODEIVPCommon<modelClassName>::varNumber;i++)	
		temp[i]=tempDeri[i]=K1[i]=K2[i]=K3[i]=K4[i]=K5[i]=K6[i]=0.;
	(static_cast<modelClassName*>(this)->*(ODEIVPCommon<modelClassName>::ODEs))(tempDeri,var);
	RKAccumulate(K1,K1,tempDeri,ODEIVPCommon<modelClassName>::ht);
	RKAccumulate(temp,var,K1,RKB21);
	(static_cast<modelClassName*>(this)->*(ODEIVPCommon<modelClassName>::ODEs))(tempDeri,temp);
	RKAccumulate(K2,K2,tempDeri,ODEIVPCommon<modelClassName>::ht);
	for (int i=0;i<ODEIVPCommon<modelClassName>::varNumber;i++)
		temp[i]=var[i]+RKB31*K1[i]+RKB32*
			K2[i];
	(static_cast<modelClassName*>(this)->*(ODEIVPCommon<modelClassName>::ODEs))(tempDeri,temp);
	RKAccumulate(K3,K3,tempDeri,ODEIVPCommon<modelClassName>::ht);
	for (int i=0;i<ODEIVPCommon<modelClassName>::varNumber;i++)
		temp[i]=var[i]+RKB41*K1[i]+RKB42*
			K2[i]+RKB43*K3[i];
	(static_cast<modelClassName*>(this)->*(ODEIVPCommon<modelClassName>::ODEs))(tempDeri,temp);
	RKAccumulate(K4,K4,tempDeri,ODEIVPCommon<modelClassName>::ht);
	for (int i=0;i<ODEIVPCommon<modelClassName>::varNumber;i++)
		temp[i]=var[i]+RKB51*K1[i]+RKB52*
			K2[i]+RKB53*K3[i]+RKB54*K4[i];
	(static_cast<modelClassName*>(this)->*(ODEIVPCommon<modelClassName>::ODEs))(tempDeri,temp);
	RKAccumulate(K5,K5,tempDeri,ODEIVPCommon<modelClassName>::ht);
	for (int i=0;i<ODEIVPCommon<modelClassName>::varNumber;i++)
		temp[i]=var[i]+RKB61*K1[i]+RKB62*
			K2[i]+RKB63*K3[i]+RKB64*K4[i]+RKB65*K5[i];
	(static_cast<modelClassName*>(this)->*(ODEIVPCommon<modelClassName>::ODEs))(tempDeri,temp);
	RKAccumulate(K6,K6,tempDeri,ODEIVPCommon<modelClassName>::ht);
	
	for (int i=0;i<ODEIVPCommon<modelClassName>::varNumber;i++)
		var[i]=var[i]+
			RKC01*K1[i]+RKC02*K2[i]+RKC03*K3[i]+RKC04*K4[i]+
			RKC05*K5[i]+RKC06*K6[i];
	deltaTime = ODEIVPCommon<modelClassName>::ht;

//deciding next stepsize
	for	(int i=0;i<ODEIVPCommon<modelClassName>::varNumber;i++)
	{
		tempDelta=abs(((RKC01-RKC11)*K1[i]+(RKC02-RKC12)*K2[i]+
				(RKC03-RKC13)*K3[i]+(RKC04-RKC14)*K4[i]+
				(RKC05-RKC15*K5[i])+(RKC06-RKC16)*K6[i])/
				(var[i]*ACCURACY));
		if (maxDelta==0) maxDelta=tempDelta;
		else if (tempDelta>maxDelta) maxDelta=tempDelta;
	}
	if (maxDelta>1.) ODEIVPCommon<modelClassName>::ht=ODEIVPCommon<modelClassName>::ht*pow(maxDelta,-0.2);
	else ODEIVPCommon<modelClassName>::ht=ODEIVPCommon<modelClassName>::ht*pow(maxDelta,-0.25);
	if (ODEIVPCommon<modelClassName>::ht>hMax && hMax!=0) 	ODEIVPCommon<modelClassName>::ht=hMax;

//Normalize result
	ODEIVPCommon<modelClassName>::Normalizer(var);

	return deltaTime;
}

#endif 	// __rungeKutta_H__ defined
