#ifndef __rungeKutta_H__
#define __rungeKutta_H__

#include <fstream>
#include <string>
#include <algorithm>
#include <iostream>
#include "basicDef.h"
#include "ODECommon.h"

//// coefficients for Multiple-Value method
//#define MV4R1 9./24.
//#define MV4R3 11./12.
//#define MV4R4 1./3.
//#define MV4R5 1./24.

#define MAXITER 10000

using namespace std;



template<typename modelClassType>
class RKmethod: public ODEIVPCommon<modelClassType>
{
private:
	double * K1, * K2,* K3, * K4, * K5, * K6, * temp, * tempDeri;

	inline void RKAccumulate(double * fState, double * iState, double * deri,
		double step)		// Basic iteration from one state to the next state.
	{
		for (int i=0;i<ODEIVPCommon<modelClassType>::varNumber;i++)	fState[i]=iState[i]+step*deri[i];
	}
// constructor and destructor for the dynamically allocated memory for the algorithm
// Constructor should be called after the network is properly defined
	void constructor();
	void destructor();

protected:
// Algorithm core, read current time and value for each variable, return
// the time and
// value of the next time step.
	void iterate(double * var, double time);

public:
//constructor and destructor
	RKmethod(int sysSize,
		double initTimeStep, 
		void (modelClassType::*targetODEs)(double *, double *)) :
		ODEIVPCommon<modelClassType>::ODEIVPCommon(sysSize, initTimeStep, targetODEs)
	{
		ODEIVPCommon<modelClassType>::Normalizer=blankNormalizer;
		constructor();
	}
	RKmethod(int sysSize,
		double initTimeStep, 
		void (modelClassType::*targetODEs)(double *, double *),
		void (*targetNormalizer)(double *)):
		ODEIVPCommon<modelClassType>::ODEIVPCommon(sysSize, initTimeStep, 
				targetODEs,targetNormalizer)

	{
		constructor();
	}

	~RKmethod()
	{
		destructor();
	}
};

template<typename modelClassType>
void RKmethod<modelClassType>::constructor()
{
	if(ODEIVPCommon<modelClassType>::iteratorPrepared==0)
	{
		ODEIVPCommon<modelClassType>::iteratorPrepared=1;
		K1=new double [ODEIVPCommon<modelClassType>::varNumber];
		K2=new double [ODEIVPCommon<modelClassType>::varNumber];
		K3=new double [ODEIVPCommon<modelClassType>::varNumber];
		K4=new double [ODEIVPCommon<modelClassType>::varNumber];
		K5=new double [ODEIVPCommon<modelClassType>::varNumber];
		K6=new double [ODEIVPCommon<modelClassType>::varNumber];
		temp=new double [ODEIVPCommon<modelClassType>::varNumber];
		tempDeri=new double [ODEIVPCommon<modelClassType>::varNumber];
	}
}

template<typename modelClassType>
void RKmethod<modelClassType>::destructor()
{
	if(ODEIVPCommon<modelClassType>::iteratorPrepared==1)
	{
		ODEIVPCommon<modelClassType>::iteratorPrepared=0;
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

template<typename modelClassType>
void RKmethod<modelClassType>::iterate(double * var, double time)
{
	double maxDelta=0,tempDelta;

	for (int i=0;i<ODEIVPCommon<modelClassType>::varNumber;i++)	
		temp[i]=tempDeri[i]=K1[i]=K2[i]=K3[i]=K4[i]=K5[i]=K6[i]=0.;
	(static_cast<modelClassType*>(this)->*(ODEIVPCommon<modelClassType>::ODEs))(tempDeri,var);
	RKAccumulate(K1,K1,tempDeri,ODEIVPCommon<modelClassType>::ht);
	RKAccumulate(temp,var,K1,RKB21);
	(static_cast<modelClassType*>(this)->*(ODEIVPCommon<modelClassType>::ODEs))(tempDeri,temp);
	RKAccumulate(K2,K2,tempDeri,ODEIVPCommon<modelClassType>::ht);
	for (int i=0;i<ODEIVPCommon<modelClassType>::varNumber;i++)
		temp[i]=var[i]+RKB31*K1[i]+RKB32*
			K2[i];
	(static_cast<modelClassType*>(this)->*(ODEIVPCommon<modelClassType>::ODEs))(tempDeri,temp);
	RKAccumulate(K3,K3,tempDeri,ODEIVPCommon<modelClassType>::ht);
	for (int i=0;i<ODEIVPCommon<modelClassType>::varNumber;i++)
		temp[i]=var[i]+RKB41*K1[i]+RKB42*
			K2[i]+RKB43*K3[i];
	(static_cast<modelClassType*>(this)->*(ODEIVPCommon<modelClassType>::ODEs))(tempDeri,temp);
	RKAccumulate(K4,K4,tempDeri,ODEIVPCommon<modelClassType>::ht);
	for (int i=0;i<ODEIVPCommon<modelClassType>::varNumber;i++)
		temp[i]=var[i]+RKB51*K1[i]+RKB52*
			K2[i]+RKB53*K3[i]+RKB54*K4[i];
	(static_cast<modelClassType*>(this)->*(ODEIVPCommon<modelClassType>::ODEs))(tempDeri,temp);
	RKAccumulate(K5,K5,tempDeri,ODEIVPCommon<modelClassType>::ht);
	for (int i=0;i<ODEIVPCommon<modelClassType>::varNumber;i++)
		temp[i]=var[i]+RKB61*K1[i]+RKB62*
			K2[i]+RKB63*K3[i]+RKB64*K4[i]+RKB65*K5[i];
	(static_cast<modelClassType*>(this)->*(ODEIVPCommon<modelClassType>::ODEs))(tempDeri,temp);
	RKAccumulate(K6,K6,tempDeri,ODEIVPCommon<modelClassType>::ht);
	
	for (int i=0;i<ODEIVPCommon<modelClassType>::varNumber;i++)
		var[i]=var[i]+
			RKC01*K1[i]+RKC02*K2[i]+RKC03*K3[i]+RKC04*K4[i]+
			RKC05*K5[i]+RKC06*K6[i];
	time+=ODEIVPCommon<modelClassType>::ht;

//deciding next stepsize
	for	(int i=0;i<ODEIVPCommon<modelClassType>::varNumber;i++)
	{
		tempDelta=abs(((RKC01-RKC11)*K1[i]+(RKC02-RKC12)*K2[i]+
				(RKC03-RKC13)*K3[i]+(RKC04-RKC14)*K4[i]+
				(RKC05-RKC15*K5[i])+(RKC06-RKC16)*K6[i])/
				(var[i]*ACCURACY));
		if (maxDelta==0) maxDelta=tempDelta;
		else if (tempDelta>maxDelta) maxDelta=tempDelta;
	}
	if (maxDelta>1.) ODEIVPCommon<modelClassType>::ht=ODEIVPCommon<modelClassType>::ht*pow(maxDelta,-0.2);
	else ODEIVPCommon<modelClassType>::ht=ODEIVPCommon<modelClassType>::ht*pow(maxDelta,-0.25);

//Normalize result
	ODEIVPCommon<modelClassType>::Normalizer(var);
}

#endif 	// __rungeKutta_H__ defined
