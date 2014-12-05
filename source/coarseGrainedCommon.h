#ifndef __COARSEGRAINEDCOMMON_H_INCLUDED__
#define __COARSEGRAINEDCOMMON_H_INCLUDED__

#include <string>
#include <fstream>
#include <stack>
#include <cmath>
#include "basicDef.h"

using namespace std;

// Data structure defined for storage of information for each reaction.
struct reaction
{
private:
	void assign(const reaction & dummy);

public:
	int code;
	double rate[MAXRATECOEF];
	int dependency[MAXDEPENDENCY];
	int updateNumber;
	int updateSet[MAXUPDSET];
	double updateOrder[MAXUPDSET];

	bool loadReaction(ifstream & reactionFile);

// constructor
	reaction();

// copy constructor
	reaction(const reaction & dummy);

// operators
	reaction & operator=(const reaction & dummy);
	void scale(double eta_c_alias, double eta_t_alias);
	void erase();
};

template<typename compType, typename updRateType>
class coarseGrainedModel
{
/* This module and other objects based on this module is designed for coarse grained modelling.
 * Other than elementary-reaction-based modeling, this type of model falls into undetermined
 * number of stereotypes: elementary reaction types, MM kinetics, hill-equation, etc. based
 * on ways of simplification. Thus for these modules, definition of its members are also more
 * open to change. 
 */
private:
// Flags
	bool initState;

// backup comp, for reset use
	compType * compBackup;

// allocate and free memory, only used when constructing/ destructing the model.
	void generateModel();
	void eraseModel();

public:
// reactants and reactions;
	int nComp;
	compType * comp;
	int nReact;
	reaction * react;

// timing related definition;
	double time; 				// current time
	double lastSavedTime; 		// last saved time point
	double stoppingTime; 		// time when simulation ends
	double saveTimeInterval; 	// time interval between two saving points

// Model loading
	int modelLoading(string & modelName);

// reaction rate determination
	int rateDetermine(double * rate, compType * comp_alias);
// in specific modules, other reaction types might be defined with additional routines.
// This function is used to add those reaction routines in.
	virtual double extendedCases(reaction currentReaction, compType * comp_alias){
		return 0;
	}

// in case the reactant level used to calculate rate is the current reactant level.
	inline int rateDetermine(double * rate);
	
// update comp based on a given rate
// 		* in deterministic case, use updRateType=double;
// 		* in gillespie algorithm, use updRateType=int;
	int reactantUpdate(compType * comp_alias, updRateType * rate);
	int reactantUpdate(updRateType * rate);
	
// constructor & destructor
	coarseGrainedModel()
	{
		initState=0;
	}
	coarseGrainedModel(string & modelName);
	coarseGrainedModel(const coarseGrainedModel<compType, updRateType> & dummy);
	~coarseGrainedModel();

//reset
	void reset();

//operator assign, duplicate all model information from model dummy;
	void assign(const coarseGrainedModel<compType, updRateType> & dummy);
	coarseGrainedModel & operator=(const coarseGrainedModel<compType, updRateType> & dummy);
};

void reaction::assign(const reaction & dummy)
{
	code=dummy.code;
	for (int i=0;i<MAXRATECOEF;i++) 	rate[i]=dummy.rate[i];
	for (int i=0;i<MAXDEPENDENCY;i++) 	dependency[i]=dummy.dependency[i];
	updateNumber=dummy.updateNumber;
	for (int i=0;i<MAXUPDSET;i++)
	{
		updateSet[i]=dummy.updateSet[i];
		updateOrder[i]=dummy.updateOrder[i];
	}
}

bool reaction::loadReaction(ifstream & reactionFile)
{
	reactionFile>>code;
	for (int i=0;i<MAXRATECOEF;i++) 	reactionFile>>rate[i];
	for (int i=0;i<MAXDEPENDENCY;i++) 	reactionFile>>dependency[i];
	reactionFile>>updateNumber;
	for (int i=0;i<MAXUPDSET;i++)
	{
		reactionFile>>updateSet[i];
		reactionFile>>updateOrder[i];
	}
	return reactionFile.eof();
}

void reaction::erase()
{
	code=0;
	for (int i=0;i<MAXRATECOEF;i++) 	rate[i]=0;
	for (int i=0;i<MAXDEPENDENCY;i++) 	dependency[i]=0;
	updateNumber=0;
	for (int i=0;i<MAXUPDSET;i++)
	{
		updateSet[i]=0;
		updateOrder[i]=0;
	}
}

reaction::reaction()
	{
		erase();
	}

reaction::reaction(const reaction & dummy)
	{
		assign(dummy);
	}

reaction & reaction::operator=(const reaction & dummy)
	{
		if (this!=&dummy) 	assign(dummy);
		return *this;
	}

// scaling the reaction network. As shown in the document
// for rate[0]
//	Type 	[c] 	[t]
// 	0 		1 		-1
// 	1 		0 		-1
// 	2 		-1 		-1
// 	3 		0 		-1
// 	4 		0 		-1
// for rate[1]
//	Type 	[c] 	[t]
// 	0 		0 		0
// 	1 		0 		0
// 	2 		0 		0
// 	3 		1 		0
// 	4 		1 		0
void reaction::scale(double eta_c_alias, double eta_t_alias)
	{
		rate[0]*=pow(eta_c_alias,(code==0?1:(code==2?-1:0)))*pow(eta_t_alias,-1);
		rate[1]*=pow(eta_c_alias,((code==3||code==4)?1:0));
	}


template<typename compType, typename updRateType>
void coarseGrainedModel<compType, updRateType>::assign
(const coarseGrainedModel<compType,updRateType> & dummy)
{
// duplicate the structure of current model
	nComp=dummy.nComp;
	nReact=dummy.nReact;
	generateModel();
	for (int i=0;i<nComp;i++) 	comp[i]=dummy.comp[i];
	for (int i=0;i<nReact;i++) 	react[i]=dummy.react[i];

// duplicate time information
	time=dummy.time;
	lastSavedTime=dummy.lastSavedTime;
	stoppingTime=dummy.stoppingTime;
	saveTimeInterval=dummy.saveTimeInterval;
}

template<typename compType, typename updRateType>
void coarseGrainedModel<compType,updRateType>::generateModel()
{
	if (initState) 	eraseModel();
	comp=new compType[nComp];
	compBackup=new compType[nComp]; 	//for the purpose of restoring to init condition.
	react=new reaction[nReact]; 	
	initState=1;
}

template<typename compType, typename updRateType>
void coarseGrainedModel<compType,updRateType>::eraseModel()
{
	if (initState) 
	{
		delete [] comp;
		delete [] compBackup;
		delete [] react;
	}
	initState=0;
}

template<typename compType, typename updRateType>
int coarseGrainedModel<compType,updRateType>::rateDetermine(double * rate, compType * comp_alias)
{
	for (int i=0;i<nReact;i++)
	{
		switch(react[i].code)
		{
			case 0 : 
				rate[i]=react[i].rate[0];
				break;
			case 1 :
				rate[i]=react[i].rate[0]*comp_alias[react[i].dependency[0]];
				break;
			case 2 :
				rate[i]=react[i].rate[0]*comp_alias[react[i].dependency[0]]*
					comp_alias[react[i].dependency[1]];
				break;
			case 3 :
				rate[i]=react[i].rate[0]*comp_alias[react[i].dependency[0]]*
					(pow(comp_alias[react[i].dependency[1]],react[i].rate[2]))/
					(pow(react[i].rate[1],react[i].rate[2])+
					 pow(comp_alias[react[i].dependency[1]],react[i].rate[2]));
				break;
			case 4 :
				rate[i]=react[i].rate[0]*comp_alias[react[i].dependency[0]]*
					comp_alias[react[i].dependency[1]]/
					(comp_alias[react[i].dependency[0]]+react[i].rate[1]);
				break;
			default :
				rate[i]=extendedCases(react[i],comp_alias);
		}
	}
	return 0;
}

template<typename compType, typename updRateType>
int coarseGrainedModel<compType,updRateType>::rateDetermine(double * rate)
	{
		return rateDetermine(rate, comp);	
	}

template<typename compType, typename updRateType>
int coarseGrainedModel<compType,updRateType>::reactantUpdate
(compType * comp_alias, updRateType * rate)
{
	for (int i=0;i<nReact;i++)
	{
		for (int j=0;j<react[i].updateNumber;j++)
			comp_alias[react[i].updateSet[j]]+=compType(rate[i]*react[i].updateOrder[j]);
	}
	return 0;
}

template<typename compType, typename updRateType>
int coarseGrainedModel<compType,updRateType>::reactantUpdate
(updRateType * rate)
{
	return reactantUpdate(comp, rate);
}

template<typename compType, typename updRateType>
int coarseGrainedModel<compType,updRateType>::modelLoading(string & modelName)
{
// load the name of model, prepare to load parameters
	string folderName=modelName+"/";
	ifstream sourceFile;
	reaction tempReaction;
	compType tempComp;
	
	stack <reaction> reactionBuffer;
	stack <compType> reactantBuffer;

// read reaction data
	string fileName=folderName+REACTIONFILE;
	sourceFile.open(fileName.c_str());
	while(!tempReaction.loadReaction(sourceFile)) 	reactionBuffer.push(tempReaction);
	sourceFile.close();
	nReact=reactionBuffer.size();
	fileName=folderName+REACTANTFILE;
	sourceFile.open(fileName.c_str());
	while(!sourceFile.eof())
	{
		sourceFile>>tempComp;
		if(!sourceFile.eof()) 	reactantBuffer.push(tempComp);
	}
	sourceFile.close();
	nComp=reactantBuffer.size();
	generateModel();
	for (int i=nReact-1;i>=0;i--)
	{
		react[i]=reactionBuffer.top();
		reactionBuffer.pop();
	}
	for (int i=nComp-1;i>=0;i--)
	{
		comp[i]=compBackup[i]=reactantBuffer.top();
		reactantBuffer.pop();
	}
	return 0;
}

//template<typename compType, typename updRateType>
//coarseGrainedModel<compType,updRateType>::coarseGrainedModel()
//{
//	initState=0;
//}

template<typename compType, typename updRateType>
coarseGrainedModel<compType,updRateType>::coarseGrainedModel(string & modelName)
{
	initState=0;
	modelLoading(modelName);
}

template<typename compType, typename updRateType>
coarseGrainedModel<compType,updRateType>::coarseGrainedModel
(const coarseGrainedModel<compType, updRateType> & dummy)
{
	initState=0;
	assign(dummy);
}

template<typename compType, typename updRateType>
coarseGrainedModel<compType,updRateType>::~coarseGrainedModel()
{
	eraseModel();
}

template<typename compType, typename updRateType>
void coarseGrainedModel<compType,updRateType>::reset()
{
	for (int i=0;i<nComp;i++) 	comp[i]=compBackup[i];
	time=lastSavedTime=0;
}

template<typename compType, typename updRateType>
coarseGrainedModel<compType,updRateType> & coarseGrainedModel<compType,updRateType>::
operator=(const coarseGrainedModel<compType, updRateType> & dummy)
{
	if (this == &dummy) return *this;
	else
	{
		assign(dummy);
		return *this;
	}
}


#endif 	//__COARSEGRAINEDCOMMON_H_INCLUDED__
