#ifndef __DISTRIBUTION_H_INCLUDED__
#define __DISTRIBUTION_H_INCLUDED__
#include <algorithm>
#include <cstdio>
#include <fstream>
#include<unistd.h>
#include<sys/stat.h>

#include "basicDef.h"

#define DISTRIBUTION_SIZECAP 2000000000
#define DIST_ELEMENT_DESCRIPTION "/element_description"
#define DIST_COUNTING "/counting"
#define MAX_VARIABLE 10

using namespace std;

// It works, but code need refine later.
class distribution
{
public:
	int nObserver;
	int observerRange[MAX_VARIABLE]; 	// maximum-minimum, with unit interval as 1.
	int highBound[MAX_VARIABLE];
	int lowBound[MAX_VARIABLE];
	int compressLevel[MAX_VARIABLE]; 	// a compress level for each observer.
	string resultFolder;
	long int allocatedSize;
	long * observer;
	bool flag_Exceed;
//	bool distributionDefined;
	bool sampleTaken;

	~distribution()
	{
		delete [] observer;
	}
	distribution(int nObserver_alias, int * observerRange_alias, string & resultFolder_alias)
	{
		struct stat sb; 	// check existance of folder
		long tempSize=1, tempCompress=1;
		int compressBaseline=1;
		int rangeMax=0;
		fstream filePointer;
		resultFolder=resultFolder_alias;

		sampleTaken=0; 	// so that when the first sample is taken, description can be complete

		nObserver=nObserver_alias;
		
		for (int i=0; i<nObserver; i++)
		{
			observerRange[i]=observerRange_alias[i];
			if (observerRange[i]>rangeMax)
				rangeMax=observerRange[i];
			tempSize*=observerRange[i];
			compressLevel[i]=1;
			tempCompress*=compressLevel[i];
		}
		// determine compress level based on the memory allowed to uptake
		while (tempSize/tempCompress>DISTRIBUTION_SIZECAP)
		{
			compressBaseline*=2;
			tempCompress=1.;
			for (int i=0; i<nObserver; i++)
			{
				compressLevel[i]=max(1, compressBaseline*observerRange[i]/rangeMax);
				while (observerRange[i]%compressLevel[i]!=0)
					compressLevel[i]-=1;
				tempCompress*=compressLevel[i];
			}
		}

		// write the rule to file
		filePointer.open((resultFolder+DIST_ELEMENT_DESCRIPTION).c_str(), ios::out | ios::trunc);
		filePointer<<"# This file contains the description for the distribution.\n";
		filePointer<<"# Observer Number: \n\t"<<nObserver<<endl;
		filePointer<<"# Range:\n";
		for (int i=0;i<nObserver;i++)
			filePointer<<observerRange[i]<<'\t';
		filePointer<<"\n# Compress level: \n";
		for (int i=0;i<nObserver;i++)
			filePointer<< compressLevel[i]<<'\t';
		filePointer<<endl;
		filePointer.close();

		// allocate memory for the distribution
		allocatedSize=tempSize/tempCompress;
		// initialize observer
		observer=new long [allocatedSize];
		for (long i=0; i<allocatedSize; i++) 	observer[i]=0;

		if (stat(resultFolder.c_str(), &sb) != 0) 	mkdir(resultFolder.c_str(),0755);
		flag_Exceed=0;	
	}
	long IDextraction(int * element_alias)
	{
		long tempID=(element_alias[0]%observerRange[0])/compressLevel[0];
		for (int i=1;i<nObserver;i++)
		{
			tempID=tempID*(observerRange[i-1]/compressLevel[i-1])+(element_alias[i]%observerRange[i])/compressLevel[i];
		}
		return tempID;
	}
	inline long insertCounting(int * element_alias)
	{
		checkBoundary(element_alias);
		return observer[IDextraction(element_alias)]++;
	}
	void checkBoundary(int * element_alias)
	{
		int tempNum;
		if (sampleTaken==0) 	// when it's the first data taken
		{
			fstream elementDescription;
			elementDescription.open((resultFolder+DIST_ELEMENT_DESCRIPTION).c_str(), ios::out | ios::app);
			elementDescription<<" sample element:\n\t";
			for (int i=0;i<nObserver;i++)
			{
				highBound[i]=lowBound[i]=element_alias[i];	// initiate highBound and lowBound value;
				elementDescription<<element_alias[i]<<'\t';
			}
			elementDescription<<endl;
			elementDescription.close();
				
			sampleTaken=1;
		}
		else
		{
			
			for (int i=0;i<nObserver;i++)
			{
				tempNum=element_alias[i];
				if (tempNum>highBound[i])
				{
					highBound[i]=tempNum;
					if (highBound[i]-lowBound[i]>=observerRange[i]-1) 	flag_Exceed=1;
				}
				if (tempNum<lowBound[i])
				{
					lowBound[i]=tempNum;
					if (highBound[i]-lowBound[i]>=observerRange[i]-1) 	flag_Exceed=1;
				}
			}
		}
	}
	void saveDistribution()
	{
		FILE * distFile;
		distFile=fopen((resultFolder+DIST_COUNTING).c_str(), "w");
		for (long i=0;i<allocatedSize;i++)
		{
			if (observer[i]!=0)
				fprintf(distFile, "%15ld\t%20ld\n",i,observer[i]);
		}
		fclose(distFile);
	}
};

#endif //__DISTRIBUTION_H_INCLUDED__
