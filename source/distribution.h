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
#define DIST_EXCEPTION "/exception"

using namespace std;

struct element_count
{
	long ID;
	long count;
	inline void assign(element_count & dummy)
	{
		ID=dummy.ID;
		count=dummy.count;
	}
	element_count & operator=( element_count & dummy)
	{
		if (this!=&dummy) 	assign(dummy);
		return *this;
	}
};


// It works, but code need refine later.
class distribution
{
public:
	int nObserver;
	int * observerRange; 	// maximum-minimum, with unit interval as 1.
	int * compressLevel; 	// a compress level for each observer.
	string resultFolder;
	long int allocatedSize;
	element_count * observer;
	long int nElement_filled;
	bool flag_Exceed;
//	bool distributionDefined;
	bool sampleTaken;

	~distribution()
	{
		delete [] compressLevel;
		delete [] observerRange;
		delete [] observer;
	}
	distribution(int nObserver_alias, int * observerRange_alias, string & resultFolder_alias)
	{
		struct stat sb; 	// check existance of folder
		double tempSize=1, tempCompress=1;
		int compressBaseline=1;
		int rangeMax=0;
		fstream filePointer;

		sampleTaken=0; 	// so that when the first sample is taken, description can be complete

		nObserver=nObserver_alias;
		compressLevel=new int [nObserver];
		observerRange=new int [nObserver];
		
		for (int i=0; i<nObserver; i++)
		{
			observerRange[i]=observerRange_alias[i];
			if (observerRange[i]>rangeMax)
				rangeMax=observerRange[i];
			tempSize*=double(observerRange[i]);
			compressLevel[i]=1;
			tempCompress*=double(compressLevel[i]);
		}
		// determine compress level based on the memory allowed to uptake
		while (tempSize/tempCompress>DISTRIBUTION_SIZECAP)
		{
			compressBaseline*=2;
			tempCompress=1.;
			for (int i=0; i<nObserver; i++)
			{
				compressLevel[i]=max(1, compressBaseline*observerRange[i]/rangeMax);
				tempCompress*=double(compressLevel[i]);
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
		allocatedSize=(long) (tempSize/tempCompress);
		// initialize observer
		observer=new element_count [allocatedSize];
		for (long i=0; i<allocatedSize; i++) 	observer[i].ID=observer[i].count=0;
		nElement_filled=0;

		resultFolder=resultFolder_alias;
		if (stat(resultFolder.c_str(), &sb) != 0) 	mkdir(resultFolder.c_str(),0755);
	// initiate necessary files;
//		filePointer.open((resultFolder+DIST_ELEMENT_DESCRIPTION).c_str(), ios::out | ios::trunc);
//		filePointer.close();
		filePointer.open((resultFolder+DIST_EXCEPTION).c_str(), ios::out | ios::trunc);
		filePointer.close();
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
	long insertCounting(int * element_alias)
	{
		long ll=0, rl=nElement_filled, templ;
		long tempID=IDextraction(element_alias);
		if (nElement_filled==0)
		{
			insertElement(element_alias,0);
			return 1;
		}
		else if (tempID>=observer[ll].ID)
		{
			while (rl-ll>1)
			{
				templ=(ll+rl)/2;
				if (tempID<observer[templ].ID) rl=templ;
				else 	ll=templ;
				if (tempID==observer[templ].ID) break;
			}
			if (observer[ll].ID==tempID) 	
			{
				return observer[ll].count++;
			}
			else 
			{
				insertElement(element_alias,ll+1);
				return 1;
			}
		}
		else
		{
			insertElement(element_alias,0);
			return 1;
		}
	}
	void insertElement(int * element_alias, long position)
	{
		fstream elementDescription;
		long ID=IDextraction(element_alias);
		if (nElement_filled<allocatedSize)	
		{
			// record the first sample element in element_description file
			if (sampleTaken==0)
			{
				elementDescription.open((resultFolder+DIST_ELEMENT_DESCRIPTION).c_str(), ios::out | ios::app);
				elementDescription<<" sample element:\n\t";
				for (int i=0;i<nObserver;i++)
					elementDescription<<element_alias[i]<<'\t';
				elementDescription<<endl;
				elementDescription.close();
				sampleTaken=1;
			}

			// insert this new element in the observer 
			for (long i = nElement_filled-1; i>=position;i--) 	observer[i+1]=observer[i];
			observer[position].ID=ID;
			observer[position].count=1;
			nElement_filled++;
		}
		else
		{
			flag_Exceed=1;
			elementDescription.open((resultFolder+DIST_EXCEPTION).c_str(), ios::out | ios::app);
			for (int i=0;i<nObserver;i++) 	elementDescription<<element_alias[i]<<'\t';
			elementDescription<<endl;
			elementDescription.close();
		}
	}
	void saveDistribution()
	{
		FILE * distFile;
		distFile=fopen((resultFolder+DIST_COUNTING).c_str(), "w");
		for (long i=0;i<nElement_filled;i++)
		{
			fprintf(distFile, "%15ld\t%20ld\n",observer[i].ID,observer[i].count);
		}
		fclose(distFile);
	}
};

#endif //__DISTRIBUTION_H_INCLUDED__
