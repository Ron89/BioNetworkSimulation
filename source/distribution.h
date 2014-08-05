#ifndef __DISTRIBUTION_H_INCLUDED__
#define __DISTRIBUTION_H_INCLUDED__
#include <algorithm>
#include <cstdio>
#include <fstream>
#include<unistd.h>
#include<sys/stat.h>
#include "basicDef"

#define DISTRIBUTION_SIZECAP 1000000000
#define DIST_ELEMENT_DESCRIPTION "/element_description"
#define DIST_COUNTING "/counting"
#define DIST_EXCEPTION "/exception"

using namespace std;

struct element_count
{
	unsigned long ID;
	unsigned long count;
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
	unsigned long int allocatedSize;
	element_count * observer;
	unsigned long int nElement_filled;
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
		double tempSize=1, tempCompress=1;
		int compressBaseline=1;
		int rangeMax_arg=0;
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
			{
				rangeMax=observerRange[i];
				rangeMax_arg=i;
			}
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
				compressLevel[i]*=max(1, compressBaseline*observerRange[i]/rangeMax);
				tempCompress*=double(compressLevel[i]);
			}
		}

		// write the rule to file
		filePointer.open((resultFolder+DIST_ELEMENT_DESCRIPTION).c_str(), ios::out | ios::trunc);
		fprintf(filePointer, "# This file contains the description for the distribution.\n");
		fprintf(filePointer, "# Observer Number: \n\t%d\n",nObserver);
		fprintf(filePointer, "# Range:\n");
		for (int i=0;i<nObserver;i++)
			fprintf(filePointer, "%d\t", observerRange[i]);
		fprintf(filePointer, "\n# Compress level: \n")
		for (int i=0;i<nObserver;i++)
			fprintf(filePointer, "%d\t", compressLevel[i]);
		filePointer<<endl;
		filePointer.close();

		// allocate memory for the distribution
		allocatedSize=unsigned long (tempSize/tempCompress);
		// initialize observer
		observer=new element_count [allocatedSize];
		for (int i=0; i<allocatedSize; i++) 	observer[i].ID=observer[i].count=0;
		nElement_filled=0;

		resultFolder=resultFolder_alias;
		if (stat(resultFolder.c_str(), &sb) != 0) 	mkdir(saveName.c_str(),0755);
	// initiate necessary files;
		filePointer.open((resultFolder+DIST_ELEMENT_DESCRIPTION).c_str(), ios::out | ios::trunc);
		filePointer.close();
		filePointer.open((resultFolder+DIST_EXCEPTION).c_str(), ios::out | ios::trunc);
		filePointer.close();
		bool flag_Exceed=0;	
	}
	unsigned long IDextraction(int * element_alias)
	{
		unsigned long tempID=(element_alias[0]%observerRange[0])/compressLevel[0];
		for (int i=1;i<nObserver;i++)
		{
			tempID=tempID*(observerRange[i-1]/compressLevel[i-1])+(element_alias[i]%observerRange[i])/compressLevel[i];
		}
		return tempID;
	}
	unsigned long insertCounting(int * element_alias)
	{
		unsigned long ll=0, rl=nElement_filled, templ;
		unsigned long tempID=IDextraction(element_alias);
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
	void insertElement(int * element_alias, unsigned long position)
	{
		fstream elementDescription;
		unsigned long ID=IDextraction(element_alias);
		if (nElement_filled<allocatedSize)	
		{
			// record the first sample element in element_description file
			if (sampleTaken==0)
			{
				element_description.open((resultFolder+DIST_ELEMENT_DESCRIPTION).c_str(), ios::out | ios::app);
				fprintf(element_description,"# sample element:\n\t");
				for (int i=0;i<nObserver;i++)
					fprintf(element_description,"%d\t",element_alias[i]);
				element_description.close();
				sampleTaken=1;
			}

			// insert this new element in the observer 
			for (int i = nElement_filled-1; i>=position;i--) 	observer[i+1]=observer[i];
			observer[position].ID=ID;
			observer[position].count=1;
			nElement_filled++;
		}
		else
		{
			flag_Exceed=1;
			elementDescription.open((resultFolder+DIST_EXCEPTION).c_str(), ios::out | ios::app);
			for (int i=0;i<nObserver;i++) 	element_description<<element_alias[i]<<'\t';
			element_description<<endl;
			element_description.close();
		}
	}
	void saveDistribution()
	{
		fstream distFile;
		distFile.open((resultFolder+DIST_COUNTING).c_str(), ios::out | ios::trunc);
		for (int i==0;i<nElement_filled;i++)
		{
			fprintf(distFile, "%u\t%u\n",observer[i].ID,observer[i].count);
		}
		distFile.close();
	}
};

#endif //__DISTRIBUTION_H_INCLUDED__
