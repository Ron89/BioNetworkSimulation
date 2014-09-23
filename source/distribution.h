#ifndef __DISTRIBUTION_H_INCLUDED__
#define __DISTRIBUTION_H_INCLUDED__
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <sstream>
#include<unistd.h>
#include<sys/stat.h>

#include "basicDef.h"

#define DISTRIBUTION_SIZECAP 4000000000
#define DIST_ELEMENT_DESCRIPTION "/element_description"
#define DIST_COUNTING "/counting"
#define MAX_VARIABLE 10

using namespace std;

// It works, but code need refine later.
// To ensure the efficiency of file writing, I use the C routine rather than the mostly 
// used c++ routine.
class distribution
{
public:
	int nObserver;
	int observerRange[MAX_VARIABLE]; 	// maximum-minimum, with unit interval as 1.
	int observerZero[MAX_VARIABLE]; 	// indicate where the zero value is.
	int highBound[MAX_VARIABLE];
	int lowBound[MAX_VARIABLE];
	int compressLevel[MAX_VARIABLE]; 	// a compress level for each observer.
	int compressLevel_optimized[MAX_VARIABLE]; 	// a compress level for each observer.
	string resultFolder;
	long int allocatedSize;
	long int allocatedSize_optimized;
	long * observer;
	bool flag_Exceed;
	bool flag_Storage; 	// storage allocated

	bool sampleTaken;

	~distribution()
	{
		delStorage();
	}

// storage manipulation
	void newStorage()
	{
		if (flag_Storage) 	delStorage();
		observer=new long [allocatedSize];
		flag_Storage=1;
	}
	void delStorage()
	{
		if (!flag_Storage) 	return;
		delete [] observer;
		flag_Storage=0;
	}
// add an element into the observer(distribution)
	inline long insertCounting(int * element_alias)
	{
		checkBoundary(element_alias);
		return observer[IDextraction(element_alias)]++;
	}

// class constructor for simulation based evaluation
	distribution(int nObserver_alias, int * observerRange_alias, string & resultFolder_alias)
	{
		struct stat sb; 	// check existance of folder
		flag_Storage=0;

		nObserver=nObserver_alias;
		for (int i=0; i<nObserver; i++)
		{
			observerRange[i]=observerRange_alias[i];
		}
		
		initiateScheme(); 	// initiating storage scheme
		optimizeScheme(); 	// optimize compress level of the distribution for memory

		fstream filePointer;
		resultFolder=resultFolder_alias;

		sampleTaken=0; 	// so that when the first sample is taken, description can be complete

		// generate folder if not exist
		if (stat(resultFolder.c_str(), &sb) != 0) 	mkdir(resultFolder.c_str(),0755);
		flag_Exceed=0;	
		// write the storage scheme to file
		filePointer.open((resultFolder+DIST_ELEMENT_DESCRIPTION).c_str(), ios::out | ios::trunc);
		filePointer<<"# This file contains the description for the distribution.\n";
		filePointer<<"# Observer Number: \n\t"<<nObserver<<endl;
		filePointer<<"# Range:\n\t";
		for (int i=0;i<nObserver;i++)
			filePointer<<observerRange[i]<<'\t';
		filePointer<<"\n# Compress level: \n\t";
		for (int i=0;i<nObserver;i++)
			filePointer<< compressLevel[i]<<'\t';
		filePointer<<endl;
		filePointer.close();

		// allocate memory for the distribution
		// initialize observer
		newStorage();
		for (long i=0; i<allocatedSize; i++) 	observer[i]=0;
	}

// class constructor for loading from file
	distribution(int nObserver_alias, int)<++>

// extract ID from element input.
	long IDextraction(int * element_alias)
	{
		long tempID=(element_alias[0]%observerRange[0])/compressLevel[0];
		for (int i=1;i<nObserver;i++)
		{
			tempID=tempID*(observerRange[i]/compressLevel[i])+
			(element_alias[i]%observerRange[i])/compressLevel[i];
		}
		return tempID;
	}

// read ID from ID_alias, put the output into element_alias
	void ID2Index(int * element_alias, long ID_alias)
	{
		long tempID=ID_alias;
		for(int i=nObserver-1;i>=0;i--) 	
		{
			element_alias[i]=observerZero[i]+
			tempID%(observerRange[i]/compressLevel[i])*compressLevel[i]+compressLevel[i]/2;
			tempID/=(observerRange[i]/compressLevel[i]);
		}
	}

// check if the predefined boundary is reached.(before insertCounting)
	void checkBoundary(int * element_alias)
	{
		int tempNum;
		if (sampleTaken==0) 	// when it's the first data taken
		{
			fstream elementDescription;
			elementDescription.open((resultFolder+DIST_ELEMENT_DESCRIPTION).c_str(), ios::out | ios::app);
			elementDescription<<"# sample element:\n\t";
			for (int i=0;i<nObserver;i++)
			{
				highBound[i]=lowBound[i]=element_alias[i];	// initiate highBound and lowBound value;
				elementDescription<<element_alias[i]<<'\t';
				observerZero[i]=element_alias[i]/observerRange[i];
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

// save the current distribution into file
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

// export the distribution in a format of: 	ob1 ob2 ... obn [count]
	void exportDistribution()
	
	
// initialize storage scheme(used when creating new distribution)
	void initiateScheme()
	{
		for (int i=0; i<nObserver; i++)
		{
			compressLevel[i]=1;
		}
	}

// optimize compress level of the distribution for memory
// based on DISTRIBUTION_SIZECAP
// hope in the future can use real memory size to determine it.
// Note that compression is based on loss of information, and lost information can
// never be regen. Thus the compress level can only be raised higher.(currently
// on base 2 basis.) 
	void optimizeScheme(bool replaceCurrentScheme=1)
	{
		int rangeMax=0;
		long tempSize=1, tempCompress=1; 	// 
		int compressBaseline=1; 
		// load current range, resource requirement and compress level.
		for (int i=0; i<nObserver; i++)
		{
			compressLevel_optimized[i]=compressLevel[i]; 	// backup current compress level;

			if (observerRange[i]>rangeMax)	rangeMax=observerRange[i];
			tempSize*=observerRange[i];
			tempCompress*=compressLevel[i];
			if (compressLevel[i]>compressBaseline) 	compressBaseline=compressLevel[i];
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
		if(replaceCurrentScheme==0)
		{
			int temp=1;
			for (int i=0; i<nObserver; i++)
			{
				temp=compressLevel_optimized[i];
				compressLevel_optimized[i]=compressLevel[i];
				compressLevel[i]=temp;
			}
			allocatedSize_optimized=tempSize/tempCompress;
		}
		else
		{
			for (int i=0; i<nObserver; i++)
			{
				compressLevel_optimized[i]=compressLevel[i];
			}
			allocatedSize=allocatedSize_optimized=tempSize/tempCompress;
		}
	}

// load distribution description from file
	void loadScheme()
	{
		fstream elementDescription;
		string func_buffer;
		elementDescription.open((resultFolder+DIST_ELEMENT_DESCRIPTION).c_str(), ios::in);
		while(getline(elementDescription, func_buffer) && func_buffer[0]!='#')
		{
			stringstream ss(func_buffer);
			ss>>Observer;
		}
		while(getline(elementDescription, func_buffer) && func_buffer[0]!='#')
		{
			stringstream ss(func_buffer);
			for (i=0;i<nObserver;i++) 	ss>>observerRange[i];
		}
		while(getline(elementDescription, func_buffer) && func_buffer[0]!='#')
		{
			stringstream ss(func_buffer);
			for (i=0;i<nObserver;i++) 	ss>>compressLevel[i];
		}
		while(getline(elementDescription, func_buffer) && func_buffer[0]!='#')
		{
			stringstream ss(func_buffer);
			for (i=0;i<nObserver;i++) 
			{
				ss>>observerZero[i];
				observerZero[i]=(observerZero[i]/observerRange[i])*observerRange[i];
			}
		}
		elementDescription.close();
	}
};

#endif //__DISTRIBUTION_H_INCLUDED__
