#ifndef __DISTRIBUTION_H_INCLUDED__
#define __DISTRIBUTION_H_INCLUDED__
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include<unistd.h>
#include<sys/stat.h>

#include "basicDef.h"

#define DISTRIBUTION_SIZECAP 4000000000
#define DIST_ELEMENT_DESCRIPTION "/element_description"
#define DIST_BOUNDARY_DESCRIPTION "/boundary_description"
#define DIST_COUNTING "/counting"
#define RESCUE_DISTRIBUTION_NAME "temp_dist"

using namespace std;

// It works, but code need refine later.
// To ensure the efficiency of file writing, I use the C routine rather than the mostly 
// used c++ routine.
class distribution
{
public:
	int nObserver;
	int observerRange[MAX_VARIABLE]; 	// expected maximum-minimum, with unit interval as 1.
	int observerZero[MAX_VARIABLE]; 	// indicate where the zero value is.
	int highBound[MAX_VARIABLE]; 		// real maximum obtained from data.
	int lowBound[MAX_VARIABLE]; 		// real minimum obtained from data.
// * 	highBound and lowBound are either loaded from file, initialized by data or
// 		manually determined. Thus no additional initialization is required.
	int compressLevel[MAX_VARIABLE]; 	// a compress level for each observer.
	int compressLevel_optimized[MAX_VARIABLE]; 	// a compress level for each observer.
	string resultFolder;
	fstream countFile;
	long int allocatedSize;
	long int allocatedSize_optimized;
	long * observer;
	bool flag_Exceed;

	bool flag_Storage; 		// storage allocated
	bool countFile_open; 	// count file opened

	bool sampleTaken;

	~distribution()
	{
		delStorage(); 			// free the memory used for the storage.
		unloadDistribution(); 	// unload count file.
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

// constructors
// create an empty distribution(default constructor)
	distribution()
	{
		// initiate flags;
		flag_Storage=0;	
		flag_Exceed=0;
		countFile_open=0;
		sampleTaken=0; 	// so that when the first sample is taken, description can be complete
	}

// class constructor for simulation based evaluation
	distribution(int nObserver_alias, int * observerRange_alias, string & resultFolder_alias)
	{
// initiate flags
		flag_Storage=0;
		countFile_open=0;
		flag_Exceed=0;	
		sampleTaken=0; 	// so that when the first sample is taken, description can be complete

		nObserver=nObserver_alias;
		for (int i=0; i<nObserver; i++)
		{
			observerRange[i]=observerRange_alias[i];
		}
		
		initiateScheme(); 	// initiating storage scheme
		optimizeScheme(); 	// optimize compress level of the distribution for memory
	// set result folder to the user specified one.
		resultFolder=resultFolder_alias;

		// allocate memory for the distribution
		// initialize observer
		newStorage();
		for (long i=0; i<allocatedSize; i++) 	observer[i]=0;
	}

// class constructor for loading from file
// for memory space reservation, we won't load the distribution right after
// the constructor is loaded. Users may use 
// 		* pop_element(ID_alias, count_alias)
// 		* load_distribution()
// to load the distribution for different purpose(no observer array created)
	distribution(string & resultFolder_alias)
	{
		// initiate flags;
		flag_Storage=0;	
		flag_Exceed=0;
		countFile_open=0;
		sampleTaken=0; 	// so that when the first sample is taken, description can be complete

		resultFolder=resultFolder_alias;
		openDistribution();
	}

// assign a distribution with certain scheme.
	void assign_scheme(int nObserver_alias, int * observerRange_alias, int * observerZero,
	int * highBound_alias,
	int * lowBound_alias, int * compressLevel_alias)
	{
		nObserver=nObserver_alias;
		for (int i=0; i<nObserver; i++)
		{
			observerRange[i]=observerRange_alias[i];
			observerZero[i]=observerRange_alias[i];
			highBound[i]=highBound_alias[i];
			lowBound[i]=lowBound_alias[i];
			compressLevel[i]=compressLevel_alias[i];
		}
	}

// open distribution for reading and writing.
// Note: 	write=0 uses void loadScheme() to load saving scheme. And open counting file
// 			as fstream countFile for future processing.
// Note: 	so far, only write=0 mode is used. Thus write=1 scheme is not complete. 
// 			Writing distribution is currently using 
// 			void saveDistribution() 	save as databulk.
// 			void exportDistribution() 	export into recognizable file.
// 				09.30.14
	void openDistribution(bool write=0)
	{
		if(countFile_open) 	unloadDistribution();
		if(resultFolder.length()==0)
		{
			if (write==1)
			{
				cout<<"\n Warning: Undefined distribution name. Rescue name:"<<RESCUE_DISTRIBUTION_NAME<<endl;
				resultFolder=RESCUE_DISTRIBUTION_NAME;
			}
			else
			{
				cout<<"\n Error: Undefined distribution name, unable to load. \n";
				exit(1);
			}
		}
		if(write==0)
		{
			loadScheme();
			countFile.open(resultFolder+DIST_COUNTING, ios::in);
		}
		if(write==1)
		{
			countFile.open(resultFolder+DIST_COUNTING, ios::out | ios::trunc);	
		}
		countFile_open=1;
		sampleTaken=1; 	// correctly saved distribution always have samples stored within it.
	}
	void unloadDistribution()
	{
		if(countFile_open) 	countFile.close();
		countFile_open=0;
		sampleTaken=0; 	
		initiateScheme(); 	// erase earlier defined scheme when processing model.
	}

// method of reading counting file(DO NOT read manually.)
// if file is not open return 1, as well as when EOF is triggered.
	bool readCounting(long * index_readout, long * count_readout)
	{
		if (countFile_open==1)
		{
			countFile>>*index_readout;
			countFile>>*count_readout;
			return countFile.eof();
		}
		else return 1;
	}
// all flags will be restored to normal, fstream will go to the begining of the file
	void returnFileBegining()
	{
		countFile.clear();
		countFile.seekg(0, ios::beg); 	// essentially the begining, equivalent to seekg(0);
	}

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
			if(element_alias[i]<lowBound[i]) 	element_alias[i]+=observerRange[i];
		}
	}

// check if the predefined boundary is reached.(before insertCounting)
	void checkBoundary(int * element_alias)
	{
		int tempNum;
		if (sampleTaken==0) 	// when it's the first data taken
		{
			for (int i=0;i<nObserver;i++)
			{
				highBound[i]=lowBound[i]=element_alias[i];	// initiate highBound and lowBound value;
			}
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
		struct stat sb; 	// check existance of folder
		FILE * distFile;
		fstream filePointer;

	// generate folder if not exist
		if (stat(resultFolder.c_str(), &sb) != 0) 	mkdir(resultFolder.c_str(),0755);
	// save description
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
		filePointer<<"# boundary description:\n";
		for (int i=0;i<nObserver;i++)
		{
			filePointer<<"#\tobserver_"<<i<<endl<<"\t\t"
			<<lowBound[i]<<"\t"<<highBound[i]<<endl;
		}
		filePointer.close();
	// save distribution;
		distFile=fopen((resultFolder+DIST_COUNTING).c_str(), "w");
		for (long i=0;i<allocatedSize;i++)
		{
			if (observer[i]!=0)
				fprintf(distFile, "%15ld\t%20ld\n",i,observer[i]);
		}
		fclose(distFile);
	}

// export the distribution in a format of: 	ob1 ob2 ... obn [count]
// currently only support exporting from in memory distribution
// Result is the same as saveDistribution, but the ID is replaced by
// the amount of each observers.
	void exportDistribution(string & export_file_name)
	{
		int index_temp[MAX_VARIABLE];
		FILE * exportFilePointer;

		exportFilePointer=fopen(export_file_name.c_str(), "w");
		for (long i=0;i<allocatedSize;i++)
		{
			if(observer[i]!=0)
			{
				ID2Index(index_temp,i);
				for (int j=0;j<nObserver;j++)
				{
					fprintf(exportFilePointer, "%8d\t",index_temp[j]);	
				}
				fprintf(exportFilePointer, "%20ld\n",observer[i]);
			}
		}
		fclose(exportFilePointer);
	}
	
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
		stringstream ss;
		fstream elementDescription;
		string func_buffer;
		int tempSize=1;
		int tempCompress=1;

		elementDescription.open((resultFolder+DIST_ELEMENT_DESCRIPTION).c_str(), ios::in);
		do
		{
			if(!(getline(elementDescription, func_buffer))) 	exit(1);
		}
		while(func_buffer[0]=='#');
		ss<<func_buffer;
		ss>>nObserver;
		ss.str("");
		ss.clear();
		do
		{
			if(!(getline(elementDescription, func_buffer))) 	exit(1);
		}
		while(func_buffer[0]=='#');
		ss<<func_buffer;
		for (int i=0;i<nObserver;i++)
		{
			ss>>observerRange[i];
			tempSize*=observerRange[i];
		}
		ss.str("");
		ss.clear();
		do
		{
			if(!(getline(elementDescription, func_buffer))) 	exit(1);
		}
		while(func_buffer[0]=='#');
		ss<<func_buffer;
		for (int i=0;i<nObserver;i++)
		{
			ss>>compressLevel[i];
			tempCompress*=compressLevel[i];
		}
		ss.str("");
		ss.clear();
		for (int i=0;i<nObserver;i++)
		{
			do
			{
				if(!(getline(elementDescription, func_buffer))) 	exit(1);
			}
			while(func_buffer[0]=='#');
			ss<<func_buffer;
			{
				ss>>lowBound[i];
				ss>>highBound[i];
			}
			ss.str("");
			ss.clear();
		}
		elementDescription.close();
	// determine baseline with lowBound[i];
		for (int i=0; i<nObserver;i++)
		{
			observerZero[i]=(lowBound[i]/observerRange[i])*observerRange[i];
		}
	// determine the memory size need to be allocated.
		allocatedSize=tempSize/tempCompress;
	}
};

#endif //__DISTRIBUTION_H_INCLUDED__
