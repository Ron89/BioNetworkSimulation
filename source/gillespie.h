#ifndef __GILLESPIE_H_INCLUDED__
#define __GILLESPIE_H_INCLUDED__

#include<cstdlib>
#include<cmath>
#include<ctime>
#include<random> 	// Mersenne Twister engine works only for c++11 or above.

#define RDALGORITHM_MT 1 	// 1 is assigned to MT19937 for rd generator
#define SEEDER_DEFAULT 10000

using namespace std;

template<typename modelClassName>
class gillespie
{
private:
// random generator
//	dsfmt_t dsfmt;
	mt19937_64 randomMTwister;
	inline double popRandom()
	{
		if (randomAlgorithm==RDALGORITHM_MT)
			return double(randomMTwister()-randomMTwister.min())/
				double(randomMTwister.max()-randomMTwister.min());
		else
			return drand48();
//		return dsfmt_genrand_close_open(&dsfmt);
//		return randGen.operator();
//
	}

// Required parameter 
	int reactionNumber;
	int (modelClassName::*rateDetermine)(double *, int *); 
								// provided by model, first argument is the address for
								// output reaction rates of each reaction.
	int (modelClassName::*reactantUpdate)(int *, double *); 	
								// second would be number of incicences for
								// each reactions. First would be the output

public:

// random reseeder
	bool randomAlgorithm; 		// 1. use MT19937
								// 0. use rand48
	inline void reseedRandom(int seed)
	{
		if (randomAlgorithm==RDALGORITHM_MT)
			randomMTwister.seed(seed);
		else
			srand48(seed);
//		randGen.seed(seed);
//		dsfmt_init_gen_rand(&dsfmt, seed);
	}
	inline void reseedRandom(int seed1, int seed2) 	// used only with MTwister
	{
		seed_seq sseq{seed1, seed2};
		randomMTwister.seed(sseq);
	}

// algorithm functional parts
	double iterate(int * comp_alias);

// constructor
	gillespie(){
		randomAlgorithm=RDALGORITHM_MT; 
		reseedRandom(SEEDER_DEFAULT);
	}

	gillespie(int reactionNumber_alias, int
			(modelClassName::*rateDetermine_alias)(double *, int *), int
			(modelClassName::*reactantUpdate_alias)(int *, double *), bool
			useMTwister=RDALGORITHM_MT)	
	{
		randomAlgorithm=useMTwister; 
		reactionNumber=reactionNumber_alias;
		rateDetermine=rateDetermine_alias;
		reactantUpdate=reactantUpdate_alias;
		reseedRandom(SEEDER_DEFAULT);
	}
	void assign(const gillespie<modelClassName> & dummy)
	{
		reactionNumber=dummy.reactionNumber;
	}
};

template<typename modelClassName>
double gillespie<modelClassName>::iterate(int * comp_alias)
{
	double r1, r2;
	
	double dt;
	double lambda;
	double * rate=new double[reactionNumber];
	double * lambdaSig=new double[reactionNumber+1];
	int ll,rl,templ;

	r1=popRandom();
	lambdaSig[0]=0;

	(static_cast<modelClassName*>(this)->*rateDetermine)(rate, comp_alias);
	
	for (int i=1;i<=reactionNumber;i++)
	{
		lambdaSig[i]=lambdaSig[i-1]+rate[i-1];
		rate[i-1]=0;
	}
	lambda=lambdaSig[reactionNumber];
	for (int i=1;i<=reactionNumber;i++) 	lambdaSig[i]/=lambda;
	dt=-log(r1)/lambda;

// ensure that the random number generated didn't coincide with any boundary.
	lambda=0;
	while(lambda==0)
	{
		lambda=1;
		r2=popRandom();
		for (int i=0;i<=reactionNumber;i++) 	lambda*=(r2-lambdaSig[i]);
	}

	ll=0;
	rl=reactionNumber;

	while(rl-ll>1)
	{
		templ=(ll+rl)/2;
		if(r2<lambdaSig[templ]) 	rl=templ;
		else 	ll=templ;
	}
	rate[ll]=1;
	(static_cast<modelClassName*>(this)->*reactantUpdate)(comp_alias, rate);

	delete [] rate;
	delete [] lambdaSig;

	return dt;
}

#endif 	//__GILLESPIE_H_INCLUDED__
