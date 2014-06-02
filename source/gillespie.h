#ifndef __GILLESPIE_H_INCLUDED__
#define __GILLESPIE_H_INCLUDED__

#include<cstdlib>
#include<cmath>
#include<ctime>

using namespace std;

template<typename modelClassName>
class gillespie
{
private:
// random generator
//	dsfmt_t dsfmt;
	inline void reseedRandom(int seed)
	{
		srand48(seed);
//		randGen.seed(seed);
//		dsfmt_init_gen_rand(&dsfmt, seed);
	}

	inline double popRandom()
	{
		return drand48();
//		return dsfmt_genrand_close_open(&dsfmt);
//		return randGen.operator();
	}

// Required parameter 
	int reactionNumber;
	int (modelClassName::*rateDetermine)(double *); 
								// provided by model, first argument is the address for
								// output reaction rates of each reaction.
	int (modelClassName::*reactantUpdate)(const double *); 	
								// second would be number of incicences for
								// each reactions. First would be the output

public:
// constructor
	gillespie(){}

	gillespie(int reactionNumber_alias, int (modelClassName::*rateDetermine_alias)(double *),
			int (modelClassName::*reactantUpdate_alias)(const double *))	
	{
		reactionNumber=reactionNumber_alias;
		rateDetermine=rateDetermine_alias;
		reactantUpdate=reactantUpdate_alias;
	}

// algorithm functional parts
	double iterate();
};

template<typename modelClassName>
double gillespie<modelClassName>::iterate()
{
	double r1, r2;
	
	double dt;
	double lambda;
	double * rate=new double[reactionNumber];
	double * lambdaSig=new double[reactionNumber+1];
	int ll,rl,templ;

	r1=popRandom();
	lambdaSig[0]=0;

	(static_cast<modelClassName*>(this)->*rateDetermine)(rate);
	
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
	(static_cast<modelClassName*>(this)->*reactantUpdate)(rate);
	delete [] rate;
	delete [] lambdaSig;

	return dt;
}

#endif 	//__GILLESPIE_H_INCLUDED__
