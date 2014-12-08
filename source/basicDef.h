#ifndef __BASICDEF_H__
#define __BASICDEF_H__
//====================
//network source file name
#define RATEFILE "rate"
#define INITCONDFILE "initCond"
#define RATEMATRIXFILE "rMatrix"
#define UPDATEMATRIXFILE "updMatrix"
// simulation result
#define RESULTFOLDER "./result/"
#define RESULTNAME "sim"

// ideal system size(only for historical reason)
#define MAXNCOMP 100
#define MAXNRATE 200
// typename abbreviation
#define TYPENAME compType,rateType,powerType

// coefficients for Runge-Kutta method
#define RKA2 1./5.
#define RKA3 3./10.
#define RKA4 3./5.
#define RKA5 1.
#define RKA6 7./8.
#define RKB21 1./5.
#define RKB31 3./40.
#define RKB32 9./40.
#define RKB41 3./10.
#define RKB42 -9./10.
#define RKB43 6./5.
#define RKB51 -11./54. 
#define RKB52 5./2.
#define RKB53 -70./27.
#define RKB54 35./27.
#define RKB61 1631./55296.
#define RKB62 175./512.
#define RKB63 575./13824.
#define RKB64 44275./110592.
#define RKB65 253./4096.
#define RKC01 37./378.
#define RKC02 0.
#define RKC03 250./621.
#define RKC04 125./594.
#define RKC05 0.
#define RKC06 512./1771.
#define RKC11 2825./27648.
#define RKC12 0.
#define RKC13 18575./48384.
#define RKC14 13525./55296.
#define RKC15 277./14336.
#define RKC16 1./4.

// accuracy requirement used in step size adjustment
#define ACCURACY 0.000001

// coarse grained modelling
#define MAXRATECOEF 4
#define MAXDEPENDENCY 3
#define MAXUPDSET 3

// model defining file
#define REACTIONFILE "reaction"
#define REACTANTFILE "reactant"

// distribution analysis(due to efficiency concern and
// memory limitation), used mainly for distribution.h
// and some related applications
#define MAX_VARIABLE 10

#endif 	//__BASICDEF_H__
