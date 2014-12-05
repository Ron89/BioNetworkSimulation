# Documentary for BCNetwork simulation related Codes.



Purpose of this article is to create an index for the classes and
modules defined in the projects. And guide the future development based
on the work.

Note that the css file used for this README article is one of the css
files hosted on [Github](https://github.com/jasonm23/markdown-css-themes
"jasonm23/markdown-css-themes").

## Biochemical models storage and conversion

### BCNetwork Class Template

defined in "BCNetwork.h"

```cpp
template<typename compType, typename rateType, typename powerType>
class BCNetwork
```

This module is the basic module to store the information necessary to
describe a biochemical network. It is supposed to be able to handle
all types of biochemical networks(from elemental to huge, both discrete
and continuous.)

typename | meaning
--- | ---
`compType` | data type of the amount of reactants
`rateType` | data type of the reactant constants
`powerType` | data type of both rate matrix and update matrix

#### Public variable members

 `type variable_name` |	meaning
------- | ----------------------
 `double t` | current time of the network 
`int nComp` | number of reactants
`compType * comp` | pointer to datablock that store amount of each reactants
`compType * compBackup` | same as `comp`, a back up datablock for resetting purpose
`int nRate` |  number of rate constants
`rateType * rate` | pointer to datablock that store each rate constants
`powerType * rateMatrix` | pointer to datablock for rate matrix. Format: `i*nComp+j`
`powerType * updateMatrix` | pointer to datablock for update matrix. Format, as above.


#### Public function members

```cpp
int loadParameter(vector<compType> &, vector<rateType> &, powerType *, powerType *);
```

Load parameter from a certain model. Four types of parameters required to fully 
describe a model are required. the four inputs are respectively : amount of reactants 
as a form of
vector; amount of rate constants as a form of vector; rate matrix as a datablock of
same format as described in variable members; update matrix as a datablock of same 
requirement.

loadParameter will import all the information obtained from the
input into the inbuilt dataset. So the destruction of input source after the 
import will not disturb the funtion of this modules or any modules that inherit from 
it. 

Importing a second network will overwrite the previous one stored in. The previous 
network will be
automatically completely erased before a second network is written in, so backup
 is suggested before this action.

---

```cpp
int fileOpen(string & condition);
int fileClose();
```

Open and close a data saving file. `fileOpen(string &)` will open a data file with
name "./result/sim$(condition)". New file will be created if it is not exist. And
older data willed be wiped in the file with a same name. Note that new folder won't
be created if not exist. So a folder "result" created under the working directory
before running the code is suggested.

Note that the folder name "./result/" is defined in "basicDef.h" as 

	#define RESULTFOLDER "./result/"

---

```cpp
virtual void saveData();
```

Essentially this function will toss a copy of the current time and amount of each 
reactants into the data file opened by `int fileOpen(string &)`. Data will be saved
in a file named "./result/simEmergency" if no file is opened. But open a file before
saving the data is still suggested.

Note that this member function covers only the basic saving feature and changable in 
classes inherit this class template.

---

```cpp
virtual void reset();
```

Restore the datablock that `* comp` pointing to with the datablock pointed by 
`* compBackup`. And again, this member function is preliminary and changable in
classed derived from this class template.

#### Constructor and Destructor

```cpp
BCNetwork();
BCNetwork(vector<compType> & initComp, vector<rateType> & inputRate, 
             powerType * inputRateMatrix, powerType * inputUpdateMatrix);
~BCNetwork();
```

`BCNetwork()` only setting private members to correct initial value for the class
to function correctly.
`BCNetwork(vector<compType> &, vector<rateType> &, powerType *, powerType *)` will
call `int loadParameter(...)` to load a model when initiating the class.
`~BCNetwork()` will free all memory that allocated from this class template before
its destruction.

### modelLoadder Class Template

defined in "modelLoader.h"

```cpp
template <typename compType, typename rateType, typename powerType >
class modelLoader
```

This class template is merely used to load a model from file. The use of this class
is temporary before it's function incorporated into `BCNetwork` class template.

#### public function members

```cpp
void loadParameter(string & mName);
```

This is the only functional member. It's used to load model parameters stored in folder
"./mName". The format of model parameters will be specified in a model creation/modification
project. To completely load a model, four files will be read. The name of them are 
specified as macros stored in "basicDef.h":

```cpp
#define RATEFILE "rate" 	// reaction rate constants, first line is the number of reactions
#define INITCONDFILE "initCond" 	//initial condition for each reactants, first line is  
									//the number of reactants
#define RATEMATRIXFILE "rMatrix" 	//rate matrix, same format as in datablock
#define UPDATEMATRIXFILE "updMatrix" 	//update matrix, same format as in datablock
```

#### Constructor/Destructor

```cpp
modelLoader(string & mName);
~modelLoader();
```

`modelLoader(string &)` calls `loadParameter(string &)`. This is the only valid way of
initiating modelLoader class. Other ways of initiating this class may involve unexpected
behavior.

### ODENetwork Class

defined in "ODECommon.h"

```cpp
class ODENetwork : public BCNetwork<double, double, int>
```

This is the basic class for biochemical networks described by ODEs. In addition to
the features already defined in `BCNetwork`class template, the ODEs required
for simulation and time variable will be defined in this class. But neither simulation
nor analysis feature is defined in this class. 

These chemical networks are defined in continuous space and are deterministic. Thus 
both reactant amount datatype and reaction rate datatype are double.

Note, the two matrices in most cases are integer, thus the powerType is defined so.
For specific cases that require non-integer powerType, change the type accordingly.

#### Protected member

`type variable_name` | meaning
--- | ---
`double h0` | initial size of time step

#### Public function members

```cpp
void ODETimeDeri(double * timeDeri, double * component);
```

This is the ODEs generated from the biochemical network. The usage of the ODEs is
to feed 2 pointers to the function. `* component` point to the datablock of 
current reactant amounts; `*timeDeri` point to the datablock that stores the time
derivations of each reactants' amount.

Note that the datablock pointed by timeDeri will be changed after this function is
called. Backup is suggested if there is anything you don't want to be erased there.

#### Constructor/ Destructor

```cpp
ODENetwork(vector<double> & initComp, vector<double> & inputRate, 
            int * inputRateMatrix, int * inputUpdateMatrix, double                              
            initTimeStep) :   
            BCNetwork<double,double,int>
            (initComp, inputRate, inputRateMatrix, inputUpdateMatrix);
```

Well, basically it's self explined. It does two things. It calls the constructor that
loads a biochemical network, and give `h0` an initial value given by user.

## Coarse Grained Model module

### Reaction network format

Reaction network will be stored as a folder named by the network name. In the folder, there
will be 2 files, `reaction` and `reactant`. 
In `reactant`, initial values of all relative reactants is to be stored respectively. Number
of values will be treated as the number of reactants.
In `reaction`, all reactions in the network is to be stored respectively. Format of reactions
is as followed:

`code` `rate_0` `rate_1` `rate_2` `rate_3` `dep_0` `dep_1` `dep_2` `upd_num` `upd_0` `upd_1` `upd_2`

where each `xxx` represent a value/group of values, either integer or decimal. Meaning of them are as following:

name | meaning
--- | ---
`code`| `int`, used to identify reaction type 
`rate_0~3`| `float`, store reaction constants. Depending on reaction type, not all 4 of them are necessarily meaningful.
`dep_0~2` | `int`, store reaction rate dependency on reactants. Each int represents the serial of the respective reactant. Depending on reaction type, not all 3 of them are necessarily meaningful.
`upd_num` | `int` ranged from `0~3`, number of reactants to be updated. 
`upd_0~2` | `int <tab> double`, store the serial number and number of  molecules change by each instance of reaction.  Dependingon `upd_num`, not all 3 of them are necessarily meaningful.


Though not necessary, 
it is recommended that one line is to store the information of one reaction. 

In coarse grained model, we describe the reactions into 5 types. Those reaction
types should be able to describe almost all reactions seen in reaction networks, the 
following table is the storage format of each reaction type. 

type | `code` | `rate_0~3` | `dep_0~2` | `upd_num` | `upd_0~2`
--- | --- | --- | --- | --- | ---
`*->..`|0|`k rand rand rand`|`rand rand rand`|`r. spec.` | `r. spec.`
`A->..`|1|`k rand rand rand`|`A rand rand`|`r. spec.` | `r. spec.`
`A+B->..`|2|`k rand rand rand`|`A B rand`|`r. spec.` | `r. spec.`
Hill equa.|3|`k K n rand`|`P factor rand` | `r. spec.` | `r. spec.`
MM kinetics|4|`k K rand rand`|`S E rand`|`r. spec.` | `r. spec.`

Note:

	1. 	Update reactants are case by case specific, thus not described here.

	2. 	The four cases defined here are the four standard cases defined in 
		`class coarseGrainedModel<compType, updRateType>`.

 	3. 	Hill equation usually follows some other reaction, here we define only the case
		that combined protein being able to trigger some type `1` reaction with rate `k`.
		The equilibrium constant of factor combining is `K`.

	4. 	New reaction types can be added, as described in
		`class coarseGrainedModel<compType, updRateType>`.

### reaction structure

defined in "coarseGrainedCommon.h"

```cpp
struct reaction
```

This structure is to store the information used to describe one reaction in coarse 
grained models. 

#### Public variable members

Note: It is not suggested for reactions to be modified after loaded from file, but on
necessary occations, the following will serve a guide on the format information stored
in the structure.

`type variable_name` | meaning
--- | ---
`int code` | code to identify the reaction type
`double rate[MAXRATECOEF]` | rate coefficients of the reaction
`int dependency[MAXDEPENDENCY]` | reactant dependency of the reaction rate
`int updateNumber` | number of reactant to be updated 
`int updateSet[MAXUPDSET]` | the set of reactants to be updated by the reaction
`double updateOrder[MAXUPDSET]` | `updateOrder[i]` gives the amount of reactant change  for reactant `updateSet[i]`.

Note: `updateOrder[i], updateOrder[i]` together is the ith instance of `upd_0~2` mentioned
in section [Reaction network format](#reaction-network-format).

#### Public function members

```cpp
reaction & operator=(const reaction & dummy);
```

The function will copy all values stored in dummy to the current `reaction`
instance and return the pointer of current instance.

---

```cpp
bool loadReaction(ifstream & reactionFile);
```

The function member will read 12 values from `reactionFile`. And store the information in
the current `reaction` instance in a format described in section 
[Reaction network format](#reaction-network-format). Returned value is the EOF flag
of file pointer `reactionFile`.

---

```cpp
void erase()
```

reset all values in this `reaction` instance to `0`.

#### constructors

```cpp
reaction()
reaction(const reaction & dummy)
```

The default constructor creates a all zero valued instance of `reaction`. The copy
constructor copies all values stored in dummy to the current `reaction` instance.

### trajectory structure template

defined in "coarseGrainedOperation.h"

```cpp
template<typename compType>
struct trajectory
```

This structure stores the trajectory of a certain network within a period of time.
It is defined before or during a simulation. A saving member is also defined so that
the data can be written into file whenever needed.

#### Public variable members

`type variable_name` | meaning
--- | ---
`int nComp` | number of reactants in the network.
`double * time` | time points on the trajectory;
`compType comp` | storage for number of each reactant at each time point. Format: `comp[t*nComp+i]`.
`long trajectoryPointer` | writing pointer, indicate the location to write  when using `void assign(int, long)` to write new data or how many data values to save when using  `void save(string, bool)` to write data to disk.

#### Public function members

```cpp
void append(double time_alias, compType * comp_alias);
```

Though it's OK for user to write trajectory data by themselves, this function will append
data to the end of the last written data block, and push `trajectoryPointer` one unit further.

---

```cpp
void save( string & outputFileName, bool append=1);
```

Saving the data currently written to the memory block to the point where `trajectoryPointer`
indicates to a file, named with `outputFileName`. If `append=1` as default, the data will
be appended at the end of the file. While if indicating `append=0`, the original file will
be replaced without warning.

---

```cpp
void erase();
```

This funtion will erase the current assigned memory if exists. `trajectoryPointer` will also
point to 0.

---

```cpp
void reallocate(int nComp_alias, unsigned long trajectorySize);
```

Erase(with `erase()`) currently assigned memory if exists, and reallocate memory for `*comp` 
with size`nComp_alias*trajectorySize`,
and `*time` with size `trajectorySize`. All values of newly assigned memory will be zero.



#### Constructors

```cpp
trajectory();
trajectory(int nComp_alias, unsigned long trajectorySize_alias);
trajectory(trajectory<compType> & dummy);
```

Default constructor create an empty `trajectory` structure. To use, user must first 
use `void reallocate(int nComp_alias, unsigned long trajectorySize)` defined as a public function
member to allocate memory for the trajectory.

`trajectory(int nComp_alias, unsigned long trajectorySize_alias);` automatically calls
`void rellocate(int, unsigned long)`, thus the datablock will be ready for all operations
when this constructor is used.

Copy constructor will create an exact duplicate of dummy `trajectory`.

### coarse grained model class template.

defined in "coarseGrainedCommon.h"

```cpp
template<typename compType, typename updRateType>
class coarseGrainedModel;
```

This class template defines a basic structure to store information of a coarse grained
biochemical network. It also provides basic reaction rate determination method and 
reactant update method. `compType` is the type for reactant amounts, can be both `int`
or `double`. `updRateType` is the type of `* rate` used when determining the amount of
reactant to be updated. It's suggested to be `double`. But other variable type is still
OK.

#### Public variable members

`type variable_name` | meaning
--- | ---
`int nComp` | number of different reactants in the network
`compType * comp` | datablock to store amount of each reactant of current time
`int nReact` | number of reactions in the network.
`reaction * react` | an array to store all reactions in the network.
`double time`| current time
`double lastSavedTime`| last saved time point
`double stoppingTime`| time when simulation ends
`double saveTimeInterval`| time interval between two saving points

#### Public function members

```cpp
virtual int modelLoading(string & modelName);
```

Load model from a folder named `modelName`. 

---

```cpp
	int rateDetermine(double * rate, compType * comp_alias);
	int rateDetermine(double * rate);
```

`int rateDetermine(double *, compType *)` use the reactant information stored in
`* comp_alias` to calculate the reaction rate. Calculated reaction rate is released
in array pointed by `* rate`. For `int rateDetermine(double *)`, reactant 
information used is the one pointed by `* comp`.

---

```cpp
int reactantUpdate(compType * comp_alias, const updRateType * rate);
int reactantUpdate(const updRateType * rate);
```

`int reactantUpdate(compType *, const updRateType *)` will use `*rate` times the
update information stored in each reaction to update reactant amounts stored in
`* comp_alias`. In `int reactantUpdate(const updRateType *)`, the target reactant
amounts storage would be `* comp`.

---

```cpp
void reset()
```
Reset data pointed by `*comp` with the condition loaded from model file. And setting
`time` and `lastSavedTime` to `0`.

---

```cpp
void assign(const coarseGrainedModel & dummy)
```

Duplicate EVERYTHING from the dummy into the current `coarseGrainedModel`.

#### constructors

```cpp
coarseGrainedModel();
coarseGrainedModel(string & modelName);
coarseGrainedModel(const coarseGrainedModel<compType, updRateType> & dummy);
```

Default constructor won't create a workable instance. User need to load a model from file
with `int modelLoading(string &)` before any other operation. 

`coarseGrainedModel(string &)` loads a model from
file automatically and creates a working instance. 

Copy constructor duplicate the dummy `coarseGrainedModel` into a new instance.

### Coarse Grained Stochastic class

defined in "coarseGrainedOperation.h"

```cpp
class coarseGrainedStochastic: public coarseGrainedModel<int,double> , 
							   public gillespie<coarseGrainedStochastic>
```

This class is based on coarseGrainedModel, using gillespie module for simulation. 

#### public function member

```cpp
void reset();
void assign(const coarseGrainedStochastic & dummy);
```

These two function members are inherited from `coarseGrainedModel`.

---

```cpp
void simulate(string fileName);
```

Simulate with gillespie module. The resulting data will be put into a file named by
`fileName`. This is a temporary design, and will be replaced by other modules when
multi-thread simulation steps in.

#### Constructor

```cpp
coarseGrainedStochastic(string & modelName, double stoppingTime_alias,
		double saveTimeInterval_alias):
	coarseGrainedModel<int,double>(modelName), gillespie<coarseGrainedStochastic>
											   (nReact,
												&coarseGrainedModel::rateDetermine,
												&coarseGrainedModel::reactantUpdate);
coarseGrainedStochastic(coarseGrainedStochastic & dummy):
	coarseGrainedModel<int,double>(dummy),
	gillespie<coarseGrainedStochastic>(nReact,
			&coarseGrainedModel::rateDetermine,
			&coarseGrainedModel::reactantUpdate)
```

Aside of the behavior inherited from the base classes, the constructors shown above will
take `stoppingTime` and `saveTimeInterval` information from either user input or `dummy`.

## Algorithms

The following classes are defined as algorithms used for biochemical network related
simulation.

### Outdated Gillespie Algorithm(gillespieStandAlone)

defined in "gillespieStandAlone.h"
Note: This module used to be called `gillespie`. However, the stand alone design is outdated
and is already replaced by the new `gillespie` class. Old version is thus changed into 
`gillespieStandAlone` for reference purpose.

```cpp
class gillespieStandAlone: public BCNetwork<int, double, int>
```

This class defines both the algorithm and the simulation methods. Gillespie algorithm 
is an exact simulator for biochemical reactions. Continuum 
approximation won't work here. Thus the datatype of reactant amounts are integer. So far
I didn't see the necessity of seperating them. But when there is a need, they might 
be seperated in the future. 


#### Public variable members

`type variable_name` | meaning
--- | ---
`double stoppingTime` | the stopping signal for the simulation
`long long int savePointInterval` | the interval of saving points by count of steps
`long long int nOfNewStep` | No. of steps since the last saving point
`double saveTimeInterval` | the interval of saving points measured by time interval
`double lastSavedTime` | the time point in which last saving action has taken places.
`bool saveMethod`| 0=save every `savePointInterval` steps; 1=save every saveTimeInterval of time;
`bool noSave` | If noSave signal is 1, no data saving is allowed

#### Public function members

```cpp
void simulate();
```

This function will iterate the network with Gillespie algorithm till `BCNetwork::t`
(defined in BCNetwork class template) surpasses `stoppingTime`. During the simulation,
once the saving condition is satisfied, the current time and the current datablock that
`BCNetwork:: *comp` points to will be saved into file specified by 
`BCNetwork::fileOpen(string &)`.

---

```cpp
inline void reset();
```

This reset function reserves the function of `BCNetwork::reset()`. It also restore
`BCNetwork::t`, `lastSavedTime` and `nOfNewStep` to 0.

#### Constructor/Destructor

```cpp
gillespieStandAlone(vector<int> & initComp, vector<double> & inputRate,
        int * inputRateMatrix, int * inputUpdateMatrix,
        double runTime ,bool sMethod=0, double saveInterval=1)
```

It calls `BCNetwork::loadParameter(...)` to load a biochemical network. 
`stoppingTime` will be determined by `runTime`, `saveMethod` will be determined by
`sMethod`. `saveInterval` will determine the value of either `savePointInterval` or
`saveTimeInterval` considering which `saveMethod` user choose.

#### Further note
The random number generator used by this class is rand48 given in gcc library.
Involved function members are

```cpp
inline void reseedRandom(int seed);
inline double popRandom();
```

Should other random number generator be used instead the current one. These two 
member functions can be changed accordingly.

### Gillespie Algorithm class template

defined in "gillespie.h"

```cpp
template<typename modelClassName>
class gillespie
```

This module is built as an stripped version of `gillespieStandAlone` class.
It works the same way as `rungeKutta` class template. A model using this module
must inherit the class template, feed the class template an rate determining 
function member to `int modelClassName::*rateDetermine_alias(double *)`. And a reactant
updating function member to `int (modelClassName::*reactantUpdate_alias)(const double *)`.
The algorithm performs a single iteration of Gillespie Algorithm, deciding the time
step size of the current iteration and which reaction to take place.

Note that

Member function name | Requirements
--- | ---
`int modelClassName::*rateDetermine_alias(double *)`| supplied by model, using model's own reactant information to calculate the reaction rate for each reaction. The rate will be pass to the algorithm as a double type pointer given as the function argument.
`int (modelClassName::*reactantUpdate_alias)(const double *)` | supplied by model, algorithm will pass the number(1 or 0)of each reaction taking places  during the time interval as a constant double type pointer. The function member will use the array multiplying the updating part of the reaction. Resulting value will be added onto the current reactant amount stored in model.

#### Public function member

```cpp
double iterate();
```

The function will do one iteration of algorithm. Reactant amounts will be 
updated by `int reactantUpdate(const double *)`, while the time step is to
be returned as a double type value.

#### Constructors

```cpp
gillespie();
gillespie(int reactionNumber_alias, int (modelClassName::*rateDetermine_alias)(double *),
		int (modelClassName::*reactantUpdate_alias)(const double *))	
```

Default constructor does nothing. Not used in most scenario.

The second constructor is what is mostly used when loading with Gillespie algorithm,
reaction number and the two function member pointer introduced earlier must be provided
when inheriting this algorithm module.


### ODEIVPCommon Class Template

defined in "ODECommon.h"

```cpp
template<typename modelClassType>
class ODEIVPCommon
```

This class is a basic class template for all algorithms created for ODE's
Initial Value Problem simulation. Currently only Runge-Kutta based simulator is
created from it.

`modelClassType` is a typename specified by the models that algorithms are applied
upon. It is used by member function pointers of the algorithms for proper scope 
names.

#### public variable members

`type variable_name` | meaning
--- | ---
`double ht` | current stepsize. For adaptive stepsize algorithms, this value is usually different from `ODENetwork::h0`
`int varNumber` | number of variables/ODEs.
`void (modelClassType::*ODEs)(double *, double *)` | member function pointer. Used to point to ODEs function provided by the specific model. First `double*` point to the returned time derivation, second `double*` points to where the reactants are stored.
`int (*Normalizer)(double *)` | a normalizer, usually need to be specified for each different problem. Thus by definition it's not bound to the member class. Note that this is a function pointer, only static function can be pointed by it.

#### Public function members

N/A

#### Constructor/Destructor

```cpp
ODEIVPCommon();
ODEIVPCommon(int sysSize, double initTimeStep,
		void (modelClassType::*targetODEs)(double *, double *));
ODEIVPCommon(int sysSize, double initTimeStep, 
		void (modelClassType::*targetODEs)(double *, double *),
		int (*targetNormalizer)(double *));
```

Three constructors are defined for this class template. 

Default constructor does nothing.

`ODEIVPCommon(int , double, void (modelClassType::*)(double*, double*));` is
used when no normalizer is required. This constructor will
load the first parameter as `varNumber`, second parameter as the initial value for `ht`,
and third parameter as the address function pointer `void (modelClassType::*ODEs)` pointing
to. Normalizer will be specified to `int blankNormalizer(double *)` defined in "ODECommon.h".
Normalizer can be otherwisely specified to other functions defined by user afterwards.

`ODEIVPCommon(int , double, void (modelClassType::*)(double*, double*), int (*)(double *));`
is used when normalizer is required. First 3 parameters have the same meaning. And the last
parameter is used to specify the pointer pointing to the normalizer function.

### RKmethod Class Template

defined in "rungeKutta.h"

```cpp
template<typename modelClassType>
class RKmethod: public ODEIVPCommon<modelClassType>;
```

This is an algorithm class template using adaptive stepsize Runge-Kutta method for ODEs
simulation of initial value problems. This class will do one thing, and one thing alone.
Given current values of all components, their value of next time step
and the time interval of current iteration will be returned. Users can implement this 
algorithm in their model for simulation and do any analysis they want with the returned
data.

Note that due to the concern of reliable use of this algorithm, this class is not supposed 
to be directly modified by any means, thus most functions are kept private. 

#### Public member

`type variable_name` | meaning
--- | ---
`double hMax` | maximum stepsize of iteration, `hMax=0` deactivates the limit. 


`hMax` is added for situation where a cap of timestep length is necessary. By default 
`hMax=0`, meaning there is no cap value for time step. For any value other than `0`, the
maximum time step will be set to that value. User can change the 
value of `hMax` during or after defining an instance of this class template.


#### Public function members

```cpp
double iterator(double * var);
```

This function takes in the datablock storying current value of the variables pointed by
`* var`. New value after one step of iteration will be directly written into the datablock 
`* var` pointing to. Time interval for this step of iteration will also be returned as a
double value. User can use the time interval to determine the exact time value after this
iteration.

#### Constructor/Destructor

```cpp
RKmethod(maxTimeStep=0);
RKmethod(int sysSize, double initTimeStep,
		void (modelClassType::*targetODEs)(double *, double *), maxTimeStep=0) :
		ODEIVPCommon<modelClassType>::ODEIVPCommon(sysSize, initTimeStep, targetODEs);
RKmethod(int sysSize, double initTimeStep, 
		void (modelClassType::*targetODEs)(double *, double *),  
		void (*targetNormalizer)(double *), maxTimeStep=0): 
		ODEIVPCommon<modelClassType>::ODEIVPCommon(sysSize, initTimeStep,
		targetODEs,targetNormalizer);
```

The use of the three constructors given by this class template mirror the three constructors
of `ODEIVPCommon`. The first one only initiate the private variables to their working states;
Second constructor is used for the case that no normalizer is required; third constructor is
used when a normalizer is required for a successful simulation. If a limit of maximum timestep
is required, users can change the value of `maxTimeStep`, whose value will be passed on to
`hMax` during the class initialization.

## User facing modules for Biochemical Network simulation

Class `gillespie` described in algorithm classes also falls to this catagory. 

### ODESimulate Class

defined in "ODEOperation.h"

```cpp
class ODESimulate : public ODENetwork , public RKmethod<ODESimulate>
```

This class use Runge-Kutta algorithm to do simulation on IVP of ODEs defined by ODENetwork.

#### Public variable members

`type variable_name` | meaning
--- | ---
`double stoppingTime` | the stopping signal for the simulation
`double saveTimeInterval` | the interval of saving points measured by time interval.
`double lastSavedTime` | the time point in which last saving action has taken places.

#### Public function members

```cpp
void simulate(string & identifier);
```

This function will simulate the system and store the result in a file specified by
`& identifier`(it automatically call the function `BCNetwork::fileOpen(string &)`).

---

```cpp
void reset();
```

This function retains the function of `BCNetwork::reset()`. It also restore 
`lastSavedTime` to 0. And restore `ODEIVPCommon::ht` to the value specified by
`ODENetwork::h0`. 

#### Constructor and Destructor

```cpp
ODESimulate(vector<double> & initComp, vector<double> & inputRate,
		int * inputRateMatrix, int * inputUpdateMatrix, 
		double initTimeStep, double runTime, double saveInterval) :
		ODENetwork(initComp, inputRate, inputRateMatrix, inputUpdateMatrix,
		initTimeStep), RKmethod<ODESimulate>(nComp, h0, &ODESimulate::ODETimeDeri);
```

Notice that no default constructor is defined. This constructor will load a typical
biochemical network, and call the constructor of `RKmethod` that doesn't use
normalizer. Two additional information is loaded to specify `stoppingTime` and
`saveTimeInterval`.

```cpp
	stoppingTime=runTime;
	saveTimeInterval=saveInterval;
```

# Model Modification Tools
defined in `modelTool.py` with alias `modelCreator.py`. 

## Class reactionNetwork:

Defines a structure that stores the information of a reaction network, and the 
operations that applies to it.

member name | type | member function
--- | --- | ---
`reactant` | dict | structure to store reactants and the amount of each of them. Indeces are the hash of corresponding reactant's name.
`reaction` | dict | structure to store reactions, their types and coefficient. Indeces are randomly generated number which is ensured to be unique within the network.
`reactant_const` | dict | indeces(hash of reactant name string) of reactants that are kept constant during scaling
`reaction_const` | dict | indeces(randomly generated, unique number) of reactions that are kept constant during scaling
`insert_new_reaction` | func | insert a new reaction to the network, if there are hashes not exists in the network, prompt to add.
`insert_new_reactant` | func | insert a new reactant to the network. If initial value is not provided, prompt to insert.
`list_reaction` | func | list the reactions with their indeces in human understandable format. If `list_target` is not provided, list the whole network.
`list_reactant` | func | list the reactants, their hash, reactant name and initial value.
`load` | func | load a network from earlier saved file. (`_dev` is automatically added, don't add it when specify the modelname)
`save` | func | save the current reaction network into a folder for future loading. (`_dev` will be automatically added after a user specified model name.)
`export` | func | export the current reaction network into a folder that is recognizable by `coarseGrainedxx` simulation modules.
`scaling` | func | scale the network by a certain spacial(`eta_c_alias`) and time(`eta_t_alias`) ratio. Reactions and reactants that are marked constant will not be changed. Reaction rate of non-const-reactions that involves the constant reactants will be adjust accordingly.
`knockout_reactant` | func | knockout the reactants listed in the argument. Involved reactions will be prompted to be removed. Since reactions with undefined reactants are not valid. Not confirming removing the involved reactions will cancel the operation altogether.
`duplicate_reaction` | func | create a duplicate of a currently existing reaction and insert it into the network. Return the index of the newly created reaction.

Note:
1. Since removing entries from an existing dictionary is well defined in python, we do not build additional methods for reaction knocking out. 
