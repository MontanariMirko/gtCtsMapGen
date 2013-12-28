module Astro
{

sequence<float> agileEvt;
const int agileEvtSize = 6;

sequence<double> agileLog;
const int agileLogSize = 13;

struct agileEvtValue
{
float theta;
float phi;
float ra;
float dec;
int energy;
float phearth;
byte evstatus;
};

sequence<agileEvtValue> agileEvtValueSeq;

struct agileEvtKey
{
double time;
float ra;
float dec;
float energy;
float theta;
float evstatus;
};

struct agileLogKey
{
double time;
double livetime;
double logStatus;
double mode;
double phase;
};

sequence<int> SimpleSeq;
sequence<SimpleSeq> Matrix;

sequence<float> Ra;
sequence<float> Dec;

sequence<agileEvtKey> SeqEvtKey;

struct IntervalTime {
double tstart;
double tstop;
};

sequence<IntervalTime> Intervals;

struct AGILECtsMapGenParams {
double mdim; 
double mres; 
long mxdim; 
double la; 
double ba; 
double emin; //ENERGY
double emax; //ENERGY
double fovradmin; //THETA
double fovradmax; //THETA
double albrad; 
double lonpole; 
int phasecode; //PHASE
int filtercode; //EV_STATUS
string projection;
Intervals paramIntervals; //non ricordo se è case sensitive
};

interface AgileCtsMapGen
{

Matrix calculateMapKey(AGILECtsMapGenParams param);
void shutdown();

};

};
