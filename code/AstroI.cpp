
#include <AstroI.h>
#include <sstream>
#include <iostream>
#include <vector>

#include <cstring>

#include <GenmapParams.h>
#include "MathUtils.h"

#include "DBAstro.h"
#include "Astro.h"
#include <Freeze/Freeze.h>

using namespace std;
using namespace Astro;

#define SIMPLE_KEY
//#define COMPOSITE_KEY

::Astro::Matrix
Astro::AgileCtsMapGenI::calculateMapKey(const ::Astro::AGILECtsMapGenParams& param,
                                        const Ice::Current& current)
{
	int numcol = 0;
	long nrows = 0;
	int status = 0;
	double l = 0, b = 0;
	double x = 0, y = 0;
	int i = 0, ii = 0;
	double the = 0;
	long mxdim= param.mxdim; // dimension (in pixels) of the map
	//unsigned short A[mxdim][mxdim]; //counts map
	Astro::Matrix A(mxdim, SimpleSeq(mxdim));
	int selectedEvents = 0;

	for (i = 0; i < mxdim; i++)
	{   for (ii = 0; ii < mxdim; ii++)
		{
			A[i][ii] = 0;
		}
	}

	double baa = param.ba * DEG2RAD;
	double laa = param.la * DEG2RAD;

	int bitpix   =  USHORT_IMG; /* 16-bit unsigned short pixel values       */
	long naxis    =   2;  /* 2-dimensional image                            */
	long naxes[2] = { mxdim, mxdim };   /* image is 300 pixels wide by 200 rows */

	const double obtlimit = 104407200.;

//	if (param.tmin < obtlimit)
//		status = 1005;
//	else
		//status = addfile(evtFits, params);
	std::cout << "AG_ctsmapgen0...................................addfile exiting STATUS : "<< status<< std::endl << std::endl ;

	vector<float> raVector, decVector;
	LogEvtString(param, raVector, decVector);

	nrows = raVector.size();

	cout << nrows << endl;
	//double ra, dec;
	//switch (param.projection) {
	//case ARC:
	if(param.projection == "ARC"){
		for (long k = 0; k<nrows; ++k) {
			//sostituire tutto questo con l'accesso al singolo elemento k del vector<double>
			//fits_get_colnum(evtFits, 1, "RA", &numcol, &status);
			//fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &ra, NULL, &status);
			//fits_get_colnum(evtFits, 1, "DEC", &numcol, &status);
			//fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &dec, NULL, &status);

			Euler(raVector[k], decVector[k], &l, &b, 1);
				//Euler(ra, dec, &l, &b, 1);

			l*=DEG2RAD;
			b*=DEG2RAD;
			the = sin(b)*sin(baa)+cos(b)*cos(baa)*cos(l-laa);
			if (the < -1.0)
				the = M_PI;
			else if (the > 1.0)
				the = 0.0;
			else
				the = acos(the);
			x = RAD2DEG/Sinaa(the) * cos(b)*sin(l-laa);
			y = RAD2DEG/Sinaa(the) * (sin(b)*cos(baa) - cos(b)*sin(baa)*cos(l-laa));

			i = (int)floor(((-x+(param.mdim/2.))/param.mres));
			ii = (int)floor(((y+(param.mdim/2.))/param.mres));

			//if (param.inmap(i,ii)) {
			if(true){
				A[ii][i]+=1;
				++selectedEvents;
			}
		}
	}
	//case AIT:
	if(param.projection == "AIT"){
		for (long k = 0; k<nrows; ++k) {
			//fits_get_colnum(evtFits, 1, "RA", &numcol, &status);
			//fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &ra, NULL, &status);
			//fits_get_colnum(evtFits, 1, "DEC", &numcol, &status);
			//fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &dec, NULL, &status);
			//Euler(ra, dec, &l, &b, 1);

		Euler(raVector[k], decVector[k], &l, &b, 1);

			l*=DEG2RAD;
			b*=DEG2RAD;
			the = sin(b)*sin(baa)+cos(b)*cos(baa)*cos(l-laa);
			if (the < -1.0)
				the = M_PI;
			else if (the > 1.0)
				the = 0.0;
			else
				the = acos(the);
			l=l-laa;

			if ( l < M_PI  )
				l=-l;
			else
				l=2*M_PI -l;

			x=RAD2DEG*(sqrt(2.0)*2.0*cos(b)*sin(l/2.0))/sqrt(1.0 + cos(b)*cos(l/2.0) ) ;
			y=RAD2DEG*(sqrt(2.0)*sin(b))/sqrt(1.0 + cos(b)*cos(l/2.0) );

			i=(int)floor(((x+(param.mdim/2.))/param.mres));
			ii=(int)floor(((y+(param.mdim/2.))/param.mres));

			//if (param.inmap(i,ii)) {
			if(true){
				A[ii][i]+=1;
				++selectedEvents;
			}
		}
		//break;
	}

	return A;
}

void
Astro::AgileCtsMapGenI::shutdown(const Ice::Current& current)
{
	cout << "Shutting down..." << endl;
	current.adapter->getCommunicator()->shutdown();
}

bool LogEvtString(AGILECtsMapGenParams params, vector<float> &ra, vector<float> &dec)
{
	//se Count() == 0 non si calcola la mappa
    if (params.paramIntervals.size() <= 0)
        return false;

    Ice::InitializationData initData;
		initData.properties = Ice::createProperties();
		initData.properties->load("config");

		// Initialize the Communicator.
		Ice::CommunicatorPtr communicator = Ice::initialize(initData);

		// Create a Freeze database connection.
		Freeze::ConnectionPtr connection = Freeze::createConnection(communicator, "db");

		//The map
		DBAgileEvt DBEvt(connection,"AgileEvtMap");
		//The iterator
		DBAgileEvt::iterator it;

	#ifdef SIMPLE_KEY
		//The evt vector
		Astro::agileEvt agileEvt;

		for(it=DBEvt.begin(); it != DBEvt.end(); ++it){
			double time = it->first;
			if (params.paramIntervals.size() == 1) {
				if (time >= params.paramIntervals[0].tstart && time <= params.paramIntervals[0].tstop) {
					agileEvt = it->second;
					if ((agileEvt[4] >= params.emin && agileEvt[4] <= params.emax) &&
							(agileEvt[3] > params.albrad) &&
							(agileEvt[2] < params.fovradmax && agileEvt[2] >= params.fovradmin) &&
							(agileEvt[1] == params.phasecode) &&
							(agileEvt[0] == params.filtercode)) {

						//Condizioni soddisfatte
						ra.push_back(agileEvt[6]);
						dec.push_back(agileEvt[5]);

					}
				}
			} else {
				bool inTime = false;
				for (int i=0; i<params.paramIntervals.size(); ++i) {
					if (time >= params.paramIntervals[i].tstart && time <= params.paramIntervals[i].tstop) {
						inTime = true;
						break;
					}
				}
				if (inTime) {
					agileEvt = it->second;
					if ((agileEvt[4] >= params.emin && agileEvt[4] <= params.emax) &&
							(agileEvt[3] > params.albrad) &&
							(agileEvt[2] < params.fovradmax && agileEvt[2] >= params.fovradmin) &&
							(agileEvt[1] == params.phasecode) &&
							(agileEvt[0] == params.filtercode)) {

						//Condizioni soddisfatte
						ra.push_back(agileEvt[6]);
						dec.push_back(agileEvt[5]);

					}
				}
			}

		}

	#endif

	#ifdef COMPOSITE_KEY

		Astro::agileEvtKey key;
		Astro::agileEvt evt;

		for(it=DBEvt.begin(); it != DBEvt.end(); ++it){
			key = it->first;
			if (params.paramIntervals.size() == 1) {
				if (key.time >= params->paramIntervals[0].tstart && key.time <= params->paramIntervals[0].tstop) {
					evt = it->second;
					if ((key.energy >= params->emin && key.energy <= params->emax) &&
							(evt[3] > params->albrad) &&
							(key.theta < params->fovradmax && agileEvt[2] >= key.theta->fovradmin) &&
							(evt[1] == params->phasecode) &&
							(evt[0] == params->filtercode)) {

						//Condizioni soddisfatte
						ra.push_back(key.ra);
						dec.push_back(key.dec);

					}
				}
			} else {
				bool inTime = false;
				for (int i=0; i<params->paramIntervals.size(); ++i) {
					if (key.time >= params->paramIntervals[i].tstart && key.time <= params->paramIntervals[i].tstop) {
						inTime = true;
						break;
					}
				}
				if (inTime) {
					evt = it->second;
					if ((key.energy >= params->emin && key.energy <= params->emax) &&
							(evt[3] > params->albrad) &&
							(key.theta < params->fovradmax && agileEvt[2] >= key.theta->fovradmin) &&
							(evt[1] == params->phasecode) &&
							(evt[0] == params->filtercode)) {

						//Condizioni soddisfatte
						ra.push_back(key.ra);
						dec.push_back(key.dec);

					}
				}
			}
		}

	#endif

		return true;
}
