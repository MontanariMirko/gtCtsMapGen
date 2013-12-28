/***************************************************************************
    begin                : Mar 16 2005
    copyright            : (C) 2005-2013 A. Chen, A. Trois, T. Contessi
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software for non commercial purpose              *
 *   and for public research institutes; you can redistribute it and/or    *
 *   modify it under the terms of the GNU General Public License.          *
 *   For commercial purpose see appropriate license terms                  *
 *                                                                         *
 ***************************************************************************/

#include <sstream>
#include <iostream>

#include <cstring>

#include <GenmapParams.h>
#include "MathUtils.h"


using namespace std;

//ignorare
string String(const Interval& intv)
{
    stringstream str(ios_base::out);
    str.precision(6);
    str << fixed << intv.Start() << ".." << intv.Stop();
    return str.str();
}

//ignorare
string String(const Intervals& intvs)
{
    stringstream str(ios_base::out);
    const char* sep = "";
    for (int i=0; i<intvs.Count(); ++i) {
        str << sep << String(intvs[i]);
        sep = ", ";
    }
    return str.str();
}

//questa è la funziona che contiene la query sui dati. Importante.
//E' da rifare accedendo al BDB
static string LogEvtString(const Intervals& intvs, double emin, double emax, double albrad, double fovradmax, double fovradmin, int phasecode, int filtercode)
{
	//se Count() == 0 non si calcola la mappa
    if (intvs.Count()<=0)
        return string("");
    stringstream str(ios_base::out);
    str.precision(6);
    //se Count() > 0 si calcola la mappa.
    //se Count() > 1 allora abbiamo una lista di tempi Start() e Stop().
    //per ora noi possiamo ignorare questo e considerare per semplicità solo un tstart ed un tstop.
    
    if (intvs.Count()==1)
        str << "TIME >= " << fixed << intvs[0].Start() << " && TIME <= " << intvs[0].Stop();
    else {
        str << "( ";
        const char* sep = "";
        for (int i=0; i<intvs.Count(); ++i) {
            str << sep << " ( TIME >= " << fixed << intvs[i].Start() << " && TIME <= " << intvs[i].Stop() << " )";
            sep = " ||";
        }
        str << " )";
    }

    str << " && ENERGY >= " << emin;
    str << " && ENERGY <= " << emax;
    str << " && PH_EARTH > " << albrad;
    str << " && THETA < " << fovradmax;
    str << " && THETA >= " << fovradmin;

    if ((phasecode & 1) == 1)
        str << " && PHASE .NE. 0";
    if ((phasecode & 2) == 2)
        str << " && PHASE .NE. 1";
    if ((phasecode & 4) == 4)
        str << " && PHASE .NE. 2";
    if ((phasecode & 8) == 8)
        str << " && PHASE .NE. 3";
    if ((phasecode & 16) == 16)
        str << " && PHASE .NE. 4";

    if ((filtercode & 1) == 1)
        str << " && EVSTATUS .NE. 'L'";
    if ((filtercode & 2) == 2)
        str << " && EVSTATUS .NE. 'G'";
    if ((filtercode & 4) == 4)
        str << " && EVSTATUS .NE. 'S'";

    return str.str();
}

//AB: questa funzione non serve. Noi facciamo accesso diretto al DB. Da ignorare
static int addfile(fitsfile* iFile, const CtsGenParams& params)
{

    std::cout << "params.evtfile: " << params.evtfile << std::endl;

    char buffer[1024];
    int status = 0;
    FILE *fp = fopen(params.evtfile, "r");
    if (!fp) {
        cerr << "Error opening file " << params.evtfile << endl;
        return 104;
    }

    cout << "Intervals: " << String(params.intervals) << endl;

    int hdutype = 0;
    bool firstLog = true;
    while (fgets(buffer , 40960, fp)) {
        char name[FLEN_FILENAME];
        double t1 = 0, t2 = 0;
        sscanf(buffer, "%s %lf %lf", name, &t1, &t2);
        Interval logIntv(t1, t2);
        Intervals logSel = Intersection(params.intervals, logIntv);

        if (logSel.Count())
            cout << "Selecting from: " << name << endl;

        for (int i=0; i<logSel.Count(); ++i) {
            Intervals logSelInt;
            logSelInt.Add(logSel[i]);
            string exprStr = LogEvtString(logSelInt, params.emin, params.emax, params.albrad, params.fovradmax, params.fovradmin, params.phasecode, params.filtercode);
            char expr[1024];
            strcpy(expr, exprStr.c_str());
            fitsfile* tempFits;
            if (firstLog) {
                if (fits_open_file(&tempFits, name, READONLY, &status) != 0 ) {
                    cerr << "FITS Error " << status << " opening file " << name << endl;
                    return status;
                }
                fits_copy_file(tempFits, iFile, 1, 1, 1, &status);
                fits_movabs_hdu(iFile, 2, &hdutype, &status);
                fits_select_rows(iFile, iFile, expr, &status);
            }
            else {
                if (fits_open_file(&tempFits, name, READONLY, &status) != 0 ) {
                    cerr << "FITS Error " << status << " opening file " << name << endl;
                    return status;
                }
                fits_movabs_hdu(tempFits, 2, &hdutype, &status);
                fits_select_rows(tempFits, iFile, expr, &status);
            }
            fits_close_file(tempFits, &status);
            firstLog = false;
        }
    }
    fclose(fp);
    if (firstLog)
        return 1005;
    return status;
}


//questa è la funzione che calcola le mappe. 
int calculateMaps(CtsGenParams & params)
{
    std::cout << "params.evtfile: " << params.evtfile << std::endl;

    int numcol = 0;
    long nrows = 0;
    int status = 0;
    double l = 0, b = 0;
    double x = 0, y = 0;
    int i = 0, ii = 0;
    double the = 0;
    long mxdim= params.mxdim; // dimension (in pixels) of the map
    unsigned short A[mxdim][mxdim]; //counts map

    int selectedEvents = 0;

    for (i = 0; i < mxdim; i++)
    {   for (ii = 0; ii < mxdim; ii++)
        {
            A[i][ii] = 0;
        }
    }

    double baa = params.ba * DEG2RAD;
    double laa = params.la * DEG2RAD;

    int bitpix   =  USHORT_IMG; /* 16-bit unsigned short pixel values       */
    long naxis    =   2;  /* 2-dimensional image                            */
    long naxes[2] = { mxdim, mxdim };   /* image is 300 pixels wide by 200 rows */

	//accesso ai file: da sostituire qui con l'accesso al DBD
    fitsfile * evtFits;
    char tempname[FLEN_FILENAME];
    strcpy(tempname, tmpnam(NULL));
    if ( fits_create_file(&evtFits, tempname, &status) != 0 ) {
        printf("Errore in apertura file %s\n",tempname);
        return status;
    }

    const double obtlimit = 104407200.;
    if (params.tmin<obtlimit)
        status = 1005;
    else
        status = addfile(evtFits, params);
    std::cout << "AG_ctsmapgen0...................................addfile exiting STATUS : "<< status<< std::endl << std::endl ;

	//creazione della mappa. Questo lo lasciamo così
    fitsfile * mapFits;
    if ( fits_create_file(&mapFits, params.outfile, &status) != 0 ) {
        printf("Errore in apertura file '%s'\n",params.outfile);
        return status;
    }

	//tutte le parti che iniziano per fits_ vanno tolte.
	//si fa qui la query a BDB che restituisce un vector<double> con il dato RA e un 
	//vector<double> con il dato DEC
    fits_movabs_hdu(evtFits, 2, NULL, &status);
    fits_get_num_rows(evtFits, &nrows, &status);
    cout << nrows << endl;
    double ra, dec;
    switch (params.projection) {
    case ARC:
        for (long k = 0; k<nrows; ++k) {
        	//sostituire tutto questo con l'accesso al singolo elemento k del vector<double>
            fits_get_colnum(evtFits, 1, "RA", &numcol, &status);
            fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &ra, NULL, &status);
            fits_get_colnum(evtFits, 1, "DEC", &numcol, &status);
            fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &dec, NULL, &status);
            Euler(ra, dec, &l, &b, 1);
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

            i = (int)floor(((-x+(params.mdim/2.))/params.mres));
            ii = (int)floor(((y+(params.mdim/2.))/params.mres));

            if (params.inmap(i,ii)) {
                A[ii][i]+=1;
                ++selectedEvents;
            }
        }
        break;
    case AIT:
        for (long k = 0; k<nrows; ++k) {
            fits_get_colnum(evtFits, 1, "RA", &numcol, &status);
            fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &ra, NULL, &status);
            fits_get_colnum(evtFits, 1, "DEC", &numcol, &status);
            fits_read_col(evtFits, TDOUBLE, numcol, k+1, 1, 1, NULL, &dec, NULL, &status);
            Euler(ra, dec, &l, &b, 1);
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

            i=(int)floor(((x+(params.mdim/2.))/params.mres));
            ii=(int)floor(((y+(params.mdim/2.))/params.mres));

            if (params.inmap(i,ii)) {
                A[ii][i]+=1;
                ++selectedEvents;
            }
        }
        break;
    }

    long nelement =  naxes[0] * naxes[1];
    std::cout<< "creating Counts Map...................................." << std::endl;
    fits_create_img(mapFits, bitpix, naxis, naxes, &status);
    std::cout<< "writinig Counts Map with " << selectedEvents << " events" << std::endl;
    fits_write_img(mapFits, bitpix, 1, nelement, A, &status);
    std::cout<< "writing header........................................" << std::endl<< std::endl;

    params.write_fits_header(mapFits, params.projection, status);

    fits_delete_file(evtFits, &status);
    fits_close_file(mapFits, &status);
    return status;
}




int main(int argC,char* argV[])
{
    std::cout << " " << std::endl;
    std::cout << " " << std::endl;
    std::cout << "#################################################################" << std::endl;
    std::cout << "####### AG_ctsmapgen0  V 1.0 - 17/11/11 - A.C., A.T., T.C. ######" << std::endl;
    std::cout << "#################################################################" << std::endl;
    std::cout << "#################################################################" << std::endl;
    std::cout << " " << std::endl;

	//questa la possiamo tenere così per il caricamento dei parametri
    CtsGenParams params;
    if (!params.Load(argC, argV))
        return -1;

    std::cout << "params.evtfile: " << params.evtfile << std::endl;
    std::cout << "AG_ctsmapgen...............................starting"<< std::endl;
    int status = calculateMaps(params);
    std::cout << "AG_ctsmapgen............................... exiting"<< std::endl;


    if (status) {
        if (status != 105) {
            remove(params.outfile);
        }
        printf("AG_ctsmapgen..................... exiting AG_ctsmapgen ERROR:");
        fits_report_error(stdout, status);
        return status;
    }
    else {
        printf("\n\n\n###################################################################\n");
        printf("#########  Task AG_ctsmapgen0.......... exiting #################\n");
        printf("#################################################################\n\n\n");
    }

    return status;
}

