/*
 * Client.cpp
 *
 *  Created on: Dec 21, 2013
 *      Author: Alberto Paganelli
 */

#include <Ice/Ice.h>
#include <Freeze/Freeze.h>

#include "GenmapParams.h"

#include "Astro.h"
#include <DBAstro.h>

using namespace std;
using namespace Astro;

//#define SIMPLE_KEY
//#define COMPOSITE_KEY

class AstroClient : public Ice::Application {
public:

	AstroClient();
	virtual int run(int, char*[]);

private:

	void menu();
};

int main(int argc, char **argv) {
	AstroClient app;
	return app.main(argc, argv, "config.client");
}

AstroClient::AstroClient() : Ice::Application(Ice::NoSignalHandling) {}

int AstroClient::run(int argc, char* argv[]){

//	if (argc > 1 ) {
//		cerr << appName() << ": too many arguments" << endl;
//		return EXIT_FAILURE;
//	}

	AGILECtsMapGenParams params;
//	if (!params.Load(argc, argv)) {
//		return -1;
//	}

	AgileCtsMapGenPrx twoway = AgileCtsMapGenPrx::checkedCast(communicator()->propertyToProxy("Astro.Proxy")->ice_twoway()->ice_timeout(-1)->ice_secure(false));

	if (!twoway)
	{
		cerr << argv[0] << ": invalid proxy" << endl;
		return EXIT_FAILURE;
	}
	AgileCtsMapGenPrx oneway = twoway->ice_oneway();
	AgileCtsMapGenPrx batchOneway = twoway->ice_batchOneway();
	AgileCtsMapGenPrx datagram = twoway->ice_datagram();
	AgileCtsMapGenPrx batchDatagram = twoway->ice_batchDatagram();

//	Ice::InitializationData initData;
//	initData.properties = Ice::createProperties();
//	initData.properties->load("config");
//
//	// Initialize the Communicator.
//	Ice::CommunicatorPtr communicator = Ice::initialize(initData);
//
//	// Create a Freeze database connection.
//	Freeze::ConnectionPtr connection = Freeze::createConnection(communicator, "db");
//
//
//	//The map
//	DBAstro DBEvt(connection,"AgileEvtMap");
//	//The iterator
//	DBAstro::iterator it;

#ifdef SIMPLE_KEY
	//The evt vector
	Astro::agileEvt agileEvt;

//	vector<double> ra, dec;
	Astro::Ra ra;
	Astro::Dec dec;

	for(it=DBEvt.begin(); it != DBEvt.end(); ++it){
		agileEvt = it->second;
		ra.push_back(agileEvt[6]);
		dec.push_back(agileEvt[5]);
	}
#endif

#ifdef COMPOSITE_KEY

	Astro::SeqEvtKey evtKeys;
	Astro::AgileEvtKey key;

	for(it=DBEvt.begin(); it != DBEvt.end(); ++it){
		key = it->first;
		evtKeys.push_back(key);
	}

#endif

	menu();

	char c;

	bool secure = false;
	int timeout = -1;

	do {
		try {
			cout << "==> ";
			cin >> c;
			if(c == 't')
			{

				Astro::Matrix retv = twoway->calculateMapKey(params);

#ifdef SIMPLE_KEY
				Astro::Matrix retv = twoway->calculateMapVector(ra, dec);
#endif

#ifdef COMPOSITE_KEY
				Astro::Matrix retv = twoway->calculateMapKey(evtKeys);
#endif

				cout << "Received back the matrix" << endl;
				//Demo::FloatSeq retv = twoway->update(1, 2.0, "ciao", v);
//				cout << "Received back the vector [ ";
//				for(unsigned int i=0; i<retv.size(); i++)
//					cout << retv[i] << " ";
//				cout << "]" << endl;
			}
			else if(c == 'o')
			{
				oneway->calculateMapKey(params);
			}
			else if(c == 'O')
			{
				batchOneway->calculateMapKey(params);
			}
			else if(c == 'd')
			{
				if(secure)
				{
					cout << "secure datagrams are not supported" << endl;
				}
				else
				{
					datagram->calculateMapKey(params);
				}
			}
			else if(c == 'D')
			{
				if(secure)
				{
					cout << "secure datagrams are not supported" << endl;
				}
				else
				{
					batchDatagram->calculateMapKey(params);
				}
			}
			else if(c == 'f')
			{
				Ice::Application::communicator()->flushBatchRequests();
			}
			else if(c == 'T')
			{
				if(timeout == -1)
				{
					timeout = 2000;
				}
				else
				{
					timeout = -1;
				}

				twoway = twoway->ice_timeout(timeout);
				oneway = oneway->ice_timeout(timeout);
				batchOneway = batchOneway->ice_timeout(timeout);

				if(timeout == -1)
				{
					cout << "timeout is now switched off" << endl;
				}
				else
				{
					cout << "timeout is now set to 2000ms" << endl;
				}
			}
			else if(c == 'S')
			{
				secure = !secure;

				twoway = twoway->ice_secure(secure);
				oneway = oneway->ice_secure(secure);
				batchOneway = batchOneway->ice_secure(secure);
				datagram = datagram->ice_secure(secure);
				batchDatagram = batchDatagram->ice_secure(secure);

				if(secure)
				{
					cout << "secure mode is now on" << endl;
				}
				else
				{
					cout << "secure mode is now off" << endl;
				}
			}
			else if(c == 's')
			{
				twoway->shutdown();
			}
			else if(c == 'x')
			{
				// Nothing to do
			}
			else if(c == '?')
			{
				menu();
			}
			else
			{
				cout << "unknown command `" << c << "'" << endl;
				menu();
			}
		} catch (const Ice::Exception& e) {
			cerr << e << endl;
		}
	} while (cin.good() && c != 'x');


}

void
AstroClient::menu()
{
    cout <<
        "usage:\n"
        "t: send a vector as twoway\n"
        "o: send a vector one way oneway\n"
        "O: send greeting as batch oneway\n"
        "d: send greeting as datagram\n"
        "D: send greeting as batch datagram\n"
        "f: flush all batch requests\n"
        "T: set a timeout\n"
        "S: switch secure mode on/off\n"
        "s: shutdown server\n"
        "x: exit\n"
        "?: help\n";
}
