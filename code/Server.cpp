/*
 * Server.cpp
 *
 *  Created on: Dec 22, 2013
 *      Author: Alberto Paganelli
 */

#include <Ice/Ice.h>
#include "AstroI.h"

using namespace std;

class AstroServer : public Ice::Application
{

	virtual int run(int, char*[]);

};


int main(int argc, char **argv) {
	AstroServer app;
	return app.main(argc, argv, "config.server");
}

int AstroServer::run(int argc, char*[]){
	if(argc > 1)
	{
		cerr << appName() << ": too many arguments" << endl;
		return EXIT_FAILURE;
	}

	Ice::ObjectAdapterPtr adapter = communicator()->createObjectAdapter("Astro");
	Astro::AgileCtsMapGenPtr astro = new Astro::AgileCtsMapGenI;
	adapter->add(astro, communicator()->stringToIdentity("astro"));
	adapter->activate();
	communicator()->waitForShutdown();
	return EXIT_SUCCESS;
}
