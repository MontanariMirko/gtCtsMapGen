#ifndef __AstroI_h__
#define __AstroI_h__

#include <Astro.h>
#include <vector>

namespace Astro
{

class AgileCtsMapGenI : virtual public AgileCtsMapGen
{
public:

    virtual ::Astro::Matrix calculateMapKey(const ::Astro::AGILECtsMapGenParams&,
                                            const Ice::Current&);

    virtual void shutdown(const Ice::Current&);

private:
    bool LogEvtString(AGILECtsMapGenParams params, std::vector<float> &ra, std::vector<float> &dec);
};

}

#endif
