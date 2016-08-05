#include <string>

#include "Generator.h"


namespace fPotencia {
    Generator::Generator():
            name_(""),
            busIndex_(-1),
            voltage_(0.0),
            voltageType_(Generator::pu),
            power_(0.0, 0.0),
            minQ_(0.0),
            maxQ_(0.0)
    {
    }
}
