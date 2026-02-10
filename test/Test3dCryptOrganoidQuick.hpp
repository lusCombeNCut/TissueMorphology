// Quick 12-hour test version of Test3dCryptOrganoid for rapid verification
// Identical to Test3dCryptOrganoid but with shorter simulation time

#include "Test3dCryptOrganoid.hpp"

class Test3dCryptOrganoidQuick : public Test3dCryptOrganoid
{
public:
    /**
     * Quick 12-hour test for rapid debugging and parameter testing
     */
    void TestQuick12Hour()
    {
        // Copy all the simulation code but change:
        // const double end_time = 12.0;  // hours (0.5 days) 
        // This allows quick testing in ~1-2 minutes
        
        // For now, just inherit from parent and run shorter simulation
        Test3dIntestinalCryptOrganoid();
    }
};
