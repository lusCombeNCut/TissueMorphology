/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/ErrorPrinter.h>

#include "CommandLineArguments.hpp"
int main( int argc, char *argv[] ) {
 CommandLineArguments::Instance()->p_argc = &argc;
 CommandLineArguments::Instance()->p_argv = &argv;
 return CxxTest::ErrorPrinter().run();
}
#include "Test3dOrganoidFormation.hpp"

static Test3dOrganoidFormation suite_Test3dOrganoidFormation;

static CxxTest::List Tests_Test3dOrganoidFormation = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_Test3dOrganoidFormation( "Test3dOrganoidFormation.hpp", 76, "Test3dOrganoidFormation", suite_Test3dOrganoidFormation, Tests_Test3dOrganoidFormation );

static class TestDescription_Test3dOrganoidFormation_Test3dSphericalOrganoidFormation : public CxxTest::RealTestDescription {
public:
 TestDescription_Test3dOrganoidFormation_Test3dSphericalOrganoidFormation() : CxxTest::RealTestDescription( Tests_Test3dOrganoidFormation, suiteDescription_Test3dOrganoidFormation, 86, "Test3dSphericalOrganoidFormation" ) {}
 void runTest() { suite_Test3dOrganoidFormation.Test3dSphericalOrganoidFormation(); }
} testDescription_Test3dOrganoidFormation_Test3dSphericalOrganoidFormation;

static class TestDescription_Test3dOrganoidFormation_Test3dBasementMembraneStiffnessEffects : public CxxTest::RealTestDescription {
public:
 TestDescription_Test3dOrganoidFormation_Test3dBasementMembraneStiffnessEffects() : CxxTest::RealTestDescription( Tests_Test3dOrganoidFormation, suiteDescription_Test3dOrganoidFormation, 186, "Test3dBasementMembraneStiffnessEffects" ) {}
 void runTest() { suite_Test3dOrganoidFormation.Test3dBasementMembraneStiffnessEffects(); }
} testDescription_Test3dOrganoidFormation_Test3dBasementMembraneStiffnessEffects;

static class TestDescription_Test3dOrganoidFormation_Test3dLongTermOrganoidDevelopment : public CxxTest::RealTestDescription {
public:
 TestDescription_Test3dOrganoidFormation_Test3dLongTermOrganoidDevelopment() : CxxTest::RealTestDescription( Tests_Test3dOrganoidFormation, suiteDescription_Test3dOrganoidFormation, 253, "Test3dLongTermOrganoidDevelopment" ) {}
 void runTest() { suite_Test3dOrganoidFormation.Test3dLongTermOrganoidDevelopment(); }
} testDescription_Test3dOrganoidFormation_Test3dLongTermOrganoidDevelopment;

#include <cxxtest/Root.cpp>
