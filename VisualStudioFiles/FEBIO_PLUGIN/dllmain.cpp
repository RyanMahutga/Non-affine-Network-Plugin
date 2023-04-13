// This file defines the required functions for registering the plugin classes
// with the FEBio framework. 
#include <FECore/FECoreKernel.h>
#include "FENetwork.h"

//-----------------------------------------------------------------------------
// This required function returns the version of the FEBio SDK that is being
// used by the plugin. This version number will be checked by FEBio to make
// sure that the FEBio executable is compatible with this SDK. Usually this
// function just returns the predefined macro FE_SDK_VERSION.
FECORE_EXPORT unsigned int GetSDKVersion()
{
	return FE_SDK_VERSION;
}

//-----------------------------------------------------------------------------
// If FEBio is successful in loading the plugin, this is the first function
// that will be called and can be used to initialize any resources that are
// needed by the plugin. FEBio passes a reference to the FECoreKernel as a parameter,
// which must be used to initialize the FECoreKernel.
FECORE_EXPORT void PluginInitialize(FECoreKernel& febio)
{
	// Set the kernel's instance to the same instance as used by febio. 
	// This is to ensure that the plugin features are registered in FEBio's kernel. 
	FECoreKernel::SetInstance(&febio);

	// This macro registers the new feature and assign a string to it that
	// can be used to reference this class in the FEBio input file.
	REGISTER_FECORE_CLASS(FENetwork, "PeriodicNetwork");
}

//-----------------------------------------------------------------------------
// This function is called when FEBio exits and gives the plugin a chance to close
// any resource that have been allocated. 
// This is an optional function.
FECORE_EXPORT void PluginCleanup()
{

}
