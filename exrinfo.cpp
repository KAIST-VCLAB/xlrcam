#include <iostream>

#include "ImfRgbaFile.h"
#include "ImfStringAttribute.h"
#include "ImfMatrixAttribute.h"
#include "ImfChannelList.h"
#include "ImfPixelType.h"
#include "Iex.h"

#include "mex.h" 

using namespace Imf;
using namespace Imath;
using namespace Iex;

using std::cout;
using std::endl;
using std::flush;


/*
 * Check inputs:
 * only one input argument that is a string (row vector of chars)
 * one or two output arguments
 * 
 * These checks were copied from the MATLAB example file revord.c
 */
void checkInputs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	if (nrhs != 1)
		mexErrMsgTxt("Incorrect number of input arguments");

	if (nlhs > 1)
		mexErrMsgTxt("Incorrect number of output arguments");

	if (mxIsChar(prhs[0]) != 1)
		mexErrMsgTxt("Input must be a string");

	if (mxGetM(prhs[0]) != 1)
		mexErrMsgTxt("Input must be a row vector.");

	return;
}

/*
 * Read the header of an EXR file.
 * Code follows examples from ReadingAndWritingImageFiles.pdf, found
 * here:
 * http://www.openexr.com/ReadingAndWritingImageFiles.pdf
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { 

	checkInputs(nlhs, plhs, nrhs, prhs);
	char *filename = mxArrayToString(prhs[0]);

	try {
		RgbaInputFile file(filename);

		const ChannelList &ch = file.header().channels();
		int ix = 1;
		for (ChannelList::ConstIterator i = ch.begin(); i != ch.end(); ++i) {
			const Channel &channel = i.channel(); 
			const char* n = i.name(); // by Min
			PixelType type = channel.type;
			const char* t = (type == UINT) ? "uint" : 
					((type == HALF) ? "half" : "float");

			cout << "channel " << ix++ << ": ";
			cout << t << ": " << n << endl;
		}
		const StringAttribute *comments =
		file.header().findTypedAttribute <StringAttribute> ("comments");

		const M44fAttribute *cameraTransform = 
		file.header().findTypedAttribute <M44fAttribute> ("cameraTransform");

		if (comments)
			cout << "comments\n   " << comments->value() << endl;

		if (cameraTransform)
			cout << "cameraTransform\n" << cameraTransform->value() << flush;

	} catch (const std::exception &exc) {
		mexErrMsgTxt(exc.what());
	}

	// Free the memory for the string
	mxFree(filename);

	return;
} 



