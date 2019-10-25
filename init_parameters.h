#ifndef INIT_PARAMETERS
#define INIT_PARAMETERS
#include "const_parameters.h"

extern std::vector< pair<int, char*> > imageTypeStringMapping;

const char* depthCode2String(int depth);

void init();

bool readParam(string&flag, string&param, bool&bFloatValid, float&value, char*arg);

bool readAndConvertImages(cv::Mat imgGray[2], cv::Mat img[2]/*char** argv*/);

void help();

bool ensureGrayImageWithDepth(cv::Mat& imgGray, int target_depth, cv::Mat& img,
	                          char*image_name = NULL, bool bPrintFormatInfo = true);

//////////////////////////////////////////////////
bool deal_init_paramers(const int& argc);

#endif
