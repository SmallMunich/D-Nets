#ifndef IMAGE_PYRAMID
#define IMAGE_PYRAMID
#include "const_parameters.h"

extern int nLayers;

class ImagePyramid
{
public:
	enum ScaleMode { UNDEFINED, DIF_UNIFORM, LOG_UNIFORM };

	// UNIFORM: the difference between sucessive scales is constant
	// LOG_UNIFORM: the fraction between successive scales is constant.
	ImagePyramid();

	virtual ~ImagePyramid();

	/*
		关于函数缺省值.h与.cpp情况 一般函数默认的缺省值要写在声明的地方
		而不是写在定义的地方，即写在.h里面不写在.cpp函数定义里
	*/
	bool create(cv::Mat *pXI, int nLayers = 4, float lastRelativeScale = 0.125f,
		        bool bLogarithmicallyEquallySpaced = false,
		        bool bAlwaysScaleFromFirstImage = false);

	bool getBestLevelIndex(int& k, double relativeDownScaleFactor);

	void updateDependencies();

	// Member variables
	ScaleMode scaleMode;
	double constScaleDif;  // Only meaningful in DIF_UNIFORM mode
	double constScaleFactor;  // Only meaningfull in LOG_UNIFORM mode
	bool bDependingValuesAreValid;
	// Members below are updated using updateDependencies()
	int nPyrImages;  // ==images.size()
	double scale_start;
	double scale_end;
	double relativeScaleRange;
	double logConstScaleFactor;
	vector<cv::Mat*> images;
	vector<double> scale;
};


#endif
