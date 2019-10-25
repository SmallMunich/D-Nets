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
		���ں���ȱʡֵ.h��.cpp��� һ�㺯��Ĭ�ϵ�ȱʡֵҪд�������ĵط�
		������д�ڶ���ĵط�����д��.h���治д��.cpp����������
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
