#include "stdafx.h"
#include "image_pyramid.h"

ImagePyramid::ImagePyramid()
{
	scaleMode = UNDEFINED;
	constScaleDif = 0;
	constScaleFactor = 0;
	bDependingValuesAreValid = false;
	nPyrImages = 0;
	scale_start = 0;
	scale_end = 0;
	relativeScaleRange = 0;
}

ImagePyramid::~ImagePyramid()
{
	for (unsigned int i = 0; i<images.size(); i++)
	{
		images[i]->release();
	}
	images.clear();
}

bool ImagePyramid::create(cv::Mat *pXI, int nLayers, float lastRelativeScale,
	                      bool bLogarithmicallyEquallySpaced,
	                      bool bAlwaysScaleFromFirstImage)
{
	scaleMode = UNDEFINED;
	constScaleDif = 0;
	constScaleFactor = 0;
	cv::Mat *pXCur = pXI;
	pXCur->addref();
	images.push_back(pXCur);
	double scale0 = 1.0;
	scale.push_back(scale0);
	double lastScale = scale0 * lastRelativeScale;
	double q = 0;
	double scaleStep = 0;
	if (nLayers - 1 >= 1)
	{
		if (bLogarithmicallyEquallySpaced)
		{
			q = pow((double)lastScale / (double)scale0, 1 / (double)(nLayers - 1));
		}
		scaleStep = (scale0 - lastScale) / (double)(nLayers - 1);
		double curScale = scale0;
		for (int i = 0; i<nLayers - 1; i++)
		{
			double factor;
			if (bLogarithmicallyEquallySpaced)
			{
				if (!bAlwaysScaleFromFirstImage) {
					factor = q;
				}
				else {
					factor = pow(q, i + 1);
				}
			}
			else
			{
				double lastScale = curScale;
				curScale -= scaleStep;
				if (!bAlwaysScaleFromFirstImage) {
					factor = curScale / lastScale;
				}
				else {
					factor = curScale / scale0;
				}
			}
			cv::Mat*pXScaled = new cv::Mat();
			if (bAlwaysScaleFromFirstImage)
			{
				int dxOrg = pXI->size().width;
				int dyOrg = pXI->size().height;
				int w = ((float)dxOrg)*(float)factor;
				if ((int)(((float)(dxOrg - 1))*factor) == w) {
					w++;
				}
				int h = ((float)dyOrg)*(float)factor;
				cv::resize(*pXI, *pXScaled, cv::Size(w, h), factor, factor, CV_INTER_AREA);
				scale.push_back(factor);
			}
			else
			{
				cv::resize(*pXCur, *pXScaled, cv::Size(0, 0), factor, factor, CV_INTER_AREA);
				double lastScale = scale.back();
				scale.push_back(lastScale*factor);
			}
			pXCur = pXScaled;
			images.push_back(pXCur);
		}
	}
	if (bLogarithmicallyEquallySpaced)
	{
		scaleMode = LOG_UNIFORM;
		constScaleFactor = q;
	}
	else
	{
		scaleMode = DIF_UNIFORM;
		constScaleDif = scaleStep;
	}
	updateDependencies();
	return true;
}


bool ImagePyramid::getBestLevelIndex(int& k, double relativeDownScaleFactor)
{
	if (!nPyrImages) {
		return false;
	}
	if (nPyrImages == 1) {
		k = 0;
		return true;
	}
	switch (scaleMode) {
	case DIF_UNIFORM:
	{
		double relativeDownScaleRange = 1.0 - relativeDownScaleFactor;
		double fk = ((double)(nPyrImages - 1))*relativeDownScaleRange /
			relativeScaleRange;
		k = (int)fk;
	}
	break;
	case LOG_UNIFORM:
	{
		double fk = log(relativeDownScaleFactor) / logConstScaleFactor;
		k = (int)fk;
	}
	break;
	default:
		return false;
	}
	if (k<0) { k = 0; }
	if (k>nPyrImages - 1) { k = nPyrImages - 1; }
	return true;
}

void ImagePyramid::updateDependencies() 
{
	nPyrImages = images.size();
	if (!nPyrImages) {
		bDependingValuesAreValid = false;
		return;
	}
	scale_start = scale[0];
	scale_end = scale.back();
	relativeScaleRange = 1.0f - scale_end / scale_start;
	if (scaleMode == LOG_UNIFORM) {
		logConstScaleFactor = log(constScaleFactor);
	}
	else {
		logConstScaleFactor = 0;
	}
	bDependingValuesAreValid = true;
}
