#include "stdafx.h"
#include "feature_extractor.h"
#include "init_parameters.h"

// FAST 角点提取算子阈值...
int FAST_thresh = 80;
// extern int FAST_thresh;

// DENSE SAMPLING...
int g_dense_spacing = 10;
bool g_dense_bAddGaussianNoise = true;
float g_dense_stdDev = 3.0f;

// typedef enum { FAST, SIFT, DENSE_SAMPLING } KeypointExtractionType;
KeypointExtractionType gKeypointExtractionType = FAST;

bool feature_extractor_nodes(cv::Mat& img, KeypointExtractionType gKeypointExtractionType, std::vector<cv::KeyPoint>& v)
{
	// gKeypointExtractionType 提取特征点方式  FAST DENSE_SAMPLE SIFT
	// FAST SIFT提取角点时候需要对输入图像进行改变图像depth与channels
	// DENSE_SAMPLE直接对输入图像进行操作
	switch (gKeypointExtractionType)
	{
	case FAST:
		if (!extractNodes_FAST(img, v)) { return 0; }
		break;
	case DENSE_SAMPLING:
		if (!extractNodes_DenseSampling(img, v)) { return 0; }
		break;
	case SIFT:
		if (!extractNodes_SIFT(img, v)) { return 0; }
		break;
	}
}

bool extractNodes_FAST(cv::Mat& img, std::vector<cv::KeyPoint>& v)
{
	int threshold = FAST_thresh;
	bool nonmaxSuppression = true;
	cv::FastFeatureDetector detector(threshold, nonmaxSuppression);
	cv::Mat img_8U;
	// 按照FAST提取角点对应的图像通道及其深度对输入图像进行改变
	if (!ensureGrayImageWithDepth(img_8U, CV_8U, img, NULL, false))
	{
		cout << "failed to convert image for point detection" << endl;
		return false;
	}
	detector.detect(img_8U, v);
	return true;
}

bool extractNodes_SIFT(cv::Mat&img, std::vector<cv::KeyPoint>&v)
{
	return true;
}
/*
bool extractNodes_SIFT(cv::Mat&img,std::vector<cv::KeyPoint>&v) 
{
	cv::Mat img_8U;
	if (!ensureGrayImageWithDepth(img_8U, CV_8U, img, NULL, false)) 
	{
	    cout << "failed to convert image for point detection" << endl;
	    return false;
	}

	cv::SIFT::CommonParams p_common;
	cv::SIFT::DetectorParams p_detector;

	p_common.nOctaves = 4;
	p_common.nOctaveLayers = 3;
	p_common.firstOctave = -1;
	p_common.angleMode = 0;
	p_detector.threshold = 0.04;
	p_detector.edgeThreshold = 10;
	int border = 0;
	cv::SiftFeatureDetector detector(p_detector, p_common);
	if (border == 0) {
	    detector.detect(img_8U, v);
	} else 
	{
	    std::vector<cv::KeyPoint> tmp;
	    detector.detect(img, tmp);
	    int n=tmp.size();
	    v.reserve(n);
	    int nLevels = p_common.nOctaves - p_common.firstOctave;
	    float*b=new float[nLevels];
	    float scalePerLevel = 2.0f;
	    for (int i=0; i<nLevels; i++) 
		{
	        b[i] = (float)border*std::pow(scalePerLevel,i);
	    }

	    float dxb = img.size().width - 1.0f;
	    float dyb = img.size().height - 1.0f;
	    for (int i=0; i<n; i++) 
		{
	      cv::KeyPoint&k = tmp[i];
	      cv::Point2f&p = k.pt;
	      int iL = k.octave - p_common.firstOctave;
	      float biL = b[iL];
	      if (p.x-biL < 0) 
		  { 
			  continue; 
		  }
	      if (p.y-biL < 0)
		  {
			  continue;
		  }
	      if (p.x+biL > dxb) 
		  { 
			  continue; 
		  }
	      if (p.y+biL > dyb) 
		  { 
			  continue;
		  }
	      v.push_back(k);
	   }
	    delete[] b;
	}

	return true;
}
*/


bool extractNodes_DenseSampling(cv::Mat&img, std::vector<cv::KeyPoint>&v)
{
	if (!(g_dense_spacing > 0))
	{
		printf("invalid spacing\n");
		return false;
	}
	int dx = img.size().width;
	int dy = img.size().height;
	if (dx<1 || dy<1)
	{
		printf("image must have at least 1 pixel\n");
		return false;
	}
	int nX = (int)((float)(dx - 1) / g_dense_spacing) + 1;
	int nY = (int)((float)(dy - 1) / g_dense_spacing) + 1;
	int n = nX*nY;

	float y = 0.0f;
	float var = g_dense_stdDev*g_dense_stdDev;
	int i = 0;
	for (int iy = 0; iy<nY; iy++)
	{
		float x = 0.0f;
		for (int ix = 0; ix<nX; ix++, i++)
		{
			if (g_dense_bAddGaussianNoise)
			{
				float gx = rndNormal((float)x, var);
				float gy = rndNormal((float)y, var);
				// Keep old if outside
				if (gx<0 || gx >= dx) { gx = x; }
				if (gy<0 || gy >= dy) { gy = y; }
				v.push_back(cv::KeyPoint(gx, gy, 0));
			}
			else
			{
				v.push_back(cv::KeyPoint((float)x, (float)y, 0));
			}
			x += g_dense_spacing;
		}
		y += g_dense_spacing;
	}
	return true;
}