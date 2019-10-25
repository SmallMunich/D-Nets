#ifndef FEATURE_EXTRACTOR
#define FEATURE_EXTRACTOR
#include "const_parameters.h"

extern int FAST_thresh;
extern int g_dense_spacing;
extern bool g_dense_bAddGaussianNoise;
extern float g_dense_stdDev;
typedef enum { FAST, SIFT, DENSE_SAMPLING } KeypointExtractionType;
extern KeypointExtractionType gKeypointExtractionType;

bool feature_extractor_nodes(cv::Mat& img, KeypointExtractionType gKeypointExtractionType, std::vector<cv::KeyPoint>& v);

/*****************************************************************************/
/*
    extractNodes_FAST() 提取FAST角点存储按标准opencv数据矢量
*/
bool extractNodes_FAST(cv::Mat& img, std::vector<cv::KeyPoint>& v);

/*****************************************************************************/
bool extractNodes_SIFT(cv::Mat&img, std::vector<cv::KeyPoint>&v);

/*****************************************************************************/

bool extractNodes_DenseSampling(cv::Mat&img, std::vector<cv::KeyPoint>&v);

/*****************************************************************************/


#endif
