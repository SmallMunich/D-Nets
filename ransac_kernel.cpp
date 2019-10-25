#include "stdafx.h"
#include "ransac_kernel.h"

float g_geometricVerification_threshDist = 16;

TransformationModelType gTransformationModel = HOMOGRAPHY; // AFFINE  HOMOGRAPHY

RansacKernel::RansacKernel()
{
	pCachedParts = NULL;
	nCachedParts = 0;
	nPartsForOneModel = 0;
	nBytesForOnePart = 0;
	nBytesForOneModel = 0;
}

RansacKernel::~RansacKernel()
{
	releasePartsMemory();
}

bool RansacKernel::cacheParts(std::vector<cv::KeyPoint> kp[2],
	                          std::vector< cv::DMatch >& candidate_matches,
	                          bool bUseGlobalFrame, bool bUseUndistortedData)
{
	return false;
}

void RansacKernel::releasePartsMemory()
{
	if (pCachedParts)
	{
		free(pCachedParts);
		pCachedParts = NULL;
	}
}

void* RansacKernel::cachedParts() 
{ 
	return pCachedParts; 
}

int RansacKernel::getNumCachedParts() 
{ 
	return nCachedParts; 
}

int RansacKernel::getNumPartsForOneModel() 
{ 
	return nPartsForOneModel; 
}

int RansacKernel::partSize() 
{
	return nBytesForOnePart; 
}

int RansacKernel::modelSize() 
{ 
	return nBytesForOneModel;
}

bool RansacKernel::getTransformation2D(class Transformation2D*& pTransformation2D, void*pModel)
{
	return false;
}


bool RansacKernel_2DPointTransformation::cacheParts(std::vector<cv::KeyPoint> kp[2],
	                                                std::vector< cv::DMatch >& candidate_matches,
	                                                bool bUseGlobalFrame, bool bUseUndistortedData)
{
	int nR = candidate_matches.size();
	releasePartsMemory();
	pCachedParts = malloc(nBytesForOnePart * nR);
	pCached2V2 = (TwoV2*)pCachedParts;
	V2 pointA, pointB;
	for (int i = 0; i<nR; i++)
	{
		int ia = candidate_matches[i].queryIdx;
		int ib = candidate_matches[i].trainIdx;
		V2&pointA = pCached2V2[i].a;
		V2&pointB = pCached2V2[i].b;
		pointA.x = kp[0].at(ia).pt.x;
		pointA.y = kp[0].at(ia).pt.y;
		pointB.x = kp[1].at(ib).pt.x;
		pointB.y = kp[1].at(ib).pt.y;
	}
	nCachedParts = nR;
	return true;
}


bool RansacKernel_AT2::construct(void* pModel, int* partIndices)
{
	AT2* m = (AT2*)pModel;
	TwoV2& pair0 = pCached2V2[partIndices[0]];
	TwoV2& pair1 = pCached2V2[partIndices[1]];
	TwoV2& pair2 = pCached2V2[partIndices[2]];
	return m->defineFromTriangles(pair0.a, pair1.a, pair2.a,
		pair0.b, pair1.b, pair2.b);
}

bool RansacKernel_AT2::accept(void* pModel, int iPart)
{
	AT2* m = (AT2*)pModel;
	TwoV2& pair = pCached2V2[iPart];
	V2 dif = pair.b - m->out(pair.a);
	float error = dif.len();
	return (error < thresh_dist);
}

bool RansacKernel_AT2::getTransformation2D(class Transformation2D*& pTransformation2D, void*pModel)
{
	const AT2* pAT2 = (AT2*)pModel;
	ATrans2* pATrans2 = new ATrans2();
	pATrans2->a[0] = pAT2->a[0];
	pATrans2->a[1] = pAT2->a[1];
	pATrans2->a[2] = pAT2->a[2];
	pATrans2->a[3] = pAT2->a[3];
	pATrans2->t = pAT2->t;
	pTransformation2D = pATrans2;
	return true;
}

bool RansacKernel_Homography::construct(void* pModel, int* partIndices)
{
	// Placement constructor (not well known),
	// constructes an object at a given memory address
	Homography*h = new (pModel) Homography();
	TwoV2&pair0 = pCached2V2[partIndices[0]];
	TwoV2&pair1 = pCached2V2[partIndices[1]];
	TwoV2&pair2 = pCached2V2[partIndices[2]];
	TwoV2&pair3 = pCached2V2[partIndices[3]];
	DV2 a[4];
	a[0] = pair0.a;
	a[1] = pair1.a;
	a[2] = pair2.a;
	a[3] = pair3.a;
	DV2 b[4];
	b[0] = pair0.b;
	b[1] = pair1.b;
	b[2] = pair2.b;
	b[3] = pair3.b;
	return h->defineFromPointCorrespondences(a, b, false);
}

void RansacKernel_Homography::destruct(void* pModel)
{
	((Homography*)pModel)->Homography::~Homography();
	//this is because we did the inplace construction in (construct).
}

bool RansacKernel_Homography::accept(void* pModel, int iPart)
{
	Homography*m = (Homography*)pModel;
	TwoV2&pair = pCached2V2[iPart];
	V2 dif = pair.b - m->out(pair.a);
	float error = dif.len();
	return (error<thresh_dist);
}

bool RansacKernel_Homography::getTransformation2D(class Transformation2D*& pTransformation2D, void*pModel)
{
	if (!pModel)
	{
		return false;
	}
	pTransformation2D = new Homography(*(Homography*)pModel);
	return true;
}