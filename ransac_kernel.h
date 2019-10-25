#ifndef RANSAC_KERNEL
#define RANSAC_KERNEL

#include "const_parameters.h"
#include "transformation2d.h"

extern float g_geometricVerification_threshDist;

typedef enum { NONE, AFFINE, HOMOGRAPHY } TransformationModelType;
extern TransformationModelType gTransformationModel;

/*****************************************************************************/
class RansacKernel
{
public:
	RansacKernel();

	virtual ~RansacKernel();

	virtual bool cacheParts(std::vector<cv::KeyPoint> kp[2],
		                    std::vector< cv::DMatch >& candidate_matches,
		                    bool bUseGlobalFrame = false,
		                    bool bUseUndistortedData = false);
	/**
	* @fn construct constructs a model from n parts. (n=nPartsForOneModel())
	* The parts are provides in terms of an array of indices into the pParts array.
	* The constructed model is stored to *pModel
	* return true, if construction was sucessfull.
	*/
	virtual bool construct(void* pModel, int* partIndices) = 0;
	/**
	*called fore each model before memory for all models is released.
	*Does not need to be implemented in many cases, but is usefull
	*for destructing (in place constructed objects).
	*/
	virtual void destruct(void* pModel) {};
	/**
	* return whether part iPart in the cached data is in consensus with model pModel
	*/
	virtual bool accept(void* pModel, int iPart) = 0;

	void releasePartsMemory();
	/**
	*gets a pointe to the parts (basically the pParts pointer)
	*/
	void* cachedParts();
	/**
	* gets the number of parts cachedn
	*/
	int getNumCachedParts();
	/**
	* return the number of parts required to construct one model.
	* returns the member variable nPartsForModel basically.
	*/
	int getNumPartsForOneModel();
	/*
	* number of bytes required to store one part.
	* basically returns the member variable nBytesForOneModel
	*/
	int partSize();
	/*
	* number of bytes required to store one model;
	* basically return the member variable nBytesForOneModel
	*/
	int modelSize();
	/*
	* returns succces
	* pTransformation2D: a 2D transformation representing the model
	* the caller is responsible for the returned objekt
	* pModel: the model
	*/
	virtual bool getTransformation2D(class Transformation2D*& pTransformation2D, void*pModel);

protected:
	void* pCachedParts;  // Pointer to cached parts
						 // NULL if no parts cached or memory released.
	int nCachedParts;       // Number of parts cached one after another in
							// pCachedParts, each having size partSize()
	int nPartsForOneModel;  // Number of parts required to construct one model.
	int nBytesForOnePart;
	int nBytesForOneModel;
};
/*****************************************************************************/
class RansacKernel_2DPointTransformation : public RansacKernel
{
public:
	RansacKernel_2DPointTransformation(float thresh_dist = 4.0f) :
		thresh_dist(thresh_dist)
	{
		nBytesForOnePart = sizeof(TwoV2);
		nBytesForOneModel = 0;
		nPartsForOneModel = 0;  // Overwrite the last two in derived class
	}
	virtual bool cacheParts(std::vector<cv::KeyPoint> kp[2],
		                    std::vector< cv::DMatch >& candidate_matches,
		                    bool bUseGlobalFrame = false,
		                    bool bUseUndistortedData = false);
	// virtual bool construct(void* pModel, int* partIndices) { return true; }
	// virtual bool accept(void* pModel, int iPart) { return true; }
	TwoV2* pCached2V2;  //==pCachedParts
	float thresh_dist;
};
/*****************************************************************************/
class RansacKernel_AT2 : public RansacKernel_2DPointTransformation
{
public:
	RansacKernel_AT2(float thresh_dist = 4.0f) :
		RansacKernel_2DPointTransformation(thresh_dist)
	{
		nPartsForOneModel = 3;
		nBytesForOneModel = sizeof(AT2);
	}
	virtual bool construct(void* pModel, int* partIndices);

	virtual bool accept(void* pModel, int iPart);

	virtual bool getTransformation2D(class Transformation2D*& pTransformation2D, void*pModel);

};

/*****************************************************************************/
class RansacKernel_Homography : public RansacKernel_2DPointTransformation
{
public:
	RansacKernel_Homography(float thresh_dist = 4.0f) :
		RansacKernel_2DPointTransformation(thresh_dist)
	{
		nPartsForOneModel = 4;
		nBytesForOneModel = sizeof(Homography);
	}
	virtual bool construct(void* pModel, int* partIndices);

	virtual void destruct(void* pModel);

	virtual bool accept(void* pModel, int iPart);

	virtual bool getTransformation2D(class Transformation2D*& pTransformation2D, void*pModel);

};

#endif
