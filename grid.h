#ifndef GRID_H
#define GRID_H

#include "const_parameters.h"
#include "image_pyramid.h"

extern float q0;
extern float q1;
extern int nSections;
extern int pyramidAccessFactor;
extern int nValuesPerSubSection;
extern float floatNValuesPerSubSection;
extern int bitsPerSection;

extern int nL;
extern int qualityMode;
extern bool bDiscardMultipleHits;
extern int nExtractOnlyNBest;

struct IndexMapInfo 
{
	int index;
	int i1stBestMatch;
	int i2ndBestMatch;
	float p1;
	float p2;
	float power;
	float entropy;
	float quality;
};

bool pix(float& f, float x, float y, float* data, cv::Size& s);

bool getStrip(bool&bAllValid, float*f, int nF, float ax,
	          float ay, float bx, float by, float*data, cv::Size s);

bool extract_dtoken(uint64& dtoken, ImagePyramid& pyr,
	                cv::KeyPoint& a, cv::KeyPoint&b,
	                float* tmpIntensities, float* tmpCount, float* avg);

///////////////////////////////////////////////////////////////////////////////
template <class REAL> void resizeVectorByAveraging(REAL*dst, int nDst, REAL*src, int nSrc, REAL*countBuf);


/*****************************************************************************/
template <class REAL>
REAL sumVector(REAL* v, int n);

/*****************************************************************************/
template <class REAL>
bool getMinAndMaxOfVector(REAL& min, REAL& max, REAL* v, int n);

/*****************************************************************************/
template <class REAL> bool normalizeVectorMinMaxZeroOne(REAL* v, int n);

/*****************************************************************************/
// Use this only if you know that all the elements of v are non-negative (>= 0)
template <class REAL> REAL normalizeVector_noabs(REAL* v, int n);

/*****************************************************************************/
int IndexMapInfo_Quality_Cmp(const void*a, const void*b);

class DList 
{
public:
	DList() : nBuckets(0), nL(0), nPassivePairingsOfBucket(NULL), pairings(NULL) 
	{
	}
	~DList();

	// nL is the maximum length of a list of a bucket
	bool create(uint64 nBuckets_p, unsigned int nL_p);

	bool tryInsert(uint64 dtoken, int i, int j);

	// Returns size of list plus pointer to start via p,
	// plus the number of passive hits to the bucket
	unsigned int getListOfBucket(pair<unsigned int, unsigned int>*& p,
		                         unsigned int& nPassive, uint64 i);

	// Member variables
	uint64 nBuckets;  //number of buckets (= number of possible d-tokens)
	unsigned int nL;

	// nPairingsOfBuckets[i] is the number of pairings contained in bucket i.
	unsigned int* nPairingsOfBucket;

	// nPassivePairingsOfBuckets[i] is the total number of pairings that we
	// tried to insert into bucket i (but were unable to, because of the limited
	// size of each list). It always holds that
	// nPassivePairingsOfBuckets[i] >= nPairingsOfBuckets[i]
	unsigned int* nPassivePairingsOfBucket;

	// pairings[i*nL],...pairings[i*nL+nPairingsOfBucket[i]-1]
	// is the list of pairings of bucket i
	pair<unsigned int, unsigned int>* pairings;
};

/*****************************************************************************/
class Grid 
{
public:
	Grid(int nS, int nT);

	~Grid();

	bool vote(DList&s, DList&t);

	bool extractCorrespondenceHypotheses(std::vector<cv::DMatch>& matches,
		                                 int qualityMode, int nBest = -1,
		                                 bool bDiscardMultipleHits = true);

	bool visualize();

	IndexMapInfo* sourceIndexMapping;

	// Member variables
	// Consecutive values in memory correspond to to the same source node,
	// but to different target nodes
	float*votes;
	float*tmp;
	int nS;  // number of nodes in image A (source)
	int nT;  // number of nodes in image B (target)
	int n;   // number of cells
};


/*****************************************************************************/
/*!
//see http://compprog.wordpress.com/2007/10/08/generating-permutations-2/
Generates the next permutation of the vector v of length n.

@return true, if there are no more permutations to be generated
@return false, otherwise
*/
bool permutsDone(int v[], int n);

#endif
