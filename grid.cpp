#include "stdafx.h"
#include "grid.h"

float q0 = 0.1f;
float q1 = 0.8f;
int nSections = 9;
int pyramidAccessFactor = 8;
int nValuesPerSubSection;
float floatNValuesPerSubSection;
int bitsPerSection = 2;

int nL = 20;
int qualityMode = 1;

bool bDiscardMultipleHits = true;
int nExtractOnlyNBest = -1;

int IndexMapInfo_Quality_Cmp(const void*a, const void*b) 
{
	IndexMapInfo* s1 = (struct IndexMapInfo*) a;
	IndexMapInfo* s2 = (struct IndexMapInfo*) b;
	if (s1->quality < s2->quality) 
	{ 
		return 1; 
	}
	if (s1->quality > s2->quality) 
	{ 
		return -1;
	}
	return(0);
}

template <class REAL>
REAL sumVector(REAL* v, int n)
{
	REAL sum = 0;
	for (int i = 0; i<n; i++)
	{
		sum += v[i];
	}
	return sum;
}

template <class REAL>
bool getMinAndMaxOfVector(REAL& min, REAL& max, REAL* v, int n)
{
	min = FLT_MAX;
	max = -FLT_MAX;
	for (int i = 0; i<n; i++)
	{
		REAL a = v[i];
		if (a < min)
		{
			min = a;
		}
		if (a > max)
		{
			max = a;
		}
	}
	return n>0;
}

template <class REAL> bool normalizeVectorMinMaxZeroOne(REAL* v, int n)
{
	REAL min, max;
	if (!getMinAndMaxOfVector(min, max, v, n))
	{
		return false;
	}
	REAL d = max - min;
	if (d == 0)
	{
		return false;
	}
	REAL factor = (REAL)1.0 / d;
	for (int i = 0; i<n; i++)
	{
		REAL a = v[i];
		v[i] = (a - min)*factor;
	}
	return true;
}

template <class REAL> REAL normalizeVector_noabs(REAL* v, int n)
{
	REAL sum = sumVector(v, n);
	if (sum == 0)
	{
		return sum;
	}
	REAL fac = ((REAL)1) / sum;
	for (int i = 0; i < n; i++)
	{
		v[i] *= fac;
	}
	return sum;
}


bool pix(float& f, float x, float y, float* data, cv::Size& s) 
{
#define px(xx,yy) data[(yy)*s.width+xx]
	int iy = (int)y;
	int ix = (int)x;
	if (iy < 0 || iy >= s.height || ix < 0 || ix >= s.width) { return false; }
	if (y >= s.height) { y = (float)s.height - 1; }
	if (x >= s.width) { x = (float)s.width - 1; }
	float r = 1.0f - (y - (float)iy);
	float c = 1.0f - (x - (float)ix);
	float omc = 1.0f - c;
	float omr = 1.0f - r;
	float px_ixiy = px(ix, iy);
	float y1 = 0;
	float y2 = 0;
	y1 = c<1 ? c * px_ixiy + omc * px(ix + 1, iy) : px_ixiy;
	if (r < 1) {
		y2 = c < 1 ? c * px(ix, iy + 1) + omc * px(ix + 1, iy + 1) : px(ix, iy + 1);
	}
	f = r * y1 + omr * y2;
	return true;
}

/*****************************************************************************/
bool getStrip(bool&bAllValid, float*f, int nF, float ax, 
	          float ay, float bx, float by, float*data, cv::Size s) 
{
	if (!f) { return false; }
	if (nF <= 0) { return false; }
	float vx = bx - ax;
	float vy = by - ay;
	float px = ax;
	float py = ay;
	float vstepx, vstepy;
	if (nF == 1) {  //sample center (if nF==1)
		px = ax + vx*0.5f;
		py = by + vy*0.5f;
		vstepx = 0;
		vstepy = 0;
	}
	else {
		float factor = 1.0f / (float)(nF - 1);
		vstepx = vx*factor;
		vstepy = vy*factor;
	}
	bAllValid = true;
	for (int i = 0; i<nF; i++) {
		if (!pix(f[i], px, py, data, s)) {
			bAllValid = false;
			return false;
		}
		px += vstepx;
		py += vstepy;
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////
template <class REAL> void resizeVectorByAveraging(REAL*dst, int nDst, REAL*src, int nSrc, REAL*countBuf) 
{
	memset(dst, 0, sizeof(*dst)*nDst);
	memset(countBuf, 0, sizeof(*countBuf)*nDst);
	if (nSrc >= nDst) 
	{
		double indexStep = (double)(nDst - 1) / (double)(nSrc - 1);
		double k = 0;
		for (int i = 0; i<nSrc; i++) 
		{
			int ik = (int)(k + 0.001);
			dst[ik] += src[i];
			countBuf[ik] += 1.0;
			k += indexStep;
		}
	}
	else {
		double indexStep = (double)(nSrc - 1) / (double)(nDst - 1);
		double k = 0;
		for (int i = 0; i<nDst; i++) 
		{
			int ik = (int)(k + 0.001);
			dst[i] += src[ik];
			countBuf[i] += 1.0;
			k += indexStep;
		}
	}
	for (int i = 0; i<nDst; i++) 
	{
		dst[i] /= countBuf[i];
	}
}


/*****************************************************************************/
bool extract_dtoken(uint64& dtoken, ImagePyramid& pyr,
	                cv::KeyPoint& a, cv::KeyPoint&b,
	                float* tmpIntensities, float* tmpCount, float* avg)
{
	// Determine which pyramid layer to use, depeding on the length of line [ab]
	float dxab = b.pt.x - a.pt.x;
	float dyab = b.pt.y - a.pt.y;
	float l2ab = dxab*dxab + dyab*dyab;
	float lab = sqrt(l2ab);
	float relativeDownScaleFactor = ((float)nSections*pyramidAccessFactor) / lab;
	int k;
	if (!pyr.getBestLevelIndex(k, relativeDownScaleFactor))
	{
		cout << "failed to access image pyramid\n";
		return false;
	}

	cv::Mat*pXI = pyr.images[k];
	float*data = (float*)pXI->data;
	cv::Size s = pXI->size();
	//project point a and b onto the respective pyramid level
	float scale = (float)pyr.scale[k];
	float ax_pyr = a.pt.x*scale;
	float ay_pyr = a.pt.y*scale;
	float bx_pyr = b.pt.x*scale;
	float by_pyr = b.pt.y*scale;
	float l_pyr = lab*scale;
	float vx_pyr = bx_pyr - ax_pyr;
	float vy_pyr = by_pyr - ay_pyr;
	float sx = ax_pyr + vx_pyr*q0;
	float sy = ay_pyr + vy_pyr*q0;
	float ex = ax_pyr + vx_pyr*q1;
	float ey = ay_pyr + vy_pyr*q1;
	int nSamples = (int)l_pyr;

	if (nSamples<nSections)
	{
		nSamples = nSections;
	}
	bool bAllValid;
	if (!getStrip(bAllValid, tmpIntensities, nSamples, sx, sy, ex, ey, data, s))
	{
		cout << "failed to get strip\n";
		return false;
	}
	resizeVectorByAveraging(avg, nSections, tmpIntensities, nSamples, tmpCount);
	if (!normalizeVectorMinMaxZeroOne(avg, nSections))
	{
		for (int i = 0; i<nSections; i++)
		{
			avg[i] = 0.5f;
		}
	}

	// Quantize each section by descretizing its normalized average
	// intensity and concatenate the bits
	dtoken = 0;
	for (int i = 0; i<nSections; i++)
	{
		int value = (int)((avg[i] * floatNValuesPerSubSection));
		if (value<0)
		{
			value = 0;
		}
		if (value >= nValuesPerSubSection)
		{
			value = nValuesPerSubSection - 1;
		}
		dtoken <<= bitsPerSection;
		dtoken |= value;
	}
	return true;
}
/*****************************************************************************/


DList::~DList() 
{
	if (nPairingsOfBucket) 
	{ 
		delete[] nPairingsOfBucket; 
	}
	if (pairings) 
	{ 
		delete[] pairings;
	}
};

// nL is the maximum length of a list of a bucket
bool DList::create(uint64 nBuckets_p, unsigned int nL_p) {
#ifdef _MSV_VER
#pragma warning( push)
#pragma warning( disable : 4127)
#endif
	if (sizeof(size_t)<8 && nBuckets_p >= 0x100000000ULL) {
		cout << "Must run on a 64-bit machine with current parameter settings"
			<< endl;
		return false;
	}
#ifdef _MSV_VER
#pragma warning( pop)
#endif
	nPairingsOfBucket = new unsigned int[(size_t)nBuckets_p];
	if (!nPairingsOfBucket) {
		cout << "insufficient memory" << endl;
		return false;
	}
	nPassivePairingsOfBucket = new unsigned int[(size_t)nBuckets_p];
	if (!nPassivePairingsOfBucket) {
		cout << "insufficient memory" << endl;
		delete[] nPairingsOfBucket;
		return false;
	}
	memset(nPairingsOfBucket, 0, sizeof(*nPairingsOfBucket)*(size_t)nBuckets_p);
	memset(nPassivePairingsOfBucket, 0,
		sizeof(*nPassivePairingsOfBucket)*(size_t)nBuckets_p);
	pairings = new pair<unsigned int, unsigned int>[(size_t)nBuckets_p*nL_p];
	if (!pairings) {
		cout << "insufficient memory" << endl;
		delete[] nPairingsOfBucket;
		delete[] nPassivePairingsOfBucket;
		nPairingsOfBucket = NULL;
		return false;
	}
	nBuckets = nBuckets_p;
	nL = nL_p;
	return true;
}

bool DList::tryInsert(uint64 dtoken, int i, int j) 
{
	nPassivePairingsOfBucket[dtoken]++;
	if (nPairingsOfBucket[dtoken] >= nL) 
	{ 
		return false; 
	}
	pair<unsigned int, unsigned int>&p =
		pairings[dtoken*nL + nPairingsOfBucket[dtoken]++];
	p.first = i;
	p.second = j;
	return true;
}

unsigned int DList::getListOfBucket(pair<unsigned int, unsigned int>*& p, unsigned int& nPassive, uint64 i) 
{
	p = pairings + i*nL;
	nPassive = nPassivePairingsOfBucket[i];
	return nPairingsOfBucket[i];
}


Grid::Grid(int nS, int nT) : nS(nS), nT(nT) 
{
	n = nS*nT;
	votes = new float[n];
	memset(votes, 0, sizeof(*votes)*n);
	sourceIndexMapping = new IndexMapInfo[nS];
	tmp = new float[nT];
}

Grid::~Grid() 
{
	delete[] tmp;
	delete[] sourceIndexMapping;
	delete[] votes;
}

bool Grid::vote(DList&s, DList&t) 
{
	if (s.nBuckets != t.nBuckets) 
	{
		cout << "lists are not compatible" << endl;
		return false;
	}
	cout << endl;
	for (unsigned int i = 0; i<s.nBuckets; i++) 
	{
		pair<unsigned int, unsigned int>* p[2];
		unsigned int nPassivePairs[2];
		unsigned int nPairs[2] = {
			s.getListOfBucket(p[0],nPassivePairs[0],i),
			t.getListOfBucket(p[1],nPassivePairs[1],i)
		};
		unsigned int nv = nPassivePairs[0] * nPassivePairs[1];
		if (nv>0) 
		{
			float voting_power = 1.0f / (float)nv;
			for (unsigned int k = 0; k<nPairs[0]; k++) 
			{
				unsigned int is0 = p[0][k].first;
				unsigned int is1 = p[0][k].second;
				for (unsigned int j = 0; j<nPairs[1]; j++) 
				{
					unsigned int it0 = p[1][j].first;
					unsigned int it1 = p[1][j].second;
					votes[is0*nT + it0] += voting_power;
					votes[is1*nT + it1] += voting_power;
				}
			}
		}
		if (i % (s.nBuckets / 333) == 0 || i == s.nBuckets - 1) 
		{
			double percentage = (double)i / (double)(s.nBuckets - 1);
			cout << '\r' << fixed << showpoint
				<< setprecision(2) << percentage*100.0 << " % of votes cast";
		}
	}
	return true;
}

bool Grid::extractCorrespondenceHypotheses(std::vector<cv::DMatch>& matches,
	                                       int qualityMode, int nBest, bool bDiscardMultipleHits) 
{
	for (int is = 0; is<nS; is++) 
	{
		memcpy(tmp, votes + is*nT, sizeof(*votes)*nT);
		float power = normalizeVector_noabs(tmp, nT);
		float sum = 0;
		float log2 = log(2.0f);
		int iBest = -1;
		int iBest2 = -1;
		float pbest = 0;
		float pbest2 = 0;
		for (int i = 0; i<nT; i++) 
		{
			float p = tmp[i];
			if (p>pbest)
			{
				if (iBest != -1)
				{
					if (p>pbest2) {
						pbest2 = pbest;
						iBest2 = iBest;
					}
				}
				pbest = p;
				iBest = i;
			}
			else if (p>pbest2) 
			{
				pbest2 = p;
				iBest2 = i;
			}
			if (p != 0) {
				// See http://en.wikipedia.org/wiki/Entropy_%28information_theory%29
				// p*log(p) is taken to be zero, which is consistent with
				// the limit \limes_{p->0+}{p log(p)}=0
				sum += p*log(p);
			}
		}
		IndexMapInfo&m = sourceIndexMapping[is];
		m.index = is;
		m.entropy = -sum / log2;
		m.i1stBestMatch = iBest;
		m.i2ndBestMatch = iBest2;
		m.p1 = pbest;
		m.p2 = pbest2;
		m.power = power;
		switch (qualityMode) {
		case 0:
			if (m.entropy != 0) {
				m.quality = m.power*m.p1 / m.entropy;
			}
			else {
				m.quality = m.power*m.p1 / 0.1f;
			}
			break;
		case 1:
			m.quality = m.power*m.p1;
			break;
		case 2:
			m.quality = -m.entropy;
			break;
		case 3:
			m.quality = m.entropy;
			break;
		case 4:
			m.quality = m.p1;
			break;
		case 5:
			if (pbest2 != 0) {
				m.quality = m.p1 / m.p2;
			}
			else {
				m.quality = 0;
			}
			break;
		case 6:
			m.quality = m.p1 - m.p2;
			break;
		default:
			m.quality = m.p1;
		}
	}
	qsort(sourceIndexMapping, nS, sizeof(*sourceIndexMapping),
		IndexMapInfo_Quality_Cmp);
	bool* covered = NULL;
	if (bDiscardMultipleHits)
	{
		covered = new bool[nT];
		memset(covered, 0, sizeof(*covered)*nT);
	}
	int count = 0;
	for (int i = 0; i<nS; i++) 
	{
		IndexMapInfo&m = sourceIndexMapping[i];
		if (m.i1stBestMatch != -1) 
		{
			count++;
			if (nBest != -1 && count >= nBest) 
			{ 
				break;
			}
			if (bDiscardMultipleHits&&covered[m.i1stBestMatch]) 
			{ 
				continue; 
			}
			matches.push_back(cv::DMatch(m.index, m.i1stBestMatch, 1.0f / m.quality));
			if (bDiscardMultipleHits)
			{ 
				covered[m.i1stBestMatch] = true; 
			}
		}
	}
	if (covered) 
	{
		delete[] covered; 
	}
	return true;
}

bool Grid::visualize()
{
	if (!votes) 
	{ 
		return false; 
	}
	float min, max;
	if (!getMinAndMaxOfVector(min, max, votes, n)) 
	{ 
		return false; 
	}
	cv::Mat I(nS, nT, CV_32FC1, votes);
	cv::Mat I2;
	I.convertTo(I2, CV_32FC1, 1.0f / (max - min), -min);
	cout << "Voting Grid Window Define the grid.cpp rows: 539 ......" << endl;
	imshow("Voting Grid...", I2);
	return true;
}

bool permutsDone(int v[], int n)
{
	/* P2 */
	/* Find the largest i */
	int t;
	int i = n - 2;
	while ((i >= 0) && (v[i] > v[i + 1])) { --i; }
	/* If i is smaller than 0, then there are no more permutations. */
	if (i < 0) { return true; }
	/* Find the largest element after vi but not larger than vi */
	int k = n - 1;
	while (v[i] > v[k]) { --k; }
	t = v[i]; v[i] = v[k]; v[k] = t;

	/* Swap the last n - i elements. */
	int j;
	k = 0;
	for (j = i + 1; j < (n + i) / 2 + 1; ++j, ++k) {
		t = v[j]; v[j] = v[n - k - 1]; v[n - k - 1] = t;
	}

	return false;
}