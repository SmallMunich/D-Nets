#ifndef DRAW_MATCHES
#define DRAW_MATCHES
#include "const_parameters.h"
#include "triangulator.h"
#include "transformation2d.h"
#include "ransac_kernel.h"

extern int g_displayMeshLineWidth;
extern cv::Scalar color1;
extern cv::Scalar color2;
// Use cv::Scalar::all(-1) below if you prefer random colors.
extern cv::Scalar g_MatchColor;
extern cv::Scalar g_singlePointColor1;
extern cv::Scalar g_singlePointColor2;

extern int g_nChunkSize;
extern bool g_bChunkedCreation;

typedef enum { LINES, MESHES } VisualizationType;
extern VisualizationType gVisualizationType;

//internal function
template <class A, class B>
int divide2arrays(A* a, B* b, int left, int right,
	bool(*a1_lessorequal_a2)(A& a1, A& a2)) 
{
	int i = left;
	int j = right - 1;
	A pivot = a[right];
	do {
		while (a1_lessorequal_a2(a[i], pivot) && i<right) { i++; }
		while (!a1_lessorequal_a2(a[j], pivot) && j>left) { j--; }
		if (i < j) {
			swap(a[i], a[j]);
			swap(b[i], b[j]);
		}
	} while (i<j);
	if (a1_lessorequal_a2(pivot, a[i])) 
	{
		swap(a[i], a[right]);
		swap(b[i], b[right]);
	}
	return i;
}
/*****************************************************************************/
template <class A, class B>
void sort2ArraysInternal(A* a, B* b, int left, int right,
	bool(*a1_lessorequal_a2)(A& a1, A& a2)) 
{
	while (right > left) 
	{
		int iDivide = divide2arrays(a, b, left, right, a1_lessorequal_a2);
		if (right - iDivide > iDivide - left) 
		{
			sort2ArraysInternal(a, b, left, iDivide - 1, a1_lessorequal_a2);
			left = iDivide + 1;
		}
		else 
		{
			sort2ArraysInternal(a, b, iDivide + 1, right, a1_lessorequal_a2);
			right = iDivide - 1;
		}
	}
}
/*****************************************************************************/
/* sorts two arrays of the same size simultaneously where the sorting is
* induced by the first array.
* Assume that you have two arrays
* A a[n];
* B b[n];
* both having size n.
* Assume further that you want to sort A based on a function that is
* able to compare to elements of a (a1_lessorequal_a2).
* The function sorts a according to this function.
* but at the same time permutates the elements in array b, accordingly.
*/
template <class A, class B>
void sort2Arrays(A* a, B* b, int n,
	bool(*a1_lessorequal_a2)(A& a1, A& a2)) 
{
	sort2ArraysInternal(a, b, 0, n - 1, a1_lessorequal_a2);
}

class KeyPointLinkAtUniquePos : public CvPoint2D32f
{
public:
	KeyPointLinkAtUniquePos(float x_p, float y_p, int idx,
		                    KeyPointLinkAtUniquePos*next = NULL) : idx(idx), next(next)
	{
		x = x_p; y = y_p;
	}
	int idx;//index of original keypoin
	class KeyPointLinkAtUniquePos* next;
};

/*****************************************************************************/

class Mesh
{
public:
	Mesh();
	Mesh(CvRect imageRect, const vector<cv::KeyPoint>& keypoints_p,
		const vector<cv::DMatch>& matches1to2) : imageRect(imageRect) {
		bDirectIndicesToOrgPoints = false;
		org_keypoints = keypoints_p;
		create_unique_point_lists(keypoints_p, matches1to2);
		triangulate();
	}

	void create_unique_point_lists(const vector<cv::KeyPoint>& keypoints_p,
		const vector<cv::DMatch>& matches1to2);

	void triangulate(bool bBiDirectional = false);

	void transfer(Mesh&m2, const vector<cv::DMatch>& matches1to2,
		const vector<cv::KeyPoint>& keypoints2);

	void getEdgeIndices(std::vector<R2>&edges_p);

	void drawIntoImage(cv::Mat&outImg, CvScalar color,
		int thickness, int lineType, int shift);

	virtual ~Mesh();
protected:
	CvRect imageRect;
	vector<cv::KeyPoint> org_keypoints;
	std::vector<KeyPointLinkAtUniquePos*> unique_keypoint_lists;
	std::vector<R2> edges;  // Indexes into unique_keypoint_lists or
							// org_keypoints
							// (if bDirectIndicesToOrgPoints is true)
	bool bDirectIndicesToOrgPoints;
};
/*****************************************************************************/
void prepareImgAndDrawKeypointsForMeshes(const cv::Mat& img1,
	const vector<cv::KeyPoint>& keypoints1,
	const cv::Mat& img2,
	const vector<cv::KeyPoint>& keypoints2,
	cv::Mat& outImg,
	cv::Mat& outImg1,
	cv::Mat& outImg2,
	const cv::Scalar& singlePointColor1,
	const cv::Scalar& singlePointColor2,
	int flags);

/*****************************************************************************/
static inline void drawKeypoint_xx(cv::Mat& img, const cv::KeyPoint& p,
	const cv::Scalar& color, int flags,
	int draw_shift_bits)
{
	int draw_multiplier = 1 << draw_shift_bits;
	CV_Assert(!img.empty());
	cv::Point center(cvRound(p.pt.x * draw_multiplier),
		cvRound(p.pt.y * draw_multiplier));
	// draw center with R=3
	int radius = 1 * draw_multiplier;
	cv::circle(img, center, radius, color, 1, CV_AA, draw_shift_bits);
}
/*****************************************************************************/

void drawMatchesThroughMeshes(const cv::Mat& img1,
	const vector<cv::KeyPoint>& keypoints1,
	const cv::Mat& img2,
	const vector<cv::KeyPoint>& keypoints2,
	const vector<cv::DMatch>& matches1to2,
	cv::Mat& outImg,
	const cv::Scalar& matchColor = cv::Scalar::all(-1),
	const cv::Scalar& singlePointColor1 = cv::Scalar::all(-1),
	const cv::Scalar& singlePointColor2 = cv::Scalar::all(-1),
	const vector<char>& matchesMask = vector<char>(),
	int flags = cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);

/*****************************************************************************/
void getMappingOfIndex(std::vector<int>& ibs,
	const vector<cv::DMatch>& matches1to2, int ia);

bool total_order_2D(cv::KeyPoint& a, cv::KeyPoint& b);


bool factorial(int& result, int n);

bool draw_k_from_n(int&result, int k, int n);

int* create_draw_index_set(int& nIndexSets, int k, int outof_n);

bool consensus_better(int& a, int& b);

bool geometrically_verify_matches(Transformation2D*& pTransformation2D,
	std::vector< cv::DMatch >& final_matches,
	std::vector< cv::KeyPoint > kp[2],
	std::vector< cv::DMatch >& grid_matches,
	RansacKernel* pKernel,
	int nChunkSize, bool bChunkedCreation);


#endif