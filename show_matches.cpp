#include "stdafx.h"
#include "show_matches.h"

// Visualization ø… ”ªØ...
cv::Scalar color1 = cv::Scalar(0, 128, 255);
cv::Scalar color2 = cv::Scalar(0, 255, 0);
int g_displayMeshLineWidth = 2;
// Use cv::Scalar::all(-1) below if you prefer random colors.
cv::Scalar g_MatchColor = cv::Scalar(0, 0xff, 0xff);
cv::Scalar g_singlePointColor1 = cv::Scalar(0, 128, 255);
cv::Scalar g_singlePointColor2 = cv::Scalar(0, 255, 0);

int g_nChunkSize = 10;
bool g_bChunkedCreation = true;


VisualizationType gVisualizationType = MESHES; // LINES MESHES

Mesh::Mesh()
{
	bDirectIndicesToOrgPoints = true;
}

void Mesh::create_unique_point_lists(const vector<cv::KeyPoint>& keypoints_p,
	const vector<cv::DMatch>& matches1to2) {
	/*
	* Triangulation methods typically have problems with multiple points
	* at identical locations. Hence, as a preprocessing step, created
	* linked lists of points with identical locations, such that each list,
	* possibly with only 1 element, can be considered as a unique point
	* for triangulation.
	* Edges than apply for all members of the linked lists of the source
	* and target of an edge, respectively.
	*/
	vector<int> idx;  // To keep track of index permutations
	for (int i = 0; i<keypoints_p.size(); i++) {  //make a local copy first
		std::vector<int> map_i;
		getMappingOfIndex(map_i, matches1to2, i);
		if (map_i.size()) {
			idx.push_back(i);
		}
	}
	int n = idx.size();
	cv::KeyPoint* keypoints = new cv::KeyPoint[n];
	for (int i = 0; i<n; i++) {
		keypoints[i] = keypoints_p[idx[i]];
	}

	// Sorting them with total_order_2d will let points with identical
	// coordinates be at subsequent positions after sorting.
	sort2Arrays(keypoints, (int*)idx.data(), n, total_order_2D);

	unique_keypoint_lists.reserve(n);  // upper bound
	for (int i = 0; i<n; i++) {
		float xroot = keypoints[i].pt.x;
		float yroot = keypoints[i].pt.y;
		KeyPointLinkAtUniquePos* pRoot =
			new KeyPointLinkAtUniquePos(xroot, yroot, idx[i]);
		KeyPointLinkAtUniquePos* pTail = pRoot;
		unique_keypoint_lists.push_back(pRoot);
		i++;
		while (i<n) {
			float xnext = keypoints[i].pt.x;
			float ynext = keypoints[i].pt.y;
			if (xnext != xroot || ynext != yroot) { break; }
			KeyPointLinkAtUniquePos* pNext =
				new KeyPointLinkAtUniquePos(xnext, ynext, idx[i], pTail);
			pTail = pNext;
			unique_keypoint_lists.push_back(pNext);
			i++;
		}
	}
	delete[] keypoints;
}

void Mesh::triangulate(bool bBiDirectional)
{
	int nP = unique_keypoint_lists.size();
	V2* p = new V2[nP];
	for (int i = 0; i<nP; i++) {
		p[i].x = unique_keypoint_lists[i]->x;
		p[i].y = unique_keypoint_lists[i]->y;
	}
	Triangulator d2(false, true, p, nP);
	fx_set<R2> tmp_edges;
	for (int i = 0; i<d2.final_triangles.size(); i++) {
		I3& i3 = d2.final_triangles[i];
		if (bBiDirectional) {
			tmp_edges.insert(R2(i3.ia, i3.ib));
			tmp_edges.insert(R2(i3.ib, i3.ia));
			tmp_edges.insert(R2(i3.ia, i3.ic));
			tmp_edges.insert(R2(i3.ic, i3.ia));
			tmp_edges.insert(R2(i3.ib, i3.ic));
			tmp_edges.insert(R2(i3.ic, i3.ib));
		}
		else {
			if (!tmp_edges.contains(R2(i3.ia, i3.ib))) {
				tmp_edges.insert(R2(i3.ib, i3.ia));
			}
			if (!tmp_edges.contains(R2(i3.ia, i3.ic))) {
				tmp_edges.insert(R2(i3.ic, i3.ia));
			}
			if (!tmp_edges.contains(R2(i3.ib, i3.ic))) {
				tmp_edges.insert(R2(i3.ic, i3.ib));
			}
		}
	}
	for (fx_set<R2>::iterator i = tmp_edges.begin();
		i != tmp_edges.end();
		i++) {
		edges.push_back(*i);
	}
	delete[] p;
}

void Mesh::transfer(Mesh&m2, const vector<cv::DMatch>& matches1to2,
	                const vector<cv::KeyPoint>& keypoints2 ) 
{
	m2.org_keypoints = keypoints2;
	//for each edge in this mesh
	//use the mapping matches1to2 to map the indices of the points
	//and add a respective edges to mesh m2
	fx_set<R2> tmp_edges; //for m2
	for (int i = 0; i<edges.size(); i++) {
		int iLA = edges[i].ia;
		int iLB = edges[i].ib;
		KeyPointLinkAtUniquePos* pa = unique_keypoint_lists[iLA];
		while (pa) {
			KeyPointLinkAtUniquePos* pb = unique_keypoint_lists[iLB];
			while (pb) {
				//ids of the original keypoints...
				int idxa = pa->idx;
				int idxb = pb->idx;
				std::vector<int> idxa_maps;
				std::vector<int> idxb_maps;
				getMappingOfIndex(idxa_maps, matches1to2, idxa);
				getMappingOfIndex(idxb_maps, matches1to2, idxb);
				for (int k = 0; k<idxa_maps.size(); k++) {
					int ia_cur_map = idxa_maps[k];
					for (int j = 0; j<idxb_maps.size(); j++) {
						int ib_cur_map = idxb_maps[j];
						tmp_edges.insert(R2(ia_cur_map, ib_cur_map));
					}
				}
				pb = pb->next;
			}
			pa = pa->next;
		}
	}
	for (fx_set<R2>::iterator i = tmp_edges.begin();
		i != tmp_edges.end();
		i++) {
		m2.edges.push_back(*i);
	}
}

void Mesh::getEdgeIndices(std::vector<R2>&edges_p)
{
	//returned indices point into unique_keypoint_lists
	edges_p = edges;
}

void Mesh::drawIntoImage(cv::Mat&outImg, CvScalar color, int thickness, int lineType, int shift) 
{
	std::vector<R2 > edges;
	getEdgeIndices(edges);
	const int draw_multiplier = 1 << shift;
	for (int i = 0; i<edges.size(); i++) 
	{
		cv::Point2f a;
		cv::Point2f b;
		if (bDirectIndicesToOrgPoints)
		{
			int ia = edges[i].ia;
			int ib = edges[i].ib;
			a.x = org_keypoints[ia].pt.x*draw_multiplier;
			a.y = org_keypoints[ia].pt.y*draw_multiplier;
			b.x = org_keypoints[ib].pt.x*draw_multiplier;
			b.y = org_keypoints[ib].pt.y*draw_multiplier;
		}
		else 
		{
			int ila = edges[i].ia;
			int ilb = edges[i].ib;
			KeyPointLinkAtUniquePos* pLA = unique_keypoint_lists[ila];
			KeyPointLinkAtUniquePos* pLB = unique_keypoint_lists[ilb];
			a = cv::Point2f(pLA->x*draw_multiplier, pLA->y*draw_multiplier);
			b = cv::Point2f(pLB->x*draw_multiplier, pLB->y*draw_multiplier);
		}
		line(outImg, a, b, color, thickness, lineType, shift);
	}
	/*Point2f pt1 = kp1.pt,
	pt2 = kp2.pt,
	dpt2 = Point2f( std::min(pt2.x+outImg1.cols, float(outImg.cols-1)), pt2.y );

	line( outImg,
	Point(cvRound(pt1.x*draw_multiplier), cvRound(pt1.y*draw_multiplier)),
	Point(cvRound(dpt2.x*draw_multiplier), cvRound(dpt2.y*draw_multiplier)),
	color, 1, CV_AA, draw_shift_bits );*/
}

Mesh:: ~Mesh() 
{
	for (unsigned int i = 0; i < unique_keypoint_lists.size(); i++) 
	{
		KeyPointLinkAtUniquePos* tmp = unique_keypoint_lists[i];
		KeyPointLinkAtUniquePos* nxt;
		do {
			nxt = tmp->next;
			delete tmp;
			tmp = nxt;
		} while (tmp);
		unique_keypoint_lists.clear();
	}
}


void prepareImgAndDrawKeypointsForMeshes(const cv::Mat& img1,
	const vector<cv::KeyPoint>& keypoints1,
	const cv::Mat& img2,
	const vector<cv::KeyPoint>& keypoints2,
	cv::Mat& outImg,
	cv::Mat& outImg1,
	cv::Mat& outImg2,
	const cv::Scalar& singlePointColor1,
	const cv::Scalar& singlePointColor2,
	int flags)
{
	cv::Size size(img1.cols + img2.cols, MAX(img1.rows, img2.rows));
	if (flags & cv::DrawMatchesFlags::DRAW_OVER_OUTIMG)
	{
		if (size.width > outImg.cols || size.height > outImg.rows)
		{
			CV_Error(CV_StsBadSize,
				"outImg size less than needed to draw img1 and img2");
		}
		outImg1 = outImg(cv::Rect(0, 0, img1.cols, img1.rows));
		outImg2 = outImg(cv::Rect(img1.cols, 0, img2.cols, img2.rows));
	}
	else
	{
		outImg.create(size, CV_MAKETYPE(img1.depth(), 3));
		outImg1 = outImg(cv::Rect(0, 0, img1.cols, img1.rows));
		outImg2 = outImg(cv::Rect(img1.cols, 0, img2.cols, img2.rows));
		if (img1.type() == CV_8U) {
			cvtColor(img1, outImg1, CV_GRAY2BGR);
		}
		else {
			img1.copyTo(outImg1);
		}
		if (img2.type() == CV_8U) {
			cvtColor(img2, outImg2, CV_GRAY2BGR);
		}
		else {
			img2.copyTo(outImg2);
		}
	}

	// draw keypoints
	if (!(flags & cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS))
	{
		cv::Mat outImg1 = outImg(cv::Rect(0, 0, img1.cols, img1.rows));
		drawKeypoints(outImg1, keypoints1, outImg1, singlePointColor1,
			flags | cv::DrawMatchesFlags::DRAW_OVER_OUTIMG);
		cv::Mat outImg2 = outImg(cv::Rect(img1.cols, 0, img2.cols, img2.rows));
		drawKeypoints(outImg2, keypoints2, outImg2, singlePointColor2,
			flags | cv::DrawMatchesFlags::DRAW_OVER_OUTIMG);
	}
}


void drawMatchesThroughMeshes(const cv::Mat& img1,
	const vector<cv::KeyPoint>& keypoints1,
	const cv::Mat& img2,
	const vector<cv::KeyPoint>& keypoints2,
	const vector<cv::DMatch>& matches1to2,
	cv::Mat& outImg,
	const cv::Scalar& matchColor,
	const cv::Scalar& singlePointColor1,
	const cv::Scalar& singlePointColor2,
	const vector<char>& matchesMask,
	int flags)
{
	if (!matchesMask.empty() && matchesMask.size() != matches1to2.size())
	{
		CV_Error(CV_StsBadSize, "matchesMask must have the same size as matches1to2");
	}
	cv::Mat outImg1, outImg2;
	prepareImgAndDrawKeypointsForMeshes(img1, keypoints1, img2, keypoints2,
		                                outImg, outImg1, outImg2,
		                                singlePointColor1, singlePointColor2,
		                                flags);
	//drawMatches(img1,keypoints1,img2,keypoints2,matches1to2,outImg);

	CvRect r1;
	r1.x = 0;
	r1.y = 0;
	r1.width = img1.size().width;
	r1.height = img1.size().height;
	Mesh m1(r1, keypoints1, matches1to2);
	CvRect r2;
	r2.x = 0;
	r2.y = 0;
	r2.width = img2.size().width;
	r2.height = img2.size().height;

	Mesh m2;
	m1.transfer(m2, matches1to2, keypoints2);
	int lineWidth = g_displayMeshLineWidth;
	m1.drawIntoImage(outImg1, color1, lineWidth, CV_AA, 4);
	m2.drawIntoImage(outImg2, color2, lineWidth, CV_AA, 4);

	cv::RNG& rng = cv::theRNG();
	bool isRandMatchColor = matchColor == cv::Scalar::all(-1);

	for (int i = 0; i<matches1to2.size(); i++) {
		cv::Scalar color = isRandMatchColor ?
			cv::Scalar(rng(256), rng(256), rng(256)) :
			matchColor;
		const cv::KeyPoint&kp1 = keypoints1[matches1to2[i].queryIdx];
		const cv::KeyPoint&kp2 = keypoints2[matches1to2[i].trainIdx];
		drawKeypoint_xx(outImg1, kp1, color, flags, 4);
		drawKeypoint_xx(outImg2, kp2, color, flags, 4);
	}
}


void getMappingOfIndex(std::vector<int>& ibs,
	const vector<cv::DMatch>& matches1to2, int ia)
{
	for (int i = 0; i<matches1to2.size(); i++)
	{
		if (matches1to2[i].queryIdx == ia)
		{
			ibs.push_back(matches1to2[i].trainIdx);
		}
	}
}

bool total_order_2D(cv::KeyPoint& a, cv::KeyPoint& b) 
{
	if (a.pt.y < b.pt.y) 
	{ 
		return true; 
	}
	if (a.pt.y > b.pt.y) 
	{
		return false; 
	}
	if (a.pt.x < b.pt.x) 
	{ 
		return true; 
	}
	if (a.pt.x > b.pt.x) 
	{
		return false; 
	}
	return true;
}


bool factorial(int& result, int n)
{
	if (n<0)
	{
		return false;
	}
	if (n == 0)
	{
		result = 1;
		return true;
	}
	if (n>12)
	{
		printf("result would exceed interger limit\n");
		return false;
	}
	int m = n;
	result = 1;
	for (int i = 0; i<n - 1; i++)
	{
		result *= m;
		m--;
	}
	return true;
}

bool draw_k_from_n(int&result, int k, int n)
{
	if (n<0)
	{
		return false;
	}
	if (n == 0)
	{
		result = 0;
		return true;
	}
	if (!k)
	{
		result = 1;
		return true;
	}
	if (k<0 || k>n)
	{
		result = 0;
		return true;
	}
	int denom = 1;
	int mul = n;
	for (int i = 0; i<k; i++)
	{
		denom *= mul;
		mul--;
	}
	int fac_k;
	if (!factorial(fac_k, k))
	{
		return false;
	}
	result = denom / fac_k;
	return true;
}

int* create_draw_index_set(int& nIndexSets, int k, int outof_n)
{
	if (k == outof_n)
	{
		int* x_1 = new int[k];
		for (int i = 0; i<k; i++)
		{
			x_1[i] = i;
		}
		nIndexSets = 1;
		return x_1;
	}
	if (!draw_k_from_n(nIndexSets, k, outof_n))
	{
		return NULL;
	}
	if (nIndexSets == 0)
	{
		return NULL;
	}
	int* maxOfPos = new int[k];
	int* x = new int[nIndexSets*k];
#ifndef NDEBUG
	memset(x, 0xff, sizeof(int)*nIndexSets*k);//just for debugging, to see the memory operations better
#endif
	int* v = x;
	//setup the first row and initialize some variables
	for (int i = 0; i<k; i++)
	{
		maxOfPos[i] = outof_n - k + i;
		v[i] = i;
	}

	int j = 1;
	int* vn;
	int curpos = k - 1;
	do {
		vn = v + k;
		// copy the last line, at the same time detect the first value
		// that has reached its maximum
		int i = 0;
		bool bMaxReached = false;
		do {
			int a = v[i];
			if (a == maxOfPos[i])
			{
				bMaxReached = true;
				break;
			}
			vn[i] = a;
			i++;
		} while (i<k);
		if (bMaxReached)
		{
			//line feed
			int b = v[i - 1] + 1;
			vn[i - 1] = b;
			for (int r = i; r<k; r++)
			{
				vn[r] = ++b;
			}
			curpos = k - 1;
		}
		else
		{
			vn[curpos]++;
		}
		j++;
		v = vn;
	} while (j<nIndexSets);

	delete[] maxOfPos;
	return x;
}

///////////////////////////////////////////////////////////////////////////
bool consensus_better(int& a, int& b)
{
	return (a>b);
}

bool geometrically_verify_matches(Transformation2D*& pTransformation2D,
	                              std::vector< cv::DMatch >& final_matches,
	                              std::vector< cv::KeyPoint > kp[2],
	                              std::vector< cv::DMatch >& grid_matches,
	                              RansacKernel* pKernel,
	                              int nChunkSize, bool bChunkedCreation)
{
	if (!pKernel) { return false; }
	if (!pKernel->cacheParts(kp, grid_matches, false, false))
	{
		return false;
	}
	int nPartsForOneModel = pKernel->getNumPartsForOneModel();
	int nParts = pKernel->getNumCachedParts();
	if (!nParts)
	{
		printf("no parts?\n");
		return false;
	}
	int nEffectiveChunkSize;

	if (!bChunkedCreation)
	{
		nEffectiveChunkSize = nParts;
	}
	else
	{
		nEffectiveChunkSize = nChunkSize;
	}

	int nCompleteChunks = nParts / nEffectiveChunkSize;
	int nPartsOfLastIncompleteChunk = nParts % nEffectiveChunkSize;

	int nModelsPerCompleteChunk;
	if (!draw_k_from_n(nModelsPerCompleteChunk,
		               nPartsForOneModel, nEffectiveChunkSize))
	{
		return false;
	}
	int nModelsForLastChunk = 0;
	if (!draw_k_from_n(nModelsForLastChunk,
		               nPartsForOneModel, nPartsOfLastIncompleteChunk))
	{
		return false;
	}
	int nExpectedModels = nCompleteChunks * nModelsPerCompleteChunk + nModelsForLastChunk;
	if (!nExpectedModels)
	{
		printf("nExpectedModels==0\n");
		return false;
	}
	int* creationIndices = new int[nExpectedModels*nPartsForOneModel];
	int nModels = 0;
	int nCurModels;
	int* curCreationIndices = create_draw_index_set(nCurModels,
		                                            nPartsForOneModel, nEffectiveChunkSize);
	if (!curCreationIndices)
	{
		printf("failed to create creation index array\n");
		assert(false);
		return false;
	}
	if (!nCurModels)
	{
		printf("nCurModels==0\n");
		assert(false);
		delete[] curCreationIndices;
		return false;
	}
	for (int i = 0; i<nCompleteChunks; i++)
	{
		//transfer indices to overall array
		int nCurIndices = nCurModels * nPartsForOneModel;
		int iStartIndex = nModels * nPartsForOneModel;
		int indexShift = i * nEffectiveChunkSize;
		for (int k = 0; k < nCurIndices; k++)
		{
			creationIndices[iStartIndex + k] = curCreationIndices[k] + indexShift;
		}
		nModels += nCurModels;
	}
	delete[] curCreationIndices;
	if (nModelsForLastChunk > 0)
	{
		curCreationIndices = create_draw_index_set(nCurModels, nPartsForOneModel,
			                                       nPartsOfLastIncompleteChunk);
		if (!curCreationIndices)
		{
			printf("failed to create creation index array\n");
			assert(false);
			return false;
		}
		if (!nCurModels)
		{
			printf("nCurModels==0\n");
			assert(false);
			delete[] curCreationIndices;
			return false;
		}
		// Transfer indices to overall array
		int nCurIndices = nCurModels*nPartsForOneModel;
		int iStartIndex = nModels*nPartsForOneModel;
		int indexShift = nCompleteChunks*nEffectiveChunkSize;
		for (int k = 0; k < nCurIndices; k++)
		{
			creationIndices[iStartIndex + k] = curCreationIndices[k] + indexShift;
		}
		delete[] curCreationIndices;
		nModels += nCurModels;
	}
	assert(nModels == nExpectedModels);
	// Allocate memory for models and create the models
	int nBytesPerModel = pKernel->modelSize();
	char* modelMem = (char*)malloc(nModels*nBytesPerModel);
	bool* modelValid = new bool[nModels];
	char* m = modelMem;
	int* ic = creationIndices;
	for (int i = 0; i<nModels; i++)
	{
		modelValid[i] = pKernel->construct(m, ic);
		ic += nPartsForOneModel;
		m += nBytesPerModel;
	}
	// Test each model
	int* nConsens = new int[nModels];
	memset(nConsens, 0, sizeof(*nConsens)*nModels);

	for (int i = 0; i<nModels; i++)
	{
		if (!modelValid[i])
		{
			continue;
		}
		for (int k = 0; k<nParts; k++)
		{
			bool bPartWasUsedForCreation = false;
			for (int j = 0; j<nPartsForOneModel; j++)
			{
				if (k == creationIndices[i*nPartsForOneModel + j])
				{
					bPartWasUsedForCreation = true;
				}
			}
			if (bPartWasUsedForCreation) { continue; }
			if (pKernel->accept(modelMem + i*nBytesPerModel, k))
			{
				nConsens[i]++;
			}
		}
	}
	int* modelIndex = new int[nModels];
	for (int i = 0; i<nModels; i++)
	{
		modelIndex[i] = i;
	}
	sort2Arrays(nConsens, modelIndex, nModels, consensus_better);
	// The best model is:
	int iBestModel = modelIndex[0];
	if (modelValid[iBestModel])
	{
		//determine again all parts that consens for this model to create the output
		for (int k = 0; k<nParts; k++)
		{
			if (pKernel->accept(modelMem + iBestModel*nBytesPerModel, k))
			{
				cv::DMatch&m = grid_matches[k];
				int ia = m.queryIdx;
				int ib = m.trainIdx;
				//pXC->getIndices(ia,ib,k);
				float v = m.distance;
				/*if (bHoldsVotesOrDistances&&!pXC->getVoteOrDistance(v,k)){
				deep_err("fields says that it holds votes or distance, but getting them failed for a structure");
				}
				pXR->addRelation(ia,ib,v);*/
				final_matches.push_back(cv::DMatch(ia, ib, v));
			}
		}
		if (pKernel->getTransformation2D(pTransformation2D,
			                             modelMem + iBestModel*nBytesPerModel))
		{
			//pOutputPinXTransformation2D->setData(pTransformation2D);
			//pTransformation2D->Release(_ev_);
		}
	}
	else
	{
		printf("no valid model was constructed\n");
	}
	for (int i = 0; i<nModels; i++)
	{
		pKernel->destruct(modelMem + i*nBytesPerModel);
	}
	delete[] modelIndex;
	delete[] nConsens;
	delete[] modelValid;
	free(modelMem);
	delete[] creationIndices;
	return true;
}

/*...END SOURCES FOR GEOMETRIC VERIFICATION*/