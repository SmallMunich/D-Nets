// D-Nets.cpp : �������̨Ӧ�ó������ڵ㡣
//
#include "stdafx.h"
#include "const_parameters.h"
#include "init_parameters.h"
#include "image_pyramid.h"
#include "feature_extractor.h"
#include "triangulator.h"
#include "transformation2d.h"
#include "ransac_kernel.h"
#include "show_matches.h"
#include "grid.h"

#include <fstream>


void writeMatrixToFile(cv::Mat& matrix, const char* filename);

//--------------------------- ��ͨ��ͼ��������ȡ------------------------------------------------------//
bool compute_twoViews_Coarse_Match(cv::Mat img[2], cv::Mat img_split[2],/*char** filename,*/ \
	                               std::vector<cv::KeyPoint> v[2], std::vector<cv::DMatch >& grid_matches);


//----------------------------Color Image �����ϲ�--------------------------------------------------//
bool merge_twoViews_Kpts_Match(std::vector<cv::KeyPoint> v_r[2], std::vector<cv::KeyPoint> v_g[2], std::vector<cv::KeyPoint> v_b[2], \
	                           std::vector<cv::DMatch > grid_matches_r, std::vector<cv::DMatch > grid_matches_g, std::vector<cv::DMatch > grid_matches_b, \
	                           std::vector<cv::KeyPoint> v[2], std::vector<cv::DMatch >& grid_matches);


bool sort_keypoints(cv::KeyPoint kpt1, cv::KeyPoint kpt2);
bool unique_kpts(cv::KeyPoint kpt1, cv::KeyPoint kpt2);

bool compute_Inliers_Homography(std::vector<cv::KeyPoint>& left_key_point, std::vector<cv::KeyPoint>& right_key_point, \
	                            std::vector<cv::DMatch >& grid_matches, cv::Mat& Homo);

/************************************* Color-D-Nets ****************************************/
int main(void) 
{
    //cv::initModule_nonfree();  // opencv nonfreeģ��  sift�㷨��ȡ������

	bool bWait = true;  // ��ʾ�ȴ�ƥ����...  true��ʾ  false����ʾ
	string outputFilename_VisualizedMatches; // ���ƥ��ͼ��·������

	//------ ��ʼ�� ------//
	init();
	//------ img[2] Ϊ����������ƥ��ͼ��
	cv::Mat img[2];
	char* filename[5];

	//------ ��ȡ�ͳ�ʼ��ͼ�� ǰ4�в���
	filename[1] = "..\\D-Nets\\image\\img1.png";
	//filename[1] = "..\\D-Nets\\image\\gold_medal_1.bmp" ;// tea_service_1 gold_medal_1
	filename[3] = "..\\D-Nets\\image\\img6.png";
	//filename[3] = "..\\D-Nets\\image\\gold_medal_2.bmp";// tea_service_2  gold_medal_2  
	filename[2] = "";    filename[4] = "";
	//------ ��ȡ�����ɫͼ�� ------------//
	cv::Mat img1 = cv::imread(filename[1], 1);
	cv::Mat img2 = cv::imread(filename[3], 1);

	//cv::cvtColor(img1, img1, CV_BGR2HSV);
	//cv::cvtColor(img2, img2, CV_BGR2HSV);

	if (img1.channels() < 3 || img2.channels() < 3)
	{
		printf("read input img1 or img2 is not color images......\n");
		system("pause");
		return 0;
	}
	//------ ��ȡ��������Ϊ�ڵ� v[2] �ֱ�Ϊ�洢����ͼ���������ʸ��
	std::vector<cv::KeyPoint> v_r[2];  // Rͨ��ͼ���ȡ��kpts  img1_r  img2_r
	std::vector<cv::KeyPoint> v_g[2];  // Gͨ��ͼ���ȡ��kpts
	std::vector<cv::KeyPoint> v_b[2];  // Bͨ��ͼ���ȡ��kpts 
	std::vector<cv::DMatch > grid_matches_r; // Rͨ��DMatchƥ��������
	std::vector<cv::DMatch > grid_matches_g; // Gͨ��DMatchƥ��������
	std::vector<cv::DMatch > grid_matches_b; // Bͨ��DMatchƥ��������
	//-------------------- Color Image Message Together ----------------//
	std::vector<cv::KeyPoint> v[2];
	std::vector<cv::DMatch > grid_matches;

	//------ Color Image Split Channels R-G-B -----//
	cv::Mat img1_rgb[3];   cv::Mat img2_rgb[3];
	//------------ R G B Single Channel Image -----------------//
	cv::Mat img_r[2];   cv::Mat img_g[2];  cv::Mat img_b[2];
	//------ ��ɫͼ��ͨ������R-G-B -------//
	cv::split(img1, img1_rgb);
	cv::split(img2, img2_rgb);
	//---------- Rͨ��ͼ�� -----------//
	img1_rgb[0].copyTo(img_r[0]);
	img2_rgb[0].copyTo(img_r[1]);
	//---------- Gͨ��ͼ�� -----------//
	img1_rgb[1].copyTo(img_g[0]);
	img2_rgb[1].copyTo(img_g[1]);
	//---------- Rͨ��ͼ�� -----------//
	img1_rgb[2].copyTo(img_b[0]);
	img2_rgb[2].copyTo(img_b[1]);

	//-------------------------- ����ͼ��ƥ��ģ��  ���Coarseƥ�� -------------------------//
	compute_twoViews_Coarse_Match(img, img_r,/*filename,*/ v_r, grid_matches_r);
	compute_twoViews_Coarse_Match(img, img_g,/*filename,*/ v_g, grid_matches_g);
	compute_twoViews_Coarse_Match(img, img_b,/*filename,*/ v_b, grid_matches_b);

	//-------------------------- ����Kpts�ϲ�--------------------------------//
	bool mark = merge_twoViews_Kpts_Match(v_r, v_g, v_b, grid_matches_r, grid_matches_g, grid_matches_b, v, grid_matches);

	if (!mark) {
		printf("Color Image Feature Extracts Failed!\n");
		return 0;
	}

	//--------------------- �����Ҫ��Ӧ��һ��������֤���裬�Ի�ø��õĽ��... -----------------//
	std::vector<cv::DMatch > final_matches; 
	// �����任ģ�Ͳ������ƽ���ƥ���Ե��ᴿ
	if (gTransformationModel != NONE) 
	{
		RansacKernel* pKernel = NULL;
		switch (gTransformationModel) 
		{
			case AFFINE: // ����任ģ�͹���
				pKernel = new RansacKernel_AT2(g_geometricVerification_threshDist);
				cout << "geometric verification with AFFINE model" << endl;
				break;
			case HOMOGRAPHY: // ͸��任ģ�͹���
				pKernel = new RansacKernel_Homography(g_geometricVerification_threshDist);
				cout << "geometric verification with HOMOGRAPHY model" << endl;
				break;
		}

		Transformation2D* pTransformation2D = NULL;
		geometrically_verify_matches(pTransformation2D, final_matches, v, grid_matches,
			                         pKernel, g_nChunkSize, g_bChunkedCreation);
		if (pTransformation2D) 
		{
			delete pTransformation2D;
		}
		if (pKernel) {
			delete pKernel;
			pKernel = NULL;
		}
	}
	else {
		cout << "no geometric verification" << endl;
		final_matches = grid_matches;
	}

	cout << endl;
	cout << "Matches extracted..." << endl;


	cv::Mat Homo;
	compute_Inliers_Homography(v[0], v[1], final_matches, Homo);
	//---- ���ͬ�ʾ���......
	for(int i=0; i < 3; ++i)	
	{
		 for(int j=0; j < 3; ++j)	
		 {
			 cout << " " << Homo.at<double>(i,j);
		 }
		 cout << endl;
	}
	//------ ��ͬ�ʾ���д���ļ���......
	const char* homoName = "..\\D-Nets\\image\\homo.txt";
	writeMatrixToFile(Homo, homoName);


	cv::Mat img_matches;
	//----- gVisualizationType ���ӻ���ʽ  ��������ƥ���Ի���������ʾ
	switch (gVisualizationType) 
	{
		case LINES:  
			drawMatches(img1, v[0], img2, v[1], final_matches, img_matches);
			//drawMatches(img1, v[0], img2, v[1], final_matches, img_matches, cv::Scalar(0, 0, 255), \
				cv::Scalar(0, 255, 0), vector<char>(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
			break;
		case MESHES:  // draw single points is not default 
			drawMatchesThroughMeshes(img1, v[0], img2, v[1], final_matches,
			                         img_matches, g_MatchColor, g_singlePointColor1, g_singlePointColor2); 
			break;
	}
	imshow("Matches", img_matches);

	//----- ���ƥ�������Ϊ�� ���ɹ� ��ƥ���ͼ�洢д��ͼ�񱣴� --------//
	//----- outputFilename_VisualizedMatches ƥ��ͼƬ����·�� ---------//
	outputFilename_VisualizedMatches = "..\\D-Nets\\image\\img_matches.png";
	if (!outputFilename_VisualizedMatches.empty()) 
	{
		cv::imwrite(outputFilename_VisualizedMatches.c_str(), img_matches);
		cout << endl;
		cout << "Visualization written to: "
			 << outputFilename_VisualizedMatches << endl;
	}
	cout << endl;
	if (bWait)
	{
		cout << endl
			 << "Press ESC to quit (while focus is in the top-level window.)" << endl;
		for (;;) 
		{
			int c = cv::waitKey(0);
			if ((char)c == 27) { 
				break; 
			}
		}
	}

	system("pause");
	return 0;
}


void writeMatrixToFile(cv::Mat& matrix, const char* filename)
{
	std::ofstream fout(filename);
	if (!fout)
	{
		std::cout << "file not open ..." << std::endl;
		return;
	}

	for (int i = 0; i < matrix.rows; ++i)
	{
		for (int j = 0; j < matrix.cols; ++j) 
		{
			fout << matrix.at<double>(i, j) << "\t";
		}
		fout << std::endl;
	}

	fout.close();
}

bool compute_twoViews_Coarse_Match(cv::Mat img[2], cv::Mat img_split[2], \
	                               std::vector<cv::KeyPoint> v[2], std::vector<cv::DMatch >& grid_matches)
{
	//------ pyr[2] ����ԭʼͼ��Ľ�����ͼ��
	ImagePyramid pyr[2];
	float sigma = 1.0;  // �߶Ȳ���
	//------ ����ͼ��·������ת��Ϊָ����ʽ ------//
	if (!readAndConvertImages(img, img_split/*filename*/)) {
		return 0;
	}

	//------ ����ͼ����ʾ(�ڳ�ʼ���߶Ȳ����뽫ͼ��ת��Ϊ������)
	//imshow("image 0", img[0]);
	//imshow("image 1", img[1]);

	int nMaxDTokens = 1 << (nSections*bitsPerSection);
	nValuesPerSubSection = 1 << bitsPerSection;
	floatNValuesPerSubSection = (float)nValuesPerSubSection;

	cout << endl;

	for (int i = 0; i<2; ++i)
	{
		// ������ȡ���� img[i]Ϊ����ͼ��  gKeypointExtractionTypeΪ��ȡ����  v[i]Ϊ��ȡ�������㼯
		feature_extractor_nodes(img[i], gKeypointExtractionType, v[i]);
		cout << v[i].size() << " interest points extracted." << endl;
	}
	cv::Mat smoothed_img[2];
	//------ ��ͼ����и�˹�˲�...
	if (sigma > 0)
	{
		//----- ���������sigma���������˴�С���и�˹�˲�...
		int kernel_size = (int)max(3.0f, 8 * sigma + 1.0f);
		for (int i = 0; i < 2; ++i)
		{
			GaussianBlur(img[i], smoothed_img[i], cv::Size(kernel_size, kernel_size),
				sigma, sigma);
		}
	}

	//------ �Ը�˹�˲����ͼ��smoothed_img[i]���н���������...
	for (int i = 0; i < 2; ++i)
	{
		// nLayers Ϊ�������flag == 'L'  nLayers = value
		pyr[i].create(&smoothed_img[i], nLayers, 1.0f / (float)nLayers, true, true);
	}
	// ��ʱ�洢ǿ�ȵı��û�����...
	int nMaxTmpIntensities = 0;
	for (int i = 0; i < 2; ++i)
	{
		cv::Size s = img[0].size();
		int nMaxCur = max(s.width, s.height) * 2; //(overestimate)
		if (nMaxCur > nMaxTmpIntensities)
		{
			nMaxTmpIntensities = nMaxCur;
		}
	}
	float* tmpIntensities = new float[nMaxTmpIntensities];
	float* tmpCount = new float[nMaxTmpIntensities];
	float* avg = new float[nSections];

	// Ϊ��ϣ�����ڴ�  ��ϣ���������벿����ɣ�����һ�����һ��ͼ����б�
	DList dlist[2];
	for (int i = 0; i < 2; ++i)
	{
		if (!dlist[i].create(nMaxDTokens, nL))
		{
			cout << "failed to create hash table, you need more memory" << endl;
			delete[] avg;
			delete[] tmpCount;
			delete[] tmpIntensities;
			return -1;
		}
	}

	//----- ����ͼ�ıߣ����ڵ�ĶԱ�(����ֻ��һ����ȫ�ġ����������Ĺ�ϵ)
	//----- ȷ�����Ǹ��Ե�d-token��Ϊɢ�м�����������d���ƽ�ÿ��ƥ����뵽��ϣ����
	//----- �����ʵ���У�����Ϊÿ����ϣ���һ�붼��һ�����ݽṹDList��ÿһ�붼����һ��ͼ��ķ����б�(ÿ��Ͱ)
	for (int i = 0; i < 2; ++i)
	{   //������ƥ��ͼ��...
		int nv = v[i].size(); // number of nodes
							  // total irreflexive relation
		uint64 nPairings = ((nv - 1)*nv);
		printf("\n");
		if (nPairings < 1000000) {
			printf("number of strips to extract for image %d: %lld\n", i, nPairings);
		}
		else {
			printf("number of strips to extract for image %d: %.1f MegaStrips",
				i, (double)nPairings / (double)1000000.0);
		}
		uint64 count = 0;
		cout << endl;
		for (int j = 0; j<nv; ++j)
		{
			for (int k = 0; k<nv; ++k)
			{
				if (k != j)
				{
					//------ �������������...
					if (count % 10000 == 0 || count == nPairings - 1)
					{
						double percentage = (double)count / (double)nPairings;
						cout << '\r' << fixed << showpoint
							<< setprecision(2) << percentage*100.0
							<< " % of d-tokens of image "
							<< i << " extracted";
					}
					count++;
					//------ ��ȡd-token...
					uint64 dtoken;
					if (!extract_dtoken(dtoken, pyr[i], v[i].at(j), v[i].at(k),
						tmpIntensities, tmpCount, avg))
					{
						cout << "failed to extract d-token" << endl;
						continue;
					}
					dlist[i].tryInsert(dtoken, j, k);
				}
			}
		}
	}
	//--------------------------------- Grid ---------------------------------------//
	Grid g(v[0].size(), v[1].size());
	if (!g.vote(dlist[0], dlist[1]))
	{
		cout << "voting failed." << endl;
	}
	// Voting Ʊ���㷨��ʾͼ��...
	g.visualize();
	// ��ȡһ���Լ���...
	vector<pair<unsigned int, unsigned int> > correspondence_hypotheses;
	vector<float> quality;
	//std::vector<cv::DMatch > grid_matches;
	if (!g.extractCorrespondenceHypotheses(grid_matches, qualityMode,
		nExtractOnlyNBest, bDiscardMultipleHits))
	{
		cout << "failed ot extract correspondence hypotheses." << endl;
	}
	cout << endl;
	cout << "preparing visualization..." << endl;

	delete[] avg;  delete[] tmpCount;
	delete[] tmpIntensities;
}



bool merge_twoViews_Kpts_Match(std::vector<cv::KeyPoint> v_r[2], std::vector<cv::KeyPoint> v_g[2], std::vector<cv::KeyPoint> v_b[2], \
	std::vector<cv::DMatch > grid_matches_r, std::vector<cv::DMatch > grid_matches_g, std::vector<cv::DMatch > grid_matches_b, \
	std::vector<cv::KeyPoint> v[2], std::vector<cv::DMatch >& grid_matches)
{
	if ((0 == v_r[0].size()) && (0 == v_g[0].size()) && (0 == v_b[0].size()) && \
		(0 == v_r[1].size()) && (0 == v_g[1].size()) && (0 == v_b[1].size())) 
	{
		return false;
	}
	//------ ��������DMatch���ݽṹ
	if ((0 == grid_matches_r.size()) && (0 == grid_matches_b.size()) && (0 == grid_matches_b.size()))
	{
		return false;
	}

	//---------- R G B Single Channel Image Data Struct Merge ---------//
	// ============ Keypoint Merge: v_r ��v_g ��v_b ==================//
	// ============ DMatch Merge: grid_matches_r�� grid_matches_g�� grid_matches_b =================//
	v[0].insert(v[0].end(), v_r[0].begin(), v_r[0].end());
	v[0].insert(v[0].end(), v_g[0].begin(), v_g[0].end());
	v[0].insert(v[0].end(), v_b[0].begin(), v_b[0].end());
	v[1].insert(v[1].end(), v_r[1].begin(), v_r[1].end());
	v[1].insert(v[1].end(), v_g[1].begin(), v_g[1].end());
	v[1].insert(v[1].end(), v_b[1].begin(), v_b[1].end());
	// ----------- Free Memory ----------------------// 
	v_r[0].clear(); v_g[0].clear(); v_b[0].clear();
	v_r[1].clear(); v_g[1].clear(); v_b[1].clear();
	grid_matches.insert(grid_matches.end(), grid_matches_r.begin(), grid_matches_r.end());
	grid_matches.insert(grid_matches.end(), grid_matches_g.begin(), grid_matches_g.end());
	grid_matches.insert(grid_matches.end(), grid_matches_b.begin(), grid_matches_b.end());
	//------------ Free Memory ---------------------//
	grid_matches_r.clear();
	grid_matches_g.clear();
	grid_matches_b.clear();
	//------------- Delete Same Keypoints -----------------//
	//sort(v[0].begin(), v[0].end(), sort_keypoints);
	//v[0].erase(unique(v[0].begin(), v[0].end(), unique_kpts), v[0].end());
	//sort(v[1].begin(), v[1].end(), sort_keypoints);
	//v[1].erase(unique(v[1].begin(), v[1].end(), unique_kpts), v[1].end());
	//------------- Delete Same DMatch -----------------//

	return true;
}


bool sort_keypoints(cv::KeyPoint kpt1, cv::KeyPoint kpt2)
{
	if (kpt1.pt.x < kpt2.pt.x)
		return true;
	else
		return false;
}

bool unique_kpts(cv::KeyPoint kpt1, cv::KeyPoint kpt2)
{
	if ((kpt1.pt.x == kpt2.pt.x) && (kpt1.pt.y == kpt2.pt.y)) {
		return true;
	}
	else
		return false;
}

bool compute_Inliers_Homography(std::vector<cv::KeyPoint>& left_key_point, std::vector<cv::KeyPoint>& right_key_point, \
	                            std::vector<cv::DMatch >& grid_matches, cv::Mat& Homo)
{
	if (0 == grid_matches.size())
		return false;

	std::vector<cv::Point2f> p_left_point;
	std::vector<cv::Point2f> p_right_point;

	//std::vector<cv::KeyPoint> left_key_point;
	//std::vector<cv::KeyPoint> right_key_point;
	cv::Mat H = cv::Mat::zeros(3, 3, CV_32F);
	int size = grid_matches.size();

	for (int i = 0; i < size; ++i)
	{
		p_left_point.push_back(left_key_point[grid_matches[i].queryIdx].pt);
		p_right_point.push_back(right_key_point[grid_matches[i].trainIdx].pt);
	}
	int npoints = grid_matches.size() / 2;
	cv::Mat status = cv::Mat::zeros(npoints, 1, CV_8UC3);
	const float MIN_H_ERROR = 2.50f; 
	H = cv::findHomography(p_left_point, p_right_point, CV_RANSAC, MIN_H_ERROR, status);

	Homo = H;

	return true;
}






