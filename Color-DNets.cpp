// D-Nets.cpp : 定义控制台应用程序的入口点。
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

//--------------------------- 单通道图像特征提取------------------------------------------------------//
bool compute_twoViews_Coarse_Match(cv::Mat img[2], cv::Mat img_split[2],/*char** filename,*/ \
	                               std::vector<cv::KeyPoint> v[2], std::vector<cv::DMatch >& grid_matches);


//----------------------------Color Image 特征合并--------------------------------------------------//
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
    //cv::initModule_nonfree();  // opencv nonfree模块  sift算法提取特征点

	bool bWait = true;  // 显示等待匹配结果...  true显示  false不显示
	string outputFilename_VisualizedMatches; // 输出匹配图像路径设置

	//------ 初始化 ------//
	init();
	//------ img[2] 为输入两幅待匹配图像
	cv::Mat img[2];
	char* filename[5];

	//------ 读取和初始化图像 前4行参数
	filename[1] = "..\\D-Nets\\image\\img1.png";
	//filename[1] = "..\\D-Nets\\image\\gold_medal_1.bmp" ;// tea_service_1 gold_medal_1
	filename[3] = "..\\D-Nets\\image\\img6.png";
	//filename[3] = "..\\D-Nets\\image\\gold_medal_2.bmp";// tea_service_2  gold_medal_2  
	filename[2] = "";    filename[4] = "";
	//------ 读取输入彩色图像 ------------//
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
	//------ 提取特征点作为节点 v[2] 分别为存储两幅图像的特征点矢量
	std::vector<cv::KeyPoint> v_r[2];  // R通道图像获取的kpts  img1_r  img2_r
	std::vector<cv::KeyPoint> v_g[2];  // G通道图像获取的kpts
	std::vector<cv::KeyPoint> v_b[2];  // B通道图像获取的kpts 
	std::vector<cv::DMatch > grid_matches_r; // R通道DMatch匹配描述子
	std::vector<cv::DMatch > grid_matches_g; // G通道DMatch匹配描述子
	std::vector<cv::DMatch > grid_matches_b; // B通道DMatch匹配描述子
	//-------------------- Color Image Message Together ----------------//
	std::vector<cv::KeyPoint> v[2];
	std::vector<cv::DMatch > grid_matches;

	//------ Color Image Split Channels R-G-B -----//
	cv::Mat img1_rgb[3];   cv::Mat img2_rgb[3];
	//------------ R G B Single Channel Image -----------------//
	cv::Mat img_r[2];   cv::Mat img_g[2];  cv::Mat img_b[2];
	//------ 彩色图像通道分离R-G-B -------//
	cv::split(img1, img1_rgb);
	cv::split(img2, img2_rgb);
	//---------- R通道图像 -----------//
	img1_rgb[0].copyTo(img_r[0]);
	img2_rgb[0].copyTo(img_r[1]);
	//---------- G通道图像 -----------//
	img1_rgb[1].copyTo(img_g[0]);
	img2_rgb[1].copyTo(img_g[1]);
	//---------- R通道图像 -----------//
	img1_rgb[2].copyTo(img_b[0]);
	img2_rgb[2].copyTo(img_b[1]);

	//-------------------------- 计算图像匹配模块  输出Coarse匹配 -------------------------//
	compute_twoViews_Coarse_Match(img, img_r,/*filename,*/ v_r, grid_matches_r);
	compute_twoViews_Coarse_Match(img, img_g,/*filename,*/ v_g, grid_matches_g);
	compute_twoViews_Coarse_Match(img, img_b,/*filename,*/ v_b, grid_matches_b);

	//-------------------------- 特征Kpts合并--------------------------------//
	bool mark = merge_twoViews_Kpts_Match(v_r, v_g, v_b, grid_matches_r, grid_matches_g, grid_matches_b, v, grid_matches);

	if (!mark) {
		printf("Color Image Feature Extracts Failed!\n");
		return 0;
	}

	//--------------------- 如果需要，应用一个几何验证步骤，以获得更好的结果... -----------------//
	std::vector<cv::DMatch > final_matches; 
	// 经过变换模型参数估计进行匹配点对的提纯
	if (gTransformationModel != NONE) 
	{
		RansacKernel* pKernel = NULL;
		switch (gTransformationModel) 
		{
			case AFFINE: // 仿射变换模型估计
				pKernel = new RansacKernel_AT2(g_geometricVerification_threshDist);
				cout << "geometric verification with AFFINE model" << endl;
				break;
			case HOMOGRAPHY: // 透射变换模型估计
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
	//---- 输出同质矩阵......
	for(int i=0; i < 3; ++i)	
	{
		 for(int j=0; j < 3; ++j)	
		 {
			 cout << " " << Homo.at<double>(i,j);
		 }
		 cout << endl;
	}
	//------ 将同质矩阵写入文件中......
	const char* homoName = "..\\D-Nets\\image\\homo.txt";
	writeMatrixToFile(Homo, homoName);


	cv::Mat img_matches;
	//----- gVisualizationType 可视化形式  点线连接匹配点对还是网格显示
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

	//----- 如果匹配点数不为空 即成功 将匹配后图存储写入图像保存 --------//
	//----- outputFilename_VisualizedMatches 匹配图片保存路径 ---------//
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
	//------ pyr[2] 输入原始图像的金字塔图像
	ImagePyramid pyr[2];
	float sigma = 1.0;  // 尺度参数
	//------ 输入图像路径将其转换为指定格式 ------//
	if (!readAndConvertImages(img, img_split/*filename*/)) {
		return 0;
	}

	//------ 输入图像显示(在初始化尺度参数与将图像转换为浮点型)
	//imshow("image 0", img[0]);
	//imshow("image 1", img[1]);

	int nMaxDTokens = 1 << (nSections*bitsPerSection);
	nValuesPerSubSection = 1 << bitsPerSection;
	floatNValuesPerSubSection = (float)nValuesPerSubSection;

	cout << endl;

	for (int i = 0; i<2; ++i)
	{
		// 特征提取函数 img[i]为输入图像  gKeypointExtractionType为提取类型  v[i]为提取的特征点集
		feature_extractor_nodes(img[i], gKeypointExtractionType, v[i]);
		cout << v[i].size() << " interest points extracted." << endl;
	}
	cv::Mat smoothed_img[2];
	//------ 对图像进行高斯滤波...
	if (sigma > 0)
	{
		//----- 根据输入的sigma来计算卷积核大小进行高斯滤波...
		int kernel_size = (int)max(3.0f, 8 * sigma + 1.0f);
		for (int i = 0; i < 2; ++i)
		{
			GaussianBlur(img[i], smoothed_img[i], cv::Size(kernel_size, kernel_size),
				sigma, sigma);
		}
	}

	//------ 对高斯滤波后的图像smoothed_img[i]进行金字塔构建...
	for (int i = 0; i < 2; ++i)
	{
		// nLayers 为输入参数flag == 'L'  nLayers = value
		pyr[i].create(&smoothed_img[i], nLayers, 1.0f / (float)nLayers, true, true);
	}
	// 暂时存储强度的备用缓冲区...
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

	// 为哈希表保存内存  哈希表由两个半部分组成，其中一半包含一个图像的列表
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

	//----- 构建图的边，即节点的对边(这里只是一个完全的、不可伸缩的关系)
	//----- 确定它们各自的d-token作为散列键，并根据其d令牌将每个匹配插入到哈希表中
	//----- 在这个实现中，我们为每个哈希表的一半都有一个数据结构DList，每一半都保存一个图像的反向列表(每个桶)
	for (int i = 0; i < 2; ++i)
	{   //两幅待匹配图像...
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
					//------ 计算与输出过程...
					if (count % 10000 == 0 || count == nPairings - 1)
					{
						double percentage = (double)count / (double)nPairings;
						cout << '\r' << fixed << showpoint
							<< setprecision(2) << percentage*100.0
							<< " % of d-tokens of image "
							<< i << " extracted";
					}
					count++;
					//------ 提取d-token...
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
	// Voting 票决算法显示图像...
	g.visualize();
	// 提取一致性假设...
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
	//------ 初步描述DMatch数据结构
	if ((0 == grid_matches_r.size()) && (0 == grid_matches_b.size()) && (0 == grid_matches_b.size()))
	{
		return false;
	}

	//---------- R G B Single Channel Image Data Struct Merge ---------//
	// ============ Keypoint Merge: v_r 、v_g 、v_b ==================//
	// ============ DMatch Merge: grid_matches_r、 grid_matches_g、 grid_matches_b =================//
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






