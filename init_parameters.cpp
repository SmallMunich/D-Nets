#include "stdafx.h"
#include "init_parameters.h"

int nLayers = 8;
std::vector< pair<int, char*> > imageTypeStringMapping;

/*****************************************************************************/
void init()
{
	imageTypeStringMapping.push_back(pair<int, char*>(CV_8U, (char*)"CV_8U"));
	imageTypeStringMapping.push_back(pair<int, char*>(CV_8S, (char*)"CV_8S"));
	imageTypeStringMapping.push_back(pair<int, char*>(CV_16U, (char*)"CV_16U"));
	imageTypeStringMapping.push_back(pair<int, char*>(CV_16S, (char*)"CV_16S"));
	imageTypeStringMapping.push_back(pair<int, char*>(CV_32S, (char*)"CV_32S"));
	imageTypeStringMapping.push_back(pair<int, char*>(CV_32F, (char*)"CV_32F"));
	imageTypeStringMapping.push_back(pair<int, char*>(CV_64F, (char*)"CV_64F"));
}


/*
    readParam() �������ú���
	    arg :  resizeString[0/1]
*/
bool readParam(string& flag, string& param, bool& bFloatValid, float& value, char* arg)
{
	// remove optional preceeding '-'
	while (*arg == '-')
	{
		arg++;
	}
	int l = strlen(arg);
	if (l < 3)
	{
		return false;
	}
	// search for '='
	int ie = 1;
	while (ie<l && arg[ie] != '=')
	{
		ie++;
	}
	if (ie == l)
	{
		return false;
	}
	// Failed to find '='
	flag  = string(arg, ie);
	param = string(arg + ie + 1, l - ie);
	bFloatValid = (sscanf_s(param.c_str(), "%f", &value) == 1);
	return true;
}

/*****************************************************************************/
/*
readAndConvertImages() ��ȡ����ͬʱת��ͼ��Ϊ������
	parameters:
		 imgGray[2]  ��Ӧ��Main  img
		 argc :   // argc���ⲿ��������ĸ���
		 argv :   // ��Ÿ�����
				  ����argv[1] argv[3] ��ԭʼͼ���ļ�·��
					  argv[2] argv[4] �ǲ���resizeString[0] resizeString[1]����
*/
bool readAndConvertImages(cv::Mat imgGray[2], cv::Mat img[2]/*char** argv*/)
{
	bool bFixedWidth[2] = { true,true };
	int fixedWidth[2];     // �޸ĵĿ��
	float scaleFactor[2];  // �߶Ȳ���
	char* filename[2];     // ͼ��·��
	char* resizeString[2]; // 

	//filename[0]     = argv[1];  // ��ƥ��ͼ��1�ļ�·��
	//resizeString[0] = argv[2];
	//filename[1]     = argv[3];  // ��ƥ��ͼ��2�ļ�·��
	//resizeString[1] = argv[4];
	filename[0] = "../D-Nets/image/.....";
	filename[1] = "../D-Nets/image/.....";
	resizeString[0] = "";
	resizeString[1] = "";
	//cv::Mat img[2]; // ��ƥ��ͼ�����
	for (int i = 0; i<2; i++)
	{   // img[i]��ȡ����ԭʼͼ���ͨ������
		//img[i] = cv::imread(filename[i], 1); 
		//if (img[i].empty()) 
		//{   // ͼ���Ƿ��ȡ�ɹ�
		//	cout << "Failed to open image \"" << filename[i] << "\"." << endl;
		//	return false;
		//}
		// ȷ���Ҷ�ͼ������
		if (!ensureGrayImageWithDepth(imgGray[i], CV_32F, img[i], filename[i], true))
		{
			cout << "Failed to convert image \"" << filename[i]
				 << "\" to a gray image with float values." << endl;
			return false;
		}
		// �޸�ͼ���С...
		string flag;
		float value;
		bool bValueOk;
		string param;
		// readParam() ��ȡ�����ӿں���
		// �������: resizeString[0] ��resizeString[1]  �������� flag�� param�� bValueOk�� value
		if (readParam(flag, param, bValueOk, value, resizeString[i]) && (flag == "w" || flag == "s") && bValueOk)
		{
			if (flag == "w")
			{
				bFixedWidth[i] = true;
				fixedWidth[i] = (int)value;
			}
			else if (flag == "s")
			{
				bFixedWidth[i] = false;
				scaleFactor[i] = value;
			}
			cv::Size sz = imgGray[i].size();
			if (bFixedWidth[i])
			{
				scaleFactor[i] = (float)fixedWidth[i] / (float)sz.width;
			}
			if (scaleFactor[i] != 1.)
			{
				cv::Mat tmp;
				cv::resize(imgGray[i], tmp, cv::Size(0, 0), scaleFactor[i], scaleFactor[i], CV_INTER_AREA);
				tmp.copyTo(imgGray[i]);
			}
		}
	}
	return true;
}

/*****************************************************************************/

/*****************************************************************************/
void help()
{
	cout << "\
This program demonstrates D-Nets for image matching. The code implements  \n\
the exhausitve version (Clique D-Nets) of our CVPR2012 paper:             \n\
                                                                          \n\
              D-Nets: Beyond Patch-Based Image Descriptors                \n\
                                                                          \n\
              Felix von Hundelshausen       Rahul Sukthankar              \n\
              felix.v.hundelshausen@live.de rahuls@cs.cmu.edu             \n\
                                                                          \n\
 IEEE International Conference on Computer Vision and Pattern Recognition \n\
              June 18-20, 2012, Providence, Rhode Island, USA             \n\
                                                                          \n\
                                                                          \n\
The program matches two images using FAST interest points as nodes.       \n\
                                                                          \n\
We recommend running the program on a 64-bit architecture with 16 GB      \n\
memory. If the program is executed on a 32-bit architecture, or with      \n\
much less memory, set parameter s<=10 (for b=2) through the command       \n\
line parameters. However the difficult matching cases with a large        \n\
scale change and many interest points will require values up to s=13.     \n\
Those cases can only be sucessfully run on a 64-bit machine with enough   \n\
memory.                                                                   \n\
Command Line Parameters                                                   \n\
-----------------------                                                   \n\
The first 4 command line parameters need to be supplied always and in a   \n\
fixed order, e.g.,:                                                       \n\
                1           2          3          4                       \n\
           -------------- ----- --------------- -----                     \n\
   ./dnets boats/img1.pgm s=1.0 boats/img2.pgm  w=640                     \n\
                                                                          \n\
The first 4 command line parameters consist of two groups of 2 parameters.\n\
Within each group, the first parameter is a filename of an input image,   \n\
the second parameter determines how the respective image should be scaled \n\
initially.                                                                \n\
Conventional image formats (such as *.png, *.jpg, *.bmp) can be read.     \n\
But color images will be converted to grayscale images internally.        \n\
The initial scaling can be specified in two ways:                         \n\
    1. The first alternative,  e.g. \"s=0.5\" specified that the image    \n\
       is to be scaled to have a resulting scale of 50% of the original   \n\
       image (in this example). That is, to leave the original image      \n\
       unscaled \"s=1\" can be specified.                                 \n\
    2. The second alternative, e.g. \"w=640\" specifies that the image    \n\
       should be scaled, such that after scaling it, the width should be  \n\
       640 pixels (in this example).                                      \n\
                                                                          \n\
Some examples for valid calls with the first 4 parameters are:            \n\
   ./dnets boats/img1.png s=1 /boats/img2.png s=1                         \n\
   ./dnets boats/img1.png w=800 /boats/img2.png s=0.4                     \n\
   ./dnets boats/img1.png w=300 /boats/img2.png w=640                     \n\
That is, the flag \"s=\" or \"w=\" always relates to the preceeding image.\n\
The fist example will leave both images at their original scale.          \n\
After the first 4 parameters in fixed order, there can be further optional\n\
parameters with no required specific order.                               \n\
They are:                                                                 \n\
   sigma=[1.0]        Smooth the initial image with a Gaussian kernel with\n\
                      standard deviation 1.0.                             \n\
   kx=[FAST]          Choose FAST as node extractor                       \n\
     =SIFT            Choose SIFT as node extractor                       \n\
     =DENSE_SAMPLING  Choose DENSE_SAMPLING add node extractor             \n\
   ds_spacing=[10]    Spacing of points in dense sampling                 \n\
   gm=NONE            Disable geometric verification                      \n\
     =AFFINE          Affine geometric verification                       \n\
     =[HOMOGRAPHY]    Homography geometric verification                   \n\
   gv_dist=[16]       Distance for geometrice verification (if enabled)   \n\
   FAST_thresh=[80]   Use a threshold of 80 for FAST-keypoint extraction. \n\
   L=[8]              Create an image pyramid with 8 levels.              \n\
   q0=0.1             Encode a strip starting at 0.1% of the strip.       \n\
   q1=0.8             Encode strip up to 0.8% of the strip.               \n\
   nS=[9]             Divide a strip into nS sections                     \n\
   b=[2]              Use 2 bits to encode each section.                  \n\
   nL=[20]            Limit the size of the lists in the hash table to 20 \n\
                      pairings.                                           \n\
   om=matches.jpg     Create file \"matches.jpg\" showing the matches     \n\
                      between both images based on OpenCV's drawMatches   \n\
                      fuction.                                            \n\
   nM=40              Only extract the best 40 matches. If not supplied   \n\
                      all matches are extractd by default.                \n\
   vis=[LINES]        Visualize final correspondences using lines.       \n\
      =MESHES         Visualize the correspondences via meshes.           \n\
   wait=false         Do not wait for ESC at the end but exit immediately.\n\
                                                                          \n\
   -h                 print this command line info and exit.              \n\
   --h                print this command line info and exit.              \n\
   (defaults are indicated by [] but the braces are not part of the input)\n\
                                                                          \n\
   Calling dnets without any argmuents will print this help, too.\n";
}

/*****************************************************************************/
const char* depthCode2String(int depth)
{
	for (unsigned int i = 0; i<imageTypeStringMapping.size(); i++)
	{
		if (imageTypeStringMapping[i].first == depth)
		{
			return imageTypeStringMapping[i].second;
		}
	}
	return "unknown depth code";
}

/*****************************************************************************/
bool ensureGrayImageWithDepth(cv::Mat& imgGray, int target_depth, cv::Mat& img,
	                          char*image_name, bool bPrintFormatInfo)
{
	cv::Size sz = img.size(); // img.sizeͼ��Ĵ�С
	int nC = img.channels();  // ͼ���ͨ����
	int depth = img.depth();  // ͼ������
	if (bPrintFormatInfo)
	{
		cout << "input image \"" << (image_name ? image_name : "") << "\" has "
			<< nC << " channels of format " << depthCode2String(depth)
			<< " with " << sz.width << "x" << sz.height << " pixels." << endl;
	}
	cv::Mat grayButWithSourceDepth;
	switch (nC) // nCΪ����ͼ���ͨ����
	{
	case 1:
		img.copyTo(grayButWithSourceDepth);//ͨ����Ϊ1ֱ�ӿ���ͼ��������grayButWithSourceDepth����
		break;
	case 3: switch (depth) // ͨ����Ϊ3��ת��Ϊ�Ҷ�ͼ����grayButWithSourceDepth����  ��֧�����Ϊ8 16 32
	{
	case CV_8U: case CV_16U: case CV_32F:
		cvtColor(img, grayButWithSourceDepth, CV_BGR2GRAY);
		break;
	default:
		return false;
	}
			break;
			// ����4ͨ��ͼ��BGRA
	case 4: switch (depth) // ͨ����Ϊ4ʱ��ת��Ϊ�Ҷ�ͼ��洢�ھ���grayButWithSourceDepth
	{
	case CV_8U: case CV_16U: case CV_32F:
		cvtColor(img, grayButWithSourceDepth, CV_BGRA2GRAY);
		break;
	default:
		return false;
	}
			break;
	default:
		return false;
	}
	// inter_nC inter_depth �ֱ�Ϊ�Ҷ�ͼ�����grayButWithSourceDepth��ͨ���������
	int inter_nC = grayButWithSourceDepth.channels();
	int inter_depth = grayButWithSourceDepth.depth();
	// �о�������ǵ�ͨ���򷵻ش���
	if (inter_nC != 1)
	{
		cout << "only one channel expected after first step of conversion" << endl;
		return false;
	}
	// Ŀ��ͼ������ת��  ����target_depth ������grayButWithSourceDepthת����target_depth�趨��ͼ��depth  ͬʱ�����ݿ�����imgGray
	switch (target_depth)
	{
	case CV_8U:
		switch (inter_depth)
		{
		case CV_8U:
			grayButWithSourceDepth.copyTo(imgGray);
			break;
		case CV_16U:
			grayButWithSourceDepth.convertTo(imgGray, CV_8U, 255.0 / 65535.0);
			break;
		case CV_32F:
			grayButWithSourceDepth.convertTo(imgGray, CV_8U, 255.0);
			break;
		}
		break;
	case CV_32F:
		switch (inter_depth)
		{
		case CV_8U:
			grayButWithSourceDepth.convertTo(imgGray, CV_32F, 1.0 / 255.0);
			break;
		case CV_16U:
			grayButWithSourceDepth.convertTo(imgGray, CV_32F, 1.0 / 65535.0);
			break;
		case CV_32F:
			grayButWithSourceDepth.copyTo(imgGray);
			break;
		}
		break;
	}
	int final_nC = imgGray.channels(); // ͼ��ͨ����
	int final_depth = imgGray.depth(); // ͼ������
	cv::Size final_sz = imgGray.size();// ͼ��Ĵ�С
									   // ���յ�ͼ���趨��֤
									   // final_sz(��ͨ��ͼ��Ĵ�С)Ҫ��sz����ͼ��Ĵ�Сһ�� 
									   // final_nCҪΪ��ͨ��   final_depthҪ���趨ͼ�����һ��
	if (final_sz != sz || final_nC != 1 || final_depth != target_depth)
	{
		return false;
	}
	return true;
}


bool deal_init_paramers(const int& argc)
{
	if (argc < 6)
	{
		return false;
	}

	return true;
}
