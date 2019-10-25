#include "stdafx.h"
#include "transformation2d.h"

Transformation2D::Transformation2D()
{
	nDefiningCorrespondences = 0;
}

Transformation2D::Transformation2D(const Transformation2D& s) 
{
	nDefiningCorrespondences = s.nDefiningCorrespondences;
	bDefined = s.bDefined;
	bInversionValid = s.bInversionValid;
}

Transformation2D::~Transformation2D()
{ 
   
}

void Transformation2D::setTo(const Transformation2D* pT)
{
	// Makes this transformation being identical to *pT in any respect.
	bDefined = pT->bDefined;
	bInversionValid = pT->bInversionValid;
}

ATrans2::ATrans2(void) 
{
	bDefined = false;
	nDefiningCorrespondences = 3;
	h[0] = 0.1;
	h[1] = 0.1;
	h[2] = 0.1;
	h[3] = 0.1;
	h[4] = 0.1;
	h[5] = 0.1;
}

void ATrans2::setTo(const Transformation2D* pT)
{
	Transformation2D::setTo(pT);
	ATrans2&at = *(ATrans2*)pT;
	a[0] = at.a[0];
	a[1] = at.a[1];
	a[2] = at.a[2];
	a[3] = at.a[3];
	t.x = at.t.x;
	t.y = at.t.y;
}

bool ATrans2::DefineFromTriangles(DV2 triFrom[3], DV2 triTo[3]) 
{
	bDefined = defineFromTriangles(triFrom, triTo);
	return bDefined;
}

bool ATrans2::DefineFromTrianglesV2(V2 triFrom[3], V2 triTo[3])
{
	bDefined = false;
	float ax = triFrom[0].x;
	float ay = triFrom[0].y;
	float bx = triFrom[1].x;
	float by = triFrom[1].y;
	float cx = triFrom[2].x;
	float cy = triFrom[2].y;
	float A = ax*(by - cy) + bx*(cy - ay) + cx*(ay - by);//TODO:optimize see a**
	if (A == 0) { return false; };
	float Ai = 1.0f / A;

	float axs = triTo[0].x;
	float ays = triTo[0].y;
	float bxs = triTo[1].x;
	float bys = triTo[1].y;
	float cxs = triTo[2].x;
	float cys = triTo[2].y;

	DV3 t_bary(bx*cy - cx*by, cx*ay - ax*cy, ax*by - ay*bx);

	DV3 Q1(by - cy, cy - ay, ay - by);//a**
	DV3 Q2(cx - bx, ax - cx, bx - ax);
	DV3 v1(axs, bxs, cxs);
	DV3 v2(ays, bys, cys);
	v1 *= Ai;
	v2 *= Ai;

	a[0] = v1*Q1;
	a[1] = v1*Q2;
	a[2] = v2*Q1;
	a[3] = v2*Q2;

	t.x = v1*t_bary;
	t.y = v2*t_bary;
	bDefined = true;
	return true;
}

bool ATrans2::defineFromPointCorrespondences(DV2* a, DV2* b, bool bCalcInverse)
{
	return DefineFromTriangles(a, b);
}

int ATrans2::getVariables(double*& x)
{
	x = a;
	return 6;
}

void ATrans2::DefineAsIdentity(void) 
{
	a[0] = 1.0f;
	a[1] = 0.0f;
	a[2] = 0.0f;
	a[3] = 1.0f;
	t.x = 0;
	t.y = 0;
	bDefined = true;
}


bool ATrans2::CalcInverseTransformation(ATrans2& at)
{
	at.bDefined = calcInverseTransformation(at);
	return at.bDefined;
}

double ATrans2::phiOut(double phi)
{
	return DAT2::phiOut(phi);
}

// DV2phi v2phiOut(DV2phi p);
float ATrans2::CalcError(ATrans2& c)
{
	return DAT2::calcError(c);
}

void ATrans2::DefineAsRigidTransformation(float phi, DV2 nt)
{
	defineAsRigidTransformation(phi, nt);
	bDefined = true;
}

void ATrans2::DefineAsScaledRigidTransformation(float scale, float phi, DV2 nt)
{
	defineAsScaledRigidTransformation(scale, phi, nt);
	bDefined = true;
}

bool ATrans2::DefineAsScaledRigidTransformation(DV2 fromA, DV2 fromB, DV2 toA, DV2 toB)
{
	bDefined = defineAsScaledRigidTransformation(fromA, fromB, toA, toB);
	return bDefined;
}

float ATrans2::getAngularDeviationFromBeingOrthogonal()
{
	return DAT2::getAngularDeviationFromBeingOrthogonal();
}

float ATrans2::GetRigidity(void)
{
	if (!bDefined)
	{
		return 0.0f;
	}
	return getRigidity();
}

bool ATrans2::getLengthOfAxes(double& lu, double& lv)
{
	if (!bDefined) { return false; }
	getLengthOfAxes(lu, lv);
	return true;
}

int ATrans2::getHStepPattern(double*& h_) 
{
	h_ = h;
	return 6;
}

Transformation2D* ATrans2::createInverse() 
{
	if (!bDefined) 
	{
		printf("affine transformation that is to be inverted is undefined");
		return NULL;
	}
	ATrans2*pAI = new ATrans2();
	if (!CalcInverseTransformation(*pAI)) 
	{
		printf("failed to determine inverse transformation");
		return NULL;
	}
	return pAI;
}


Homography::~Homography() 
{
	delete pA;
	delete pAi;
	delete pX;
	delete pB;
	delete pH;
	delete pHi;
	delete p3;
	delete r3;
}


bool Homography::defineFromPointCorrespondences(DV2*morg, DV2*morg_, bool bCalcInversion)
{
	if (bAvoidPlaneFlippers&&flipsPlane(morg, morg_)) { return false; }
	// See homography.pdf for a documenation of the names of the variables
	// and a mathematical derivation of the underlying equations
	// to improve numerical precision define local coordinate systems in
	// the source and target space
	localOrigin.x = 0;
	localOrigin.y = 0;
	localOrigin_.x = 0;
	localOrigin_.y = 0;

	//localOrigin=0.25*(morg[0]+morg[1]+morg[2]+morg[3]);
	//localOrigin_=0.25*(morg_[0]+morg_[1]+morg_[2]+morg_[3]);
	DV2 m[4];
	DV2 m_[4];
	for (int i = 0; i<4; i++) {
		m[i] = morg[i] - localOrigin;
		m_[i] = morg_[i] - localOrigin_;
	}
	bInversionValid = false;
	bDefined = false;
#define u0 m[0].x
#define v0 m[0].y
#define u1 m[1].x
#define v1 m[1].y
#define u2 m[2].x
#define v2 m[2].y
#define u3 m[3].x
#define v3 m[3].y
#define u0_ m_[0].x
#define v0_ m_[0].y
#define u1_ m_[1].x
#define v1_ m_[1].y
#define u2_ m_[2].x
#define v2_ m_[2].y
#define u3_ m_[3].x
#define v3_ m_[3].y

	double* aa = a;
	*aa++ = 0; *aa++ = 0; *aa++ = 0; *aa++ = -u0; *aa++ = -v0; *aa++ = -1.0; *aa++ = u0*v0_; *aa++ = v0*v0_;
	*aa++ = u0; *aa++ = v0; *aa++ = 1.0; *aa++ = 0; *aa++ = 0; *aa++ = 0; *aa++ = -u0*u0_; *aa++ = -v0*u0_;
	*aa++ = 0; *aa++ = 0; *aa++ = 0; *aa++ = -u1; *aa++ = -v1; *aa++ = -1.0; *aa++ = u1*v1_; *aa++ = v1*v1_;
	*aa++ = u1; *aa++ = v1; *aa++ = 1.0; *aa++ = 0; *aa++ = 0; *aa++ = 0; *aa++ = -u1*u1_; *aa++ = -v1*u1_;
	*aa++ = 0; *aa++ = 0; *aa++ = 0; *aa++ = -u2; *aa++ = -v2; *aa++ = -1.0; *aa++ = u2*v2_; *aa++ = v2*v2_;
	*aa++ = u2; *aa++ = v2; *aa++ = 1.0; *aa++ = 0; *aa++ = 0; *aa++ = 0; *aa++ = -u2*u2_; *aa++ = -v2*u2_;
	*aa++ = 0; *aa++ = 0; *aa++ = 0; *aa++ = -u3; *aa++ = -v3; *aa++ = -1.0; *aa++ = u3*v3_; *aa++ = v3*v3_;
	*aa++ = u3; *aa++ = v3; *aa++ = 1.0; *aa++ = 0; *aa++ = 0; *aa++ = 0; *aa++ = -u3*u3_; *aa = -v3*u3_;

	b[0] = -v0_; b[1] = u0_; b[2] = -v1_; b[3] = u1_; b[4] = -v2_; b[5] = u2_; b[6] = -v3_; b[7] = u3_;

	double det = cvInvert(pA, pAi, CV_LU);
	if (fabs(det) <= tinyDet) {
		return false;
	}
	cvMatMul(pAi, pB, pX);
	bDefined = true;  // Set flag to allow mapping temporarily
					  // Verify whether the original points are mapped correctly
	double worstDist = 0;
	for (int i = 0; i<4; i++) {
		DV2 po = out(morg[i]);
		double curDist = (po - morg_[i]).len();
		if (curDist>worstDist) {
			worstDist = curDist;
		}
	}
	if (worstDist>thresh_allowedDefinitionDist) {
#ifndef NDEBUG
		printf("homography defined too imprecisely\n");
#endif
		bDefined = false;
		return false;
	}

	/*#ifndef NDEBUG
	double mem_check_is_b[8];
	CvMat cvCheck_is_b=cvMat(8,1,CV_64FC1,mem_check_is_b);
	cvMatMul(pA,pX,&cvCheck_is_b);
	#endif*/
	if (bCalcInversion) {
		tryCalcInversion();
	}
	bDefined = true;
	return true;
}


void Homography::DefineAsIdentity(void) 
{
	memset(x, 0, 8 * sizeof(*x));
	x[0] = 1.0;
	x[4] = 1.0;
	//x[8] is already 1.0
	bDefined = true;
	memset(hi, 0, 9 * sizeof(*hi));
	hi[0] = 1.0;
	hi[4] = 1.0;
	hi[8] = 1.0;
	bInversionValid = true;
}

void Homography::ini() 
{
	nDefiningCorrespondences = 4;
	bDefined = false;
	bInversionValid = false;

	// The new constructs here only allocate the headers not the data.
	// Allocation is done in this way to avoid the need to include
	// respective heades in the header of this file.
	pA = new CvMat();
	pAi = new CvMat();
	pX = new CvMat();
	pB = new CvMat();
	pH = new CvMat();
	pHi = new CvMat();
	p3 = new CvMat();
	r3 = new CvMat();

	// The initialization of the matrices is somehow unconventional, but is done so on purpsoe.
	// The tricky part of the code is that the memory block double x[9] (see the header file)
	// is used as underlying memory for two different matrices.
	// First only the first 8 of the 9 elemenets are used for the 8x1 vector pX which receives the result of a matrix multiplication as a solution of a system of linear equations.
	// Second the same data with an aditional element initialized to 1.0 is used as a 3x3 matrix to acutally project points.
	*pA = cvMat(8, 8, CV_64FC1, a);
	*pAi = cvMat(8, 8, CV_64FC1, ai);
	*pX = cvMat(8, 1, CV_64FC1, x);
	*pB = cvMat(8, 1, CV_64FC1, b);
	*pH = cvMat(3, 3, CV_64FC1, x);//duplicate use of x intentionally!
	*pHi = cvMat(3, 3, CV_64FC1, hi);
	*p3 = cvMat(3, 1, CV_64FC1, p);
	*r3 = cvMat(3, 1, CV_64FC1, r);
	x[8] = 1.0;
	p[2] = 1.0;
	r[2] = 1.0;

	hstep[0] = 0.01;
	hstep[1] = 0.01;
	hstep[2] = 0.01;
	hstep[3] = 0.01;
	hstep[4] = 0.01;
	hstep[5] = 0.01;
	hstep[6] = 0.000001;
	// @TODO: I believe that there is an error in how the hstp pattern is set,
	// because I think that the 9 components x of the homography when seen as
	// a matrix, are stored rowwise and the translation vector would then be
	// at indices [2] and [5]
	hstep[7] = 0.000001;
}

void Homography::set(double h00, double h01, double h02,
	                 double h10, double h11, double h12,
	                 double h20, double h21, double h22,
	                 bool bCalcInversion) 
{
	x[0] = h00;
	x[1] = h01;
	x[2] = h02;
	x[3] = h10;
	x[4] = h11;
	x[5] = h12;
	x[6] = h20;
	x[7] = h21;
	x[8] = h22;
	if (bCalcInversion) {
		tryCalcInversion();
	}
	bDefined = true;
}

void Homography::tryCalcInversion() 
{
	double detHi = cvInvert(pH, pHi, CV_LU);
	if (fabs(detHi)>tinyDet) 
	{
		bInversionValid = true;
	}
}

bool Homography::flipsPlane(DV2* a, DV2* b)
{
	for (int i = 0; i<4; i++)
	{
		DV2 pa = a[i];
		DV2 pb = b[i];
		for (int k = i + 1; k<4; k++) 
		{
			DV2 va = a[k] - pa;
			DV2 vb = b[k] - pb;
			DV2 na = ~va;
			DV2 nb = ~vb;
			for (int j = 0; j<4; j++)
			{
				if (j != i && j != k)
				{
					DV2 ta = a[j] - pa;
					DV2 tb = b[j] - pb;
					double da = ta*na;
					double db = tb*nb;
					if (da*db < 0) 
					{
						return true;
					}
				}
			}
		}
	}
	return false;
}

int Homography::getVariables(double*& x_)
{
	x_ = x;
	return 8;
}

DV2 Homography::out(DV2 m) 
{
	assert(bDefined);
	p[0] = m.x - localOrigin.x;
	p[1] = m.y - localOrigin.y;  //last homogeneous coordinate p[2]==1.0 (see constructor)
	cvMatMul(pH, p3, r3);
	double w = r[2];
	return localOrigin_ + DV2(r[0] / w, r[1] / w);
}

DV2 Homography::in(DV2 m_) 
{
	assert(bDefined);
	assert(bInversionValid);
	p[0] = m_.x - localOrigin_.x;
	p[1] = m_.y - localOrigin_.y;  //last homogeneous coordinate p[2]==1.0 (see constructor)
	cvMatMul(pHi, p3, r3);
	double w = r[2];
	return localOrigin + DV2(r[0] / w, r[1] / w);
}

int Homography::getHStepPattern(double*& h_) 
{
	h_ = hstep;
	return 8;
}

void Homography::setTo(const Transformation2D* pT)
{
	// Makes this transformation being identical to *pT in any respect.
	Transformation2D::setTo(pT);
	Homography&ht = *(Homography*)pT;
	tinyDet = ht.tinyDet;
	thresh_allowedDefinitionDist = ht.thresh_allowedDefinitionDist;
	bAvoidPlaneFlippers = ht.bAvoidPlaneFlippers;
	localOrigin = ht.localOrigin;
	localOrigin_ = ht.localOrigin_;
	memcpy(a, ht.a, sizeof(a));
	memcpy(ai, ht.ai, sizeof(ai));
	memcpy(b, ht.b, sizeof(b));
	memcpy(x, ht.x, sizeof(x));
	memcpy(hi, ht.hi, sizeof(hi));
	memcpy(hstep, ht.hstep, sizeof(hstep));
}

Transformation2D* Homography::createInverse()
{
	// Caller is responsible for cleanup, NULL if impossible
	if (!bDefined) 
	{
		printf("homography that is to be inverted is undefined\n");
		return NULL;
	}
	if (!bInversionValid) 
	{ 
		tryCalcInversion();
	}
	if (!bInversionValid) 
	{
		printf("failed to invert homography\n");
		return NULL;
	}
	Homography*pHI = new Homography(tinyDet, thresh_allowedDefinitionDist, bAvoidPlaneFlippers);
	bool bCalcInversion = true;  //(inversion of inversion)
	pHI->set(hi[0], hi[1], hi[2], hi[3], hi[4], hi[5], hi[6], hi[7], hi[8], true);
	return pHI;
}

DV2 Homography::vout(DV2 v) 
{ 
	printf("not implemented\n");
	return DV2();
}
DV2 Homography::vin(DV2 v) 
{ 
	printf("not implemented\n"); 
	return DV2();
}

