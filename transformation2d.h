#ifndef TRANSFORMATION2D
#define TRANSFORMATION2D
#include "const_parameters.h"

/*****************************************************************************/
/* BEGIN SOURCES FOR GEOMETRIIC VERIFICATION...*/
/*****************************************************************************/
/*****************************************************************************/
class Transformation2D 
{
public:
	Transformation2D();

	Transformation2D(const Transformation2D& s);

	Transformation2D& operator=(const Transformation2D&s) {
		nDefiningCorrespondences = s.nDefiningCorrespondences;
		bDefined = s.bDefined;
		bInversionValid = s.bInversionValid;
		return *this;
	}

	~Transformation2D();

	virtual void DefineAsIdentity(void) = 0;

	virtual bool defineFromPointCorrespondences(DV2* a, DV2* b, bool bCalcInverse = true) = 0;

	virtual void setTo(const Transformation2D* pT);

	virtual DV2 out(DV2 p) = 0;
	virtual DV2 in(DV2 p) = 0;
	virtual DV2 vout(DV2 v) = 0;
	virtual DV2 vin(DV2 v) = 0;
	virtual int getVariables(double*& x) = 0;
	virtual int getHStepPattern(double*& h) = 0;
	// For below, caller is responsible for cleanup, NULL if impossible
	virtual Transformation2D*createInverse() = 0;

protected:

	int nDefiningCorrespondences;
	bool bDefined;
	bool bInversionValid;
};

/*****************************************************************************/

class ATrans2 : public Transformation2D, public DAT2 
{
public:
	ATrans2(void);

	ATrans2(const ATrans2& s) : Transformation2D(s), DAT2(s) 
	{
		h[0] = s.h[0];
		h[1] = s.h[1];
		h[2] = s.h[2];
		h[3] = s.h[3];
		h[4] = s.h[4];
		h[5] = s.h[5];
	}

	ATrans2& operator=(const ATrans2& s) 
	{
		(*(Transformation2D*)this) = s;
		a[0] = s.a[0];
		a[1] = s.a[1];
		a[2] = s.a[2];
		a[3] = s.a[3];
		t = s.t;
		h[0] = s.h[0];
		h[1] = s.h[1];
		h[2] = s.h[2];
		h[3] = s.h[3];
		h[4] = s.h[4];
		h[5] = s.h[5];
		return *this;
	}
	virtual void setTo(const Transformation2D* pT);

	bool operator==(const ATrans2& other) const {
		return DAT2::operator==(other);
	}
	bool DefineFromTriangles(DV2 triFrom[3], DV2 triTo[3]);

	bool DefineFromTrianglesV2(V2 triFrom[3], V2 triTo[3]);

	bool defineFromPointCorrespondences(DV2* a, DV2* b, bool bCalcInverse = true);

	virtual int getVariables(double*& x);

	virtual void DefineAsIdentity(void);
	/***************************************************************/
	virtual DV2 out(DV2 p) {
		return DAT2::out(p);
	}
	DV2 vout(DV2 v) {
		return DAT2::vout(v);
	}
	virtual DV2 in(DV2 p) {
		printf("not implemented");
		return DV2();
	}
	DV2 vin(DV2 v) {
		printf("not implemented");
		return DV2();
	}
	/******************************************************************/
	bool CalcInverseTransformation(ATrans2& at);

	double phiOut(double phi);

	// DV2phi v2phiOut(DV2phi p);
	float CalcError(ATrans2& c);

	void DefineAsRigidTransformation(float phi, DV2 nt);

	void DefineAsScaledRigidTransformation(float scale, float phi, DV2 nt);

	bool DefineAsScaledRigidTransformation(DV2 fromA, DV2 fromB, DV2 toA, DV2 toB);
	/*
	void DefineAsRigidTransformationAround(DV2 center, float dphi, DV2 dt) {
	defineAsRigidTransformationAround(center,dphi,dt);
	bDefined = true;
	}
	*/
	float getAngularDeviationFromBeingOrthogonal();

	float GetRigidity(void);

	bool getLengthOfAxes(double& lu, double& lv);

	virtual int getHStepPattern(double*& h_);

	virtual Transformation2D*createInverse();

	double h[6];  // Do not add member variables before h
};

/*****************************************************************************/

class Homography : public Transformation2D 
{
public:
	Homography(const Homography& s) :
		tinyDet(s.tinyDet),
		thresh_allowedDefinitionDist(s.thresh_allowedDefinitionDist),
		bAvoidPlaneFlippers(s.bAvoidPlaneFlippers) {
		ini();
		set(s.x[0], s.x[1], s.x[2], s.x[3], s.x[4], s.x[5], s.x[6], s.x[7], s.x[8],
			s.bInversionValid);
	}

	Homography(double tinyDet = 1e-5, double thresh_allowedDefinitionDist = 0.01,
			   bool bAvoidPlaneFlippers = true) :
		tinyDet(tinyDet),
		thresh_allowedDefinitionDist(thresh_allowedDefinitionDist),
		bAvoidPlaneFlippers(bAvoidPlaneFlippers) {
			ini();
		}

	Homography(double h00, double h01, double h02,
		       double h10, double h11, double h12,
		       double h20, double h21, double h22,
		       bool bCalcInversion = true, double tinyDet = 1e-5,
		       double thresh_allowedDefinitionDist = 0.01,
		       bool bAvoidPlaneFlippers = true) :
		tinyDet(tinyDet),
		thresh_allowedDefinitionDist(thresh_allowedDefinitionDist),
		bAvoidPlaneFlippers(bAvoidPlaneFlippers)
	{
		ini();
		set(h00, h01, h02, h10, h11, h12, h20, h21, h22, bCalcInversion);
	}

	~Homography();

	virtual bool defineFromPointCorrespondences(DV2*morg, DV2*morg_, bool bCalcInversion = true);

	virtual void DefineAsIdentity(void);

	void ini();

	void set(double h00, double h01, double h02,
		double h10, double h11, double h12,
		double h20, double h21, double h22,
		bool bCalcInversion = true);

	void tryCalcInversion();

	bool flipsPlane(DV2* a, DV2* b);

	virtual int getVariables(double*& x_);

	virtual DV2 out(DV2 m);

	virtual DV2 in(DV2 m_);

	virtual int getHStepPattern(double*& h_);

	virtual void setTo(const Transformation2D* pT);

	virtual Transformation2D*createInverse();

	virtual DV2 vout(DV2 v);

	virtual DV2 vin(DV2 v);

	DV2 localOrigin;
	DV2 localOrigin_;
	double a[8 * 8];
	double ai[8 * 8];
	double b[8];
	double x[9];  // Intentionally 1 larger! 9 components of homography-matrix
	double hi[9];
	double p[3];
	double r[3];
	double hstep[8];
	class CvMat*pA;
	class CvMat*pAi;
	class CvMat*pX;
	class CvMat*pB;
	class CvMat*pH;
	class CvMat*pHi;
	class CvMat*p3;
	class CvMat*r3;
	double tinyDet;
	double thresh_allowedDefinitionDist;
	bool bAvoidPlaneFlippers;
};

#endif
