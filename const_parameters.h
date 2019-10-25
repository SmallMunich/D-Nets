#ifndef CONST_PARAMETERS
#define CONST_PARAMETERS

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/legacy/legacy.hpp"
//#include <opencv2\nonfree\nonfree.hpp>
//#include <opencv2\nonfree\features2d.hpp>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <set>

using namespace std;
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif


// Variables for specific node extractors...
// FAST 角点提取算子阈值...
// int FAST_thresh = 80;

// DENSE SAMPLING...
//int g_dense_spacing = 10;
//bool g_dense_bAddGaussianNoise = true;
//float g_dense_stdDev = 3.0f;

template <class REAL> class TV2;
typedef class TV2<float> V2;
typedef class TV2<double> DV2;

template <class REAL> class TV2 
{
public:
	TV2() { x = (REAL)0; y = (REAL)0; }
	TV2(REAL x, REAL y) : x(x), y(y) {};
	TV2(const V2& v) : x((REAL)v.x), y((REAL)v.y) {}
	TV2(const DV2& v) : x((REAL)v.x), y((REAL)v.y) {}
	inline TV2<REAL> operator-(const TV2<REAL> &v) const {
		return TV2(x - v.x, y - v.y);
	}
	inline REAL len() const { return (REAL)sqrt(x*x + y*y); }
	inline REAL operator!() {
		REAL p = atan2(this->y, this->x);
		if (p < 0) { p += 2.0*M_PI; }
		return p;
	}
	inline REAL angleTo(TV2 v) { return angleFromTo(!*this, !v); }
	inline REAL toOne() {
		REAL l = len();
		if (exactly_equal(l, (REAL)0)) { return 0; }
		x /= l;
		y /= l;
		return l;
	}
	inline TV2<REAL> asOne() {
		REAL l = len();
		if (exactly_equal(l, (REAL)0)) { return *this; }
		return TV2<REAL>(x / l, y / l);
	}
	inline REAL operator*(const TV2& v) const {
		return x*v.x + y*v.y;
	}
	inline TV2<REAL> operator~() const {
		return TV2((-1)*y, x);
	}
	inline TV2<REAL> operator+(const TV2<REAL> &v) const {
		return TV2(x + v.x, y + v.y);
	}
	REAL x;
	REAL y;
};

/*****************************************************************************/
struct TwoV2 
{
	V2 a;
	V2 b;
};
/*****************************************************************************/
template <class REAL> class TV3
{
public:
	TV3() { x = y = z = 0; }
	TV3(REAL x, REAL y, REAL z) : x(x), y(y), z(z) {}
	TV3(const TV3& v) : x(v.x), y(v.y), z(v.z) {}
	TV3(const TV2<REAL>&s) : x(s.x), y(s.y), z(0) {}

	//inline TV3<REAL> operator-(const TV3&v) const;
	//inline TV3<REAL> operator+(const TV3&v) const;
	//inline TV3<REAL> operator*(REAL m) const;
	//inline TV3<REAL> operator/(REAL d) const;
	//inline TV3<REAL> operator-()const;
	inline REAL operator*(const TV3&v) const { return x*v.x + y*v.y + z*v.z; }
	//inline TV3<REAL> operator |(const TV3&v);
	//inline void operator /=(REAL f);
	inline void operator*=(REAL f) { x *= f; y *= f; z *= f; };
	//inline void operator +=(const TV3&v);
	//inline void operator -=(const TV3&v);
	//inline bool operator ==(const TV3&v) const;
	//inline bool operator !=(const TV3&v) const;
	//inline TV3<REAL>& operator=(const TV2<REAL> &v){x=v.x;y=v.y;z=0;return *this;};
	//inline operator TV2<REAL>() const{return TV2<REAL>(x,y);};
	REAL x, y, z;
};

typedef  class TV3<float> V3;
typedef  class TV3<double> DV3;
/*****************************************************************************/
template <class REAL> class TAT2 
{
public:
	TAT2() {};
	//////////////////////////////////////////
	TAT2(const TAT2<REAL>& s) 
	{
		a[0] = s.a[0];
		a[1] = s.a[1];
		a[2] = s.a[2];
		a[3] = s.a[3];
		t = s.t;
	};
	//////////////////////////////////////////
	~TAT2() {};
	//////////////////////////////////////////
	TAT2<REAL>& operator=(const TAT2<REAL>&s) 
	{
		a[0] = s.a[0];
		a[1] = s.a[1];
		a[2] = s.a[2];
		a[3] = s.a[3];
		t = s.t;
		return *this;
	}
	//////////////////////////////////////////
	bool operator==(const TAT2<REAL>&other) const 
	{
		if (t.x != other.t.x) { return false; }
		if (t.y != other.t.y) { return false; }
		if (a[0] != other.a[0]) { return false; }
		if (a[1] != other.a[1]) { return false; }
		if (a[2] != other.a[2]) { return false; }
		if (a[3] != other.a[3]) { return false; }
		return true;
	}
	//////////////////////////////////////////
	bool defineFromTriangles(TV2<REAL> triFrom[3], TV2<REAL> triTo[3]) 
	{
		REAL ax = triFrom[0].x;
		REAL ay = triFrom[0].y;
		REAL bx = triFrom[1].x;
		REAL by = triFrom[1].y;
		REAL cx = triFrom[2].x;
		REAL cy = triFrom[2].y;
		REAL A = ax*(by - cy) + bx*(cy - ay) + cx*(ay - by);//TODO:optimize see a**
		if (A == 0) { return false; }
		REAL Ai = 1.0f / A;
		REAL axs = triTo[0].x;
		REAL ays = triTo[0].y;
		REAL bxs = triTo[1].x;
		REAL bys = triTo[1].y;
		REAL cxs = triTo[2].x;
		REAL cys = triTo[2].y;
		TV3<REAL> t_bary(bx*cy - cx*by, cx*ay - ax*cy, ax*by - ay*bx);
		TV3<REAL> Q1(by - cy, cy - ay, ay - by);//a**
		TV3<REAL> Q2(cx - bx, ax - cx, bx - ax);
		TV3<REAL> v1(axs, bxs, cxs);
		TV3<REAL> v2(ays, bys, cys);
		v1 *= Ai;
		v2 *= Ai;
		a[0] = v1*Q1;
		a[1] = v1*Q2;
		a[2] = v2*Q1;
		a[3] = v2*Q2;
		t.x = v1*t_bary;
		t.y = v2*t_bary;
		return true;
	}
	//////////////////////////////////////////
	bool defineFromTriangles(TV2<REAL>& from0,
		TV2<REAL>& from1,
		TV2<REAL>& from2,
		TV2<REAL>& to0,
		TV2<REAL>& to1,
		TV2<REAL>& to2) 
	{
		REAL A = from0.x*(from1.y - from2.y) +
			from1.x*(from2.y - from0.y) +
			from2.x*(from0.y - from1.y);  //TODO:optimize see a**
		if (A == 0) { return false; }
		REAL Ai = 1.0f / A;
		TV3<REAL> t_bary(from1.x*from2.y - from2.x*from1.y,
			from2.x*from0.y - from0.x*from2.y,
			from0.x*from1.y - from0.y*from1.x);
		TV3<REAL> Q1(from1.y - from2.y, from2.y - from0.y, from0.y - from1.y);//a**
		TV3<REAL> Q2(from2.x - from1.x, from0.x - from2.x, from1.x - from0.x);
		TV3<REAL> v1(to0.x, to1.x, to2.x);
		TV3<REAL> v2(to0.y, to1.y, to2.y);
		v1 *= Ai;
		v2 *= Ai;
		a[0] = v1*Q1;
		a[1] = v1*Q2;
		a[2] = v2*Q1;
		a[3] = v2*Q2;
		REAL det = a[0] * a[3] - a[1] * a[2];
		if (det == 0) { return false; }
		t.x = v1*t_bary;
		t.y = v2*t_bary;
		return true;
	}
	//////////////////////////////////////////
	void defineAsIdentity() {
		a[0] = 1.0f;
		a[1] = 0.0f;
		a[2] = 0.0f;
		a[3] = 1.0f;
		t.x = 0;
		t.y = 0;
	}
	//////////////////////////////////////////
	TV2<REAL> out(TV2<REAL> p) {
		TV2<REAL> ps;
		ps.x = a[0] * p.x + a[1] * p.y + t.x;
		ps.y = a[2] * p.x + a[3] * p.y + t.y;
		return ps;
	}
	//////////////////////////////////////////
	TV2<REAL> vout(TV2<REAL> v) {
		TV2<REAL> ps;
		ps.x = a[0] * v.x + a[1] * v.y;
		ps.y = a[2] * v.x + a[3] * v.y;
		return ps;
	}
	//////////////////////////////////////////
	bool calcInverseTransformation(TAT2<REAL>& at) {
		REAL det = a[1] * a[2] - a[0] * a[3];
		if (det == 0) { return false; }
		at.a[0] = -a[3] / det;
		at.a[1] = a[1] / det;
		at.a[2] = a[2] / det;
		at.a[3] = -a[0] / det;
		at.t.x = -at.a[0] * t.x - at.a[1] * t.y;
		at.t.y = -at.a[2] * t.x - at.a[3] * t.y;
		return true;

	}
	//////////////////////////////////////////
	REAL phiOut(REAL phi) {
		TV2<REAL> v(cos(phi), sin(phi));
		TV2<REAL> vt = vout(v);
		float phit = !vt;
		return phit;
	}
	/*   //////////////////////////////////////////
	TV2phi<REAL> v2phiOut(TV2phi<REAL> p)
	{
	TV2phi<REAL> v2phit;
	TV2<REAL> pt=out(TV2<REAL>(p.x,p.y));
	float phiout=phiOut(p.phi);
	return TV2phi<REAL>(pt.x,pt.y,phiout);
	}*/
	//////////////////////////////////////////
	REAL calcError(TAT2<REAL>& c) {
		TV2<REAL> ref[3];
		ref[0].x = 0.0f;
		ref[0].y = 0.0f;
		ref[1].x = 1.0f;
		ref[1].y = 0.0f;
		ref[2].x = 0.0f;
		ref[2].y = 1.0f;
		TV2<REAL> t1[3];
		TV2<REAL> t2[3];
		REAL sumDist = 0;
		for (int i = 0; i<3; i++) {
			t1[i] = out(ref[i]);
			t2[i] = c.out(ref[i]);
			TV2<REAL> v = t1[i] - t2[i];
			REAL l = v.len();
			sumDist += l;
		}
		return sumDist;
	}
	//////////////////////////////////////////
	void defineAsRigidTransformation(REAL phi, TV2<REAL> nt) {
		t.x = nt.x;
		t.y = nt.y;
		REAL c = cos(phi);
		REAL s = sin(phi);
		a[0] = c;
		a[2] = s;
		a[1] = -s;
		a[3] = c;
	}
	//////////////////////////////////////////
	void defineAsScaledRigidTransformation(REAL scale, REAL phi, TV2<REAL> nt) {
		t.x = nt.x;
		t.y = nt.y;
		REAL c = scale*cosf(phi);
		REAL s = scale*sinf(phi);
		a[0] = c;
		a[2] = s;
		a[1] = -s;
		a[3] = c;
	}
	//////////////////////////////////////////
	bool defineAsScaledRigidTransformation(TV2<REAL> fromA,
		TV2<REAL> fromB,
		TV2<REAL> toA,
		TV2<REAL> toB) {
		TV2<REAL> vSource = fromB - fromA;
		TV2<REAL> vTarget = toB - toA;
		REAL lSource = vSource.len();
		REAL lTarget = vTarget.len();
		if (lSource == 0) { return false; }
		REAL scale = lTarget / lSource;
		defineAsScaledRigidTransformation(scale,
			angleFromTo(!vSource, !vTarget),
			TV2<REAL>());
		t = toA - out(fromA);
		return true;
	}
	//////////////////////////////////////////
	/* bool defineAsRigidTransformationAround(TV2<REAL> center, REAL dphi, TV2<REAL> dt)
	{
	TV2<REAL> p[3],ps[3];
	p[0]=TV2<REAL>(0.0f,0.0f);
	p[1]=TV2<REAL>(1.0f,0.0f);
	p[2]=TV2<REAL>(0.0f,1.0f);
	for (int i=0;i<3;i++){
	TV2<REAL> v=p[i]-center;
	v.rotate(dphi);
	ps[i]=center+v+dt;
	}
	return defineFromTriangles(p,ps);
	}*/
	//////////////////////////////////////////
	// Return the absolute delta angle in [rad] by which the
	// absolute angle between the column vectors of the a matrix differs
	// from 90 degrees.
	REAL getAngularDeviationFromBeingOrthogonal() {
		TV2<REAL> u(a[0], a[2]);
		TV2<REAL> v(a[1], a[3]);
		REAL dphi_from_ortho = fabs(angleFromTo((REAL)(M_PI / 2.0),
			(REAL)fabs(u.angleTo(v))));
		return dphi_from_ortho;
	}
	//////////////////////////////////////////
	REAL getRigidity() {
		TV2<REAL> va(a[0], a[2]);
		TV2<REAL> vb(a[1], a[3]);
		va.toOne();
		vb.toOne();
		REAL rigidity = va*vb;
		return 1.0 - fabs(rigidity);
	}
	//////////////////////////////////////////
	void getLengthOfAxes(REAL&lu, REAL&lv) {
		TV2<REAL> u(a[0], a[2]);
		TV2<REAL> v(a[1], a[3]);
		lu = u.len();
		lv = v.len();
	}
	//////////////////////////////////////////
	REAL a[4];
	TV2<REAL> t;
};
//////////////////////////////////////////////////////////////////////////////
typedef class TAT2<float> AT2;
typedef class TAT2<double> DAT2;  //(D)ouble (A)ffine (T)ransformation (2)d

/*****************************************************************************/
template <class REAL> inline REAL rnd() 
{
	return (REAL)rand() / (REAL)RAND_MAX;
}

/*****************************************************************************/
template <class REAL> REAL rndNormal(REAL mean, REAL var) 
{
	int l = 0;
	REAL s, v1;
	do {
		REAL u1 = rnd<REAL>();
		REAL u2 = rnd<REAL>();
		v1 = 2 * u1 - 1;
		REAL v2 = 2 * u2 - 1;
		s = v1*v1 + v2*v2;
		if (s<1) { l = 1; }
	} while (!l);
	REAL z = sqrt(-2 * log(s) / s)*v1;
	REAL x = z*sqrt(var) + mean;
	return x;
}

template <class DataType> class fx_set : public std::set<DataType> 
{
public:
	bool contains(DataType d) 
	{
		return std::find(this->begin(), this->end(), d) != this->end();
	}
	void remove(DataType d) 
	{
		erase(std::find(this->begin(), this->end(), d));
	}
};


template<class REAL> void swapvalues(REAL&a, REAL&b)
{
	REAL c = a;
	a = b;
	b = c;
}

/*****************************************************************************/
template <class REAL> inline REAL TAbs(REAL a) {
	return (a < 0.0) ? -a : a;
};

/*****************************************************************************/
template <class REAL> inline REAL wrapAngle(REAL phiFrom, REAL phiTo, REAL len)
{
	REAL a = phiTo - phiFrom;
	REAL b = len - a;
	if (TAbs(a) <= TAbs(b)) {
		return a;
	}
	else {
		return -b;
	}
}
/*****************************************************************************/
template <class REAL> inline void NormalizePhi(REAL& phi) {
	if (phi < 0.0) { phi += ((REAL)2.0*(REAL)M_PI); }
	if (phi >((REAL)2.0*(REAL)M_PI)) { phi -= ((REAL)2.0*(REAL)M_PI); }
}

/*****************************************************************************/
template <class REAL> inline REAL angleFromTo(REAL phiFrom, REAL phiTo)
{
	NormalizePhi(phiFrom);
	NormalizePhi(phiTo);
	REAL sig = 1;
	if (phiFrom>phiTo) {
		sig = -1;
		swapvalues(phiFrom, phiTo);
	}
	REAL dphi = wrapAngle(phiFrom, phiTo, ((REAL)2.0*(REAL)M_PI));
	return dphi*sig;
}


/*****************************************************************************/
template <class REAL> inline bool exactly_equal(REAL a, REAL b) {
	return !(a<b) && !(b<a);
}


#endif