#ifndef TRIANGULATOR
#define TRIANGULATOR
#include "const_parameters.h"

/*BEGIN DELAUNEY*/
////////////////////////////////////////
#define DELETED -2
#define le 0
#define re 1
////////////////////////////////////////
struct Freenode {
	struct Freenode *nextfree;
};
////////////////////////////////////////
struct Freelist {
	std::vector<char*> blocks;
	struct Freenode *head;
	int nodesize;
};
////////////////////////////////////////
struct Point {
	float x, y;
};
////////////////////////////////////////
struct Site { // structure used both for sites and for vertices
	struct Point	coord;
	int sitenbr;
	int refcnt;
};
////////////////////////////////////////
struct Edge {
	float a, b, c;
	struct Site *ep[2];
	struct Site *reg[2];
	int edgenbr;
};
////////////////////////////////////////
struct Halfedge {
	struct Halfedge *ELleft, *ELright;
	struct Edge *ELedge;
	int ELrefcnt;
	char ELpm;
	struct Site *vertex;
	float ystar;
	struct Halfedge *PQnext;
};

////////////////////////////////////////
struct I2 {
	int ia;
	int ib;
};
////////////////////////////////////////
struct I3 {
	int ia;
	int ib;
	int ic;
};
////////////////////////////////////////

class R2 {
public:
	R2() {}
	R2(int ia, int ib) :ia(ia), ib(ib) {}
	bool operator==(const R2&b) const {
		return ia == b.ia && ib == b.ib;
	}
	bool operator<(const R2&b) const {
		if (ia == b.ia) { return ib < b.ib; }
		return ia < b.ia;
	}
	int ia;
	int ib;
};

////////////////////////////////////////
//The following class Delauney2 uses/modifies code of Steve Fortune's.
//The original disclaimer for that code is:
/*
* The author of this software is Steven Fortune.  Copyright (c) 1994 by AT&T
* Bell Laboratories.
* Permission to use, copy, modify, and distribute this software for any
* purpose without fee is hereby granted, provided that this entire notice
* is included in all copies of any software which is or includes a copy
* or modification of this software and in all copies of the supporting
* documentation for such software.
* THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
* WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
* REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
* OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
*/
/////////////////////////////////////////////////

class Triangulator {
public:
	Triangulator(bool bOutEdges, bool bOutTriangles, V2*p, int nP,
		int triangulate = 1, int sorted = 0, int plot = 0, int debug = 0);
	virtual ~Triangulator();
	void freeinit(struct	Freelist *fl, int size);
	char *getfree(struct	Freelist *fl);
	void makefree(struct Freenode *curr, struct Freelist *fl);
	void releaseList(struct Freelist*fl);
	char *myalloc(unsigned n);
	void geominit();
	void plotinit();
	void iniSites(V2*p, int nP);
	struct Site *next();
	void PQinitialize();
	void PQinsert(struct Halfedge *he, struct Site *v, float 	offset);
	void PQdelete(struct Halfedge *he);
	int PQbucket(struct Halfedge *he);
	struct Halfedge *PQextractmin();

	int PQempty();
	struct Point PQ_min();
	void ELinitialize();
	void ELinsert(struct	Halfedge *lb, struct	Halfedge * newone);
	void ELdelete(struct Halfedge *he);
	struct Halfedge * ELright(struct Halfedge *he);
	struct Halfedge * ELleft(struct Halfedge *he);
	struct Halfedge * ELgethash(int b);
	struct Halfedge * ELleftbnd(struct Point *p);

	struct Halfedge *HEcreate(struct Edge *e, int pm);
	void out_site(struct Site *s);
	void out_vertex(struct Site *v);
	void out_bisector(struct Edge *e);
	void out_ep(struct Edge *e);
	void out_triple(struct Site *s1, struct Site *s2, struct Site * s3);
	void clip_line(struct Edge *e);
	void line(int x1, int y1, int x2, int y2);
	int right_of(struct Halfedge *el, struct Point *p);
	struct Site *leftreg(struct Halfedge *he);
	struct Site *rightreg(struct Halfedge *he);
	struct Edge *bisect(struct	Site *s1, struct Site*s2);
	struct Site *intersect(struct Halfedge *el1, struct Halfedge *el2, struct Point *p = NULL);
	float dist(struct Site *s, struct Site *t);
	void deref(struct	Site *v);
	void ref(struct Site *v);
	void makevertex(struct Site *v);
	void endpoint(struct Edge *e, int lr, struct Site *s);

	void voronoi();

	struct Site *sites;
	int nsites;
	int siteidx;
	int sqrt_nsites;
	int nvertices;
	struct Freelist sfl;
	struct Site *bottomsite;
	int total_alloc;
	int triangulate;
	int sorted;
	int plot;
	int debug;
	float xmin, xmax, ymin, ymax, deltax, deltay;
	struct Halfedge *PQhash;
	int PQhashsize;
	int PQcount;
	int PQmin;
	int ntry;
	int totalsearch;
	struct Freelist hfl;
	struct Halfedge *ELleftend, *ELrightend;
	int 	ELhashsize;
	struct Halfedge **ELhash;
	float pxmin, pxmax, pymin, pymax, cradius;
	int nedges;
	struct Freelist efl;
	std::vector<I3> final_triangles;
	fx_set<R2> final_edges;
	bool bOutEdges;
	bool bOutTriangles;
	int nMallocs;//for debugging
};

#endif