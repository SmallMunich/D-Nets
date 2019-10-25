#include "stdafx.h"
#include "triangulator.h"

Triangulator::Triangulator(bool bOutEdges, bool bOutTriangles, V2* p, int nP,
	int triangulate, int sorted, int plot, int debug)
	: triangulate(triangulate), sorted(sorted), plot(plot),
	debug(debug), bOutEdges(bOutEdges), bOutTriangles(bOutTriangles) {
	nMallocs = 0;
	PQhash = NULL;
	ELhash = NULL;
	total_alloc = 0;
	sites = NULL;
	freeinit(&sfl, sizeof(*sites));
	iniSites(p, nP);
	siteidx = 0;
	geominit();
	if (plot) { plotinit(); }
	voronoi();
}
/////////////////////////////////////////////////
Triangulator::~Triangulator() {
	if (sites) {
		free(sites);
		sites = NULL;
	}
	if (ELhash) {
		free(ELhash);
		ELhash = NULL;
	}
	if (PQhash) {
		free(PQhash);
		PQhash = NULL;
	}
	releaseList(&sfl);
	releaseList(&efl);
	releaseList(&hfl);
}
/////////////////////////////////////////////////
void Triangulator::freeinit(struct	Freelist *fl, int size) {
	fl->head = (struct Freenode *) NULL;
	fl->nodesize = size;
}
/////////////////////////////////////////////////
void Triangulator::releaseList(struct Freelist*fl) {
	for (int i = 0; i < fl->blocks.size(); i++) {
		free(fl->blocks[i]);
	}
	fl->blocks.clear();
	fl->head = NULL;
	fl->nodesize = 0;
}
/////////////////////////////////////////////////
char *Triangulator::getfree(struct	Freelist *fl)
{
	int i;
	struct Freenode *t;
	if (fl->head == (struct Freenode *) NULL) {
		t = (struct Freenode *) myalloc(sqrt_nsites * fl->nodesize);
		fl->blocks.push_back((char*)t);
		for (i = 0; i<sqrt_nsites; i++) {
			makefree((struct Freenode *)((char *)t + i*fl->nodesize), fl);
		}
	};
	t = fl->head;
	fl->head = fl->head->nextfree;
	return (char *)t;
}
/////////////////////////////////////////////////
void Triangulator::makefree(struct Freenode *curr, struct Freelist *fl) {
	curr->nextfree = fl->head;
	fl->head = curr;
}
/////////////////////////////////////////////////
char *Triangulator::myalloc(unsigned n) {
	char *t = (char*)malloc(n);
	if (!t) {
		fprintf(stderr, "Insufficient memory site %d (%d bytes in use)\n",
			siteidx, total_alloc);
		return NULL;
	};
	nMallocs++;
	total_alloc += n;
	return t;
}
/////////////////////////////////////////////////
void Triangulator::geominit() {
	struct Edge e;
	float sn;
	freeinit(&efl, sizeof e);
	nvertices = 0;
	nedges = 0;
	sn = nsites + 4;
	sqrt_nsites = sqrt(sn);
	deltay = ymax - ymin;
	deltax = xmax - xmin;
}
/////////////////////////////////////////////////
void Triangulator::plotinit() {
	float dx, dy, d;

	dy = ymax - ymin;
	dx = xmax - xmin;
	d = (dx > dy ? dx : dy) * 1.1;
	pxmin = xmin - (d - dx) / 2.0;
	pxmax = xmax + (d - dx) / 2.0;
	pymin = ymin - (d - dy) / 2.0;
	pymax = ymax + (d - dy) / 2.0;
	cradius = (pxmax - pxmin) / 350.0;
	//openpl();
	//range(pxmin, pymin, pxmax, pymax);
}
/////////////////////////////////////////////////
/* sort sites on y, then x, coord */
int scomp(const void*a, const void*b) {
	struct Point*s1 = (struct Point*)a;
	struct Point*s2 = (struct Point*)b;
	if (s1->y < s2->y) { return(-1); }
	if (s1->y > s2->y) { return(1); }
	if (s1->x < s2->x) { return(-1); }
	if (s1->x > s2->x) { return(1); }
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
void Triangulator::iniSites(V2*p, int nP) {
	sites = (struct Site*) myalloc(nP * sizeof(*sites));
	for (int i = 0; i<nP; i++) {
		sites[i].coord.x = p[i].x;
		sites[i].coord.y = p[i].y;
		sites[i].sitenbr = i;
		sites[i].refcnt = 0;
	}
	nsites = nP;

	qsort(sites, nsites, sizeof(*sites), scomp);
	xmin = sites[0].coord.x;
	xmax = sites[0].coord.x;
	for (int i = 1; i<nsites; i += 1) {
		if (sites[i].coord.x < xmin) { xmin = sites[i].coord.x; }
		if (sites[i].coord.x > xmax) { xmax = sites[i].coord.x; }
	}
	ymin = sites[0].coord.y;
	ymax = sites[nsites - 1].coord.y;
}
/////////////////////////////////////////////////
struct Site *Triangulator::next() {
	struct Site *s;
	if (siteidx < nsites) {
		s = &sites[siteidx];
		siteidx += 1;
		return s;
	}
	else {
		return (struct Site *)NULL;
	}
}
/////////////////////////////////////////////////
void Triangulator::PQinitialize() {
	int i;
	struct Point *s;
	PQcount = 0;
	PQmin = 0;
	PQhashsize = 4 * sqrt_nsites;
	PQhash = (struct Halfedge *) myalloc(PQhashsize * sizeof *PQhash);
	for (i = 0; i<PQhashsize; i++) {
		PQhash[i].PQnext = (struct Halfedge *)NULL;
	}
}
/////////////////////////////////////////////////
struct Halfedge *Triangulator::PQextractmin() {
	struct Halfedge *curr;
	curr = PQhash[PQmin].PQnext;
	PQhash[PQmin].PQnext = curr->PQnext;
	PQcount--;
	return curr;
}
/////////////////////////////////////////////////
void Triangulator::out_site(struct Site *s) {
	if (!triangulate & plot & !debug) {
		//circle (s->coord.x, s->coord.y, cradius);
		printf("circle visualization not implemented\n");
	}
	if (!triangulate & !plot & !debug) {
		printf("s %f %f\n", s->coord.x, s->coord.y);
	}
	if (debug) {
		printf("site (%d) at %f %f\n", s->sitenbr, s->coord.x, s->coord.y);
	}
}
/////////////////////////////////////////////////
void Triangulator::out_vertex(struct Site *v) {
	if (!triangulate & !plot & !debug) {
		printf("v %f %f\n", v->coord.x, v->coord.y);
	}
	if (debug) {
		printf("vertex(%d) at %f %f\n", v->sitenbr, v->coord.x, v->coord.y);
	}
}
/////////////////////////////////////////////////
void Triangulator::out_bisector(struct Edge *e)
{
	if (triangulate & plot & !debug) {
		line(e->reg[0]->coord.x, e->reg[0]->coord.y,
			e->reg[1]->coord.x, e->reg[1]->coord.y);
	}
	if (!triangulate & !plot & !debug) {
		printf("l %f %f %f", e->a, e->b, e->c);
	}
	if (debug) {
		printf("line(%d) %gx+%gy=%g, bisecting %d %d\n",
			e->edgenbr, e->a, e->b, e->c,
			e->reg[le]->sitenbr, e->reg[re]->sitenbr);
	}
}
/////////////////////////////////////////////////
void Triangulator::out_triple(struct Site *s1, struct Site *s2, struct Site * s3) {
	if (bOutTriangles) {
		I3 i3;
		i3.ia = s1->sitenbr;
		i3.ib = s2->sitenbr;
		i3.ic = s3->sitenbr;
		final_triangles.push_back(i3);
	}
	if (bOutEdges) {
		R2 e1(s1->sitenbr, s2->sitenbr);
		R2 e2(s2->sitenbr, s3->sitenbr);
		R2 e3(s3->sitenbr, s1->sitenbr);
		R2 e1m(s2->sitenbr, s1->sitenbr);
		R2 e2m(s3->sitenbr, s2->sitenbr);
		R2 e3m(s1->sitenbr, s3->sitenbr);
		if (!final_edges.contains(e1) && !final_edges.contains(e1m)) {
			final_edges.insert(e1);
		}
		if (!final_edges.contains(e2) && !final_edges.contains(e2m)) {
			final_edges.insert(e2);
		}
		if (!final_edges.contains(e3) && !final_edges.contains(e3m)) {
			final_edges.insert(e3);
		}
	}
	/*
	if(triangulate & !plot &!debug)
	printf("%d %d %d\n", s1->sitenbr, s2->sitenbr, s3->sitenbr);
	if(debug)
	printf("circle through left=%d right=%d bottom=%d\n",
	s1->sitenbr, s2->sitenbr, s3->sitenbr);
	*/
}
/////////////////////////////////////////////////
void Triangulator::out_ep(struct Edge *e) {
	/*
	if(!triangulate & plot) {
	clip_line(e);
	}
	if(!triangulate & !plot)
	{	printf("e %d", e->edgenbr);
	printf(" %d ", e->ep[le] != (struct Site *)NULL ? e->ep[le]->sitenbr : -1);
	printf("%d\n", e->ep[re] != (struct Site *)NULL ? e->ep[re]->sitenbr : -1);
	};
	*/
}
/////////////////////////////////////////////////
void Triangulator::clip_line(struct Edge *e) {
	struct Site *s1, *s2;
	float x1, x2, y1, y2;

	if (e->a == 1.0 && e->b >= 0.0) {
		s1 = e->ep[1];
		s2 = e->ep[0];
	}
	else {
		s1 = e->ep[0];
		s2 = e->ep[1];
	}

	if (e->a == 1.0) {
		y1 = pymin;
		if (s1 != (struct Site *)NULL && s1->coord.y > pymin) { y1 = s1->coord.y; }
		if (y1>pymax) { return; }
		x1 = e->c - e->b * y1;
		y2 = pymax;
		if (s2 != (struct Site *)NULL && s2->coord.y < pymax) { y2 = s2->coord.y; }
		if (y2<pymin) { return; }
		x2 = e->c - e->b * y2;

		if ((x1 > pxmax && x2 >pxmax) || (x1 < pxmin && x2<pxmin)) { return; }
		if (x1 > pxmax) { x1 = pxmax; y1 = (e->c - x1) / e->b; }
		if (x1 < pxmin) { x1 = pxmin; y1 = (e->c - x1) / e->b; }
		if (x2 > pxmax) { x2 = pxmax; y2 = (e->c - x2) / e->b; }
		if (x2 < pxmin) { x2 = pxmin; y2 = (e->c - x2) / e->b; }
	}
	else {
		x1 = pxmin;
		if (s1 && (s1->coord.x > pxmin)) { x1 = s1->coord.x; }
		if (x1 > pxmax) { return; }
		y1 = e->c - e->a * x1;
		x2 = pxmax;
		if (s2 && (s2->coord.x < pxmax)) { x2 = s2->coord.x; }
		if (x2 < pxmin) { return; }
		y2 = e->c - e->a * x2;
		if ((y1 > pymax && y2 > pymax) || (y1 < pymin && y2<pymin)) { return; }
		if (y1 > pymax) { y1 = pymax; x1 = (e->c - y1) / e->a; }
		if (y1 < pymin) { y1 = pymin; x1 = (e->c - y1) / e->a; }
		if (y2 > pymax) { y2 = pymax; x2 = (e->c - y2) / e->a; }
		if (y2 < pymin) { y2 = pymin; x2 = (e->c - y2) / e->a; }
	};
	line(x1, y1, x2, y2);
}
/////////////////////////////////////////////////
void Triangulator::line(int x1, int y1, int x2, int y2) {
	printf("visualization of line not implemented");
}
/////////////////////////////////////////////////
void Triangulator::ELinitialize() {
	int i;
	freeinit(&hfl, sizeof **ELhash);
	ELhashsize = 2 * sqrt_nsites;
	ELhash = (struct Halfedge **) myalloc(sizeof *ELhash * ELhashsize);
	for (i = 0; i < ELhashsize; i++) { ELhash[i] = (struct Halfedge *)NULL; }
	ELleftend = HEcreate((struct Edge *)NULL, 0);
	ELrightend = HEcreate((struct Edge *)NULL, 0);
	ELleftend->ELleft = (struct Halfedge *)NULL;
	ELleftend->ELright = ELrightend;
	ELrightend->ELleft = ELleftend;
	ELrightend->ELright = (struct Halfedge *)NULL;
	ELhash[0] = ELleftend;
	ELhash[ELhashsize - 1] = ELrightend;
}
/////////////////////////////////////////////////
void Triangulator::ELinsert(struct	Halfedge *lb, struct Halfedge *newone) {
	newone->ELleft = lb;
	newone->ELright = lb->ELright;
	lb->ELright->ELleft = newone;
	lb->ELright = newone;
}
/////////////////////////////////////////////////
struct Halfedge *Triangulator::HEcreate(struct Edge *e, int pm) {
	struct Halfedge *answer;
	answer = (struct Halfedge *) getfree(&hfl);
	answer->ELedge = e;
	answer->ELpm = pm;
	answer->PQnext = (struct Halfedge *) NULL;
	answer->vertex = (struct Site *) NULL;
	answer->ELrefcnt = 0;
	return answer;
}
/////////////////////////////////////////////////
/* This delete routine can't reclaim node, since pointers from hash
table may be present.   */
void Triangulator::ELdelete(struct Halfedge *he) {
	he->ELleft->ELright = he->ELright;
	he->ELright->ELleft = he->ELleft;
	he->ELedge = (struct Edge *)DELETED;
}
/////////////////////////////////////////////////
struct Halfedge	* Triangulator::ELright(struct Halfedge *he) {
	return he->ELright;
}
/////////////////////////////////////////////////
struct Halfedge	* Triangulator::ELleft(struct Halfedge *he) {
	return he->ELleft;
}
/////////////////////////////////////////////////
void Triangulator::PQinsert(struct Halfedge *he, struct Site *v, float offset) {
	struct Halfedge *last, *next;

	he->vertex = v;
	ref(v);
	he->ystar = v->coord.y + offset;
	last = &PQhash[PQbucket(he)];
	while ((next = last->PQnext) != (struct Halfedge *) NULL &&
		(he->ystar > next->ystar ||
		(he->ystar == next->ystar && v->coord.x > next->vertex->coord.x))) {
		last = next;
	}
	he->PQnext = last->PQnext;
	last->PQnext = he;
	PQcount++;
}
/////////////////////////////////////////////////
void Triangulator::PQdelete(struct Halfedge *he) {
	struct Halfedge *last;
	if (he->vertex) {
		last = &PQhash[PQbucket(he)];
		while (last->PQnext != he) { last = last->PQnext; }
		last->PQnext = he->PQnext;
		PQcount--;
		deref(he->vertex);
		he->vertex = (struct Site *) NULL;
	}
}
/////////////////////////////////////////////////
int Triangulator::PQbucket(struct Halfedge *he) {
	int bucket;
	bucket = (he->ystar - ymin) / deltay * PQhashsize;
	if (bucket < 0) { bucket = 0; }
	if (bucket >= PQhashsize) { bucket = PQhashsize - 1; }
	if (bucket < PQmin) { PQmin = bucket; }
	return bucket;
}
/////////////////////////////////////////////////
int Triangulator::PQempty() {
	return PQcount == 0;
}
/////////////////////////////////////////////////
struct Point Triangulator::PQ_min() {
	struct Point answer;
	while (!PQhash[PQmin].PQnext) { PQmin++; }
	answer.x = PQhash[PQmin].PQnext->vertex->coord.x;
	answer.y = PQhash[PQmin].PQnext->ystar;
	return answer;
}
/////////////////////////////////////////////////
/* Get entry from hash table, pruning any deleted nodes */
struct Halfedge *Triangulator::ELgethash(int b) {
	struct Halfedge *he;
	if (b < 0 || b >= ELhashsize) { return (struct Halfedge *) NULL; }
	he = ELhash[b];
	if (he == (struct Halfedge *) NULL ||
		he->ELedge != (struct Edge *) DELETED) {
		return he;
	}

	/* Hash table points to deleted half edge.  Patch as necessary. */
	ELhash[b] = (struct Halfedge *) NULL;
	if (--(he->ELrefcnt) == 0) { makefree((Freenode*)he, &hfl); }
	return (struct Halfedge *) NULL;
}
/////////////////////////////////////////////////
/* returns 1 if p is to right of halfedge e */
int Triangulator::right_of(struct Halfedge *el, struct Point *p) {
	struct Edge *e;
	struct Site *topsite;
	int right_of_site, above, fast;
	float dxp, dyp, dxs, t1, t2, t3, yl;

	e = el->ELedge;
	topsite = e->reg[1];
	right_of_site = p->x > topsite->coord.x;
	if (right_of_site && el->ELpm == le) { return(1); }
	if (!right_of_site && el->ELpm == re) { return (0); }

	if (e->a == 1.0) {
		dyp = p->y - topsite->coord.y;
		dxp = p->x - topsite->coord.x;
		fast = 0;
		if ((!right_of_site & e->b<0.0) | (right_of_site & e->b >= 0.0)) {
			above = dyp >= e->b*dxp;
			fast = above;
		}
		else {
			above = p->x + p->y*e->b > e->c;
			if (e->b < 0.0) { above = !above; }
			if (!above) { fast = 1; }
		}
		if (!fast) {
			dxs = topsite->coord.x - (e->reg[0])->coord.x;
			above = e->b * (dxp*dxp - dyp*dyp) <
				dxs*dyp*(1.0 + 2.0*dxp / dxs + e->b*e->b);
			if (e->b < 0.0) { above = !above; }
		}
	}
	else { // e->b==1.0
		yl = e->c - e->a*p->x;
		t1 = p->y - yl;
		t2 = p->x - topsite->coord.x;
		t3 = yl - topsite->coord.y;
		above = t1*t1 > t2*t2 + t3*t3;
	}
	return el->ELpm == le ? above : !above;
}
/////////////////////////////////////////////////
struct Halfedge *Triangulator::ELleftbnd(struct Point *p) {
	int i, bucket;
	struct Halfedge *he;

	/* Use hash table to get close to desired halfedge */
	bucket = (p->x - xmin) / deltax * ELhashsize;
	if (bucket<0) { bucket = 0; }
	if (bucket >= ELhashsize) { bucket = ELhashsize - 1; }
	he = ELgethash(bucket);
	if (he == (struct Halfedge *) NULL) {
		for (i = 1; 1; i++) {
			if ((he = ELgethash(bucket - i)) != (struct Halfedge *) NULL) { break; }
			if ((he = ELgethash(bucket + i)) != (struct Halfedge *) NULL) { break; }
		}
		totalsearch += i;
	}
	ntry++;
	/* Now search linear list of halfedges for the corect one */
	if (he == ELleftend || (he != ELrightend && right_of(he, p))) {
		do { he = he->ELright; } while (he != ELrightend && right_of(he, p));
		he = he->ELleft;
	}
	else {
		do { he = he->ELleft; } while (he != ELleftend && !right_of(he, p));
	}

	/* Update hash table and reference counts */
	if (bucket > 0 && bucket <ELhashsize - 1) {
		if (ELhash[bucket]) { --(ELhash[bucket]->ELrefcnt); }
		ELhash[bucket] = he;
		++(ELhash[bucket]->ELrefcnt);
	};
	return he;
}

/////////////////////////////////////////////////
struct Site *Triangulator::leftreg(struct Halfedge *he) {
	if (he->ELedge == (struct Edge *)NULL) { return(bottomsite); }
	return he->ELpm == le ? he->ELedge->reg[le] : he->ELedge->reg[re];
}
/////////////////////////////////////////////////
struct Site *Triangulator::rightreg(struct Halfedge *he) {
	if (he->ELedge == (struct Edge *)NULL) { return(bottomsite); }
	return he->ELpm == le ? he->ELedge->reg[re] : he->ELedge->reg[le];
}
/////////////////////////////////////////////////
void Triangulator::deref(struct Site *v) {
	--(v->refcnt);
	if (v->refcnt == 0) { makefree((Freenode*)v, &sfl); }
}
/////////////////////////////////////////////////
void Triangulator::ref(struct	Site *v) {
	++(v->refcnt);
}
/////////////////////////////////////////////////
struct Edge *Triangulator::bisect(struct Site *s1, struct Site*s2) {
	float dx, dy, adx, ady;
	struct Edge *newedge;
	newedge = (struct Edge *) getfree(&efl);

	newedge->reg[0] = s1;
	newedge->reg[1] = s2;
	ref(s1);
	ref(s2);
	newedge->ep[0] = (struct Site *) NULL;
	newedge->ep[1] = (struct Site *) NULL;

	dx = s2->coord.x - s1->coord.x;
	dy = s2->coord.y - s1->coord.y;
	adx = dx>0 ? dx : -dx;
	ady = dy>0 ? dy : -dy;
	newedge->c = s1->coord.x * dx + s1->coord.y * dy + (dx*dx + dy*dy)*0.5;
	if (adx>ady) {
		newedge->a = 1.0; newedge->b = dy / dx; newedge->c /= dx;
	}
	else {
		newedge->b = 1.0; newedge->a = dx / dy; newedge->c /= dy;
	}

	newedge->edgenbr = nedges;
	out_bisector(newedge);
	nedges++;
	return newedge;
}
/////////////////////////////////////////////////
struct Site *Triangulator::intersect(struct Halfedge *el1,
	struct Halfedge *el2, struct Point *p) {
	struct Edge *e1, *e2, *e;
	struct Halfedge *el;
	float d, xint, yint;
	int right_of_site;
	struct Site *v;

	e1 = el1->ELedge;
	e2 = el2->ELedge;
	if (!e1 || !e2) { return NULL; }
	if (e1->reg[1] == e2->reg[1]) { return NULL; }

	d = e1->a * e2->b - e1->b * e2->a;
	if (-1.0e-10 < d && d < 1.0e-10) { return NULL; }

	xint = (e1->c*e2->b - e2->c*e1->b) / d;
	yint = (e2->c*e1->a - e1->c*e2->a) / d;

	if ((e1->reg[1]->coord.y < e2->reg[1]->coord.y) ||
		(e1->reg[1]->coord.y == e2->reg[1]->coord.y &&
			e1->reg[1]->coord.x < e2->reg[1]->coord.x)) {
		el = el1;
		e = e1;
	}
	else {
		el = el2;
		e = e2;
	}
	right_of_site = xint >= e->reg[1]->coord.x;
	if ((right_of_site && el->ELpm == le) ||
		(!right_of_site && el->ELpm == re)) {
		return NULL;
	}

	v = (struct Site *) getfree(&sfl);
	v->refcnt = 0;
	v->coord.x = xint;
	v->coord.y = yint;
	return v;
}
/////////////////////////////////////////////////
float Triangulator::dist(struct Site *s, struct Site *t) {
	float dx, dy;
	dx = s->coord.x - t->coord.x;
	dy = s->coord.y - t->coord.y;
	return sqrt(dx*dx + dy*dy);
}
/////////////////////////////////////////////////
void Triangulator::makevertex(struct Site *v) {
	v->sitenbr = nvertices;
	++nvertices;
	out_vertex(v);
}
/////////////////////////////////////////////////
void Triangulator::endpoint(struct Edge *e, int lr, struct Site *s) {
	e->ep[lr] = s;
	ref(s);
	if (!(e->ep[re - lr])) { return; }
	out_ep(e);
	deref(e->reg[le]);
	deref(e->reg[re]);
	makefree((Freenode*)e, &efl);
}
/////////////////////////////////////////////////

void Triangulator::voronoi() {
	struct Site *newsite, *bot, *top, *temp, *p;
	struct Site *v;
	struct Point newintstar;
	int pm;
	struct Halfedge *lbnd, *rbnd, *llbnd, *rrbnd, *bisector;
	struct Edge *e;

	PQinitialize();
	bottomsite = next();
	out_site(bottomsite);
	ELinitialize();

	newsite = next();
	while (1) {
		if (!PQempty()) newintstar = PQ_min();
		if (newsite &&
			(PQempty() ||
				newsite->coord.y < newintstar.y ||
				(newsite->coord.y == newintstar.y &&
					newsite->coord.x < newintstar.x))) { // new site is smallest

			out_site(newsite);
			lbnd = ELleftbnd(&(newsite->coord));
			rbnd = ELright(lbnd);
			bot = rightreg(lbnd);
			e = bisect(bot, newsite);
			bisector = HEcreate(e, le);
			ELinsert(lbnd, bisector);
			if ((p = intersect(lbnd, bisector)) != NULL) {
				PQdelete(lbnd);
				PQinsert(lbnd, p, dist(p, newsite));
			}
			lbnd = bisector;
			bisector = HEcreate(e, re);
			ELinsert(lbnd, bisector);
			if ((p = intersect(bisector, rbnd)) != NULL) {
				PQinsert(bisector, p, dist(p, newsite));
			};

			newsite = next();
		}
		else if (!PQempty()) {  // intersection is smallest
			lbnd = PQextractmin();
			llbnd = ELleft(lbnd);
			rbnd = ELright(lbnd);
			rrbnd = ELright(rbnd);
			bot = leftreg(lbnd);
			top = rightreg(rbnd);

			out_triple(bot, top, rightreg(lbnd));
			v = lbnd->vertex;
			makevertex(v);
			endpoint(lbnd->ELedge, lbnd->ELpm, v);
			endpoint(rbnd->ELedge, rbnd->ELpm, v);
			ELdelete(lbnd);
			PQdelete(rbnd);
			ELdelete(rbnd);
			pm = le;
			if (bot->coord.y > top->coord.y) {
				temp = bot;
				bot = top;
				top = temp;
				pm = re;
			}
			e = bisect(bot, top);
			bisector = HEcreate(e, pm);
			ELinsert(llbnd, bisector);
			endpoint(e, re - pm, v);
			deref(v);
			if ((p = intersect(llbnd, bisector)) != NULL) {
				PQdelete(llbnd);
				PQinsert(llbnd, p, dist(p, bot));
			}
			if ((p = intersect(bisector, rrbnd)) != NULL) {
				PQinsert(bisector, p, dist(p, bot));
			};
		}
		else {
			break;
		}
	}

	for (lbnd = ELright(ELleftend); lbnd != ELrightend; lbnd = ELright(lbnd)) {
		e = lbnd->ELedge;
		out_ep(e);
	}
}

/*END DELAUNEY*/