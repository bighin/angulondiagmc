#ifndef __DIAGRAMS_H__
#define __DIAGRAMS_H__

/*
	The following structure defines an arc, i.e. a phonon line
*/

struct arc_t
{
	double k;
	int lambda,mu;

	int startmidpoint,endmidpoint;
	double starttau,endtau;
};

/*
	A free propagator
*/

struct g0_t
{
	int j,m;

	int startmidpoint,endmidpoint;
	double starttau,endtau;
	
	/*
		The first and last line in a diagram have immutable j and m.
	*/
};

/*
	A structure vertex_info_t is defined at every vertex, containing
	a pointer to the left free line, the right free line and the phonon line.
*/

struct vertex_info_t
{
	struct g0_t *left,*right;
	struct arc_t *phononline;
};

int vertex_get_j1(struct vertex_info_t *vif);
int vertex_get_m1(struct vertex_info_t *vif);
int vertex_get_j2(struct vertex_info_t *vif);
int vertex_get_m2(struct vertex_info_t *vif);
int vertex_get_lambda(struct vertex_info_t *vif);
int vertex_get_mu(struct vertex_info_t *vif);

/*
	This structure defines a diagram
*/

struct diagram_t
{
	double mintau,maxtau,chempot;

	struct vlist_t *phonons;
	struct vlist_t *midpoints;
	struct vlist_t *free_propagators;
	struct vlist_t *vertices;
};

struct diagram_t *init_diagram(double maxtau,int j,int m,double chempot);
void fini_diagram(struct diagram_t *dgr);

struct arc_t *get_phonon_line(struct diagram_t *dgr,int c);
double get_midpoint(struct diagram_t *dgr,int c);
struct g0_t *get_free_propagator(struct diagram_t *dgr,int c);
struct vertex_info_t *get_vertex(struct diagram_t *dgr,int c);

int get_nr_phonons(struct diagram_t *dgr);
int get_nr_midpoints(struct diagram_t *dgr);
int get_nr_free_propagators(struct diagram_t *dgr);
int get_nr_vertices(struct diagram_t *dgr);

struct g0_t *get_left_neighbour(struct diagram_t *dgr,int midpoint);
struct g0_t *get_right_neighbour(struct diagram_t *dgr,int midpoint);
void delete_left_neighbour(struct diagram_t *dgr,int midpoint);
void delete_right_neighbour(struct diagram_t *dgr,int midpoint);
struct arc_t *get_phonon_line_after_propagator(struct diagram_t *dgr,int c);

void diagram_check_consistency_of_times(struct diagram_t *dgr,double tau,int c);
void diagram_check_consistency(struct diagram_t *dgr);
double diagram_weight(struct diagram_t *dgr);
void print_diagram(struct diagram_t *dgr);

#endif //__DIAGRAMS_H__

