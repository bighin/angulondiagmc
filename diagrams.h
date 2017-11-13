#ifndef __DIAGRAMS_H__
#define __DIAGRAMS_H__

#include <stdbool.h>
#include <gsl/gsl_rng.h>

/*
    GCC does not define M_PI in C11 mode
*/

#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif

/*
	The following structure defines an arc, i.e. a phonon line
*/

struct arc_t
{
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

	int arcs_over_me;
};

/*
	A structure vertex_t is defined at every vertex, containing
	a pointer to the left free line, the right free line and the phonon line.
*/

struct vertex_t
{
	struct g0_t *left,*right;
	struct arc_t *phononline;
};

int vertex_get_j1(struct vertex_t *vif);
int vertex_get_m1(struct vertex_t *vif);
int vertex_get_j2(struct vertex_t *vif);
int vertex_get_m2(struct vertex_t *vif);
int vertex_get_lambda(struct vertex_t *vif);
int vertex_get_mu(struct vertex_t *vif);

struct diagram_parameters_t
{
	/*
		General configuration parameters, to be specified
		when creating a diagram:
	
		- Angular momenta of the initial and final line
		- Initial length in imaginary time
		- Chemical potential
	*/

	int j;
	double endtau;
	double chempot;

	/*
		Parameters defining the interaction potential
	*/
	
	double n;
	
	/*
		The maximum allowed lenght
	*/
	
	double maxtau;
};

/*
	This structure defines a diagram
*/

struct diagram_t
{
	double mintau,endtau,chempot;

	struct vlist_t *phonons;
	struct vlist_t *midpoints;
	struct vlist_t *free_propagators;
	struct vlist_t *vertices;

	/*
		Parameters defining the interaction potential:

		- n is the adimensional density
	*/

	double n;

	/*
		A GSL random number generator context
	*/

	gsl_rng *rng_ctx;
	
	/*
		The context defining the phonon PDFs
	*/

	struct phonon_ctx_t *phonon_ctx;
};

struct diagram_t *init_diagram(struct diagram_parameters_t *dpars,bool verbose);
void fini_diagram(struct diagram_t *dgr);

struct arc_t *get_phonon_line(struct diagram_t *dgr,int c);
double get_midpoint(struct diagram_t *dgr,int c);
struct g0_t *get_free_propagator(struct diagram_t *dgr,int c);
struct vertex_t *get_vertex(struct diagram_t *dgr,int c);

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

bool is_last_propagator(struct diagram_t *dgr,struct g0_t *g0);
double calculate_free_propagator_weight(struct diagram_t *dgr,struct g0_t *g0);
double calculate_arc_weight(struct diagram_t *dgr,struct arc_t *arc);
double calculate_vertex_weight(struct diagram_t *dgr,int index);
double diagram_weight(struct diagram_t *dgr);

void diagram_copy(struct diagram_t *src,struct diagram_t *dst);
struct diagram_t *diagram_clone(struct diagram_t *src);

#define PRINT_TOPOLOGY		(0x01)
#define PRINT_PROPAGATORS	(0x02)
#define PRINT_DRYRUN		(0x04)
#define PRINT_INFO0		(0x08)
#define PRINT_ARCS_OVER_ME	(0x10)

int print_diagram(struct diagram_t *dgr,int flags);

bool check_triangle_condition_and_parity(struct diagram_t *dgr,struct vertex_t *thisvertex);
bool check_triangle_condition_and_parity_from_js(int j1,int j2,int j3);

#endif //__DIAGRAMS_H__
