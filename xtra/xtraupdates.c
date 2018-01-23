#warning Probably truncated_phonon_dist() has numerical errors!

int update_shift_vertex(struct diagram_t *dgr,struct configuration_t *cfg)
{
	int target;
	double min,max,minlength,maxlength,oldlength,newlength,oldtau,newtau,weightratio,acceptance_ratio;

#ifndef NDEBUG
	double oldweight=diagram_weight(dgr);
#endif

	struct vertex_t *vertex;
	struct g0_t *left,*right;
	struct arc_t *arc;

	bool is_accepted,startshere;

	if(get_nr_vertices(dgr)<=0)
		return UPDATE_UNPHYSICAL;

	/*
		We select a random vertex...
	*/

	target=gsl_rng_uniform_int(dgr->rng_ctx,get_nr_vertices(dgr));
	vertex=get_vertex(dgr,target);

	/*
		...and we load the informations about the left and right neighbours.
	*/

	left=vertex->left;
	right=vertex->right;
	arc=vertex->phononline;

	assert((left!=NULL)&&(right!=NULL));
	assert(left->endtau==right->starttau);
	assert(left->endtau==get_midpoint(dgr,target));

	/*
		We load the position of the previous vertex, of the next vertex
		and of the selected vertex, which we will try to move.
	*/

	min=left->starttau;
	max=right->endtau;
	oldtau=get_midpoint(dgr,target);
	oldlength=arc->endtau-arc->starttau;

	assert((min<=oldtau)&&(oldtau<=max));
	assert((arc->starttau==get_midpoint(dgr,target))||(arc->endtau==get_midpoint(dgr,target)));

	weightratio=1.0f;
	weightratio/=calculate_free_propagator_weight(dgr,left);
	weightratio/=calculate_free_propagator_weight(dgr,right);
	weightratio/=calculate_arc_weight(dgr,arc);

	if(arc->starttau==get_midpoint(dgr,target))
	{
		/*
			In this case the phonon *starts* at the vertex we selected

 		                 _____....
			        |
			..._____|_____....
		*/

		minlength=arc->endtau-max;
		maxlength=arc->endtau-min;
		
		assert(minlength>=0.0f);
		assert(maxlength>=0.0f);
		assert(maxlength>minlength);

		/*
			We have to select oldtau in the interval [min,max]
			according a (truncated) phonon distribution function.
		*/

		newtau=arc->endtau-phonon_truncated_dist(dgr->rng_ctx,dgr->phonon_ctx,arc->lambda,minlength,maxlength);

		left->endtau=newtau;
		right->starttau=newtau;
		arc->starttau=newtau;
		*((double *)(vlist_get_element(dgr->midpoints,target)))=newtau;

		startshere=true;
	}
	else
	{
		/*
			In this case the phonon *ends* at the vertex we selected

		        ..._____
	        		|
			..._____|_____....
		*/

		minlength=min-arc->starttau;
		maxlength=max-arc->starttau;

		assert(minlength>=0.0f);
		assert(maxlength>=0.0f);
		assert(maxlength>minlength);

		newtau=arc->starttau+phonon_truncated_dist(dgr->rng_ctx,dgr->phonon_ctx,arc->lambda,minlength,maxlength);

		left->endtau=newtau;
		right->starttau=newtau;
		arc->endtau=newtau;
		*((double *)(vlist_get_element(dgr->midpoints,target)))=newtau;


		startshere=false;
	}

	weightratio*=calculate_free_propagator_weight(dgr,left);
	weightratio*=calculate_free_propagator_weight(dgr,right);
	weightratio*=calculate_arc_weight(dgr,arc);

	assert(weightratio>=0.0f);

#ifndef NDEBUG
	diagram_check_consistency(dgr);
#endif

	/*
		Finally we calculate the acceptance ratio for the update.
	*/

	newlength=arc->endtau-arc->starttau;

	acceptance_ratio=weightratio;
	acceptance_ratio*=phonon_truncated_pdf(dgr->phonon_ctx,arc->lambda,minlength,maxlength,oldlength);
	acceptance_ratio/=phonon_truncated_pdf(dgr->phonon_ctx,arc->lambda,minlength,maxlength,newlength);

	is_accepted=(gsl_rng_uniform(dgr->rng_ctx)<acceptance_ratio)?(true):(false);

	if(is_accepted==false)
	{
		left->endtau=oldtau;
		right->starttau=oldtau;
		*((double *)(vlist_get_element(dgr->midpoints,target)))=oldtau;
		
		if(startshere==true)
			arc->starttau=oldtau;
		else
			arc->endtau=oldtau;

#ifndef NDEBUG
		if(isfinite(oldweight)&&isfinite(diagram_weight(dgr)))
			assert(almost_same_float(oldweight,diagram_weight(dgr))==true);
#endif

		return UPDATE_REJECTED;
	}

	return UPDATE_ACCEPTED;
}

int update_swap_deltajs(struct diagram_t *dgr,struct configuration_t *cfg)
{
	int target,deltaj1,deltaj2;
	double acceptance_ratio,weightratio;
	bool is_accepted;

	struct arc_t *arc;

#ifndef NDEBUG
	double oldweight=diagram_weight(dgr);
#endif

	if(get_nr_phonons(dgr)<=0)
		return UPDATE_UNPHYSICAL;

	target=gsl_rng_uniform_int(dgr->rng_ctx,get_nr_phonons(dgr));
	arc=get_phonon_line(dgr,target);
	
	deltaj1=deltaj(dgr,arc->startmidpoint);
	deltaj2=deltaj(dgr,arc->endmidpoint);

	/*
		We swap the Delta_j's
	*/


	weightratio=1.0f/calculate_propagators_and_vertices(dgr,0,get_nr_vertices(dgr)-1);

	change_deltaj(dgr,arc->endmidpoint,deltaj1);
	change_deltaj(dgr,arc->startmidpoint,deltaj2);

	if(propagators_are_allowed(dgr)==false)
	{	
		change_deltaj(dgr,arc->endmidpoint,deltaj2);
		change_deltaj(dgr,arc->startmidpoint,deltaj1);

#ifndef NDEBUG
		if(isfinite(oldweight)&&isfinite(diagram_weight(dgr)))
			assert(almost_same_float(oldweight,diagram_weight(dgr))==true);
#endif

		return UPDATE_REJECTED;
	}

	weightratio*=calculate_propagators_and_vertices(dgr,0,get_nr_vertices(dgr)-1);

	/*
		Acceptance ratio is just the weight ratio.
	*/

	acceptance_ratio=fabs(weightratio);

	is_accepted=(gsl_rng_uniform(dgr->rng_ctx)<acceptance_ratio)?(true):(false);

	if(is_accepted==false)
	{
		change_deltaj(dgr,arc->endmidpoint,deltaj2);
		change_deltaj(dgr,arc->startmidpoint,deltaj1);

#ifndef NDEBUG
		if(isfinite(oldweight)&&isfinite(diagram_weight(dgr)))
			assert(almost_same_float(oldweight,diagram_weight(dgr))==true);
#endif

		return UPDATE_REJECTED;
	}

	if(weightratio<0.0f)
		dgr->sign*=-1;

	return UPDATE_ACCEPTED;
}

int update_change_mu(struct diagram_t *dgr,struct configuration_t *cfg)
{
	int c,target,oldmu,newmu;
	double acceptance_ratio,weightratio;
	bool is_accepted;

	struct arc_t *arc;

#ifndef NDEBUG
	double oldweight=diagram_weight(dgr);
#endif

	if(get_nr_phonons(dgr)<=0)
		return UPDATE_UNPHYSICAL;

#warning It would be better to select only arcs with \lambda != 0

	target=gsl_rng_uniform_int(dgr->rng_ctx,get_nr_phonons(dgr));
	arc=get_phonon_line(dgr,target);
	oldmu=arc->mu;
	newmu=gsl_rng_uniform_int(dgr->rng_ctx,1+2*arc->lambda)-arc->lambda;

	weightratio=1.0f/calculate_propagators_and_vertices(dgr,0,get_nr_vertices(dgr)-1);

	for(c=arc->startmidpoint+1;c<=arc->endmidpoint;c++)
		get_free_propagator(dgr,c)->m-=newmu-oldmu;

	if(propagators_are_allowed(dgr)==false)
	{
		arc->mu=oldmu;

		for(c=arc->startmidpoint+1;c<=arc->endmidpoint;c++)
			get_free_propagator(dgr,c)->m+=newmu-oldmu;

#ifndef NDEBUG
		if(isfinite(oldweight)&&isfinite(diagram_weight(dgr)))
			assert(almost_same_float(oldweight,diagram_weight(dgr))==true);
#endif

		return UPDATE_REJECTED;
	}

	weightratio*=calculate_propagators_and_vertices(dgr,0,get_nr_vertices(dgr)-1);

	/*
		Acceptance ratio is just the weight ratio.
	*/

	acceptance_ratio=fabs(weightratio);

	is_accepted=(gsl_rng_uniform(dgr->rng_ctx)<acceptance_ratio)?(true):(false);

	if(is_accepted==false)
	{
		arc->mu=oldmu;

		for(c=arc->startmidpoint+1;c<=arc->endmidpoint;c++)
			get_free_propagator(dgr,c)->m+=newmu-oldmu;

#ifndef NDEBUG
		if(isfinite(oldweight)&&isfinite(diagram_weight(dgr)))
			assert(almost_same_float(oldweight,diagram_weight(dgr))==true);
#endif

		return UPDATE_REJECTED;
	}

	if(weightratio<0.0f)
		dgr->sign*=-1;

	return UPDATE_ACCEPTED;
}
