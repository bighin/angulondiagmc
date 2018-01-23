double calculate_propagators_and_vertices_ratio(struct diagram_t *numerator,int startmidpoint1,int endmidpoint1,
						struct diagram_t *denominator,int startmidpoint2,int endmidpoint2)
{
	assert(endmidpoint1>=startmidpoint1);
	assert(endmidpoint2>=startmidpoint2);

	assert((endmidpoint1-startmidpoint1)==(endmidpoint2-startmidpoint2));

	int length=endmidpoint1-startmidpoint1;
	double ret=0.0f;

	for(c=0;c<=length;c++)
	{

#warning This ratio can be optimized (taken in a single step)

		ret*=calculate_free_propagator_weight(numerator,get_left_neighbour(numerator,startmidpoint1+c));
		ret/=calculate_free_propagator_weight(denominator,get_left_neighbour(denominator,startmidpoint2+c));

		ret*=calculate_vertex_weight(numerator,startmidpoint1+c);
		ret/=calculate_vertex_weight(denominator,startmidpoint2+c);
	}

	assert((startmidpoint1+c)==endmidpoint1);
	assert((startmidpoint2+c)==endmidpoint2);

#warning This ratio can be optimized (taken in a single step)

	ret*=calculate_free_propagator_weight(numerator,get_right_neighbour(numerator,endmidpoint1));
	ret/=calculate_free_propagator_weight(denominator,get_right_neighbour(denominator,endmidpoint2));

	return ret;
}
