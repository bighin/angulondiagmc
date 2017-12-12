#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "config.h"
#include "inih/ini.h"

int configuration_handler(void *user,const char *section,const char *name,const char *value)
{
	struct configuration_t *pconfig=(struct configuration_t *)(user);

#define MATCH(s,n) ((strcmp(section,s)==0)&&(strcmp(name,n)==0))

	if(MATCH("general","prefix"))
	{
		pconfig->prefix=strdup(value);
	}
	else if(MATCH("general","seedrng"))
	{
		if(!strcmp(value,"true"))
			pconfig->seedrng=true;
		else if(!strcmp(value,"false"))
			pconfig->seedrng=false;
		else
			return 0;
	}
	else if(MATCH("general","progressbar"))
	{
		if(!strcmp(value,"true"))
			pconfig->progressbar=true;
		else if(!strcmp(value,"false"))
			pconfig->progressbar=false;
		else
			return 0;
	}
	else if(MATCH("general","animate"))
	{
		if(!strcmp(value,"true"))
			pconfig->animate=true;
		else if(!strcmp(value,"false"))
			pconfig->animate=false;
		else
			return 0;
	}
	else if(MATCH("general","liveplot"))
	{
		if(!strcmp(value,"true"))
			pconfig->liveplot=true;
		else if(!strcmp(value,"false"))
			pconfig->liveplot=false;
		else
			return 0;
	}
	else if(MATCH("general","hashtable"))
	{
		if(!strcmp(value,"true"))
			pconfig->use_hashtable=true;
		else if(!strcmp(value,"false"))
			pconfig->use_hashtable=false;
		else
			return 0;
	}
	else if(MATCH("parameters","j"))
	{
		pconfig->j=atoi(value);
	}
	else if(MATCH("parameters","endtau"))
	{
		pconfig->endtau=atof(value);
	}
	else if(MATCH("parameters","chempot"))
	{
		pconfig->chempot=atof(value);
	}
	else if(MATCH("parameters","maxtau"))
	{
		pconfig->maxtau=atof(value);
	}
	else if(MATCH("parameters","maxorder"))
	{
		pconfig->maxorder=atoi(value);
	}
	else if(MATCH("potential","logn"))
	{
		pconfig->n=exp(atof(value));
	}
	else if(MATCH("sampling","bins"))
	{
		pconfig->bins=atoi(value);
	}
	else if(MATCH("sampling","thermalization"))
	{
		pconfig->thermalization=atoi(value);
	}

	else if(MATCH("sampling","iterations"))
	{
		pconfig->iterations=atoi(value);
	}
	else if(MATCH("sampling","timelimit"))
	{
		pconfig->timelimit=atof(value);
	}
	else if(MATCH("sampling","width"))
	{	
		pconfig->width=atof(value);
	}
	else if(MATCH("parallel","parallel"))
	{
		if(!strcmp(value,"true"))
			pconfig->parallel=true;
		else if(!strcmp(value,"false"))
			pconfig->parallel=false;
		else
			return 0;
	}
	else if(MATCH("parallel","nthreads"))
	{		
		if(!strcmp(value,"auto"))
		{
#if defined(_OPENMP)
			pconfig->nthreads=omp_get_max_threads();
#else
			pconfig->nthreads=1;
#endif
		}
		else
			pconfig->nthreads=atoi(value);
	}
	else
	{
		/* Unknown section/name, error */
		return 0;
	}

	return 1;
}

void load_config_defaults(struct configuration_t *config)
{
	config->prefix=strdup("default");
	config->seedrng=false;
	config->progressbar=true;
	config->animate=false;
	config->liveplot=false;

	config->j=1;
	config->endtau=1.0f;
	config->chempot=0.95f;
	config->maxtau=100.0f;
	config->maxorder=1024;
	
	config->n=exp(-10.0f);

	config->iterations=10000000;
	config->thermalization=config->iterations/100;
	config->timelimit=0.0f;
	config->bins=50;
	config->width=0.25;

	config->parallel=false;
	config->nthreads=1;
	
	config->ininame=NULL;
}

bool load_configuration(char *configfile,struct configuration_t *config)
{
	load_config_defaults(config);

	if(ini_parse(configfile,configuration_handler,config)<0)
	{
		fprintf(stderr,"Couldn't read or parse '%s'\n",configfile);
		return false;
	}

	if((config->animate==true)&&(config->progressbar==true))
	{
		fprintf(stderr,"The options 'animate' and 'progressbar' cannot be set at the same time!\n");
		return false;
	}

	if((config->animate==true)&&(config->parallel==true))
	{
		fprintf(stderr,"The options 'animate' and 'parallel' cannot be set at the same time!\n");
		return false;
	}

	if(config->j*(config->j+1)<=config->chempot)
	{
		fprintf(stderr,"Warning: the chemical potential should be set below j(j+1).\n");
	}

	if((config->parallel==true)&&(config->use_hashtable==true))
	{
		fprintf(stderr,"Warning: the hashtable is incompatible with the 'parallel' option. The hashtable will be switched off.\n");
		config->use_hashtable=false;
	}

	if(config->parallel==false)
		config->nthreads=1;

	config->ininame=strdup(configfile);

	fprintf(stderr,"Loaded '%s'\n",configfile);
	return true;
}

