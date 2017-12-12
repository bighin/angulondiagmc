#include <iostream>

#include <alps/accumulators.hpp>
#include <alps/mc/api.hpp>
#include <alps/mc/mcbase.hpp>
#include <alps/mc/stop_callback.hpp>

extern "C"
{
#include "../diagrams.h"
#include "../config.h"
#include "../mc.h"
#include "../physics.h"
}

/*
	The angulon diagrammatic Monte Carlo simulation class.

	We extend alps::mcbase, which is the base class of all Monte Carlo simulations.
*/

class angulon_diagmc : public alps::mcbase
{
	private:

	/*
		Global variables and statistics
	*/

	unsigned long sweeps,total_sweeps,thermalization_sweeps;
	size_t nloop,nhist,nbins;
	double binwidth;

	/*
		Acceptance/rejection statistics for each update
	*/

	std::vector<int> accepted;
	std::vector<int> rejected;

	/*
		Number of physical/unphysical updates, order by order
	*/

	std::vector<long int> physical_updates;
	std::vector<long int> unphysical_updates;

	/*
		The diagram state
	*/

	diagram_t *dgr;
	configuration_t config;

	typedef int (*update_function)(struct diagram_t *dgr,struct configuration_t *cfg);

	std::vector<std::pair<update_function,const char *>> updates;

	/*
		The values calculated or updated at the current MC step
	*/

	std::vector<double> hist_g;
	std::vector<double> hist_g0;
	std::vector<double> hist_g1;

	std::vector<double> order_frequencies;
	
	double intG0,hist_order;

	public:

	angulon_diagmc(parameters_type const &params, std::size_t seed_offset);
	~angulon_diagmc();
	
	void update(void);
	void fill_histograms(void);
	void reset_histograms(void);
	void measure(void);
	double fraction_completed(void) const;

	double histogram_bin_center(int index);	
	double histogram_index_time(double time);

	void print_data(const alps::results_type<angulon_diagmc>::type &results);
};

/*
	The constructor for our simulation.

	We always need the parameters and the seed as we need to pass it to the alps::mcbase constructor.
*/

angulon_diagmc::angulon_diagmc(parameters_type const &params, std::size_t seed_offset=42) : alps::mcbase(params, seed_offset)
{
	sweeps=0;
	total_sweeps=params["mc.totalsweeps"];
	thermalization_sweeps=params["mc.thermalization_sweeps"];

	nbins=params["mc.nbins"];
	nhist=params["mc.nhist"];
	nloop=params["mc.nloop"];

	config.j=params["parameters.j"];
	config.endtau=params["parameters.endtau"];
	config.chempot=params["parameters.chempot"];
	config.maxtau=params["parameters.maxtau"];
	config.maxorder=params["parameters.maxorder"];
	config.n=exp(params["potential.logn"]);

	binwidth=config.maxtau/nbins;

	measurements << alps::accumulators::FullBinningAccumulator<std::vector<double>>("G");
	measurements << alps::accumulators::FullBinningAccumulator<std::vector<double>>("G0");
	measurements << alps::accumulators::FullBinningAccumulator<std::vector<double>>("G1");
	measurements << alps::accumulators::FullBinningAccumulator<double>("intG0");
	measurements << alps::accumulators::FullBinningAccumulator<double>("order");
        measurements << alps::accumulators::MeanAccumulator<std::vector<double>>("order_frequencies");

	updates.push_back(std::make_pair(update_length,"UpdateLength"));
	updates.push_back(std::make_pair(update_add_phonon_line,"AddPhononLine"));
	updates.push_back(std::make_pair(update_remove_phonon_line,"RemovePhononLine"));
	updates.push_back(std::make_pair(update_shift_vertex,"ShiftVertex"));
	updates.push_back(std::make_pair(update_swap_deltajs,"SwapDeltajs"));
	updates.push_back(std::make_pair(update_change_mu,"ChangeMu"));

	int nr_updates=updates.size();
	
	accepted.resize(nr_updates);
	rejected.resize(nr_updates);
	
	physical_updates.resize(1+config.maxorder);
	unphysical_updates.resize(1+config.maxorder);

	order_frequencies.resize(1+config.maxorder);
	hist_g.resize(nbins);
	hist_g0.resize(nbins);
	hist_g1.resize(nbins);

	{
		struct diagram_parameters_t dpars;
		
	        dpars.j=config.j;
	        dpars.endtau=config.endtau;
	        dpars.maxtau=config.maxtau;
	        dpars.chempot=config.chempot;
	        dpars.n=config.n;

	        dgr=init_diagram(&dpars,true);
		
		if(dgr==NULL)
		{
			std::cout << "Couldn't initialize diagram!" << std::endl;
			throw std::exception();
		}

		dgr->rng_ctx=gsl_rng_alloc(gsl_rng_mt19937);

		if(dgr->rng_ctx==NULL)
		{
			std::cout << "Couldn't initialize internal (GSL) RNG!" << std::endl;
			throw std::exception();
		}
	}
}

angulon_diagmc::~angulon_diagmc()
{
	fini_diagram(dgr);
}


/*
	 This performs the actual calculation at each MC step.
*/

void angulon_diagmc::update(void)
{
	reset_histograms();

	for(int j=0;j<nloop;j++)
	{
        	for(int i=0;i<nhist;i++)
		{
			int update_type,status;
		
			update_type=int(3*random())%3;

			status=updates[update_type].first(dgr,&config);

			switch(status)
			{
				case UPDATE_ACCEPTED:
				accepted[update_type]++;
				break;

				case UPDATE_UNPHYSICAL:
				case UPDATE_REJECTED:
				rejected[update_type]++;
				break;

				case UPDATE_ERROR:
				break;
			}

			int diagram_order=get_nr_phonons(dgr);

			if(configuration_is_physical(dgr)==true)
			{
				physical_updates[diagram_order]++;
			}
			else
			{
				unphysical_updates[diagram_order]++;
			}
		}

		fill_histograms();
	}
}

void angulon_diagmc::fill_histograms()
{
	int itime=histogram_index_time(dgr->endtau);

	hist_g[itime] += dgr->sign;

	int diagram_order=get_nr_phonons(dgr);

	hist_order += diagram_order;
	order_frequencies[diagram_order] += 1.;
	
	if(diagram_order==0)
	{
		intG0 += 1.;
		hist_g0[itime] += dgr->sign;
	}
	
	if(diagram_order==1)
	{
		hist_g1[itime] += dgr->sign;
	}
}

void angulon_diagmc::reset_histograms()
{
	for(int i=0;i<nbins;i++)
	{
		hist_g[i] = 0.0f;
		hist_g0[i] = 0.0f;
		hist_g1[i] = 0.0f;
	}

	std::fill(order_frequencies.begin(),order_frequencies.end(),0.0f);

	intG0 = 0.0f;
	hist_order = 0.0f;
}

/*
	This collects the measurements at each MC step.
*/

void angulon_diagmc::measure(void)
{
	sweeps++;

	if (sweeps<thermalization_sweeps)
		return;

	std::vector<double> GF(nbins, 0.);
	std::vector<double> GF0(nbins, 0.);
	std::vector<double> GF1(nbins, 0.);

	for (size_t j=0; j < nbins; j++)
	{
		GF[j] = hist_g[j]/binwidth;
		GF0[j] = hist_g0[j]/binwidth;
		GF1[j] = hist_g1[j]/binwidth;
	}

	for (auto it = order_frequencies.begin(); it != order_frequencies.end(); it++)
        	*it /= nloop;

	measurements["G"] << GF;
	measurements["G0"] << GF0;
	measurements["G1"] << GF1;
	measurements["intG0"] << intG0;
	measurements["order"] << hist_order / nloop;
	measurements["order_frequencies"] << order_frequencies;
}

/*
	A number in the range from 0.0 to 1.0 that says how much
	of the simulation has been completed
*/

double angulon_diagmc::fraction_completed() const
{	
	std::cout << "pct " << ((sweeps > thermalization_sweeps) ? (sweeps - thermalization_sweeps) / double(total_sweeps) : 0.) * 100 << std::endl;
	
	return ((sweeps > thermalization_sweeps) ? (sweeps - thermalization_sweeps) / double(total_sweeps) : 0.);
}

double angulon_diagmc::histogram_bin_center(int index)
{
	if(index>=nbins)
		index=nbins-1;

	if(index<0)
		index=0;

	return binwidth*(0.5f+index);
}

double angulon_diagmc::histogram_index_time(double time)
{
	int index=int(time/binwidth);

	if(index<0)
		return 0;
	
	if(index>=nbins)
		return nbins-1;
	
	return index;
}

void angulon_diagmc::print_data(const alps::results_type<angulon_diagmc>::type &results)
{
	double Ej=config.j*(config.j+1.0f);
	double I0=(1.0f-exp(-(Ej-config.chempot)*config.maxtau))/(Ej-config.chempot);

        const alps::accumulators::result_wrapper& G0norm=I0*results["G0"]/results["intG0"];
        const alps::accumulators::result_wrapper& G1norm=I0*results["G1"]/results["intG0"];
        const alps::accumulators::result_wrapper& Gnorm=I0*results["G"]/results["intG0"];
        const alps::accumulators::result_wrapper& logGnorm=log(Gnorm);

	for(int c=0;c<nbins/10;c++)
	{
		double bincenter=histogram_bin_center(c);
		double free_rotor_energy=config.j*(config.j+1.0f);

		std::cout << histogram_bin_center(c) << " ";
		std::cout << Gnorm.mean<std::vector<double>>()[c] << " ";
		std::cout << Gnorm.error<std::vector<double>>()[c] << " ";
		std::cout << exp(-(free_rotor_energy-config.chempot)*bincenter) << std::endl;
	}
}

/*
	Simulation entry point
*/

int main(int argc, char **argv)
{
	/*
		Creates the parameters for the simulation
	*/

	std::cout << "Initializing parameters..." << std::endl;
	alps::parameters_type<angulon_diagmc>::type params(argc, argv);

	/*
		Define the parameters for our simulation
	*/

	params.define<unsigned long>("mc.totalsweeps", 100000, "Total number of sweeps");
	params.define<unsigned long>("mc.thermalization_sweeps", 50, "Number of thermalization sweeps");
	params.define<size_t>("mc.nloop", 1000, "Number of histogram samples gathered before measuring");
	params.define<size_t>("mc.nhist", 1,"Number of updates until counting configuration towards histograms.");
	params.define<size_t>("mc.nbins", 10000, "Number of bins");

	params.define<std::string>("general.outputfile", "test.output", "Output file name");

	params.define<int>("parameters.j", 0, "Angular momentum of the external line");
	params.define<double>("parameters.endtau", 1.0f, "Initial diagram length");
	params.define<double>("parameters.chempot", -4.90f, "Chemical potential");
	params.define<double>("parameters.maxtau", 100.0f, "Maximum diagram length");
	params.define<int>("parameters.maxorder", 100, "Maximum diagram order");

	params.define<double>("potential.logn", 5.0f, "Logarithm of the density (interaction parameter)");

	angulon_diagmc::define_parameters(params);

	if(params.help_requested(std::cout))
	{
		return 0;
	}
   
	/*
		Create and run the simulation
	*/
	
	std::cout << "Running simulation..." << std::endl;

	angulon_diagmc *mysim;

	try
	{	
		mysim=new angulon_diagmc(params);
		mysim->run(alps::stop_callback(0));
		//mysim->run([]{return false;});
	}
	catch(const std::exception&)
	{
		return EXIT_FAILURE;
	}

	/*
		Collect the results from the simulation
	*/

	std::cout << "Collecting and printing results..." << std::endl;

	alps::results_type<angulon_diagmc>::type results = alps::collect_results(*mysim);
	mysim->print_data(results);

	/*
		Saving to the output file
	*/
	
	std::cout << "Collecting results..." << std::endl;

	std::string output_file = params["general.outputfile"];
	alps::hdf5::archive ar(output_file, "w");
	ar["/parameters"] << params;
	ar["/simulation/results"] << results;
	
	delete mysim;
	
	return EXIT_SUCCESS;
}
