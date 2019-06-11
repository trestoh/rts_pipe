// oms_tide_2013_tests.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <ctime>
#include <map>
#include <vector>
#include <gflags/gflags.h>
#include <ratio>
#include <windows.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#define BOOST_ALL_NO_LIB
#define CRUX_VERSION "3.1"

#include <util/GlobalParams.h>
#include <util/WinCrux.h>

#include <app/tide/abspath.h>
//#include <app/tide/records_to_vector-inl.h>
#include <app/tide/records.h>

#include <io/carp.h>
#include <parameter.h>
#include <io/SpectrumRecordWriter.h>
#include <app/TideIndexApplication.h>
#include <app/ParamMedicApplication.h>
#include <app/PSMConvertApplication.h>
#include <app/tide/mass_constants.h>
#include <app/TideMatchSet.h>
#include <app/TideSearchApplication.h>
#include <util/Params.h>
#include <util/FileUtils.h>
#include <util/StringUtils.h>
#include <util/ArgParser.h>
#include <io/SpectrumCollectionFactory.h>

//#include "peptides.pb.h"
//#include "spectrum.pb.h"
//#include "app/tide/theoretical_peak_set.h"
//#include "app/tide/max_mz.h"

#include "inmemindex.h"
#include "mem_peptide_queue.h"
#include "rt_sequencer.h"

//const double XCORR_SCALING = 100000000.0;
const double CYSTEINE_DEFAULT = 57.021464;

//global vars for timing
double PCFreq = 0.0;
__int64 CounterStart = 0;

void StartCounter()
{
	LARGE_INTEGER li;
	if (!QueryPerformanceFrequency(&li))
		std::cout << "QueryPerformanceFrequency failed!\n";

	PCFreq = double(li.QuadPart) / 1000.0;

	QueryPerformanceCounter(&li);
	CounterStart = li.QuadPart;
}
double GetCounter()
{
	LARGE_INTEGER li;
	QueryPerformanceCounter(&li);
	return double(li.QuadPart - CounterStart) / PCFreq;
}


/* This constant is used to put the refactored XCorr back into the
* same range as the original XCorr score.  It is the XCorr "magic
* number" (10000) divided by the EVIDENCE_SCALE_INT (defined in
* tide/spectrum_preprocess2.cc). */
const double RESCALE_FACTOR = 20.0;

struct InputFile {
	std::string OriginalName;
	std::string SpectrumRecords;
	bool Keep;
	InputFile(const std::string& name,
		const std::string& spectrumrecords,
		bool keep) :
		OriginalName(name), SpectrumRecords(spectrumrecords), Keep(keep) {}
};


struct specStruct
{
	double key;
	int index;

	specStruct(const double& mz, int ind): key(mz), index(ind) {}

	bool operator < (const specStruct& str) const
	{
		return (key < str.key);
	}
};


int main(int argc, char * argv[])
{
#ifdef _MSC_VER
	// Turn off auto-tranlation of line-feed to 
	// carriage-return/line-feed
	_set_fmode(_O_BINARY);
#endif 

	set_verbosity_level(CARP_DEBUG);
	
	string output_file_name_;
	
	bool HAS_DECOYS = false;

	int spectrum_it = 0;

	//variables for parameters and such
	const string appName = "rts-pipe";
	const std::vector<std::string> appArgs;
	const std::vector<std::string> appOptions;

	//Implemented in parameter.cpp, I am not entirely sure what this function does
	initialize_parameters();

	//below parses what arguments we do have
	ArgParser argParser;
	try {

		argParser.Parse(argc, argv, appArgs);

		const map<string, string>& options = argParser.GetOptions();

		// Read parameter file if specified
		string parameter_file = argParser.GetOption("parameter-file");
		if (!parameter_file.empty()) {
			parse_parameter_file(parameter_file.c_str());
		}

		// Process command line options
		for (map<string, string>::const_iterator i = options.begin(); i != options.end(); i++) {
			if (i->first != "oms-src" && i->first != "oms-index" && i->first != "xc1" && i->first != "xc2" && i->first != "xc3" && i->first != "deltacn" )
			{
				try {
					Params::Set(i->first, i->second);
				}
				catch (const runtime_error& e) {
					std::cout << "Failed option add" << std::endl;
					throw ArgParserException(e.what());
				}
			}
		}

		// Process command line arguments
		const map< string, vector<string> >& args = argParser.GetArgs();
		for (map< string, vector<string> >::const_iterator i = args.begin(); i != args.end(); i++) {
			for (vector<string>::const_iterator j = i->second.begin(); j != i->second.end(); j++) {
				try {
					Params::AddArgValue(i->first, *j);
				}
				catch (const runtime_error& e) {
					std::cout << "Failed arg add" << std::endl;
					throw ArgParserException(e.what());
				}
			}
		}
	}

	catch (const runtime_error& e) {
		carp(CARP_FATAL, "%s", e.what());
	}

	set_verbosity_level(Params::GetInt("verbosity"));

	// Set seed for random number generation 
	if (Params::GetString("seed") == "time") {
		time_t seconds; // use current time to seed
		time(&seconds); // Get value from sys clock and set seconds variable.
		mysrandom((unsigned)seconds); // Convert seconds to a unsigned int
	}
	else {
		mysrandom(StringUtils::FromString<unsigned>(Params::GetString("seed")));
	}

	//re parse arguments looking for the custom options we've added
	argParser.Parse(argc, argv, appArgs);
	const map<string, string> options = argParser.GetOptions();
	map<string, string>::const_iterator oms_arg = options.find("oms-src");
	map<string, string>::const_iterator oms_index = options.find("oms-index");
	std::string oms = oms_arg->second;
	std::string index = oms_index->second;

	double xcorr1, xcorr2, xcorr3, deltacn;

	map<string, string>::const_iterator oms_xc1 = options.find("xc1");
	map<string, string>::const_iterator oms_xc2 = options.find("xc2");
	map<string, string>::const_iterator oms_xc3 = options.find("xc3");
	map<string, string>::const_iterator oms_deltacn = options.find("deltacn");

	if (oms_xc1 != options.end())
		xcorr1 = atof(oms_xc1->second.c_str());
	else
		xcorr1 = 1.8;

	if (oms_xc2 != options.end())
		xcorr2 = atof(oms_xc2->second.c_str());
	else
		xcorr2 = 2.2;

	if (oms_xc3 != options.end())
		xcorr3 = atof(oms_xc3->second.c_str());
	else
		xcorr3 = 3.5;

	if (oms_deltacn != options.end())
		deltacn = atof(oms_deltacn->second.c_str());
	else
		deltacn = 0.08;

	// Create output directory if appliation needs it.
	// Create output directory 
	string output_folder = Params::GetString("output-dir");
	if (create_output_directory(output_folder, Params::GetBool("overwrite")) == -1) {
		carp(CARP_FATAL, "Unable to create output directory %s.", output_folder.c_str());
	}

	// Open the log file to record carp messages 
	open_log_file(appName + ".log.txt");

	// Store the host name, start date and time, version number, and command line.
	carp(CARP_INFO, "CPU: %s", hostname());
	carp(CARP_INFO, "Crux version: %s", CRUX_VERSION);
	carp(CARP_INFO, date_and_time());

	// Write the parameter file
	string paramFile = make_file_path(appName + ".params.txt");
	ofstream* file = FileUtils::GetWriteStream(paramFile, Params::GetBool("overwrite"));
	if (file == NULL) {
		throw runtime_error("Could not open " + paramFile + " for writing");
	}

	Params::Write(file);
	delete file;

	/*
	 *	This is the point where we begin to emulate the behavior of the tide search app main and manipulate it to use non-file data
	 */

	carp(CARP_INFO, "Running tide-search...");

	/*
	
		This section allows for threading, which we are not currently using but probably could? (set up index in one space
		of memory and have threads as querying entities. Have to refactor things
	
	*/

	// prevent different output formats from using threading
	//if (Params::GetBool("peptide-centric-search") == true) {
	//	NUM_THREADS = 1;
	//	carp(CARP_INFO, "Threading for peptide-centric formats are not supported yet");
	//}
	//else {
		//NUM_THREADS = Params::GetInt("num-threads");
	//	NUM_THREADS = 1;
	//}
	//if (NUM_THREADS < 1) {
	//	NUM_THREADS = boost::thread::hardware_concurrency(); // MINIMUM # = 1.
		// (Meaning just main thread) Do not make this value below 1.
		// make sure that number of threads are reasonable, e.g. user did not specify millions of threads...
	//}
	//else if (NUM_THREADS > 64) {
	//	carp(CARP_FATAL, "Requested more than 64 threads.");
	//}
	//carp(CARP_INFO, "Number of Threads: %d", NUM_THREADS); // prints the number of threads

	/*
			END THREADING SETUP
	*/

	// Check concat parameter
	bool concat = Params::GetBool("concat");

	// Check compute-sp parameter
	bool compute_sp = Params::GetBool("compute-sp");
	if (Params::GetBool("sqt-output") && !compute_sp) {
		compute_sp = true;
		carp(CARP_INFO, "Enabling parameter compute-sp since SQT output is enabled "
			"(this will increase runtime).");
	}
	
	//Everything below here is an aspect of "main" and will not be used in the library version of the program

	ofstream* target_file = NULL;
	ofstream* decoy_file = NULL;

	carp(CARP_DEBUG, "Using TideMatchSet to write matches");
	bool overwrite = Params::GetBool("overwrite");
	if (!concat) {
		string target_file_name = make_file_path("tide-search-junk.target.txt");
		target_file = create_stream_in_path(target_file_name.c_str(), NULL, overwrite);
		output_file_name_ = target_file_name;
		if (HAS_DECOYS) {
			string decoy_file_name = make_file_path("tide-search-junk.decoy.txt");
			decoy_file = create_stream_in_path(decoy_file_name.c_str(), NULL, overwrite);
		}
	}
	else {
		string concat_file_name = make_file_path("tide-search-junk.txt");
		target_file = create_stream_in_path(concat_file_name.c_str(), NULL, overwrite);
		output_file_name_ = concat_file_name;
	}

	if (target_file) {
		carp(CARP_INFO, "Target file successfully created");
		TideMatchSet::writeHeaders(target_file, false, compute_sp);
		TideMatchSet::writeHeaders(decoy_file, true, compute_sp);
	}

	//inittialize sequencer, may want to have initalize as a separate function

	RT_Sequencer* sequencer = new RT_Sequencer();
	//sequencer->init(index.c_str(), xcorr1, xcorr2, xcorr3, deltacn);
	//sequencer->init(index.c_str(), NULL, NULL, 1.8, 2.2, 3.5, 0.08, 3, 20.0);
	sequencer->init(index.c_str(), NULL, NULL, 1.8, 1.3, 1.5, 0.08, 3, 20.0);


	carp(CARP_INFO, "Re-parsed succefully");

	OpenMS::MzMLFile mzMLDataFileProfile;
	OpenMS::MSExperiment msExperimentProfile;
	std::string write_file = "";
	std::string junk;
	std::string scan;

	try {
		mzMLDataFileProfile.load(oms.c_str(), msExperimentProfile);
	}
	catch (std::exception& e) {
		std::cout << e.what() << std::endl;
		//usage();
		return 1;
	}

	OpenMS::MSSpectrum<> s = msExperimentProfile.getSpectrum(spectrum_it);
	std::istringstream natID = std::istringstream(s.getNativeID());

	carp(CARP_INFO, "Open MS experiment loaded");

	double spectrum_min_mz = Params::GetDouble("spectrum-min-mz");
	double spectrum_max_mz = Params::GetDouble("spectrum-max-mz");
	int min_peaks = Params::GetInt("min-peaks");
	int top_matches = Params::GetInt("top-match");
	double highest_mz = 0.0f;

	// Start the timer.
	wall_clock();
	
	//various counting variables
	int db_hits = 0;
	int ms2_count = 0;
	int ms1_count = 0;
	int ms3_count = 0;
	int decoy_count = 0;
	int target_count = 0;


	ofstream* target_hits = NULL;
	string target_hit_file_name = make_file_path("target_hit_seq.txt");
	target_hits = create_stream_in_path(target_hit_file_name.c_str(), NULL, overwrite);

	string fd_summary_file_name = make_file_path("fd_rates.txt");
	ofstream fd_rates(fd_summary_file_name, std::ios_base::app);

	string times_out_name = make_file_path("search_times.txt");
	ofstream times_out(times_out_name, std::ios_base::app);

	//outputtting the parameters
	(*target_hits) << xcorr1 << " " << xcorr2 << " " << xcorr3 << " " << deltacn << std::endl;

	double max_time = 0.0;
	std::vector<double> times;
	std::vector<std::string> pep_hits;

	carp(CARP_INFO, "Hitting main loop");

	for (int i = spectrum_it; i < msExperimentProfile.getNrSpectra(); i++)
	{
		//non-sorted version
		s = msExperimentProfile.getSpectrum(i);
	
		//sorted version
		//s = msExperimentProfile.getSpectrum(specs[i].index);
		if (s.getMSLevel() == 1)
			ms1_count++;

		if (s.getMSLevel() == 3)
			ms3_count++;

		if (s.getMSLevel() == 2)
		{
			StartCounter();

			++ms2_count;

			natID = std::istringstream(s.getNativeID());
			natID >> junk >> junk >> scan;
			scan = scan.substr(scan.find('=') + 1);

			if (scan == "10663")
			{
				std::cout << "Scan of interest here" << std::endl;
			}

			Spectrum sspec = Spectrum(atoi(scan.c_str()), s.getPrecursors()[0].getMZ());

			sequencer->makeSpec(atoi(scan.c_str()), s.getPrecursors()[0].getMZ());


			for (OpenMS::MSSpectrum<>::ConstIterator it = s.begin(); it != s.end(); ++it)
			{
				sspec.AddPeak(it->getPos(), it->getIntensity());

				sequencer->addPeak(it->getPos(), it->getIntensity());

				//if (atoi(scan.c_str()) == 16613)
				//{
				//	carp(CARP_INFO, "testSpec.AddPeak(%f,%f);", it->getPos(), it->getIntensity());
				//}
			}			

			highest_mz = (s.end() - 1)->getPos();

			std::string pep_seq;
			bool isDecoy;

			

			//bool db_match = sequencer->is_match(sspec, highest_mz, s.getPrecursors()[0].getUnchargedMass(), s.getPrecursors()[0].getCharge(), pep_seq, isDecoy);
			bool db_match = sequencer->is_match(highest_mz, s.getPrecursors()[0].getUnchargedMass(), s.getPrecursors()[0].getCharge());


			if (db_match)
			{

				//carp(CARP_INFO, "Hit: Scan No: %d, prec mz: %f, prec uncharged mass: %f prec charge: %d mighest mz: %f", atoi(scan.c_str()), s.getPrecursors()[0].getMZ(), s.getPrecursors()[0].getUnchargedMass(), s.getPrecursors()[0].getCharge(), highest_mz);

				db_hits++;
				/*
				if (isDecoy)
					++decoy_count;
				else {
					++target_count;
					pep_hits.push_back(pep_seq);
				}
				*/
			}

			//more file non-sense
			PSMConvertApplication converter;

			// convert tab delimited to other file formats.
			if (!concat) {
				string target_file_name = make_file_path("tide-search.target.txt");
				if (Params::GetBool("pin-output")) {
					converter.convertFile("tsv", "pin", target_file_name, "tide-search.target.", Params::GetString("protein-database"), true);
				}
				if (Params::GetBool("pepxml-output")) {
					converter.convertFile("tsv", "pepxml", target_file_name, "tide-search.target.", Params::GetString("protein-database"), true);
				}
				if (Params::GetBool("mzid-output")) {
					converter.convertFile("tsv", "mzidentml", target_file_name, "tide-search.target.", Params::GetString("protein-database"), true);
				}
				if (Params::GetBool("sqt-output")) {
					converter.convertFile("tsv", "sqt", target_file_name, "tide-search.target.", Params::GetString("protein-database"), true);
				}

				if (HAS_DECOYS) {
					string decoy_file_name = make_file_path("tide-search.decoy.txt");
					if (Params::GetBool("pin-output")) {
						converter.convertFile("tsv", "pin", decoy_file_name, "tide-search.decoy.", Params::GetString("protein-database"), true);
					}
					if (Params::GetBool("pepxml-output")) {
						converter.convertFile("tsv", "pepxml", decoy_file_name, "tide-search.decoy.", Params::GetString("protein-database"), true);
					}
					if (Params::GetBool("mzid-output")) {
						converter.convertFile("tsv", "mzidentml", decoy_file_name, "tide-search.decoy.", Params::GetString("protein-database"), true);
					}
					if (Params::GetBool("sqt-output")) {
						converter.convertFile("tsv", "sqt", decoy_file_name, "tide-search.decoy.", Params::GetString("protein-database"), true);
					}
				}
			}
			else {
				string concat_file_name = make_file_path("tide-search.txt");
				if (Params::GetBool("pin-output")) {
					converter.convertFile("tsv", "pin", concat_file_name, "tide-search.", Params::GetString("protein-database"), true);
				}
				if (Params::GetBool("pepxml-output")) {
					converter.convertFile("tsv", "pepxml", concat_file_name, "tide-search.", Params::GetString("protein-database"), true);
				}
				if (Params::GetBool("mzid-output")) {
					converter.convertFile("tsv", "mzidentml", concat_file_name, "tide-search.", Params::GetString("protein-database"), true);
				}
				if (Params::GetBool("sqt-output")) {
					converter.convertFile("tsv", "sqt", concat_file_name, "tide-search.", Params::GetString("protein-database"), true);
				}
			}

			double elapsed = GetCounter();
			times.push_back(elapsed);

		} //end search for spectra that are MS2

	} //end individual spectra loop

	double false_discovery = ((double)decoy_count) / db_hits;

	//various information output
	carp(CARP_INFO, "Elapsed total time: %.3g s", wall_clock() / 1e6);
	carp(CARP_INFO, "Elapsed time per MS2 spectrum conversion: %.3g s", wall_clock() / (1e6*ms2_count) );
	carp(CARP_INFO, "There are %d MS1 spectra in the input file.", ms1_count);
	carp(CARP_INFO, "There are %d MS2 spectra in the input file.", ms2_count);
	carp(CARP_INFO, "There are %d MS3 spectra in the input file.", ms3_count);
	carp(CARP_INFO, "Of those spectra, %d were hits in the search for the current criteria", db_hits);
	carp(CARP_INFO, "Of the hits %d were decoys and %d were targets", decoy_count, target_count);
	carp(CARP_INFO, "False Discovery rate: %f", false_discovery);

	fd_rates << "For Parameters: (" << xcorr1 << ", " << xcorr2 << ", " << xcorr3 << ", " << deltacn << ") FD rate is " << false_discovery << std::endl;

	for (int i = 0; i < pep_hits.size(); i++)
	{
		(*target_hits) << pep_hits[i] << std::endl;
	}

	for (int i = 0; i < times.size(); i++)
	{
		times_out << times[i] << std::endl;
	}

	//end various outputs

	if (target_file) {
		delete target_file;
		if (decoy_file) {
			delete decoy_file;
		}
	}
	
	delete sequencer;

	return 0;
}
