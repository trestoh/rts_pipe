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
#include <app/tide/records_to_vector-inl.h>
#include <app/tide/records.h>

//#include <app/TideSearchApplication.h>

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

#include "peptides.pb.h"
#include "spectrum.pb.h"
#include "app/tide/theoretical_peak_set.h"
#include "app/tide/max_mz.h"

#include "inmemindex.h"
#include "mem_peptide_queue.h"

//TideSearchApplication::HAS_DECOYS = true;

const double XCORR_SCALING = 100000000.0;

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

void collectScoresCompiled(
	//ActivePeptideQueue* active_peptide_queue,
	ActivePeptideQueue2* active_peptide_queue,
	const Spectrum* spectrum,
	const ObservedPeakSet& observed,
	TideMatchSet::Arr2* match_arr,
	int queue_size,
	int charge);

void computeWindow(
	const SpectrumCollection::SpecCharge& sc,
	WINDOW_TYPE_T window_type,
	double precursor_window,
	int max_charge,
	vector<int>* negative_isotope_errors,
	vector<double>* out_min,
	vector<double>* out_max,
	double* min_range,
	double* max_range
	);

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

	/*
	OpenMS::MzMLFile mzMLDataFileProfile;
	OpenMS::MSExperiment msExperimentProfile;
	try {
	mzMLDataFileProfile.load(argv[1], msExperimentProfile);
	}
	catch (std::exception& e) {
	std::cout << e.what() << std::endl;
	//usage();
	return 1;
	}
	*/

	set_verbosity_level(CARP_DEBUG);

	//std::ofstream out(argv[2]);
	//std::ofstream calc_out(argv[3]);

	//ofstream for when we are writing to spectrum for tide to use
	std::ofstream out("./temp/spec_0.ms2");
	out.close();

	std::string remove_index_ = "";
	string output_file_name_;
	const double CYSTEINE_DEFAULT = 57.021464;
	bool HAS_DECOYS = false;
	bool PROTEIN_LEVEL_DECOYS = false;

	int spectrum_it = 0;

	//test our ability to parse arguments
	unsigned int NUM_THREADS;

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
		//options = argParser.GetOptions();
		for (map<string, string>::const_iterator i = options.begin(); i != options.end(); i++) {
			if (i->first != "oms-src" && i->first != "oms-index")
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

	//std::cout << "Trying to make a Search App" << std::endl;
	//TideSearchApplication pass;
	//std::cout << "Made the search app" << std::endl;

	//PROBLEM AREA, TEAR OUT THE OLD FUNCTION, GET RID OF THE APP??
	/*
	*
	*
	* PROCESS PARAMS FOR TIDE SEARCH BY HAND
	* Can't get this to work so it has been eliminated for now. I think the main replicates some of this functionality anyway?
	*
	*/

	//
	//	EXTREMELY DUMB CODE: Index hard coded, needs to be changed
	//

	argParser.Parse(argc, argv, appArgs);
	const map<string, string> options = argParser.GetOptions();
	map<string, string>::const_iterator oms_arg = options.find("oms-src");
	map<string, string>::const_iterator oms_index = options.find("oms-index");
	std::string oms = oms_arg->second;
	std::string index = oms_index->second;

	//string index = "C:\\dev\\real_time_seq\\oms_tide_2013_tests\\oms_tide_2013_tests\\human-index";
	string peptides_file = FileUtils::Join(index, "pepix");
	string proteins_file = FileUtils::Join(index, "protix");
	string auxlocs_file = FileUtils::Join(index, "auxlocs");


	/*
	const string index = "C:\\dev\\real_time_seq\\oms_tide_2013_tests\\oms_tide_2013_tests\\human-index";
	if (!FileUtils::Exists(index)) {
	carp(CARP_FATAL, "'%s' does not exist", index.c_str());
	}
	*/
	/*
	*
	* Ignore case where database is not already an index for now, this involves invoking ANOTHER
	* application, this may be fine down the road but for now it is not worth worrying about
	*
	*
	else if (FileUtils::IsRegularFile(index)) {
	// Index is FASTA file
	carp(CARP_INFO, "Creating index from '%s'", index.c_str());
	string targetIndexName = Params::GetString("store-index");
	if (targetIndexName.empty()) {
	targetIndexName = FileUtils::Join(Params::GetString("output-dir"),
	"tide-search.tempindex");
	remove_index_ = targetIndexName;
	}

	string default_cysteine = "C+" + StringUtils::ToString(CYSTEINE_DEFAULT);
	string mods_spec = Params::GetString("mods-spec");
	if (mods_spec.find('C') == string::npos) {
	mods_spec = mods_spec.empty() ?
	default_cysteine : default_cysteine + ',' + mods_spec;
	carp(CARP_DEBUG, "Using default cysteine mod '%s' ('%s')",
	default_cysteine.c_str(), mods_spec.c_str());
	}
	Params::Set("mods-spec", mods_spec);

	if (indexApp.main(index, targetIndexName) != 0) {
	carp(CARP_FATAL, "tide-index failed.");
	}
	Params::Set("tide database", targetIndexName);
	}
	*/

	//problems here, peptide stuff is failing


	// Index is Tide index directory
	pb::Header pepts_header;

	HeadedRecordReader pept_reader(peptides_file, &pepts_header);
	if ((pepts_header.file_type() != pb::Header::PEPTIDES) ||
		!pepts_header.has_peptides_header()) {
		carp(CARP_FATAL, "Error reading index (%s)", peptides_file.c_str());
	}

	const pb::Header::PeptidesHeader& pepsHeader = pepts_header.peptides_header();

	Params::Set("enzyme", pepsHeader.enzyme());
	const char* digestString = digest_type_to_string(pepsHeader.full_digestion() ? FULL_DIGEST : PARTIAL_DIGEST);
	Params::Set("digestion", digestString);
	Params::Set("isotopic-mass", pepsHeader.monoisotopic_precursor() ? "mono" : "average");

	// run param-medic? 
	// May need param-medic later but for now it seems like no
	/*
	const string autoPrecursor = Params::GetString("auto-precursor-window");
	const string autoFragment = Params::GetString("auto-mz-bin-width");
	if (autoPrecursor != "false" || autoFragment != "false") {
	if (autoPrecursor != "false" && Params::GetString("precursor-window-type") != "ppm") {
	carp(CARP_FATAL, "Automatic peptide mass tolerance detection is only supported with ppm "
	"units. Please rerun with either auto-precursor-window set to 'false' or "
	"precursor-window-type set to 'ppm'.");
	}
	ParamMedicErrorCalculator errCalc;
	errCalc.processFiles(Params::GetStrings("tide spectra file"));
	string precursorFailure, fragmentFailure;
	double precursorSigmaPpm = 0;
	double fragmentSigmaPpm = 0;
	double fragmentSigmaTh = 0;
	double precursorPredictionPpm = 0;
	double fragmentPredictionPpm = 0;
	double fragmentPredictionTh = 0;
	errCalc.calcMassErrorDist(&precursorFailure, &fragmentFailure,
	&precursorSigmaPpm, &fragmentSigmaPpm,
	&precursorPredictionPpm, &fragmentPredictionTh);

	if (autoPrecursor != "false") {
	if (precursorFailure.empty()) {
	carp(CARP_INFO, "precursor ppm standard deviation: %f", precursorSigmaPpm);
	carp(CARP_INFO, "Precursor error estimate (ppm): %.2f", precursorPredictionPpm);
	Params::Set("precursor-window", precursorPredictionPpm);
	}
	else {
	carp(autoPrecursor == "fail" ? CARP_FATAL : CARP_ERROR,
	"failed to calculate precursor error: %s", precursorFailure.c_str());
	}
	}
	if (autoFragment != "false") {
	if (fragmentFailure.empty()) {
	carp(CARP_INFO, "fragment standard deviation (ppm): %f", fragmentSigmaPpm);
	carp(CARP_INFO, "Fragment bin size estimate (Th): %.4f", fragmentPredictionTh);
	Params::Set("mz-bin-width", fragmentPredictionTh);
	}
	else {
	carp(autoFragment == "fail" ? CARP_FATAL : CARP_ERROR,
	"failed to calculate fragment error: %s", fragmentFailure.c_str());
	}
	}
	}
	*/

	//pass.processParams();
	Params::Finalize();
	GlobalParams::set();

	//open_log_file( "tide-search.log.txt");
	//carp(CARP_DEBUG, "CPU: %s", hostname());

	// Store the host name, start date and time, version number, and command line.
	//carp(CARP_DEBUG, "CPU: %s", hostname());
	//carp(CARP_DEBUG, "Crux version: %s", CRUX_VERSION);
	//carp(CARP_DEBUG, date_and_time());
	//log_command_line(argc, argv);

	//carp(CARP_DEBUG, "Running tide-search...");

	carp(CARP_INFO, "Beginning %s.", appName.c_str());

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
	//log_command_line(argc, argv);
	//std::cout << "Passed log_command_line..." << std::endl;

	// Write the parameter file
	string paramFile = make_file_path(appName + ".params.txt");
	ofstream* file = FileUtils::GetWriteStream(paramFile, Params::GetBool("overwrite"));
	if (file == NULL) {
		throw runtime_error("Could not open " + paramFile + " for writing");
	}

	Params::Write(file);
	delete file;

	/*
	*
	*
	*
	*
	*	This is the point where we begin to emulate the behavior of the tide search app main and manipulate it to use non-file data
	*
	*
	*
	*
	*
	*/

	carp(CARP_INFO, "Running tide-search...");

	// prevent different output formats from using threading
	if (Params::GetBool("peptide-centric-search") == true) {
		NUM_THREADS = 1;
		carp(CARP_INFO, "Threading for peptide-centric formats are not supported yet");
	}
	else {
		//NUM_THREADS = Params::GetInt("num-threads");
		NUM_THREADS = 1;
	}
	if (NUM_THREADS < 1) {
		NUM_THREADS = boost::thread::hardware_concurrency(); // MINIMUM # = 1.
		// (Meaning just main thread) Do not make this value below 1.
		// make sure that number of threads are reasonable, e.g. user did not specify millions of threads...
	}
	else if (NUM_THREADS > 64) {
		carp(CARP_FATAL, "Requested more than 64 threads.");
	}
	carp(CARP_INFO, "Number of Threads: %d", NUM_THREADS); // prints the number of threads


	double window = Params::GetDouble("precursor-window");
	WINDOW_TYPE_T window_type = string_to_window_type(Params::GetString("precursor-window-type"));

	// Check spectrum-charge parameter
	string charge_string = Params::GetString("spectrum-charge");
	int charge_to_search;
	if (charge_string == "all") {
		carp(CARP_DEBUG, "Searching all charge states");
		charge_to_search = 0;
	}
	else {
		charge_to_search = atoi(charge_string.c_str());
		if (charge_to_search < 1 || charge_to_search > 6) {
			carp(CARP_FATAL, "Invalid spectrum-charge value %s", charge_string.c_str());
		}
		carp(CARP_DEBUG, "Searching charge state %d", charge_to_search);
	}

	// Check scan-number parameter
	string scan_range = Params::GetString("scan-number");
	int min_scan, max_scan;
	if (scan_range.empty()) {
		min_scan = 0;
		max_scan = BILLION;
		carp(CARP_DEBUG, "Searching all scans");
	}
	else if (scan_range.find('-') == string::npos) {
		// Single scan
		min_scan = max_scan = atoi(scan_range.c_str());
		carp(CARP_DEBUG, "Searching single scan %d", min_scan);
	}
	else {
		if (!get_range_from_string(scan_range.c_str(), min_scan, max_scan)) {
			carp(CARP_FATAL, "The scan number range '%s' is invalid. "
				"Must be of the form <first>-<last>", scan_range.c_str());
		}
		else {
			if (min_scan > max_scan) {
				int tmp_scan = min_scan;
				min_scan = max_scan;
				max_scan = tmp_scan;
				carp(CARP_DEBUG, "Switched scan range min and max");
			}
			carp(CARP_DEBUG, "Searching scan range %d-%d", min_scan, max_scan);
		}
	}

	//check to compute exact p-value 
	bool exact_pval_search = Params::GetBool("exact-p-value");
	double bin_width = Params::GetDouble("mz-bin-width");
	double bin_offset = Params::GetDouble("mz-bin-offset");

	// for now don't allow XCorr p-value searches with variable bin width
	if (exact_pval_search && !Params::IsDefault("mz-bin-width")) {
		carp(CARP_FATAL, "Tide-search with XCorr p-values and variable bin width "
			"is not allowed in this version of Crux.");
	}

	// Check concat parameter
	bool concat = Params::GetBool("concat");

	// Check compute-sp parameter
	bool compute_sp = Params::GetBool("compute-sp");
	if (Params::GetBool("sqt-output") && !compute_sp) {
		compute_sp = true;
		carp(CARP_INFO, "Enabling parameter compute-sp since SQT output is enabled "
			"(this will increase runtime).");
	}

	// Check isotope-error parameter
	string isotope_errors_string = Params::GetString("isotope-error");
	if (isotope_errors_string[0] == ',') {
		carp(CARP_FATAL, "Error in isotope_error parameter formatting: (%s)", isotope_errors_string.c_str());
	}
	for (int i = 0; isotope_errors_string[i] != '\0'; ++i) {
		if (isotope_errors_string[i] == ',' && (isotope_errors_string[i + 1] == ',' || isotope_errors_string[i + 1] == '\0')) {
			carp(CARP_FATAL, "Error in isotope_error parameter formatting: (%s) ", isotope_errors_string.c_str());
		}
	}

	vector<int>* negative_isotope_errors = new vector<int>();
	negative_isotope_errors->push_back(0);

	if (isotope_errors_string != "") {
		vector<int> isotope_errors = StringUtils::Split<int>(Params::GetString("isotope-error"), ',');
		for (vector<int>::iterator it = isotope_errors.begin(); it != isotope_errors.end(); ++it) {
			if (*it < 0) {
				carp(CARP_FATAL, "Found a negative isotope error: %d. There should not be any legitimate reasons to use negative isotope errors. Try modeling with a modification instead.", *it);
			}
			else if (find(negative_isotope_errors->begin(), negative_isotope_errors->end(), -1 * *it) != negative_isotope_errors->end()) {
				carp(CARP_FATAL, "Found duplicate when parsing isotope_error parameter: %d", *it);
			}
			negative_isotope_errors->push_back(-1 * *it);
		}
	}

	sort(negative_isotope_errors->begin(), negative_isotope_errors->end());

	//everything above here honestly seems good. So I'm getting worried now lol
	//this is where the "Process Params" started to get into trouble

	carp(CARP_INFO, "Reading index %s", index.c_str());

	//std::cout << "Reading Index: " << index.c_str() << std::endl;

	

	// Read proteins index file
	ProteinVec proteins;
	pb::Header protein_header;
	if (!ReadRecordsToVector<pb::Protein, const pb::Protein>(&proteins,
		proteins_file, &protein_header)) {
		carp(CARP_FATAL, "Error reading index (%s)", proteins_file.c_str());
	}

	double* aaFreqN = NULL;
	double* aaFreqI = NULL;
	double* aaFreqC = NULL;
	int* aaMass = NULL;
	int nAA = 0;

	// Read peptides index file
	pb::Header peptides_header;

	std::vector<HeadedRecordReader*> peptide_reader;
	for (int i = 0; i < NUM_THREADS; i++) {
		peptide_reader.push_back(new HeadedRecordReader(peptides_file, &peptides_header));
	}

	if ((peptides_header.file_type() != pb::Header::PEPTIDES) ||
		!peptides_header.has_peptides_header()) {
		carp(CARP_FATAL, "Error reading index (%s)", peptides_file.c_str());
	}

	const pb::Header::PeptidesHeader& pepHeader = peptides_header.peptides_header();
	DECOY_TYPE_T headerDecoyType = (DECOY_TYPE_T)pepHeader.decoys();
	if (headerDecoyType != NO_DECOYS) {
		HAS_DECOYS = true;
		if (headerDecoyType == PROTEIN_REVERSE_DECOYS) {
			PROTEIN_LEVEL_DECOYS = true;
		}
	}



	
	
	

	if (exact_pval_search) {
		pb::Header aaf_peptides_header;
		HeadedRecordReader aaf_peptide_reader(peptides_file, &aaf_peptides_header);

		if ((aaf_peptides_header.file_type() != pb::Header::PEPTIDES) ||
			!aaf_peptides_header.has_peptides_header()) {
			carp(CARP_FATAL, "Error reading index (%s)", peptides_file.c_str());
		}
		MassConstants::Init(&aaf_peptides_header.peptides_header().mods(),
			&aaf_peptides_header.peptides_header().nterm_mods(),
			&aaf_peptides_header.peptides_header().cterm_mods(),
			bin_width, bin_offset);
		ActivePeptideQueue* active_peptide_queue = new ActivePeptideQueue(aaf_peptide_reader.Reader(), proteins);
		//ActivePeptideQueue2* active_peptide_queue = new ActivePeptideQueue2(test_index, aaf_peptide_reader.Reader(), proteins);
		nAA = active_peptide_queue->CountAAFrequency(bin_width, bin_offset,
			&aaFreqN, &aaFreqI, &aaFreqC, &aaMass);
		delete active_peptide_queue;
	} // End calculation AA frequencies

	// Read auxlocs index file
	vector<const pb::AuxLocation*> locations;
	if (!ReadRecordsToVector<pb::AuxLocation>(&locations, auxlocs_file)) {
		carp(CARP_FATAL, "Error reading index (%s)", auxlocs_file.c_str());
	}
	carp(CARP_INFO, "Read %d auxlocs", locations.size());

	MassConstants::Init(&pepHeader.mods(), &pepHeader.nterm_mods(),
		&pepHeader.cterm_mods(), bin_width, bin_offset);
	ModificationDefinition::ClearAll();
	TideMatchSet::initModMap(pepHeader.mods(), ANY);
	TideMatchSet::initModMap(pepHeader.nterm_mods(), PEPTIDE_N);
	TideMatchSet::initModMap(pepHeader.cterm_mods(), PEPTIDE_C);

	MaxBin::SetGlobalMax(2000.0);

	//reading "test_index"
	InMemIndex* test_index = new InMemIndex(peptide_reader[0]->Reader(), proteins);

	/*
	for (int i = 0; i < test_index->pep_array.size(); i++)
	{
		Peptide* temp_temp = (test_index->pep_array)[i];
		carp(CARP_INFO, "Index has peptide %s with mass %f and programs at %d and %d", temp_temp->Seq(), temp_temp->Mass(), temp_temp->Prog(1), temp_temp->Prog(3));
	}
	*/
	//
	//	All below here I believe is involved with writing stuff to disk, meaning we will not want to use it and
	//	will instead want to reimplement with everything kept in memory
	//

	ofstream* target_file = NULL;
	ofstream* decoy_file = NULL;

	carp(CARP_DEBUG, "Using TideMatchSet to write matches");
	bool overwrite = Params::GetBool("overwrite");
	stringstream ss;
	ss << Params::GetString("enzyme") << '-' << Params::GetString("digestion");
	TideMatchSet::CleavageType = ss.str();
	if (!concat) {
		string target_file_name = make_file_path("tide-search.target.txt");
		target_file = create_stream_in_path(target_file_name.c_str(), NULL, overwrite);
		output_file_name_ = target_file_name;
		if (HAS_DECOYS) {
			string decoy_file_name = make_file_path("tide-search.decoy.txt");
			decoy_file = create_stream_in_path(decoy_file_name.c_str(), NULL, overwrite);
		}
	}
	else {
		string concat_file_name = make_file_path("tide-search.txt");
		target_file = create_stream_in_path(concat_file_name.c_str(), NULL, overwrite);
		output_file_name_ = concat_file_name;
	}

	if (target_file) {
		carp(CARP_INFO, "Target file successfully created");
		TideMatchSet::writeHeaders(target_file, false, compute_sp);
		TideMatchSet::writeHeaders(decoy_file, true, compute_sp);
	}

	// Try to read all spectrum files as spectrumrecords, convert those that fail
	// AKA: we must by hand convert spectrum into the format wanted. In the ultimate
	// version of the program this will be basically sitting around waiting for a trigger
	// but for now we can build a loop using our OpenMS stuff. Ultimately we need to figure
	// out how the spectrumrecord format works and how to build it up using OpenMS or maybe
	// something external.
	//
	//There's a bunch of file convert/clean shit that we also don't need
	/*
	vector<string> input_files;

	input_files.push_back("C:\\dev\\real_time_seq\\oms_tide_2013_tests\\oms_tide_2013_tests\\HELA_2017-04-27_294.mzML");

	vector<InputFile> input_sr;

	//single file version
	vector<string>::const_iterator f = input_files.begin();
	SpectrumCollection spectra;
	pb::Header spectrum_header;
	string spectrumrecords = *f;
	bool keepSpectrumrecords = true;
	if (!spectra.ReadSpectrumRecords(spectrumrecords, &spectrum_header)) {
	// Failed, try converting to spectrumrecords file
	carp(CARP_INFO, "Converting %s to spectrumrecords format", f->c_str());
	carp(CARP_INFO, "Elapsed time starting conversion: %.3g s", wall_clock() / 1e6);
	spectrumrecords = Params::GetString("store-spectra");
	keepSpectrumrecords = !spectrumrecords.empty();
	if (!keepSpectrumrecords) {
	spectrumrecords = make_file_path(FileUtils::BaseName(*f) + ".spectrumrecords.tmp");
	}
	else if (input_files.size() > 1) {
	carp(CARP_FATAL, "Cannot use store-spectra option with multiple input "
	"spectrum files");
	}
	carp(CARP_DEBUG, "New spectrumrecords filename: %s", spectrumrecords.c_str());

	//this is where convert/read needs to go
	//auto_ptr<Crux::SpectrumCollection> spectra(SpectrumCollectionFactory::create(f->c_str()));

	carp(CARP_DEBUG, "Reading converted spectra file %s", spectrumrecords.c_str());
	// Re-read converted file as spectrumrecords file
	if (!spectra.ReadSpectrumRecords(spectrumrecords, &spectrum_header)) {
	carp(CARP_DEBUG, "Deleting %s", spectrumrecords.c_str());
	remove(spectrumrecords.c_str());
	carp(CARP_FATAL, "Error reading spectra file %s", spectrumrecords.c_str());
	}
	}
	input_sr.push_back(InputFile(*f, spectrumrecords, keepSpectrumrecords));
	*/


	//Set up active peptide queue. May need to set up per spectra? Or is each one grabbing a spectra? Figure it out
	//in the code for search

	//problem with the APQ here?

	/*
	vector<ActivePeptideQueue*> active_peptide_queue;
	for (int i = 0; i < NUM_THREADS; i++) {
		active_peptide_queue.push_back(new ActivePeptideQueue(peptide_reader[i]->Reader(), proteins));
		active_peptide_queue[i]->SetBinSize(bin_width, bin_offset);
	}
	*/

	carp(CARP_INFO, "Set up peptide queue");

	//can we add to this collection by hand? Or not?
	//SpectrumCollection spectra;

	// read oms-src, for now (we will use similar stuff in the future but not yet). build up a single spectra spectrum collection
	// from that, modeled after functions called heretofore, then do our version of search on that singular collection. Or single
	// spectra if need be (pulling that OUT of search function)

	//build spectra and then do search
	

	carp(CARP_INFO, "Re-parsed succefully");

	//printf("OMS eating spectra from %s", oms.c_str());

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
	double precursor_window = window;

	if (!peptide_reader[0]) {
		for (int j = 0; j < NUM_THREADS; j++) {
			peptide_reader[j] = new HeadedRecordReader(peptides_file, &peptides_header);
		}
	}

	
	
	//for (int k = 0; k < test_index->pep_array.size(); k++)
	//	carp(CARP_INFO, "Mass for index member %d is %f", k , test_index->pep_array[k].mass());


	//delete test_index;

	
	if (!peptide_reader[0]) {
		for (int j = 0; j < NUM_THREADS; j++) {
			peptide_reader[j] = new HeadedRecordReader(peptides_file, &peptides_header);
		}
	}

	vector<ActivePeptideQueue2*> active_peptide_queue;
	//vector<ActivePeptideQueue*> active_peptide_queue;
	for (int j = 0; j < NUM_THREADS; j++) {
		//active_peptide_queue.push_back(new ActivePeptideQueue(peptide_reader[j]->Reader(), proteins));
		active_peptide_queue.push_back(new ActivePeptideQueue2(test_index, peptide_reader[j]->Reader(), proteins));
		active_peptide_queue[j]->SetBinSize(bin_width, bin_offset);
	}

	// Start the timer.
	wall_clock();
	

	//stuff for sorting the spectra
	
	//vector<OpenMS::MSSpectrum<>> testa;
	/*
	vector<specStruct> specs;

	for (int i = 0; i < msExperimentProfile.getNrSpectra(); i++)
	{
		s = msExperimentProfile.getSpectrum(i);
		if (s.getMSLevel() == 2)
		{
			specs.push_back(specStruct((double)s.getPrecursors()[0].getUnchargedMass(), i));
		}
		//testa.push_back(msExperimentProfile.getSpectrum(i));
	}

	//for (int i = 0; i < testa.size(); i++)
	//	specs.push_back(&(testa[i]));

	std::sort(specs.begin(), specs.end());
	*/

	//for (int i = 0; i < specs.size(); i++)

	int db_hits = 0;
	int ms2_count = 0;
	int decoy_count = 0;
	int target_count = 0;

	for (int i = spectrum_it; i < msExperimentProfile.getNrSpectra(); i++)
	{

		//non-sorted version
		s = msExperimentProfile.getSpectrum(i);
		
		//sorted version
		//s = msExperimentProfile.getSpectrum(specs[i].index);

		if (s.getMSLevel() == 2)
		{
			++ms2_count;

			/*
			if (!peptide_reader[0]) {
				for (int j = 0; j < NUM_THREADS; j++) {
					peptide_reader[j] = new HeadedRecordReader(peptides_file, &peptides_header);
				}
			}

			vector<ActivePeptideQueue*> active_peptide_queue;
			for (int j = 0; j < NUM_THREADS; j++) {
				active_peptide_queue.push_back(new ActivePeptideQueue(peptide_reader[j]->Reader(), proteins));
				//active_peptide_queue.push_back(new ActivePeptideQueue2(test_index, peptide_reader[j]->Reader(), proteins));
				active_peptide_queue[j]->SetBinSize(bin_width, bin_offset);
			}
			*/

			natID = std::istringstream(s.getNativeID());
			natID >> junk >> junk >> scan;
			scan = scan.substr(scan.find('=') + 1);

			Spectrum sspec = Spectrum(atoi(scan.c_str()), s.getPrecursors()[0].getMZ());
			for (OpenMS::MSSpectrum<>::ConstIterator it = s.begin(); it != s.end(); ++it)
			{
				sspec.AddPeak(it->getPos(), it->getIntensity());
			}

			highest_mz = (s.end() - 1)->getPos();
			//MaxBin::SetGlobalMax(highest_mz);
			vector<SpectrumCollection::SpecCharge> spec;
			vector<SpectrumCollection::SpecCharge>* spec_charges = &spec;

			spec_charges->push_back(SpectrumCollection::SpecCharge(s.getPrecursors()[0].getUnchargedMass(), s.getPrecursors()[0].getCharge(), &sspec, 0));

			//need to do something different if we want a per spectrum value for time
			//carp(CARP_INFO, "Elapsed time per spectrum %d conversion: %.3g s", i, wall_clock() /  (1e6*msExperimentProfile.getNrSpectra()) );

			//spectrum is set up. do the search and figure out how to get results in memory nice and clean
			int elution_window = Params::GetInt("elution-window-size");
			bool peptide_centric = Params::GetBool("peptide-centric-search");
			bool use_neutral_loss_peaks = Params::GetBool("use-neutral-loss-peaks");
			bool use_flanking_peaks = Params::GetBool("use-flanking-peaks");
			int max_charge = Params::GetInt("max-precursor-charge");
			ObservedPeakSet observed(bin_width, bin_offset, use_neutral_loss_peaks, use_flanking_peaks);

			if (peptide_centric == false) {
				elution_window = 0;
			}

			for (int j = 0; j < NUM_THREADS; j++) {
				active_peptide_queue[j]->setElutionWindow(elution_window);
				active_peptide_queue[j]->setPeptideCentric(peptide_centric);
			}

			if (elution_window > 0 && elution_window % 2 == 0) {
				for (int j = 0; j < NUM_THREADS; j++) {
					active_peptide_queue[j]->setElutionWindow(elution_window + 1);
				}
			}

			if (!peptide_centric || !exact_pval_search) {
				for (int j = 0; j < NUM_THREADS; j++) {
					active_peptide_queue[j]->setElutionWindow(0);
				}
			}

			//why does the queue need outputs? will this be used later?
			for (int j = 0; j < NUM_THREADS; j++) {
				active_peptide_queue[j]->SetOutputs(NULL, &locations, top_matches, compute_sp, target_file, decoy_file, highest_mz);
			}

			//carp(CARP_INFO, "Iterating on spec_charges, with size %d", spec_charges-> size());
			//carp(CARP_INFO, "Spec goes from %d to %d", spec_charges->begin(), spec_charges->end() );

			/*for (int i = 0; i < 10; i++)
			while (sc < spec_charges->end())
			{
			carp(CARP_INFO, "Spec iteration test at: %d", sc);
			sc++;
			}
			*/

			for (vector<SpectrumCollection::SpecCharge>::const_iterator sc = spec_charges->begin();
				sc < spec_charges->end();
				sc++)
			{

				Spectrum* spectrum = sc->spectrum;
				double precursor_mz = spectrum->PrecursorMZ();
				int charge = sc->charge;
				int scan_num = spectrum->SpectrumNumber();

				//there's some shit about spectrum flags here for some reason
				if (precursor_mz < spectrum_min_mz || precursor_mz > spectrum_max_mz ||
					scan_num < min_scan || scan_num > max_scan ||
					spectrum->Size() < min_peaks ||
					(charge_to_search != 0 && charge != charge_to_search) || charge > max_charge) {
					carp(CARP_INFO, "Exited search early on spectrum %d", i);
					if (spectrum->Size() < min_peaks)
						carp(CARP_INFO, "Not Enough Peaks");
					continue;
				}

				vector<double>* min_mass = new vector<double>();
				vector<double>* max_mass = new vector<double>();
				vector<bool>* candidatePeptideStatus = new vector<bool>();
				double min_range, max_range;
				double *min_ran = &min_range;
				double *max_ran = &max_range;

				computeWindow(*sc, window_type, precursor_window, max_charge, negative_isotope_errors, min_mass, max_mass, &min_range, &max_range);

				if (!exact_pval_search) {  //execute original tide-search program

					//verbose output of search ranges
					/*
					for (vector<double>::const_iterator di = min_mass->begin(); di != min_mass->end(); di++)
					{
						carp(CARP_INFO, "min mass at %f", *di);
					}

					for (vector<double>::const_iterator di = max_mass->begin(); di != max_mass->end(); di++)
					{
						carp(CARP_INFO, "max mass at %f", *di);
					}

					carp(CARP_INFO, "min range at %f", min_range);
					carp(CARP_INFO, "max range at %f", min_range);
					*/




					// Normalize the observed spectrum and compute the cache of
					// frequently-needed values for taking dot products with theoretical
					// spectra.
					observed.PreprocessSpectrum(*spectrum, charge);
					int nCandPeptide = active_peptide_queue[0]->SetActiveRange(min_mass, max_mass, min_range, max_range, candidatePeptideStatus);
					if (nCandPeptide == 0) {
						delete min_mass;
						delete max_mass;
						delete candidatePeptideStatus;
						//carp(CARP_INFO, "SetActiveRange on APQ failed, stopping search");
						continue;
					}

					int candidatePeptideStatusSize = candidatePeptideStatus->size();
					TideMatchSet::Arr2 match_arr2(candidatePeptideStatusSize); // Scored peptides will go here.


					for (int j = 1; j <= candidatePeptideStatusSize; j++)
					{
						const Peptide* tempo = active_peptide_queue[0]->GetPeptide(j);
						double cand_mass = tempo->Mass();
						std::string cand_name = tempo->Seq();

						//carp(CARP_INFO, "Candidate at %d, %s has Mass %f and program at %d" , j, cand_name.c_str(), cand_mass, tempo->Prog(charge) );
					}

					//carp(CARP_INFO, "Peaks for observed spectra:");
					//observed.ShowPeaks();
					//observed.ShowCache();

					// Programs for taking the dot-product with the observed spectrum are laid
					// out in memory managed by the active_peptide_queue, one program for each
					// candidate peptide. The programs will store the results directly into
					// match_arr. We now pass control to those programs.

					//replicate collect scores compiled.  active queue will be active_peptide_queue[0]
					//will use match_arr2 rather than match_arr, others stay same name. candidatePeptidestatussize is called queue_size

					collectScoresCompiled(active_peptide_queue[0], spectrum, observed, &match_arr2,
						candidatePeptideStatusSize, charge);


					// matches will arrange the results in a heap by score, return the top
					// few, and recover the association between counter and peptide. We output
					// the top matches.
					if (peptide_centric) {
						carp(CARP_INFO, "Peptide Centric Search, iterating matches");
						deque<Peptide*>::const_iterator iter_ = active_peptide_queue[0]->iter_;
						TideMatchSet::Arr2::iterator it = match_arr2.begin();
						for (; it != match_arr2.end(); ++iter_, ++it) {
							int peptide_idx = candidatePeptideStatusSize - (it->second);
							if ((*candidatePeptideStatus)[peptide_idx]) {
								(*iter_)->AddHit(spectrum, it->first, 0.0, it->second, charge);
							}
						}
					}
					else {  //spectrum centric match report.
						TideMatchSet::Arr match_arr(nCandPeptide);
						for (TideMatchSet::Arr2::iterator it = match_arr2.begin();
							it != match_arr2.end();
							++it) {
							int peptide_idx = candidatePeptideStatusSize - (it->second);
							if ((*candidatePeptideStatus)[peptide_idx]) {
								TideMatchSet::Scores curScore;
								curScore.xcorr_score = (double)(it->first / XCORR_SCALING);
								curScore.rank = it->second;
								match_arr.push_back(curScore);
								//if ((double)(it->first / XCORR_SCALING) > 1.9f )
								//carp(CARP_INFO, "Candidate with rank %d has XCORR %f", it->second, (double)(it->first / XCORR_SCALING));
							}
						}


						//
						//
						//	Here down is match set stuff
						//
						//

						//TideMatchSet::Arr* matches = &match_arr;

						//TideMatchSet matches(&match_arr, highest_mz);
						//matches.exact_pval_search_ = exact_pval_search;

						//true is highScoreBest, exact_pval_search_ is same
						//top_n is top_matches, peptides is active peptide queue[0] probably
						// matches_ = matches , max_mz_ = highest_mz , exact_pval_search_ = exact_pval_search, elution_window_ = elution_window
						//files stuff needs to be changed/fixed rather than writing out

						/*
						*
						*	full ass all wires pulled out version
						*
						*/

						//set up App and do everything but the write to file

						//fix this later, right now other stuff is broken
						
						double max_corr = 0.0;
						int max_corr_rank = 0;
						int precision = Params::GetInt("precision");
						bool decoy = false;

						carp(CARP_INFO, "XCORR Results for spectrum %d", spectrum->SpectrumNumber());
						//carp(CARP_INFO, "Tide MatchSet reporting %d matches", match_arr.size());

						for (TideMatchSet::Arr::iterator it = match_arr.begin();
							it != match_arr.end();
							++it)
						{

							//carp(CARP_INFO, "Candidate with rank %d has XCORR: %f", it->rank, it->xcorr_score);

							if (it->xcorr_score > max_corr)
							{
								max_corr = it->xcorr_score;
								max_corr_rank = it->rank;
								
							}

							//const Peptide& peptid = *(active_peptide_queue[0]->GetPeptide(max_corr_rank ));
							//if (peptid.IsDecoy())
								decoy = true;

						}

						//const Peptide* temp_pep = active_peptide_queue[0]->GetPeptide(max_corr_rank);
						//carp(CARP_INFO, "Top candidate peptide has sequence: %s", temp_pep->Seq() );
						//if (temp_pep->IsDecoy())
						//	decoy = true;

						//carp(CARP_INFO, "XCORR Results for spectrum %d", spectrum->SpectrumNumber());
						//carp(CARP_INFO, "Tide MatchSet reporting top %d of %d matches", top_matches, match_arr.size());
						double threshold;
						switch (charge)
						{
						case 1:
							threshold = 1.9;
							break;
						case 2:
							threshold = 2.2;
							break;
						default:
							threshold = 3.75;
							break;

						}

						if (max_corr >= threshold)
						{
							carp(CARP_INFO, "Candidate with rank %d has max XCORR, and a hit, with: %f", max_corr_rank, max_corr);
							db_hits++;
							if (decoy)
								++decoy_count;
							else
								++target_count;
						}

						TideMatchSet::Arr* matches_ = &match_arr;

						

						TideMatchSet matches(&match_arr, highest_mz);
						matches.exact_pval_search_ = exact_pval_search;

					

						//report will not work with ActivePeptideQueue2
						//*********************************************
						//*********************************************
						//*********************************************
						//*********************************************
						/*
						boost::mutex * rwlock = new boost::mutex();

						const string spectrum_filename = "HELA_2017-04-27_294.mzML";

						matches.report(target_file, decoy_file, top_matches, spectrum_filename,
							spectrum, charge, active_peptide_queue[0], proteins,
							locations, compute_sp, true, rwlock);

						delete rwlock;
						*/

						if (matches_->size() == 0) {
							carp(CARP_INFO, "Matches empty, moving to next spectra");

							delete min_mass;
							delete max_mass;
							delete candidatePeptideStatus;
							continue;
						}
						
						//verbose reports per spectra

						//carp(CARP_INFO, "Searching spectrum %d", spectrum->SpectrumNumber());

						//carp(CARP_INFO, "Tide MatchSet reporting top %d of %d matches",
						//	top_matches, matches_->size());



						//matches.report(target_file, decoy_file, top_matches, spectrum_filename,
						//spectrum, charge, active_peptide_queue, proteins,
						//locations, compute_sp, true, locks_array[0]);
						

						//
						//
						//	END Match set report stuff
						//
						//


					}  //end peptide_centric == true
				}	//end !exact_pval_search

				/*


				Extremeley long code for exact PVAL search would go here
				WE DO NOT DO EXACT P VAL



				*/


				delete min_mass;
				delete max_mass;
				delete candidatePeptideStatus;

			} //end inner "search loop" over spec charge pairs

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

			
			//delete is here for APQ1
			/*
			for (int i = 0; i < NUM_THREADS; i++) {
				delete active_peptide_queue[i];
				delete peptide_reader[i];
				peptide_reader[i] = NULL;
			}
			*/

		} //end search for spectra that are MS2

	} //end individual spectra loop

	//delete is here for APQ2
	for (int i = 0; i < NUM_THREADS; i++) {
		delete active_peptide_queue[i];
		delete peptide_reader[i];
		peptide_reader[i] = NULL;
	}

	carp(CARP_INFO, "Elapsed time per spectrum conversion: %.3g s", wall_clock() / (1e6*msExperimentProfile.getNrSpectra()));
	carp(CARP_INFO, "There are %d MS2 spectra in the input file.", ms2_count);
	carp(CARP_INFO, "Of those spectra, %d were hits in the search for the current criteria", db_hits);
	carp(CARP_INFO, "Of the hits %d were decoys and %d were targets", decoy_count, target_count);

	//this is the multifile system, we will be modifying to work on either single file or single spectra
	/*
	for (vector<string>::const_iterator f = input_files.begin(); f != input_files.end(); f++)
	{
	SpectrumCollection spectra;
	pb::Header spectrum_header;
	string spectrumrecords = *f;
	bool keepSpectrumrecords = true;
	if (!spectra.ReadSpectrumRecords(spectrumrecords, &spectrum_header)) {
	// Failed, try converting to spectrumrecords file
	carp(CARP_INFO, "Converting %s to spectrumrecords format", f->c_str());
	carp(CARP_INFO, "Elapsed time starting conversion: %.3g s", wall_clock() / 1e6);
	spectrumrecords = Params::GetString("store-spectra");
	keepSpectrumrecords = !spectrumrecords.empty();
	if (!keepSpectrumrecords) {
	spectrumrecords = make_file_path(FileUtils::BaseName(*f) + ".spectrumrecords.tmp");
	}
	else if (input_files.size() > 1) {
	carp(CARP_FATAL, "Cannot use store-spectra option with multiple input "
	"spectrum files");
	}
	carp(CARP_DEBUG, "New spectrumrecords filename: %s", spectrumrecords.c_str());

	//
	//
	//This SpectrumRecordWriter::convert function is the one we must replicate the behavior of
	//
	//
	if (!SpectrumRecordWriter::convert(*f, spectrumrecords)) {
	carp(CARP_FATAL, "Error converting %s to spectrumrecords format", f->c_str());
	}
	carp(CARP_DEBUG, "Reading converted spectra file %s", spectrumrecords.c_str());
	// Re-read converted file as spectrumrecords file
	if (!spectra.ReadSpectrumRecords(spectrumrecords, &spectrum_header)) {
	carp(CARP_DEBUG, "Deleting %s", spectrumrecords.c_str());
	remove(spectrumrecords.c_str());
	carp(CARP_FATAL, "Error reading spectra file %s", spectrumrecords.c_str());
	}
	}
	input_sr.push_back(InputFile(*f, spectrumrecords, keepSpectrumrecords));
	}
	*/
	/*
	for (vector<InputFile>::const_iterator f = input_sr.begin();
	f != input_sr.end();
	f++)
	{
	string spectra_file = f->SpectrumRecords;
	if (!f->Keep) {
	carp(CARP_DEBUG, "Deleting %s", spectra_file.c_str());
	remove(spectra_file.c_str());
	}
	}
	*/
	//next is loop through the spectra, including the search itself. In the final version this will be a function
	//that is waiting and will be called whenever an MS2 takes place



	//
	//
	//
	//clean up
	//
	//
	//

	/*
	for (int i = 0; i < NUM_THREADS; i++) {
		delete active_peptide_queue[i];
		delete peptide_reader[i];
		peptide_reader[i] = NULL;
	}
	*/

	delete test_index;
	delete negative_isotope_errors;

	for (ProteinVec::iterator i = proteins.begin(); i != proteins.end(); ++i) {
		delete *i;
	}
	if (target_file) {
		delete target_file;
		if (decoy_file) {
			delete decoy_file;
		}
	}
	delete[] aaFreqN;
	delete[] aaFreqI;
	delete[] aaFreqC;
	delete[] aaMass;

	return 0;
}

void collectScoresCompiled(
	//ActivePeptideQueue* active_peptide_queue,
	ActivePeptideQueue2* active_peptide_queue,
	const Spectrum* spectrum,
	const ObservedPeakSet& observed,
	TideMatchSet::Arr2* match_arr,
	int queue_size,
	int charge
	) {
	if (!active_peptide_queue->HasNext()) {
		return;
	}
	// prog gets the address of the dot-product program for the first peptide
	// in the active queue.
	const void* prog = active_peptide_queue->NextPeptide()->Prog(charge);
	const int* cache = observed.GetCache();
	// results will get (score, counter) pairs, where score is the dot product
	// of the observed peak set with a candidate peptide. The candidate
	// peptide is given by counter, which refers to the index within the
	// ActivePeptideQueue, counting from the back. This complication
	// simplifies the generated programs, which now simply dump the counter.
	pair<int, int>* results = match_arr->data();

	// See compiler.h for a description of the programs beginning at prog and
	// how they are generated. Here we initialize certain registers to the
	// values expected by the programs and call the first one (*prog).
	//
	// See gnu assembler format for more on this format. We tell the compiler
	// to set these registers:
	// edx/rdx points to the cache.
	// eax/rax points to the first program.
	// ecx/rcx is the counter and gets the size of the active queue.
	// edi/rdi points to the results buffer.
	//
	// The push and pop operations are a workaround for a compiler that
	// doesn't understand that %ecx and %edi (or %rcx and %rdi) get
	// clobbered. Since they're already input registers, they can't be
	// included in the clobber list.

#ifdef _MSC_VER
#ifdef _WIN64
	DWORD64 rcx;
	DWORD64 rdi;

	volatile bool restored = false;
	CONTEXT context;
	RtlCaptureContext(&context);
	if (!restored) {
		rcx = context.Rcx;
		rdi = context.Rdi;

		context.Rdx = (DWORD64)cache;
		context.Rax = (DWORD64)prog;
		context.Rcx = (DWORD64)queue_size;
		context.Rdi = (DWORD64)results;

		restored = true;
		RtlRestoreContext(&context, NULL);
	}
	else {
		((void(*)(void))prog)();
	}

	restored = false;
	RtlCaptureContext(&context);
	if (!restored) {
		context.Rcx = rcx;
		context.Rdi = rdi;

		restored = true;
		RtlRestoreContext(&context, NULL);
	}
#else
	__asm {
		cld
			push ecx
			push edi
			mov edx, cache
			mov eax, prog
			mov ecx, queue_size
			mov edi, results
			call eax
			pop edi
			pop ecx
	}
#endif
#else
	__asm__ __volatile__("cld\n" // stos operations increment edi
#ifdef __x86_64__
		"push %%rcx\n"
		"push %%rdi\n"
		"call *%%rax\n"
		"pop %%rdi\n"
		"pop %%rcx\n"
#else
		"push %%ecx\n"
		"push %%edi\n"
		"call *%%eax\n"
		"pop %%edi\n"
		"pop %%ecx\n"
#endif
		: // no outputs
	: "d" (cache),
		"a" (prog),
		"c" (queue_size),
		"D" (results)
		);
#endif

	// match_arr is filled by the compiled programs, not by calls to
	// push_back(). We have to set the final size explicitly.
	match_arr->set_size(queue_size);
}

void computeWindow(
	const SpectrumCollection::SpecCharge& sc,
	WINDOW_TYPE_T window_type,
	double precursor_window,
	int max_charge,
	vector<int>* negative_isotope_errors,
	vector<double>* out_min,
	vector<double>* out_max,
	double* min_range,
	double* max_range
	) {

	double unit_dalton = BIN_WIDTH;

	switch (window_type) {
	case WINDOW_MASS:
		for (vector<int>::const_iterator ie = negative_isotope_errors->begin(); ie != negative_isotope_errors->end(); ++ie) {
			out_min->push_back(sc.neutral_mass + (*ie * unit_dalton) - precursor_window);
			out_max->push_back(sc.neutral_mass + (*ie * unit_dalton) + precursor_window);
		}
		*min_range = (sc.neutral_mass + (negative_isotope_errors->front() * unit_dalton)) - precursor_window;
		*max_range = (sc.neutral_mass + (negative_isotope_errors->back() * unit_dalton)) + precursor_window;
		break;
	case WINDOW_MZ: {
		double mz_minus_proton = sc.spectrum->PrecursorMZ() - MASS_PROTON;
		for (vector<int>::const_iterator ie = negative_isotope_errors->begin(); ie != negative_isotope_errors->end(); ++ie) {
			out_min->push_back((mz_minus_proton - precursor_window) * sc.charge + (*ie * unit_dalton));
			out_max->push_back((mz_minus_proton + precursor_window) * sc.charge + (*ie * unit_dalton));
		}
		*min_range = (mz_minus_proton*sc.charge + (negative_isotope_errors->front() * unit_dalton)) - precursor_window*max_charge;
		*max_range = (mz_minus_proton*sc.charge + (negative_isotope_errors->back() * unit_dalton)) + precursor_window*max_charge;
		break;
	}
	case WINDOW_PPM: {
		double tiny_precursor = precursor_window * 1e-6;
		for (vector<int>::const_iterator ie = negative_isotope_errors->begin(); ie != negative_isotope_errors->end(); ++ie) {
			out_min->push_back((sc.neutral_mass + (*ie * unit_dalton)) * (1.0 - tiny_precursor));
			out_max->push_back((sc.neutral_mass + (*ie * unit_dalton)) * (1.0 + tiny_precursor));
		}
		*min_range = (sc.neutral_mass + (negative_isotope_errors->front() * unit_dalton)) * (1.0 - tiny_precursor);
		*max_range = (sc.neutral_mass + (negative_isotope_errors->back() * unit_dalton)) * (1.0 + tiny_precursor);
		break;
	}
	default:
		carp(CARP_FATAL, "Invalid window type");
	}
	carp(CARP_DETAILED_DEBUG, "Scan %d.%d mass window is [%f, %f]",
		sc.spectrum->SpectrumNumber(), sc.charge, (*out_min)[0], (*out_max)[0]);
}

	/*
	*
	*
	*
	*	HERE DOWN IS CUSTOM CODE
	*
	*
	*
	*/


	//carp(CARP_DEBUG, "Number of Threads: %d", NUM_THREADS); // prints the number of threads
	//printf("Number of Threads: %d\n", NUM_THREADS);
	
	//we should open and deal with the index files here before doing anything else. Human index is only about 100MB


		/*
		*
		*
		*	Everything below is a useful model for going through individual spectra, which we will eventually pass to the above
		*
		*
		*

		OpenMS::MSSpectrum<> s = msExperimentProfile.getSpectrum(spectrum);
		std::string junk;
		std::string scan;

		std::istringstream natID = std::istringstream( s.getNativeID() );

		natID >> junk >> junk >> scan;

		scan = scan.substr(scan.find('=') + 1);

		while (s.getMSLevel() == 1)
		{
		spectrum++;
		s = msExperimentProfile.getSpectrum(spectrum);
		}

		out << "H	CreationDate	10 / 4 / 2017 3:02 : 18 PM" << std::endl;
		out << "H	Extractor	TSMS2" << std::endl;
		out << "H	ExtractorVersion	1.0" << std::endl;
		out << "H	Comments	TSMS2 written by Trent Stohrer, 2017" << std::endl;
		out << "H	ExtractorOptions	MS2 / MS1" << std::endl;
		out << "S \t" << scan << "\t" << scan << "\t" << s.getPrecursors()[0].getMZ() << std::endl;
		out << "Z \t" << s.getPrecursors()[0].getCharge() << "\t" << (s.getPrecursors()[0].getUnchargedMass() + OpenMS::Constants::PROTON_MASS_U) << std::endl;
		for (OpenMS::MSSpectrum<>::ConstIterator it = s.begin(); it != s.end(); ++it)
		{
		out << it->getPos() << "\t" << it->getIntensity() << std::endl;
		}

		out.close();
		std::string write_file = "";

		for (int i = spectrum + 1; i < msExperimentProfile.getNrSpectra(); i++)
		//for (int i = 0; i <= spectrum; i++)
		{
		s = msExperimentProfile.getSpectrum(i);
		if (s.getMSLevel() == 2)
		{
		write_file = "./temp/spec_" + std::to_string(i) + ".ms2";
		out.open(write_file);
		natID = std::istringstream(s.getNativeID());
		natID >> junk >> junk >> scan;
		scan = scan.substr(scan.find('=') + 1);
		out << "S \t" << scan << "\t" << scan << "\t" << s.getPrecursors()[0].getMZ() << std::endl;
		out << "Z \t" << s.getPrecursors()[0].getCharge() << "\t" << (s.getPrecursors()[0].getUnchargedMass() + OpenMS::Constants::PROTON_MASS_U) << std::endl;
		for (OpenMS::MSSpectrum<>::ConstIterator it = s.begin(); it != s.end(); ++it)
		{
		out << it->getPos() << "\t" << it->getIntensity() << std::endl;
		}

		out.close();
		}
		//std::cout << "The " << i << " spectrum has " << s.size() << " peaks and is at MS level: " << s.getMSLevel() << std::endl;
		}

		out.close();

		*/

		//std::cout << "There are " << msExperimentProfile.getNrSpectra() << " spectra in the input file." << std::endl;

		//TideSearchApplication pass;

		//pass.TideSearchApplication::main(argc, argv);

		/*
		 *
		 *
		 *
		 *
		 *
		 *	DENNIS'S STUFF PROVIDING EXAMPLE CODE BELOW HERE
		 *
		 *
		 *
		 *
		 *
		 *
		 */



		//out << "isotope.range" << "\t" << "ion.index" << "\t" << "ion.name" << "\t" << "mz" << "\t" << "int" << std::endl;
		//calc_out << "isotope.range" << "\t" << "ion.index" << "\t" << "ion.name" << "\t" << "mz" << "\t" << "int" << "\t" << "method" << std::endl;

		//const Ion precursorIon = Ion(OpenMS::AASequence::fromString("[-18.010565]ELYENKPRRPYIL"), OpenMS::Residue::Full, 3);

		//create list of b and y ions
		//std::vector<Ion> ionList;
		//Ion::generateFragmentIons(ionList, precursorIon.sequence, precursorIon.charge);

		//double isotopeStep = OpenMS::Constants::NEUTRON_MASS_U / precursorIon.charge;

		//Loop through all spectra
		//for (int specIndex = 0; specIndex < msExperimentCentroid.getNrSpectra(); ++specIndex) {
		//for (int specIndex = 0; specIndex < 10; ++specIndex)
		//{

		//	OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrumCentroid = msExperimentCentroid.getSpectrum(specIndex);
		//	OpenMS::MSSpectrum<OpenMS::Peak1D> currentSpectrumProfile = msExperimentProfile.getSpectrum(specIndex);

		//	if (currentSpectrumCentroid.getMSLevel() == 1) continue;

		//	currentSpectrumCentroid.sortByPosition();

		//const OpenMS::Precursor precursorInfo = currentSpectrumCentroid.getPrecursors()[0];
		//	std::vector<OpenMS::UInt> precursorIsotopes;

		//}

//old code

//compute window replication
/*
double unit_dalton = BIN_WIDTH;

switch (window_type) {
case WINDOW_MASS:
for (vector<int>::const_iterator ie = negative_isotope_errors->begin(); ie != negative_isotope_errors->end(); ++ie) {
min_mass->push_back(sc->neutral_mass + (*ie * unit_dalton) - precursor_window);
max_mass->push_back(sc->neutral_mass + (*ie * unit_dalton) + precursor_window);
}
*min_ran = (sc->neutral_mass + (negative_isotope_errors->front() * unit_dalton)) - precursor_window;
*max_ran = (sc->neutral_mass + (negative_isotope_errors->back() * unit_dalton)) + precursor_window;
break;
case WINDOW_MZ: {
double mz_minus_proton = sc->spectrum->PrecursorMZ() - MASS_PROTON;
for (vector<int>::const_iterator ie = negative_isotope_errors->begin(); ie != negative_isotope_errors->end(); ++ie) {
min_mass->push_back((mz_minus_proton - precursor_window) * sc->charge + (*ie * unit_dalton));
max_mass->push_back((mz_minus_proton + precursor_window) * sc->charge + (*ie * unit_dalton));
}
*min_ran = (mz_minus_proton*sc->charge + (negative_isotope_errors->front() * unit_dalton)) - precursor_window*max_charge;
*max_ran = (mz_minus_proton*sc->charge + (negative_isotope_errors->back() * unit_dalton)) + precursor_window*max_charge;
break;
}
case WINDOW_PPM: {
double tiny_precursor = precursor_window * 1e-6;
for (vector<int>::const_iterator ie = negative_isotope_errors->begin(); ie != negative_isotope_errors->end(); ++ie) {
min_mass->push_back((sc->neutral_mass + (*ie * unit_dalton)) * (1.0 - tiny_precursor));
max_mass->push_back((sc->neutral_mass + (*ie * unit_dalton)) * (1.0 + tiny_precursor));
}
*min_ran = (sc->neutral_mass + (negative_isotope_errors->front() * unit_dalton)) * (1.0 - tiny_precursor);
*max_ran = (sc->neutral_mass + (negative_isotope_errors->back() * unit_dalton)) * (1.0 + tiny_precursor);
break;
}
default:
carp(CARP_FATAL, "Invalid window type");
}
carp(CARP_DETAILED_DEBUG, "Scan %d.%d mass window is [%f, %f]",
sc->spectrum->SpectrumNumber(), sc->charge, (*min_mass)[0], (*max_mass)[0]);
//end compute window replication
*/