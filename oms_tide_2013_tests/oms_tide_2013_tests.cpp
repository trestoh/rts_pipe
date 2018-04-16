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
#include <app/tide/records_to_vector-inl.h>
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

#include "peptides.pb.h"
#include "spectrum.pb.h"
#include "app/tide/theoretical_peak_set.h"
#include "app/tide/max_mz.h"

#include "inmemindex.h"
#include "mem_peptide_queue.h"
//#include "rt_sequencer.h"

const double XCORR_SCALING = 100000000.0;
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

void collectScoresCompiled(
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

class RT_Sequencer
{
public:
	RT_Sequencer();
	void init(const char* index_name = NULL, double xc1 = 1.8, double xc2 = 2.2, double xc3 = 3.5, double dcn = 0.08);
	~RT_Sequencer();
	bool is_match(Spectrum& sspec, double high_mz, double precursor_mass, OpenMS::Int precursor_charge, string& peptide_hit, bool& decoy);

private:
	ActivePeptideQueue2* peptide_queue;
	InMemIndex* pep_index;
	double window;
	WINDOW_TYPE_T window_type;
	double bin_width;
	double bin_offset;
	bool compute_sp;
	vector<int>* negative_isotope_errors;
	double* aaFreqN;
	double* aaFreqI;
	double* aaFreqC;
	int* aaMass;
	int nAA;
	bool HAS_DECOYS;
	bool PROTEIN_LEVEL_DECOYS;
	bool overwrite;
	bool exact_pval_search;
	double spectrum_min_mz;
	double spectrum_max_mz;
	int min_peaks;
	int top_matches;
	double highest_mz;
	double precursor_window;
	int min_scan;
	int max_scan;
	int charge_to_search;
	double xcorr1;
	double xcorr2;
	double xcorr3;
	double deltacn;

};

RT_Sequencer::RT_Sequencer()
{
	pep_index = NULL;
	peptide_queue = NULL;
	negative_isotope_errors = NULL;
	aaFreqN = NULL;
	aaFreqI = NULL;
	aaFreqC = NULL;
	aaMass = NULL;
}

void RT_Sequencer::init(const char* index_name, double xc1, double xc2, double xc3, double dcn)
{
	xcorr1 = xc1;
	xcorr2 = xc2;
	xcorr3 = xc3;
	deltacn = dcn;

	HAS_DECOYS = false;
	PROTEIN_LEVEL_DECOYS = false;

	if (index_name == NULL)
	{
		carp(CARP_FATAL, "Error Reading Index, no name provided");
		exit(1);
	}

	std::string peptides_file = FileUtils::Join(index_name, "pepix");
	std::string proteins_file = FileUtils::Join(index_name, "protix");
	std::string auxlocs_file = FileUtils::Join(index_name, "auxlocs");

	pb::Header pepts_header;
	HeadedRecordReader pept_reader(peptides_file, &pepts_header);

	if ((pepts_header.file_type() != pb::Header::PEPTIDES) || !pepts_header.has_peptides_header())
	{
		carp(CARP_FATAL, "Error Reading Index, Peptide Header not found");
		exit(1);
	}

	const pb::Header::PeptidesHeader& pepsHeader = pepts_header.peptides_header();
	Params::Set("enzyme", pepsHeader.enzyme());
	const char* digestString = digest_type_to_string(pepsHeader.full_digestion() ? FULL_DIGEST : PARTIAL_DIGEST);
	Params::Set("digestion", digestString);
	Params::Set("isotopic-mass", pepsHeader.monoisotopic_precursor() ? "mono" : "average");

	Params::Finalize();
	GlobalParams::set();

	window = Params::GetDouble("precursor-window");
	window_type = string_to_window_type(Params::GetString("precursor-window-type"));

	string charge_string = Params::GetString("spectrum-charge");
	if (charge_string == "all"){
		carp(CARP_DEBUG, "Searching all charge states");
		charge_to_search = 0;
	}

	else {
		charge_to_search = atoi(charge_string.c_str());
		if (charge_to_search < 1 || charge_to_search > 6) {
			carp(CARP_FATAL, "Invalid  spec-charge state");
		}
		carp(CARP_DEBUG, "Searching charge state %d", charge_to_search);
	}

	string scan_range = Params::GetString("scan-number");
	if (scan_range.empty()) {
		min_scan = 0;
		max_scan = BILLION;
		carp(CARP_DEBUG, "Searching all scans");
	}
	else if (scan_range.find('-') == string::npos) {
		min_scan = max_scan = atoi(scan_range.c_str());
		carp(CARP_DEBUG, "Searching single scan %d", min_scan);
	}
	else {
		if (!get_range_from_string(scan_range.c_str(), min_scan, max_scan))
			carp(CARP_FATAL, "Scan range invalid, must be of form <first>-<last>");
		else {
			if (min_scan > max_scan) {
				int tmp_scan = min_scan;
				min_scan = max_scan;
				max_scan = tmp_scan;
			}
			carp(CARP_DEBUG, "Searching scan range %d-%d", min_scan, max_scan);
		}
	}

	exact_pval_search = Params::GetBool("exact-p-value");
	bin_width = Params::GetDouble("mz-bin-width");
	bin_offset = Params::GetDouble("mz-bin-offset");

	bool concat = Params::GetBool("concat");

	compute_sp = Params::GetBool("compute-sp");
	if (Params::GetBool("sqt-output") && !compute_sp) {
		compute_sp = true;
		carp(CARP_INFO, "Enabling parameter compute-sp since SQT output is enabled, will increase runtime");
	}

	string isotope_errors_string = Params::GetString("isotope-error");
	if (isotope_errors_string[0] == ',') {
		carp(CARP_FATAL, "Error in isotope_error parameter formatting: (%s) ", isotope_errors_string.c_str());
	}
	for (int i = 0; isotope_errors_string[i] != '\0'; ++i){
		if (isotope_errors_string[i] == ',' && (isotope_errors_string[i + 1] == ',' || isotope_errors_string[i + 1] == '\0'))
			carp(CARP_FATAL, "Error in isotope_error parameter formatting: (%s) ", isotope_errors_string.c_str());
	}

	negative_isotope_errors = new vector<int>();
	negative_isotope_errors->push_back(0);

	if (isotope_errors_string != "") {
		vector<int> isotope_errors = StringUtils::Split<int>(Params::GetString("isotope-error"), ',');
		for (vector<int>::iterator it = isotope_errors.begin(); it != isotope_errors.end(); ++it){
			if (*it < 0)
				carp(CARP_FATAL, "Found a negative isotope error, don't do this, model with a modification.");
			else if (find(negative_isotope_errors->begin(), negative_isotope_errors->end(), -1 * *it) != negative_isotope_errors->end())
				carp(CARP_FATAL, "Found duplicate when parsing isotope_error parameter");
			negative_isotope_errors->push_back(-1 * *it);
		}
	}

	sort(negative_isotope_errors->begin(), negative_isotope_errors->end());

	ProteinVec proteins;
	pb::Header protein_header;

	if (!ReadRecordsToVector<pb::Protein, const pb::Protein>(&proteins, proteins_file, &protein_header))
		carp(CARP_FATAL, "Error Reading Index!");

	aaFreqN = NULL;
	aaFreqI = NULL;
	aaFreqC = NULL;
	aaMass = NULL;
	nAA = 0;

	pb::Header peptides_header;

	HeadedRecordReader* peptide_reader = new HeadedRecordReader(peptides_file, &peptides_header);

	if ((peptides_header.file_type() != pb::Header::PEPTIDES) || !peptides_header.has_peptides_header())
		carp(CARP_FATAL, "Error Reading Index");

	const pb::Header::PeptidesHeader& pepHeader = peptides_header.peptides_header();
	DECOY_TYPE_T headerDecoyType = (DECOY_TYPE_T)pepHeader.decoys();
	if (headerDecoyType != NO_DECOYS){
		HAS_DECOYS = true;
		if (headerDecoyType == PROTEIN_REVERSE_DECOYS)
			PROTEIN_LEVEL_DECOYS = true;
	}

	std::vector<const pb::AuxLocation*> locations;
	if (!ReadRecordsToVector<pb::AuxLocation>(&locations, auxlocs_file))
		carp(CARP_FATAL, "Error Reading Index");
	carp(CARP_INFO, "Read %d auxlocs", locations.size());

	MassConstants::Init(&pepHeader.mods(), &pepHeader.nterm_mods(), &pepHeader.cterm_mods(), bin_width, bin_offset);
	ModificationDefinition::ClearAll();
	TideMatchSet::initModMap(pepHeader.mods(), ANY);
	TideMatchSet::initModMap(pepHeader.nterm_mods(), PEPTIDE_N);
	TideMatchSet::initModMap(pepHeader.cterm_mods(), PEPTIDE_C);

	MaxBin::SetGlobalMax(2000.0);

	pep_index = new InMemIndex(peptide_reader->Reader(), proteins);

	overwrite = Params::GetBool("overwrite");
	stringstream ss;
	ss << Params::GetString("enzyme") << '-' << Params::GetString("digestion");
	TideMatchSet::CleavageType = ss.str();

	spectrum_min_mz = Params::GetDouble("spectrum-min-mz");
	spectrum_max_mz = Params::GetDouble("spectrum-max-mz");
	min_peaks = Params::GetInt("min-peaks");
	top_matches = Params::GetInt("top-match");
	highest_mz = 0.0f;
	precursor_window = window;

	peptide_queue = new ActivePeptideQueue2(pep_index, peptide_reader->Reader(), proteins);
	peptide_queue->SetBinSize(bin_width, bin_offset);

}

RT_Sequencer::~RT_Sequencer()
{
	delete pep_index;
	delete peptide_queue;
	delete[] aaFreqN;
	delete[] aaFreqI;
	delete[] aaFreqC;
	delete[] aaMass;
	delete negative_isotope_errors;
}

bool RT_Sequencer::is_match(Spectrum& sspec, double high_mz, double precursor_mass, OpenMS::Int precursor_charge, std::string& peptide_hit, bool& decoy)
{
	highest_mz = high_mz;
	vector<SpectrumCollection::SpecCharge> spec;
	vector<SpectrumCollection::SpecCharge>* spec_charges = &spec;

	spec_charges->push_back(SpectrumCollection::SpecCharge(precursor_mass, precursor_charge, &sspec, 0));
	int elution_window = Params::GetInt("elution-window-size");
	bool peptide_centric = Params::GetBool("peptide-centric-search");
	bool use_neutral_loss_peaks = Params::GetBool("use-neutral-loss-peaks");
	bool use_flanking_peaks = Params::GetBool("use-flanking-peaks");
	int max_charge = Params::GetInt("max-precursor-charge");
	ObservedPeakSet observed(bin_width, bin_offset, use_neutral_loss_peaks, use_flanking_peaks);

	if (peptide_centric == false) {
		elution_window = 0;
	}

	peptide_queue->setElutionWindow(elution_window);
	peptide_queue->setPeptideCentric(peptide_centric);

	if (elution_window > 0 && elution_window % 2 == 0)
		peptide_queue->setElutionWindow(elution_window + 1);

	if (!peptide_centric || !exact_pval_search)
		peptide_queue->setElutionWindow(0);

	for (vector<SpectrumCollection::SpecCharge>::const_iterator sc = spec_charges->begin(); sc < spec_charges->end(); sc++)
	{
		Spectrum* spectrum = sc->spectrum;
		double precursor_mz = spectrum->PrecursorMZ();
		int charge = sc->charge;
		int scan_num = spectrum->SpectrumNumber();

		if (precursor_mz < spectrum_min_mz || precursor_mz > spectrum_max_mz || scan_num < min_scan || scan_num > max_scan
			|| spectrum->Size() < min_peaks || (charge_to_search != 0 && charge != charge_to_search) || charge > max_charge) {
			carp(CARP_INFO, "Exited search early");
			if (spectrum->Size() < min_peaks)
				carp(CARP_INFO, "Not Enough Peaks.");
			return false;
		}

		vector<double>* min_mass = new vector<double>();
		vector<double>* max_mass = new vector<double>();
		vector<bool>* candidatePeptideStatus = new vector<bool>();
		double min_range, max_range;
		double* min_ran = &min_range;
		double* max_ran = &max_range;

		computeWindow(*sc, window_type, precursor_window, max_charge, negative_isotope_errors, min_mass, max_mass, &min_range, &max_range);

		if (!exact_pval_search) {
			observed.PreprocessSpectrum(*spectrum, charge);
			int nCandPeptide = peptide_queue->SetActiveRange(min_mass, max_mass, min_range, max_range, candidatePeptideStatus);
			if (nCandPeptide == 0) {
				delete min_mass;
				delete max_mass;
				delete candidatePeptideStatus;
				return false;
			}

			int candidatePeptideStatusSize = candidatePeptideStatus->size();
			TideMatchSet::Arr2 match_arr2(candidatePeptideStatusSize);

			collectScoresCompiled(peptide_queue, spectrum, observed, &match_arr2, candidatePeptideStatusSize, charge);

			if (peptide_centric) {
				carp(CARP_INFO, "Peptide Centric Search, iterating matches");
				deque<Peptide*>::const_iterator iter_ = peptide_queue->iter_;
				TideMatchSet::Arr2::iterator it = match_arr2.begin();
				for (; it != match_arr2.end(); ++iter_, ++it) {
					int peptide_idx = candidatePeptideStatusSize - (it->second);
					if ((*candidatePeptideStatus)[peptide_idx]) {
						(*iter_)->AddHit(spectrum, it->first, 0.0, it->second, charge);
					}
				}
			} //end peptide centric

			else {
				TideMatchSet::Arr match_arr(nCandPeptide);
				for (TideMatchSet::Arr2::iterator it = match_arr2.begin(); it != match_arr2.end(); ++it){
					int peptide_idx = candidatePeptideStatusSize - (it->second);
					if ((*candidatePeptideStatus)[peptide_idx]) {
						TideMatchSet::Scores curScore;
						curScore.xcorr_score = (double)(it->first / XCORR_SCALING);
						curScore.rank = it->second;
						match_arr.push_back(curScore);
					}
				}

				double max_corr = match_arr.begin()->xcorr_score;
				double second_corr = match_arr.begin()->xcorr_score;
				int max_corr_rank = 0;
				int precision = Params::GetInt("precision");

				for (TideMatchSet::Arr::iterator it = match_arr.begin(); it != match_arr.end(); ++it) {
					if (it->xcorr_score > max_corr) {
						second_corr = max_corr;
						max_corr = it->xcorr_score;
						max_corr_rank = it->rank;
					}
				}

				double delta_corr = (max_corr - second_corr) / second_corr;
				double threshold;
				switch (charge) {
				case 1:
					threshold = xcorr1;
					break;
				case 2:
					threshold = xcorr2;
					break;
				default:
					threshold = xcorr3;
					break;
				}

				double delta_threshold = deltacn;

				const Peptide& peptid = *(peptide_queue->GetPeptide(max_corr_rank));

				if (peptid.IsDecoy())
					decoy = true;
				else
					decoy = false;

				vector<Crux::Modification> modVector;
				string seq(peptid.Seq());
				const ModCoder::Mod* mods;
				int pep_mods = peptid.Mods(&mods);
				for (int i = 0; i < pep_mods; i++){
					int mod_index;
					double mod_delta;
					MassConstants::DecodeMod(mods[i], &mod_index, &mod_delta);
					const ModificationDefinition* modDef = ModificationDefinition::Find(mod_delta, false);
					if (modDef == NULL)
						continue;
					modVector.push_back(Crux::Modification(modDef, mod_index));
				}

				Crux::Peptide cruxPep = Crux::Peptide(peptid.Seq(), modVector);
					
				peptide_hit = cruxPep.getModifiedSequenceWithMasses();

				if (max_corr >= threshold && delta_corr >= delta_threshold) {
					delete min_mass;
					delete max_mass;
					delete candidatePeptideStatus;
					return true;
				}

			} //end spectrum centric search

		}//end NOT exact p value version of search

		delete min_mass;
		delete max_mass;
		delete candidatePeptideStatus;

	}//end loop through spectrum charges

	return false;

}

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
	RT_Sequencer sequencer;
	sequencer.init(index.c_str(), xcorr1, xcorr2, xcorr3, deltacn);

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

			Spectrum sspec = Spectrum(atoi(scan.c_str()), s.getPrecursors()[0].getMZ());
			for (OpenMS::MSSpectrum<>::ConstIterator it = s.begin(); it != s.end(); ++it)
			{
				sspec.AddPeak(it->getPos(), it->getIntensity());
			}			

			highest_mz = (s.end() - 1)->getPos();

			std::string pep_seq;
			bool isDecoy;

			bool db_match = sequencer.is_match(sspec, highest_mz, s.getPrecursors()[0].getUnchargedMass(), s.getPrecursors()[0].getCharge(), pep_seq, isDecoy);

			if (db_match)
			{
				db_hits++;
				if (isDecoy)
					++decoy_count;
				else {
					++target_count;
					pep_hits.push_back(pep_seq);
				}
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
		try{
			RtlRestoreContext(&context, NULL);
		}

		catch (...)
		{
			carp(CARP_INFO, "Bad Bad Not Good...");
		}


	}
	else {
		try{
			((void(*)(void))prog)();
		}

		catch (...)
		{
			//carp(CARP_INFO, "Bad Bad Not Good...");
			carp(CARP_INFO, "Queue Dump");
			active_peptide_queue->dumpQueue(charge);
		}

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