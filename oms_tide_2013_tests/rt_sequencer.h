#ifndef MY_RT_SEQUENCER_H
#define MY_RT_SEQUENCER_H

#define BOOST_ALL_NO_LIB

#include "stdafx.h"

#include <io.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
//#include <thread> //maybe not needed
#include <ctime>
#include <map>
#include <vector>
#include <gflags/gflags.h>
#include <ratio>
#include <windows.h>
#include <queue>
#include <functional>

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

#include <OpenMS/CHEMISTRY/AASequence.h>

#include "peptides.pb.h"
#include "spectrum.pb.h"
#include "app/tide/theoretical_peak_set.h"
#include "app/tide/max_mz.h"

#include "mem_peptide_queue.h"
#include "inmemindex.h"

const double XCORR_SCALING = 100000000.0;

class RT_Sequencer
{
public:
	RT_Sequencer();
	void init(const char* index_name = NULL, const char* param_name = NULL, const char* inclusion_name = NULL, double xc1 = 1.8, double xc2 = 2.2, double xc3 = 3.5, double dcn = 0.08, int window_type_index = 1, double window_size = 3.0);
	~RT_Sequencer();
	bool is_match(double high_mz, double precursor_mass, int precursor_charge);
	void addPeak(double mz, double intensity);
	void makeSpec(int scan_num, double prec_mz);
	double score();
	std::string sequence() { return open_ms_seq.toString(); }
	std::vector<double> get_sps() { return sps_targets; }

private:
	bool protInclude(std::string prot_name);
	void computeFragments(int charge);
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

private:
	ActivePeptideQueue2* peptide_queue;
	InMemIndex* pep_index;
	Spectrum* curr_spectrum;
	ProteinVec proteins;
	WINDOW_TYPE_T window_type;
	double bin_width;
	double bin_offset;
	bool compute_sp;
	vector<int>* negative_isotope_errors;
	bool HAS_DECOYS;
	bool PROTEIN_LEVEL_DECOYS;
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
	double last_score;
	std::vector<double> last_fragments;
	std::vector<double> sps_targets;
	std::map<std::string, std::string> inclusion_list;
	std::string last_seq;
	OpenMS::AASequence open_ms_seq;
	bool static_mods;
	std::string static_delta;

	int low_peak;

};

class fragmentNode
{
public:
	double mz;
	int list;

	friend bool operator>(const fragmentNode& l, const fragmentNode& r) { return ((l.mz > r.mz) ? true : false); }
	fragmentNode(double mymz, int mylist) : mz(mymz), list(mylist) {}

};

class peakNode
{
public:
	double mz;
	double intensity;

	friend bool operator>(const peakNode& l, const peakNode& r) { return ((l.intensity > r.intensity) ? true : false); }
	peakNode(double mymz, double myintense) : mz(mymz), intensity(myintense) {}
};

#endif

