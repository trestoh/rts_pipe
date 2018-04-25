#ifndef MY_RT_SEQUENCER_H
#define MY_RT_SEQUENCER_H

#define BOOST_ALL_NO_LIB

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

#include "mem_peptide_queue.h"
#include "inmemindex.h"

const double XCORR_SCALING = 100000000.0;

class RT_Sequencer
{
public:
	RT_Sequencer();
	void init(const char* index_name = NULL, double xc1 = 1.8, double xc2 = 2.2, double xc3 = 3.5, double dcn = 0.08);
	~RT_Sequencer();
	bool is_match(Spectrum& sspec, double high_mz, double precursor_mass, OpenMS::Int precursor_charge, string& peptide_hit, bool& decoy);
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

#endif

