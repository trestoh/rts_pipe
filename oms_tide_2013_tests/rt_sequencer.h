#ifndef MY_RT_SEQUENCER_H
#define MY_RT_SEQUENCER_H

#include "mem_peptide_queue.h"
#include "inmemindex.h"

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

#endif

