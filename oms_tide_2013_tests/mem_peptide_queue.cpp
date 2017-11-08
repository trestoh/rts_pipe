// original author: Benjamin Diament
// subsequently modified by Attila Kertesz-Farkas, Jeff Howbert
#include "stdafx.h"

#define BOOST_ALL_NO_LIB

#include <deque>
#include <gflags/gflags.h>
#include "app/tide/records.h"
#include "peptides.pb.h"
#include "app/tide/peptide.h"
#include "mem_peptide_queue.h"
#include "app/tide/records_to_vector-inl.h"
#include "app/tide/theoretical_peak_set.h"
#include "app/tide/compiler.h"
#include "app/TideMatchSet.h"
#include "app/tide/active_peptide_queue.h"
#define CHECK(x) GOOGLE_CHECK((x))


ActivePeptideQueue2::ActivePeptideQueue2(InMemIndex* index, RecordReader* reader, const vector<const pb::Protein*>& proteins)
	: reader_(reader),
	proteins_(proteins),
	theoretical_peak_set_(2000),   // probably overkill, but no harm
	theoretical_b_peak_set_(200),  // probably overkill, but no harm
	active_targets_(0), active_decoys_(0),
	fifo_alloc_peptides_(fifo_page_size << 20),
	fifo_alloc_prog1_(fifo_page_size << 20),
	fifo_alloc_prog2_(fifo_page_size << 20),
	all_peptides(index){
	CHECK(reader_->OK());
	compiler_prog1_ = new TheoreticalPeakCompiler(&fifo_alloc_prog1_);
	compiler_prog2_ = new TheoreticalPeakCompiler(&fifo_alloc_prog2_);
	peptide_centric_ = false;
	elution_window_ = 0;
}

ActivePeptideQueue2::~ActivePeptideQueue2() {
	deque<Peptide*>::iterator i = queue_.begin();
	// for (; i != queue_.end(); ++i)
	//   delete (*i)->PB();
	fifo_alloc_peptides_.ReleaseAll();
	fifo_alloc_prog1_.ReleaseAll();
	fifo_alloc_prog2_.ReleaseAll();

	delete compiler_prog1_;
	delete compiler_prog2_;
}

// Compute the theoretical peaks of the peptide in the "back" of the queue
// (i.e. the one most recently read from disk -- the heaviest).
void ActivePeptideQueue2::ComputeTheoreticalPeaksBack() {
	theoretical_peak_set_.Clear();
	Peptide* peptide = queue_.back();
	peptide->ComputeTheoreticalPeaks(&theoretical_peak_set_, current_pb_peptide_,
		compiler_prog1_, compiler_prog2_);
}

bool ActivePeptideQueue2::isWithinIsotope(vector<double>* min_mass, vector<double>* max_mass, double mass, int* isotope_idx) {
	for (int i = *isotope_idx; i < min_mass->size(); ++i) {
		if (mass >= (*min_mass)[i] && mass <= (*max_mass)[i]) {
			if (i > *isotope_idx) {
				*isotope_idx = i;
			}
			return true;
		}
	}
	return false;
}

int ActivePeptideQueue2::SetActiveRange(vector<double>* min_mass, vector<double>* max_mass, double min_range, double max_range, vector<bool>* candidatePeptideStatus) {
	//min_range and max_range have been introduced to fix a bug 
	//introduced by m/z selection. see #222 in sourceforge
	//this has to be true:
	// min_range <= min_mass <= max_mass <= max_range

	// queue front() is lightest; back() is heaviest

	// for now delete everything, otherwise we'd have to do some wonky stuff with the enqueing
	//we can add that back later but given how "random access" our use of the index is I don't
	//know how much that optimization would actually help
	while (!queue_.empty()) 
	{
		// would delete peptide's underlying pb::Peptide;
		queue_.pop_front();
		//    delete peptide;
	}
	
	//clean out queue
	fifo_alloc_peptides_.ReleaseAll();
	fifo_alloc_prog1_.ReleaseAll();
	fifo_alloc_prog2_.ReleaseAll();

	// Enqueue all peptides that are not yet queued but are lighter than
	// max_range. For each new enqueued peptide compute the corresponding
	// theoretical peaks. Data associated with each peptide is allocated by
	// fifo_alloc_peptides_.
	bool done = false;
	/*
	if (queue_.empty() || queue_.back()->Mass() <= max_range) {
		if (!queue_.empty()) {
			ComputeTheoreticalPeaksBack();
		}
		while (!(done = reader_->Done())) {
			// read all peptides lighter than max_range
			reader_->Read(&current_pb_peptide_);
			if (current_pb_peptide_.mass() < min_range) {
				// we would delete current_pb_peptide_;
				continue; // skip peptides that fall below min_range
			}
			Peptide* peptide = new(&fifo_alloc_peptides_)
				Peptide(current_pb_peptide_, proteins_, &fifo_alloc_peptides_);
			queue_.push_back(peptide);
			if (peptide->Mass() > max_range) {
				break;
			}
			ComputeTheoreticalPeaksBack();
		}
	}
	*/
	
	auto pep_it = all_peptides->lowerBound(min_range);
	if (pep_it == all_peptides->pep_array.cend())
		done = true;

	else
	{
		while ((*pep_it)->Mass() <= max_range)
		{
			queue_.push_back(*pep_it);
			if ((*pep_it)->Mass() > max_range) {
				break;
			}
			pep_it++;
		}
	}
	


	//I believe all the actual work is above here

	// by now, if not EOF, then the last (and only the last) enqueued
	// peptide is too heavy
	assert(!queue_.empty() || done);

	// Set up iterator for use with HasNext(),
	// GetPeptide(), and NextPeptide(). Return the number of enqueued peptides.
	if (queue_.empty()) {
		return 0;
	}

	iter_ = queue_.begin();
	while (iter_ != queue_.end() && (*iter_)->Mass() < min_mass->front()){
		++iter_;
	}

	int* isotope_idx = new int(0);
	end_ = iter_;
	int active = 0;
	active_targets_ = active_decoys_ = 0;
	while (end_ != queue_.end() && (*end_)->Mass() < max_mass->back()){
		if (isWithinIsotope(min_mass, max_mass, (*end_)->Mass(), isotope_idx)) {
			++active;
			candidatePeptideStatus->push_back(true);
			if (!(*end_)->IsDecoy()) {
				++active_targets_;
			}
			else {
				++active_decoys_;
			}
		}
		else {
			candidatePeptideStatus->push_back(false);
		}
		++end_;
	}
	delete isotope_idx;
	if (active == 0) {
		return 0;
	}

	return active;

} 

//END SET ACTIVE RANGE, HERE IS WHERE THE WORK IS TO BE DONE



// Compute the b ion only theoretical peaks of the peptide in the "back" of the queue
// (i.e. the one most recently read from disk -- the heaviest).
void ActivePeptideQueue2::ComputeBTheoreticalPeaksBack() {
	theoretical_b_peak_set_.Clear();
	Peptide* peptide = queue_.back();
	peptide->ComputeBTheoreticalPeaks(&theoretical_b_peak_set_);
	b_ion_queue_.push_back(theoretical_b_peak_set_);
}

int ActivePeptideQueue2::SetActiveRangeBIons(vector<double>* min_mass, vector<double>* max_mass, double min_range, double max_range, vector<bool>* candidatePeptideStatus) {
	exact_pval_search_ = true;
	// queue front() is lightest; back() is heaviest

	// delete anything already loaded that falls below min_range
	while (!queue_.empty() && queue_.front()->Mass() < min_range) {
		Peptide* peptide = queue_.front();
		// would delete peptide's underlying pb::Peptide;
		ReportPeptideHits(peptide);
		peptide->spectrum_matches_array.clear();
		vector<Peptide::spectrum_matches>().swap(peptide->spectrum_matches_array);
		queue_.pop_front();
		b_ion_queue_.pop_front();
		//    delete peptide;
	}
	if (queue_.empty()) {
		fifo_alloc_peptides_.ReleaseAll();
	}
	else {
		Peptide* peptide = queue_.front();
		// Free all peptides up to, but not including peptide.
		fifo_alloc_peptides_.Release(peptide);
	}

	// Enqueue all peptides that are not yet queued but are lighter than
	// max_range. For each new enqueued peptide compute the corresponding
	// theoretical peaks. Data associated with each peptide is allocated by
	// fifo_alloc_peptides_.
	bool done;
	if (queue_.empty() || queue_.back()->Mass() <= max_range) {
		while (!(done = reader_->Done())) {
			// read all peptides lighter than max_range
			reader_->Read(&current_pb_peptide_);
			if (current_pb_peptide_.mass() < min_range) {
				// we would delete current_pb_peptide_;
				continue; // skip peptides that fall below min_range
			}
			Peptide* peptide = new(&fifo_alloc_peptides_)
				Peptide(current_pb_peptide_, proteins_, &fifo_alloc_peptides_);
			queue_.push_back(peptide);
			ComputeBTheoreticalPeaksBack();
			if (peptide->Mass() > max_range) {
				break;
			}
		}
	}
	// by now, if not EOF, then the last (and only the last) enqueued
	// peptide is too heavy
	assert(!queue_.empty() || done);

	iter1_ = b_ion_queue_.begin();
	iter_ = queue_.begin();
	while (iter_ != queue_.end() && (*iter_)->Mass() < min_mass->front()){
		++iter_;
		++iter1_;
	}

	int* isotope_idx = new int(0);
	end_ = iter_;
	end1_ = iter1_;
	int active = 0;
	active_targets_ = active_decoys_ = 0;
	while (end_ != queue_.end() && (*end_)->Mass() < max_mass->back()){
		if (isWithinIsotope(min_mass, max_mass, (*end_)->Mass(), isotope_idx)) {
			++active;
			candidatePeptideStatus->push_back(true);
			if (!(*end_)->IsDecoy()) {
				++active_targets_;
			}
			else {
				++active_decoys_;
			}
		}
		else {
			candidatePeptideStatus->push_back(false);
		}
		++end_;
		++end1_;
	}
	delete isotope_idx;
	if (active == 0) {
		return 0;
	}

	return active;
}

int ActivePeptideQueue2::CountAAFrequency(
	double binWidth,
	double binOffset,
	double** dAAFreqN,
	double** dAAFreqI,
	double** dAAFreqC,
	int** dAAMass
	) {

	unsigned int i = 0;
	unsigned int cntTerm = 0;
	unsigned int cntInside = 0;
	const unsigned int MaxModifiedAAMassBin = 2000 / binWidth;   //2000 is the maximum size of a modified amino acid 
	unsigned int* nvAAMassCounterN = new unsigned int[MaxModifiedAAMassBin];   //N-terminal amino acids
	unsigned int* nvAAMassCounterC = new unsigned int[MaxModifiedAAMassBin];   //C-terminal amino acids
	unsigned int* nvAAMassCounterI = new unsigned int[MaxModifiedAAMassBin];   //inner amino acids in the peptides
	memset(nvAAMassCounterN, 0, MaxModifiedAAMassBin * sizeof(unsigned int));
	memset(nvAAMassCounterC, 0, MaxModifiedAAMassBin * sizeof(unsigned int));
	memset(nvAAMassCounterI, 0, MaxModifiedAAMassBin * sizeof(unsigned int));

	while (!(reader_->Done())) { // read all peptides in index
		reader_->Read(&current_pb_peptide_);
		Peptide* peptide = new(&fifo_alloc_peptides_) Peptide(current_pb_peptide_, proteins_, &fifo_alloc_peptides_);

		double* dAAResidueMass = peptide->getAAMasses(); //retrieves the amino acid masses, modifications included

		int nLen = peptide->Len(); //peptide length
		++nvAAMassCounterN[(unsigned int)(dAAResidueMass[0] / binWidth + 1.0 - binOffset)];
		for (i = 1; i < nLen - 1; ++i) {
			++nvAAMassCounterI[(unsigned int)(dAAResidueMass[i] / binWidth + 1.0 - binOffset)];
			++cntInside;
		}
		++nvAAMassCounterC[(unsigned int)(dAAResidueMass[nLen - 1] / binWidth + 1.0 - binOffset)];
		++cntTerm;

		delete[] dAAResidueMass;
		fifo_alloc_peptides_.ReleaseAll();
	}

	//calculate the unique masses
	unsigned int uiUniqueMasses = 0;
	for (i = 0; i < MaxModifiedAAMassBin; ++i) {
		if (nvAAMassCounterN[i] || nvAAMassCounterI[i] || nvAAMassCounterC[i]) {
			++uiUniqueMasses;
		}
	}

	//calculate the unique amino acid masses
	*dAAMass = new int[uiUniqueMasses];     //a vector for the unique (integerized) amino acid masses present in the sample
	*dAAFreqN = new double[uiUniqueMasses]; //a vector for the amino acid frequencies at the N-terminus
	*dAAFreqI = new double[uiUniqueMasses]; //a vector for the amino acid frequencies inside the peptide
	*dAAFreqC = new double[uiUniqueMasses]; //a vector for the amino acid frequencies at the C-terminus
	unsigned int cnt = 0;
	for (i = 0; i < MaxModifiedAAMassBin; ++i) {
		if (nvAAMassCounterN[i] || nvAAMassCounterI[i] || nvAAMassCounterC[i]) {
			(*dAAFreqN)[cnt] = (double)nvAAMassCounterN[i] / cntTerm;
			(*dAAFreqI)[cnt] = (double)nvAAMassCounterI[i] / cntInside;
			(*dAAFreqC)[cnt] = (double)nvAAMassCounterC[i] / cntTerm;
			(*dAAMass)[cnt] = i;
			cnt++;
		}
	}

	delete[] nvAAMassCounterN;
	delete[] nvAAMassCounterI;
	delete[] nvAAMassCounterC;
	return uiUniqueMasses;
}

void ActivePeptideQueue2::ReportPeptideHits(Peptide* peptide) {
	if (!peptide_centric_) {
		return;
	}

	current_peptide_ = peptide;
	TideMatchSet matches(peptide, highest_mz_);
	matches.exact_pval_search_ = exact_pval_search_;
	matches.elution_window_ = elution_window_;

	/*if (!output_files_) { //only tab-delimited output is supported
		matches.report(target_file_, decoy_file_, top_matches_,
			this, proteins_, *locations_, compute_sp_);
	}*/
}

