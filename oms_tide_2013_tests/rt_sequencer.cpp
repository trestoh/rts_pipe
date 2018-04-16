#include "rt_sequencer.h"

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