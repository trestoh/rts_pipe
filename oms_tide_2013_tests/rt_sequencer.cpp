#include "stdafx.h"

#include "rt_sequencer.h"

RT_Sequencer::RT_Sequencer()
{
	pep_index = NULL;
	peptide_queue = NULL;
	negative_isotope_errors = NULL;
	curr_spectrum = NULL;
}

bool RT_Sequencer::protInclude(std::string prot_name)
{
	if (inclusion_list.empty())
		return true;

	if (inclusion_list.count(prot_name))
		return true;

	else
		return false;

}

void RT_Sequencer::addPeak(double mz, double intensity)
{
	if (curr_spectrum == NULL)
		return;

	curr_spectrum->AddPeak(mz, intensity);

}

void RT_Sequencer::makeSpec(int scan_num, double prec_mz)
{
	if (curr_spectrum != NULL)
		delete curr_spectrum;

	curr_spectrum = new Spectrum(scan_num, prec_mz);
}

double RT_Sequencer::score()
{
	return last_score;
}

void RT_Sequencer::init(const char* index_name, const char* param_name, const char* inclusion_name, double xc1, double xc2, double xc3, double dcn, int window_type_index, double window_size)
{
	//wchar_t* msgbuf = (wchar_t*) index_name;
	set_verbosity_level(30);
	//sprintf(msgbuf, "Index is %s", index_name);

	OutputDebugStringA(index_name);

	xcorr1 = xc1;
	xcorr2 = xc2;
	xcorr3 = xc3;
	deltacn = dcn;

	low_peak = 0;

	last_score = 0.0;
	last_seq = "";

	//open_ms_seq = OpenMS::AASequence::fromString("EHQEALHQQR");
	//computeFragments(3);
	//for (int j = 0; j < last_fragments.size(); j++)
	//carp(CARP_INFO, "Theoretical Fragment of %f", last_fragments.at(j));
	//	std::cout << "Theoretical Fragment of " << last_fragments.at(j) << std::endl;

	HAS_DECOYS = false;
	PROTEIN_LEVEL_DECOYS = false;

	if (index_name == NULL)
	{
		OutputDebugString(L"Error Reading Index, name misread");
		carp(CARP_FATAL, "Error Reading Index, no name provided");
		exit(2);
	}

	if (param_name != NULL)
	{
		parse_parameter_file(param_name);
		precursor_window = Params::GetDouble("precursor-window");
		window_type = string_to_window_type(Params::GetString("precursor-window-type"));
	}

	else
	{
		precursor_window = window_size;
		window_type = WINDOW_TYPE_T(window_type_index);
	}

	if (inclusion_name != NULL)
	{
		std::ifstream inc_list(inclusion_name);
		if (inc_list)
		{
			std::string line;
			std::string prot;
			std::string kin_family;
			while (std::getline(inc_list, line))
			{
				std::istringstream iss(line);
				iss >> prot >> kin_family;
				inclusion_list.insert({ prot, kin_family });
			}
		}
	}

	std::string peptides_file = FileUtils::Join(index_name, "pepix");
	std::string proteins_file = FileUtils::Join(index_name, "protix");
	std::string auxlocs_file = FileUtils::Join(index_name, "auxlocs");

	//OutputDebugStringA(peptides_file.c_str());

	pb::Header pepts_header;
	carp(CARP_INFO, "Declaring Headed Record Reader");
	HeadedRecordReader pept_reader(peptides_file, &pepts_header);
	carp(CARP_INFO, "Declared Headed Record Reader");

	if ((pepts_header.file_type() != pb::Header::PEPTIDES) || !pepts_header.has_peptides_header())
	{
		OutputDebugString(L"Error Reading Index, Peptide Header not found");
		Sleep(5000);
		carp(CARP_FATAL, "Error Reading Index, Peptide Header not found");
		//exit(2);
	}

	const pb::Header::PeptidesHeader& pepsHeader = pepts_header.peptides_header();
	Params::Set("enzyme", pepsHeader.enzyme());
	const char* digestString = digest_type_to_string(pepsHeader.full_digestion() ? FULL_DIGEST : PARTIAL_DIGEST);
	Params::Set("digestion", digestString);
	Params::Set("isotopic-mass", pepsHeader.monoisotopic_precursor() ? "mono" : "average");

	Params::Finalize();
	GlobalParams::set();

	const pb::ModTable mods = pepsHeader.mods();

	static_mods = false;

	auto s_mods = mods.static_mod();

	auto s_walker = s_mods.begin();
	s_walker++;

	if (s_walker == s_mods.end())
		static_delta = "";

	while (s_walker != s_mods.end())
	{
		pb::Modification ok = (*s_walker);
		//double is_tmt = ok.delta() - 229.162932;
		//if (is_tmt <= 0.000001)
		//	static_mods = true;

		static_mods = true;
		static_delta = to_string(ok.delta());

		s_walker++;
	}

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
		OutputDebugString(L"Error in isotope_error parameter formatting");
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

	pb::Header protein_header;

	if (!ReadRecordsToVector<pb::Protein, const pb::Protein>(&proteins, proteins_file, &protein_header))
		carp(CARP_FATAL, "Error Reading Index!");


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

	stringstream ss;
	ss << Params::GetString("enzyme") << '-' << Params::GetString("digestion");
	TideMatchSet::CleavageType = ss.str();

	spectrum_min_mz = Params::GetDouble("spectrum-min-mz");
	spectrum_max_mz = Params::GetDouble("spectrum-max-mz");
	min_peaks = Params::GetInt("min-peaks");
	top_matches = Params::GetInt("top-match");
	highest_mz = 0.0f;

	peptide_queue = new ActivePeptideQueue2(pep_index, peptide_reader->Reader(), proteins);
	peptide_queue->SetBinSize(bin_width, bin_offset);

	OutputDebugString(L"Sequencer Initialized");

}

RT_Sequencer::~RT_Sequencer()
{
	delete pep_index;
	delete peptide_queue;
	delete negative_isotope_errors;
	delete curr_spectrum;
}

bool RT_Sequencer::is_match(double high_mz, double precursor_mass, int precursor_charge)
{
	highest_mz = high_mz;
	vector<SpectrumCollection::SpecCharge> spec;
	vector<SpectrumCollection::SpecCharge>* spec_charges = &spec;

	spec_charges->push_back(SpectrumCollection::SpecCharge(precursor_mass, precursor_charge, curr_spectrum, 0));
	int elution_window = Params::GetInt("elution-window-size");
	bool use_neutral_loss_peaks = Params::GetBool("use-neutral-loss-peaks");
	bool use_flanking_peaks = Params::GetBool("use-flanking-peaks");
	int max_charge = Params::GetInt("max-precursor-charge");
	ObservedPeakSet observed(bin_width, bin_offset, use_neutral_loss_peaks, use_flanking_peaks);

	peptide_queue->setElutionWindow(elution_window);
	peptide_queue->setPeptideCentric(false);

	if (elution_window > 0 && elution_window % 2 == 0)
		peptide_queue->setElutionWindow(elution_window + 1);

	peptide_queue->setElutionWindow(0);

	for (vector<SpectrumCollection::SpecCharge>::const_iterator sc = spec_charges->begin(); sc < spec_charges->end(); sc++)
	{
		Spectrum* spectrum = sc->spectrum;
		double precursor_mz = spectrum->PrecursorMZ();
		int charge = sc->charge;
		int scan_num = spectrum->SpectrumNumber();

		OutputDebugString(L"In Spec Collec loop, peaks for spectrum is: ");
		char spec_size[100];
		sprintf(spec_size, "%d", spectrum->Size());
		OutputDebugStringA(spec_size);


		if (precursor_mz < spectrum_min_mz || precursor_mz > spectrum_max_mz || scan_num < min_scan || scan_num > max_scan
			|| spectrum->Size() < min_peaks || (charge_to_search != 0 && charge != charge_to_search) || charge > max_charge) {
			carp(CARP_INFO, "Exited search early");
			if (spectrum->Size() < min_peaks)
			{
				carp(CARP_INFO, "Not Enough Peaks.");
				low_peak++;
			}
			delete curr_spectrum;
			curr_spectrum = NULL;
			return false;
		}

		vector<double>* min_mass = new vector<double>();
		vector<double>* max_mass = new vector<double>();
		vector<bool>* candidatePeptideStatus = new vector<bool>();
		double min_range, max_range;
		double* min_ran = &min_range;
		double* max_ran = &max_range;

		computeWindow(*sc, window_type, precursor_window, max_charge, negative_isotope_errors, min_mass, max_mass, &min_range, &max_range);

		
		observed.PreprocessSpectrum(*spectrum, charge);
		int nCandPeptide = peptide_queue->SetActiveRange(min_mass, max_mass, min_range, max_range, candidatePeptideStatus);
		if (nCandPeptide == 0) {
			delete min_mass;
			delete max_mass;
			delete candidatePeptideStatus;
			delete curr_spectrum;
			curr_spectrum = NULL;
			return false;
		}

		int candidatePeptideStatusSize = candidatePeptideStatus->size();
		TideMatchSet::Arr2 match_arr2(candidatePeptideStatusSize);

		OutputDebugString(L"Hitting collect scores compiled");
		this->collectScoresCompiled(peptide_queue, spectrum, observed, &match_arr2, candidatePeptideStatusSize, charge);
		OutputDebugString(L"Passed collect scores compiled");

			
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
		int max_corr_rank = match_arr.size();
		int precision = Params::GetInt("precision");

		for (TideMatchSet::Arr::iterator it = match_arr.begin(); it != match_arr.end(); ++it) {
			if (it->xcorr_score > max_corr) {
				second_corr = max_corr;
				max_corr = it->xcorr_score;
				max_corr_rank = it->rank;
			}
		}

		last_score = max_corr;

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
		/*
		const Peptide& peptid = *(peptide_queue->GetPeptide(max_corr_rank));

		//if (peptid.IsDecoy())
		//	decoy = true;
		//else
		//	decoy = false;

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
		*/
		//peptide_hit = cruxPep.getModifiedSequenceWithMasses();

		//carp(CARP_DETAILED_DEBUG, "Max XCorr for Scan is: %f and dCn is: %f", max_corr, delta_corr);
		const Peptide& peptid01 = *(peptide_queue->GetPeptide(max_corr_rank));
		const pb::Protein* protein = proteins[peptid01.FirstLocProteinId()];
		int pos = peptid01.FirstLocPos();
		string proteinName = protein->name();//getProteinName(*protein,
		//(!protein->has_target_pos()) ? pos : protein->target_pos());
		int start_pos = 0;
		while (proteinName[start_pos] != '|')
			start_pos++;
		start_pos++;
		int end_pos = start_pos;
		while (proteinName[end_pos] != '-' && proteinName[end_pos] != '|')
			end_pos++;
		proteinName = proteinName.substr(start_pos, end_pos - start_pos);


		if (max_corr >= threshold && delta_corr >= delta_threshold && protInclude(proteinName)) {
			const Peptide& peptid = *(peptide_queue->GetPeptide(max_corr_rank));
			last_seq = peptid.SeqWithMods();
			//if (last_seq.find('C') != std::string::npos)
				//std::cout << last_seq << std::endl;
			for (int i = 0; i < last_seq.size(); i++)
			{
				if (last_seq.at(i) == '[')
					if (last_seq.substr(i, 8) == "[+229.2]")
					{
						//carp(CARP_INFO, "Found a dang TMT");
						//carp(CARP_INFO, "Sequence started as %s", last_seq.c_str());
						//last_seq.replace(i, 8, "(TMT6plex)");
						if (i == 1)
						{
							last_seq.erase(i, 8);
							last_seq = "(TMT6plex)" + last_seq;
							i = 10;
						}
						else
						{
							last_seq.replace(i, 8, "(TMT6plex)");
							i += 9;
						}
						//carp(CARP_INFO, "Changed sequence to %s", last_seq.c_str());
					}
			}

			if (static_mods)
			{
				//last_seq.insert(0, ".(TMT6plex)");
				last_seq.insert(0, (".[+" + static_delta + "]"));
				if (last_seq.at(last_seq.size() - 1) == 'K')
					//last_seq.append("(TMT6plex)");
					last_seq.append(("[+" + static_delta + "]"));
			}

			open_ms_seq = OpenMS::AASequence::fromString(last_seq);

			computeFragments(precursor_charge);
			//for (int j = 0; j < last_fragments.size(); j++)
			//carp(CARP_INFO, "Theoretical Fragment of %f", last_fragments.at(j));
			//	std::cout << "Theoretical Fragment of " << last_fragments.at(j) << std::endl;

			//for (int k = 0; k < curr_spectrum->Size(); k++)
			//	std::cout << "Observed Peak MZ: " << curr_spectrum->M_Z(k) << " Intensity: " << curr_spectrum->Intensity(k) << std::endl;

			//for (int m = 0; m < sps_targets.size(); m++)
			//	std::cout << "Sps Target: " << sps_targets.at(m) << std::endl;

			//carp(CARP_INFO, "Search was at charge state %d", precursor_charge);
			delete min_mass;
			delete max_mass;
			delete candidatePeptideStatus;
			delete curr_spectrum;
			curr_spectrum = NULL;
			return true;
		}


		

		delete min_mass;
		delete max_mass;
		delete candidatePeptideStatus;

	}//end loop through spectrum charges

	delete curr_spectrum;
	curr_spectrum = NULL;
	return false;

}

void RT_Sequencer::computeFragments(int charge)
{
	last_fragments.clear();
	sps_targets.clear();

	std::vector<std::vector<double>> temp_fragment_lists;
	temp_fragment_lists.push_back(std::vector<double>());
	temp_fragment_lists.push_back(std::vector<double>());
	temp_fragment_lists.push_back(std::vector<double>());
	temp_fragment_lists.push_back(std::vector<double>());

	bool nterm = open_ms_seq.hasNTerminalModification() ? ((open_ms_seq.getNTerminalModificationName() == "TMT6plex") ? true : false) : false;

	nterm = open_ms_seq.hasNTerminalModification();

	int second_tmt = open_ms_seq.size();
	bool second_mod = false;

	OpenMS::AASequence::Iterator it = --open_ms_seq.end();
	while (it != open_ms_seq.begin())
	{
		if (it->isModified())
		{
			std::string mod = it->getModificationName();
			carp(CARP_INFO, "%s", mod);
			if (mod == "TMT6plex")
			{
				second_mod = true;
				break;
			}

			//second_mod = true;
			//break;

		}
		--second_tmt;
		--it;
	}

	int right_index = open_ms_seq.size() - second_tmt;

	//if (second_mod && right_index != 0)
	//	bool temp = true;

	for (int i = 1; i < open_ms_seq.size(); i++)
	{
		if (nterm)
		{
			temp_fragment_lists.at(0).push_back(open_ms_seq.getPrefix(i).getMonoWeight(OpenMS::Residue::ResidueType::BIon, 1) - 18.0091422);
			temp_fragment_lists.at(0).push_back(open_ms_seq.getPrefix(i).getMonoWeight(OpenMS::Residue::ResidueType::BIon, 1) - 17.0086343);
			temp_fragment_lists.at(0).push_back(open_ms_seq.getPrefix(i).getMonoWeight(OpenMS::Residue::ResidueType::BIon, 1));
		}

		if (second_mod)
		{
			if (i > right_index)
			{
				temp_fragment_lists.at(1).push_back(open_ms_seq.getSuffix(i).getMonoWeight(OpenMS::Residue::ResidueType::YIon, 1) - 18.0091422);
				temp_fragment_lists.at(1).push_back(open_ms_seq.getSuffix(i).getMonoWeight(OpenMS::Residue::ResidueType::YIon, 1) - 17.0086343);
				temp_fragment_lists.at(1).push_back(open_ms_seq.getSuffix(i).getMonoWeight(OpenMS::Residue::ResidueType::YIon, 1));
			}
		}

		if (charge > 2)
		{
			if (nterm)
				temp_fragment_lists.at(2).push_back(open_ms_seq.getPrefix(i).getMonoWeight(OpenMS::Residue::ResidueType::BIon, 2) / 2);

			if (second_mod)
				if (i > right_index)
					temp_fragment_lists.at(3).push_back(open_ms_seq.getSuffix(i).getMonoWeight(OpenMS::Residue::ResidueType::YIon, 2) / 2);
		}
	}

	std::vector<double>::iterator iter[4];
	iter[0] = temp_fragment_lists.at(0).begin();
	iter[1] = temp_fragment_lists.at(1).begin();
	iter[2] = temp_fragment_lists.at(2).begin();
	iter[3] = temp_fragment_lists.at(3).begin();

	std::priority_queue<fragmentNode, std::vector<fragmentNode>, std::greater<fragmentNode>> queue;
	for (int i = 0; i < 4; i++)
	{
		if (iter[i] != temp_fragment_lists.at(i).end())
		{
			queue.push(fragmentNode(*(iter[i]++), i));
		}
	}

	fragmentNode temp(0.0, 0);

	while (!queue.empty())
	{
		temp = queue.top();
		//std::cout << "Queue Test: "  << temp.mz << std::endl;
		last_fragments.push_back(temp.mz);
		queue.pop();
		if (iter[temp.list] != temp_fragment_lists.at(temp.list).end())
		{
			queue.push(fragmentNode(*(iter[temp.list]++), temp.list));
		}
	}

	std::priority_queue<peakNode, std::vector<peakNode>, std::greater<peakNode>> peak_queue;
	int peaks = 0;
	int frags = 0;
	while (peaks < curr_spectrum->Size() && frags < last_fragments.size())
	{
		double mz_diff = curr_spectrum->M_Z(peaks) - last_fragments.at(frags);

		if (abs(mz_diff) < 0.5)
		{
			if (peak_queue.size() != 10)
				peak_queue.push(peakNode(curr_spectrum->M_Z(peaks), curr_spectrum->Intensity(peaks)));
			else if (curr_spectrum->Intensity(peaks) > peak_queue.top().intensity)
			{
				peak_queue.pop();
				peak_queue.push(peakNode(curr_spectrum->M_Z(peaks), curr_spectrum->Intensity(peaks)));
			}

			peaks++;

		}

		else if (mz_diff >= 0.0)
			frags++;
		else
			peaks++;
	}

	while (!peak_queue.empty())
	{
		sps_targets.push_back(peak_queue.top().mz);
		peak_queue.pop();
	}


}

void RT_Sequencer::collectScoresCompiled(
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
		//try{
		((void(*)(void))prog)();
		//}

		//catch (...)
		//{
		//carp(CARP_INFO, "Bad Bad Not Good...");
		//carp(CARP_INFO, "Queue Dump");
		//active_peptide_queue->dumpQueue(charge);
		//}

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

void RT_Sequencer::computeWindow(
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