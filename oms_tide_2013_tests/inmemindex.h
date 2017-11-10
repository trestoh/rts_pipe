

#include <gflags/gflags.h>
#include <app/tide/records.h>
#include <peptides.pb.h>
#include <app/tide/peptide.h>
#include <app/tide/theoretical_peak_set.h>
#include <app/tide/compiler.h>
#include <app/TideMatchSet.h>
#include <vector>
#include <algorithm>
#include <app/tide/active_peptide_queue.h>
#include <io/carp.h>

#ifndef IN_MEM_INDEX_H
#define IN_MEM_INDEX_H

DECLARE_int32(fifo_page_size);

class InMemIndex {

public:
	InMemIndex(RecordReader* reader, const vector<const pb::Protein*>& proteins) : fifo_alloc_peptides(FLAGS_fifo_page_size << 20), fifo_alloc_prog1(FLAGS_fifo_page_size << 20),
		fifo_alloc_prog2(FLAGS_fifo_page_size << 20), theoretical_peak_set(2000)
	{

		compiler_prog1 = new TheoreticalPeakCompiler(&fifo_alloc_prog1);
		compiler_prog2 = new TheoreticalPeakCompiler(&fifo_alloc_prog2);

		bool done = false;
		pb::Peptide current_pb_peptide;

		while (!(done = reader->Done())) {
			// read all peptides
			theoretical_peak_set.Clear();
			reader->Read(&current_pb_peptide);
			Peptide* peptide = new(&fifo_alloc_peptides) Peptide(current_pb_peptide, proteins, &fifo_alloc_peptides);

			peptide->ComputeTheoreticalPeaks(&theoretical_peak_set, current_pb_peptide,
				compiler_prog1, compiler_prog2);

			pep_array.push_back(peptide);

		}
	}

	~InMemIndex()
	{
		vector<Peptide*>::iterator i = pep_array.begin();
		// for (; i != queue_.end(); ++i)
		//   delete (*i)->PB();
		fifo_alloc_peptides.ReleaseAll();
		fifo_alloc_prog1.ReleaseAll();
		fifo_alloc_prog2.ReleaseAll();
	}

	vector<Peptide*> pep_array;

	vector<Peptide*>::iterator lowerBound(const double& mass)
	{

		vector<Peptide*>::iterator first = pep_array.begin();
		vector<Peptide*>::iterator last = pep_array.end();

		vector<Peptide*>::iterator it;
		std::iterator_traits<vector<Peptide*>::iterator>::difference_type count, step;
		count = std::distance(first, last);

		while (count > 0) {
			it = first;
			step = count / 2;
			std::advance(it, step);
			if ((*it)->Mass() < mass) {
				first = ++it;
				count -= step + 1;
			}
			else
				count = step;
		}
		return first;
	}


private:
	
	FifoAllocator fifo_alloc_peptides;
	ST_TheoreticalPeakSet theoretical_peak_set;
	FifoAllocator fifo_alloc_prog1;
	FifoAllocator fifo_alloc_prog2;
	TheoreticalPeakCompiler* compiler_prog1;
	TheoreticalPeakCompiler* compiler_prog2;

};


#endif //IN_MEM_INDEX_H