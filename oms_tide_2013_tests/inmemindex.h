

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

class InMemIndex {

public:
	InMemIndex(RecordReader* reader, const vector<const pb::Protein*>& proteins)
	{
		bool done = false;
		pb::Peptide current_pb_peptide;

		while (!(done = reader->Done())) {
			// read all peptides
			reader->Read(&current_pb_peptide);
			//Peptide* peptide = new(&fifo_alloc_peptides) Peptide(current_pb_peptide, proteins, &fifo_alloc_peptides);

			pep_array.push_back(current_pb_peptide);

		}
	}

	vector<Peptide*> pep_array;

	vector<pb::Peptide>::iterator lowerBound(const double& mass)
	{

		vector<pb::Peptide>::iterator first = pep_array.begin();
		vector<pb::Peptide>::iterator last = pep_array.end();

		vector<pb::Peptide>::iterator it;
		std::iterator_traits<vector<pb::Peptide>::iterator>::difference_type count, step;
		count = std::distance(first, last);

		while (count > 0) {
			it = first;
			step = count / 2;
			std::advance(it, step);
			if (it->mass() < mass) {
				first = ++it;
				count -= step + 1;
			}
			else
				count = step;
		}
		return first;
	}


private:
	
	//FifoAllocator fifo_alloc_peptides;

};


#endif //IN_MEM_INDEX_H