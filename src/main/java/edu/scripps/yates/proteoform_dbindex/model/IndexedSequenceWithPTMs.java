package edu.scripps.yates.proteoform_dbindex.model;

import java.util.HashSet;
import java.util.Set;

import edu.scripps.yates.utilities.fasta.dbindex.IndexedSequence;

public class IndexedSequenceWithPTMs extends IndexedSequence {
	private final Set<PTM> ptms = new HashSet<PTM>();

	public IndexedSequenceWithPTMs(long id, double mass, String sequence, String resLeft, String resRight) {
		super(id, mass, sequence, resLeft, resRight);
	}

	public void addPTM(PTM ptm) {
		ptms.add(ptm);
	}

	public int getBySize() {
		final int numPtms = ptms != null ? ptms.size() : 0;
		final int ret = 8 + 2 + 2 + (3 * numPtms) + 1 + (4 * getProteinIds().size()) + 4;
		return ret;
	}

}
