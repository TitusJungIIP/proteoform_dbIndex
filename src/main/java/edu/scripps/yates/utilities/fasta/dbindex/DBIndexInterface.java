package edu.scripps.yates.utilities.fasta.dbindex;

import java.util.List;
import java.util.Set;

public interface DBIndexInterface {

	public List<IndexedSequence> getSequences(double precursorMass, double massTolerance) throws DBIndexStoreException;

	public List<edu.scripps.yates.utilities.fasta.dbindex.IndexedSequence> getSequences(List<MassRange> massRanges)
			throws DBIndexStoreException;

	public List<IndexedProtein> getProteins(IndexedSequence seq)
			throws edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;

	public Set<IndexedProtein> getProteins(String seq) throws DBIndexStoreException;

	public IndexedProtein getIndexedProteinById(int proteinId) throws DBIndexStoreException;

	public String getProteinSequenceById(int proteinId) throws DBIndexStoreException;
}
