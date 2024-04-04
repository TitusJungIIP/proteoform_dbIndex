package edu.scripps.yates.utilities.fasta.dbindex;

public abstract class PeptideFilter {
	public abstract boolean isValid(String peptideSequence);

	@Override
	public abstract String toString();
}
