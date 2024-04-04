package edu.scripps.yates.utilities.proteomicsmodel;

public interface HasPeptide {
	public boolean setPeptide(Peptide peptide, boolean recursively);

	public Peptide getPeptide();
}
