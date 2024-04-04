package edu.scripps.yates.utilities.proteomicsmodel;

import java.util.Set;

public interface HasPeptides {
	public boolean addPeptide(Peptide peptide, boolean recursively);

	public Set<Peptide> getPeptides();
}
