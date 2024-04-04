package edu.scripps.yates.utilities.proteomicsmodel.factories;

import java.io.Serializable;

import edu.scripps.yates.utilities.proteomicsmodel.AbstractProtein;
import edu.scripps.yates.utilities.proteomicsmodel.Accession;
import edu.scripps.yates.utilities.proteomicsmodel.Organism;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AccessionType;

public class ProteinEx extends AbstractProtein implements Serializable {
	/**
	 *
	 */
	private static final long serialVersionUID = -1435542806814270031L;

	public ProteinEx(AccessionType accessionType, String accession, Organism organism) {
		final AccessionEx accessionEx = new AccessionEx(accession, accessionType);
		setPrimaryAccession(accessionEx);
		setOrganism(organism);
		setKey(getAccession());
	}

	public ProteinEx(String accession) {
		setPrimaryAccession(accession);
		setKey(getAccession());
	}

	public ProteinEx(AccessionType accessionType, String accession) {
		this(accessionType, accession, null);
	}

	public ProteinEx(Accession accession) {
		setPrimaryAccession(accession);
		setKey(getAccession());
	}

}
