package edu.scripps.yates.utilities.proteomicsmodel;

import java.util.List;

import edu.scripps.yates.utilities.proteomicsmodel.enums.AccessionType;

public interface Accession {

	public AccessionType getAccessionType();

	public String getAccession();

	public String getDescription();

	public boolean equals(Accession accession);

	public List<String> getAlternativeNames();

	public void setDescription(String description);

	public void setAccession(String accession);

	public void setAccessionType(AccessionType accessionType);

}
