package edu.scripps.yates.utilities.proteomicsmodel.factories;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang.builder.HashCodeBuilder;

import edu.scripps.yates.utilities.proteomicsmodel.Accession;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AccessionType;
import edu.scripps.yates.utilities.util.Pair;

public class AccessionEx implements Accession, Serializable {

	private static final long serialVersionUID = 285204540497978045L;
	private AccessionType accessionType;
	private String accession;
	private String description;
	private List<String> alternativeNames;

	public AccessionEx() {

	}

	public AccessionEx(String accession, AccessionType accessionType) {
		this.accession = accession;
		this.accessionType = accessionType;
	}

	public AccessionEx(Pair<String, AccessionType> acc) {
		accession = acc.getFirstelement();
		accessionType = acc.getSecondElement();
	}

	/**
	 * @param description
	 *            the description to set
	 */
	@Override
	public void setDescription(String description) {
		this.description = description;
	}

	@Override
	public AccessionType getAccessionType() {
		return accessionType;
	}

	@Override
	public String getAccession() {
		return accession;
	}

	@Override
	public String getDescription() {
		return description;
	}

	@Override
	public boolean equals(Accession accession) {
		if (accession.getAccessionType().equals(getAccessionType()))
			if (accession.getAccession().equalsIgnoreCase(getAccession()))
				return true;
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		if (obj instanceof Accession)
			return this.equals((Accession) obj);
		return super.equals(obj);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return getAccession();
	}

	/**
	 * @return the alternativeNames
	 */
	@Override
	public List<String> getAlternativeNames() {
		return alternativeNames;
	}

	/**
	 * @param alternativeNames
	 *            the alternativeNames to set
	 */
	public void setAlternativeNames(List<String> alternativeNames) {
		this.alternativeNames = alternativeNames;
	}

	public void addAlternativeName(String alternativeName) {
		if (alternativeNames == null)
			alternativeNames = new ArrayList<String>();
		alternativeNames.add(alternativeName);
	}

	@Override
	public int hashCode() {
		return new HashCodeBuilder().append(getAccessionType()).append(getAccession()).toHashCode();
	}

	@Override
	public void setAccession(String accession) {
		this.accession = accession;
	}

	@Override
	public void setAccessionType(AccessionType accessionType) {
		this.accessionType = accessionType;
	}
}
