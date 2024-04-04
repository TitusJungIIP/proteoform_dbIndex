package edu.scripps.yates.utilities.proteomicsmodel;

public class PrincipalInvestigator {
	private final String name;
	private String email;
	private String institution;
	private String country;

	public PrincipalInvestigator(String name) {
		this.name = name;
	}

	public String getEmail() {
		return email;
	}

	public void setEmail(String email) {
		this.email = email;
	}

	public String getInstitution() {
		return institution;
	}

	public void setInstitution(String institution) {
		this.institution = institution;
	}

	public String getCountry() {
		return country;
	}

	public void setCountry(String country) {
		this.country = country;
	}

	public String getName() {
		return name;
	}

}
