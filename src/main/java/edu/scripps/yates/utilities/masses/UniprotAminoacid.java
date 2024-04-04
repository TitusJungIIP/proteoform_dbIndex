package edu.scripps.yates.utilities.masses;

public class UniprotAminoacid {
	private final char oneLetterCode;
	private final String threeLetterCode;
	private final String aminoacidName;
	private Double monoMass;

	public UniprotAminoacid(char oneLetterCode, String threeLetterCode, String aminoacidName) {
		super();
		this.oneLetterCode = oneLetterCode;
		this.threeLetterCode = threeLetterCode;
		this.aminoacidName = aminoacidName;
	}

	public char getOneLetterCode() {
		return oneLetterCode;
	}

	public String getThreeLetterCode() {
		return threeLetterCode;
	}

	public String getAminoacidName() {
		return aminoacidName;
	}

	public Double getMonoMass() {
		if (monoMass == null) {
			AssignMass.getInstance(true);
			monoMass = AssignMass.getMass(oneLetterCode);
		}
		return monoMass;
	}
}
