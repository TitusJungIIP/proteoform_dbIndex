package edu.scripps.yates.utilities.annotations.uniprot;

public enum ProteinExistence {
	PROTEIN_LEVEL(1, "Evidence at protein level"), //
	TRANSCRIPT_LEVEL(2, "Evidence at transcript level"), //
	INFERRED_FROM_HOMOLOGY(3, "Inferred from homology"), //
	PROTEIN_PREDCTED(4, "Predicted"), //
	PROTEIN_UNCERTAIN(5, "Uncertain");
	private final int num;
	private final String description;

	private ProteinExistence(int num, String description) {
		this.num = num;
		this.description = description;
	}

	public int getNum() {
		return num;
	}

	public String getDescription() {
		return description;
	}
}
