package edu.scripps.yates.utilities.grouping;

/**
 *
 * @author gorka
 */
public enum ProteinEvidence {
	/**
	 * At least one unique peptides
	 */
	CONCLUSIVE,
	/**
	 * Shared Discriminating peptides
	 */
	AMBIGUOUSGROUP,
	/**
	 * Same peptides and at least one Discriminating
	 */
	INDISTINGUISHABLE,
	/**
	 * Only NonDiscrimitating peptides
	 */
	NONCONCLUSIVE,
	/**
	 * No peptides
	 */
	FILTERED
}
