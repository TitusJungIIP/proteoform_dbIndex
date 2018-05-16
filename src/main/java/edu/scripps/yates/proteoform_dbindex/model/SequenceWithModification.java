package edu.scripps.yates.proteoform_dbindex.model;

/**
 * Represents a peptide or protein sequence with a modification, whether it is a
 * PTM or a sequence variation
 * 
 * @author salvador
 *
 */
public class SequenceWithModification {
	private final String sequenceWithModification;
	private String sequenceAfterModification;
	private final String sequenceBeforeModification;
	private Double sequenceMassAfterModification;
	private final String proteinSequence;
	private final ExtendedAssignMass extendedAssignMass;

	public SequenceWithModification(String originalsequenceSequence, String sequenceWithModification,
			ExtendedAssignMass extendedAssignMass, String proteinSequence) {
		sequenceBeforeModification = originalsequenceSequence;
		this.sequenceWithModification = sequenceWithModification;
		this.extendedAssignMass = extendedAssignMass;
		this.proteinSequence = proteinSequence;

	}

	public String getSequenceAfterModification() {
		if (sequenceAfterModification == null) {
			final StringBuilder sb = new StringBuilder();
			StringBuilder change = null;
			for (int i = 0; i < sequenceWithModification.length(); i++) {
				final char aa = sequenceWithModification.charAt(i);
				if (aa == '[') {
					change = new StringBuilder();
					continue;
				}
				if (aa == ']') {
					if (change != null) {
						if (change.toString().contains("->")) {
							final String[] split = change.toString().split("->");
							if (split.length > 1) {
								sb.append(split[1]);
							}
						}
						change = null;
					}
				} else {
					if (change != null) {
						change.append(aa);
					} else {
						sb.append(aa);
					}
				}
			}
			sequenceAfterModification = sb.toString();
		}
		return sequenceAfterModification;
	}

	public String getSequenceWithModification() {
		return sequenceWithModification;
	}

	public String getSequenceBeforeModification() {
		return sequenceBeforeModification;
	}

	@Override
	public String toString() {
		final StringBuilder sb = new StringBuilder();
		sb.append("Before: '" + sequenceBeforeModification + "', '" + sequenceWithModification + "' and after '"
				+ getSequenceAfterModification() + "'");
		return sb.toString();
	}

	public double getSequenceMassAfterModification() {
		if (sequenceMassAfterModification == null) {
			sequenceMassAfterModification = extendedAssignMass.calculateMass(sequenceWithModification);
		}
		return sequenceMassAfterModification;
	}

	public String getProteinSequence() {
		return proteinSequence;
	}

	@Override
	public int hashCode() {
		return getSequenceWithModification().hashCode() + getSequenceAfterModification().hashCode();
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof SequenceWithModification) {
			return toString().equals(((SequenceWithModification) obj).toString());
		}
		return super.equals(obj);
	}
}
