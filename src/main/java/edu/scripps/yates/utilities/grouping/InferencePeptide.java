package edu.scripps.yates.utilities.grouping;

import java.util.ArrayList;
import java.util.List;

public class InferencePeptide {
	private final List<InferenceProtein> inferenceProteins = new ArrayList<InferenceProtein>();
	private List<GroupablePeptide> mergedPeptides = new ArrayList<GroupablePeptide>();
	private PeptideRelation relation;
	private final String id;
	private final String sequence;

	public InferencePeptide(GroupablePeptide pept) {
		this(pept, PeptideRelation.NONDISCRIMINATING);
	}

	public InferencePeptide(GroupablePeptide pept, PeptideRelation r) {

		relation = r;
		mergedPeptides.add(pept);
		id = pept.getIdentifier();
		sequence = pept.getSequence();
	}

	public void addPeptide(GroupablePeptide p) {
		mergedPeptides.add(p);
	}

	public List<GroupablePeptide> getPeptidesMerged() {
		return mergedPeptides;
	}

	@Override
	public String toString() {
		switch (relation) {
		case DISCRIMINATING:
			return getId() + "*";
		case NONDISCRIMINATING:
			return getId() + "**";
		}
		return new Integer(getId()).toString();
	}

	private String getId() {
		return id;
	}

	public List<InferenceProtein> getInferenceProteins() {
		return inferenceProteins;
	}

	public List<GroupablePeptide> getMergedPeptides() {
		return mergedPeptides;
	}

	public void setMergedPeptides(List<GroupablePeptide> mergedPeptides) {
		this.mergedPeptides = mergedPeptides;
	}

	public PeptideRelation getRelation() {
		return relation;
	}

	public void setRelation(PeptideRelation relation) {
		for (final GroupablePeptide p : mergedPeptides) {
			p.setRelation(relation);
		}
		this.relation = relation;
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof InferencePeptide))
			return super.equals(obj);
		else {
			final InferencePeptide peptide = (InferencePeptide) obj;
			if (peptide.getSequence() == null) {
				return false;
			}
			return peptide.getSequence().equals(getSequence());
		}
	}

	private String getSequence() {
		return sequence;
	}
}
