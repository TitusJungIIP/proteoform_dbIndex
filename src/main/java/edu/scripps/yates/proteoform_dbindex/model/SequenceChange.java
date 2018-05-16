package edu.scripps.yates.proteoform_dbindex.model;

import java.text.DecimalFormat;

import org.apache.commons.lang.math.IntRange;

import edu.scripps.yates.annotations.uniprot.proteoform.ProteoformType;

/**
 * Class representing a change in the peptide sequence, stating in which
 * position (or range of positions), the sequence change (a substitution or a
 * PTM)
 * 
 * @author salvador
 *
 */
public class SequenceChange implements Comparable<SequenceChange> {
	private final IntRange positionInPeptide;
	private final String original;
	private final String change;
	private final boolean ptm;
	private final Double massChange;
	private String proteoformID;
	private final ProteoformType proteoformType;
	public final static DecimalFormat df = new DecimalFormat("+#.####;-#.####");

	public SequenceChange(String proteoformID, ProteoformType proteoformType, int startOfChange, String original,
			String change, Double massChange) {
		this.proteoformID = proteoformID;
		if (original != null) {
			positionInPeptide = new IntRange(startOfChange, startOfChange + original.length() - 1);
		} else {
			positionInPeptide = new IntRange(startOfChange);
		}
		this.original = original;
		this.change = change;
		ptm = massChange != null;
		this.massChange = massChange;
		this.proteoformType = proteoformType;

	}

	public SequenceChange(String proteoformID, ProteoformType proteoformType, IntRange positionInPeptide,
			String original, String change, Double massChange) {
		this.proteoformID = proteoformID;
		this.positionInPeptide = positionInPeptide;
		this.original = original;
		this.change = change;
		ptm = massChange != null;
		this.massChange = massChange;
		this.proteoformType = proteoformType;
		if (original == null && change == null && massChange == null) {
			System.out.println(this);
		}

	}

	public IntRange getPositionsInPeptide() {
		return positionInPeptide;
	}

	public String getOriginal() {
		return original;
	}

	public String getChange() {
		return change;
	}

	public boolean isPtm() {
		return ptm;
	}

	public Double getMassChange() {
		return massChange;
	}

	public String getProteoformID() {
		return proteoformID;
	}

	public void setProteoformID(String proteoformID) {
		this.proteoformID = proteoformID;
	}

	public ProteoformType getProteoformType() {
		return proteoformType;
	}

	public String getKey() {
		if (massChange != null) {
			return df.format(massChange);
		} else {
			return original + "->" + change;
		}
	}

	@Override
	public int compareTo(SequenceChange o) {
		int ret = Integer.compare(positionInPeptide.getMinimumInteger(), o.positionInPeptide.getMinimumInteger());
		if (ret != 0) {
			return ret;
		}
		if (ret == 0) {
			ret = Integer.compare(positionInPeptide.getMaximumInteger(), o.positionInPeptide.getMaximumInteger());
		}
		if (ret != 0) {
			return ret;
		}
		return proteoformID.compareTo(o.proteoformID);

	}

	public int getFirstPositionOfChangeInPeptide() {
		return positionInPeptide.getMinimumInteger();
	}

}
