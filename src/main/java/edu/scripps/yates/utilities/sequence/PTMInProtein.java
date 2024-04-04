package edu.scripps.yates.utilities.sequence;

import java.text.DecimalFormat;

/**
 * A class to represent a PTM in a protein, including the deltaMass
 * 
 * @author salvador
 *
 */
public class PTMInProtein extends PositionInProtein {

	private final Double deltaMass;
	private static final DecimalFormat formatter = new DecimalFormat("+#.###;-#.###");

	public PTMInProtein(int position, char aa, String proteinACC, Double deltaMass) {
		super(position, aa, proteinACC);
		this.deltaMass = deltaMass;

	}

	public Double getDeltaMass() {
		return deltaMass;
	}

	@Override
	public String toString() {
		String string = getProteinACC() + SEPARATOR + getAa() + getPosition();
		if (deltaMass != null) {
			string += "(" + formatter.format(deltaMass) + ")";
		}
		return string;
	}

}
