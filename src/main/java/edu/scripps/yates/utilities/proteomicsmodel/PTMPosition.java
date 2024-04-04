package edu.scripps.yates.utilities.proteomicsmodel;

public enum PTMPosition {

	CTERM, NTERM, PCTERM, PNTERM, NONE;
	public static PTMPosition parsePositon(String s) {
		if (s != null) {
			if (s.compareToIgnoreCase("N-Term") == 0)
				return PTMPosition.NTERM;
			if (s.compareToIgnoreCase("C-Term") == 0)
				return PTMPosition.CTERM;
			if (s.compareToIgnoreCase("Protein N-term") == 0)
				return PTMPosition.PNTERM;
			if (s.compareToIgnoreCase("Protein C-Term") == 0)
				return PTMPosition.PCTERM;
		}
		return PTMPosition.NONE;
	}

	public static PTMPosition getPTMPositionFromSequence(String sequence, int position) {
		if (position == 0) {
			return PTMPosition.NTERM;
		}
		if (position == sequence.length() + 1) {
			return PTMPosition.CTERM;
		}
		return PTMPosition.NONE;
	}
}
