package edu.scripps.yates.proteoform_dbindex.model;

public class PTMCodeObj {
	public static final String FILE_NAME = "ptmCodes.txt";
	private final String description;
	private final short ptmCode;
	private final Double ptmMassDiff;

	public PTMCodeObj(String description, short ptmCode, Double ptmMassDiff) {
		super();
		this.description = description;
		this.ptmCode = ptmCode;
		this.ptmMassDiff = ptmMassDiff;
	}

	public short getPtmCode() {
		return ptmCode;
	}

	public Double getPtmMassDiff() {
		return ptmMassDiff;
	}

	public String getDescription() {
		return description;
	}

	@Override
	public String toString() {

		String string = description + "\t" + ptmCode;
		if (ptmMassDiff != null) {
			string += "\t" + ptmMassDiff;
		}
		return string;
	}

	public static PTMCodeObj getFromLine(String line) {
		final String[] split = line.split("\t");
		final String description = split[0];
		final short code = Short.valueOf(split[1]);
		Double massDiff = null;
		if (split.length > 2) {
			massDiff = Double.valueOf(split[2]);
		}
		final PTMCodeObj ret = new PTMCodeObj(description, code, massDiff);
		return ret;

	}
}
