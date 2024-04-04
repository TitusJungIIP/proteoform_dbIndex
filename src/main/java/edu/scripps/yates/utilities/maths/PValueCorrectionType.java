package edu.scripps.yates.utilities.maths;

public enum PValueCorrectionType {
	BH("Yoav Benjamini, Yosef Hochberg \"Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing\", Journal of the Royal Statistical Society. Series B, Vol. 57, No. 1 (1995), pp. 289-300, JSTOR:2346101"), //
	BY("Yoav Benjamini, Daniel Yekutieli, \"The control of the false discovery rate in multiple testing under dependency\", Ann. Statist., Vol. 29, No. 4 (2001), pp. 1165-1188, DOI:10.1214/aos/1013699998 JSTOR:2674075"), //
	BONFERRONI("Reference not available"), //
	HOCHBERG(
			"Yosef Hochberg, \"A sharper Bonferroni procedure for multiple tests of significance\", Biometrika, Vol. 75, No. 4 (1988), pp 800–802, DOI:10.1093/biomet/75.4.800 JSTOR:2336325"), //
	HOLM("Sture Holm, \"A Simple Sequentially Rejective Multiple Test Procedure\", Scandinavian Journal of Statistics, Vol. 6, No. 2 (1979), pp. 65-70, JSTOR:4615733"), //
	HOMMEL("Gerhard Hommel, \"A stagewise rejective multiple test procedure based on a modified Bonferroni test\", Biometrika, Vol. 75, No. 2 (1988), pp 383–386, DOI:10.1093/biomet/75.2.383 JSTOR:2336190");

	private final String reference;

	private PValueCorrectionType(String reference) {
		this.reference = reference;
	}

	public String getReference() {
		return reference;
	}

	public static String getValuesString() {
		final StringBuilder sb = new StringBuilder();
		for (final PValueCorrectionType pValueCorrectionType : PValueCorrectionType.values()) {
			if (!"".equals(sb.toString())) {
				sb.append(",");
			}
			sb.append(pValueCorrectionType.name());
		}
		return sb.toString();
	}
}
