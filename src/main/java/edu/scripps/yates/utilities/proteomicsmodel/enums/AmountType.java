package edu.scripps.yates.utilities.proteomicsmodel.enums;

public enum AmountType {

	INTENSITY, NORMALIZED_INTENSITY, AREA, XIC, SPC, NSAF, dNSAF, NSAF_NORM, EMPAI, EMPAI_COV, OTHER, REGRESSION_FACTOR,
	CORRIONINJECTION_INTENSITY, CCS;

	public static AmountType translateStringToAmountType(String amountTypeString) {

		final AmountType[] values = AmountType.values();
		for (final AmountType amountType : values) {
			if (amountType.name().equalsIgnoreCase(amountTypeString))
				return amountType;
		}
		return null;

	}

	public static String getValuesString() {
		final StringBuilder sb = new StringBuilder();
		final AmountType[] values = AmountType.values();
		for (final AmountType amountType : values) {
			if (!"".equals(sb.toString())) {
				sb.append(",");
			}
			sb.append(amountType.name());
		}
		return sb.toString();
	}

}
