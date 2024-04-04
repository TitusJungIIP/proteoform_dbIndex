package edu.scripps.yates.utilities.matrix;

public class DefaultIntegerMatrix extends Matrix<String, Integer> {

	public DefaultIntegerMatrix(Object emptyValue) {
		super(emptyValue);
	}

	@Override
	public String printValue(Integer t) {
		return String.valueOf(t);
	}

}
