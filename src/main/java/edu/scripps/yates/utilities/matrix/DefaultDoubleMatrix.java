package edu.scripps.yates.utilities.matrix;

public class DefaultDoubleMatrix extends Matrix<String, Double> {

	public DefaultDoubleMatrix(Object emptyValue) {
		super(emptyValue);

	}

	@Override
	public String printValue(Double t) {
		return String.valueOf(t);
	}

}
