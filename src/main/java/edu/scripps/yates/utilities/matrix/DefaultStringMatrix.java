package edu.scripps.yates.utilities.matrix;

public class DefaultStringMatrix extends Matrix<String, String> {

	public DefaultStringMatrix(Object emptyValue) {
		super(emptyValue);

	}

	@Override
	public String printValue(String t) {
		return t;
	}

}
