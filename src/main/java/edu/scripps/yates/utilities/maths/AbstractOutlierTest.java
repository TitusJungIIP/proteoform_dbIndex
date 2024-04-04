package edu.scripps.yates.utilities.maths;

import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

public abstract class AbstractOutlierTest {
	protected static final Logger log = Logger.getLogger(AbstractOutlierTest.class);

	protected double[] array;
	protected final double valueToTest;
	protected final List<Double> populationValues;

	public abstract OutlierTestResult test();

	protected AbstractOutlierTest(double valueToTest, List<Double> populationValues) {
		List<Double> list = filteroutNonValidValues(populationValues);
		if (list.isEmpty()) {
			throw new IllegalArgumentException("Non valid values to test");
		}
		if (!populationValues.contains(valueToTest)) {
			populationValues.add(0, valueToTest);
		}
		array = new double[list.size()];
		int i = 0;
		for (Double d : list) {
			array[i++] = d;
		}
		this.valueToTest = valueToTest;
		this.populationValues = populationValues;

	}

	protected List<Double> filteroutNonValidValues(List<Double> populationValues) {
		List<Double> ret = new ArrayList<Double>();
		for (Double value : populationValues) {
			if (value == null || value.isInfinite() || value.isNaN()) {
				continue;
			}
			ret.add(value);
		}
		return ret;
	}

	public String getPopulationValuesString() {
		StringBuilder sb = new StringBuilder();
		for (double d : array) {
			if (!"".equals(sb.toString())) {
				sb.append(", ");
			}
			sb.append(d);
		}
		return sb.toString();
	}
}
