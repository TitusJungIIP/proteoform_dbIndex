package edu.scripps.yates.utilities.maths;

import java.util.List;

public class OutlierTestResult {
	private boolean outlier = false;

	public OutlierTestResult(double valueToTest, List<Double> populationValues, List<Object> outlierResult) {

		final boolean isOutlier = (boolean) outlierResult.get(0);
		if (isOutlier) {
			final int[] indexOfOutliers = (int[]) outlierResult.get(2);
			for (int index : indexOfOutliers) {
				if (index == 0) {
					outlier = true;
				}
			}
		}
	}

	/**
	 * @return the outlier
	 */
	public boolean isOutlier() {
		return outlier;
	}

}
