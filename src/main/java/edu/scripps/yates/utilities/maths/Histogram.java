package edu.scripps.yates.utilities.maths;

import gnu.trove.list.TDoubleList;

public class Histogram {

	public static double[][] calcHistogram(double[] data, double min, double max, int numBins) {
		final double[] breaks = smile.math.Histogram.breaks(min, max, numBins);
		return smile.math.Histogram.histogram(data, breaks);
	}

	public static double[][] calcHistogram(float[] data, float min, float max, int numBins) {
		final double[] breaks = smile.math.Histogram.breaks(min, max, numBins);
		final float[] breaksFloat = new float[breaks.length];
		for (int i = 0; i < breaks.length; i++) {
			breaksFloat[i] = Double.valueOf(breaks[i]).floatValue();
		}
		final double[][] histogram = smile.math.Histogram.histogram(data, breaksFloat);
		return histogram;

	}

	public static float getBinSize(float min, float max, int numBins) {
		return (max - min) / numBins;
	}

	public static double getBinSize(double min, double max, int numBins) {
		return (max - min) / numBins;
	}

	public static double getBinSize(TDoubleList data, int numBins) {
		return (data.max() - data.min()) / numBins;
	}

	/**
	 * This returns the Sturgis rule which determines the number of bins in an
	 * histogram as 1+3.3Log10n.<br>
	 * https://www.me.psu.edu/cimbala/me345/Lectures/Histograms.pdf#:~:text=an%20integer%20number%20of%20bins%2C%20we%20would%20use,10.9%20bins%2C%20which%20we%20round%20to%2011%20bins.
	 * 
	 * @param n the number of data points in the histogram
	 * @return
	 */
	public static int getSturgisRuleForHistogramBins(int n) {
		final int numberOfBins = Long.valueOf(Math.round(1 + 3.3 * Math.log10(n))).intValue();
		return numberOfBins;
	}

	/**
	 * This returns the Rice rule which determines the number of bins in an
	 * histogram as 2*n^(1/3)<br>
	 * https://www.me.psu.edu/cimbala/me345/Lectures/Histograms.pdf#:~:text=an%20integer%20number%20of%20bins%2C%20we%20would%20use,10.9%20bins%2C%20which%20we%20round%20to%2011%20bins.
	 * 
	 * @param n the number of data points in the histogram
	 * @return
	 */
	public static int getRiceRuleForHistogramBins(int n) {
		final int numberOfBins = Long.valueOf(Math.round(2 * Math.pow(n, 1.0 / 3.0))).intValue();
		return numberOfBins;
	}
}
