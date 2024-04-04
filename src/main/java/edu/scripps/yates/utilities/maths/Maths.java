package edu.scripps.yates.utilities.maths;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Median;

import gnu.trove.list.TDoubleList;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.TLongList;
import gnu.trove.list.array.TDoubleArrayList;

public class Maths {
	private Maths() {
	}

	/**
	 * Returns the maximum value in the array a[], -infinity if no such value.
	 */
	public static double max(double[] a) {
		double max = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < a.length; i++) {
			if (Double.isNaN(a[i]))
				return Double.NaN;
			if (a[i] > max)
				max = a[i];
		}
		return max;
	}

	/**
	 * Returns the maximum value in the subarray a[lo..hi], -infinity if no such
	 * value.
	 */
	public static double max(double[] a, int lo, int hi) {
		if (lo < 0 || hi >= a.length || lo > hi)
			throw new RuntimeException("Subarray indices out of bounds");
		double max = Double.NEGATIVE_INFINITY;
		for (int i = lo; i <= hi; i++) {
			if (Double.isNaN(a[i]))
				return Double.NaN;
			if (a[i] > max)
				max = a[i];
		}
		return max;
	}

	public static double max(Double... doubles) {
		double max = -Double.MAX_VALUE;
		for (final Double double1 : doubles) {
			if (double1 != null) {
				if (max < double1) {
					max = double1;
				}
			}
		}
		return max;
	}

	public static float max(Float... numbers) {
		float max = -Float.MAX_VALUE;
		for (final Float number : numbers) {
			if (number != null) {
				if (max < number) {
					max = number;
				}
			}
		}
		return max;
	}

	public static double min(Double... doubles) {
		double min = Double.MAX_VALUE;
		for (final Double double1 : doubles) {
			if (double1 != null) {
				if (min > double1) {
					min = double1;
				}
			}
		}
		return min;
	}

	public static float min(float... numbers) {
		float min = Float.MAX_VALUE;
		for (final Float number : numbers) {
			if (min > number) {
				min = number;
			}
		}
		return min;
	}

	/**
	 * Returns the maximum value in the array a[], Integer.MIN_VALUE if no such
	 * value.
	 */
	public static int max(int[] a) {
		int max = Integer.MIN_VALUE;
		for (int i = 0; i < a.length; i++) {
			if (a[i] > max)
				max = a[i];
		}
		return max;
	}

	/**
	 * Returns the minimum value in the array a[], +infinity if no such value.
	 */
	public static double min(double[] a) {
		double min = Double.POSITIVE_INFINITY;
		for (int i = 0; i < a.length; i++) {
			if (Double.isNaN(a[i]))
				return Double.NaN;
			if (a[i] < min)
				min = a[i];
		}
		return min;
	}

	/**
	 * Returns the minimum value in the subarray a[lo..hi], +infinity if no such
	 * value.
	 */
	public static double min(double[] a, int lo, int hi) {
		if (lo < 0 || hi >= a.length || lo > hi)
			throw new RuntimeException("Subarray indices out of bounds");
		double min = Double.POSITIVE_INFINITY;
		for (int i = lo; i <= hi; i++) {
			if (Double.isNaN(a[i]))
				return Double.NaN;
			if (a[i] < min)
				min = a[i];
		}
		return min;
	}

	/**
	 * Returns the minimum value in the array a[], Integer.MAX_VALUE if no such
	 * value.
	 */
	public static int min(int[] a) {
		int min = Integer.MAX_VALUE;
		for (int i = 0; i < a.length; i++) {
			if (a[i] < min)
				min = a[i];
		}
		return min;
	}

	/**
	 * Returns the average value in the array a[], NaN if no such value.
	 */
	public static double mean(double[] a) {
		if (a.length == 0)
			return Double.NaN;
		final double sum = sum(a);
		return sum / a.length;
	}

	/**
	 * Returns the average value in the array a[], NaN if no such value.
	 */
	public static float mean(float[] a) {
		if (a.length == 0)
			return Float.NaN;
		final float sum = sum(a);
		return sum / a.length;
	}

	/**
	 * Returns the average value in the array a[], NaN if no such value.
	 */
	public static double mean(Double[] a) {
		if (a.length == 0)
			return Double.NaN;
		final double sum = sum(a);
		return sum / a.length;
	}

	/**
	 * Returns the average value in the array a[], NaN if no such value.
	 */
	public static double mean(TDoubleList a) {
		if (a.isEmpty())
			return Double.NaN;

		return a.sum() / a.size();
	}

	/**
	 * Returns the average value in the array a[], NaN if no such value.
	 */
	public static double mean(TIntList a) {
		if (a.isEmpty())
			return Double.NaN;

		return 1.0 * a.sum() / a.size();
	}

	/**
	 * Returns the average value in the array a[], NaN if no such value.
	 */
	public static double mean(TLongList a) {
		if (a.isEmpty())
			return Double.NaN;

		return 1.0 * a.sum() / a.size();
	}

	/**
	 * Returns the average value in the array a[], NaN if no such value.
	 */
	public static float mean(TFloatList a) {
		if (a.isEmpty())
			return Float.NaN;

		return a.sum() / a.size();
	}

	/**
	 * Returns the average value in the subarray a[lo..hi], NaN if no such value.
	 */
	public static double mean(double[] a, int lo, int hi) {
		final int length = hi - lo + 1;
		if (lo < 0 || hi >= a.length || lo > hi)
			throw new RuntimeException("Subarray indices out of bounds");
		if (length == 0)
			return Double.NaN;
		final double sum = sum(a, lo, hi);
		return sum / length;
	}

	/**
	 * Returns the average value in the array a[], NaN if no such value.
	 */
	public static double mean(int[] a) {
		if (a.length == 0)
			return Double.NaN;
		double sum = 0.0;
		for (int i = 0; i < a.length; i++) {
			sum = sum + a[i];
		}
		return sum / a.length;
	}

	/**
	 * Returns the sample variance in the array a[], NaN if no such value.
	 */
	public static double var(double[] a) {
		if (a.length == 0)
			return Double.NaN;
		final double avg = mean(a);
		double sum = 0.0;
		for (int i = 0; i < a.length; i++) {
			sum += (a[i] - avg) * (a[i] - avg);
		}
		return sum / (a.length - 1);
	}

	/**
	 * Returns the sample variance in the array a[], NaN if no such value.
	 */
	public static float var(float[] a) {
		if (a.length == 0)
			return Float.NaN;
		final float avg = mean(a);
		float sum = 0.0f;
		for (int i = 0; i < a.length; i++) {
			sum += (a[i] - avg) * (a[i] - avg);
		}
		return sum / (a.length - 1);
	}

	/**
	 * Returns the sample variance in the subarray a[lo..hi], NaN if no such value.
	 */
	public static double var(double[] a, int lo, int hi) {
		final int length = hi - lo + 1;
		if (lo < 0 || hi >= a.length || lo > hi)
			throw new RuntimeException("Subarray indices out of bounds");
		if (length == 0)
			return Double.NaN;
		final double avg = mean(a, lo, hi);
		double sum = 0.0;
		for (int i = lo; i <= hi; i++) {
			sum += (a[i] - avg) * (a[i] - avg);
		}
		return sum / (length - 1);
	}

	/**
	 * Returns the sample variance in the array a[], NaN if no such value.
	 */
	public static double var(int[] a) {
		if (a.length == 0)
			return Double.NaN;
		final double avg = mean(a);
		double sum = 0.0;
		for (int i = 0; i < a.length; i++) {
			sum += (a[i] - avg) * (a[i] - avg);
		}
		return sum / (a.length - 1);
	}

	/**
	 * Returns the population variance in the array a[], NaN if no such value.
	 */
	public static double varp(double[] a) {
		if (a.length == 0)
			return Double.NaN;
		final double avg = mean(a);
		double sum = 0.0;
		for (int i = 0; i < a.length; i++) {
			sum += (a[i] - avg) * (a[i] - avg);
		}
		return sum / a.length;
	}

	/**
	 * Returns the population variance in the subarray a[lo..hi], NaN if no such
	 * value.
	 */
	public static double varp(double[] a, int lo, int hi) {
		final int length = hi - lo + 1;
		if (lo < 0 || hi >= a.length || lo > hi)
			throw new RuntimeException("Subarray indices out of bounds");
		if (length == 0)
			return Double.NaN;
		final double avg = mean(a, lo, hi);
		double sum = 0.0;
		for (int i = lo; i <= hi; i++) {
			sum += (a[i] - avg) * (a[i] - avg);
		}
		return sum / length;
	}

	/**
	 * Returns the sample standard deviation in the array a[], NaN if no such value.
	 */
	public static double stddev(double[] a) {
		return Math.sqrt(var(a));
	}

	/**
	 * Returns the sample standard deviation in the array a[], NaN if no such value.
	 */
	public static double stddev(float[] a) {
		return Math.sqrt(var(a));
	}

	/**
	 * Returns the sample standard deviation in the array a[], NaN if no such value.
	 */
	public static double stddev(TDoubleList a) {
		return stddev(a.toArray());
	}

	/**
	 * Returns the sample standard deviation in the array a[], NaN if no such value.
	 */
	public static double stddev(TFloatList a) {
		return stddev(a.toArray());
	}

	/**
	 * Returns the sample standard deviation in the subarray a[lo..hi], NaN if no
	 * such value.
	 */
	public static double stddev(double[] a, int lo, int hi) {
		return Math.sqrt(var(a, lo, hi));
	}

	/**
	 * Returns the sample standard deviation in the array a[], NaN if no such value.
	 */
	public static double stddev(int[] a) {
		return Math.sqrt(var(a));
	}

	/**
	 * Returns the sample standard deviation in the array a[], NaN if no such value.
	 */
	public static double stddev(TIntList a) {
		return stddev(a.toArray());
	}

	/**
	 * Returns the population standard deviation in the array a[], NaN if no such
	 * value.
	 */
	public static double stddevp(double[] a) {
		return Math.sqrt(varp(a));
	}

	/**
	 * Returns the population standard deviation in the subarray a[lo..hi], NaN if
	 * no such value.
	 */
	public static double stddevp(double[] a, int lo, int hi) {
		return Math.sqrt(varp(a, lo, hi));
	}

	/**
	 * Returns the sum of all values in the {@link TDoubleList}
	 */
	public static double sum(TDoubleList a) {
		return a.sum();
	}

	/**
	 * Returns the sum of all values in the {@link TFloatList}
	 */
	public static double sum(TFloatList a) {
		return a.sum();
	}

	/**
	 * Returns the sum of all values in the array a[].
	 */
	public static double sum(double[] a) {
		double sum = 0.0;
		for (int i = 0; i < a.length; i++) {
			sum += a[i];
		}
		return sum;
	}

	/**
	 * Returns the sum of all values in the array a[].
	 */
	public static float sum(float[] a) {
		float sum = 0.0f;
		for (int i = 0; i < a.length; i++) {
			sum += a[i];
		}
		return sum;
	}

	public static double sum(Double[] a) {
		double sum = 0.0;
		for (int i = 0; i < a.length; i++) {
			if (a[i] != null) {
				sum += a[i];
			}
		}
		return sum;
	}

	/**
	 * Returns the sum of all values in the subarray a[lo..hi].
	 */
	public static double sum(double[] a, int lo, int hi) {
		if (lo < 0 || hi >= a.length || lo > hi)
			throw new RuntimeException("Subarray indices out of bounds");
		double sum = 0.0;
		for (int i = lo; i <= hi; i++) {
			sum += a[i];
		}
		return sum;
	}

	/**
	 * Returns the sum of all values in the array a[].
	 */
	public static int sum(int[] a) {
		int sum = 0;
		for (int i = 0; i < a.length; i++) {
			sum += a[i];
		}
		return sum;
	}

	/**
	 * Checks if the value is equals to Double.MAX_VALUE or Double.MIN_VALUE or the
	 * negative of these numbers
	 *
	 * @param value
	 * @return
	 */
	public static boolean isMaxOrMinValue(double value) {
		if (Double.compare(Double.MAX_VALUE, value) == 0 || Double.compare(Double.MIN_VALUE, value) == 0)
			return true;
		if (Double.compare(Double.MAX_VALUE, -value) == 0 || Double.compare(Double.MIN_VALUE, -value) == 0)
			return true;
		return false;
	}

	/**
	 * Calculates de median absolute deviation (MAD)
	 *
	 * @param values
	 * @return
	 */
	public static double mad(double[] values) {
		final Median medianCalculator = new Median();
		final double median = medianCalculator.evaluate(values);

		final double[] tmp = new double[values.length];
		for (int i = 0; i < values.length; i++) {
			final double value = values[i];
			tmp[i] = Math.abs(value - median);
		}
		final double ret = medianCalculator.evaluate(tmp);
		return ret;
	}

	/**
	 * This is an outliers test: Mi=0.6745(xi−x~)/MAD where MAD is the median
	 * absolute deviation x~ is the median of the population and xi is the value to
	 * test wether is an outlier or not.<br>
	 *
	 *
	 * @param valueToTest
	 * @param populationValues
	 * @return if the return value is greater than 3.5, then it should be considered
	 *         as a potential outlier.
	 */
	public static double iglewiczHoaglinTest(double valueToTest, double[] populationValues) {
		final double factor = 0.6745;
		final Median medianCalculator = new Median();
		final double populationMedian = medianCalculator.evaluate(populationValues);
		final double mad = mad(populationValues);

		final double ret = Math.abs(factor * (valueToTest - populationMedian) / mad);
		return ret;

	}

	//

	/**
	 * Function to calculate sample error or Standard Error of Measurement, that is
	 * the standard deviation divided by the squared root of the size of the sample
	 * 
	 * @param arr
	 * @return
	 */
	public static double sem(double arr[]) {
		if (arr.length == 0) {
			return 0.0;
		}
		// Formula to find sample error.
		return stddev(arr) / Math.sqrt(1.0 * arr.length);
	}

	/**
	 * Function to calculate sample error or Standard Error of Measurement, that is
	 * the standard deviation divided by the squared root of the size of the sample
	 * 
	 * @param values
	 * @return
	 */
	public static double sem(TDoubleList values) {
		return sem(values.toArray());
	}

	/**
	 * Z-score calculation: (value-mean)/stdev a z-score greater than 3 could be
	 * considered as an outlier
	 *
	 * @param valueToTest
	 * @param populationValues
	 * @return
	 */
	public static double zScore(double valueToTest, double[] populationValues) {

		final double populationMean = new Mean().evaluate(populationValues);
		final double populationSTD = new StandardDeviation().evaluate(populationValues);

		final double ret = (valueToTest - populationMean) / populationSTD;
		return ret;

	}

	public static void main(String[] args) {
		final double[] data = { 1.58, 0, 1.73 };

		System.out.println(iglewiczHoaglinTest(0, data));

	}

	public static int factorial(int n) {
		if (n >= 13) {
			throw new IllegalArgumentException("n cannot be greater than 13, otherwise the integer overflows");
		}
		if (n == 0) {
			return 1;
		} else {
			return n * factorial(n - 1);
		}
	}

	public static double log(double x, int base) {
		return Math.log(x) / Math.log(base);
	}

	public static double correlationCoefficient(TDoubleList x, TDoubleList y) {
		// double sum_X = 0, sum_Y = 0, sum_XY = 0;
		// double squareSum_X = 0, squareSum_Y = 0;
		// final int n = x.size();
		// for (int i = 0; i < n; i++) {
		// // sum of elements of array X.
		// sum_X = sum_X + x.get(i);
		//
		// // sum of elements of array Y.
		// sum_Y = sum_Y + y.get(i);
		//
		// // sum of X[i] * Y[i].
		// sum_XY = sum_XY + x.get(i) * y.get(i);
		//
		// // sum of square of array elements.
		// squareSum_X = squareSum_X + x.get(i) * x.get(i);
		// squareSum_Y = squareSum_Y + y.get(i) * y.get(i);
		// }
		// // use formula for calculating correlation
		// // coefficient.
		// final double corr = (n * sum_XY - sum_X * sum_Y)
		// / Math.sqrt((n * squareSum_X - sum_X * sum_X) * (n * squareSum_Y -
		// sum_Y * sum_Y));
		//
		// return corr;

		// replaced by Maths3 from apache
		return new PearsonsCorrelation().correlation(x.toArray(), y.toArray());
	}

	/**
	 * Coefficient of variation: 100 * stdev / mean
	 * 
	 * @param values
	 * @return
	 */
	public static double cv(TDoubleList values) {
		final double stddev = stddev(values);
		final double mean = mean(values);
		final double cv = 100.0 * stddev / mean;
		return cv;
	}

	public static double median(TDoubleList values) {
		final TDoubleList sortedList = new TDoubleArrayList();
		sortedList.addAll(values);
		sortedList.sort();
		final int middle = sortedList.size() / 2;
		final double median = sortedList.get(middle);
		return median;
	}

	public static double dotProduct(TDoubleList values1, TDoubleList values2) {
		return dotProduct(values1.toArray(), values2.toArray());
	}

	public static double dotProduct(double[] values1, double[] values2) {
		if (values1.length != values2.length) {
			throw new IllegalArgumentException("Both arrays have to have the same size");
		}
		double product = 0;

		// Loop for calculate cot product
		for (int i = 0; i < values1.length; i++)
			product = product + values1[i] * values2[i];
		return product;
	}

	/**
	 * It normalizes the vector of numbers by dividing all of them by the maximum
	 * 
	 * @param numbers
	 * @return
	 */
	public static TDoubleList normalize(TDoubleList numbers) {
		final double max = numbers.max();
		final TDoubleList ret = new TDoubleArrayList();
		for (final double number : numbers.toArray()) {
			ret.add(number / max);
		}
		return ret;
	}

	/**
	 * Calculates the cosine similarity between two vectors of the same size
	 * 
	 * @param vectorA
	 * @param vectorB
	 * @return
	 */
	public static double cosineSimilarity(double[] vectorA, double[] vectorB) {
		if (vectorA.length != vectorB.length) {
			throw new IllegalArgumentException("Vectors must have same length");
		}
		if (vectorA.length == 0) {
			throw new IllegalArgumentException("vector A is empty");
		}
		if (vectorB.length == 0) {
			throw new IllegalArgumentException("vector B is empty");
		}
		double dotProduct = 0.0;
		double normA = 0.0;
		double normB = 0.0;
		for (int i = 0; i < vectorA.length; i++) {
			dotProduct += vectorA[i] * vectorB[i];
			normA += Math.pow(vectorA[i], 2);
			normB += Math.pow(vectorB[i], 2);
		}
		return dotProduct / (Math.sqrt(normA) * Math.sqrt(normB));
	}
}
