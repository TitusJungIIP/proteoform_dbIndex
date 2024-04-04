package edu.scripps.yates.utilities.maths;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import gnu.trove.set.hash.TDoubleHashSet;

public class DixonOutlierQTest {
	/**
	 * Perform a Dixon's qTest over a {@link Collection} of double numbers
	 *
	 * @param numbers
	 * @return
	 */
	public boolean qTest(Collection<Double> numbers) {
		if (numbers.size() < 2) {
			return false;
		}
		List<Double> numberList = new ArrayList<Double>();
		numberList.addAll(numbers);
		TDoubleHashSet numberSet = new TDoubleHashSet();
		numberSet.addAll(numbers);
		if (numberSet.size() < 3) {
			return false;
		}

		// sorts the ratios to get the largest and smallest
		Collections.sort(numberList);

		Double largest = numberList.get(numberList.size() - 1);
		Double smallest = numberList.get(0);

		// finds the nearest number to the largest number
		int substractor = 2;
		Double nearestNumLarge = numberList.get(numbers.size() - substractor);
		// takes into account for duplicates
		while (Double.compare(nearestNumLarge, largest) == 0) {
			substractor++;
			if (numbers.size() <= substractor) {
				return false;
			}
			nearestNumLarge = numberList.get(numbers.size() - substractor);
		}

		// finds the nearest number to the smallest number
		int summator = 1;
		Double nearestNumSmall = numberList.get(summator);
		// takes into account for duplicates
		while (Double.compare(nearestNumSmall, smallest) == 0) {
			summator++;
			if (numbers.size() - 2 == summator) {
				return false;
			}
			nearestNumSmall = numberList.get(summator);
		}

		// qTest calculation
		Double qResult1 = Math.abs(smallest - nearestNumSmall) / Math.abs(largest - smallest);
		Double qResult2 = Math.abs(largest - nearestNumLarge) / Math.abs(largest - smallest);

		// sees if the qtest result falls below or above the threshold
		boolean flag = qTestThreshold(qResult1, qResult2, numberList);

		// true or false, is it inconsistent of consistent
		return flag;
	}

	// qTest threshold comparer
	private static boolean qTestThreshold(Double qResult1, Double qResult2, List<Double> ratios) {
		if (ratios.size() < 3) {
			return false;
		}

		if (ratios.size() == 3) {
			if (qResult1 > 0.97 || qResult2 > 0.97) {
				return true;
			}
		}
		return false;
	}
}
