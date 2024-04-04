package edu.scripps.yates.utilities.maths;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;

public class PValuesCollection<T> {
	private final TObjectDoubleMap<T> pValues = new TObjectDoubleHashMap<T>();
	private List<T> sortedKeysByPValue;

	public PValuesCollection(TObjectDoubleMap<T> pValues) {
		this.pValues.putAll(pValues);
		process();
	}

	private void process() {
		sortedKeysByPValue = new ArrayList<T>();
		sortedKeysByPValue.addAll(pValues.keySet());
		// sort pvalues
		Collections.sort(sortedKeysByPValue, new Comparator<T>() {

			@Override
			public int compare(T o1, T o2) {

				final Double d1 = pValues.get(o1);
				final Double d2 = pValues.get(o2);
				final double num1 = d1 != null ? d1 : Double.MAX_VALUE;
				final double num2 = d2 != null ? d2 : Double.MAX_VALUE;
				return Double.compare(num1, num2);
			}
		});

	}

	public List<T> getSortedKeysByPValue() {
		return sortedKeysByPValue;
	}

	public int size() {
		return pValues.size();
	}

	public Double getPValue(T key) {
		if (pValues.containsKey(key)) {
			return pValues.get(key);
		} else {
			return null;
		}
	}

	public TObjectDoubleMap<T> getPValues() {
		return pValues;
	}

}
