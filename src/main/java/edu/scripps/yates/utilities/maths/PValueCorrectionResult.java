package edu.scripps.yates.utilities.maths;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;

public class PValueCorrectionResult<T> {
	private PValuesCollection<T> originalPValues;
	private PValuesCollection<T> correctedPValues;

	/**
	 * @return the originalPValues
	 */
	public PValuesCollection<T> getOriginalPValues() {
		return originalPValues;
	}

	/**
	 * @param originalPValues the originalPValues to set
	 */
	public void setOriginalPValues(PValuesCollection<T> originalPValues) {
		this.originalPValues = originalPValues;
	}

	/**
	 * @param originalPValues the originalPValues to set
	 */
	public void setOriginalPValues(TObjectDoubleMap<T> originalPValues) {
		setOriginalPValues(new PValuesCollection<T>(originalPValues));
	}

	/**
	 * @return the correctedPValues
	 */
	public PValuesCollection<T> getCorrectedPValues() {
		return correctedPValues;
	}

	/**
	 * @param correctedPValues the correctedPValues to set
	 */
	public void setCorrectedPValues(PValuesCollection<T> correctedPValues) {
		this.correctedPValues = correctedPValues;
	}

	/**
	 * @param correctedPValues the correctedPValues to set
	 */
	public void setCorrectedPValues(TObjectDoubleHashMap<T> pValues) {
		setCorrectedPValues(new PValuesCollection<T>(pValues));
	}

	public List<T> getSortedKeysByCorrectedPValue() {
		return getCorrectedPValues().getSortedKeysByPValue();
	}

	public List<T> getSortedKeysByOriginalPValue() {
		return getOriginalPValues().getSortedKeysByPValue();
	}

	public List<T> getSortedKeysByKeys(Comparator<T> comparator) {
		final List<T> keys = getSortedKeysByCorrectedPValue();
		Collections.sort(keys, comparator);
		return keys;
	}

}
