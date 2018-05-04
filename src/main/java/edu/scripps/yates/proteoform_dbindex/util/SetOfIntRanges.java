package edu.scripps.yates.proteoform_dbindex.util;

import java.util.HashSet;

import org.apache.commons.lang.math.IntRange;

public class SetOfIntRanges extends HashSet<IntRange> {

	/**
	 * 
	 */
	private static final long serialVersionUID = -4586056702636461047L;

	@Override
	public boolean contains(Object o) {
		if (o instanceof IntRange) {
			final IntRange range = (IntRange) o;

			for (final IntRange range2 : this) {
				if (range2.getMinimumInteger() >= range.getMinimumInteger()
						&& range2.getMinimumInteger() <= range.getMaximumInteger()) {
					return true;
				}
				if (range2.getMaximumInteger() >= range.getMinimumInteger()
						&& range2.getMaximumInteger() <= range.getMaximumInteger()) {
					return true;
				}
			}
			return false;
		}

		return super.contains(o);
	}
}
