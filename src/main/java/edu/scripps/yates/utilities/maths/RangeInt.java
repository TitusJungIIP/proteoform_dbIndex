package edu.scripps.yates.utilities.maths;

import java.util.ArrayList;
import java.util.List;

import gnu.trove.TIntCollection;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;

public class RangeInt {
	private final int minInclusive;
	private final int maxInclusive;
	private TIntArrayList tintList;
	private ArrayList<Integer> list;

	public RangeInt(int minInclusive, int maxInclusibe) {
		super();
		this.minInclusive = minInclusive;
		this.maxInclusive = maxInclusibe;
	}

	public int getMinimumInclusive() {
		return minInclusive;
	}

	public int getMaximumInclusive() {
		return maxInclusive;
	}

	public List<Integer> getListOfTotalRangeIntegers() {
		if (list == null) {
			list = new ArrayList<Integer>();
			for (int i = minInclusive; i <= maxInclusive; i++) {
				list.add(i);
			}
		}
		return list;
	}

	public TIntList getTIntListOfTotalRangeIntegers() {
		if (tintList == null) {
			tintList = new TIntArrayList();
			for (int i = minInclusive; i <= maxInclusive; i++) {
				tintList.add(i);
			}
		}
		return tintList;
	}

	public static List<RangeInt> createListOfRangeInts(TIntCollection tintCollection) {
		final TIntList list = new TIntArrayList();
		list.addAll(tintCollection);
		list.sort();
		final List<RangeInt> ret = new ArrayList<RangeInt>();
		final TIntList numbers = new TIntArrayList();
		for (int i = 0; i < list.size(); i++) {
			final int number = list.get(i);
			if (numbers.isEmpty()) {
				numbers.add(number);
			} else {
				// check if this number is equal to the previous one +1
				if (number != numbers.get(numbers.size() - 1) + 1) {
					final RangeInt range = new RangeInt(numbers.min(), numbers.max());
					ret.add(range);
					numbers.clear();
				}
				numbers.add(number);
			}
		}
		if (!numbers.isEmpty()) {
			final RangeInt range = new RangeInt(numbers.min(), numbers.max());
			ret.add(range);
		}
		return ret;
	}
}
