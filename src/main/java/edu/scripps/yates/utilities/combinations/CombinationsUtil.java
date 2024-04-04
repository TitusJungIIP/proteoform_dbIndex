package edu.scripps.yates.utilities.combinations;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.util.CombinatoricsUtils;

public class CombinationsUtil<T> {

	/**
	 * Returns all the combinations of N elements from the list of objects
	 * 
	 * @param n
	 * @param objs
	 * @return
	 */
	public List<List<T>> getCollectionsOfNElements(int n, List<T> objs) {
		final List<List<T>> ret = new ArrayList<List<T>>();
		final Iterator<int[]> combinationsIterator = CombinatoricsUtils.combinationsIterator(objs.size(), n);
		while (combinationsIterator.hasNext()) {
			final int[] combinationsIndexes = combinationsIterator.next();
			final List<T> list = new ArrayList<T>();
			for (final int index : combinationsIndexes) {
				list.add(objs.get(index));
			}
			ret.add(list);
		}
		return ret;
	}

}
