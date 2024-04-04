package edu.scripps.yates.utilities.index;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.utilities.util.Pair;

public interface FileIndex<T> {
	public T getItem(String key);

	public List<T> getItems(Collection<String> keys);

	public Map<String, Pair<Long, Long>> addItem(T item, Set<String> keys);

	public boolean isEmpty();
}
