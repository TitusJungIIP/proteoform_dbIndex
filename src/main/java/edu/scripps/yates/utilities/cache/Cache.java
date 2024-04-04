package edu.scripps.yates.utilities.cache;

import java.util.Collection;
import java.util.Map;
import java.util.Set;

public interface Cache<T, V> {

	public void addtoCache(T t, V key);

	public void addtoCache(Map<V, T> map);

	public T getFromCache(V key);

	public Set<T> getFromCache(Collection<V> keys);

	public boolean contains(V key);

	public T removeFromCache(V key);

	boolean containsAll(Collection<V> keys);

	public V processKey(V key);

	public void clearCache();

	public boolean isEmpty();
}
