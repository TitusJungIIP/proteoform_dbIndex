package edu.scripps.yates.utilities.collections;

import java.util.Collection;

public abstract class CollectionObserver<T> {

	public abstract boolean add(T obj);

	public abstract void clear();

	public abstract boolean remove(Object obj);

	public abstract boolean addAll(Collection<? extends T> collection);

	public abstract boolean removeAll(Collection<?> collection);

}
