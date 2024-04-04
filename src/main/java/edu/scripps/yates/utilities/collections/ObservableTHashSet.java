package edu.scripps.yates.utilities.collections;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import gnu.trove.iterator.hash.TObjectHashIterator;
import gnu.trove.set.hash.THashSet;

public class ObservableTHashSet<T> extends THashSet<T> {
	private final Set<CollectionObserver<T>> observers = new HashSet<CollectionObserver<T>>();

	public ObservableTHashSet() {

	}

	public ObservableTHashSet(CollectionObserver<T> observer) {
		addCollectionObserver(observer);
	}

	public void addCollectionObserver(CollectionObserver<T> observer) {
		observers.add(observer);
	}

	@Override
	public boolean add(T obj) {
		for (final CollectionObserver<T> observer : observers) {
			observer.add(obj);
		}
		return super.add(obj);
	}

	@Override
	public void clear() {
		for (final CollectionObserver<T> observer : observers) {
			observer.clear();
		}
		super.clear();
	}

	@Override
	public boolean remove(Object obj) {
		for (final CollectionObserver<T> observer : observers) {
			observer.remove(obj);
		}
		return super.remove(obj);
	}

	@Override
	public TObjectHashIterator<T> iterator() {
		return new ObservableTObjectHashIterator<T>(this, observers);
	}

	@Override
	public boolean addAll(Collection<? extends T> collection) {
		for (final CollectionObserver<T> observer : observers) {
			observer.addAll(collection);
		}

		return super.addAll(collection);
	}

	@Override
	public boolean removeAll(Collection<?> collection) {
		for (final CollectionObserver<T> observer : observers) {
			observer.removeAll(collection);
		}

		return super.removeAll(collection);
	}

}
