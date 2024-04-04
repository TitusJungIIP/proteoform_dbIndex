package edu.scripps.yates.utilities.collections;

import java.util.HashSet;
import java.util.Set;

import gnu.trove.impl.hash.TObjectHash;
import gnu.trove.iterator.hash.TObjectHashIterator;

public class ObservableTObjectHashIterator<T> extends TObjectHashIterator<T> {
	private final Set<CollectionObserver<T>> collectionObservers = new HashSet<CollectionObserver<T>>();

	public ObservableTObjectHashIterator(TObjectHash<T> hash, CollectionObserver<T> collectionObserver) {
		super(hash);
		collectionObservers.add(collectionObserver);
	}

	public ObservableTObjectHashIterator(TObjectHash<T> hash, Set<CollectionObserver<T>> collectionObservers) {
		super(hash);
		this.collectionObservers.addAll(collectionObservers);
	}

	@Override
	public void remove() {
		for (final CollectionObserver<T> collectionListener : collectionObservers) {
			final T objectAtIndex = super.objectAtIndex(_index);
			collectionListener.remove(objectAtIndex);
		}
		super.remove();
	}

}
