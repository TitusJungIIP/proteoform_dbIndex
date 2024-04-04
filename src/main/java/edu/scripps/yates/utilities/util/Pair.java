package edu.scripps.yates.utilities.util;

public class Pair<T, V> {
	private T firstElement;
	private V secondElement;

	public Pair() {

	}

	public Pair(T firstElement, V secondElement) {
		this.firstElement = firstElement;
		this.secondElement = secondElement;
	}

	public T getFirstelement() {
		return this.firstElement;
	}

	public V getSecondElement() {
		return this.secondElement;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "[" + firstElement + "," + secondElement + "]";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#hashCode()
	 */
	@Override
	public int hashCode() {
		int hash = 17;
		hash = hash * 31 + firstElement.hashCode();
		hash = hash * 31 + secondElement.hashCode();
		return hash;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof Pair) {
			final Pair pair2 = (Pair) obj;
			if (this.firstElement.equals(pair2.firstElement) && this.secondElement.equals(pair2.secondElement)) {
				return true;
			}
			return false;
		}
		return super.equals(obj);
	}

}
