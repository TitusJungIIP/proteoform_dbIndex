package edu.scripps.yates.utilities.venndata;

import java.util.Map;
import java.util.Set;

import gnu.trove.set.hash.THashSet;

public class VennDataUtils<T extends ContainsMultipleKeys> {

	protected VennDataUtils() {

	}

	protected Set<T> getObjectsByKeys(Set<String> keys, Map<String, T>... hashes) {
		Set<T> ret = new THashSet<T>();

		for (String key : keys) {
			for (Map<String, T> hash : hashes) {
				if (hash.containsKey(key)) {
					ret.add(hash.get(key));
					continue;
				}
			}
		}
		return ret;
	}

	public Set<T> getObjects(Map<String, T> hash) {
		Set<T> ret = new THashSet<T>();

		for (T obj : hash.values()) {
			ret.add(obj);
		}

		return ret;
	}

}
