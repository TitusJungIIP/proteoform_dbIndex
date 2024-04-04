package edu.scripps.yates.utilities.staticstorage;

import java.util.Collections;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.locks.ReentrantReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock.ReadLock;
import java.util.concurrent.locks.ReentrantReadWriteLock.WriteLock;

import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.THashSet;

public class ItemStorage<T> {
	private final Map<String, Set<T>> itemByRunID = new THashMap<String, Set<T>>();
	private final ReentrantReadWriteLock itemByRunIDLock = new ReentrantReadWriteLock();

	private final Map<String, Set<T>> itemByConditionID = new THashMap<String, Set<T>>();
	private final ReentrantReadWriteLock itemByConditionIDLock = new ReentrantReadWriteLock();

	private final TIntObjectHashMap<Set<T>> itemByExcelRow = new TIntObjectHashMap<Set<T>>();
	private final ReentrantReadWriteLock itemByExcelRowLock = new ReentrantReadWriteLock();

	private final Map<String, Set<T>> itemByKey = new THashMap<String, Set<T>>();
	private final ReentrantReadWriteLock itemByKeyLock = new ReentrantReadWriteLock();

	public void clearData() {
		WriteLock wl = itemByConditionIDLock.writeLock();
		try {
			wl.lock();
			itemByConditionID.clear();
		} finally {
			wl.unlock();
		}
		wl = itemByExcelRowLock.writeLock();
		try {
			wl.lock();
			itemByExcelRow.clear();
		} finally {
			wl.unlock();
		}
		wl = itemByRunIDLock.writeLock();
		try {
			wl.lock();
			itemByRunID.clear();
		} finally {
			wl.unlock();
		}
		wl = itemByKeyLock.writeLock();
		try {
			wl.lock();
			itemByKey.clear();
		} finally {
			wl.unlock();
		}
	}

	/**
	 *
	 * @param item
	 * @param msRunID
	 * @param conditionID
	 * @param excelRowIndex
	 * @param key
	 *            it can be a sequence, or a protein accession
	 */
	public void add(T item, String msRunID, String conditionID, int excelRowIndex, String key) {
		if (item != null) {
			if (msRunID != null && !"".equals(msRunID)) {
				WriteLock writeLock = itemByRunIDLock.writeLock();
				try {
					writeLock.lock();
					addToMap(item, itemByRunID, msRunID);
				} finally {
					writeLock.unlock();
				}
			}
			if (conditionID != null && !"".equals(conditionID)) {
				WriteLock writeLock = itemByConditionIDLock.writeLock();
				try {
					writeLock.lock();
					addToMap(item, itemByConditionID, conditionID);
				} finally {
					writeLock.unlock();
				}
			}
			if (excelRowIndex > -1) {
				WriteLock writeLock = itemByExcelRowLock.writeLock();
				try {
					writeLock.lock();
					addToMap(item, itemByExcelRow, excelRowIndex);
				} finally {
					writeLock.unlock();
				}
			}
			if (key != null && !"".equals(key)) {
				WriteLock writeLock = itemByKeyLock.writeLock();
				try {
					writeLock.lock();
					addToMap(item, itemByKey, key);
				} finally {
					writeLock.unlock();
				}
			}
		}
	}

	/**
	 * Gets an item from the storage. Any non null input parameter will be
	 * checked and if not found, an empty set will be returned.<br>
	 * >If several non null parameters are submitted, it will return the
	 * intersection of the items retrieved from each of them.
	 *
	 * @param msRunID
	 * @param conditionID
	 * @param excelRowIndex
	 * @param key
	 * @return
	 */
	public Set<T> get(String msRunID, String conditionID, int excelRowIndex, String key) {
		Set<T> ret = new THashSet<T>();
		if (key != null && !"".equals(key)) {
			ReadLock readLock = itemByKeyLock.readLock();
			try {
				readLock.lock();
				if (itemByKey.containsKey(key)) {
					if (ret.isEmpty()) {
						ret.addAll(itemByKey.get(key));
					} else {
						ret = getIntersection(ret, itemByKey.get(key));
					}
				} else {
					return Collections.EMPTY_SET;
				}
			} finally {
				readLock.unlock();
			}
			if (ret.isEmpty()) {
				return Collections.EMPTY_SET;
			}
		}
		if (excelRowIndex > -1) {
			ReadLock readLock = itemByExcelRowLock.readLock();
			try {
				readLock.lock();
				if (itemByExcelRow.containsKey(excelRowIndex)) {
					if (ret.isEmpty()) {
						ret.addAll(itemByExcelRow.get(excelRowIndex));
					} else {
						ret = getIntersection(ret, itemByExcelRow.get(excelRowIndex));
					}
				} else {
					return Collections.EMPTY_SET;
				}
			} finally {
				readLock.unlock();
			}
			if (ret.isEmpty()) {
				return Collections.EMPTY_SET;
			}
		}

		if (msRunID != null && !"".equals(msRunID)) {
			ReadLock readLock = itemByRunIDLock.readLock();
			try {
				readLock.lock();
				if (itemByRunID.containsKey(msRunID)) {
					if (ret.isEmpty()) {
						ret.addAll(itemByRunID.get(msRunID));
					} else {
						ret = getIntersection(ret, itemByRunID.get(msRunID));
					}
				} else {
					return Collections.EMPTY_SET;
				}
			} finally {
				readLock.unlock();
			}
			if (ret.isEmpty()) {
				return Collections.EMPTY_SET;
			}
		}
		// NEW in 4th Nov 2016
		// only in case of not having runID, we look for condition
		// this is because if someone looks for a protein with a certain runID,
		// it doesnt matter the condition
		// only in case of someone looking for all proteins in a condition (no
		// matter the msrun) this would be useful
		if (msRunID == null) {
			if (conditionID != null && !"".equals(conditionID)) {
				ReadLock readLock = itemByConditionIDLock.readLock();
				try {
					readLock.lock();
					if (itemByConditionID.containsKey(conditionID)) {
						if (ret.isEmpty()) {
							ret.addAll(itemByConditionID.get(conditionID));
						} else {
							ret = getIntersection(ret, itemByConditionID.get(conditionID));
						}
					} else {
						return Collections.EMPTY_SET;
					}
				} finally {
					readLock.unlock();
				}
				if (ret.isEmpty()) {
					return Collections.EMPTY_SET;
				}
			}
		}
		return ret;
	}

	/**
	 * Gets a set containing the elements present in both sets
	 *
	 * @param ret
	 * @param set
	 * @return
	 */
	private Set<T> getIntersection(Set<T> set1, Set<T> set2) {
		Set<T> ret = new THashSet<T>();
		Set<T> smallerSet = set1;
		Set<T> biggerSet = set2;
		if (set2.size() < set1.size()) {
			smallerSet = set2;
			biggerSet = set1;
		}
		for (T t1 : smallerSet) {
			if (biggerSet.contains(t1)) {
				ret.add(t1);
			}
		}
		return ret;
	}

	private void addToMap(T item, Map<String, Set<T>> map, String key) {
		if (map.containsKey(key)) {
			map.get(key).add(item);
		} else {
			Set<T> set = new THashSet<T>();
			set.add(item);
			map.put(key, set);
		}
	}

	private void addToMap(T item, TIntObjectHashMap<Set<T>> map, int key) {
		if (map.containsKey(key)) {
			map.get(key).add(item);
		} else {
			Set<T> set = new THashSet<T>();
			set.add(item);
			map.put(key, set);
		}
	}

	public boolean contains(String msRunID, String conditionID, int excelRowIndex, String key) {
		return !get(msRunID, conditionID, excelRowIndex, key).isEmpty();
	}

	public boolean isEmpty() {
		ReadLock readLock = itemByConditionIDLock.readLock();
		ReadLock readLock2 = itemByExcelRowLock.readLock();
		ReadLock readLock3 = itemByKeyLock.readLock();
		ReadLock readLock4 = itemByRunIDLock.readLock();
		try {
			readLock.lock();
			readLock2.lock();
			readLock3.lock();
			readLock4.lock();
			return itemByConditionID.isEmpty() || itemByExcelRow.isEmpty() || itemByKey.isEmpty()
					|| itemByRunID.isEmpty();
		} finally {
			readLock.unlock();
			readLock2.unlock();
			readLock3.unlock();
			readLock4.unlock();
		}
	}

	public int sizeByKeys() {
		ReadLock readLock = itemByKeyLock.readLock();
		try {
			readLock.lock();
			return itemByKey.size();
		} finally {
			readLock.unlock();
		}
	}
}
