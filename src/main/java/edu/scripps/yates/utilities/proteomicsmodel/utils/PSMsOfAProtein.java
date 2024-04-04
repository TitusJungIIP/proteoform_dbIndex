package edu.scripps.yates.utilities.proteomicsmodel.utils;

import java.util.Collection;
import java.util.Set;

import edu.scripps.yates.utilities.collections.ObservableArrayList;
import edu.scripps.yates.utilities.proteomicsmodel.PSM;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import gnu.trove.set.hash.THashSet;

public class PSMsOfAProtein extends ObservableArrayList<PSM> {
	/**
	 * 
	 */
	private static final long serialVersionUID = 5984394901348424280L;

	private final Set<PSM> set = new THashSet<PSM>();

	public PSMsOfAProtein(Protein protein) {
		super(new PSMsOfAProteinCollectionObserver(protein));
	}

	@Override
	public boolean add(PSM obj) {
		set.add(obj);
		return super.add(obj);
	}

	@Override
	public void clear() {
		set.clear();
		super.clear();
	}

	@Override
	public boolean remove(Object obj) {
		set.remove(obj);
		return super.remove(obj);
	}

	@Override
	public boolean addAll(Collection<? extends PSM> collection) {
		set.addAll(collection);
		return super.addAll(collection);
	}

	@Override
	public boolean removeAll(Collection<?> collection) {
		set.removeAll(collection);
		return super.removeAll(collection);
	}

	@Override
	public boolean contains(Object o) {
		return set.contains(o);
	}

	@Override
	public boolean containsAll(Collection<?> c) {
		return set.containsAll(c);
	}
}
