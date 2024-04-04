package edu.scripps.yates.utilities.proteomicsmodel.utils;

import java.util.Collection;

import edu.scripps.yates.utilities.collections.CollectionObserver;
import edu.scripps.yates.utilities.proteomicsmodel.PSM;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;

/**
 * This class ensures that every time a PSM is added or deleted to a protein,
 * the corresponding peptide is added
 * 
 * @author salvador
 *
 */
public class PSMsOfAProteinCollectionObserver extends CollectionObserver<PSM> {
	private final Protein protein;

	public PSMsOfAProteinCollectionObserver(Protein protein) {
		this.protein = protein;
	}

	@Override
	public boolean add(PSM psm) {
		final boolean ret = protein.addPeptide(psm.getPeptide(), false);
		protein.addMSRun(psm.getMSRun());
		return ret;
	}

	@Override
	public void clear() {
		protein.getPeptides().clear();
	}

	@Override
	public boolean remove(Object obj) {
		if (obj instanceof PSM) {
			final PSM psm = (PSM) obj;
			// only remove peptide from protein if that peptide doesn't have
			// more psms
			psm.getPeptide().getPSMs().remove(psm);
			if (psm.getPeptide().getPSMs().isEmpty()) {
				return protein.getPeptides().remove(psm.getPeptide());
			}
		}
		return false;
	}

	@Override
	public boolean addAll(Collection<? extends PSM> collection) {
		boolean ret = false;
		for (final PSM psm : collection) {
			final boolean b = protein.addPeptide(psm.getPeptide(), false);
			if (b) {
				ret = true;
			}
		}
		return ret;
	}

	@Override
	public boolean removeAll(Collection<?> collection) {
		boolean ret = true;
		for (final Object obj : collection) {
			if (obj instanceof PSM) {
				final PSM psm = (PSM) obj;
				psm.getPeptide().getPSMs().remove(psm);
				if (psm.getPeptide().getPSMs().isEmpty()) {
					final boolean b = protein.getPeptides().remove(psm.getPeptide());
					if (!b) {
						ret = false;
					}
				}
			}
		}
		return ret;
	}

}
