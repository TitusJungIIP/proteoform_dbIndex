package edu.scripps.yates.utilities.proteomicsmodel.staticstorage;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.proteomicsmodel.Accession;
import edu.scripps.yates.utilities.proteomicsmodel.Condition;
import edu.scripps.yates.utilities.proteomicsmodel.MSRun;
import edu.scripps.yates.utilities.proteomicsmodel.PSM;
import edu.scripps.yates.utilities.proteomicsmodel.Peptide;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.proteomicsmodel.utils.KeyUtils;
import edu.scripps.yates.utilities.staticstorage.ItemStorage;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

/**
 * This class is intended to storage {@link Protein}s {@link Peptide}s and
 * {@link PSM}s that are being created during import process, in a way that can
 * be retrieved by {@link Condition}, {@link MSRun} or rowIndex
 *
 * @author Salva
 *
 */
public class StaticProteomicsModelStorage {
	private static final Logger log = Logger.getLogger(StaticProteomicsModelStorage.class);
	private final static ItemStorage<Protein> proteinStorage = new ItemStorage<Protein>();
	private final static ItemStorage<Peptide> peptideStorage = new ItemStorage<Peptide>();
	private final static ItemStorage<PSM> psmStorage = new ItemStorage<PSM>();
	private final static Map<String, MSRun> msRunsByID = new THashMap<String, MSRun>();

	public static void clearData() {
		if (!proteinStorage.isEmpty()) {
			log.info("Clearing static protein storage with " + proteinStorage.sizeByKeys());
		}
		proteinStorage.clearData();
		if (!peptideStorage.isEmpty()) {
			log.info("Clearing static peptide storage with " + peptideStorage.sizeByKeys());
		}
		peptideStorage.clearData();
		if (!psmStorage.isEmpty()) {
			log.info("Clearing static psm storage with " + psmStorage.sizeByKeys());
		}
		psmStorage.clearData();
		msRunsByID.clear();
	}

	public static void addMSRun(MSRun msRun) {
		msRunsByID.put(msRun.getRunId(), msRun);
	}

	public static boolean containsMSRun(String msRunID) {
		return msRunsByID.containsKey(msRunID);
	}

	public static MSRun getMSRun(String msRunID) {
		return msRunsByID.get(msRunID);
	}

	public static void addProtein(Protein protein, String msRunID, String conditionID) {
		addProtein(protein, msRunID, conditionID, -1);
	}

	public static void addProtein(Protein protein, MSRun msRun, String conditionID) {
		String msRunID = null;
		if (msRun != null) {
			msRunID = msRun.getRunId();
		}
		addProtein(protein, msRunID, conditionID);
	}

	public static void addProtein(Protein protein, Collection<MSRun> msRuns, String conditionID) {
		for (final MSRun msRun2 : msRuns) {
			addProtein(protein, msRun2, conditionID);
		}
	}

	public static void addProtein(Protein protein, Collection<MSRun> msRuns, String conditionID, int excelRowIndex) {
		for (final MSRun msRun : msRuns) {
			addProtein(protein, msRun.getRunId(), conditionID, excelRowIndex);
		}
	}

	public static void addProtein(Protein protein, String msRunID, String conditionID, int excelRowIndex) {
		proteinStorage.add(protein, msRunID, conditionID, excelRowIndex, protein.getAccession());
		if (protein.getSecondaryAccessions() != null) {
			for (final Accession secondaryAccession : protein.getSecondaryAccessions()) {
				proteinStorage.add(protein, msRunID, conditionID, excelRowIndex, secondaryAccession.getAccession());
			}
		}
	}

	public static void addPeptide(Peptide peptide, String msRunID, String conditionID,
			boolean distinguishModifiedPeptides, boolean chargeStateSensible) {
		addPeptide(peptide, msRunID, conditionID, -1, distinguishModifiedPeptides, chargeStateSensible);
	}

	public static void addPeptide(Peptide peptide, String msRunID, String conditionID, String peptideKey) {
		addPeptide(peptide, msRunID, conditionID, -1, peptideKey);
	}

	public static void addPeptide(Peptide peptide, MSRun msRun, String conditionID, boolean distinguishModifiedPeptides,
			boolean chargeStateSensible) {
		addPeptide(peptide, msRun.getRunId(), conditionID, distinguishModifiedPeptides, chargeStateSensible);
	}

	public static void addPeptide(Peptide peptide, MSRun msRun, String conditionID, String peptideKey) {
		if (msRun != null) {
			addPeptide(peptide, msRun.getRunId(), conditionID, peptideKey);
		} else {
			final String runID = null;
			addPeptide(peptide, runID, conditionID, peptideKey);
		}
	}

	public static void addPeptide(Peptide peptide, Collection<MSRun> msRuns, String conditionID,
			boolean distinguishModifiedPeptides, boolean chargeStateSensible) {
		addPeptide(peptide, msRuns, conditionID, -1, distinguishModifiedPeptides, chargeStateSensible);
	}

	public static void addPeptide(Peptide peptide, Collection<MSRun> msRuns, String conditionID, String peptideKey) {
		addPeptide(peptide, msRuns, conditionID, -1, peptideKey);
	}

	public static void addPeptide(Peptide peptide, Collection<MSRun> msRuns, String conditionID, int excelRowIndex,
			boolean distinguishModifiedPeptides, boolean chargeStateSensible) {
		for (final MSRun msRun2 : msRuns) {
			addPeptide(peptide, msRun2.getRunId(), conditionID, excelRowIndex, distinguishModifiedPeptides,
					chargeStateSensible);
		}
	}

	public static void addPeptide(Peptide peptide, Collection<MSRun> msRuns, String conditionID, int excelRowIndex,
			String peptideKey) {
		for (final MSRun msRun2 : msRuns) {
			addPeptide(peptide, msRun2.getRunId(), conditionID, excelRowIndex, peptideKey);
		}
	}

	public static void addPeptide(Peptide peptide, String msRunID, String conditionID, int excelRowIndex,
			String peptideKey) {
		if (conditionID == null) {
			// log.info("condition is null for peptide");
		}

		peptideStorage.add(peptide, msRunID, conditionID, excelRowIndex, peptideKey);
	}

	public static void addPeptide(Peptide peptide, String msRunID, String conditionID, int excelRowIndex,
			boolean distinguishModifiedPeptides, boolean chargeStateSensible) {
		if (conditionID == null) {
			// log.info("condition is null for peptide");
		}
		final String key = KeyUtils.getInstance().getSequenceChargeKey(peptide, distinguishModifiedPeptides,
				chargeStateSensible);
		peptideStorage.add(peptide, msRunID, conditionID, excelRowIndex, key);
	}

	public static void addPSM(PSM psm, String runID, String conditionID) {
		addPSM(psm, runID, conditionID, -1);
	}

	public static void addPSM(PSM psm, MSRun msRun, String conditionID) {
		addPSM(psm, msRun.getRunId(), conditionID);
	}

	public static void addPSM(PSM psm, String runID, String conditionID, int excelRowIndex) {
		if (conditionID == null) {
			// log.info("condition is null for psm");
		}
		psmStorage.add(psm, runID, conditionID, excelRowIndex, psm.getIdentifier());
	}

	public static boolean containsProtein(Collection<MSRun> msRuns, String conditionID, String accession) {
		for (final MSRun msRun : msRuns) {
			final boolean b = containsProtein(msRun, conditionID, accession);
			if (b) {
				return b;
			}
		}
		return false;
	}

	public static boolean containsProtein(MSRun msRun, String conditionID, String accession) {
		String msRunID = null;
		if (msRun != null) {
			msRunID = msRun.getRunId();
		}
		return containsProtein(msRunID, conditionID, -1, accession);
	}

	public static boolean containsProtein(String msRunID, String conditionID, String accession) {
		return containsProtein(msRunID, conditionID, -1, accession);
	}

	public static boolean containsProtein(String msRunID, String conditionID, int excelRowIndex, String accession) {
		return proteinStorage.contains(msRunID, conditionID, excelRowIndex, accession);
	}

	public static boolean containsPSM(String msRun, String conditionID, String psmID) {
		return containsPSM(msRun, conditionID, -1, psmID);
	}

	public static boolean containsPSM(MSRun msRun, String conditionID, String psmID) {
		return containsPSM(msRun.getRunId(), conditionID, psmID);
	}

	public static boolean containsPSM(String msRun, String conditionID, int excelRowIndex, String psmID) {
		String psmId2 = psmID;
		if (psmId2 == null) {
			psmId2 = (excelRowIndex + 1) + "-" + msRun;
		}
		return psmStorage.contains(msRun, conditionID, excelRowIndex, psmId2);
	}

	public static boolean containsPeptide(String msRunID, String conditionID, String key) {
		return containsPeptide(msRunID, conditionID, -1, key);
	}

	public static boolean containsPeptide(MSRun msRun, String conditionID, String key) {
		if (msRun != null) {
			return containsPeptide(msRun.getRunId(), conditionID, key);
		} else {
			final String runID = null;
			return containsPeptide(runID, conditionID, key);
		}
	}

	public static boolean containsPeptide(Collection<MSRun> msRuns, String conditionID, String key) {
		for (final MSRun msRun2 : msRuns) {
			final boolean b = containsPeptide(msRun2, conditionID, key);
			if (b) {
				return b;
			}
		}
		return false;
	}

	public static boolean containsPeptide(String msRunID, String conditionID, int excelRowIndex, String key) {
		return peptideStorage.contains(msRunID, conditionID, excelRowIndex, key);
	}

	public static boolean containsPeptide(MSRun msRun, String conditionID, int excelRowIndex, String key) {
		return containsPeptide(msRun.getRunId(), conditionID, excelRowIndex, key);
	}

	public static Set<Protein> getProtein(String conditionID, String accession, Collection<MSRun> msRuns) {
		final Set<Protein> ret = new THashSet<Protein>();
		for (final MSRun msRun : msRuns) {
			ret.addAll(getProtein(msRun, conditionID, accession));
		}
		return ret;
	}

	public static Protein getSingleProtein(String conditionID, String accession, Collection<MSRun> msRuns) {
		final Set<Protein> proteins = getProtein(conditionID, accession, msRuns);
		if (proteins == null || proteins.isEmpty()) {
			return null;
		}
		if (proteins.size() > 1) {
			log.warn("Retrieved proteins are multiple!");
		}
		return proteins.iterator().next();
	}

	public static Set<Protein> getProtein(MSRun msRun, String conditionID, String accession) {
		String msRunID = null;
		if (msRun != null) {
			msRunID = msRun.getRunId();
		}
		return getProtein(msRunID, conditionID, accession);
	}

	public static Set<Protein> getProtein(String msRunID, String conditionID, String accession) {
		return getProtein(msRunID, conditionID, -1, accession);
	}

	public static Set<Protein> getProtein(Collection<String> msRunIDs, String conditionID, String accession) {
		final Set<Protein> set = new THashSet<Protein>();
		for (final String msRunID : msRunIDs) {
			set.addAll(getProtein(msRunID, conditionID, -1, accession));
		}
		return set;
	}

	public static Set<Protein> getProtein(Collection<String> msRunIDs, String conditionID, int excelRowIndex,
			String accession) {
		final Set<Protein> set = new THashSet<Protein>();
		for (final String msRunID : msRunIDs) {
			set.addAll(proteinStorage.get(msRunID, conditionID, excelRowIndex, accession));
		}
		return set;
	}

	public static Set<Protein> getProtein(String msRunID, String conditionID, int excelRowIndex, String accession) {
		return proteinStorage.get(msRunID, conditionID, excelRowIndex, accession);
	}

	public static Set<PSM> getPSM(String runID, String conditionID, String psmID) {
		return getPSM(runID, conditionID, -1, psmID);
	}

	public static Set<PSM> getPSM(MSRun msRun, String conditionID, String psmID) {
		return getPSM(msRun.getRunId(), conditionID, psmID);
	}

	public static PSM getSinglePSM(MSRun msRun, String conditionID, String psmID) {
		final Set<PSM> psms = getPSM(msRun, conditionID, psmID);
		if (psms == null || psms.isEmpty()) {
			return null;
		}
		if (psms.size() > 1) {
			log.warn("Retrieved psms are multiple!");
		}
		return psms.iterator().next();
	}

	public static Set<PSM> getPSM(String runID, String conditionID, int excelRowIndex, String psmID) {
		String psmId2 = psmID;
		if (psmId2 == null) {
			psmId2 = (excelRowIndex + 1) + "-" + runID;
		}
		return psmStorage.get(runID, conditionID, excelRowIndex, psmId2);
	}

	public static Set<Peptide> getPeptide(String msRunID, String conditionID, String sequence) {
		return getPeptide(msRunID, conditionID, -1, sequence);
	}

	public static Set<Peptide> getPeptide(String conditionID, String key, Collection<MSRun> msRuns) {
		final Set<Peptide> ret = new THashSet<Peptide>();
		for (final MSRun msRun : msRuns) {
			final Set<Peptide> ret2 = getPeptide(msRun, conditionID, key);
			if (ret2 != null) {
				ret.addAll(ret2);
			}
		}
		return ret;
	}

	public static Set<Peptide> getPeptide(MSRun msRun, String conditionID, String key) {
		if (msRun != null) {
			return getPeptide(msRun.getRunId(), conditionID, key);
		} else {
			final String runID = null;
			return getPeptide(runID, conditionID, key);
		}
	}

	public static Peptide getSinglePeptide(Collection<MSRun> msRuns, String conditionID, String key) {
		final Set<Peptide> peptides = getPeptide(conditionID, key, msRuns);
		if (peptides == null || peptides.isEmpty()) {
			return null;
		}
		if (peptides.size() > 1) {
			log.warn("Retrieved peptides are multiple!");
		}
		return peptides.iterator().next();
	}

	public static Peptide getSinglePeptide(String msRunID, String conditionID, String sequence) {
		final Set<Peptide> peptides = getPeptide(msRunID, conditionID, sequence);
		if (peptides == null || peptides.isEmpty()) {
			return null;
		}
		if (peptides.size() > 1) {
			log.warn("Retrieved peptides are multiple!");
		}
		return peptides.iterator().next();
	}

	public static Peptide getSinglePeptide(MSRun msRun, String conditionID, String sequence) {
		final Set<Peptide> peptides = getPeptide(msRun, conditionID, sequence);
		if (peptides == null || peptides.isEmpty()) {
			return null;
		}
		if (peptides.size() > 1) {
			log.warn("Retrieved peptides are multiple!");
		}
		return peptides.iterator().next();
	}

	public static Protein getSingleProtein(String msRunID, String conditionID, String accession) {
		final Set<Protein> proteins = getProtein(msRunID, conditionID, accession);
		if (proteins == null || proteins.isEmpty()) {
			return null;
		}
		if (proteins.size() > 1) {
			log.warn("Retrieved proteins are multiple!");
		}
		return proteins.iterator().next();
	}

	public static Set<Peptide> getPeptide(String msRunID, String conditionID, int excelRowIndex, String sequence) {
		return peptideStorage.get(msRunID, conditionID, excelRowIndex, sequence);
	}

	public static Set<Peptide> getPeptide(Collection<String> msRunIDs, String conditionID, String sequence) {
		final Set<Peptide> set = new THashSet<Peptide>();
		for (final String msRunID : msRunIDs) {
			set.addAll(getPeptide(msRunID, conditionID, -1, sequence));
		}
		return set;
	}

	public static Set<Peptide> getPeptide(Collection<String> msRunIDs, String conditionID, int excelRowIndex,
			String sequence) {
		final Set<Peptide> set = new THashSet<Peptide>();
		for (final String msRunID : msRunIDs) {
			set.addAll(peptideStorage.get(msRunID, conditionID, excelRowIndex, sequence));
		}
		return set;

	}

	public static Set<Peptide> getPeptide(String conditionID, int excelRowIndex, String sequence,
			Collection<MSRun> msRuns) {
		final Set<Peptide> set = new THashSet<Peptide>();
		for (final MSRun msRun : msRuns) {
			set.addAll(peptideStorage.get(msRun.getRunId(), conditionID, excelRowIndex, sequence));
		}
		return set;

	}

	public static void addProtein(Protein protein, List<String> msRunIDs, String conditionID) {
		for (final String msRunID : msRunIDs) {
			addProtein(protein, msRunID, conditionID);
		}
	}

	public static void addPeptide(Peptide peptide, List<String> msRunIDs, String conditionID, String peptideKey) {
		for (final String msRunID : msRunIDs) {
			addPeptide(peptide, msRunID, conditionID, -1, peptideKey);
		}
	}

	public static void addPeptide(Peptide peptide, List<String> msRunIDs, String conditionID,
			boolean distinguishModifiedPeptides, boolean chargeStateSensible) {
		for (final String msRunID : msRunIDs) {
			addPeptide(peptide, msRunID, conditionID, -1, distinguishModifiedPeptides, chargeStateSensible);
		}
	}
}
