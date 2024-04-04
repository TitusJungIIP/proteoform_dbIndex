package edu.scripps.yates.utilities.proteomicsmodel.utils;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.proteomicsmodel.Accession;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.proteomicsmodel.AnnotationType;
import edu.scripps.yates.utilities.proteomicsmodel.HasAmounts;
import edu.scripps.yates.utilities.proteomicsmodel.HasRatios;
import edu.scripps.yates.utilities.proteomicsmodel.HasScores;
import edu.scripps.yates.utilities.proteomicsmodel.Organism;
import edu.scripps.yates.utilities.proteomicsmodel.PSM;
import edu.scripps.yates.utilities.proteomicsmodel.PTM;
import edu.scripps.yates.utilities.proteomicsmodel.PTMSite;
import edu.scripps.yates.utilities.proteomicsmodel.Peptide;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.proteomicsmodel.ProteinAnnotation;
import edu.scripps.yates.utilities.proteomicsmodel.Ratio;
import edu.scripps.yates.utilities.proteomicsmodel.Score;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AccessionType;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AmountType;
import edu.scripps.yates.utilities.proteomicsmodel.factories.OrganismEx;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.THashSet;

public class ModelUtils {
	private static final Logger log = Logger.getLogger(ModelUtils.class);
	// private static final HashMap<MSRun, Map<String, Peptide>> peptidesByMSRun
	// = new THashMap<MSRun, Map<String, Peptide>>();
	private static Organism ORGANISM_CONTAMINANT;

	public static Organism getOrganismContaminant() {
		if (ORGANISM_CONTAMINANT == null) {
			ORGANISM_CONTAMINANT = new OrganismEx("000000");
			((OrganismEx) ORGANISM_CONTAMINANT).setName(FastaParser.CONTAMINANT_PREFIX);
		}
		return ORGANISM_CONTAMINANT;
	}

	public static Set<Ratio> getProteinRatiosBetweenTwoConditions(Protein protein, String condition1Name,
			String condition2Name) {
		final Set<Ratio> ret = new THashSet<Ratio>();
		if (protein.getRatios() != null) {
			for (final Ratio ratio : protein.getRatios()) {

				if (ratio.getCondition1().getName().equals(condition1Name)) {
					if (ratio.getCondition2().getName().equals(condition2Name)) {
						ret.add(ratio);
					}
				} else if (ratio.getCondition2().getName().equals(condition1Name)) {
					if (ratio.getCondition1().getName().equals(condition2Name)) {
						ret.add(ratio);
					}

				}
			}
		}
		return ret;
	}

	/**
	 * Return a list of protein as a result of the union of the two input list of
	 * proteins
	 *
	 * @param list1
	 * @param list2
	 * @return
	 */
	public static Collection<Protein> proteinUnion(Collection<Protein> list1, Collection<Protein> list2) {

		if (list1.isEmpty() && !list2.isEmpty())
			return list2;
		if (list2.isEmpty() && !list1.isEmpty())
			return list1;
		final Set<String> proteinPrimaryAccs = new THashSet<String>();
		final List<Protein> ret = new ArrayList<Protein>();

		for (final Protein protein : list1) {
			if (!proteinPrimaryAccs.contains(protein.getPrimaryAccession().getAccession())) {
				proteinPrimaryAccs.add(protein.getPrimaryAccession().getAccession());
				ret.add(protein);
			}
		}
		for (final Protein protein : list2) {
			if (!proteinPrimaryAccs.contains(protein.getPrimaryAccession().getAccession())) {
				proteinPrimaryAccs.add(protein.getPrimaryAccession().getAccession());
				ret.add(protein);
			}
		}
		return ret;
	}

	public static Collection<PSM> psmUnion(Collection<PSM> list1, Collection<PSM> list2) {
		final Set<String> ids = new THashSet<String>();
		final List<PSM> ret = new ArrayList<PSM>();

		for (final PSM psm : list1) {
			if (!ids.contains(psm.getIdentifier())) {
				ids.add(psm.getIdentifier());
				ret.add(psm);
			}
		}
		for (final PSM psm : list2) {
			if (!ids.contains(psm.getIdentifier())) {
				ids.add(psm.getIdentifier());
				ret.add(psm);
			}
		}
		return ret;
	}

	public static Set<PSM> getPSMIntersection(Collection<PSM> set1, Collection<PSM> set2) {
		final Set<PSM> ret = new THashSet<PSM>();
		for (final PSM t1 : set1) {
			if (set2.contains(t1)) {
				ret.add(t1);
			}
		}
		return ret;
	}

	public static Set<Peptide> getPeptideIntersection(Collection<Peptide> set1, Collection<Peptide> set2) {
		final Set<Peptide> ret = new THashSet<Peptide>();
		for (final Peptide t1 : set1) {
			if (set2.contains(t1)) {
				ret.add(t1);
			}
		}
		return ret;
	}

	public static Set<Protein> getProteinIntersection(Collection<Protein> set1, Collection<Protein> set2) {
		final Set<Protein> ret = new THashSet<Protein>();
		for (final Protein p1 : set1) {
			for (final Protein p2 : set2) {
				if (p1.getPrimaryAccession().getAccession().equals(p2.getPrimaryAccession().getAccession())) {
					ret.add(p1);
					ret.add(p2);
				}
			}
		}
		return ret;
	}

	/**
	 * Gets all accessions, primary or not from the protein that are of a certain
	 * {@link AccessionType}
	 *
	 * @param prot
	 * @param accType
	 * @return
	 */
	public static List<Accession> getAccessions(Protein prot, AccessionType accType) {
		return getAccessions(prot, accType.name());
	}

	/**
	 * Gets all accessions, primeray or not from a set of proteins that are of a
	 * certain {@link AccessionType}
	 *
	 * @param prot
	 * @param accType
	 * @return
	 */
	public static List<Accession> getAccessions(Protein prot, String accType) {
		final List<Accession> ret = new ArrayList<Accession>();
		final Set<Accession> proteinAccessions = prot.getSecondaryAccessions();
		if (proteinAccessions != null) {
			for (final Accession proteinAccession : proteinAccessions) {
				if (proteinAccession.getAccessionType().name().equalsIgnoreCase(accType))
					ret.add(proteinAccession);
			}
		}
		if (prot.getPrimaryAccession().getAccessionType().name().equalsIgnoreCase(accType))
			ret.add(prot.getPrimaryAccession());
		return ret;
	}

	public static List<Protein> getProteinsFromPsms(Collection<PSM> psms) {
		final List<Protein> ret = new ArrayList<Protein>();
		final Set<String> ids = new THashSet<String>();
		for (final PSM psm : psms) {
			final Set<Protein> proteins = psm.getProteins();
			for (final Protein protein : proteins) {
				if (!ids.contains(protein.getPrimaryAccession().getAccession())) {
					ret.add(protein);
					ids.add(protein.getPrimaryAccession().getAccession());
				}
			}
		}
		return ret;
	}

	public static Set<ProteinAnnotation> getProteinAnnotations(Collection<ProteinAnnotation> proteinAnnotations,
			AnnotationType annotationType) {
		final Set<ProteinAnnotation> ret = new THashSet<ProteinAnnotation>();
		for (final ProteinAnnotation proteinAnnotation : proteinAnnotations) {
			if (proteinAnnotation.getAnnotationType().getKey().equals(annotationType.getKey())) {
				ret.add(proteinAnnotation);
			}
		}
		return ret;
	}

	public static Collection<String> getPrimaryAccessions(Collection<Protein> proteinList) {
		final Set<String> ret = new THashSet<String>();

		for (final Protein protein : proteinList) {
			ret.add(protein.getPrimaryAccession().getAccession());
		}

		return ret;
	}

	public static Collection<String> getPrimaryAccessions(Collection<Protein> proteinList, String accType) {
		final Set<String> ret = new THashSet<String>();

		for (final Protein protein : proteinList) {
			final Accession primaryAccession = protein.getPrimaryAccession();
			if (primaryAccession.getAccessionType().name().equals(accType))
				ret.add(primaryAccession.getAccession());
		}

		return ret;
	}

	public static Collection<String> getPrimaryAccessions(Collection<Protein> proteinList, AccessionType accType) {
		final Set<String> ret = new THashSet<String>();

		for (final Protein protein : proteinList) {
			final Accession primaryAccession = protein.getPrimaryAccession();
			if (primaryAccession.getAccessionType().name().equals(accType.name()))
				ret.add(primaryAccession.getAccession());
		}

		return ret;
	}

	public static Score getScore(HasScores obj, String scoreName) {

		final Set<Score> scores = obj.getScores();
		for (final Score score : scores) {
			if (score.getScoreName().equalsIgnoreCase(scoreName))
				return score;
		}

		return null;
	}

	public static List<Amount> getAmounts(HasAmounts obj, String conditionName) {
		final List<Amount> ret = new ArrayList<Amount>();

		final Set<Amount> amounts = obj.getAmounts();
		for (final Amount amount : amounts) {
			if (amount.getCondition().getName().equalsIgnoreCase(conditionName))
				ret.add(amount);
		}

		return ret;
	}

	public static List<Amount> getAmounts(HasAmounts obj, String conditionName, AmountType amountType) {
		final List<Amount> ret = new ArrayList<Amount>();
		final List<Amount> amounts = getAmounts(obj, conditionName);
		for (final Amount amount : amounts) {
			if (amount.getAmountType() == amountType)
				ret.add(amount);
		}

		return ret;
	}

	public static List<Ratio> getRatios(HasRatios obj, String condition1, String condition2) {
		final List<Ratio> ret = new ArrayList<Ratio>();

		final Set<? extends Ratio> ratios = obj.getRatios();
		for (final Ratio ratio : ratios) {
			if (ratio.getCondition1().getName().equalsIgnoreCase(condition1)) {
				if (ratio.getCondition2().getName().equalsIgnoreCase(condition2)) {
					ret.add(ratio);
				}
			} else
			// the other way round
			if (ratio.getCondition1().getName().equalsIgnoreCase(condition2)) {
				if (ratio.getCondition2().getName().equalsIgnoreCase(condition1)) {
					ret.add(ratio);
				}
			}
		}

		return ret;
	}

	public static Map<String, List<Protein>> mergeProteins(List<Protein> proteins) {
		final Map<String, List<Protein>> map = new THashMap<String, List<Protein>>();
		for (final Protein protein : proteins) {
			final String accession = protein.getPrimaryAccession().getAccession();
			if (map.containsKey(accession)) {
				map.get(accession).add(protein);
			} else {
				final List<Protein> list = new ArrayList<Protein>();
				list.add(protein);
				map.put(accession, list);
			}
		}
		return map;
	}

	public static void addToMap(Map<String, Set<Protein>> map, Protein protein) {
		final String primaryAcc = protein.getPrimaryAccession().getAccession();
		if (map.containsKey(primaryAcc)) {
			map.get(primaryAcc).add(protein);
		} else {
			final Set<Protein> set = new THashSet<Protein>();
			set.add(protein);
			map.put(primaryAcc, set);
		}
		if (protein.getSecondaryAccessions() != null) {
			for (final Accession acc : protein.getSecondaryAccessions()) {
				if (map.containsKey(acc.getAccession())) {
					map.get(acc.getAccession()).add(protein);
				} else {
					final Set<Protein> set = new THashSet<Protein>();
					set.add(protein);
					map.put(acc.getAccession(), set);
				}
			}
		}
	}

	public static void addToMap(Map<String, Set<Protein>> map, Collection<Protein> proteins) {
		for (final Protein protein : proteins) {
			addToMap(map, protein);
		}

	}

	public static void addToMap(Map<String, Set<Protein>> receiverMap, Map<String, Set<Protein>> donorMap) {
		for (final String key : donorMap.keySet()) {
			final Set<Protein> set = donorMap.get(key);
			if (receiverMap.containsKey(key)) {
				receiverMap.get(key).addAll(set);
			} else {
				final Set<Protein> set2 = new THashSet<Protein>();
				set2.addAll(set);
				receiverMap.put(key, set2);
			}
		}
	}

	public static List<Protein> getAllProteinsFromMap(Map<String, Set<Protein>> map) {
		final List<Protein> list = new ArrayList<Protein>();
		for (final Set<Protein> proteinSet : map.values()) {
			for (final Protein protein : proteinSet) {
				list.add(protein);
			}
		}
		return list;
	}

	public static Map<String, Set<PSM>> getPSMMapBySequence(Collection<PSM> psms, boolean chargeStateSensible,
			boolean ptmSensible) {
		final Map<String, Set<PSM>> ret = new THashMap<String, Set<PSM>>();
		if (psms != null) {
			for (final PSM psm : psms) {
				if (psm == null)
					continue;
				final String key = KeyUtils.getInstance().getSequenceChargeKey(psm, ptmSensible,
						chargeStateSensible);
				if (ret.containsKey(key)) {
					ret.get(key).add(psm);
				} else {
					final Set<PSM> set = new THashSet<PSM>();
					set.add(psm);
					ret.put(key, set);
				}
			}
		}
		return ret;
	}

	// /**
	// * Create the peptides of a protein taking into account the peptides per
	// * MSRun in an static way
	// *
	// * @param protein
	// */
	// public static void createPeptides(Protein protein) {
	// MSRun msRun = protein.getMSRun();
	// Set<PSM> psms = protein.getPSMs();
	// if (psms == null || psms.isEmpty()) {
	// log.warn("The protein has no psms, so, no peptides can be created");
	// return;
	// }
	// Map<String, Set<PSM>> psmMapBySequence =
	// ModelUtils.getPSMMapBySequence(psms);
	// Map<String, Peptide> peptideSet = null;
	// // Create peptides grouped by the MSRun
	//
	// if (peptidesByMSRun.containsKey(msRun)) {
	// peptideSet = peptidesByMSRun.get(msRun);
	// } else {
	// peptideSet = new THashMap<String, Peptide>();
	// peptidesByMSRun.put(msRun, peptideSet);
	// }
	// for (String sequence : psmMapBySequence.keySet()) {
	// final Set<PSM> psmsWithThatSequence = psmMapBySequence.get(sequence);
	// Peptide peptide = null;
	// if (peptideSet.containsKey(sequence)) {
	// peptide = peptideSet.get(sequence);
	// } else {
	// // create the peptide
	// peptide = new PeptideEx(sequence, msRun);
	// peptideSet.put(sequence, peptide);
	// }
	//
	// for (PSM psmWithThatSequence : psmsWithThatSequence) {
	// // psm-peptide relation
	// peptide.addPSM(psmWithThatSequence);
	// psmWithThatSequence.setPeptide(peptide);
	//
	// // protein-psm
	// protein.addPSM(psmWithThatSequence);
	// psmWithThatSequence.addProtein(protein);
	//
	// }
	// final Set<PSM> psms2 = peptide.getPSMs();
	// for (PSM psm : psms2) {
	// protein.addPSM(psm);
	// psm.addProtein(protein);
	// }
	// // protein-peptide relation
	// protein.addPeptide(peptide);
	// peptide.addProtein(protein);
	//
	// }
	//
	// }
	/**
	 * GEts a fullSequence from a cleanSequnce + PTMs<br>
	 * i.e. Having SEQUENCE and PTM at position 4 of 80Daltons:<br>
	 * the result would be SEQU(80)ENCE
	 * 
	 * @param cleanSequence
	 * @param ptms
	 * @return
	 */
	public static String getFullSequence(String cleanSequence, Collection<PTM> ptms) {
		try {
			if (ptms == null || ptms.isEmpty()) {
				return cleanSequence;
			}
			// sort ptms
			final TIntObjectHashMap<PTM> map = new TIntObjectHashMap<PTM>();
			for (final PTM ptm : ptms) {
				final List<PTMSite> ptmSites = ptm.getPTMSites();
				for (final PTMSite ptmSite : ptmSites) {
					map.put(ptmSite.getPosition(), ptm);
				}
			}
			final TIntArrayList positions = new TIntArrayList();
			positions.addAll(map.keys());
			positions.sort();
			final StringBuilder sb = new StringBuilder();
			// n-terminal modification is at position 0
			if (map.contains(0)) {
				final PTM ptm = map.get(0);
				final Double massShift = ptm.getMassShift();
				if (massShift != null) {
					sb.append("(");
					sb.append(getPtmFormatter().format(massShift));
					sb.append(")");
				}
			}
			for (int index = 0; index < cleanSequence.length(); index++) {
				if (map.contains(index + 1)) {
					sb.append(cleanSequence.charAt(index));
					final PTM ptm = map.get(index + 1);

					final Double massShift = ptm.getMassShift();
					if (massShift != null) {
						sb.append("(");
						sb.append(getPtmFormatter().format(massShift));
						sb.append(")");
					}

				} else {
					sb.append(cleanSequence.charAt(index));
				}
			}
			return sb.toString();
		} catch (final Exception e) {
			e.printStackTrace();
			log.error(e.getMessage());
			throw e;
		}
	}

	public static String getFullSequenceFromSequenceAndPTMs(String sequence, List<PTM> ptms) {
		final List<PTM> ptmList = new ArrayList<PTM>();
		// sort ptms by position
		ptmList.addAll(ptms);
		ptmList.sort(new Comparator<PTM>() {

			@Override
			public int compare(PTM o1, PTM o2) {
				return Integer.compare(sortPTMSites(o1.getPTMSites()).get(0).getPosition(),
						sortPTMSites(o2.getPTMSites()).get(0).getPosition());
			}
		});
		final StringBuilder sb = new StringBuilder();

		final List<PTMSite> sortedPTMSites = getAllPTMSitesSorted(ptms);
		int currentPosition = 1;
		try {
			for (final PTMSite ptmSite : sortedPTMSites) {
				final PTM ptm = getPTM(ptmSite, ptms);
				if (ptm == null) {
					log.error("This cannot happen");
				} else {
					int position = ptmSite.getPosition();
					if (position == sequence.length() + 1) {
						// it is c-term
						position = position - 1;
					}
					final String ptmString = "(" + getPtmFormatter().format(ptm.getMassShift()) + ")";
					sb.append(sequence.substring(currentPosition - 1, position)).append(ptmString);
					if (position < sequence.length() + 1) {
						currentPosition = position + 1;
					}
				}
			}
			sb.append(sequence.substring(currentPosition - 1));
		} catch (final Exception e) {
			e.printStackTrace();
		}
		return sb.toString();
	}

	public static int PTM_NUM_DECIMALS = 3;

	/**
	 * Call this before calling getPtmFormatter to set the number of decimals you
	 * want to consider in the PTM mass shifts
	 * 
	 * @param num
	 */
	public static void setPTM_NUM_DECIMALS(int num) {
		PTM_NUM_DECIMALS = num;
	}

	public static DecimalFormat getPtmFormatter() {

		String decimals = "";
		for (int i = 1; i <= PTM_NUM_DECIMALS; i++) {
			decimals += "#";
		}
		final DecimalFormat formatter = new DecimalFormat("+#." + decimals + ";-#." + decimals);
		return formatter;
	}

	private static Comparator<PTMSite> ptmSiteComparator;

	private static PTM getPTM(PTMSite ptmSite, List<PTM> ptMs2) {
		for (final PTM ptm : ptMs2) {
			if (ptm.getPTMSites().contains(ptmSite)) {
				return ptm;
			}
		}
		return null;
	}

	private static List<PTMSite> getAllPTMSitesSorted(Collection<PTM> ptms) {
		final List<PTMSite> ret = new ArrayList<PTMSite>();
		for (final PTM ptm : ptms) {
			ret.addAll(ptm.getPTMSites());
		}
		return sortPTMSites(ret);
	}

	private static List<PTMSite> sortPTMSites(List<PTMSite> ptmSites) {
		if (ptmSites.size() == 1) {
			return ptmSites;
		}
		if (ptmSiteComparator == null) {
			ptmSiteComparator = new Comparator<PTMSite>() {

				@Override
				public int compare(PTMSite o1, PTMSite o2) {
					return Integer.compare(o1.getPosition(), o2.getPosition());
				}
			};
		}
		ptmSites.sort(ptmSiteComparator);
		return ptmSites;
	}

}
