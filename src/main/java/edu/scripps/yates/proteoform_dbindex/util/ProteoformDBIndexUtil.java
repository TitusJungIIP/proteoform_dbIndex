package edu.scripps.yates.proteoform_dbindex.util;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.log4j.Logger;
import org.proteored.miapeapi.psimod.PSIModOBOPlainTextReader;

import com.compomics.util.general.UnknownElementMassException;

import edu.scripps.yates.annotations.uniprot.UniprotCVTermCode;
import edu.scripps.yates.annotations.uniprot.UniprotPTMCVReader;
import edu.scripps.yates.annotations.uniprot.UniprotPTMCVTerm;
import edu.scripps.yates.annotations.uniprot.proteoform.Proteoform;
import edu.scripps.yates.annotations.uniprot.proteoform.ProteoformType;
import edu.scripps.yates.annotations.uniprot.proteoform.UniprotPTM;
import edu.scripps.yates.dbindex.Constants;
import edu.scripps.yates.proteoform_dbindex.ProteoformDBIndexer;
import edu.scripps.yates.proteoform_dbindex.model.ExtendedAssignMass;
import edu.scripps.yates.proteoform_dbindex.model.PhosphositeDB;
import edu.scripps.yates.proteoform_dbindex.model.SequenceChange;
import edu.scripps.yates.proteoform_dbindex.model.SequenceWithModification;
import edu.scripps.yates.utilities.masses.FormulaCalculator;
import edu.scripps.yates.utilities.sequence.MyEnzyme;
import edu.scripps.yates.utilities.strings.StringUtils;
import gnu.trove.iterator.TDoubleIterator;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.THashSet;

public class ProteoformDBIndexUtil {
	private final static Logger log = Logger.getLogger(ProteoformDBIndexUtil.class);
	private final static String MISSING_AT = "_missing_at_";
	private static ProteoformDBIndexUtil instance;
	private final static FormulaCalculator calc = new FormulaCalculator();
	private final Set<String> ignoredPTMs = new HashSet<String>();
	private final PSIModOBOPlainTextReader psiModReader;
	private final Set<String> ptms = new HashSet<String>();
	private static int staticCallDepp = 0;
	private static final String AT = "at_";

	private ProteoformDBIndexUtil() {
		try {
			psiModReader = PSIModOBOPlainTextReader.getInstance();
		} catch (final IOException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(e);
		}
	}

	public static ProteoformDBIndexUtil getInstance() {
		if (instance == null) {
			instance = new ProteoformDBIndexUtil();
		}
		return instance;
	}

	private String getProteoformSequenceOriginalChange(Proteoform proteoform) {
		if (proteoform.isOriginal()) {
			return null;
		}
		// _D->N_
		final String id = proteoform.getId();
		if (id.contains("->")) {
			final String[] split = id.split("->");
			final String ret = split[0].substring(split[0].lastIndexOf("_") + 1);
			return ret;
		} else {

			if (id.contains(MISSING_AT)) {
				String missing = id.substring(id.indexOf(MISSING_AT) + MISSING_AT.length());
				if (missing.startsWith("[")) {
					missing = missing.substring(1);
				}
				if (missing.endsWith("]")) {
					missing = missing.substring(0, missing.length() - 1);
				}
				if (missing.contains("-")) {
					final int start = Integer.valueOf(missing.split("-")[0]);
					final int end = Integer.valueOf(missing.split("-")[1]);
					return proteoform.getOriginalSeq().substring(start - 1, end - 1);
				} else {
					return String.valueOf(proteoform.getOriginalSeq().charAt(Integer.valueOf(missing) - 1));
				}

			}
		}
		throw new IllegalArgumentException();
	}

	private String getProteoformSequenceChange(Proteoform proteoform) {
		if (proteoform.isOriginal()) {
			return null;
		}
		// _D->N_
		final String id = proteoform.getId();
		if (id.contains("->")) {
			final String[] split = id.split("->");
			final String ret = split[1].substring(0, split[1].indexOf("_"));
			return ret;
		} else if (id.contains("_missing_at")) {
			return "";
		}
		throw new IllegalArgumentException();
	}

	private TDoubleArrayList getProteoformPTMMassChanges(int proteoformPositionInProtein, Proteoform proteoform)
			throws UnknownElementMassException {
		if (proteoform.getPtms() != null && !proteoform.getPtms().isEmpty()) {
			final TDoubleArrayList ret = new TDoubleArrayList();
			for (final UniprotPTM uniprotPTM : proteoform.getPtms()) {
				if (uniprotPTM.getPositionInProtein() == proteoformPositionInProtein) {
					final Double massChange = getMassChange(uniprotPTM);
					if (massChange != null) {
						ret.add(massChange);
					} else {
						ignoredPTMs.add(uniprotPTM.getSingleValue(UniprotCVTermCode.AC) + "\t"
								+ uniprotPTM.getSingleValue(UniprotCVTermCode.ID));
					}
				}
			}
			return ret;
		}
		return null;
	}

	private Double getMassChange(UniprotPTM uniprotPTM) throws UnknownElementMassException {

		String formula = uniprotPTM.getSingleValue(UniprotCVTermCode.CF);
		if (formula != null) {
			formula = formula.replace(" ", "");
			return calc.calculateMass(formula);
		} else {
			final Set<String> valueSet = uniprotPTM.getValueSet(UniprotCVTermCode.DR);
			if (valueSet != null) {
				for (final String string2 : valueSet) {
					if (string2.startsWith("PSI-MOD")) {
						// as PSI-MOD; MOD:00214.
						final String psiModID = string2.split(" ")[1].replace(".", "");

						formula = psiModReader.getTermXRefDiffFormula(psiModID);
						if (formula != null) {
							if ("none".equals(formula)) {
								return null;
							}
							try {
								return calc.calculateMass(formula);
							} catch (final UnknownElementMassException e) {
								e.printStackTrace();
							}
						}
					}
				}
			}
		}
		return null;
	}

	public Map<String, List<Proteoform>> loadProteoformMapFromPhosphoSite(PhosphositeDB phosphositeDB,
			String accession) {
		final Map<String, List<Proteoform>> ret = new THashMap<String, List<Proteoform>>();
		final TIntArrayList phosphorilatedPositions = phosphositeDB.getPhosphorilatedPositions(accession);
		if (phosphorilatedPositions != null) {
			final List<Proteoform> list = new ArrayList<Proteoform>();
			ret.put(accession, list);
			final String proteinSeq = phosphositeDB.getProteinSeq(accession);

			for (final int position : phosphorilatedPositions.toArray()) {
				final char aa = proteinSeq.charAt(position - 1);
				final Proteoform proteoform = new Proteoform(accession, proteinSeq, accession, proteinSeq, "",
						"Phospho(" + aa + ")_at_" + position, null, null, null, ProteoformType.PTM, true, true);

				final UniprotPTMCVTerm uniprotPTMCVTerm = UniprotPTMCVReader.getInstance()
						.getPtmsByID(getPhosphoSitePTMID(aa));

				final UniprotPTM uniprotPTM = new UniprotPTM(position, uniprotPTMCVTerm);
				proteoform.addPTM(uniprotPTM);
				list.add(proteoform);
			}

		}
		return ret;
	}

	private String getPhosphoSitePTMID(char aa) {
		if (aa == 'y') {
			return "Phosphotyrosine";
		} else if (aa == 's') {
			return "Phosphoserine";
		} else if (aa == 't') {
			return "Phosphothreonine";

		}
		return null;
	}

	public TIntObjectHashMap<List<Proteoform>> getProteoformsByPositionInProtein(List<Proteoform> proteoforms) {
		final TIntObjectHashMap<List<Proteoform>> proteoformsByPositionInProtein = new TIntObjectHashMap<List<Proteoform>>();
		for (final Proteoform proteoform : proteoforms) {
			// ignore mutagenesis, because they are artificial mutations
			if (proteoform.getProteoformType() == ProteoformType.MUTAGENESIS_SITE) {
				continue;
			}
			// deal with isoforms
			if (proteoform.getProteoformType() == ProteoformType.ISOFORM) {
				System.out.println(proteoform.getId());
			}

			final String id = proteoform.getId();
			if (id.contains(AT)) {
				final String substring = id.substring(id.indexOf(AT) + 3, id.length());
				int position;
				if (substring.startsWith("[")) {
					position = Integer.valueOf(substring.split("-")[0].substring(1));
				} else {
					position = Integer.valueOf(substring);
				}
				if (proteoformsByPositionInProtein.containsKey(position)) {
					proteoformsByPositionInProtein.get(position).add(proteoform);
				} else {
					final List<Proteoform> list = new ArrayList<Proteoform>();
					list.add(proteoform);
					proteoformsByPositionInProtein.put(position, list);
				}
			} else {
				if (proteoform.getPtms() != null && !proteoform.getPtms().isEmpty()) {
					for (final UniprotPTM uniprotPTM : proteoform.getPtms()) {
						final int position = uniprotPTM.getPositionInProtein();
						if (proteoformsByPositionInProtein.containsKey(position)) {
							proteoformsByPositionInProtein.get(position).add(proteoform);
						} else {
							final List<Proteoform> list = new ArrayList<Proteoform>();
							list.add(proteoform);
							proteoformsByPositionInProtein.put(position, list);
						}
					}
				}
			}
		}

		return proteoformsByPositionInProtein;
	}

	/**
	 * 
	 * @param peptideSeq
	 * @param peptideInitInProtein
	 * @param nonIsoformsProteoformsByPositionInMainProtein
	 * @param isoform                                       can be null
	 * @return
	 * @throws UnknownElementMassException
	 */
	public List<SequenceChange>  getSequenceChangesInPeptide(String peptideSeq, int peptideInitInProtein,
			TIntObjectHashMap<List<Proteoform>> nonIsoformsProteoformsByPositionInMainProtein, Proteoform isoform)
			throws UnknownElementMassException {
		if (isoform == null) {
			return getSequenceChangesInPeptideTowardsMainProteinform(peptideSeq, peptideInitInProtein,
					nonIsoformsProteoformsByPositionInMainProtein);
		} else {
			return getSequenceChangesInPeptideTowardsIsoform(peptideSeq, nonIsoformsProteoformsByPositionInMainProtein,
					isoform);
		}
	}

	/**
	 * 
	 * @param peptideSeq
	 * @param peptideInitInOriginalProtein
	 * @param proteoformsByPositionInProtein
	 * @param isoform                        can be null
	 * @return
	 * @throws UnknownElementMassException
	 */
	private List<SequenceChange> getSequenceChangesInPeptideTowardsMainProteinform(String peptideSeq,
			int peptideInitInOriginalProtein, TIntObjectHashMap<List<Proteoform>> proteoformsByPositionInProtein)
			throws UnknownElementMassException {
		final List<SequenceChange> ret = new ArrayList<SequenceChange>();

		final List<Integer> positionsInOrder = new ArrayList<Integer>();
		proteoformsByPositionInProtein.forEachKey(k -> positionsInOrder.add(k));

		Collections.sort(positionsInOrder);
		final int peptideLength = peptideSeq.length();
		final int peptideEndInProtein = peptideInitInOriginalProtein + peptideLength - 1;
		for (final Integer proteoformPositionInProtein : positionsInOrder) {
			// if the position of the sequence variance is after the end of the peptide in
			// the protein, skip it
			if (peptideEndInProtein < proteoformPositionInProtein) {
				break;
			}
			// if the position of the sequence variance goes into the peptide
			if (proteoformPositionInProtein >= peptideInitInOriginalProtein
					&& proteoformPositionInProtein <= peptideEndInProtein) {
				final List<Proteoform> proteoforms = proteoformsByPositionInProtein.get(proteoformPositionInProtein);
				for (final Proteoform proteoform : proteoforms) {

					final int positionInPeptide = proteoformPositionInProtein - peptideInitInOriginalProtein + 1;
					final String original = getProteoformSequenceOriginalChange(proteoform);
					if (original != null && proteoformPositionInProtein + original.length() - 1 > peptideEndInProtein) {
						// the change is out of the peptide
						continue;
					}

					if (original != null && !peptideSeq.contains(original)) {
						// the change is out of the peptide
						continue;
					}
					final TDoubleArrayList massChanges = getProteoformPTMMassChanges(proteoformPositionInProtein,
							proteoform);
					final String change = getProteoformSequenceChange(proteoform);

					if (massChanges != null) {
						final TDoubleIterator iterator = massChanges.iterator();
						while (iterator.hasNext()) {
							final Double massChange = iterator.next();
							final SequenceChange seqChange = new SequenceChange(proteoform.getId(), ProteoformType.PTM,
									positionInPeptide, original, change, massChange);
							ret.add(seqChange);
						}
					} else {
						final SequenceChange seqChange = new SequenceChange(proteoform.getId(),
								proteoform.getProteoformType(), positionInPeptide, original, change, null);
						ret.add(seqChange);
					}
				}
			}

		}
		// because implements comparable
		Collections.sort(ret);
		return ret;

	}

	/**
	 * 
	 * @param isoformPeptideSeq
	 * @param proteoformsByPositionInMainProtein
	 * @param isoform                            is not null
	 * 
	 * @return
	 * @throws UnknownElementMassException
	 */
	private List<SequenceChange> getSequenceChangesInPeptideTowardsIsoform(String isoformPeptideSeq,
			TIntObjectHashMap<List<Proteoform>> proteoformsByPositionInMainProtein, Proteoform isoform)
			throws UnknownElementMassException {
		final List<SequenceChange> ret = new ArrayList<SequenceChange>();

		final List<Integer> positionsInOrder = new ArrayList<Integer>();
		proteoformsByPositionInMainProtein.forEachKey(k -> positionsInOrder.add(k));

		Collections.sort(positionsInOrder);
		for (final Integer proteoformPositionInMainProtein : positionsInOrder) {
			final List<Proteoform> proteoformsToApply = proteoformsByPositionInMainProtein
					.get(proteoformPositionInMainProtein);
			final int isoformPeptideLength = isoformPeptideSeq.length();
			for (final Proteoform proteoformToApply : proteoformsToApply) {
				final String mainProteinSeq = proteoformToApply.getOriginalSeq();
				final int peptideInitInMainProtein = mainProteinSeq.indexOf(isoformPeptideSeq) + 1;
				final int peptideEndInMainProtein = peptideInitInMainProtein + isoformPeptideLength - 1;
				if (peptideInitInMainProtein > 0) {// otherwise it cannot be
													// mapped
					if (peptideEndInMainProtein < proteoformPositionInMainProtein) {
						break;
					}
					if (proteoformPositionInMainProtein >= peptideInitInMainProtein
							&& proteoformPositionInMainProtein <= peptideEndInMainProtein) {

						final int positionInPeptide = proteoformPositionInMainProtein - peptideInitInMainProtein + 1;
						final String original = getProteoformSequenceOriginalChange(proteoformToApply);
						if (original != null
								&& proteoformPositionInMainProtein + original.length() - 1 > peptideEndInMainProtein) {
							// the change is out of the peptide
							continue;
						}
						if (original != null && !isoformPeptideSeq.contains(original)) {
							// the change is out of the peptide
							continue;
						}
						final TDoubleArrayList massChanges = getProteoformPTMMassChanges(
								proteoformPositionInMainProtein, proteoformToApply);

						final String change = getProteoformSequenceChange(proteoformToApply);

						if (massChanges != null) {
							final TDoubleIterator iterator = massChanges.iterator();
							while (iterator.hasNext()) {
								final Double massChange = iterator.next();
								final SequenceChange seqChange = new SequenceChange(isoform.getId(), ProteoformType.PTM,
										positionInPeptide, original, change, massChange);
								ret.add(seqChange);
							}
						} else {
							// if we have the proteoform id as
							// P12345_missint_at_23,
							// change to P12345_missint_at_23_mappedto_P12345-2
							final String id = proteoformToApply.getId() + "_mappedto_" + isoform.getId();
							final SequenceChange seqChange = new SequenceChange(id,
									proteoformToApply.getProteoformType(), positionInPeptide, original, change, null);
							ret.add(seqChange);
						}
					}
				}
			}
		}

		return ret;

	}

	public List<SequenceWithModification> getAllCombinationsForPeptide(String peptideSeq, String proteinSequence,
			int positionOfPeptideInProtein, MyEnzyme enzyme, List<SequenceChange> sequenceChanges,
			int maxNumVariationsPerPeptide, ExtendedAssignMass extendedAssignMass, Set<String> peptidesProcessed)
			throws IOException {
		staticCallDepp++;
		if (staticCallDepp > 4) {
//			log.info("deep call " + staticCallDepp);
		}
		peptidesProcessed.add(peptideSeq);
		final List<SequenceWithModification> ret = new ArrayList<SequenceWithModification>();
		// check whether the peptide has the correct number of
		// allowed misscleavages
		if (getNumMissedClavages(peptideSeq, enzyme) <= enzyme.getMiscleavages()
				&& peptideSeq.length() >= Constants.MIN_PEP_LENGTH) {
			ret.add(new SequenceWithModification(peptideSeq, peptideSeq, extendedAssignMass, proteinSequence));
		}
		final List<SequenceChange> sequenceChangesNaturalVariants = sequenceChanges.stream()
				.filter(s -> s.getProteoformType() == ProteoformType.NATURAL_VARIANT).collect(Collectors.toList());
		final List<SequenceChange> sequenceChangesOthers = sequenceChanges.stream()
				.filter(s -> s.getProteoformType() != ProteoformType.NATURAL_VARIANT).collect(Collectors.toList());

		for (int i = 0; i <= sequenceChangesNaturalVariants.size(); i++) {
			final List<SequenceChange> listOfSequenceChangesToApply = new ArrayList<SequenceChange>();
			listOfSequenceChangesToApply.addAll(sequenceChangesOthers);
			if (i < sequenceChangesNaturalVariants.size()) {
				listOfSequenceChangesToApply.add(sequenceChangesNaturalVariants.get(i));
			}

			final int maxNumberOfVariationsInPeptide = Math.min(maxNumVariationsPerPeptide,
					Math.min(listOfSequenceChangesToApply.size(), peptideSeq.length()));
			for (int numSequenceChanges = 0; numSequenceChanges <= maxNumberOfVariationsInPeptide; numSequenceChanges++) {
				final Iterator<int[]> combinationsIterator = CombinatoricsUtils
						.combinationsIterator(listOfSequenceChangesToApply.size(), numSequenceChanges);
				while (combinationsIterator.hasNext()) {

					// in order to not repeat locations
					final SetOfIntRanges setOfPositions = new SetOfIntRanges();
					int[] indexes = combinationsIterator.next();
					for (final int index : indexes) {
						final SequenceChange sequenceChange = listOfSequenceChangesToApply.get(index);
						if (!setOfPositions.contains(sequenceChange.getPositionsInPeptide())) {
							setOfPositions.add(sequenceChange.getPositionsInPeptide());
						} else {
							// if there are more than one sequence change in a
							// single position, discard this combination
							indexes = null;
							break;
						}
					}

					if (indexes != null && indexes.length > 0) {
						// iterate over them to construct the key
						final SequenceWithModification sequenceChanged = getModifiedSequence(peptideSeq,
								listOfSequenceChangesToApply, 0, indexes, extendedAssignMass, proteinSequence);
						if (sequenceChanged.getSequenceAfterModification().length() < Constants.MIN_PEP_LENGTH) {
							continue;
						}
						// check whether after the change in the sequence the
						// resulting peptide is valid (comparing it to the ones
						// coming from digesting the protein after the change)

						// first apply changes to protein
						final SequenceWithModification modifiedProteinSequence = getModifiedSequence(proteinSequence,
								listOfSequenceChangesToApply, positionOfPeptideInProtein - 1, indexes,
								extendedAssignMass, proteinSequence);
						final String newProteinSequence = modifiedProteinSequence.getSequenceAfterModification();
						final Collection<String> cleaves = enzyme.cleave(newProteinSequence,
								ProteoformDBIndexer.MIN_PEPTIDE_LENGHT, ProteoformDBIndexer.MAX_PEPTIDE_LENGTH,
								new THashSet<String>());
						// only keep the ones overlapping with the position of the
//						final Iterator<String> iterator = cleaves.iterator();
//						while (iterator.hasNext()) {
//							final String cleave = iterator.next();
//							final TIntArrayList starts = StringUtils.allPositionsOf(newProteinSequence, cleave);
//							boolean overlaps = false;
//							for (final int start : starts.toArray()) {
//								final int end = start + cleave.length();
//								final int startPositionOfOriginalPeptideOnProtein = positionOfPeptideInProtein;
//								final int endPositionOfOriginalPeptideOnProtein = startPositionOfOriginalPeptideOnProtein
//										+ peptideSeq.length();
//								if (start >= startPositionOfOriginalPeptideOnProtein
//										&& start <= endPositionOfOriginalPeptideOnProtein) {
//									overlaps = true;
//								}
//								if (start <= startPositionOfOriginalPeptideOnProtein
//										&& end >= startPositionOfOriginalPeptideOnProtein) {
//									overlaps = true;
//								}
//								if (start <= endPositionOfOriginalPeptideOnProtein
//										&& end >= startPositionOfOriginalPeptideOnProtein) {
//									overlaps = true;
//								}
//							}
//							if (!overlaps) {
//								iterator.remove();
//							}
//						}

						// then check if one of the cleavage products is this
						// peptide
						if (peptideIsInCleavagesProducts(sequenceChanged.getSequenceAfterModification(), cleaves)) {
							ret.add(sequenceChanged);
						} else {
							// if it is not valid, take peptides starting with
							// the original peptide sequence:
							final List<String> newPeptides = getCleavageProductsStartingWithPeptide(
									sequenceChanged.getSequenceAfterModification(), cleaves);
							// and for each of them, apply all the variations
							for (final String newPeptideSequence : newPeptides) {
								final String newPeptideSequenceOriginal = getOriginalPeptideSequence(
										sequenceChanged.getSequenceWithModification(), newPeptideSequence);
								if (
//										newPeptideSequenceOriginal.equals(peptideSeq) 										|| 
								peptidesProcessed.contains(newPeptideSequenceOriginal)) {
									continue;
								}
								if (!proteinSequence.contains(newPeptideSequenceOriginal)) {
									// this happens when the peptide sequenceChanged.getSequenceAfterModification()
									// is in several positions in the protein sequence. then, the extension of that
									// peptide may not match with the protein<br>
									// just do not use it
									continue;
								}
								ret.addAll(getAllCombinationsForPeptide(newPeptideSequenceOriginal, proteinSequence,
										positionOfPeptideInProtein, enzyme, listOfSequenceChangesToApply,
										maxNumVariationsPerPeptide, extendedAssignMass, peptidesProcessed));
							}
						}

					}
				}
			}
		}
		staticCallDepp--;
		return ret;

	}

	private int getNumMissedClavages(String peptideSeq, MyEnzyme enzyme) {
		int numMissedClavages = 0;
		// loop until the second last AA (we dont count the last one) to be a
		// missedclavage
		for (int i = 0; i < peptideSeq.length() - 1; i++) {
			if (isOneOfThese(enzyme.getCleavage(), peptideSeq.charAt(i))) {
				// i + 1 is always fine, because we loop until the second last
				// AA
				if (!isOneOfThese(enzyme.getRestrict(), peptideSeq.charAt(i + 1))) {
					numMissedClavages++;
				}
			}
		}
		return numMissedClavages;
	}

	private boolean isOneOfThese(char[] chars, char charAt) {
		if (chars != null) {
			for (final char c : chars) {
				if (c == charAt) {
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * Having a peptide that changed as<br>
	 * ABCD[R->C] <br>
	 * and a newPeptide that is the result of the change and something else as:<br>
	 * ABCDCEFGK<br>
	 * returns the sequence ABCDREFGK where the change has not been applied yet
	 * Takes into account that it may have multiple variants and PTM changes are
	 * ignored
	 * 
	 * @param sequenceChanged
	 * @param newPeptideSequence
	 * @return
	 */
	public static String getOriginalPeptideSequence(String sequenceWithModification, String newPeptideSequence) {
		final StringBuilder sb = new StringBuilder();
		// common beginning would be ABDC in the example

		final String commonBeginning = StringUtils.getCommonBeginning(sequenceWithModification, newPeptideSequence);

		sb.append(commonBeginning);
		int index1 = commonBeginning.length();
		int index2 = commonBeginning.length();
		boolean variantOpen = false;
		boolean variantOriginal = false;
		boolean ptmOpen = false;
		while (index2 < newPeptideSequence.length()) {
			if (index1 < sequenceWithModification.length()) {
				final char c1 = sequenceWithModification.charAt(index1);
				if (c1 == '[') {
					// check whether is a variant modification or a PTM
					final int indexOf = sequenceWithModification.indexOf("->", index1);
					if (indexOf == -1 || indexOf > sequenceWithModification.indexOf("]", index1)) {
						// it is a PTM
						ptmOpen = true;
					}
					variantOpen = true;
					variantOriginal = true;

				} else if (c1 == ']') {
					variantOpen = false;
					ptmOpen = false;
					variantOriginal = false;

				} else if (c1 == '-') {
					// do nothing
				} else if (c1 == '>') {
					variantOriginal = false;
				} else {
					if (!ptmOpen) {
						if (!variantOpen || (variantOpen && variantOriginal)) {
							sb.append(c1);
						}
						if (!variantOpen || variantOpen && !variantOriginal) {
							index2++;
						}

					}
				}
			} else {
				sb.append(newPeptideSequence.charAt(index2));
				index2++;
			}
			index1++;

		}
		return sb.toString();
//		// change string would be R->C in the example
//		final String changeString = sequenceWithModification.substring(sequenceWithModification.indexOf("[") + 1,
//				sequenceWithModification.indexOf("]"));
//		// change would be C in the example
//		String change = "";
//		if (changeString.contains("->")) {
//			// original would be R in the example
//			final String original = changeString.substring(0, changeString.indexOf("->"));
//			sb.append(original);
//			if (!changeString.endsWith("->")) {
//				change = changeString.substring(changeString.indexOf("->") + 2);
//			}
//
//		}
//		sb.append(newPeptideSequence.substring(commonBeginning.length() + change.length()));
////		final String ret = FastaParser.cleanSequenceAndNotApplySequenceVariances(sequenceWithModification);
//		return ret;
	}

	public static void main(String[] arg) {
		System.out.println(getOriginalPeptideSequence("ABCD[R->C]", "ABCDCEFGK"));
	}

	private List<String> getCleavageProductsStartingWithPeptide(String peptideSeq, Collection<String> cleaves) {
		final List<String> ret = new ArrayList<String>();
		for (final String cleavageProduct : cleaves) {
			if (cleavageProduct.startsWith(peptideSeq)) {
				ret.add(cleavageProduct);
			}
		}
		return ret;
	}

	private boolean peptideIsInCleavagesProducts(String peptide, Collection<String> cleaves) {
		return cleaves.contains(peptide);
	}

	private SequenceWithModification getModifiedSequence(String peptideSeq, List<SequenceChange> sequenceChanges,
			int offset, int[] indexes, ExtendedAssignMass extendedAssignMass, String proteinSequence)
			throws IOException {
		final StringBuilder sb = new StringBuilder();
		final List<Integer> indexesList = new ArrayList<Integer>();
		if (indexes.length > 1) {

			// sort sequenceChanges by position
			final Map<SequenceChange, Integer> map = new HashMap<SequenceChange, Integer>();
			for (final int index : indexes) {
				map.put(sequenceChanges.get(index), index);
			}
			final List<SequenceChange> list = new ArrayList<SequenceChange>();
			list.addAll(map.keySet());
			Collections.sort(list);

			for (final SequenceChange sequenceChange : list) {
				final int index = map.get(sequenceChange);
				indexesList.add(index);
			}
		} else {
			indexesList.add(indexes[0]);
		}
		int i = 0;

		SequenceChange currentSequenceChange = null;
		if (i < indexesList.size()) {
			currentSequenceChange = sequenceChanges.get(indexesList.get(i));

		}

		for (int index = 0; index < peptideSeq.length(); index++) {
			if (currentSequenceChange == null) {
				sb.append(peptideSeq.substring(index));
				break;
			}
			if (currentSequenceChange != null
					&& index + 1 < currentSequenceChange.getFirstPositionOfChangeInPeptide() + offset
					&& peptideSeq.length() > currentSequenceChange.getFirstPositionOfChangeInPeptide() + offset - 1) {

				final String substring = peptideSeq.substring(index,
						currentSequenceChange.getFirstPositionOfChangeInPeptide() + offset - 1);
				sb.append(substring);
				index = currentSequenceChange.getFirstPositionOfChangeInPeptide() + offset - 1;

			}

			if (currentSequenceChange != null
					&& currentSequenceChange.getFirstPositionOfChangeInPeptide() + offset - 1 == index) {

				if (currentSequenceChange.getChange() != null) {
					final String change = currentSequenceChange.getChange();
					final String original = currentSequenceChange.getOriginal();
					sb.append("[" + original + "->" + change + "]");
					index += original.length() - 1;

				} else {
					sb.append(peptideSeq.charAt(index));
					final String key = currentSequenceChange.getKey();
					ptms.add(key);
					sb.append("[" + SequenceChange.df.format(currentSequenceChange.getMassChange()) + "]");
				}

				if (++i < indexesList.size()) {
					currentSequenceChange = sequenceChanges.get(indexesList.get(i));
				} else {
					currentSequenceChange = null;
				}
			}
		}

		final SequenceWithModification peptideWithModification = new SequenceWithModification(peptideSeq, sb.toString(),
				extendedAssignMass, proteinSequence);
		if (peptideWithModification.getSequenceBeforeModification().equals("VLATVTKPVGGDQNGGTRVVK")) {
			log.info(peptideWithModification);
		}
		return peptideWithModification;
	}

}
