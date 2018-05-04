package proteoform_dbindex;

import static org.junit.Assert.fail;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.junit.Test;
import org.proteored.miapeapi.psimod.PSIModOBOPlainTextReader;

import com.compomics.dbtoolkit.io.EnzymeLoader;
import com.compomics.dbtoolkit.io.implementations.FASTADBLoader;
import com.compomics.util.general.UnknownElementMassException;
import com.compomics.util.protein.AASequenceImpl;
import com.compomics.util.protein.Protein;

import edu.scripps.yates.annotations.uniprot.UniprotCVTermCode;
import edu.scripps.yates.annotations.uniprot.UniprotPTMCVReader;
import edu.scripps.yates.annotations.uniprot.UniprotPTMCVTerm;
import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.annotations.uniprot.proteoform.Proteoform;
import edu.scripps.yates.annotations.uniprot.proteoform.ProteoformType;
import edu.scripps.yates.annotations.uniprot.proteoform.ProteoformUtil;
import edu.scripps.yates.annotations.uniprot.proteoform.UniprotPTM;
import edu.scripps.yates.annotations.uniprot.proteoform.xml.UniprotProteoformRetrieverFromXML;
import edu.scripps.yates.proteoform_dbindex.model.PhosphositeDB;
import edu.scripps.yates.proteoform_dbindex.model.SequenceChange;
import edu.scripps.yates.proteoform_dbindex.util.SetOfIntRanges;
import edu.scripps.yates.utilities.masses.FormulaCalculator;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;

public class ProteoformStatistics {

	private static final String AT = "at_";
	private final int nMissedCleavages = 3;
	private final int nMin = 6;
	private final int nMax = 30;
	private final String uniprotVersion = "2018_03";
	final String enzymeName = "Trypsin";
	private final int MAX_NUM_VARIATIONS_PER_PEPTIDE = 5;
	private final static FormulaCalculator calc = new FormulaCalculator();
	private PSIModOBOPlainTextReader psiModReader;
	// new CalculateDistributions(aInputFile, aLengthInterval, aMassInterval,
	// aEnzyme).
	private UniprotProteinLocalRetriever uplr;
	private final File uniprotReleasesFolder = new File("Z:\\share\\Salva\\data\\uniprotKB");
	private final boolean retrieveIsoforms = true;
	private final static Set<String> ptms = new HashSet<String>();
	private final static Set<String> ignoredPTMs = new HashSet<String>();
	private PhosphositeDB phosphositeDB;
	private final boolean useUniprot = false;
	private final boolean usePhosphosite = !useUniprot;

	@Test
	public void proteoformStatistics() {
		try {
			phosphositeDB = new PhosphositeDB("human");
			System.out.println(phosphositeDB.getTotalPhosphoSites() + " phosphosites in PhosphoSiteDB");
			psiModReader = new PSIModOBOPlainTextReader();
			uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
			final UniprotProteoformRetrieverFromXML proteoformRetriever = new UniprotProteoformRetrieverFromXML(uplr,
					uniprotVersion);
			uplr.setCacheEnabled(false);
			final com.compomics.util.protein.Enzyme enzyme = EnzymeLoader.loadEnzyme(enzymeName,
					String.valueOf(nMissedCleavages));

			proteoformRetriever.setRetrieveIsoforms(retrieveIsoforms);
			final File fastaFile = new File("D:\\Downloads\\casimir\\Uniprot_Human.fasta");
			final String cleavageKey = Arrays.toString(enzyme.getCleavage()).replace("[", "").replace("]", "")
					.replace(",", "").replace(" ", "");
			final String analysisKey = "cleav_" + cleavageKey + "_miss_" + enzyme.getMiscleavages();
			final FileWriter log = new FileWriter(new File(fastaFile.getParent() + File.separator
					+ FilenameUtils.getBaseName(fastaFile.getAbsolutePath()) + "_stats_" + analysisKey + ".txt"));
			final FASTADBLoader loader = new FASTADBLoader();

			loader.load(fastaFile.getAbsolutePath());
			Protein protein = null;
			final Set<String> uniquePeptidesWithVariationsTotal = new HashSet<String>();
			final Set<String> uniquePeptidesWithNoVariationsTotal = new HashSet<String>();
			final ProgressCounter proteinCounter = new ProgressCounter(loader.countNumberOfEntries(),
					ProgressPrintingType.PERCENTAGE_STEPS, 0);
			final List<Integer> numberOfPeptidesWithNoVariationsPerProtein = new ArrayList<Integer>();
			final List<Integer> numberOfVariationsPerPeptide = new ArrayList<Integer>();
			final List<Integer> numberOfPeptidesWithVariationsPerProtein = new ArrayList<Integer>();
			while ((protein = loader.nextProtein()) != null) {
				final Set<String> uniquePeptidesWithVariationsInProtein = new HashSet<String>();
				final Set<String> uniquePeptidesWithNoVariationsInProtein = new HashSet<String>();
				proteinCounter.increment();
				final String printIfNecessary = proteinCounter.printIfNecessary();
				if (!"".equals(printIfNecessary)) {
					System.out.println("Protein: " + printIfNecessary);
				}
				final String fullHeaderWithAddenda = protein.getHeader().getFullHeaderWithAddenda();
				if (fullHeaderWithAddenda.contains("Reverse")) {
					continue;
				}
				final Protein[] mainProteinPeptides = enzyme.cleave(protein, nMin, nMax);
				// increase number of peptides

				final String accession = protein.getHeader().getAccession();
				final String proteinSeq = protein.getSequence().getSequence();
				Map<String, List<Proteoform>> proteoformMap = null;
				if (useUniprot) {
					proteoformMap = proteoformRetriever.getProteoforms(accession);
				} else if (usePhosphosite) {
					proteoformMap = loadProteoformMapFromPhosphosite(accession);
				}
				if (proteoformMap.size() > 1) {
					fail();
				}
				// if (!proteoformMap.containsKey(accession)) {
				// // continue;
				// }
				final List<Proteoform> proteoforms = proteoformMap.get(accession);
				final List<Proteoform> nonIsoformProteoforms = ProteoformUtil
						.getProteoformsDifferentThan(ProteoformType.ISOFORM, proteoforms);
				final List<Proteoform> isoformProteoforms = ProteoformUtil.getProteoformsAs(ProteoformType.ISOFORM,
						proteoforms);

				final Map<Integer, List<Proteoform>> proteoformsByPositionInMainProtein = getProteoformsByPositionInProtein(
						nonIsoformProteoforms);
				for (int i = 0; i < isoformProteoforms.size() + 1; i++) {
					Proteoform isoform = null;
					if (i < isoformProteoforms.size()) {
						isoform = isoformProteoforms.get(i);
					}
					// here, if isoform is null, we are in the latest index 'i'
					// which is when we apply all the nonIsoforms to the main
					// protein entry
					Protein[] peptides;
					if (isoform == null) {
						peptides = new Protein[mainProteinPeptides.length];
						System.arraycopy(mainProteinPeptides, 0, peptides, 0, mainProteinPeptides.length);
					} else {
						// digest isoform sequence
						peptides = enzyme.cleave(new Protein(new AASequenceImpl(isoform.getSeq())), nMin, nMax);
					}

					// final ProgressCounter counter = new
					// ProgressCounter(peptides.length,
					// ProgressPrintingType.PERCENTAGE_STEPS, 0);
					final Set<String> peptideKeys = new HashSet<String>();
					for (final Protein peptide : peptides) {
						// counter.increment();
						// final String printIfNecessary1 =
						// counter.printIfNecessary();
						// if (!"".equals(printIfNecessary1)) {

						// System.out.println(proteinProgress + "\t" + accession
						// + "\t" + printIfNecessary1);
						// }
						final String peptideSeq = peptide.getSequence().getSequence();
						uniquePeptidesWithNoVariationsInProtein.add(peptideSeq);
						uniquePeptidesWithNoVariationsTotal.add(peptideSeq);
						int positionInProtein = -1;
						if (isoform == null) {
							positionInProtein = proteinSeq.indexOf(peptideSeq) + 1;
						} else {
							positionInProtein = isoform.getSeq().indexOf(peptideSeq) + 1;
						}
						// in order to save time:
						// peptideSeq + positionInMainProtein is a key, and we
						// avoid this to be repeated because it will return the
						// same things
						final String key = accession + "|" + peptideSeq + "|" + proteinSeq.indexOf(peptideSeq) + 1;
						if (peptideKeys.contains(key)) {
							continue;
						} else {
							peptideKeys.add(key);
						}
						final List<SequenceChange> sequenceChanges = getSequenceChangesInPeptide(peptideSeq,
								positionInProtein, proteoformsByPositionInMainProtein, isoform);
						final Set<String> combinations = getAllCombinationsForPeptide(peptideSeq, sequenceChanges);
						uniquePeptidesWithVariationsTotal.addAll(combinations);
						uniquePeptidesWithVariationsInProtein.addAll(combinations);
						numberOfVariationsPerPeptide.add(combinations.size());
						final String acc = isoform != null ? isoform.getId() : protein.getHeader().getAccession();
						log.write(peptideSeq + "\tACC=" + acc + "\tpositionInProtein=[" + positionInProtein + "-"
								+ (positionInProtein + peptideSeq.length() - 1) + "]\t# proteoforms="
								+ sequenceChanges.size());
						try {
							final int numberOfNaturalVariants = getNumNaturalVariants(sequenceChanges);
							final double numberOfCombinations = combinations.size();

							log.write("\t" + numberOfCombinations + "\t parent mass combinations " + " ("
									+ numberOfNaturalVariants + " natural variants + "
									+ (sequenceChanges.size() - numberOfNaturalVariants)
									+ " other variant annotations)");
						} catch (final IllegalArgumentException e) {

						}
						for (final SequenceChange sequenceChange : sequenceChanges) {
							log.write("\t" + sequenceChange.getProteoformID());
							if (sequenceChange.isPtm()) {
								log.write("_at_"
										+ (sequenceChange.getFirstPositionOfChangeInPeptide() + positionInProtein - 1)
										+ ":" + sequenceChange.getKey());
							}
						}

						log.write("\n");
						if (positionInProtein < 1) {
							throw new IllegalArgumentException("ERROR " + peptideSeq + " not found");
						}

					}
				}
				numberOfPeptidesWithNoVariationsPerProtein.add(uniquePeptidesWithNoVariationsInProtein.size());
				numberOfPeptidesWithVariationsPerProtein.add(uniquePeptidesWithVariationsInProtein.size());
			}
			final StringBuilder result = new StringBuilder();
			result.append("Analsysis at " + new Date() + "\tEnzyme:" + enzyme.getTitle() + "("
					+ Arrays.toString(enzyme.getCleavage()) + ", " + enzyme.getMiscleavages() + ")\n");
			result.append(
					uniquePeptidesWithNoVariationsTotal.size() + " different peptides (no counting variations)\n");
			result.append(uniquePeptidesWithVariationsTotal.size() + " unique peptides counting variations\n");

			final double avgNumberOfPeptidesPerProtein = Maths
					.mean(numberOfPeptidesWithNoVariationsPerProtein.toArray(new Integer[0]));
			final double avgNumberOfVariationPeptidesPerPeptide = Maths
					.mean(numberOfVariationsPerPeptide.toArray(new Integer[0]));
			final double avgNumberOfVariationPeptidesPerProtein = Maths
					.mean(numberOfPeptidesWithVariationsPerProtein.toArray(new Integer[0]));
			result.append("Avg number of peptides per protein " + avgNumberOfPeptidesPerProtein + "\n");
			result.append(
					"Avg number of peptide variations per peptide " + avgNumberOfVariationPeptidesPerPeptide + "\n");
			result.append(
					"Avg number of peptide variations per protein " + avgNumberOfVariationPeptidesPerProtein + "\n");
			final double ratio = uniquePeptidesWithVariationsTotal.size() * 1.0
					/ uniquePeptidesWithNoVariationsTotal.size();
			result.append("Increased by : " + ratio + "\n\n");
			result.append("Number of different PTMs: " + ptms.size() + "\n\n");
			for (final String ptm : ptms) {
				result.append(ptm + "\n");
			}

			result.append(
					"Number of different PTMs ignored because no DELTA mass was found: " + ignoredPTMs.size() + "\n\n");
			for (final String ptm : ignoredPTMs) {
				result.append("ignored\t" + ptm + "\n");
			}
			log.write(result.toString());
			System.out.println(result);
			log.close();
		} catch (final IOException e) {
			e.printStackTrace();
		} catch (final UnknownElementMassException e) {
			e.printStackTrace();
		}
	}

	private Map<String, List<Proteoform>> loadProteoformMapFromPhosphosite(String accession) {
		final Map<String, List<Proteoform>> ret = new HashMap<String, List<Proteoform>>();
		final List<Integer> phosphorilatedPositions = phosphositeDB.getPhosphorilatedPositions(accession);
		if (phosphorilatedPositions != null) {
			final List<Proteoform> list = new ArrayList<Proteoform>();
			ret.put(accession, list);
			final String proteinSeq = phosphositeDB.getProteinSeq(accession);

			for (final Integer position : phosphorilatedPositions) {
				final char aa = proteinSeq.charAt(position - 1);
				final Proteoform proteoform = new Proteoform(accession, proteinSeq, accession, proteinSeq,
						"Phospho(" + aa + ")_at_" + position, null, null, ProteoformType.PTM, true);

				final UniprotPTMCVTerm uniprotPTMCVTerm = UniprotPTMCVReader.getInstance().getPtmsByID(getPTMID(aa));

				final UniprotPTM uniprotPTM = new UniprotPTM(position, uniprotPTMCVTerm);
				proteoform.addPTM(uniprotPTM);
				list.add(proteoform);
			}

		}
		return ret;
	}

	private String getPTMID(char aa) {
		if (aa == 'y') {
			return "Phosphotyrosine";
		} else if (aa == 's') {
			return "Phosphoserine";
		} else if (aa == 't') {
			return "Phosphothreonine";

		}
		return null;
	}

	private Set<String> getAllCombinationsForPeptide(String peptideSeq, List<SequenceChange> sequenceChanges) {
		final Set<String> ret = new HashSet<String>();
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

			final int maxNumberOfVariationsInPeptide = Math.min(MAX_NUM_VARIATIONS_PER_PEPTIDE,
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
					if (indexes == null) {
						continue;
					}
					if (indexes.length > 0) {
						// iterate over them to construct the key
						final String sequenceChanged = getSequenceChange(peptideSeq, listOfSequenceChangesToApply,
								indexes);
						ret.add(sequenceChanged);
					} else {
						ret.add(peptideSeq);
					}
				}
			}
		}
		return ret;
	}

	private String getSequenceChange(String peptideSeq, List<SequenceChange> sequenceChanges, int[] indexes) {
		final StringBuilder sb = new StringBuilder();
		int i = 0;

		SequenceChange currentSequenceChange = null;
		if (i < indexes.length) {
			currentSequenceChange = sequenceChanges.get(indexes[i]);
		}
		for (int index = 0; index < peptideSeq.length(); index++) {
			if (currentSequenceChange == null) {
				sb.append(peptideSeq.substring(index));
				break;
			}
			if (currentSequenceChange != null
					&& index + 1 < currentSequenceChange.getFirstPositionOfChangeInPeptide()) {
				sb.append(peptideSeq.substring(index, currentSequenceChange.getFirstPositionOfChangeInPeptide() - 1));
				index = currentSequenceChange.getFirstPositionOfChangeInPeptide() - 1;

			}

			if (currentSequenceChange != null
					&& currentSequenceChange.getFirstPositionOfChangeInPeptide() - 1 == index) {

				if (currentSequenceChange.getChange() != null) {
					final String change = currentSequenceChange.getChange();
					sb.append(change);
				} else {
					sb.append(peptideSeq.charAt(index));
					final String key = currentSequenceChange.getKey();
					ptms.add(key);
					sb.append("[" + key + "]");
				}

				if (++i < indexes.length) {
					currentSequenceChange = sequenceChanges.get(indexes[i]);
				} else {
					currentSequenceChange = null;
				}
			}
		}
		return sb.toString();
	}

	private int getNumNaturalVariants(List<SequenceChange> sequenceChanges) {
		int ret = 0;
		for (final SequenceChange sequenceChange : sequenceChanges) {
			if (sequenceChange.getProteoformType() == ProteoformType.NATURAL_VARIANT) {
				ret++;
			}
		}
		return ret;
	}

	private Map<Integer, List<Proteoform>> getProteoformsByPositionInProtein(List<Proteoform> proteoforms) {
		final Map<Integer, List<Proteoform>> proteoformsByPositionInProtein = new HashMap<Integer, List<Proteoform>>();
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
	 * @param proteoformsByPositionInMainProtein
	 * @param isoform
	 *            can be null
	 * @return
	 * @throws UnknownElementMassException
	 */
	private List<SequenceChange> getSequenceChangesInPeptide(String peptideSeq, int peptideInitInProtein,
			Map<Integer, List<Proteoform>> proteoformsByPositionInMainProtein, Proteoform isoform)
			throws UnknownElementMassException {
		if (isoform == null) {
			return getSequenceChangesInPeptideTowardsMainProteinform(peptideSeq, peptideInitInProtein,
					proteoformsByPositionInMainProtein);
		} else {
			return getSequenceChangesInPeptideTowardsIsoform(peptideSeq, proteoformsByPositionInMainProtein, isoform);
		}
	}

	/**
	 * 
	 * @param peptideSeq
	 * @param peptideInitInOriginalProtein
	 * @param proteoformsByPositionInProtein
	 * @param isoform
	 *            can be null
	 * @return
	 * @throws UnknownElementMassException
	 */
	private List<SequenceChange> getSequenceChangesInPeptideTowardsMainProteinform(String peptideSeq,
			int peptideInitInOriginalProtein, Map<Integer, List<Proteoform>> proteoformsByPositionInProtein)
			throws UnknownElementMassException {
		final List<SequenceChange> ret = new ArrayList<SequenceChange>();

		final List<Integer> positionsInOrder = new ArrayList<Integer>();
		positionsInOrder.addAll(proteoformsByPositionInProtein.keySet());
		Collections.sort(positionsInOrder);
		for (final Integer proteoformPositionInProtein : positionsInOrder) {
			final int peptideLength = peptideSeq.length();
			final int peptideEndInProtein = peptideInitInOriginalProtein + peptideLength - 1;
			if (peptideEndInProtein < proteoformPositionInProtein) {
				break;
			}
			if (proteoformPositionInProtein >= peptideInitInOriginalProtein
					&& proteoformPositionInProtein <= peptideEndInProtein) {
				final List<Proteoform> proteoforms = proteoformsByPositionInProtein.get(proteoformPositionInProtein);
				for (final Proteoform proteoform : proteoforms) {

					final int positionInPeptide = proteoformPositionInProtein - peptideInitInOriginalProtein + 1;
					final String change = getProteoformSequenceChange(proteoform);
					if (change != null && proteoformPositionInProtein + change.length() - 1 > peptideEndInProtein) {
						// the change is out of the peptide
						continue;
					}
					final List<Double> massChanges = getProteoformPTMMassChanges(proteoformPositionInProtein,
							proteoform);

					final String original = getProteoformSequenceOriginalChange(proteoform,
							proteoformPositionInProtein);

					if (massChanges != null) {
						for (final Double massChange : massChanges) {
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
		// because implements Comparable we can sort
		Collections.sort(ret);
		return ret;

	}

	/**
	 * 
	 * @param isoformPeptideSeq
	 * @param proteoformsByPositionInMainProtein
	 * @param isoform
	 *            is not null
	 * 
	 * @return
	 * @throws UnknownElementMassException
	 */
	private List<SequenceChange> getSequenceChangesInPeptideTowardsIsoform(String isoformPeptideSeq,
			Map<Integer, List<Proteoform>> proteoformsByPositionInMainProtein, Proteoform isoform)
			throws UnknownElementMassException {
		final List<SequenceChange> ret = new ArrayList<SequenceChange>();

		final List<Integer> positionsInOrder = new ArrayList<Integer>();
		positionsInOrder.addAll(proteoformsByPositionInMainProtein.keySet());
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
						final String change = getProteoformSequenceChange(proteoformToApply);
						if (change != null
								&& proteoformPositionInMainProtein + change.length() - 1 > peptideEndInMainProtein) {
							// the change is out of the peptide
							continue;
						}
						final List<Double> massChanges = getProteoformPTMMassChanges(proteoformPositionInMainProtein,
								proteoformToApply);

						final String original = getProteoformSequenceOriginalChange(proteoformToApply,
								proteoformPositionInMainProtein);

						if (massChanges != null) {
							for (final Double massChange : massChanges) {
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

	private List<Double> getProteoformPTMMassChanges(int proteoformPositionInProtein, Proteoform proteoform)
			throws UnknownElementMassException {
		if (proteoform.getPtms() != null && !proteoform.getPtms().isEmpty()) {
			final List<Double> ret = new ArrayList<Double>();
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

	public Double getMassChange(UniprotPTM uniprotPTM) throws UnknownElementMassException {

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

	private String getProteoformSequenceOriginalChange(Proteoform proteoform, Integer proteoformPositionInProtein) {
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
			final String MISSING_AT = "_missing_at_";
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

}
