package edu.scripps.yates.proteoform_dbindex;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.compomics.util.general.UnknownElementMassException;
import com.compomics.util.protein.AASequenceImpl;
import com.compomics.util.protein.Enzyme;
import com.compomics.util.protein.Protein;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.annotations.uniprot.proteoform.Proteoform;
import edu.scripps.yates.annotations.uniprot.proteoform.ProteoformType;
import edu.scripps.yates.annotations.uniprot.proteoform.ProteoformUtil;
import edu.scripps.yates.annotations.uniprot.proteoform.xml.UniprotProteoformRetrieverFromXML;
import edu.scripps.yates.dbindex.Constants;
import edu.scripps.yates.dbindex.DBIndexStore.FilterResult;
import edu.scripps.yates.dbindex.DBIndexStoreException;
import edu.scripps.yates.dbindex.DBIndexer;
import edu.scripps.yates.dbindex.ProteinCache;
import edu.scripps.yates.dbindex.io.DBIndexSearchParams;
import edu.scripps.yates.dbindex.model.AssignMass;
import edu.scripps.yates.proteoform_dbindex.model.ExtendedAssignMass;
import edu.scripps.yates.proteoform_dbindex.model.PTMCodeObj;
import edu.scripps.yates.proteoform_dbindex.model.PhosphositeDB;
import edu.scripps.yates.proteoform_dbindex.model.SequenceChange;
import edu.scripps.yates.proteoform_dbindex.model.SequenceWithModification;
import edu.scripps.yates.proteoform_dbindex.util.ProteoformDBIndexUtil;
import edu.scripps.yates.utilities.fasta.FastaParser;

public class ProteoformDBIndexer extends DBIndexer {
	private static final org.apache.log4j.Logger logger = org.apache.log4j.Logger.getLogger(ProteoformDBIndexer.class);
	////////////////////////////////////////////
	// TODO
	// TO PUT THESE PARAMETERS IN sParam? :
	public static final int MIN_PEPTIDE_LENGHT = 6;
	public static final int MAX_PEPTIDE_LENGTH = 40;
	////////////////////////////////////////////
	private final boolean useUniprot;
	private final boolean usePhosphosite;
	private final PhosphositeDB phosphositeDB;
	private final UniprotProteoformRetrieverFromXML proteoformRetriever;
	private final int maxNumVariationsPerPeptide;
	private Enzyme enzyme;
	private final ExtendedAssignMass extendedAssignMass;

	public ProteoformDBIndexer(DBIndexSearchParams sparam, IndexerMode indexerMode, boolean useUniprot,
			boolean usePhosphosite, String phosphoSiteSpecies, File uniprotReleasesFolder, String uniprotVersion,
			int maxNumVariationsPerPeptide) throws IOException {
		super(sparam, indexerMode,
				new ProteoformDBIndexStoreSQLiteMult(sparam, false,
						ExtendedAssignMass.getInstance(sparam.isUseMonoParent(),
								new File(new File(sparam.getFullIndexFileName()).getAbsolutePath() + File.separator
										+ PTMCodeObj.FILE_NAME))));
		this.usePhosphosite = usePhosphosite;
		this.useUniprot = useUniprot;
		if (usePhosphosite) {
			if (phosphoSiteSpecies != null && !"".equals(phosphoSiteSpecies)) {
				phosphositeDB = new PhosphositeDB(phosphoSiteSpecies);
			} else {
				throw new IllegalArgumentException(
						"If usePhosphoSiteDB is true, the species parameter is required and it is null or empty");
			}
		} else {
			phosphositeDB = null;
		}
		final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
		proteoformRetriever = new UniprotProteoformRetrieverFromXML(uplr, uniprotVersion);

		this.maxNumVariationsPerPeptide = maxNumVariationsPerPeptide;
		extendedAssignMass = ExtendedAssignMass.getInstance(sparam.isUseMonoParent(), new File(
				new File(sparam.getFullIndexFileName()).getAbsolutePath() + File.separator + PTMCodeObj.FILE_NAME));
	}

	/**
	 * Cut a Fasta protein sequence according to params spec and index the
	 * protein and generated sequences
	 *
	 * Should be called once per unique Fasta sequence
	 *
	 * @param fasta
	 *            fasta to index
	 * @throws IOException
	 */
	@Override
	protected void cutSeq(final String proteinFastaHeader, String canonicalProtSeq) throws IOException {
		final String protAccession = FastaParser.getUniProtACC(proteinFastaHeader);
		if (protAccession == null) {
			throw new IllegalArgumentException(
					"Uniprot accession cannot be extracted from fasta header: '" + proteinFastaHeader + "'");
		}

		final Protein[] canonicalProteinPeptides = digestProtein(canonicalProtSeq);

		Map<String, List<Proteoform>> proteoformMap = new HashMap<String, List<Proteoform>>();
		if (useUniprot) {
			proteoformMap = proteoformRetriever.getProteoforms(protAccession);
		}
		if (usePhosphosite) {
			mergeMaps(proteoformMap,
					ProteoformDBIndexUtil.getInstance().loadProteoformMapFromPhosphoSite(phosphositeDB, protAccession));
		}
		// get proteoforms for the protein
		final List<Proteoform> proteoforms = proteoformMap.get(protAccession);
		// separate isoforms from others
		// others here
		final List<Proteoform> isoformProteoforms = ProteoformUtil.getProteoformsAs(ProteoformType.ISOFORM,
				proteoforms);
		// isoforms here
		final List<Proteoform> nonIsoformProteoforms = ProteoformUtil
				.getProteoformsDifferentThan(ProteoformType.ISOFORM, proteoforms);

		final Map<Integer, List<Proteoform>> nonIsoformsProteoformsByPositionInMainProtein = ProteoformDBIndexUtil
				.getInstance().getProteoformsByPositionInProtein(nonIsoformProteoforms);

		for (int i = 0; i < isoformProteoforms.size() + 1; i++) {
			String proteinSequence = canonicalProtSeq;
			String proteinAccession = protAccession;
			Proteoform isoform = null;
			if (i < isoformProteoforms.size()) {
				isoform = isoformProteoforms.get(i);
				proteinSequence = isoform.getSeq();
				proteinAccession = isoform.getId();
			}
			final int proteinId = proteinCache.addProtein(proteinAccession, proteinSequence);
			final int proteinLength = proteinSequence.length();
			// here, if isoform is null, we are in the latest index 'i'
			// which is when we apply all the nonIsoforms to the main
			// protein entry

			Protein[] peptides;
			if (isoform == null) {
				peptides = new Protein[canonicalProteinPeptides.length];
				System.arraycopy(canonicalProteinPeptides, 0, peptides, 0, canonicalProteinPeptides.length);
			} else {
				// digest isoform sequence
				peptides = enzyme.cleave(new Protein(new AASequenceImpl(proteinSequence)), MIN_PEPTIDE_LENGHT,
						MAX_PEPTIDE_LENGTH);
			}

			final Set<String> peptideKeys = new HashSet<String>();
			for (final Protein peptide : peptides) {
				final String peptideSequence = peptide.getSequence().getSequence();
				// at least one of this AAs has to be in the sequence:
				final char[] mandatoryInternalAAs = sparam.getMandatoryInternalAAs();

				try {
					// peptide start and end in protein
					int start = 0;
					int end = 0;

					if (peptideFilter != null && !peptideFilter.isValid(peptideSequence)) {
						logger.debug("Skipping peptide '" + peptideSequence + "' by filter");
						break;
					}

					final List<SequenceChange> sequenceChanges = ProteoformDBIndexUtil.getInstance()
							.getSequenceChangesInPeptide(peptideSequence, proteinSequence.indexOf(peptideSequence) + 1,
									nonIsoformsProteoformsByPositionInMainProtein, isoform);
					final Set<SequenceWithModification> modifiedPeptides = ProteoformDBIndexUtil.getInstance()
							.getAllCombinationsForPeptide(peptideSequence, proteinSequence,
									proteinSequence.indexOf(peptideSequence) + 1, enzyme, sequenceChanges,
									maxNumVariationsPerPeptide, extendedAssignMass);

					for (final SequenceWithModification modifiedPeptide : modifiedPeptides) {
						// TODO, change this to check the modified peptide
						final String sequenceAfterModification = modifiedPeptide.getSequenceAfterModification();
						if (sequenceAfterModification.length() < Constants.MIN_PEP_LENGTH) { // Constants.MIN_PRECURSOR
							continue;
						}
						if (mandatoryInternalAAs != null) {
							boolean found = false;
							for (final char internalAA : mandatoryInternalAAs) {
								if (sequenceAfterModification.indexOf(internalAA) >= 0) {
									found = true;
								}
							}
							if (!found) {
								continue;
							}
						}
						// Salva added 24Nov2014
						double precMass = 0.0;
						if (sparam.isH2OPlusProtonAdded())
							precMass += AssignMass.H2O_PROTON;
						precMass += AssignMass.getcTerm();
						precMass += AssignMass.getnTerm();
						// CALCULATE MASS
						precMass += modifiedPeptide.getSequenceMassAfterModification();

						if (precMass > sparam.getMaxPrecursorMass()) {
							break;
						}
						if (precMass < sparam.getMinPrecursorMass()) {
							break;
						}
						// check if index will accept it
						final FilterResult filterResult = indexStore.filterSequence(precMass,
								sequenceAfterModification);
						if (filterResult.equals(FilterResult.SKIP_PROTEIN_START)) {
							// bail out earlier as we are no longer
							// interested in this protein starting at
							// start
							break; // move to new start position
						} else if (filterResult.equals(FilterResult.INCLUDE)) {

							start = modifiedPeptide.getProteinSequence().indexOf(sequenceAfterModification);
							end = start + sequenceAfterModification.length() - 1;
							final int resLeftI = start >= Constants.MAX_INDEX_RESIDUE_LEN
									? start - Constants.MAX_INDEX_RESIDUE_LEN : 0;
							final int resLeftLen = Math.min(Constants.MAX_INDEX_RESIDUE_LEN, start);
							final StringBuilder sbLeft = new StringBuilder(Constants.MAX_INDEX_RESIDUE_LEN);
							for (int ii = 0; ii < resLeftLen; ++ii) {
								sbLeft.append(modifiedPeptide.getProteinSequence().charAt(ii + resLeftI));
							}
							final int resRightI = end + 1;
							final int resRightLen = Math.min(Constants.MAX_INDEX_RESIDUE_LEN, proteinLength - end - 1);
							final StringBuilder sbRight = new StringBuilder(Constants.MAX_INDEX_RESIDUE_LEN);
							if (resRightI < proteinLength) {
								for (int jj = 0; jj < resRightLen; ++jj) {
									sbRight.append(modifiedPeptide.getProteinSequence().charAt(jj + resRightI));
								}
							}

							// add -- markers to fill
							// Constants.MAX_INDEX_RESIDUE_LEN length
							final int lLen = sbLeft.length();
							for (int c = 0; c < Constants.MAX_INDEX_RESIDUE_LEN - lLen; ++c) {
								sbLeft.insert(0, '-');
							}
							final int rLen = sbRight.length();
							for (int c = 0; c < Constants.MAX_INDEX_RESIDUE_LEN - rLen; ++c) {
								sbRight.append('-');
							}

							final String resLeft = sbLeft.toString();
							final String resRight = sbRight.toString();

							// before adding the sequence:
							final String key = protAccession + "|" + sequenceAfterModification + "|"
									+ modifiedPeptide.getProteinSequence().indexOf(sequenceAfterModification) + 1;
							if (peptideKeys.contains(key)) {
								continue;
							} else {
								peptideKeys.add(key);
							}

							indexStore.addSequence(precMass, proteinSequence.indexOf(peptideSequence),
									peptideSequence.length(), modifiedPeptide.getSequenceWithModification(), resLeft,
									resRight, proteinId);
						}
					}
					// System.out.println("\t" +
					// peptideSeqString);

				} catch (final DBIndexStoreException e) {
					e.printStackTrace();
					logger.error("Error writing sequence to db index store, ", e);
				} catch (final UnknownElementMassException e) {
					e.printStackTrace();
					logger.error("Error writing sequence to db index store, ", e);
				}
			}
		}

	}

	private Protein[] digestProtein(String canonicalProtSeq) {
		final Enzyme enzyme = getEnzyme();
		final Protein protein = new Protein(new AASequenceImpl(canonicalProtSeq));
		final Protein[] peptides = enzyme.cleave(protein, MIN_PEPTIDE_LENGHT, MAX_PEPTIDE_LENGTH);
		return peptides;
	}

	private Enzyme getEnzyme() {
		if (enzyme == null) {
			enzyme = new Enzyme(null, sparam.getEnzymeResidues(), sparam.getEnzymeNocutResidues(), "Cterm",
					sparam.getMaxMissedCleavages());
		}
		return enzyme;
	}

	private void mergeMaps(Map<String, List<Proteoform>> target, Map<String, List<Proteoform>> toMerge) {
		for (final String key : toMerge.keySet()) {
			if (target.containsKey(key)) {
				target.get(key).addAll(toMerge.get(key));
			} else {
				target.put(key, toMerge.get(key));
			}
		}
	}

	@Override
	protected ProteinCache getProteinCache() {
		if (proteinCache == null) {
			proteinCache = new ProteoformProteinCache(extendedAssignMass,
					new File(new File(sparam.getFullIndexFileName()).getAbsolutePath() + File.separator
							+ ProteoformProteinCache.FILE_NAME));
		}
		return proteinCache;
	}
}
