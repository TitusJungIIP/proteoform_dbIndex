package edu.scripps.yates.utilities.proteomicsmodel;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.annotations.UniprotProteinLocalRetrieverInterface;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotEntryUtil;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.grouping.GroupableProtein;
import edu.scripps.yates.utilities.grouping.PeptideRelation;
import edu.scripps.yates.utilities.proteomicsmodel.adapters.PTMAdapter;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AccessionType;
import edu.scripps.yates.utilities.proteomicsmodel.utils.KeyUtils;
import edu.scripps.yates.utilities.proteomicsmodel.utils.ModelUtils;
import edu.scripps.yates.utilities.sequence.PTMInPeptide;
import edu.scripps.yates.utilities.sequence.PTMInProtein;
import edu.scripps.yates.utilities.sequence.PositionInPeptide;
import edu.scripps.yates.utilities.sequence.PositionInProtein;
import edu.scripps.yates.utilities.sequence.ProteinSequenceUtils;
import edu.scripps.yates.utilities.strings.StringUtils;
import edu.scripps.yates.utilities.util.StringPosition;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntDoubleHashMap;
import gnu.trove.set.hash.THashSet;

public abstract class AbstractPSM implements PSM {
	private final static Logger log = Logger.getLogger(AbstractPSM.class);
	private static Comparator<PTMSite> ptmSiteComparator;
	private Set<Score> scores;
	private Set<Ratio> ratios;
	private Set<Amount> amounts;
	private Set<Condition> conditions;
	private MSRun msRun;
	private String sequence;
	private String identifier;
	private PeptideRelation relation;
	private Set<Protein> proteins;
	private List<PTM> ptms;
	private Float experimentalMH;
	private Float calcMH;
	private Float massErrorPPM;
	private Float totalIntensity;
	private Integer spr;
	private Float ionProportion;
	private Float pi;
	private String fullSequence;
	private Peptide peptide;
	private int uniqueIdentifier;
	private String afterSeq;
	private String beforeSeq;
	private Integer chargeState;
	private String scanNumber;
	private Float rtInMinutes;
	private String searchEngine;
	private Float deltaCn;
	private Float xCorr;
	private Set<String> taxonomies;
	private List<PTMInPeptide> ptmsInPeptide;
	private boolean ignoreTaxonomy;
	private Map<Character, List<PositionInProtein>> positionsInProteinsByQuantifiedAA;
	private String key;
	private final boolean distinguishModifiedSequence;
	private final boolean chargeStateSensible;

	protected AbstractPSM(boolean distinguishModifiedSequence, boolean chargeStateSensible) {
		this.distinguishModifiedSequence = distinguishModifiedSequence;
		this.chargeStateSensible = chargeStateSensible;
	}

	@Override
	public Set<Score> getScores() {
		if (scores == null) {
			scores = new THashSet<Score>();
		}
		return scores;
	}

	@Override
	public boolean addScore(Score score) {
		if (scores == null) {
			scores = new THashSet<Score>();
		}
		return scores.add(score);
	}

	@Override
	public Set<Ratio> getRatios() {
		if (ratios == null) {
			ratios = new THashSet<Ratio>();
		}
		return ratios;
	}

	@Override
	public boolean addRatio(Ratio ratio) {
		if (ratios == null) {
			ratios = new THashSet<Ratio>();
		}
		return ratios.add(ratio);
	}

	@Override
	public Set<Amount> getAmounts() {
		if (amounts == null) {
			amounts = new THashSet<Amount>();
		}
		return amounts;
	}

	@Override
	public boolean addAmount(Amount amount) {
		if (amounts == null) {
			amounts = new THashSet<Amount>();
		}
		addCondition(amount.getCondition());
		return amounts.add(amount);
	}

	@Override
	public Set<Condition> getConditions() {
		if (conditions == null) {
			conditions = new THashSet<Condition>();
		}
		return conditions;
	}

	@Override
	public boolean addCondition(Condition condition) {
		if (conditions == null) {
			conditions = new THashSet<Condition>();
		}
		if (condition != null) {
			return conditions.add(condition);
		}
		return false;
	}

	@Override
	public MSRun getMSRun() {
		return msRun;
	}

	@Override
	public void setMSRun(MSRun msRun) {
		this.msRun = msRun;
	}

	@Override
	public String getSequence() {
		return sequence;
	}

	@Override
	public String getIdentifier() {
		return identifier;
	}

	@Override
	public void setRelation(PeptideRelation relation) {
		this.relation = relation;
	}

	@Override
	public PeptideRelation getRelation() {
		return relation;
	}

	@Override
	public List<GroupableProtein> getGroupableProteins() {
		final List<GroupableProtein> ret = new ArrayList<GroupableProtein>();
		ret.addAll(getProteins());
		return ret;
	}

	@Override
	public Set<Protein> getProteins() {
		if (proteins == null) {
			proteins = new THashSet<Protein>();
		}
		return proteins;
	}

	@Override
	public boolean addProtein(Protein protein, boolean recursively) {
		if (protein != null) {
			if (proteins == null) {
				proteins = new THashSet<Protein>();
			}
			final boolean ret = proteins.add(protein);
			protein.addMSRun(getMSRun());
			if (recursively) {
				protein.addPSM(this, false);
				protein.addPeptide(getPeptide(), false);
				if (getPeptide() != null) {
					getPeptide().addProtein(protein, false);
				}
			}
			return ret;
		}
		return false;
	}

	@Override
	public List<PTM> getPTMs() {
		if (ptms == null) {
			ptms = new ArrayList<PTM>();
		}
		return ptms;
	}

	@Override
	public boolean addPTM(PTM newPtm) {
		if (newPtm != null) {
			if (ptms == null) {
				ptms = new ArrayList<PTM>();
			}
			boolean found = false;
			for (final PTM ptm : ptms) {
				if (ptm.getName().equals(newPtm.getName())) {
					boolean anyPtmIsNew = false;
					for (final PTMSite newPtmSite : newPtm.getPTMSites()) {
						boolean ptmSiteFound = false;
						for (final PTMSite ptmSite : ptm.getPTMSites()) {
							if (newPtmSite.getPosition() == ptmSite.getPosition()) {
								ptmSiteFound = true;
							}
						}
						if (!ptmSiteFound) {
							anyPtmIsNew = true;
						}
					}
					if (!anyPtmIsNew) {
						found = true;
					}
				}
			}
			if (!found) {
				for (final PTMSite ptmSite : newPtm.getPTMSites()) {
					if (ptmSite.getPTMPosition() == PTMPosition.NONE) {
						// maybe is because without knowing the length it cannot
						// know that is a C-terminal
						if (ptmSite.getPosition() == getSequence().length() + 1) {
							ptmSite.setPTMPosition(PTMPosition.CTERM);
						}
					}
				}

				final boolean ret = ptms.add(newPtm);
				if (ret) {
					fullSequence = null;
				}
				return ret;
			}
		}
		return false;
	}

	@Override
	public Float getExperimentalMH() {
		return experimentalMH;
	}

	@Override
	public Float getCalcMH() {
		return calcMH;
	}

	@Override
	public Float getMassErrorPPM() {
		return massErrorPPM;
	}

	@Override
	public Float getTotalIntensity() {
		return totalIntensity;
	}

	@Override
	public Float getIonProportion() {
		return ionProportion;
	}

	@Override
	public Float getPi() {
		return pi;
	}

	@Override
	public String getFullSequence() {
		if (fullSequence == null || (getPTMs() != null && !getPTMs().isEmpty() && fullSequence != null
				&& fullSequence.equals(sequence))) {
			if (getPTMs() != null && !getPTMs().isEmpty()) {
				fullSequence = ModelUtils.getFullSequenceFromSequenceAndPTMs(getSequence(), getPTMs());

			} else {
				fullSequence = sequence;
			}
		}
		return fullSequence;

	}

	@Override
	public Peptide getPeptide() {

		return peptide;
	}

	@Override
	public boolean setPeptide(Peptide peptide, boolean recursively) {
		if (peptide != null) {
			this.peptide = peptide;
			if (recursively) {
				peptide.addPSM(this, false);
				for (final Protein protein : getProteins()) {
					protein.addPeptide(peptide, false);
					peptide.addProtein(protein, false);
				}
			}
			return true;
		}
		return false;

	}

	@Override
	public int getUniqueIdentifier() {
		return uniqueIdentifier;
	}

	@Override
	public String getAfterSeq() {
		return afterSeq;
	}

	@Override
	public String getBeforeSeq() {
		return beforeSeq;
	}

	@Override
	public Integer getChargeState() {
		return chargeState;
	}

	@Override
	public String getScanNumber() {
		return scanNumber;
	}

	public void setIdentifier(String identifier) {
		this.identifier = identifier;
	}

	public void setScanNumber(String scanNumber) {
		this.scanNumber = scanNumber;
	}

	public void setCalcMH(Float calcMH) {
		this.calcMH = calcMH;
	}

	public void setTotalIntensity(Float intensity) {
		totalIntensity = intensity;
	}

	@Override
	public Integer getSpr() {
		return spr;
	}

	public void setSpr(Integer spr) {
		this.spr = spr;
	}

	public void setFullSequence(String fullSequence) {
		this.fullSequence = FastaParser.getSequenceInBetween(fullSequence);
		sequence = FastaParser.cleanSequence(fullSequence);
		afterSeq = FastaParser.getAfterSeq(fullSequence);
		beforeSeq = FastaParser.getBeforeSeq(fullSequence);
		this.ptms = null;
		this.ptmsInPeptide = null;
		final TIntDoubleHashMap ptMsFromSequence = FastaParser.getPTMsFromSequence(fullSequence, true);
		if (!ptMsFromSequence.isEmpty()) {
			for (final int position : ptMsFromSequence.keys()) {
				final double deltaMass = ptMsFromSequence.get(position);
				String aa = null;
				if (position > 0 && position <= getSequence().length()) { // not N-term or C-term
					aa = getSequence().substring(position - 1, position);
				}
				final PTM ptm = new PTMAdapter(deltaMass, aa, position,
						PTMPosition.getPTMPositionFromSequence(getSequence(), position)).adapt();
				addPTM(ptm);
			}
		}
	}

	public void setChargeState(Integer chargeState) {
		this.chargeState = chargeState;
	}

	public void setMassErrorPPM(Float massErrorPPM) {
		this.massErrorPPM = massErrorPPM;
	}

	public void setPi(Float pi) {
		this.pi = pi;
	}

	public void setIonProportion(Float ionProportion) {
		this.ionProportion = ionProportion;
	}

	@Override
	public Float getRtInMinutes() {
		return rtInMinutes;
	}

	public void setRtInMinutes(Float rtInMinutes) {
		this.rtInMinutes = rtInMinutes;
	}

	public void setExperimentalMH(Float experimentalMH) {
		this.experimentalMH = experimentalMH;
	}

	@Override
	public String toString() {
		StringBuilder sb = null;
		if (getMSRun() != null) {
			sb = new StringBuilder();
			sb.append(" in MSrun: " + getMSRun().getRunId());
		}
		if (sb != null) {
			return getIdentifier() + sb.toString();
		}
		return getIdentifier();
	}

	@Override
	public String getSearchEngine() {
		return searchEngine;
	}

	@Override
	public void setSearchEngine(String searchEngine) {
		this.searchEngine = searchEngine;
	}

	@Override
	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

	@Override
	public Float getDeltaCn() {
		return deltaCn;
	}

	@Override
	public Float getXCorr() {
		return xCorr;
	}

	@Override
	public void setXCorr(Float xcorr2) {
		xCorr = xcorr2;
	}

	@Override
	public void setDeltaCn(Float deltaCn) {
		this.deltaCn = deltaCn;
	}

	@Override
	public final boolean equals(Object obj) {
		if (obj instanceof PSM) {
			final PSM psm = (PSM) obj;
			return psm.getKey().equals(getKey());
		}
		return super.equals(obj);

	}

	@Override
	public final int hashCode() {
		final HashCodeBuilder hashCodeBuilder = new HashCodeBuilder();
		hashCodeBuilder.append(getKey());

		return hashCodeBuilder.toHashCode();
	}

	@Override
	public Set<String> getTaxonomies() {
		if ((taxonomies == null || taxonomies.isEmpty()) && !ignoreTaxonomy) {
			for (final Protein protein : getProteins()) {
				final String fastaHeader = protein.getDescription();
				final String accession = protein.getAccession();
				addTaxonomy(FastaParser.getOrganismNameFromFastaHeader(fastaHeader, accession));
			}
		}
		return taxonomies;
	}

	@Override
	public boolean addTaxonomy(String taxonomy) {
		if (taxonomies == null) {
			taxonomies = new THashSet<String>();
		}
		return taxonomies.add(taxonomy);
	}

	/**
	 * @return the ptms
	 */
	@Override
	public List<PTMInPeptide> getPTMsInPeptide() {
		if (ptmsInPeptide == null) {
			ptmsInPeptide = new ArrayList<PTMInPeptide>();
			final List<StringPosition> tmp = FastaParser.getInside(getFullSequence());
			for (final StringPosition stringPosition : tmp) {
				Double deltaMass = null;
				try {
					deltaMass = Double.valueOf(stringPosition.string);
				} catch (final NumberFormatException e) {

				}
				char aa;
				if (stringPosition.position < 1 || stringPosition.position > getSequence().length()) {
					// it is a N-term or C-term
					aa = ' ';
				} else {
					aa = getSequence().charAt(stringPosition.position - 1);
				}
				final PTMInPeptide ptmInPeptide = new PTMInPeptide(stringPosition.position, aa, getSequence(),
						deltaMass);
				ptmsInPeptide.add(ptmInPeptide);
			}
		}
		return ptmsInPeptide;
	}

	@Override
	public List<PTMInProtein> getPTMsInProtein(UniprotProteinLocalRetrieverInterface uplr,
			Map<String, String> proteinSequences) {
		final List<PTMInProtein> ptmsInProtein = new ArrayList<PTMInProtein>();
		final List<PTMInPeptide> ptms = getPTMsInPeptide();
		if (ptms.isEmpty()) {
			return ptmsInProtein;
		}

		for (final PTMInPeptide ptmInPeptide : ptms) {
			final int positionInPeptide = ptmInPeptide.getPosition();

			final Set<Protein> proteins = getProteins();
			for (final Protein quantifiedProtein : proteins) {
				final String acc = quantifiedProtein.getAccession();
				String proteinSequence = null;
				// it is important that we look for any protein
				// CONTAINING the accession, so that we can search for
				// the isoforms and proteoforms
				if (proteinSequences != null && proteinSequences.containsKey(acc)) {
					proteinSequence = proteinSequences.get(acc);
				}
				if (proteinSequence == null && uplr != null) {
					final Map<String, Entry> annotatedProtein = uplr.getAnnotatedProtein(null, acc);
					final Entry entry = annotatedProtein.get(acc);
					if (entry != null) {
						proteinSequence = UniprotEntryUtil.getProteinSequence(entry);
					}
				}
				if (proteinSequence != null) {
					final String sequence2 = getSequence();
					final TIntArrayList positionsInProteinSequence = StringUtils.allPositionsOf(proteinSequence,
							sequence2);
					if (!positionsInProteinSequence.isEmpty()) {
						for (final int positionInProteinSequence : positionsInProteinSequence.toArray()) {
							final int positionOfSiteInProtein = positionInProteinSequence + positionInPeptide - 1;
							final PTMInProtein ptmInProtein = new PTMInProtein(positionOfSiteInProtein,
									proteinSequence.charAt(positionOfSiteInProtein - 1), acc,
									ptmInPeptide.getDeltaMass());
							if (!ptmsInProtein.contains(ptmInProtein)) {
								ptmsInProtein.add(ptmInProtein);
							}
						}
					} else {
						log.warn("Sequence not found\t" + this.getMSRun().getRunId() + "\tscan:\t"
								+ this.getScanNumber() + "\t" + getFullSequence() + "\t" + getSequence() + "\t" + acc);
					}
				} else {
					throw new IllegalArgumentException("Protein sequence from protein " + acc
							+ " not found neither in the fasta file nor in Uniprot");
				}

			}

		}

		return ptmsInProtein;
	}

	@Override
	public Map<PositionInPeptide, List<PositionInProtein>> getProteinKeysByPeptideKeysForPTMs(
			UniprotProteinLocalRetrieverInterface uplr, Map<String, String> proteinSequences) {

		final Map<PositionInPeptide, List<PositionInProtein>> ret = new THashMap<PositionInPeptide, List<PositionInProtein>>();
		final Set<String> aas = new HashSet<String>();
		aas.clear();

		for (final PTMInPeptide positionInPeptide : getPTMsInPeptide()) {

			final List<PositionInProtein> positionsInProtein = new ArrayList<PositionInProtein>();

			final Set<Protein> proteins = getProteins();
			for (final Protein quantifiedProtein : proteins) {
				final String acc = quantifiedProtein.getAccession();

				String proteinSequence = null;
				// it is important that we look for any protein
				// CONTAINING the accession, so that we can search for
				// the isoforms and proteoforms
				if (proteinSequences != null && proteinSequences.containsKey(acc)) {
					proteinSequence = proteinSequences.get(acc);
				}
				if (proteinSequence == null && uplr != null) {
					final Map<String, Entry> annotatedProtein = uplr.getAnnotatedProtein(null, acc);
					final Entry entry = annotatedProtein.get(acc);
					if (entry != null) {
						proteinSequence = UniprotEntryUtil.getProteinSequence(entry);
					}
				}
				if (proteinSequence != null) {

					final TIntArrayList positionsInProteinSequence = StringUtils.allPositionsOf(proteinSequence,
							getSequence());
					for (final int positionInProteinSequence : positionsInProteinSequence.toArray()) {
						final int positionOfSiteInProtein = positionInProteinSequence + positionInPeptide.getPosition()
								- 1;
						final PositionInProtein positionInProtein = new PositionInProtein(positionOfSiteInProtein,
								proteinSequence.charAt(positionOfSiteInProtein - 1), acc);
						if (!positionsInProtein.contains(positionInProtein)) {
							positionsInProtein.add(positionInProtein);
						}
					}
				} else {
					if (FastaParser.getACC(acc).getAccessionType() == AccessionType.UNKNOWN) {
						log.warn("Protein " + acc + " is ignored because we don't have its protein sequence");
						continue;
					}
					throw new IllegalArgumentException("Protein sequence from protein " + acc
							+ " not found neither in the fasta file nor in Uniprot");
				}

			}
			if (!positionsInProtein.isEmpty()) {
				ret.put(positionInPeptide, positionsInProtein);
			}

		}
		return ret;
	}

	@Override
	public Map<PositionInPeptide, List<PositionInProtein>> getProteinKeysByPeptideKeysForQuantifiedAAs(
			char[] quantifiedAAs, UniprotProteinLocalRetrieverInterface uplr, Map<String, String> proteinSequences) {

		final Map<PositionInPeptide, List<PositionInProtein>> ret = new THashMap<PositionInPeptide, List<PositionInProtein>>();
		final Set<String> aas = new HashSet<String>();
		aas.clear();
		for (final char c : quantifiedAAs) {
			final String aa = String.valueOf(c).toUpperCase();
			if (!aas.contains(aa)) {
				aas.add(aa);
				final TIntArrayList positionsInPeptideSequence = StringUtils.allPositionsOf(getSequence(), aa);
				for (final int positionInPeptide : positionsInPeptideSequence.toArray()) {
					final PositionInPeptide positionInPeptideObj = new PositionInPeptide(positionInPeptide,
							getSequence().charAt(positionInPeptide - 1), getSequence());
					final List<PositionInProtein> positionsInProtein = new ArrayList<PositionInProtein>();

					final Set<Protein> proteins = getProteins();
					for (final Protein quantifiedProtein : proteins) {
						final String acc = quantifiedProtein.getAccession();

						String proteinSequence = null;
						// it is important that we look for any protein
						// CONTAINING the accession, so that we can search for
						// the isoforms and proteoforms
						if (proteinSequences != null && proteinSequences.containsKey(acc)) {
							proteinSequence = proteinSequences.get(acc);
						}
						if (proteinSequence == null && uplr != null) {
							final Map<String, Entry> annotatedProtein = uplr.getAnnotatedProtein(null, acc);
							final Entry entry = annotatedProtein.get(acc);
							if (entry != null) {
								proteinSequence = UniprotEntryUtil.getProteinSequence(entry);
							}
						}
						if (proteinSequence != null) {
							final TIntArrayList positionsInProteinSequence = StringUtils.allPositionsOf(proteinSequence,
									getSequence());
							for (final int positionInProteinSequence : positionsInProteinSequence.toArray()) {
								final int positionOfSiteInProtein = positionInProteinSequence + positionInPeptide - 1;
								final PositionInProtein positionInProtein = new PositionInProtein(
										positionOfSiteInProtein, proteinSequence.charAt(positionOfSiteInProtein - 1),
										acc);
								if (!positionsInProtein.contains(positionInProtein)) {
									positionsInProtein.add(positionInProtein);
								}
							}
						} else {
							if (FastaParser.getACC(acc).getAccessionType() == AccessionType.UNKNOWN) {
								log.warn("Protein " + acc + " is ignored because we don't have its protein sequence");
								continue;
							}
							throw new IllegalArgumentException("Protein sequence from protein " + acc
									+ " not found neither in the fasta file nor in Uniprot");
						}

					}
					if (!positionsInProtein.isEmpty()) {
						ret.put(positionInPeptideObj, positionsInProtein);
					}
				}
			}
		}
		return ret;
	}

	@Override
	public Map<Character, List<PositionInProtein>> getPositionInProteinForSites(char[] quantifiedAAs,
			UniprotProteinLocalRetrieverInterface uplr, Map<String, String> proteinSequences) {
		if (positionsInProteinsByQuantifiedAA == null) {
			positionsInProteinsByQuantifiedAA = new THashMap<Character, List<PositionInProtein>>();
		}

		final Set<Protein> proteins = getProteins();
		for (final Protein quantifiedProtein : proteins) {
			for (final char quantifiedAA : quantifiedAAs) {
				if (positionsInProteinsByQuantifiedAA.containsKey(quantifiedAA)) {
					continue;
				}
				final List<PositionInProtein> list = new ArrayList<PositionInProtein>();
				final String acc = quantifiedProtein.getAccession();

				String proteinSequence = null;
				if (proteinSequences != null && proteinSequences.containsKey(acc)) {
					proteinSequence = proteinSequences.get(acc);
				}
				if (proteinSequence == null && uplr != null) {
					final Map<String, Entry> annotatedProtein = uplr.getAnnotatedProtein(null, acc);
					final Entry entry = annotatedProtein.get(acc);
					if (entry != null) {
						proteinSequence = UniprotEntryUtil.getProteinSequence(entry);
					}
				}

				if (proteinSequence != null) {
					list.addAll(ProteinSequenceUtils.getPositionsInProteinForSites(quantifiedAAs, getSequence(),
							proteinSequence, acc));
					positionsInProteinsByQuantifiedAA.put(quantifiedAA, list);
				}
			}
		}

		return positionsInProteinsByQuantifiedAA;
	}

	@Override
	public List<PositionInProtein> getStartingPositionsInProtein(String proteinACC,
			UniprotProteinLocalRetrieverInterface uplr, Map<String, String> proteinSequences) {
		final List<PositionInProtein> ret = new ArrayList<PositionInProtein>();

		String proteinSequence = null;
		if (proteinSequences != null && proteinSequences.containsKey(proteinACC)) {
			proteinSequence = proteinSequences.get(proteinACC);
		}
		if (proteinSequence == null && uplr != null) {
			final Map<String, Entry> annotatedProtein = uplr.getAnnotatedProtein(null, proteinACC);
			final Entry entry = annotatedProtein.get(proteinACC);
			if (entry != null) {
				proteinSequence = UniprotEntryUtil.getProteinSequence(entry);
			}
		}

		if (proteinSequence != null) {
			ret.addAll(ProteinSequenceUtils.getPositionsOfPeptideSequenceInProteinSequence(getSequence(),
					proteinSequence, proteinACC));
		}

		return ret;
	}

	@Override
	public boolean containsPTMs() {

		return !getPTMsInPeptide().isEmpty();
	}

	@Override
	public boolean isIgnoreTaxonomy() {
		return ignoreTaxonomy;
	}

	@Override
	public void setIgnoreTaxonomy(boolean ignoreTaxonomy) {
		this.ignoreTaxonomy = ignoreTaxonomy;
	}

	@Override
	public final String getKey() {
		if (key == null) {
			key = KeyUtils.getInstance().getSpectrumKey(this, isDistinguishModifiedSequence(), isChargeStateSensible());
		}
		return key;
	}

	public void setKey(String key) {
		this.key = key;
	}

	public boolean isDistinguishModifiedSequence() {
		return distinguishModifiedSequence;
	}

	public boolean isChargeStateSensible() {
		return chargeStateSensible;
	}
}
