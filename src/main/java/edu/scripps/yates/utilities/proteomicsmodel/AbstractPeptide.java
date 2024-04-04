package edu.scripps.yates.utilities.proteomicsmodel;

import java.util.ArrayList;
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

public abstract class AbstractPeptide implements Peptide {
	private static Logger log = Logger.getLogger(AbstractPeptide.class);
	private Set<Score> scores;
	private Set<Ratio> ratios;
	private Set<Amount> amounts;
	private Set<Condition> conditions;
	private Set<MSRun> msRuns;
	private String sequence;
	private PeptideRelation relation;
	private Set<Protein> proteins;
	private List<PSM> psms;
	private List<PTM> ptms;

	private String fullSequence;
	private String afterSeq;
	private String beforeSeq;
	private String searchEngine;
	private Integer spectrumCount;
	private int dbID;
	private List<PTMInPeptide> ptmsInPeptide;
	private Map<Character, List<PositionInProtein>> positionsInProteinsByQuantifiedAA;
	private Set<String> taxonomies;
	private boolean ignoreTaxonomy;
	private final String key;
	private boolean distinguishModifiedSequence;
	private boolean chargeSensible;

	public AbstractPeptide(String key) {
		this.key = key;
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
	public Set<MSRun> getMSRuns() {
		if (msRuns == null) {
			msRuns = new THashSet<MSRun>();
		}
		return msRuns;
	}

	@Override
	public boolean addMSRun(MSRun msRun) {
		if (msRun != null) {
			if (msRuns == null) {
				msRuns = new THashSet<MSRun>();
			}
			return msRuns.add(msRun);
		}
		return false;
	}

	@Override
	public String getSequence() {
		return sequence;
	}

	@Override
	public String getIdentifier() {
		return getFullSequence();
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
			if (recursively) {
				protein.addPeptide(this, false);
				if (getPSMs() != null) {
					for (final PSM psm : getPSMs()) {
						protein.addPSM(psm, false);
						psm.addProtein(protein, false);
					}
				}
			}
			return ret;
		}
		return false;
	}

	@Override
	public List<PTM> getPTMs() {
		return ptms;
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
	public List<PSM> getPSMs() {
		if (psms == null) {
			psms = new ArrayList<PSM>();
		}
		return psms;
	}

	@Override
	public boolean addPSM(PSM psm, boolean recursively) {
		if (psm != null) {
			if (psms == null) {
				psms = new ArrayList<PSM>();
			}
			if (!psms.contains(psm)) {
				final boolean ret = psms.add(psm);
				if (psm.getMSRun() != null) {
					addMSRun(psm.getMSRun());
				}
				if (recursively) {
					psm.setPeptide(this, false);
					for (final Protein protein : psm.getProteins()) {
						addProtein(protein, false);
					}
					for (final Protein protein : getProteins()) {
						psm.addProtein(protein, false);
					}
				}
				return ret;
			}
		}
		return false;
	}

	@Override
	public String getAfterSeq() {
		return afterSeq;
	}

	@Override
	public String getBeforeSeq() {
		return beforeSeq;
	}

	public void setFullSequence(String fullSequence) {
		this.fullSequence = fullSequence;
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

	@Override
	public String toString() {
		final StringBuilder sb = new StringBuilder();
		if (getMSRuns() != null) {
			sb.append(" in run(s): ");
			for (final MSRun msRun : getMSRuns()) {
				if (!"".equals(sb.toString())) {
					sb.append(",");
				}
				sb.append(msRun.getRunId());
			}
		}
		return getIdentifier() + sb.toString();
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
	public Integer getSpectrumCount() {
		return spectrumCount;
	}

	@Override
	public void setSpectrumCount(Integer spc) {
		spectrumCount = spc;
	}

	@Override
	public int getDBId() {
		return dbID;
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
					// add it to the PSMs too
					for (final PSM psm : getPSMs()) {
						psm.addPTM(newPtm);
					}
					fullSequence = null;
				}
				return ret;
			}
		}
		return false;
	}

	@Override
	public void setSequence(String sequence2) {
		sequence = sequence2;
	}

	@Override
	public final int hashCode() {
		return new HashCodeBuilder().append(getKey()).toHashCode();
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof Peptide) {
			return ((Peptide) obj).getKey().equals(getKey());
		}
		return super.equals(obj);
	}

	@Override
	public void mergeWithPeptide(Peptide otherPeptide) {
		for (final PSM psm : otherPeptide.getPSMs()) {
			addPSM(psm, true);
		}
		for (final Protein protein : otherPeptide.getProteins()) {
			addProtein(protein, true);
		}
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
					final TIntArrayList positionsInProteinSequence = StringUtils.allPositionsOf(proteinSequence,
							getSequence());
					for (final int positionInProteinSequence : positionsInProteinSequence.toArray()) {
						final int positionOfSiteInProtein = positionInProteinSequence + positionInPeptide - 1;
						final PTMInProtein ptmInProtein = new PTMInProtein(positionOfSiteInProtein,
								proteinSequence.charAt(positionOfSiteInProtein - 1), acc, ptmInPeptide.getDeltaMass());
						if (!ptmsInProtein.contains(ptmInProtein)) {
							ptmsInProtein.add(ptmInProtein);
						}
					}
				} else {
					throw new IllegalArgumentException("Protein sequence from protein " + acc
							+ " not found neither in the fasta file nor in Uniprot");
				}

			}

		}

		return ptmsInProtein;
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
				final int position = stringPosition.position;
				char aa = 0;
				// it can be 0 when it is n-terminal or length +1 when it is on c-term
				if (position > 0 && position <= getSequence().length()) {
					aa = getSequence().charAt(position - 1);
				} else {
//					log.debug("Modification at position " + position + ": " + getFullSequence());
				}
				Double deltaMass = null;
				try {
					deltaMass = Double.valueOf(stringPosition.string);
				} catch (final NumberFormatException e) {

				}
				final PTMInPeptide ptm = new PTMInPeptide(position, aa, getSequence(), deltaMass);
				ptmsInPeptide.add(ptm);
			}
		}

		return ptmsInPeptide;
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
			UniprotProteinLocalRetrieverInterface uplr) {
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
				final Map<String, Entry> annotatedProtein = uplr.getAnnotatedProtein(null, acc);
				final Entry entry = annotatedProtein.get(acc);
				if (entry != null) {
					final String proteinSequence = UniprotEntryUtil.getProteinSequence(entry);
					if (proteinSequence != null) {
						list.addAll(ProteinSequenceUtils.getPositionsInProteinForSites(quantifiedAAs, getSequence(),
								proteinSequence, acc));
						positionsInProteinsByQuantifiedAA.put(quantifiedAA, list);
					}
				}
			}
		}

		return positionsInProteinsByQuantifiedAA;
	}

	@Override
	public List<PositionInProtein> getStartingPositionsInProtein(String proteinACC,
			UniprotProteinLocalRetrieverInterface uplr) {
		final List<PositionInProtein> ret = new ArrayList<PositionInProtein>();

		final Map<String, Entry> annotatedProtein = uplr.getAnnotatedProtein(null, proteinACC);
		final Entry entry = annotatedProtein.get(proteinACC);
		if (entry != null) {
			final String proteinSequence = UniprotEntryUtil.getProteinSequence(entry);
			if (proteinSequence != null) {
				ret.addAll(ProteinSequenceUtils.getPositionsOfPeptideSequenceInProteinSequence(getSequence(),
						proteinSequence, proteinACC));
			}
		}

		return ret;
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
		return key;
	}

	public void setDistinguishModifiedSequence(boolean distinguishModifiedSequence) {
		this.distinguishModifiedSequence = distinguishModifiedSequence;
	}

	public void setChargeSensible(boolean chargeSensible) {
		this.chargeSensible = chargeSensible;
	}

	public boolean isDistinguishModifiedSequence() {
		return distinguishModifiedSequence;
	}

	public boolean isChargeSensible() {
		return chargeSensible;
	}

}
