package edu.scripps.yates.utilities.proteomicsmodel;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.grouping.GroupablePeptide;
import edu.scripps.yates.utilities.grouping.ProteinEvidence;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AccessionType;
import edu.scripps.yates.utilities.proteomicsmodel.factories.AccessionEx;
import edu.scripps.yates.utilities.proteomicsmodel.factories.GeneEx;
import edu.scripps.yates.utilities.proteomicsmodel.utils.ModelUtils;
import edu.scripps.yates.utilities.proteomicsmodel.utils.PSMsOfAProtein;
import gnu.trove.set.hash.THashSet;

public abstract class AbstractProtein implements Protein {
	private static final Logger log = Logger.getLogger(AbstractProtein.class);

	private Set<Score> scores;
	private Set<Ratio> ratios;
	private Set<Amount> amounts;
	private Set<Condition> conditions;
	private PSMsOfAProtein psms;
	private Integer spectrumCount;
	private Set<MSRun> msRuns;
	private ProteinGroup proteinGroup;
	private ProteinEvidence evidence;
	private Accession primaryAccession;
	private Set<Accession> secondaryAccessions;
	private Set<Gene> genes;
	private Set<ProteinAnnotation> annotations;
	private Set<Threshold> thresholds;
	private THashSet<Peptide> peptides;
	private Float coverage;
	private Float empai;
	private Float mw;
	private Float nsaf_norm;
	private Float nsaf;
	private String searchEngine;
	private Float pi;
	private Integer length;
	private Organism organism;
	private String sequence;
	private Set<String> taxonomies;
	private boolean ignoreTaxonomy;
	private boolean genesParsed = false;
	private boolean organismParsed;
	private String key;

	@Override
	public Set<Score> getScores() {
		return scores;
	}

	@Override
	public boolean addScore(Score score) {
		if (score != null) {
			if (scores == null) {
				scores = new THashSet<Score>();
			}
			return scores.add(score);
		}
		return false;
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
		if (ratio != null) {
			if (ratios == null) {
				ratios = new THashSet<Ratio>();
			}
			return ratios.add(ratio);
		}
		return false;
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
		if (amount != null) {
			if (amounts == null) {
				amounts = new THashSet<Amount>();
			}
			return amounts.add(amount);
		}
		return false;
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
		if (condition != null) {
			if (conditions == null) {
				conditions = new THashSet<Condition>();
			}
			if (condition.getSample() != null && condition.getSample().getOrganism() != null) {
				setOrganism(condition.getSample().getOrganism());
			}
			return conditions.add(condition);
		}
		return false;
	}

	@Override
	public List<PSM> getPSMs() {
		if (psms == null) {
			psms = new PSMsOfAProtein(this);
		}
		return psms;
	}

	@Override
	public boolean addPSM(PSM psm, boolean recursively) {
		if (psm != null) {
			if (psms == null) {
				psms = new PSMsOfAProtein(this);
			}
			if (!psms.contains(psm)) {
				final boolean ret = psms.add(psm);
				if (psm.getMSRun() != null) {
					addMSRun(psm.getMSRun());
				}
				if (recursively) {
					addPeptide(psm.getPeptide(), false);
					if (psm.getPeptide() != null) {
						psm.getPeptide().addProtein(this, false);
					}
					psm.addProtein(this, false);
				}
				return ret;
			}
		}
		return false;
	}

	@Override
	public Integer getSpectrumCount() {
		if (spectrumCount != null) {
			return spectrumCount;
		}
		return getPSMs().size();
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
	public List<GroupablePeptide> getGroupablePeptides() {
		final List<GroupablePeptide> ret = new ArrayList<GroupablePeptide>();
		for (final GroupablePeptide groupablePeptide : getPSMs()) {
			ret.add(groupablePeptide);
		}
		return ret;
	}

	@Override
	public ProteinGroup getProteinGroup() {
		return proteinGroup;
	}

	@Override
	public String getUniqueID() {
		return getKey();
	}

	@Override
	public String getAccession() {
		return getPrimaryAccession().getAccession();
	}

	@Override
	public void setEvidence(ProteinEvidence evidence) {
		this.evidence = evidence;
	}

	@Override
	public ProteinEvidence getEvidence() {
		return evidence;
	}

	@Override
	public void setProteinGroup(ProteinGroup proteinGroup) {
		this.proteinGroup = proteinGroup;
	}

	@Override
	public Accession getPrimaryAccession() {
		return primaryAccession;
	}

	@Override
	public Set<Accession> getSecondaryAccessions() {
		return secondaryAccessions;
	}

	@Override
	public Set<Gene> getGenes() {
		if ((genes == null || genes.isEmpty()) && !genesParsed) {
			final String geneFromFastaHeader = FastaParser.getGeneFromFastaHeader(getDescription());
			if (geneFromFastaHeader != null) {
				final GeneEx gene = new GeneEx(geneFromFastaHeader);
				addGene(gene);
			}

			genesParsed = true;
		}
		return genes;
	}

	@Override
	public Set<ProteinAnnotation> getAnnotations() {
		if (annotations == null) {
			annotations = new THashSet<ProteinAnnotation>();
		}
		return annotations;
	}

	@Override
	public Set<Threshold> getThresholds() {
		if (thresholds == null) {
			thresholds = new THashSet<Threshold>();
		}
		return thresholds;
	}

	@Override
	public boolean addThreshold(Threshold threshold) {
		if (threshold != null) {
			if (thresholds == null) {
				thresholds = new THashSet<Threshold>();
			}
			return thresholds.add(threshold);
		}
		return false;
	}

	@Override
	public Boolean passThreshold(String thresholdName) {
		if (thresholds != null) {
			for (final Threshold threshold : thresholds) {
				if (threshold.getName().equalsIgnoreCase(thresholdName)) {
					return threshold.isPassThreshold();
				}
			}
		}
		return null;
	}

	@Override
	public Set<Peptide> getPeptides() {
		if (peptides == null) {
			peptides = new THashSet<Peptide>();
		}
		return peptides;
	}

	@Override
	public Integer getLength() {
		return length;
	}

	@Override
	public Float getPi() {
		return pi;
	}

	@Override
	public String getSequence() {
		return sequence;
	}

	@Override
	public Organism getOrganism() {
		if (!organismParsed && organism == null) {
			if (FastaParser.isContaminant(getAccession())) {
				organism = ModelUtils.getOrganismContaminant();
			}
			organismParsed = true;
		}
		return organism;
	}

	@Override
	public void setOrganism(Organism organism) {
		this.organism = organism;
	}

	@Override
	public void setMw(Float mw) {
		this.mw = mw;
	}

	@Override
	public void setPi(Float pi) {
		this.pi = pi;
	}

	@Override
	public void setLength(Integer length) {
		this.length = length;
	}

	@Override
	public boolean addPeptide(Peptide peptide, boolean recursively) {
		if (peptide != null) {
			if (peptides == null) {
				peptides = new THashSet<Peptide>();
			}
			final boolean ret = peptides.add(peptide);
			if (recursively) {
				for (final PSM psm : peptide.getPSMs()) {
					addPSM(psm, false);
					psm.addProtein(this, false);
				}
				peptide.addProtein(this, false);
			}
			return ret;
		}
		return false;
	}

	@Override
	public boolean addGene(Gene gene) {
		if (gene != null) {
			if (genes == null) {
				genes = new THashSet<Gene>();
			}
			return genes.add(gene);
		}
		return false;
	}

	@Override
	public void mergeWithProtein(Protein protein) {
		if (protein == this) {
			return;
		}
		if (protein.getCoverage() != null) {
			setCoverage(protein.getCoverage());
		}
		if (getPrimaryAccession() == null) {
			setPrimaryAccession(protein.getPrimaryAccession());
		}
		final String description = protein.getDescription();
		if (protein.getDescription() != null) {
			setDescription(description);
		}
		if (protein.getEmpai() != null) {
			setEmpai(protein.getEmpai());
		}
		if (protein.getLength() != null) {
			setLength(protein.getLength());
		}

		if (getMw() == null) {
			setMw(protein.getMw());
		}
		if (protein.getNsaf() != null) {
			setNsaf(protein.getNsaf());
		}
		if (getNsaf_norm() == null) {
			setNsaf_norm(protein.getNsaf_norm());
		}
		if (getPi() == null) {
			setPi(protein.getPi());
		}
		// protein may have to be grouped again since it contains new PSMs, so
		// we reset group
		setProteinGroup(null);
		// spectrumCount may not be accurate now with new psms
		setSpectrumCount(null);

		if (protein.getPSMs() != null) {
			for (final PSM psm : protein.getPSMs()) {
				addPSM(psm, true);
			}
		}
		if (protein.getPeptides() != null) {
			for (final Peptide peptide2 : protein.getPeptides()) {
				boolean found = false;
				for (final Peptide peptide : getPeptides()) {
					if (peptide.equals(peptide2)) {
						found = true;
					}
				}
				if (!found) {
					addPeptide(peptide2, true);
				}
			}
		}

		if (getSearchEngine() == null) {
			setSearchEngine(protein.getSearchEngine());
		}
	}

	@Override
	public Float getCoverage() {
		return coverage;
	}

	@Override
	public Float getEmpai() {
		return empai;
	}

	@Override
	public String getDescription() {
		return getPrimaryAccession().getDescription();
	}

	@Override
	public void setDescription(String description) {
		final Accession primaryAccession2 = getPrimaryAccession();
		if (primaryAccession2 != null) {
			primaryAccession2.setDescription(description);
		}
	}

	@Override
	public Float getMw() {
		return mw;
	}

	@Override
	public Float getNsaf_norm() {
		return nsaf_norm;
	}

	@Override
	public Float getNsaf() {
		return nsaf;
	}

	@Override
	public String getSearchEngine() {
		return searchEngine;
	}

	@Override
	public void setCoverage(Float coverage) {
		this.coverage = coverage;
	}

	@Override
	public void setNsaf_norm(Float nsaf_norm) {
		this.nsaf_norm = nsaf_norm;
	}

	@Override
	public void setNsaf(Float nsaf) {
		this.nsaf = nsaf;
	}

	@Override
	public void setEmpai(Float empai) {
		this.empai = empai;
	}

	@Override
	public void setPrimaryAccession(Accession accession) {
		primaryAccession = accession;
	}

	@Override
	public void setPrimaryAccession(String accession) {
		setPrimaryAccession(FastaParser.getACC(accession));
	}

	@Override
	public void setSearchEngine(String searchEngine) {
		this.searchEngine = searchEngine;
	}

	@Override
	public void setSpectrumCount(Integer spc) {
		spectrumCount = spc;
	}

	@Override
	public void setPrimaryAccession(AccessionType accessionType, String accession) {
		final AccessionEx newAccession = new AccessionEx(accession, accessionType);
		setPrimaryAccession(newAccession);
	}

	@Override
	public boolean addSecondaryAccession(AccessionType accessionType, String accession) {

		if (accessionType != null && accession != null) {
			final AccessionEx newAccession = new AccessionEx(accession, accessionType);
			return addSecondaryAccession(newAccession);
		}
		return false;
	}

	@Override
	public boolean addSecondaryAccession(Accession accession) {
		if (accession != null) {
			if (secondaryAccessions == null) {
				secondaryAccessions = new THashSet<Accession>();
			}

			return secondaryAccessions.add(accession);
		}
		return false;
	}

	@Override
	public boolean addProteinAnnotation(ProteinAnnotation proteinAnnotation) {

		if (proteinAnnotation != null) {
			if (annotations == null) {
				annotations = new THashSet<ProteinAnnotation>();
			}
			return annotations.add(proteinAnnotation);
		}
		return false;
	}

	@Override
	public void addProteinAnnotations(Collection<ProteinAnnotation> proteinAnnotations) {
		for (final ProteinAnnotation annotation : proteinAnnotations) {
			addProteinAnnotation(annotation);
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
		return getAccession() + sb.toString();
	}

	@Override
	public Set<String> getTaxonomies() {
		if ((taxonomies == null || taxonomies.isEmpty()) && !ignoreTaxonomy) {
			final String fastaHeader = getDescription();
			final String accession = getAccession();
			addTaxonomy(FastaParser.getOrganismNameFromFastaHeader(fastaHeader, accession));
		}
		return taxonomies;
	}

	@Override
	public boolean addTaxonomy(String taxonomy) {

		if (taxonomies == null) {
			taxonomies = new THashSet<String>();
		}
		if (taxonomy != null) {
			return taxonomies.add(taxonomy);
		}
		return false;
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
	public final int hashCode() {
		return new HashCodeBuilder().append(getKey()).toHashCode();
	}

	@Override
	public final boolean equals(Object obj) {
		if (obj instanceof Protein) {
			return ((Protein) obj).getKey().equals(getKey());
		}
		return super.equals(obj);
	}

	@Override
	public final String getKey() {
		return key;
	}

	protected void setKey(String key) {
		this.key = key;
	}

}
