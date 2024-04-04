package edu.scripps.yates.utilities.proteomicsmodel;

import java.util.Collection;
import java.util.Set;

import edu.scripps.yates.utilities.grouping.GroupableProtein;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AccessionType;

/**
 * This class represents a Protein, that has been detected in N {@link MSRun}s
 * and be measured with an ({@link Amount}) in a Experimental condition
 * ({@link Condition}) .<br>
 * {@link Protein} may belongs to a different conditions ({@link Condition}).
 *
 * @author Salva
 *
 */
public interface Protein
		extends HasScores, HasRatios, HasAmounts, HasConditions, HasPsms, HasPeptides, HasMsRuns, GroupableProtein,
		// a protein has several taxonomies because sometimes it represents a
		// set of proteins
		HasTaxonomies, HasKey {

	public Accession getPrimaryAccession();

	public Set<Accession> getSecondaryAccessions();

	/**
	 * @return the genes
	 */
	public Set<Gene> getGenes();

	/**
	 * @return the proteinAnnotations
	 */
	public Set<ProteinAnnotation> getAnnotations();

	/**
	 * Gets all the {@link Threshold}s that has been applied to the protein.
	 *
	 * @return
	 */
	public Set<Threshold> getThresholds();

	/**
	 * If a {@link Threshold} with the name indicated in the parameter has been
	 * applied to the protein it will return yes if it has been passed or false
	 * otherwise. <br>
	 * In case of the {@link Threshold} has not been applied to the {@link Protein},
	 * it will return null value.
	 *
	 * @param thresholdName
	 * @return
	 */
	public Boolean passThreshold(String thresholdName);

	public Integer getLength();

	public Float getPi();

	public String getSequence();

	public Organism getOrganism();

	public void setOrganism(Organism organism);

	public void setLength(Integer length);

	public boolean addGene(Gene gene);

	public void mergeWithProtein(Protein otherProtein);

	public Float getCoverage();

	public Float getEmpai();

	public String getDescription();

	public void setDescription(String description);

	public Float getMw();

	public Float getNsaf_norm();

	public Float getNsaf();

	public String getSearchEngine();

	public boolean addThreshold(Threshold threshold);

	public void setEmpai(Float empai);

	public void setNsaf(Float nsaf);

	public void setNsaf_norm(Float nsaf_norm);

	public void setCoverage(Float coverage);

	public void setMw(Float mw);

	public void setPi(Float pi);

	public void setPrimaryAccession(Accession accession);

	public void setPrimaryAccession(String accession);

	public void setSearchEngine(String searchEngine);

	public void setPrimaryAccession(AccessionType accessionType, String accession);

	public boolean addSecondaryAccession(AccessionType accessionType, String accession);

	public boolean addSecondaryAccession(Accession accession);

	public boolean addProteinAnnotation(ProteinAnnotation proteinAnnotation);

	public void addProteinAnnotations(Collection<ProteinAnnotation> proteinAnnotations);

}
