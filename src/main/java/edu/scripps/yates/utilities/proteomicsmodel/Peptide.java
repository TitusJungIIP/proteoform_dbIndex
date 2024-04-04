package edu.scripps.yates.utilities.proteomicsmodel;

import java.util.List;
import java.util.Map;

import edu.scripps.yates.utilities.annotations.UniprotProteinLocalRetrieverInterface;
import edu.scripps.yates.utilities.grouping.GroupablePeptide;
import edu.scripps.yates.utilities.sequence.PTMInPeptide;
import edu.scripps.yates.utilities.sequence.PTMInProtein;
import edu.scripps.yates.utilities.sequence.PositionInPeptide;
import edu.scripps.yates.utilities.sequence.PositionInProtein;

public interface Peptide extends HasRatios, HasScores, HasAmounts, HasConditions, HasPsms, HasMsRuns, HasProteins,
		HasPTMs, GroupablePeptide, HasTaxonomies, HasKey {

	public int getDBId();

	public String getFullSequence();

	public String getAfterSeq();

	public String getBeforeSeq();

	public void setSearchEngine(String searchEngine);

	public String getSearchEngine();

	public void setSequence(String sequence);

	public void mergeWithPeptide(Peptide otherPeptide);

	public List<PTMInPeptide> getPTMsInPeptide();

	/**
	 * Get a list of {@link PositionInProtein} for each quantified site in the
	 * peptide sequence (represented as a {@link PositionInPeptide}).<br>
	 * Examples:<br>
	 * "PEPTIDE#4 {PROTEIN1#234, PROTEIN2#123}
	 * 
	 * @param quantifiedAAs
	 * @param uplr             used in order to get the protein sequence
	 * @param proteinSequences map of protein sequences
	 * @param dBIndex
	 * @return
	 */
	public Map<PositionInPeptide, List<PositionInProtein>> getProteinKeysByPeptideKeysForQuantifiedAAs(
			char[] quantifiedAAs, UniprotProteinLocalRetrieverInterface uplr, Map<String, String> proteinSequences);

	/**
	 * Get a list of {@link PositionInProtein} for each PTM site in the peptide
	 * sequence (represented as a {@link PositionInPeptide}).<br>
	 * Examples:<br>
	 * "PEPTIDE#4 {PROTEIN1#234, PROTEIN2#123}
	 * 
	 * @param quantifiedAAs
	 * @param uplr             used in order to get the protein sequence
	 * @param proteinSequences map of protein sequences
	 * @param dBIndex
	 * @return
	 */
	public Map<PositionInPeptide, List<PositionInProtein>> getProteinKeysByPeptideKeysForPTMs(
			UniprotProteinLocalRetrieverInterface uplr, Map<String, String> proteinSequences);

	public Map<Character, List<PositionInProtein>> getPositionInProteinForSites(char[] quantifiedAAs,
			UniprotProteinLocalRetrieverInterface uplr);

	public List<PositionInProtein> getStartingPositionsInProtein(String proteinACC,
			UniprotProteinLocalRetrieverInterface uplr);

	public boolean containsPTMs();

	/**
	 * Get a list of {@link PositionInProtein} for each ptm site in the peptide
	 * sequence (represented as a {@link PositionInPeptide}).<br>
	 * Examples:<br>
	 * "PEPTIDE#4 {PROTEIN1#234#238, PROTEIN2#123#127}
	 * 
	 * @param uplr             used in order to get the protein sequence
	 * @param proteinSequences map of protein sequences
	 * @param dBIndex
	 * @return
	 */
	public List<PTMInProtein> getPTMsInProtein(UniprotProteinLocalRetrieverInterface uplr,
			Map<String, String> proteinSequences);
}
