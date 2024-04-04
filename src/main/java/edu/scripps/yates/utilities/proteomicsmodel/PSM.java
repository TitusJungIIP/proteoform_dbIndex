package edu.scripps.yates.utilities.proteomicsmodel;

import java.util.List;
import java.util.Map;

import edu.scripps.yates.utilities.annotations.UniprotProteinLocalRetrieverInterface;
import edu.scripps.yates.utilities.grouping.GroupablePeptide;
import edu.scripps.yates.utilities.sequence.PTMInPeptide;
import edu.scripps.yates.utilities.sequence.PTMInProtein;
import edu.scripps.yates.utilities.sequence.PositionInPeptide;
import edu.scripps.yates.utilities.sequence.PositionInProtein;

public interface PSM extends HasScores, HasRatios, HasAmounts, HasConditions, HasMSRun, GroupablePeptide, HasProteins,
		HasPTMs, HasPeptide, HasTaxonomies, HasKey {

	public Float getExperimentalMH();

	public Float getCalcMH();

	public Float getMassErrorPPM();

	public Float getTotalIntensity();

	public Integer getSpr();

	public Float getIonProportion();

	public Float getPi();

	/**
	 * Gets the peptide sequence including PTMs Nterm and Cterm aminoacids, as it is
	 * reported by the search engine
	 *
	 * @return
	 */
	public String getFullSequence();

	public int getUniqueIdentifier();

	public String getAfterSeq();

	public String getBeforeSeq();

	public Integer getChargeState();

	public String getScanNumber();

	public void setSearchEngine(String searchEngine);

	public String getSearchEngine();

	public Float getDeltaCn();

	public void setDeltaCn(Float deltaCn);

	public Float getXCorr();

	public void setXCorr(Float xcorr);

	public Float getRtInMinutes();

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

	public boolean containsPTMs();

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
			UniprotProteinLocalRetrieverInterface uplr, Map<String, String> proteinSequences);

	public List<PositionInProtein> getStartingPositionsInProtein(String proteinACC,
			UniprotProteinLocalRetrieverInterface uplr, Map<String, String> proteinSequences);

	public void setSequence(String seq);
}
