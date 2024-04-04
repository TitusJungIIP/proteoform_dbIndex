package edu.scripps.yates.utilities.sequence;

import java.util.ArrayList;
import java.util.List;

import edu.scripps.yates.utilities.strings.StringUtils;
import gnu.trove.list.array.TIntArrayList;

public class ProteinSequenceUtils {

	/**
	 * Returns the positions in which the peptide sequence is found in the protein
	 * sequence.<br>
	 * This functions is equals to StringUtils.allPositionsOf(peptideSequence,
	 * proteinSequence);
	 * 
	 * @param peptideSequence
	 * @param proteinSequence
	 * @param proteinACC
	 * @return
	 */
	public static List<PositionInProtein> getPositionsOfPeptideSequenceInProteinSequence(String peptideSequence,
			String proteinSequence, String proteinACC) {
		final List<PositionInProtein> ret = new ArrayList<PositionInProtein>();

		final String proteinSequenceSafe = proteinSequence.toUpperCase().replace("I", "&").replace("L", "&");
		final String peptideSequenceSafe = peptideSequence.toUpperCase().replace("I", "&").replace("L", "&");
		// because in proteomics we cannot distinguish between I and L, we convert all
		// to & in the protein and in the peptide

		final TIntArrayList allPositionsOf = StringUtils.allPositionsOf(proteinSequenceSafe, peptideSequenceSafe);
		allPositionsOf.forEach(
				position -> ret.add(new PositionInProtein(position, proteinSequence.charAt(position - 1), proteinACC)));

		return ret;
	}

	/**
	 * Returns the positions in which the aminoacids in aas contained in the peptide
	 * sequence are in the protein sequence
	 * 
	 * @param aas
	 * @param peptideSequence
	 * @param proteinSequence
	 * @param proteinACC
	 * @return
	 */
	public static List<PositionInProtein> getPositionsInProteinForSites(char[] aas, String peptideSequence,
			String proteinSequence, String proteinACC) {
		final List<PositionInProtein> ret = new ArrayList<PositionInProtein>();

		// because in proteomics we cannot distinguish between I and L, we convert all
		// to & in the protein and in the peptide
		final String proteinSequenceSafe = proteinSequence.toUpperCase().replace("I", "&").replace("L", "&");
		final String peptideSequenceSafe = peptideSequence.toUpperCase().replace("I", "&").replace("L", "&");
		for (char aa : aas) {
			aa = Character.toUpperCase(aa);
			if (aa == 'I' || aa == 'L') {
				aa = '&';
			}
			final TIntArrayList positionsInPeptide = StringUtils.allPositionsOf(peptideSequenceSafe, aa);
			for (final int positionInPeptide : positionsInPeptide.toArray()) {

				final TIntArrayList positionsInProtein = StringUtils.allPositionsOf(proteinSequenceSafe,
						peptideSequenceSafe);
				for (final int positionInProtein : positionsInProtein.toArray()) {
					final int positionOfSiteInProtein = positionInProtein + positionInPeptide - 1;
					ret.add(new PositionInProtein(positionOfSiteInProtein,
							proteinSequence.charAt(positionOfSiteInProtein - 1), proteinACC));
				}
			}
		}
		return ret;
	}
}
