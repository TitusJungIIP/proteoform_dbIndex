package edu.scripps.yates.utilities.parsers.idparser;

import java.util.ArrayList;
import java.util.List;

import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.proteomicsmodel.PTM;
import edu.scripps.yates.utilities.proteomicsmodel.PTMPosition;
import edu.scripps.yates.utilities.proteomicsmodel.factories.PTMEx;
import edu.scripps.yates.utilities.staticstorage.StaticStrings;
import gnu.trove.map.hash.TIntDoubleHashMap;

public class PeptideSequence {
	private final String sequence;
	private final char beforeSeq;
	private final char afterSeq;
	private final List<PTM> modifications = new ArrayList<PTM>();
	private final String fullSequence;

	public static final char NULL_SEQ = '-';

	public PeptideSequence(String sequenceToParse, boolean ptmsAsDeltaMass) {
		final String rawSequence = StaticStrings.getUniqueInstance(sequenceToParse);
		fullSequence = StaticStrings.getUniqueInstance(FastaParser.getSequenceInBetween(rawSequence));
		sequence = StaticStrings.getUniqueInstance(FastaParser.cleanSequence(fullSequence));
		final TIntDoubleHashMap ptms = FastaParser.getPTMsFromSequence(fullSequence, ptmsAsDeltaMass);
		for (final int position : ptms.keys()) {
			final double deltaMass = ptms.get(position);
			if (position > 0 && position <= sequence.length()) {
				addModification(deltaMass, position, sequence.substring(position - 1, position));
			} else if (position == 0) {
				final PTMPosition ptmPosition = PTMPosition.NTERM;
				addModification(deltaMass, position, null, ptmPosition);
			} else if (position == sequence.length()) {
				final PTMPosition ptmPosition = PTMPosition.CTERM;
				addModification(deltaMass, position, null, ptmPosition);
			}
		}
		beforeSeq = FastaParser.getBeforeSeq(rawSequence).charAt(0);
		afterSeq = FastaParser.getAfterSeq(rawSequence).charAt(0);
	}

	/**
	 *
	 * @param modificationShift
	 * @param modPosition
	 *            position of the modification, starting by 1 at the first AA
	 * @param aa
	 */
	private void addModification(Double modificationShift, int modPosition, String aa) {
		addModification(modificationShift, modPosition, aa, PTMPosition.NONE);
	}

	private void addModification(Double modificationShift, int modPosition, String aa, PTMPosition ptmPosition) {
		modifications.add(new PTMEx(modificationShift, aa, modPosition, ptmPosition));
	}

	/**
	 * @return the sequence with no modifications
	 */
	public String getSequence() {
		return sequence.toString();
	}

	/**
	 * @return the beforeSeq
	 */
	public char getBeforeSeq() {
		return beforeSeq;
	}

	/**
	 * @return the afterSeq
	 */
	public char getAfterSeq() {
		return afterSeq;
	}

	/**
	 * @return the modifications
	 */
	public List<PTM> getModifications() {
		return modifications;
	}

	public String getFullSequence() {
		return fullSequence;
	}

}
