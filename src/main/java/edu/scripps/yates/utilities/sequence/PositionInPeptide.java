package edu.scripps.yates.utilities.sequence;

import java.util.ArrayList;
import java.util.List;

/**
 * Same class as {@link PositionInProtein}, but just to distinguish them.<br>
 * The difference is that the key here is the peptide sequence and there is a
 * validation in the constructor that checks if the position is under the length
 * of the sequence
 * 
 * @author Salva
 *
 */
public class PositionInPeptide extends PositionInProtein {

	public PositionInPeptide(int position, char aa, String peptideSequence) {
		super(position, aa, peptideSequence);
		// we allow position 0 and length+1 to refer to n-term and c-term respectivelly
		if (position < 0 || position > peptideSequence.length() + 1) {
			throw new IllegalArgumentException(
					"Position " + position + " cannot belong to peptide sequence " + peptideSequence);
		}
	}

	public static List<PositionInPeptide> parseStringToPositionInPeptide(String string, String separator) {
		final List<PositionInPeptide> ret = new ArrayList<PositionInPeptide>();
		final List<String> subStrings = new ArrayList<String>();
		if (string.contains(separator)) {
			final String[] split = string.split(separator);
			for (final String string2 : split) {
				subStrings.add(string2);
			}
		} else {
			subStrings.add(string);
		}
		for (final String string2 : subStrings) {
			if (string2.contains(SEPARATOR)) {
				try {
					final String[] split = string2.split(SEPARATOR);
					final PositionInPeptide positionInProtein = new PositionInPeptide(
							Integer.valueOf(String.valueOf(split[1].substring(1))), split[1].charAt(0), split[0]);
					ret.add(positionInProtein);
				} catch (final NumberFormatException e) {

				}
			}
		}
		return ret;
	}
}
