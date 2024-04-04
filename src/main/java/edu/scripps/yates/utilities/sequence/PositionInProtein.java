package edu.scripps.yates.utilities.sequence;

import java.util.ArrayList;
import java.util.List;

/**
 * Class that represents a position in a protein sequence, where we use the
 * protein accession as key
 * 
 * @author Salva
 *
 */
public class PositionInProtein {
	private final int position;
	private final String proteinACC;
	public final static char NULL_CHAR = '0';
	private char aa = NULL_CHAR;
	public final static String SEPARATOR = "#";

	public PositionInProtein(int position, char aa, String proteinACC) {
		this.position = position;
		this.proteinACC = proteinACC;
		this.aa = aa;
	}

	public int getPosition() {
		return position;
	}

	public String getProteinACC() {
		return proteinACC;
	}

	@Override
	public String toString() {
		return proteinACC + SEPARATOR + aa + position;
	}

	@Override
	public int hashCode() {
		int hash = 23;
		hash = hash * 31 + position;
		hash = hash * 31 + proteinACC.hashCode();
		return hash;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof PositionInProtein) {
			final PositionInProtein positionInProtein = (PositionInProtein) obj;
			if (positionInProtein.position == position && positionInProtein.proteinACC.equals(proteinACC)) {
				return true;
			}
			return false;
		}
		return super.equals(obj);
	}

	private static boolean isNumber(char c) {
		return c >= 48 && c <= 57;
	}

	/**
	 * Converts the string "Q9Y3B7#K81-Q9Y3B7-2#K55" to position 81 of Q9Y3B7 and
	 * position 55 of Q9Y3B7
	 * 
	 * @param string
	 * @param separator
	 * @return
	 */
	public static List<PositionInProtein> parseStringToPositionInProtein(String string, String separator) {
		final List<PositionInProtein> ret = new ArrayList<PositionInProtein>();
		final List<String> subStrings = new ArrayList<String>();
		if (string.contains(separator)) {
			final String[] split = string.split(separator);
			for (final String string2 : split) {
				if (isNumber(string2.charAt(0))) {
					final int lastIndex = subStrings.size() - 1;
					subStrings.set(lastIndex, subStrings.get(lastIndex) + separator + string2);
				} else {
					subStrings.add(string2);
				}
			}
		} else {
			subStrings.add(string);
		}

		for (final String string2 : subStrings) {
			if (string2.contains(SEPARATOR)) {
				try {
					final String[] split = string2.split(SEPARATOR);
					final String positionInteger = String.valueOf(split[1].substring(1));
					final PositionInProtein positionInProtein = new PositionInProtein(Integer.valueOf(positionInteger),
							split[1].charAt(0), split[0]);
					ret.add(positionInProtein);
				} catch (final NumberFormatException e) {

				}
			}
		}
		return ret;
	}

	public char getAa() {
		return aa;
	}
}
