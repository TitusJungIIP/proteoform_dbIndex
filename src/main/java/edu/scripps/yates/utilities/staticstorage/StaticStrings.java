package edu.scripps.yates.utilities.staticstorage;

import gnu.trove.map.hash.THashMap;

/**
 * Class that storage strings statically, so that we dont create the same string
 * again and again.<br>
 * This is intended to be used in strings objects that repeat a lot, such as for
 * example peptide sequences
 * 
 * @author salvador
 *
 */
public class StaticStrings {
	private final static THashMap<String, String> strings = new THashMap<String, String>();

	/**
	 * Gets a unique instance of a certain string.<br>
	 * This is intended to be used in strings objects that repeat a lot, such as
	 * for example peptide sequences.<br>
	 * BE AWARE THAT THE STRING RETURNED WILL BE SHARED BY MANY OTHER OBJECTS,
	 * SO DO NOT MODIFY IT.
	 * 
	 * @param string
	 * @return
	 */
	public synchronized static String getUniqueInstance(String string) {
		if (string == null) {
			return null;
		}
		if (strings.containsKey(string)) {
			return strings.get(string);
		} else {
			strings.put(string, string);
			return string;
		}
	}

}
