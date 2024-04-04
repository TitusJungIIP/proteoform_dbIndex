package edu.scripps.yates.utilities.hash;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

public class HashGenerator {
	public static int generateHash(Collection<String> strings) {
		final List<String> list = new ArrayList<String>();
		list.addAll(strings);
		return generateHash(list);
	}

	public static int generateHash(String... strings) {
		final List<String> list = new ArrayList<String>();
		for (final String string : strings) {
			list.add(string);
		}
		return generateHash(list);
	}

	public static int generateHash(List<String> list) {

		Collections.sort(list);
		int hash = 7;
		for (final String string : list) {
			for (int i = 0; i < string.length(); i++) {
				hash = hash * 31 + string.charAt(i);
			}
		}
		return hash;
	}
}
