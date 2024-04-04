package edu.scripps.yates.utilities.masses;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;
import org.springframework.core.io.ClassPathResource;

public class UniprotAANomenclature {
	private final static Logger log = Logger.getLogger(UniprotAANomenclature.class);
	private final static Map<Character, UniprotAminoacid> aaByOneLetterCode = new HashMap<Character, UniprotAminoacid>();

	private final static Map<String, UniprotAminoacid> aaByAminoacidName = new HashMap<String, UniprotAminoacid>();

	private final static Map<String, UniprotAminoacid> aaByThreeLetterCode = new HashMap<String, UniprotAminoacid>();
	private final static List<UniprotAminoacid> aaList = new ArrayList<UniprotAminoacid>();
	private static UniprotAANomenclature instance;

	private UniprotAANomenclature() {
		try {
			final List<String> lines = Files
					.readAllLines(Paths.get(new ClassPathResource("uniprotAANomemclature.txt").getURI()));
			// lines are 3 elements separated by tabs, where the first is the
			// letter code, then the 3 letter aa and then the aa name
			for (final String line : lines) {
				if (line.contains("\t")) {
					final String[] split = line.split("\t");
					if (split[0].trim().length() == 1) {
						final char aa = split[0].trim().charAt(0);
						final String threeLetterCode = split[1].trim();
						final String aaName = split[2].trim();
						final UniprotAminoacid uniprotAA = new UniprotAminoacid(aa, threeLetterCode, aaName);
						aaByAminoacidName.put(aaName, uniprotAA);
						aaByOneLetterCode.put(aa, uniprotAA);
						aaByThreeLetterCode.put(threeLetterCode, uniprotAA);
						aaList.add(uniprotAA);
					}
				}
			}
		} catch (final IOException e) {
			e.printStackTrace();
			log.error(e);
		}
	}

	public static UniprotAANomenclature getInstance() {
		if (instance == null) {
			instance = new UniprotAANomenclature();
		}
		return instance;
	}

	public UniprotAminoacid getByAminoacidName(String aminacidName) {
		return aaByAminoacidName.get(aminacidName);
	}

	public UniprotAminoacid getByThreeLetterCode(String threeLetterCode) {
		return aaByThreeLetterCode.get(threeLetterCode);
	}

	public UniprotAminoacid getByOneLetterCode(char aa) {
		return aaByOneLetterCode.get(aa);
	}
}
