package edu.scripps.yates.proteoform_dbindex.model;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import edu.scripps.yates.utilities.strings.StringUtils;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.THashMap;

public class PhosphositeDB {
	private final static File dbPath = new File("D:\\Downloads\\Phosphosite_PTM_seq.fasta\\Phosphosite_PTM_seq.fasta");
	private final Map<String, TIntArrayList> phosphorilatedPositionsByUniprotACC = new THashMap<String, TIntArrayList>();
	private static final char[] PHOSPHORILATED_SITES = { 'y', 's', 't' };
	private int totalPhosphoSites;
	private final Map<String, String> proteinSeqs = new HashMap<String, String>();

	public PhosphositeDB(String species) throws IOException {

		final FileReader fr = new FileReader(dbPath);
		final BufferedReader br = new BufferedReader(fr);
		String line;
		StringBuilder seq = null;
		String uniprotACC = null;
		while ((line = br.readLine()) != null) {
			if (line.startsWith(">")) {
				if (seq != null) {
					final TIntArrayList positions = new TIntArrayList();
					for (final char phosphoAA : PHOSPHORILATED_SITES) {
						positions.addAll(StringUtils.allPositionsOf(seq.toString(), phosphoAA));
					}
					positions.sort();
					totalPhosphoSites += positions.size();
					phosphorilatedPositionsByUniprotACC.put(uniprotACC, positions);
					proteinSeqs.put(uniprotACC, seq.toString());
				}
				uniprotACC = null;
				seq = null;
				if (line.contains(species)) {
					uniprotACC = line.substring(line.lastIndexOf("|") + 1);
				}
			} else if (uniprotACC != null) {
				if (seq == null) {
					seq = new StringBuilder();
				}
				seq.append(line.trim());
			}

		}
		br.close();
	}

	public TIntArrayList getPhosphorilatedPositions(String uniprotACC) {
		return phosphorilatedPositionsByUniprotACC.get(uniprotACC);
	}

	public int getTotalPhosphoSites() {
		return totalPhosphoSites;
	}

	public String getProteinSeq(String accession) {
		return proteinSeqs.get(accession);
	}

}
