package edu.scripps.yates.proteoform_dbindex.model;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.files.ZipManager;
import edu.scripps.yates.utilities.strings.StringUtils;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.THashMap;

public class PhosphositeDB {
	private final static Logger log = Logger.getLogger(PhosphositeDB.class);
	private final static File dbPath = new File("D:\\Downloads\\Phosphosite_PTM_seq.fasta\\Phosphosite_PTM_seq.fasta");
	private final Map<String, TIntArrayList> phosphorilatedPositionsByUniprotACC = new THashMap<String, TIntArrayList>();
	private static final char[] PHOSPHORILATED_SITES = { 'y', 's', 't' };
	private int totalPhosphoSites;
	private final Map<String, String> proteinSeqs = new HashMap<String, String>();
	private boolean ready = false;
	private final File dbFile;
	private final String species;

	/**
	 * Using the default location at
	 * D:\\Downloads\\Phosphosite_PTM_seq.fasta\\Phosphosite_PTM_seq.fasta. Only for
	 * testing purposes.
	 * 
	 * @param species mouse, human, rat, cow...
	 * @throws IOException
	 */
	@Deprecated
	public PhosphositeDB(String species) {
		this(species, dbPath);
	}

	/**
	 * 
	 * @param species mouse, human, rat, cow...
	 * @param dbFile  it can be the Phosphosite_PTM_seq.fasta file or the
	 *                Phosphosite_PTM_seq.fasta.gz file (it is decompressed if
	 *                necessary)
	 * @throws IOException
	 */
	public PhosphositeDB(String species, File dbFile) {

		this.dbFile = dbFile;
		this.species = species;
	}

	private void init() {
		try {
			final File dbFile = ZipManager.decompressFileIfNeccessary(this.dbFile);
			final FileReader fr = new FileReader(dbFile);
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
			ready = true;
		} catch (final IOException e) {
			e.printStackTrace();
			log.error(e.getMessage());
		}
	}

	public TIntArrayList getPhosphorilatedPositions(String uniprotACC) {
		if (!ready) {
			init();
		}
		return phosphorilatedPositionsByUniprotACC.get(uniprotACC);
	}

	public int getTotalPhosphoSites() {
		if (!ready) {
			init();
		}
		return totalPhosphoSites;
	}

	public String getProteinSeq(String accession) {
		if (!ready) {
			init();
		}
		return proteinSeqs.get(accession);
	}

}
