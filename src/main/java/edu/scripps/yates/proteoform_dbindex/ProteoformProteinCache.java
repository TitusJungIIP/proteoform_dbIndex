package edu.scripps.yates.proteoform_dbindex;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import edu.scripps.yates.dbindex.DBIndexStoreException;
import edu.scripps.yates.dbindex.ProteinCache;
import edu.scripps.yates.proteoform_dbindex.model.ExtendedAssignMass;
import edu.scripps.yates.proteoform_dbindex.model.PTM;
import edu.scripps.yates.proteoform_dbindex.model.PTMCodeObj;

public class ProteoformProteinCache extends ProteinCache {
	private final static Logger log = Logger.getLogger(ProteoformProteinCache.class);
	public static final String FILE_NAME = "proteinCache.txt";
	private final ExtendedAssignMass extendedAssignMass;
	private final File proteinCacheFile;
	private boolean loaded;
	private final static Charset CHARSET = StandardCharsets.UTF_8;

	public ProteoformProteinCache(ExtendedAssignMass extendedAssignMass, File proteinCacheFile) {
		this.extendedAssignMass = extendedAssignMass;
		this.proteinCacheFile = proteinCacheFile;
	}

	public String getPeptideSequence(int proteinId, short seqOffset, short seqLen, List<PTM> ptms)
			throws DBIndexStoreException {
		String protSeq = null;
		try {
			protSeq = getProteinSequence(proteinId);
		} catch (final Exception e) {
			e.printStackTrace();
			final String message = "Error trying to get protein from protein cache with index " + proteinId;
			log.error(message);
			log.error("Using beginIndex:" + seqOffset + " and length: " + seqLen);
			log.error(e.getMessage());
			throw new DBIndexStoreException(message, e);
		}

		try {
			final String peptideSeq = protSeq.substring(seqOffset, seqOffset + seqLen);
			return applyPTMs(peptideSeq, ptms);
		} catch (final Exception e) {
			e.printStackTrace();
			log.error("Error tryin to get substring from " + protSeq);
			log.error("Using beginIndex:" + seqOffset + " and length: " + seqLen);
			log.error(e.getMessage());
		}

		return null;
	}

	private String applyPTMs(String peptideSeq, List<PTM> ptms) {
		if (ptms == null || ptms.isEmpty()) {
			return peptideSeq;
		}
		final StringBuilder sb = new StringBuilder();
		final Map<Integer, PTM> ptmsByPosition = new HashMap<Integer, PTM>();
		for (final PTM ptm : ptms) {
			ptmsByPosition.put(ptm.getPosInPeptide(), ptm);
		}

		for (int pos = 1; pos <= peptideSeq.length(); pos++) {
			if (ptmsByPosition.containsKey(pos)) {
				final PTM ptm = ptmsByPosition.get(pos);
				final PTMCodeObj ptmCodeObj = extendedAssignMass.getPTMbyPTMCode(ptm.getPtmCode());
				if (!ptmCodeObj.getDescription().contains("->")) {
					sb.append(peptideSeq.charAt(pos - 1));
				}
				sb.append("[" + ptmCodeObj.getDescription() + "]");
				if (ptmCodeObj.getDescription().contains("->")) {
					final String original = ptmCodeObj.getDescription().split("->")[0];
					pos += original.length() - 1;
				}
			} else {
				sb.append(peptideSeq.charAt(pos - 1));
			}
		}
		if (sb.toString().equals("MEEPQS[P->S]PSVEPPLSQETFSDLWK")) {
			log.info(applyPTMs(peptideSeq, ptms));
		}
		return sb.toString();
	}

	private void load() throws IOException {
		if (loaded) {
			return;
		}
		sequences.clear();
		defs.clear();
		if (proteinCacheFile.exists()) {
			final InputStreamReader isr = new InputStreamReader(new FileInputStream(proteinCacheFile), CHARSET);
			final BufferedReader br = new BufferedReader(isr);
			String line = null;
			while ((line = br.readLine()) != null) {
				final String[] split = line.split("\t");
				defs.add(split[1]);
				if (split.length > 2) {
					sequences.add(split[2]);
				}
			}
			br.close();
		}
	}

	@Override
	public void addProtein(String def) {

		try {
			load();
			if (!proteinCacheFile.getParentFile().exists()) {
				proteinCacheFile.getParentFile().mkdirs();
			}
			final int proteinIndex = defs.size();
			final OutputStreamWriter out = new OutputStreamWriter(new FileOutputStream(proteinCacheFile, true),
					CHARSET);
			out.write(proteinIndex + "\t" + def + "\n");
			out.close();

			super.addProtein(def);

		} catch (final IOException e) {
			e.printStackTrace();
		}
	}

	@Override
	public int addProtein(String def, String protein) {
		try {
			load();
			if (!proteinCacheFile.getParentFile().exists()) {
				proteinCacheFile.getParentFile().mkdirs();
			}
			final int proteinIndex = defs.size();
			final OutputStreamWriter out = new OutputStreamWriter(new FileOutputStream(proteinCacheFile, true),
					CHARSET);
			out.write(proteinIndex + "\t" + def + "\t" + protein + "\n");
			out.close();

			return super.addProtein(def, protein);
		} catch (final FileNotFoundException e) {
		} catch (final IOException e) {
			e.printStackTrace();
		}
		return -1;
	}

}
