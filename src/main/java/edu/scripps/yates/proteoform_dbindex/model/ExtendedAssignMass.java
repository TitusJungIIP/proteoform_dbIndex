package edu.scripps.yates.proteoform_dbindex.model;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.masses.AssignMass;

public class ExtendedAssignMass extends AssignMass {
	private final static Logger log = Logger.getLogger(ExtendedAssignMass.class);
	private final static Charset CHARSET = StandardCharsets.UTF_8;

	private final Map<String, PTMCodeObj> ptmsPerKey = new HashMap<String, PTMCodeObj>();
	private final Map<Short, PTMCodeObj> ptmsByPTMCode = new HashMap<Short, PTMCodeObj>();
	private final File ptmCodesFile;
	private short shortCounter = 1;
	private static ExtendedAssignMass instance;
	private boolean fileExists = false;

	public static ExtendedAssignMass getInstance(boolean useMono, File ptmCodesFile) {
		if (instance == null) {
			instance = new ExtendedAssignMass(useMono, ptmCodesFile);
		}
		return instance;
	}

	private ExtendedAssignMass(boolean useMono, File ptmCodesFile) {
		super(useMono);
		this.ptmCodesFile = ptmCodesFile;
		if (!ptmCodesFile.exists()) {
			new File(ptmCodesFile.getParent()).mkdirs();
		} else {
			fileExists = true;
		}
		load();
	}

	public short getPTMCode(String description, Double massDiff) throws IOException {
		if (ptmsPerKey.containsKey(description)) {
			return ptmsPerKey.get(description).getPtmCode();
		} else {
			final PTMCodeObj newPTM = getNewPTMCodeObj(description, massDiff);
			return newPTM.getPtmCode();
		}
	}

	private void write(PTMCodeObj newPTM) throws IOException {
		OutputStreamWriter out = null;
		if (!fileExists) {
			out = new OutputStreamWriter(new FileOutputStream(ptmCodesFile), CHARSET);
			fileExists = true;
		} else {
			out = new OutputStreamWriter(new FileOutputStream(ptmCodesFile, true), CHARSET);
		}

		out.write(newPTM.toString() + "\n");
		out.close();
	}

	private PTMCodeObj getNewPTMCodeObj(String description, Double massDiff) throws IOException {
		final short ptmCode = getNewPTMCode();

		final PTMCodeObj ret = new PTMCodeObj(description, ptmCode, massDiff);
		write(ret);
		addPTMToMaps(ret);
		return ret;
	}

	private short getNewPTMCode() {
		return shortCounter++;
	}

	private void load() {

		try {
			if (!ptmCodesFile.exists()) {
				return;
			}
			final InputStreamReader isr = new InputStreamReader(new FileInputStream(ptmCodesFile), CHARSET);
			final BufferedReader buffReader = new BufferedReader(isr);
			String line = null;
			while ((line = buffReader.readLine()) != null) {
				final PTMCodeObj ptm = PTMCodeObj.getFromLine(line);
				addPTMToMaps(ptm);
			}
			buffReader.close();
		} catch (final IOException e) {
			e.printStackTrace();
		}
	}

	private void addPTMToMaps(PTMCodeObj ptm) {
		ptmsPerKey.put(ptm.getDescription(), ptm);
		ptmsByPTMCode.put(ptm.getPtmCode(), ptm);
	}

	public static boolean isAA(char aa) {
		return aaMassMono[aa] != 0.0;
	}

	/**
	 * Calculates the monoisotopic mass from the peptide like:<br>
	 * ASDFASDF[-12.123]SADFSAD[M->K]ASDF
	 * 
	 * @param modifiedPeptide
	 * @return
	 */

	public double calculateMass(String modifiedPeptide) {
		double mass = 0.0;
		StringBuilder change = null;
		for (int i = 0; i < modifiedPeptide.length(); i++) {
			final char aa = modifiedPeptide.charAt(i);
			if (change != null) {
				if (aa == ']') {
					// change as: ASDF->IASDFM
					// parse change
					if (change.toString().contains("->")) {
						final String[] split = change.toString().split("->");

						if (split.length == 2) {
							mass += calculateMass(split[1]);
						} else {
							if (split[0].length() > 100) {
								log.info("asdf");
							}
						}
						change = null;
					} else {
						final Double massPTM = Double.valueOf(change.toString());
						mass += massPTM;
						change = null;
					}
				} else {
					change.append(aa);
				}
			} else if (isAA(aa)) {
				mass += getMass(aa);
			} else {
				if (aa == '[') {
					change = new StringBuilder();
				} else {
					throw new IllegalArgumentException("Error parsing " + modifiedPeptide);
				}
			}
		}
		return mass;
	}

	public PTMCodeObj getPTMByDescription(String description, Double massDiff) throws IOException {
		if (ptmsPerKey.containsKey(description)) {
			return ptmsPerKey.get(description);
		} else {
			final PTMCodeObj newPTM = getNewPTMCodeObj(description, massDiff);
			return newPTM;
		}
	}

	public PTMCodeObj getPTMbyPTMCode(short ptmCode) {

		if (ptmsByPTMCode.containsKey(ptmCode)) {
			return ptmsByPTMCode.get(ptmCode);
		}
		return null;
	}

}
