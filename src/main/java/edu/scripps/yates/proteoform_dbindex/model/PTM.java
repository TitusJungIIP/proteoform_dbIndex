package edu.scripps.yates.proteoform_dbindex.model;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import edu.scripps.yates.utilities.bytes.DynByteBuffer;

public class PTM {
	private final byte posInPeptide;
	private int posInPeptideInt = -1;
	private final short ptmCode;

	public PTM(byte posInPeptide, short ptmCode) {
		this.posInPeptide = posInPeptide;
		this.ptmCode = ptmCode;
	}

	public int getPosInPeptide() {
		if (posInPeptideInt == -1) {
			posInPeptideInt = posInPeptide & (0xff);
		}
		return posInPeptideInt;
	}

	public short getPtmCode() {
		return ptmCode;
	}

	public byte getPosInPeptideByte() {
		return posInPeptide;
	}

	public static List<PTM> extractPTMsFromSequence(String modifiedSequence, ExtendedAssignMass extendedAssignMass)
			throws IOException {
		final List<PTM> ptms = new ArrayList<PTM>();
		int position = 0;
		boolean newChangeStarts = false;
		StringBuilder change = null;
		for (final char aa : modifiedSequence.toCharArray()) {
			if (change == null && ExtendedAssignMass.isAA(aa)) {
				// is an AA
				position++;
			} else {
				if (change != null) {
					if (aa == ']') {
						// parse change
						if (change.toString().contains("->")) {
							int tmp = 0;
							final String[] split = change.toString().split("->");
							if (split.length == 2) {
								tmp = split[1].length();
							}
							final int pos = position - tmp;
							final PTMCodeObj ptmObj = extendedAssignMass.getPTMByDescription(change.toString(), null);
							final PTM ptm = new PTM(DynByteBuffer.toByteArray(pos)[0], ptmObj.getPtmCode());
							ptms.add(ptm);
						} else {
							final Double massdiff = Double.valueOf(change.toString());
							final PTMCodeObj ptmObj = extendedAssignMass.getPTMByDescription(change.toString(),
									massdiff);
							final PTM ptm = new PTM(DynByteBuffer.toByteArray(position)[0], ptmObj.getPtmCode());
							ptms.add(ptm);
						}
						change = null;
						newChangeStarts = false;
					} else {
						change.append(aa);
						if (newChangeStarts) {
							position++;
						}
					}
					if (aa == '>') {
						newChangeStarts = true;
					}
				} else {
					if (aa == '[') {
						change = new StringBuilder();
					} else {
						throw new IllegalArgumentException("ERROR PARSING " + modifiedSequence);
					}
				}
			}

		}

		return ptms;
	}

	public static List<PTM> extractPTMsFromSequenceOLD(String modifiedSequence, ExtendedAssignMass extendedAssignMass)
			throws IOException {
		final List<PTM> ptms = new ArrayList<PTM>();
		int position = 0;
		StringBuilder change = null;
		for (final char aa : modifiedSequence.toCharArray()) {
			if (change == null && ExtendedAssignMass.isAA(aa)) {
				// is an AA
				position++;
			} else {
				if (change != null) {
					if (aa == ']') {
						// parse change
						if (change.toString().contains("->")) {
							position++;
							final PTMCodeObj ptmObj = extendedAssignMass.getPTMByDescription(change.toString(), null);
							final PTM ptm = new PTM(DynByteBuffer.toByteArray(position)[0], ptmObj.getPtmCode());
							ptms.add(ptm);
						} else {
							final Double massdiff = Double.valueOf(change.toString());
							final PTMCodeObj ptmObj = extendedAssignMass.getPTMByDescription(change.toString(),
									massdiff);
							final PTM ptm = new PTM(DynByteBuffer.toByteArray(position)[0], ptmObj.getPtmCode());
							ptms.add(ptm);
						}
						change = null;
					} else {
						change.append(aa);
					}
				} else {
					if (aa == '[') {
						change = new StringBuilder();
					} else {
						throw new IllegalArgumentException("ERROR PARSING " + modifiedSequence);
					}
				}
			}

		}

		return ptms;
	}

}
