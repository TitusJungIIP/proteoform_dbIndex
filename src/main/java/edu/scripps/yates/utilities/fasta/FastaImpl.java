package edu.scripps.yates.utilities.fasta;

import java.io.UnsupportedEncodingException;

public class FastaImpl implements Fasta {

	// description line of this Fasta sequence
	protected static String defline;
	// the sequence string of this Fasta
	protected byte[] sequence;
	private String seq;
	private String accession = null;
	private String sequestLikeAccession = null;
	private double mPlusH = 0;
	private final boolean proteoform;

	public FastaImpl(String defline, String sequence, boolean isProteoform) {
		FastaImpl.defline = defline;
		proteoform = isProteoform;
		// System.out.println("defline: " + this.defline);
		sequence = sequence.toUpperCase();
		// seq = sequence;
		try {
			seq = sequence;
			this.sequence = sequence.getBytes("US-ASCII");
		} catch (final UnsupportedEncodingException e) {
			System.err.println("Unknow charset");
			System.exit(1);
		}
		// System.out.println("hahaha, defline: " + this.defline);
	}

	@Override
	public String toString() {
		return "Fasta{defline=" + defline + '}';
	}

	@Override
	public void setMPlusH(double mh) {
		mPlusH = mh;
	}

	@Override
	public double getMPlusH() {
		return mPlusH;
	}

	// public double setMPlusH() {
	// return mPlusH;
	// }
	// to order the sequence from long to short
	@Override
	public int compareTo(Fasta f) {
		if (f == null)
			return -1;
		return f.getLength() - getLength();
	}

	@Override
	public String getAccessionWithNoVersion() {
		final int index = accession.indexOf(".");
		if (index == -1) {
			return accession;
		} else {
			return getAccession().substring(0, index);
		}
	}

	@Override
	public String getSequence() {
		if (seq == null) {
			seq = new String(sequence);
		}
		return seq;
		// return new String(sequence);
	}

	@Override
	public byte[] getSequenceAsBytes() {
		return sequence;
	}

	@Override
	public String getOriginalDefline() {

		// return defline.substring(1,defline.length());
		return defline;
	}

	@Override
	public String getDefline() {
		if (defline.startsWith(">")) {
			return defline.substring(1, defline.length());
		}
		return defline;
	}

	@Override
	public byte byteAt(int index) {
		return sequence[index];
	}

	@Override
	public int getLength() {
		return sequence.length;
	}

	// get accession without version
	@Override
	public String getAccession() {
		if (accession == null) {
			// System.out.println("defline: " + defline);
			accession = FastaParser.getACC(getDefline()).getAccession();

		}

		return accession;
	}

	@Override
	public String getSequestLikeAccession() {
		if (sequestLikeAccession == null) {
			sequestLikeAccession = Fasta.getSequestLikeAccession(defline.substring(1));
		}
		return sequestLikeAccession;
	}

	@Override
	public boolean isReversed() {
		return getAccession().startsWith("Re");
	}

	@Override
	public boolean isProteoform() {
		return proteoform;
	}

}
