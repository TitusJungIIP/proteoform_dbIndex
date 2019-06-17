package edu.scripps.yates.proteoform_dbindex.model;

import java.util.List;

import edu.scripps.yates.dbindex.DBIndexStoreSQLiteByteIndexMerge;
import edu.scripps.yates.utilities.bytes.DynByteBuffer;

/**
 * internal temporarly representation of indexed sequence before sequences are
 * merged into IndexedSequence
 * 
 * two IndexedSeqInternal are equal if they sequence strings are equal
 * 
 * two IndexedSeqInternal sequences are comparable by their mass
 * 
 */
public class IndexedSeqMergedWithPtms implements Comparable<IndexedSeqMergedWithPtms> {

	public IndexedSeqMergedWithPtms(double mass, int offset, short length, List<Integer> proteinIds, List<PTM> ptms) {
		setMass(mass);
		setOffset(offset);
		setLength(length);
		setProteinIds(proteinIds);
		setPtms(ptms);
	}

	private double mass;
	private int offset;
	private short length;
	private List<Integer> proteinIds;
	private List<PTM> ptms;

	@Override
	public int compareTo(IndexedSeqMergedWithPtms o) {
		if (getMass() < o.getMass()) {
			return -1;
		} else if (getMass() > o.getMass()) {
			return 1;
		} else {
			return 0;
		}
	}

	private String getPtmSring() {
		final StringBuilder sb = new StringBuilder();
		if (getPtms() != null) {
			for (final PTM ptm : getPtms()) {
				sb.append(ptm.getPosInPeptide() + "[" + ptm.getPtmCode() + "]");
			}
		}
		return sb.toString();
	}

	@Override
	public String toString() {
		return "IndexedSeqMerged{" + "mass=" + getMass() + ", offset=" + getOffset() + ", length=" + getLength()
				+ ", proteinIds=" + getProteinIds() + ", ptms=" + getPtmSring() + "}";
	}

	public double getMass() {
		return mass;
	}

	public void setMass(double mass) {
		this.mass = mass;
	}

	public int getOffset() {
		return offset;
	}

	public void setOffset(int offset) {
		this.offset = offset;
	}

	public short getLength() {
		return length;
	}

	public void setLength(short length) {
		this.length = length;
	}

	public List<PTM> getPtms() {
		return ptms;
	}

	public void setPtms(List<PTM> ptms) {
		this.ptms = ptms;
	}

	public List<Integer> getProteinIds() {
		return proteinIds;
	}

	public void setProteinIds(List<Integer> proteinIds) {
		this.proteinIds = proteinIds;
	}

	public byte[] toByteArray() {
		final DynByteBuffer ret = new DynByteBuffer();
		final byte[] seqMassB = DynByteBuffer.toByteArray(getMass());
		ret.add(seqMassB);

		final byte[] seqOffsetB = DynByteBuffer.toByteArray(getOffset());
		ret.add(seqOffsetB);

		final byte[] seqLengthB = DynByteBuffer.toByteArray(getLength());
		ret.add(seqLengthB);
		final List<PTM> ptms = getPtms();
		if (ptms != null && !ptms.isEmpty()) {
			for (final PTM ptm : ptms) {
				ret.add(ptm.getPosInPeptideByte());
				ret.add(DynByteBuffer.toByteArray(ptm.getPtmCode()));
			}
		}

		final byte b = -1;
		ret.add(b);
		for (final int proteinId : getProteinIds()) {
			final byte[] proteinIdB = DynByteBuffer.toByteArray(proteinId);
			ret.add(proteinIdB);
		}
		// add separator
		ret.add(DBIndexStoreSQLiteByteIndexMerge.SEQ_SEPARATOR);
		return ret.getData();
	}
}
