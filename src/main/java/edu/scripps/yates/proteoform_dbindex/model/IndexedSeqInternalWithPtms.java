package edu.scripps.yates.proteoform_dbindex.model;

import java.util.ArrayList;
import java.util.List;

/**
 * internal temporarly representation of indexed sequence before sequences are
 * merged into IndexedSequence
 * 
 * two IndexedSeqInternal are equal if they sequence strings are equal
 * 
 * two IndexedSeqInternal sequences are comparable by their mass
 * 
 * @author Adam
 */
public class IndexedSeqInternalWithPtms implements Comparable<IndexedSeqInternalWithPtms> {

	public IndexedSeqInternalWithPtms(double mass, short offset, short length, int proteinId, String sequence,
			String protDescription, List<PTM> ptms) {
		this.mass = mass;
		this.offset = offset;
		this.length = length;
		this.proteinId = proteinId;
		this.sequence = sequence;
		this.protDescription = protDescription;
		this.ptms = ptms;

	}

	public IndexedSeqInternalWithPtms(double mass, short offset, short length, int proteinId, String sequence,
			List<PTM> ptms) {
		this(mass, offset, length, proteinId, sequence, null, ptms);

	}

	private final double mass;
	private final short offset;
	private final short length;
	private final int proteinId;
	private final String protDescription;
	private final String sequence;
	private List<PTM> ptms;

	@Override
	public int hashCode() {
		int hash = 7;
		hash = 29 * hash + (sequence != null ? sequence.hashCode() : 0);
		return hash;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		final IndexedSeqInternalWithPtms other = (IndexedSeqInternalWithPtms) obj;
		if ((sequence == null) ? (other.sequence != null) : !sequence.equals(other.sequence)) {
			return false;
		}
		return true;
	}

	@Override
	public int compareTo(IndexedSeqInternalWithPtms o) {
		if (getMass() < o.getMass()) {
			return -1;
		} else if (getMass() > o.getMass()) {
			return 1;
		} else {
			return 0;
		}
	}

	public void addPTM(PTM ptm) {
		if (getPtms() == null) {
			ptms = new ArrayList<PTM>();
		}
		getPtms().add(ptm);
	}

	public int getProteinId() {
		return proteinId;
	}

	public double getMass() {
		return mass;
	}

	public short getOffset() {
		return offset;
	}

	public short getLength() {
		return length;
	}

	public List<PTM> getPtms() {
		return ptms;
	}

	/**
	 * Returns the size of the byte array that stores this information
	 * 
	 * @return
	 */
	public int byteSize() {
		final int numPTMs = ptms != null ? ptms.size() : 0;
		final int ret = 8 + // mass
				2 + // position in protein
				2 + // length
				(numPTMs * 3) + // ptms(1+2)
				1 + // no more ptms
				4 + // protein id
				4;// no more protein ids;
		return ret;
	}

	public String getSequence() {
		return sequence;
	}

}
