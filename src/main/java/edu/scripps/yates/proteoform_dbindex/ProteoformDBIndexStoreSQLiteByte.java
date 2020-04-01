package edu.scripps.yates.proteoform_dbindex;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import edu.scripps.yates.dbindex.Constants;
import edu.scripps.yates.dbindex.DBIndexStoreSQLiteByte;
import edu.scripps.yates.dbindex.ProteinCache;
import edu.scripps.yates.dbindex.Util;
import edu.scripps.yates.proteoform_dbindex.model.ExtendedAssignMass;
import edu.scripps.yates.proteoform_dbindex.model.IndexedSeqInternalWithPtms;
import edu.scripps.yates.proteoform_dbindex.model.PTM;
import edu.scripps.yates.proteoform_dbindex.util.ByteArrayUtil;
import edu.scripps.yates.utilities.bytes.DynByteBuffer;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexSearchParams;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedSequence;
import edu.scripps.yates.utilities.fasta.dbindex.ResidueInfo;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.TIntHashSet;

public class ProteoformDBIndexStoreSQLiteByte extends DBIndexStoreSQLiteByte {

	private final ExtendedAssignMass extendedAssignMass;

	protected ProteoformDBIndexStoreSQLiteByte(int bucketId, ProteinCache proteinCache,
			ExtendedAssignMass extendedAssignMass) {
		super(bucketId, proteinCache);
		this.extendedAssignMass = extendedAssignMass;
	}

	protected ProteoformDBIndexStoreSQLiteByte(boolean inMemory, int bucketId, ProteinCache proteinCache,
			ExtendedAssignMass extendedAssignMass) {
		super(inMemory, bucketId, proteinCache);
		this.extendedAssignMass = extendedAssignMass;
	}

	protected ProteoformDBIndexStoreSQLiteByte(DBIndexSearchParams searchParams, boolean inMemory, int bucketId,
			ProteinCache proteinCache, ExtendedAssignMass extendedAssignMass) {
		super(searchParams, inMemory, bucketId, proteinCache);
		this.extendedAssignMass = extendedAssignMass;
	}

	@Override
	public void addSequence(double precMass, int seqOffset, int seqLength, String sequence, String resLeft,
			String resRight, long proteinId) throws DBIndexStoreException {
		if (!inited) {
			throw new DBIndexStoreException("Indexer is not initialized");
		}

		totalSeqCount++;
		try {
			final List<PTM> ptms = PTM.extractPTMsFromSequenceBasedOnResultingSequence(sequence, extendedAssignMass);
			final char seqOffsetChar = DynByteBuffer.toChar(DynByteBuffer.toByteArray(seqOffset));

			updateCachedData(precMass, seqOffsetChar, (short) seqLength, (int) proteinId, ptms);

			if (true) {
				// despite commiting each sequence, commit and clear all after
				// 100mil
				// to free up memory, optimize buffers, and ensure a commit
				if (seqCount++ == FULL_CACHE_COMMIT_INTERVAL) {
					// logger.info( "Commiting all");
					seqCount = 0;
					try {
						commitCachedData();
					} catch (final SQLException ex) {
						logger.error("Error commiting cached data", ex);
					}
				}
			}
		} catch (final IOException e) {
			throw new DBIndexStoreException(e);
		}
	}

	/**
	 * Add the new sequence to in memory cached pre-commit sequence hash
	 *
	 * @param precMass
	 * @param seqOffset
	 * @param seqLength
	 * @param proteinId
	 * @param ptms
	 */

	protected void updateCachedData(double precMass, char seqOffset, short seqLength, int proteinId, List<PTM> ptms) {
		// changed by Salva 11Nov2014, using the value on IndexUtil
		final int rowId = (int) (precMass * sparam.getMassGroupFactor());

		// change by Salva 21Nov2014
		// DynByteBuffer byteBuffer = data[rowId];
		DynByteBuffer byteBuffer = dataMap.get(rowId);
		if (byteBuffer == null) {
			byteBuffer = new DynByteBuffer();
			// change by Salva 21Nov2014
			// data[rowId] = byteBuffer;
			dataMap.put(rowId, byteBuffer);
		}

		byteBuffer.add(ByteArrayUtil.toByteArray(precMass, seqOffset, seqLength, proteinId, ptms));

		// commit this mass after COMMIT_SEQUENCES sequences
		if (true) {
			if (byteBuffer.getSize() > Constants.COMMIT_SEQUENCES * Constants.BYTE_PER_SEQUENCE) {
				try {
					insertSequence(rowId, byteBuffer);
					// clear since we wrote it to db
					byteBuffer.clear();
				} catch (final SQLException e) {
					logger.error("Error commiting mass buffer");
				}
			}
		}
	}

	@Override
	protected void parseAddPeptideInfo(byte[] data, List<IndexedSequence> toInsert, double minMass, double maxMass)
			throws DBIndexStoreException {

		// to collapse multiple sequences into single one, with multiproteins
		final Map<String, List<IndexedSeqInternalWithPtms>> temp = new THashMap<String, List<IndexedSeqInternalWithPtms>>();

		final int dataLength = data.length;

		// if (dataLength % 4 != 0) {
		// throw new RuntimeException("Unexpected number of peptide items: "
		// + dataLength);
		// }
		int i = 0;
		while (i < dataLength) {

			final List<PTM> ptms = new ArrayList<PTM>();
			// System.out.println("pep: " + peptideSequence);

			// we are cheating and supplying protein id instead of peptide
			// id
			// to set it temporarily, before we merge them into a list
			final IndexedSeqInternalWithPtms tempSequence = ByteArrayUtil
					.getIndexedSeqInternalWithPtmsFromByteArray(data, i, (ProteoformProteinCache) proteinCache);
			i += tempSequence.byteSize();
			if (tempSequence.getMass() < minMass || tempSequence.getMass() > maxMass) {
				// skip the sequence, it qualified the bucket, but not the
				// actual mass
				continue;
			}
			final String peptideSequence = proteinCache.getPeptideSequence(tempSequence.getProteinId(),
					tempSequence.getOffset(), tempSequence.getLength());
			List<IndexedSeqInternalWithPtms> sequences = temp.get(peptideSequence);
			if (sequences == null) {
				sequences = new ArrayList<IndexedSeqInternalWithPtms>();
				temp.put(peptideSequence, sequences);
			}
			sequences.add(tempSequence);

		}

		// group the same peptides from many proteins into single peptide
		// with protein id list
		for (final String pepSeqKey : temp.keySet()) {

			final List<IndexedSeqInternalWithPtms> sequences = temp.get(pepSeqKey);

			final TIntHashSet proteinIds = new TIntHashSet();
			for (final IndexedSeqInternalWithPtms tempSeq : sequences) {
				proteinIds.add(tempSeq.getProteinId());
			}

			final IndexedSeqInternalWithPtms firstSeq = sequences.get(0);
			final IndexedSequence mergedSequence = new IndexedSequence(0, firstSeq.getMass(), pepSeqKey, "", "");
			final ArrayList<Integer> proteinIds2 = new ArrayList<Integer>();
			for (final int integer : proteinIds._set) {
				proteinIds2.add(integer);
			}
			mergedSequence.setProteinIds(proteinIds2);
			// set residues
			final String protSequence = proteinCache.getProteinSequence(firstSeq.getProteinId());
			final ResidueInfo residues = Util.getResidues(null, firstSeq.getOffset(), firstSeq.getLength(),
					protSequence);
			mergedSequence.setResidues(residues);
			toInsert.add(mergedSequence);
		}

	}

	@Override
	public void stopAddSeq() throws DBIndexStoreException {
		if (proteinCache instanceof ProteoformProteinCache) {
			try {
				((ProteoformProteinCache) proteinCache).writeBuffer();
			} catch (final IOException e) {
				e.printStackTrace();
				throw new DBIndexStoreException(e);
			}
		}
		super.stopAddSeq();
	}

}
