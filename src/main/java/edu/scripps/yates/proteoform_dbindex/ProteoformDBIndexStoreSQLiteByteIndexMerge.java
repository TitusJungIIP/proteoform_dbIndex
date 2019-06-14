package edu.scripps.yates.proteoform_dbindex;

import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import edu.scripps.yates.dbindex.Constants;
import edu.scripps.yates.dbindex.DBIndexStoreSQLiteByteIndexMerge;
import edu.scripps.yates.proteoform_dbindex.model.ExtendedAssignMass;
import edu.scripps.yates.proteoform_dbindex.model.IndexedSeqInternalWithPtms;
import edu.scripps.yates.proteoform_dbindex.model.IndexedSeqMergedWithPtms;
import edu.scripps.yates.proteoform_dbindex.model.IndexedSequenceWithPTMs;
import edu.scripps.yates.proteoform_dbindex.model.PTM;
import edu.scripps.yates.proteoform_dbindex.util.ByteArrayUtil;
import edu.scripps.yates.utilities.bytes.DynByteBuffer;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexSearchParams;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedSequence;
import gnu.trove.map.hash.THashMap;

public class ProteoformDBIndexStoreSQLiteByteIndexMerge extends DBIndexStoreSQLiteByteIndexMerge {
	private static final org.apache.log4j.Logger logger = org.apache.log4j.Logger
			.getLogger(ProteoformDBIndexStoreSQLiteByteIndexMerge.class);
	private final ExtendedAssignMass extendedAssignMass;

	public ProteoformDBIndexStoreSQLiteByteIndexMerge(boolean inMemory, int bucketId,
			ProteoformProteinCache proteinCache, ExtendedAssignMass extendedAssignMass) {
		super(inMemory, bucketId, proteinCache);
		this.extendedAssignMass = extendedAssignMass;

	}

	public ProteoformDBIndexStoreSQLiteByteIndexMerge(DBIndexSearchParams sParams, boolean inMemory, int bucketId,
			ProteoformProteinCache proteinCache, ExtendedAssignMass extendedAssignMass) {
		super(sParams, inMemory, bucketId, proteinCache);
		this.extendedAssignMass = extendedAssignMass;

	}

	private ProteoformProteinCache getProteoformProteinCache() {
		return (ProteoformProteinCache) proteinCache;
	}

	/**
	 * Parse the encoded sequences and return them in toInsert
	 *
	 * @param data
	 * @param toInsert
	 * @param minMass  used to eliminate not matching masses in the 1ppm row
	 * @param maxMass  used to eliminate not matching masses in the 1ppm row
	 * @throws DBIndexStoreException
	 */
	@Override
	protected void parseAddPeptideInfo(byte[] data, List<IndexedSequence> toInsert, double minMass, double maxMass)
			throws DBIndexStoreException {

		final int dataLength = data.length;

		// if (dataLength % 4 != 0) {
		// throw new RuntimeException("Unexpected number of peptide items: "
		// + dataLength);
		// }
		int i = 0;
		while (i < dataLength) {
			final IndexedSequenceWithPTMs pep = ByteArrayUtil.getIndexedPeptideFromByteArray(data, i,
					getProteoformProteinCache());
			i += pep.getBySize();
			toInsert.add(pep);

		} // end for every byte

	}

	/**
	 * Parse the encoded sequences and return them in toInsert This is version with
	 * multiple mass ranges (slightly slower as we need to check every range if
	 * sequence qualifies)
	 *
	 * @param data      the data row for all peptides grouped by a single key, i.e.
	 *                  1 ppm value
	 * @param toInsert  parsed sequences that qualify
	 * @param minMasses used to eliminate not matching masses in the 1ppm row
	 * @param maxMasses used to eliminate not matching masses in the 1ppm row
	 * @throws DBIndexStoreException
	 */
	@Override
	protected void parseAddPeptideInfo(byte[] data, List<IndexedSequence> toInsert, double[] minMasses,
			double[] maxMasses) throws DBIndexStoreException {

		final int dataLength = data.length;

		// if (dataLength % 4 != 0) {
		// throw new RuntimeException("Unexpected number of peptide items: "
		// + dataLength);
		// }

		final int numRanges = minMasses.length;
		int i = 0;
		// go over every sequence in the data row
		while (i < dataLength) {
			final IndexedSequenceWithPTMs pep = ByteArrayUtil.getIndexedPeptideFromByteArray(data, i,
					getProteoformProteinCache());
			i += pep.getBySize();
			final double seqMass = pep.getMass();
			// check how the actual sequence mass fits in all mass ranges
			// requested
			// we should bail out if pass all the ranges
			boolean greaterThanMax = true;
			// boolean lesserThanMin = true;
			boolean qualifiesRange = false; // check if qualifies in any
											// range supplied
			for (int range = 0; range < numRanges; ++range) {
				if (seqMass < maxMasses[range]) {
					greaterThanMax = false;
				}

				// if (seqMass > minMasses[range]) {
				// lesserThanMin = false;
				// }

				if (seqMass >= minMasses[range] && seqMass <= maxMasses[range]) {
					qualifiesRange = true;
				}

				// if fits, bail out of this check earlier, otherwise check
				// all ranges
				if (qualifiesRange == true) {
					break;
				}

				// check if makes any range

			}

			// since it's sorted
			// if current mass great than maxmass, break early - optimize
			if (greaterThanMax && !qualifiesRange) {
				break;
			}

			toInsert.add(pep);

		} // end for every byte

	}

	@Override
	protected void mergePeptides() throws DBIndexStoreException {
		// logger.info( "Merging sequence index start");
		ResultSet massKeysRs = null;
		try {
			con.setAutoCommit(false);

			massKeysRs = getMassDataStatement.executeQuery();

			// go over each row and merge peptides
			while (massKeysRs.next()) {
				final int massKey = massKeysRs.getInt(1);

				final byte[] seqData = massKeysRs.getBytes(2);
				if (massKey == 27662450) {
					logger.info(seqData);
				}
				final byte[] seqDataMerged = getMergedData(seqData);
				if (seqDataMerged != null) {
					updateSeqStatement.setBytes(1, seqDataMerged);
					updateSeqStatement.setInt(2, massKey);
					updateSeqStatement.execute();
				}
			}

			con.commit();

		} catch (final SQLException ex) {
			ex.printStackTrace();
			logger.error("Error merging sequences", ex);
		} finally {
			try {
				if (massKeysRs != null) {
					massKeysRs.close();
				}
				con.setAutoCommit(true);
			} catch (final SQLException ex) {
				ex.printStackTrace();
				logger.error("Error restoring autocommit", ex);
			}
		}
		// logger.info( "Merging sequence index end");
	}

	@Override
	protected byte[] getMergedData(byte[] seqData) throws DBIndexStoreException {

		// to collapse multiple sequences into single one, with multiproteins
		final Map<String, List<IndexedSeqInternalWithPtms>> temp = new THashMap<String, List<IndexedSeqInternalWithPtms>>();

		final int dataLength = seqData.length;

		// modification made by Salvador Martinez on 5/21/2014
		int i = 0;
		while (i < dataLength) {
			final IndexedSeqInternalWithPtms tempSequence = ByteArrayUtil
					.getIndexedSeqInternalWithPtmsFromByteArray(seqData, i, getProteoformProteinCache());
			i += tempSequence.byteSize();
			final String peptideSequence = tempSequence.getSequence();

			List<IndexedSeqInternalWithPtms> sequences = temp.get(peptideSequence);
			if (sequences == null) {
				sequences = new ArrayList<IndexedSeqInternalWithPtms>();
				temp.put(peptideSequence, sequences);
			}
			sequences.add(tempSequence);

		}

		// sort sequences by masses for search optimization
		final List<IndexedSeqMergedWithPtms> sortedMerged = new ArrayList<IndexedSeqMergedWithPtms>();

		// group the same peptides from many proteins into single peptide
		// with protein id list
		// sort them by mass for optimization
		for (final String pepSeqKey : temp.keySet()) {
			// for each sequence str

			final List<IndexedSeqInternalWithPtms> sequences = temp.get(pepSeqKey);

			final List<Integer> proteinIds = new ArrayList<Integer>();
			for (final IndexedSeqInternalWithPtms tempSeq : sequences) {
				proteinIds.add(tempSeq.getProteinId());
			}

			// make sure the 1st protein id is that of the first sequence
			final IndexedSeqInternalWithPtms firstSeq = sequences.get(0);

			final IndexedSeqMergedWithPtms merged = new IndexedSeqMergedWithPtms(firstSeq.getMass(),
					firstSeq.getOffset(), firstSeq.getLength(), proteinIds, firstSeq.getPtms());

			sortedMerged.add(merged);

		} // end of this sequence

		Collections.sort(sortedMerged);

		final DynByteBuffer seqDataMerged = new DynByteBuffer();
		// write out grouped merged peptides to byte buf
		for (final IndexedSeqMergedWithPtms merged : sortedMerged) {

			final byte[] array = merged.toByteArray();

			// add separator
			seqDataMerged.add(array);

		} // end of this sequence

		return seqDataMerged.getData();
	}

	@Override
	public void addSequence(double precMass, int seqOffset, int seqLength, String sequence, String resLeft,
			String resRight, long proteinId) throws DBIndexStoreException {
		if (!inited) {
			throw new DBIndexStoreException("Indexer is not initialized");
		}
		try {
			totalSeqCount++;
			final List<PTM> ptms = PTM.extractPTMsFromSequence(sequence, extendedAssignMass);
			if (seqOffset > 32767 || seqOffset < 0) {
				logger.info("track this downt to bytes");
			}
			new DynByteBuffer();
			final char seqOffsetChar = DynByteBuffer.toChar(DynByteBuffer.toByteArray(seqOffset));
			final short int1 = DynByteBuffer.toShort(DynByteBuffer.toByteArray(seqOffsetChar));
			if (int1 < 0) {
				logger.info("Some error ocurred! 2 bytes of offset are overflown");
			}
			updateCachedData(precMass, seqOffsetChar, Integer.valueOf(seqLength).shortValue(),
					Long.valueOf(proteinId).intValue(), ptms);

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
