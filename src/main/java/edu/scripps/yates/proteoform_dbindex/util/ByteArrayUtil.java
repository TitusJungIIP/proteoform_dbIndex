package edu.scripps.yates.proteoform_dbindex.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import edu.scripps.yates.dbindex.DBIndexStoreSQLiteByteIndexMerge;
import edu.scripps.yates.dbindex.Util;
import edu.scripps.yates.proteoform_dbindex.ProteoformProteinCache;
import edu.scripps.yates.proteoform_dbindex.model.IndexedSeqInternalWithPtms;
import edu.scripps.yates.proteoform_dbindex.model.IndexedSequenceWithPTMs;
import edu.scripps.yates.proteoform_dbindex.model.PTM;
import edu.scripps.yates.utilities.bytes.DynByteBuffer;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;
import edu.scripps.yates.utilities.fasta.dbindex.ResidueInfo;

public class ByteArrayUtil {
	public static IndexedSequenceWithPTMs getIndexedPeptideFromByteArray(byte[] data, int startIndex,
			ProteoformProteinCache proteinCache) throws DBIndexStoreException {
		try {
			int start = startIndex;
			int end = startIndex + 8;
			// mass (double, 8 bytes)
			byte[] slice = Arrays.copyOfRange(data, start, end);
			final double mass = DynByteBuffer.toDouble(slice);

			// offset (int, 4 bytes)
			start = end;
			end = start + 4;
			slice = Arrays.copyOfRange(data, start, end);
			final int offset = DynByteBuffer.toInt(slice);
			if (offset < 0) {
				throw new IllegalArgumentException("Offset overflown");
			}
			// length (short, 2 bytes)
			start = end;
			end = start + 2;
			slice = Arrays.copyOfRange(data, start, end);
			final short length = DynByteBuffer.toShort(slice);

			// ptms (1 byte for position in peptide and 1 byte for ptm code)
			final List<PTM> ptms = new ArrayList<PTM>();
			while (true) {
				start = end;
				end = start + 1;
				slice = Arrays.copyOfRange(data, start, end);
				final byte posInPeptide = slice[0];
				if (posInPeptide == -1) {
					break;
				}
				start = end;
				end = start + 2;
				slice = Arrays.copyOfRange(data, start, end);
				final short ptmCode = DynByteBuffer.toShort(slice);
				final PTM ptm = new PTM(posInPeptide, ptmCode);
				ptms.add(ptm);
			}

			// proteinId (integer, 4 bytes)
			final List<Integer> proteinIds = new ArrayList<Integer>();
			while (true) {
				start = end;
				end = start + 4;
				slice = Arrays.copyOfRange(data, start, end);
				final int proteinID = DynByteBuffer.toInt(slice);
				if (proteinID == DBIndexStoreSQLiteByteIndexMerge.SEQ_SEPARATOR_INT) {
					break;
				}
				proteinIds.add(proteinID);
			}

			final String sequence = proteinCache.getPeptideSequence(proteinIds.get(0), offset, length, ptms);

			final String cleanSequence = FastaParser.cleanSequenceAndApplySequenceVariances(sequence);
			// pre and post will be set below at the ResidueInfo
			final IndexedSequenceWithPTMs ret = new IndexedSequenceWithPTMs(0, mass, cleanSequence, "", "");

			ret.setIsModified(!ptms.isEmpty());
			if (!ptms.isEmpty()) {
				ret.setModSequence(sequence);
			}
			// ret.setProteinDescArray(proteinDescArray);
			ret.setProteinIds(proteinIds);
			for (final PTM ptm : ptms) {
				ret.addPTM(ptm);
			}

			// set residues
			final String protSequence = proteinCache.getProteinSequence(proteinIds.get(0));
			final ResidueInfo residues = Util.getResidues(null, offset, length, protSequence);
			ret.setResidues(residues);
			return ret;

		} catch (final ArrayIndexOutOfBoundsException e) {
			e.printStackTrace();
			throw new DBIndexStoreException(e);
		}
	}

	public static IndexedSeqInternalWithPtms getIndexedSeqInternalWithPtmsFromByteArray(byte[] data, int startIndex,
			ProteoformProteinCache proteinCache) throws DBIndexStoreException {

		int start = startIndex;
		int end = startIndex + 8;
		// mass (double, 8 bytes)
		byte[] slice = Arrays.copyOfRange(data, start, end);
		final double mass = DynByteBuffer.toDouble(slice);

		// offset (int, 4 bytes)
		start = end;
		end = start + 4;
		slice = Arrays.copyOfRange(data, start, end);
		final int offset = DynByteBuffer.toInt(slice);

		// length (short, 2 bytes)
		start = end;
		end = start + 2;
		slice = Arrays.copyOfRange(data, start, end);
		final short length = DynByteBuffer.toShort(slice);

		// ptms (1 byte for position in peptide and 1 byte for ptm code)
		final List<PTM> ptms = new ArrayList<PTM>();
		while (true) {
			start = end;
			end = start + 1;
			slice = Arrays.copyOfRange(data, start, end);
			final byte posInPeptide = slice[0];
			if (posInPeptide == -1) {
				break;
			}
			start = end;
			end = start + 2;
			slice = Arrays.copyOfRange(data, start, end);
			final short ptmCode = DynByteBuffer.toShort(slice);
			final PTM ptm = new PTM(posInPeptide, ptmCode);
			ptms.add(ptm);
		}

		// proteinId (integer, 4 bytes)
		start = end;
		end = start + 4;
		slice = Arrays.copyOfRange(data, start, end);
		final int proteinID = DynByteBuffer.toInt(slice);

		final String sequence = proteinCache.getPeptideSequence(proteinID, offset, length, ptms);
		final IndexedSeqInternalWithPtms ret = new IndexedSeqInternalWithPtms(mass, offset, length, proteinID, sequence,
				ptms);
		return ret;

	}

	public static byte[] toByteArray(double precMass, int seqOffset, short seqLength, int proteinId, List<PTM> ptms) {
		final DynByteBuffer byteBuffer = new DynByteBuffer();

		final byte[] seqMassB = DynByteBuffer.toByteArray(precMass);
		byteBuffer.add(seqMassB);
		if (seqOffset < 0) {
			System.out.println("asdf ");
		}
		final byte[] seqOffsetB = DynByteBuffer.toByteArray(seqOffset);
		byteBuffer.add(seqOffsetB);

		final byte[] seqLengthB = DynByteBuffer.toByteArray(seqLength);
		byteBuffer.add(seqLengthB);

		if (ptms != null && !ptms.isEmpty()) {
			for (final PTM ptm : ptms) {
				final byte pos = ptm.getPosInPeptideByte();
				byteBuffer.add(pos);
				final short ptmCode = ptm.getPtmCode();
				byteBuffer.add(DynByteBuffer.toByteArray(ptmCode));
			}
		}
		final byte b = -1;
		byteBuffer.add(b);

		final byte[] proteinIdB = DynByteBuffer.toByteArray(proteinId);
		byteBuffer.add(proteinIdB);

		byteBuffer.add(DBIndexStoreSQLiteByteIndexMerge.SEQ_SEPARATOR);
		return byteBuffer.getData();
	}
}
