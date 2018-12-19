package edu.scripps.yates.proteoform_dbindex;

import java.io.File;
import java.io.IOException;

import edu.scripps.yates.dbindex.Constants;
import edu.scripps.yates.dbindex.DBIndexStoreSQLiteMult;
import edu.scripps.yates.proteoform_dbindex.model.ExtendedAssignMass;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexSearchParams;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;

public class ProteoformDBIndexStoreSQLiteMult extends DBIndexStoreSQLiteMult {
	private static final org.apache.log4j.Logger logger = org.apache.log4j.Logger
			.getLogger(ProteoformDBIndexStoreSQLiteMult.class);
	private final ExtendedAssignMass extendedAssignMass;

	public ProteoformDBIndexStoreSQLiteMult(DBIndexSearchParams sparam, boolean inMemoryIndex,
			ExtendedAssignMass extendedAssignMass) {
		super(sparam, inMemoryIndex, new ProteoformDBIndexStoreSQLiteByteIndexMerge[sparam.getIndexFactor()]);
		this.extendedAssignMass = extendedAssignMass;
	}

	@Override
	public void init(String databaseID) throws DBIndexStoreException {
		if (databaseID == null || databaseID.equals("")) {
			throw new DBIndexStoreException("Index path is missing, cannot initialize the indexer.");
		}

		if (inited) {
			throw new DBIndexStoreException("Already intialized");
		}

		final File dbBase = new File(databaseID).getAbsoluteFile();
		String baseName = dbBase.getName();
		if (!baseName.endsWith(IDX_SUFFIX)) {
			baseName = baseName + IDX_SUFFIX;
		}
		final File parentDir = dbBase.getParentFile();

		final File indexDir = new File(parentDir.getAbsolutePath() + File.separator + baseName);
		if (indexDir.exists() && !indexDir.isDirectory()) {
			logger.info("Trying to delete old index file: " + indexDir.getAbsolutePath());
			indexDir.delete();

		}

		try {
			indexDir.mkdir();
		} catch (final SecurityException e) {
		}

		dbPathBase = indexDir.getAbsolutePath();

		logger.info("Using database index dir: " + dbPathBase);

		if (!indexDir.exists()) {
			throw new DBIndexStoreException(
					"Index dir does not exist: " + indexDir + ", cannot initialize the indexer.");
		}

		// initialize buckets
		for (int i = 0; i < Constants.NUM_BUCKETS; ++i) {
			// buckets[i] = new DBIndexStoreSQLiteByte(i); //version with merge
			// during search
			buckets[i] = new ProteoformDBIndexStoreSQLiteByteIndexMerge(sparam, inMemoryIndex, i,
					(ProteoformProteinCache) proteinCache, extendedAssignMass); // version
			// with
			// merge
			// during
			// index
			final String bucketId = dbPathBase + File.separator + i + IDX_SUFFIX;
			try {
				buckets[i].init(bucketId);
			} catch (final DBIndexStoreException e) {
				logger.error("Error initializing bucket " + i, e);
				throw e;
			}
		}

		inited = true;

	}

	public void setProteinCache(ProteoformProteinCache proteinCache) {
		this.proteinCache = proteinCache;
		for (int i = 0; i < Constants.NUM_BUCKETS; ++i) {
			buckets[i].setProteinCache(proteinCache);
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
