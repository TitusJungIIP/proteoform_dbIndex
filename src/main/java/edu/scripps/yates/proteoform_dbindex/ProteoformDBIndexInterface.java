package edu.scripps.yates.proteoform_dbindex;

import java.io.File;
import java.io.IOException;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.dbindex.DBIndexImpl;
import edu.scripps.yates.dbindex.DBIndexer;
import edu.scripps.yates.dbindex.DBIndexer.IndexerMode;
import edu.scripps.yates.dbindex.DBIndexerException;
import edu.scripps.yates.dbindex.SearchParams;
import edu.scripps.yates.dbindex.io.SearchParamReader;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexSearchParams;
import edu.scripps.yates.utilities.masses.AssignMass;

public class ProteoformDBIndexInterface extends DBIndexImpl {
	private final static Logger log = Logger.getLogger(ProteoformDBIndexInterface.class);

	/**
	 * 
	 * @param sParam                     the parameters object. You can create it
	 *                                   from the static method:
	 *                                   {@link DBIndexImpl}.getDefaultDBIndexParamsForProteoformAnalysis()
	 * 
	 * @param maxNumVariationsPerPeptide
	 */
	public ProteoformDBIndexInterface(DBIndexSearchParams sParam, int maxNumVariationsPerPeptide) {
		super();

		try {
			final boolean useUniprot = true;
			final boolean usePhosphosite = false;
			final String phosphoSiteSpecies = null;
			final File phosphoSiteDBFile = null;
			final Set<String> peptideInclusionList = null;

			final edu.scripps.yates.dbindex.DBIndexer.IndexerMode indexerMode = sParam.isUseIndex()
					? IndexerMode.SEARCH_INDEXED
					: IndexerMode.SEARCH_UNINDEXED;
			final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(
					sParam.getUniprotReleasesFolder(), true);
			final String uniprotVersion = sParam.getUniprotVersion();

			if (!sParam.isLookProteoforms()) {
				throw new IllegalArgumentException(
						"lookForProteoforms is FALSE. Set it to TRUE when creating DBIndexSearchParams object");
			}

			indexer = createIndexer(sParam, indexerMode, useUniprot, usePhosphosite, phosphoSiteSpecies,
					phosphoSiteDBFile, uplr, uniprotVersion, maxNumVariationsPerPeptide, peptideInclusionList);
			try {
				indexer.init();

			} catch (final DBIndexerException ex) {
				indexer = createIndexer(sParam, IndexerMode.INDEX, useUniprot, usePhosphosite, phosphoSiteSpecies,
						phosphoSiteDBFile, uplr, uniprotVersion, maxNumVariationsPerPeptide, peptideInclusionList);
				try {
					indexer.init();
					indexer.run();
				} catch (final DBIndexerException e) {
					e.printStackTrace();
					log.error("Could not initialize the indexer in search mode and init the worker thread");
				}
			}

		} catch (final IOException e) {
			e.printStackTrace();
			log.error(e.getMessage());

		} catch (final Exception e) {
			e.printStackTrace();
			log.error(e.getMessage());
		}
	}

	public ProteoformDBIndexInterface(File paramFile, UniprotProteinLocalRetriever uplr, String uniprotVersion,
			int maxNumVariationsPerPeptide) {
		this(paramFile, true, false, null, null, uplr, uniprotVersion, maxNumVariationsPerPeptide, null);
	}

	public ProteoformDBIndexInterface(File paramFile, boolean useUniprot, boolean usePhosphosite,
			String phosphoSiteSpecies, File phosphoSiteDBFile, UniprotProteinLocalRetriever uplr, String uniprotVersion,
			int maxNumVariationsPerPeptide, Set<String> peptideInclusionList) {
		super();
		final String paramFileName = SearchParamReader.DEFAULT_PARAM_FILE_NAME;

		try {

			SearchParamReader pr = null;
			if (paramFile != null) {
				pr = new SearchParamReader(paramFile);
			} else {
				final String dbIndexPath = getDBIndexPath();
				pr = new SearchParamReader(dbIndexPath, paramFileName);
			}

			final SearchParams sParam = pr.getSearchParams();
			final edu.scripps.yates.dbindex.DBIndexer.IndexerMode indexerMode = sParam.isUseIndex()
					? IndexerMode.SEARCH_INDEXED
					: IndexerMode.SEARCH_UNINDEXED;
			if (uplr != null && sParam.getUniprotReleasesFolder() == null) {
				sParam.setUniprotReleasesFolder(uplr.getUniprotReleasesFolder());
			}
			sParam.setLookProteoforms(true);
			indexer = createIndexer(sParam, indexerMode, useUniprot, usePhosphosite, phosphoSiteSpecies,
					phosphoSiteDBFile, uplr, uniprotVersion, maxNumVariationsPerPeptide, peptideInclusionList);
			try {
				indexer.init();

			} catch (final DBIndexerException ex) {
				indexer = createIndexer(sParam, IndexerMode.INDEX, useUniprot, usePhosphosite, phosphoSiteSpecies,
						phosphoSiteDBFile, uplr, uniprotVersion, maxNumVariationsPerPeptide, peptideInclusionList);
				try {
					indexer.init();
					indexer.run();
				} catch (final DBIndexerException e) {
					e.printStackTrace();
					log.error("Could not initialize the indexer in search mode and init the worker thread");
				}
			}

			dbIndexByFile.put(paramFile, this);
		} catch (final IOException e) {
			e.printStackTrace();
			log.error(e.getMessage());

		} catch (final Exception e) {
			e.printStackTrace();
			log.error(e.getMessage());
		}
	}

	private DBIndexer createIndexer(DBIndexSearchParams sParam, IndexerMode indexerMode, boolean useUniprot,
			boolean usePhosphosite, String phosphoSiteSpecies, File phosphoSiteDBFile,
			UniprotProteinLocalRetriever uplr, String uniprotVersion, int maxNumVariationsPerPeptide,
			Set<String> peptideInclusionList) throws IOException {
		return new ProteoformDBIndexer(sParam, indexerMode, useUniprot, usePhosphosite, phosphoSiteSpecies,
				phosphoSiteDBFile, uplr, uniprotVersion, maxNumVariationsPerPeptide, peptideInclusionList);
	}

	public ProteoformDBIndexInterface(DBIndexSearchParams sParam, boolean useUniprot, boolean usePhosphosite,
			String phosphoSiteSpecies, File phosphoSiteDBFile, UniprotProteinLocalRetriever uplr, String uniprotVersion,
			int maxNumVariationsPerPeptide, Set<String> peptideInclusionList) {
		super();
		try {

			// init the masses
			AssignMass.getInstance(sParam.isUseMonoParent());
			final edu.scripps.yates.dbindex.DBIndexer.IndexerMode indexerMode = sParam.isUseIndex()
					? IndexerMode.SEARCH_INDEXED
					: IndexerMode.SEARCH_UNINDEXED;
			indexer = new ProteoformDBIndexer(sParam, indexerMode, useUniprot, usePhosphosite, phosphoSiteSpecies,
					phosphoSiteDBFile, uplr, uniprotVersion, maxNumVariationsPerPeptide, peptideInclusionList);
			try {
				indexer.init();

			} catch (final DBIndexerException ex) {
				indexer = new ProteoformDBIndexer(sParam, IndexerMode.INDEX, useUniprot, usePhosphosite,
						phosphoSiteSpecies, phosphoSiteDBFile, uplr, uniprotVersion, maxNumVariationsPerPeptide,
						peptideInclusionList);
				try {
					indexer.init();
					indexer.run();
				} catch (final DBIndexerException e) {
					e.printStackTrace();
					log.error("Could not initialize the indexer in search mode and init the worker thread");
				}
			}
			dbIndexByParamKey.put(sParam.getFullIndexFileName(), this);
		} catch (final IOException e) {
			e.printStackTrace();
			log.error(e.getMessage());

		} catch (final Exception e) {
			e.printStackTrace();
			log.error(e.getMessage());

		}

	}

}
