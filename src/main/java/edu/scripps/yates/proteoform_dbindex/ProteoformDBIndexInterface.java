package edu.scripps.yates.proteoform_dbindex;

import java.io.File;
import java.io.IOException;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.dbindex.DBIndexInterface;
import edu.scripps.yates.dbindex.DBIndexer.IndexerMode;
import edu.scripps.yates.dbindex.DBIndexerException;
import edu.scripps.yates.dbindex.SearchParams;
import edu.scripps.yates.dbindex.io.SearchParamReader;
import edu.scripps.yates.dbindex.model.DBIndexSearchParams;
import edu.scripps.yates.utilities.masses.AssignMass;

public class ProteoformDBIndexInterface extends DBIndexInterface {
	private final static Logger log = Logger.getLogger(ProteoformDBIndexInterface.class);

	public ProteoformDBIndexInterface(File paramFile, boolean useUniprot, boolean usePhosphosite,
			String phosphoSiteSpecies, UniprotProteinLocalRetriever uplr, String uniprotVersion,
			int maxNumVariationsPerPeptide) {
		super();
		final String paramFileName = SearchParamReader.DEFAULT_PARAM_FILE_NAME;

		try {

			final String dbIndexPath = getDBIndexPath();
			SearchParamReader pr = null;
			if (paramFile != null) {
				pr = new SearchParamReader(paramFile);
			} else {
				pr = new SearchParamReader(dbIndexPath, paramFileName);
			}

			final SearchParams sParam = pr.getSearchParams();
			final edu.scripps.yates.dbindex.DBIndexer.IndexerMode indexerMode = sParam.isUseIndex()
					? IndexerMode.SEARCH_INDEXED : IndexerMode.SEARCH_UNINDEXED;
			indexer = new ProteoformDBIndexer(sParam, indexerMode, useUniprot, usePhosphosite, phosphoSiteSpecies, uplr,
					uniprotVersion, maxNumVariationsPerPeptide);
			try {
				indexer.init();

			} catch (final DBIndexerException ex) {
				indexer = new ProteoformDBIndexer(sParam, IndexerMode.INDEX, useUniprot, usePhosphosite,
						phosphoSiteSpecies, uplr, uniprotVersion, maxNumVariationsPerPeptide);
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

	public ProteoformDBIndexInterface(DBIndexSearchParams sParam, boolean useUniprot, boolean usePhosphosite,
			String phosphoSiteSpecies, UniprotProteinLocalRetriever uplr, String uniprotVersion,
			int maxNumVariationsPerPeptide) {
		super();
		try {

			// init the masses
			AssignMass.getInstance(sParam.isUseMonoParent());
			final edu.scripps.yates.dbindex.DBIndexer.IndexerMode indexerMode = sParam.isUseIndex()
					? IndexerMode.SEARCH_INDEXED : IndexerMode.SEARCH_UNINDEXED;
			indexer = new ProteoformDBIndexer(sParam, indexerMode, useUniprot, usePhosphosite, uniprotVersion, uplr,
					uniprotVersion, maxNumVariationsPerPeptide);
			try {
				indexer.init();

			} catch (final DBIndexerException ex) {
				indexer = new ProteoformDBIndexer(sParam, IndexerMode.INDEX, useUniprot, usePhosphosite, uniprotVersion,
						uplr, uniprotVersion, maxNumVariationsPerPeptide);
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
