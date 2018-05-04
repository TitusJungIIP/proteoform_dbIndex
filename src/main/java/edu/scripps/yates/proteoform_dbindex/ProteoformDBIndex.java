package edu.scripps.yates.proteoform_dbindex;

import java.io.File;

import edu.scripps.yates.dbindex.DBIndexInterface;
import edu.scripps.yates.dbindex.io.DBIndexSearchParams;
import edu.scripps.yates.dbindex.io.DBIndexSearchParamsImpl;

public class ProteoformDBIndex {
	public ProteoformDBIndex(File fastaFile, char[] enzymeArray, int missedCleavages, boolean semicleavage,
			File uniprotReleasesFolder) {
		final DBIndexSearchParams defaultDBIndexParams = DBIndexInterface.getDefaultDBIndexParams(fastaFile);

		((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeArr(enzymeArray, missedCleavages, semicleavage);
		((DBIndexSearchParamsImpl) defaultDBIndexParams).setSemiCleavage(semicleavage);
		((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeOffset(0);
		((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeNocutResidues("");
		((DBIndexSearchParamsImpl) defaultDBIndexParams).setH2OPlusProtonAdded(true);
		((DBIndexSearchParamsImpl) defaultDBIndexParams).setUniprotReleasesFolder(uniprotReleasesFolder);
		((DBIndexSearchParamsImpl) defaultDBIndexParams).setLookProteoforms(true);
		// if looking for proteoforms, not use in memory
		final boolean inMemoryIndex = false;
		((DBIndexSearchParamsImpl) defaultDBIndexParams).setInMemoryIndex(inMemoryIndex);
		final DBIndexInterface dbIndex = new DBIndexInterface(defaultDBIndexParams);
	}
}
