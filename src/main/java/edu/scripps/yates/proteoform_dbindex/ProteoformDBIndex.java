package edu.scripps.yates.proteoform_dbindex;

import java.io.File;
import java.io.IOException;

import edu.scripps.yates.dbindex.DBIndexImpl;
import edu.scripps.yates.dbindex.io.DBIndexSearchParamsImpl;
import edu.scripps.yates.proteoform_dbindex.model.ExtendedAssignMass;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexSearchParams;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;
import edu.scripps.yates.utilities.masses.AssignMass;
import org.jdom.JDOMException;
import org.junit.Assert;

public class ProteoformDBIndex {
	public ProteoformDBIndex(File fastaFile, char[] enzymeArray, int missedCleavages, boolean semicleavage,
			File uniprotReleasesFolder) {
		final DBIndexSearchParams defaultDBIndexParams = DBIndexImpl.getDefaultDBIndexParams(fastaFile);

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
		final DBIndexImpl dbIndex = new DBIndexImpl(defaultDBIndexParams);
	}

	public static void main(String[] args) throws IOException, DBIndexStoreException, JDOMException {
		String pathToSearchXml = args[0];
		SearchParams sp = new SearchParams(pathToSearchXml);

		// write your code here


		String name = sp.getDatabaseName();
		File f = new File(name);

		String directory = f.getParentFile().getPath();
	//	sp.getMaxInternalMisCleavage();
		// test();

        /*String directory = "/home/yateslab/projectData/prolucid/3101ProteoformWork/small_search";
        String name = "/home/yateslab/projectData/prolucid/3101ProteoformWork/small_search/" +
                "UniProt_Human_reviewed_contaminant_04-02-2019_reversed.fasta";*/

		int maxMisclevages = sp.getMaxInternalMisCleavage();
		maxMisclevages = maxMisclevages == -1 ? Integer.MAX_VALUE : maxMisclevages;
		double minPrcMass = sp.getMinPrecursorMass();
		double maxPrcMass = sp.getMaxPrecursorMass();
	String residues = sp.getResidues();
		int maxNumVariationsPerPeptide = sp.getMaxAlter();
		int offset =  sp.getOffset();
		//  String uniprotVersion = null;
		String uniprotVersion = sp.getUniprotVersion();
		boolean semiCleavage = sp.getEnzymeSpecificity() != 2;
		boolean useMono = "mono".equals(sp.getPrecursorIsotope());
		DBIndexSearchParams params =  DBIndexImpl.getDefaultDBIndexParamsForProteoformAnalysis(name,
				maxMisclevages, minPrcMass,maxPrcMass, residues, null,offset,semiCleavage,
				uniprotVersion, null, directory);
		System.out.println(params.getDiscardDecoyRegexp());
		//AssignMass.setnTerm(sp.getStaticCTermMod());
		((DBIndexSearchParamsImpl) params).setInMemoryIndex(true);

		ProteoformDBIndexInterface proteoformDBIndex = new ProteoformDBIndexInterface(params,maxNumVariationsPerPeptide);
		ExtendedAssignMass extendedAssignMass = ExtendedAssignMass.getInstance(useMono,null);
		System.out.println(">>> Done");
	}

}
