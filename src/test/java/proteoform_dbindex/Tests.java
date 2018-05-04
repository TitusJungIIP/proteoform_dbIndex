package proteoform_dbindex;

import static org.junit.Assert.fail;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.springframework.core.io.ClassPathResource;

import edu.scripps.yates.dbindex.DBIndexInterface;
import edu.scripps.yates.dbindex.IndexedProtein;
import edu.scripps.yates.dbindex.IndexedSequence;
import edu.scripps.yates.dbindex.io.DBIndexSearchParams;
import edu.scripps.yates.dbindex.io.DBIndexSearchParamsImpl;
import edu.scripps.yates.dbindex.model.AssignMass;
import edu.scripps.yates.dbindex.util.IndexUtil;
import edu.scripps.yates.proteoform_dbindex.ProteoformDBIndexInterface;
import edu.scripps.yates.proteoform_dbindex.model.ExtendedAssignMass;

public class Tests {

	@Test
	public void testNormalIndex() {
		File fastaFile;
		try {
			fastaFile = new ClassPathResource("P04637.fasta").getFile();
			final char[] enzymeArray = { 'K', 'R' };
			final int missedCleavages = 3;
			final boolean semicleavage = false;
			final File uniprotReleasesFolder = new File("Z:\\share\\Salva\\data\\uniprotKB");
			final DBIndexSearchParams defaultDBIndexParams = DBIndexInterface.getDefaultDBIndexParams(fastaFile);

			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeArr(enzymeArray, missedCleavages, semicleavage);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setSemiCleavage(semicleavage);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeOffset(0);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeNocutResidues("");
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setH2OPlusProtonAdded(true);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setUniprotReleasesFolder(uniprotReleasesFolder);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setLookProteoforms(false);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setUseMonoParent(true);

			// if looking for proteoforms, not use in memory
			final boolean inMemoryIndex = false;
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setInMemoryIndex(inMemoryIndex);
			final DBIndexInterface dbIndex = new DBIndexInterface(defaultDBIndexParams);
			final Set<IndexedProtein> proteins = dbIndex.getProteins("MEEPQSDPSVEPPLSQETFSDLWK");
			for (final IndexedProtein indexedProtein : proteins) {
				System.out.println(indexedProtein.getAccession());
			}
		} catch (final IOException e) {
			e.printStackTrace();
			fail();
		}

	}

	@Test
	public void testProteoformIndex() throws IOException {
		File fastaFile;
		// fastaFile = new ClassPathResource("P04637.fasta").getFile();
		fastaFile = new File("D:\\Downloads\\casimir\\Uniprot_Human.fasta");
		final char[] enzymeArray = { 'K', 'R' };
		final int missedCleavages = 3;
		final boolean semicleavage = false;
		final File uniprotReleasesFolder = new File("Z:\\share\\Salva\\data\\uniprotKB");
		final DBIndexSearchParams defaultDBIndexParams = DBIndexInterface.getDefaultDBIndexParams(fastaFile);
		final int numVariationsPerPeptide = 5;
		final boolean useUniprot = true;
		final boolean usePhosphoSite = false;
		((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeArr(enzymeArray, missedCleavages, semicleavage);
		((DBIndexSearchParamsImpl) defaultDBIndexParams).setSemiCleavage(semicleavage);
		((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeOffset(0);
		((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeNocutResidues("");
		((DBIndexSearchParamsImpl) defaultDBIndexParams).setH2OPlusProtonAdded(true);
		((DBIndexSearchParamsImpl) defaultDBIndexParams).setUniprotReleasesFolder(uniprotReleasesFolder);
		((DBIndexSearchParamsImpl) defaultDBIndexParams).setLookProteoforms(true);
		((DBIndexSearchParamsImpl) defaultDBIndexParams).setUseMonoParent(true);
		// if looking for proteoforms, not use in memory
		final boolean inMemoryIndex = false;
		((DBIndexSearchParamsImpl) defaultDBIndexParams).setInMemoryIndex(inMemoryIndex);

		final ProteoformDBIndexInterface dbIndex = new ProteoformDBIndexInterface(defaultDBIndexParams, useUniprot,
				usePhosphoSite, "human", uniprotReleasesFolder, null, numVariationsPerPeptide);
		final String peptideSeq = "MEEPQSDPSVEPPLSQETFSDLWK";
		final Set<IndexedProtein> proteins = dbIndex.getProteins(peptideSeq);
		Assert.assertFalse(proteins.isEmpty());
		for (final IndexedProtein indexedProtein : proteins) {
			System.out.println("Proteins of peptide " + "MEEPQSDPSVEPPLSQETFSDLWK\t" + indexedProtein.getAccession());
		}
		final double calculateMass = IndexUtil.calculateMass(peptideSeq, defaultDBIndexParams.isH2OPlusProtonAdded());
		final List<IndexedSequence> sequences = dbIndex.getSequences(calculateMass, 0.1);
		for (final IndexedSequence indexedSequence : sequences) {
			System.out.println("Peptide with mass close to " + calculateMass + "\t" + indexedSequence.getSequence()
					+ "\t" + indexedSequence.getMass());
			double x = ExtendedAssignMass.getInstance(true, null).calculateMass(indexedSequence.getSequence());
			x += AssignMass.H2O_PROTON;
			System.out.println(x + "\tvs\t" + indexedSequence.getMass());

		}
		final double precursorMass = 2766.26 + 80;
		final List<IndexedSequence> sequences2 = dbIndex.getSequences(precursorMass, 0.1);
		for (final IndexedSequence indexedSequence : sequences2) {
			System.out.println("Peptide with mass close to " + precursorMass + "\t" + indexedSequence.getSequence()
					+ "\t" + indexedSequence.getMass());
			double x = ExtendedAssignMass.getInstance(true, null).calculateMass(indexedSequence.getSequence());
			x += AssignMass.H2O_PROTON;
			System.out.println(x + "\tvs\t" + indexedSequence.getMass());
		}

	}

	@Before
	public void before() {
		new AssignMass(true);
	}

	@Test
	public void test() {
		System.out.println(AssignMass.NH3);
		System.out.println(AssignMass.H2O);
		System.out.println(AssignMass.CO);

	}

	@Test
	public void printShorts() {
		for (short i = 0; i < 1000; i++) {
			System.out.println(i);
		}
	}
}
