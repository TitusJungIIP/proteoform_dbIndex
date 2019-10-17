package proteoform_dbindex;

import java.io.File;
import java.util.List;

import org.junit.Test;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.dbindex.DBIndexImpl;
import edu.scripps.yates.proteoform_dbindex.ProteoformDBIndexInterface;
import edu.scripps.yates.proteoform_dbindex.model.ExtendedAssignMass;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexSearchParams;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedSequence;

public class TitusTest {
	private static final File uniprotReleasesFolder = new File("Z:\\share\\Salva\\data\\uniprotKB");

	@Test
	public void testTitus() throws DBIndexStoreException {
		// input parameters file
		final File paramFile = new File("d:\\titus\\blazmass_Q13523.params");

//maximum number of variations (sequence variations and PTMs) per peptide
		final int maxNumVariationsPerPeptide = 4;

//Uniprot repository version release
//null for latest version or "2019_05" for May 2019 version, for example
		final String uniprotVersion = null;

//Uniprot annotations retriever. It will retrieve the annotations to folder uniprotReleasesFolder
		final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
		// System.out.println("TESTTTING ");
//Create ProteoformDBIndexInterface instance
		final String sufix = null;
		final ProteoformDBIndexInterface proteoformDBIndex = new ProteoformDBIndexInterface(paramFile, sufix, uplr,
				uniprotVersion, maxNumVariationsPerPeptide);
		final List<IndexedSequence> sequences = proteoformDBIndex.getSequences(2275.1242204659998, 0.0001);
		System.out.println("PRINTING SEQ: ");
		for (final IndexedSequence iseq : sequences) {
			System.out.println(iseq.getSequence());
		}

	}

	@Test
	public void testTitus2() throws DBIndexStoreException {
		// input parameters file
		final File paramFile = new File("d:\\titus\\blazmass_UniProt_Human_sprot_11-08-2010_reversed.params");

//maximum number of variations (sequence variations and PTMs) per peptide
		final int maxNumVariationsPerPeptide = 0;

//Uniprot repository version release
//null for latest version or "2019_05" for May 2019 version, for example
		final String uniprotVersion = "2019_08";
//Uniprot annotations retriever. It will retrieve the annotations to folder uniprotReleasesFolder
		final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
		// System.out.println("TESTTTING ");
//Create ProteoformDBIndexInterface instance
		final String sufix = null;
		final ProteoformDBIndexInterface proteoformDBIndex = new ProteoformDBIndexInterface(paramFile, sufix, uplr,
				uniprotVersion, maxNumVariationsPerPeptide);
		final List<IndexedSequence> sequences = proteoformDBIndex.getSequences(2275.1242204659998, 0.0001);
		System.out.println("PRINTING SEQ: ");
		for (final IndexedSequence iseq : sequences) {
			System.out.println(iseq.getSequence());
		}

	}

	@Test
	public void testTitus3() throws DBIndexStoreException {

//maximum number of variations (sequence variations and PTMs) per peptide
		final int maxNumVariationsPerPeptide = 2;

//Uniprot repository version release
//null for latest version or "2019_05" for May 2019 version, for example
		final String uniprotVersion = "2019_08";
//Uniprot annotations retriever. It will retrieve the annotations to folder uniprotReleasesFolder
		final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
		final String fastaFilePath = "D:\\titus\\UniProt_Human_sprot_11-08-2010_reversed.fasta";
		final int maxMissedCleavages = 3;
		final double minPrecursorMass = 600;
		final double maxPrecursorMass = 6000;
		final String enzymeResidues = "KR";

		final String enzymeNocutResidues = null;
		final int enzymeOffset = 1;
		final String discardDecoyRegexp = null;
		final boolean semiCleavage = false;
		// System.out.println("TESTTTING ");
//Create ProteoformDBIndexInterface instance
		final DBIndexSearchParams params = DBIndexImpl.getDefaultDBIndexParamsForProteoformAnalysis(fastaFilePath,
				maxMissedCleavages, minPrecursorMass, maxPrecursorMass, enzymeResidues, enzymeNocutResidues,
				enzymeOffset, semiCleavage, uniprotVersion, discardDecoyRegexp,
				uniprotReleasesFolder.getAbsolutePath());
		final ProteoformDBIndexInterface proteoformDBIndex = new ProteoformDBIndexInterface(params,
				maxNumVariationsPerPeptide);
		final List<IndexedSequence> sequences = proteoformDBIndex.getSequences(2275.1242204659998, 0.0001);
		System.out.println("PRINTING SEQ: ");
		for (final IndexedSequence iseq : sequences) {
			System.out.println(iseq.getSequence());
		}

	}

	@Test
	public void test6() {

		final String directory = "d:\\titus";

		final int maxMisclevages = 3;
		final double minPrcMass = 600;
		final double maxPrcMass = 6000;
		final String residues = "KR";
		final int maxNumVariationsPerPeptide = 0;
		final int offset = 1;

		final String fastaFilePath = "D:\\titus\\UniProt_Human_sprot_11-08-2010_reversed.fasta";
		final String discarDecoys = null;
		final String uniprotVersion = "2019_08";
		final DBIndexSearchParams params = DBIndexImpl.getDefaultDBIndexParamsForProteoformAnalysis(fastaFilePath,
				maxMisclevages, minPrcMass, maxPrcMass, residues, null, offset, false, uniprotVersion, discarDecoys,
				uniprotReleasesFolder.getAbsolutePath());

		final ProteoformDBIndexInterface proteoformDBIndex = new ProteoformDBIndexInterface(params,
				maxNumVariationsPerPeptide);
		final Object extendedAssignMass = ExtendedAssignMass.getInstance(false, null);

	}
}
