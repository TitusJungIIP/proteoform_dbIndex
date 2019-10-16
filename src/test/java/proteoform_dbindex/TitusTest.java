package proteoform_dbindex;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.annotations.uniprot.UniprotProteinRemoteRetriever;
import edu.scripps.yates.dbindex.DBIndexImpl;
import edu.scripps.yates.proteoform_dbindex.ProteoformDBIndexInterface;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotEntryUtil;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexSearchParams;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedSequence;
import gnu.trove.set.hash.THashSet;

public class TitusTest {
	private static final File uniprotReleasesFolder = new File("Z:\\share\\Salva\\data\\uniprotKB");

	@Test
	public void testTitus() throws DBIndexStoreException {
		// input parameters file
		final File paramFile = new File("C:\\Users\\salvador\\Desktop\\tmp\\titus\\blazmass_Q13523.params");

//maximum number of variations (sequence variations and PTMs) per peptide
		final int maxNumVariationsPerPeptide = 4;

//Uniprot repository version release
//null for latest version or "2019_05" for May 2019 version, for example
		final String uniprotVersion = null;
		final File uniprotReleasesFolder = new File("C:\\Users\\salvador\\Desktop\\tmp\\titus");
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
	public void testTitus3() throws DBIndexStoreException {
		// input parameters file
		final File paramFile = new File("d:\\titus\\blazmass_UniProt_Human_sprot_11-08-2010_reversed.params");

//maximum number of variations (sequence variations and PTMs) per peptide
		final int maxNumVariationsPerPeptide = 0;

//Uniprot repository version release
//null for latest version or "2019_05" for May 2019 version, for example
		final String uniprotVersion = null;
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
	public void test4() {
		final Set<String> accessions = new THashSet<String>();
		accessions.add("Q00005-2");
		accessions.add("Q00005-3");
		final Map<String, Entry> entries = UniprotProteinRemoteRetriever.getFASTASequencesInParallel(accessions);
		System.out.println(entries.size() + " entries");
		for (final String acc : accessions) {
			System.out.println(entries.containsKey(acc) + " " + acc);
			final Entry entry = entries.get(acc);
			System.out.println("\n\n==============================\n" + acc);
			System.out.println(UniprotEntryUtil.getPrimaryAccession(entry));
			System.out.println(UniprotEntryUtil.getProteinSequence(entry));
		}
	}

	@Test
	public void test5() {
		final Set<String> accessions = new THashSet<String>();
		accessions.add("Q00005-2");
		accessions.add("Q00005-3");
		final Map<String, Entry> entries = new UniprotProteinRemoteRetriever()
				.getIsoformFASTASequencesFromUniprot(accessions, uniprotReleasesFolder, null);
		System.out.println(entries.size() + " entries");
		for (final String acc : accessions) {
			System.out.println(entries.containsKey(acc) + " " + acc);
			final Entry entry = entries.get(acc);
			System.out.println("\n\n==============================\n" + acc);
			System.out.println(UniprotEntryUtil.getPrimaryAccession(entry));
			System.out.println(UniprotEntryUtil.getProteinSequence(entry));
		}
	}
}
