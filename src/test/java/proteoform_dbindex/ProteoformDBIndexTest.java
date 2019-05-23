package proteoform_dbindex;

import static org.junit.Assert.fail;

import java.io.File;
import java.util.List;
import java.util.Set;

import org.junit.Before;
import org.junit.Test;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.dbindex.util.IndexUtil;
import edu.scripps.yates.proteoform_dbindex.ProteoformDBIndexInterface;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedProtein;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedSequence;

public class ProteoformDBIndexTest {

	private final String phosphorilatedPeptide = "ALYDFLPREPCNLALR";
	final char[] enzymeArray = { 'K', 'R' };
	final int missedCleavages = 3;
	final boolean semicleavage = false;
	private File paramFile;
	private File phosphoSiteDBFile;
	private String phosphoSiteSpecies;
	private ProteoformDBIndexInterface proteoformDBIndex;

	@Before
	public void loadResources() {
		// paramFile with input FASTA file, missedcleavages, etc.
		paramFile = new File("D:\\Salva\\git_projects\\proteoform_dbindex\\src\\main\\resources\\blazmass.params");

		// phosphoSiteDBFile
		// if you want to consider also www.phosphoSite.org phosphorilation annotations,
		// use this file.
		// this file is downloaded from https://www.phosphosite.org/staticDownloads in
		// the link named "Phosphosite_PTM_seq.fasta.gz"
		phosphoSiteDBFile = new File("Z:\\share\\Salva\\data\\PhosphositePTMseq\\Phosphosite_PTM_seq.fasta.gz");

		// phosphoSiteSpecies can be human, rat, cow, etc...
		phosphoSiteSpecies = "human";

		// uniprotReleasesFolder
		// location where uniprot annotations will be stored and indexed
		final File uniprotReleasesFolder = new File("Z:\\share\\Salva\\data\\uniprotKB");

		// uniprotVersion
		// null for latest version or "2019_03" for March 2019 version, for example
		final String uniprotVersion = null;

		// maxNumVariationsPerPeptide
		// maximum number of variations per peptide allowed
		final int maxNumVariationsPerPeptide = 4;

		// if provided, peptideInclusionList is a list of peptides to index (the rest
		// will be ignored). If not provided (null) all peptides under the parameters
		// will be indexed
		final Set<String> peptideInclusionList = null;

		// create index
		proteoformDBIndex = new ProteoformDBIndexInterface(paramFile, true, true, phosphoSiteSpecies, phosphoSiteDBFile,
				new UniprotProteinLocalRetriever(uniprotReleasesFolder, true), uniprotVersion,
				maxNumVariationsPerPeptide, peptideInclusionList);

	}

	@Test
	public void testingProteoformIndex() {

		try {
			final Set<IndexedProtein> proteins = proteoformDBIndex.getProteins(phosphorilatedPeptide);
			for (final IndexedProtein indexedProtein : proteins) {
				System.out.println(indexedProtein.getAccession());
				System.out.println(indexedProtein.getFastaDefLine());
			}
			final double parentMass = IndexUtil.calculateMass(phosphorilatedPeptide);
			final List<IndexedSequence> sequences = proteoformDBIndex.getSequences(parentMass, 0.0001);
			for (final IndexedSequence indexedSequence : sequences) {
				System.out.println(indexedSequence.getSequence() + "\t" + +indexedSequence.getMass());
			}

			final double phosphorilatedParentMass = parentMass + 79.966331;
			final List<IndexedSequence> phosphorilatedSequences = proteoformDBIndex
					.getSequences(phosphorilatedParentMass, 0.0001);
			for (final IndexedSequence indexedSequence : phosphorilatedSequences) {
				System.out.println(indexedSequence.getSequence() + "\t" + +indexedSequence.getMass());
			}

		} catch (final DBIndexStoreException e) {
			e.printStackTrace();
			fail();
		}

	}
}
