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

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.dbindex.util.IndexUtil;
import edu.scripps.yates.proteoform_dbindex.ProteoformDBIndexInterface;
import edu.scripps.yates.proteoform_dbindex.model.IndexedSequenceWithPTMs;
import edu.scripps.yates.proteoform_dbindex.model.PTM;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedProtein;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedSequence;

public class ProteoformDBIndexTest {

//	private File phosphoSiteDBFile;
//	private String phosphoSiteSpecies;
//	// uniprotReleasesFolder
	// location where uniprot annotations will be stored and indexed
	private final File uniprotReleasesFolder = new File(
			"D:\\Salva\\git_projects\\proteoform_dbindex\\src\\test\\resources\\uniprotKB");

	// uniprotVersion
	// null for latest version or "2019_05" for May 2019 version, for example
	private final String uniprotVersion = null;// "2019_05";

	// maxNumVariationsPerPeptide
	// maximum number of variations per peptide allowed
	private final int maxNumVariationsPerPeptide = 4;

	// if provided, peptideInclusionList is a list of peptides to index (the rest
	// will be ignored). If not provided (null) all peptides under the parameters
	// will be indexed
//	private final Set<String> peptideInclusionList = null;

	@Before
	public void loadResources() throws IOException {

		// phosphoSiteDBFile
		// if you want to consider also www.phosphoSite.org phosphorilation annotations,
		// use this file.
		// this file is downloaded from https://www.phosphosite.org/staticDownloads in
		// the link named "Phosphosite_PTM_seq.fasta.gz"
//		phosphoSiteDBFile = new ClassPathResource("Phosphosite_PTM_seq.fasta").getFile();

		// phosphoSiteSpecies can be human, rat, cow, etc...
//		phosphoSiteSpecies = "human";

	}

	@Test
	public void testingProteoformIndex_UniprotSwissProtHuman() throws IOException {
		// paramFile with input FASTA file, missedcleavages, etc.
		// this file is pointing to P42681.fasta file
		// change the paths correspondingly
		final File paramFile = new ClassPathResource("blazmass_uniprot_swissprot_human.params").getFile();
		// create index
		final String sufix = null;
		final ProteoformDBIndexInterface proteoformDBIndex = new ProteoformDBIndexInterface(paramFile, sufix,
				new UniprotProteinLocalRetriever(uniprotReleasesFolder, true), uniprotVersion,
				maxNumVariationsPerPeptide);
		try {
			System.out.println();
			System.out.println("Looking for the proteins from peptide " + "ALYDFLPR");
			final Set<IndexedProtein> proteins = proteoformDBIndex.getProteins("ALYDFLPR");
			for (final IndexedProtein indexedProtein : proteins) {
				System.out.println(indexedProtein.getAccession());
				System.out.println(indexedProtein.getFastaDefLine());
			}
			double parentMass = IndexUtil.calculateMass("ALYDFLPR");
			System.out.println();
			System.out.println("Looking for the mass " + parentMass + " that is from peptide " + "ALYDFLPR");
			List<IndexedSequence> sequences = proteoformDBIndex.getSequences(parentMass, 0.0001);
			for (final IndexedSequence indexedSequence : sequences) {
				System.out.println(
						indexedSequence.getSequence() + "\t" + IndexUtil.calculateMass(indexedSequence.getSequence())
								+ "\t" + indexedSequence.getMass() + "\t" + parentMass);
			}

			final double phosphorilatedParentMass = parentMass + 79.966331;
			System.out.println();
			System.out.println("Looking for the mass " + phosphorilatedParentMass + " that is from peptide "
					+ "ALYDFLPR" + " plus a phosphorilation");
			final List<IndexedSequence> phosphorilatedSequences = proteoformDBIndex
					.getSequences(phosphorilatedParentMass, 0.0001);
			Assert.assertFalse(phosphorilatedSequences.isEmpty());
			for (final IndexedSequence indexedSequence : phosphorilatedSequences) {
				if (indexedSequence instanceof IndexedSequenceWithPTMs) {
					final IndexedSequenceWithPTMs indexedSequenceWithPTM = (IndexedSequenceWithPTMs) indexedSequence;
					System.out.println(indexedSequenceWithPTM.getSequence() + "\t"
							+ IndexUtil.calculateMass(indexedSequence.getSequence()) + "\t"
							+ indexedSequenceWithPTM.getModSequence() + "\t" + indexedSequenceWithPTM.getMass() + "\t"
							+ phosphorilatedParentMass);
				} else {
					System.out.println(indexedSequence.getSequence() + "\t" + indexedSequence.getModSequence() + "\t"
							+ +indexedSequence.getMass() + "\t" + phosphorilatedParentMass);
				}
			}

			parentMass = IndexUtil.calculateMass("TQISLSTDEELPEKYTQRR");
			System.out.println();
			System.out.println("Looking for the mass " + parentMass + " that is from peptide " + "TQISLSTDEELPEKYTQRR");
			sequences = proteoformDBIndex.getSequences(parentMass, 0.0001);
			for (final IndexedSequence indexedSequence : sequences) {
				System.out.println(indexedSequence.getSequence() + "\t" + indexedSequence.getModSequence() + "\t"
						+ +indexedSequence.getMass() + "\t" + parentMass);
			}
			Assert.assertFalse(sequences.isEmpty());
			parentMass = IndexUtil.calculateMass("TQISLSTDEELPEKYTQHR");
			System.out.println();
			System.out.println("Looking for the mass " + parentMass + " that is from peptide " + "TQISLSTDEELPEKYTQHR");
			sequences = proteoformDBIndex.getSequences(parentMass, 0.0001);
			for (final IndexedSequence indexedSequence : sequences) {
				System.out.println(indexedSequence.getSequence() + "\t"
						+ IndexUtil.calculateMass(indexedSequence.getSequence()) + "\t"
						+ indexedSequence.getModSequence() + "\t" + +indexedSequence.getMass() + "\t" + parentMass);
				final List<Integer> proteinIds = indexedSequence.getProteinIds();
				for (final Integer proteinId : proteinIds) {
					final IndexedProtein indexedProtein = proteoformDBIndex.getIndexedProteinById(proteinId);
					System.out.println("protein id: " + proteinId + "\t" + indexedProtein.getAccession() + "\t"
							+ indexedProtein.getFastaDefLine());
				}
			}
			Assert.assertFalse(sequences.isEmpty());

		} catch (final Exception e) {
			e.printStackTrace();
			fail();
		} catch (final DBIndexStoreException e) {
			e.printStackTrace();
			fail();
		}

	}

	@Test
	public void testingProteoformIndex_Q13523() throws IOException {
		// paramFile with input FASTA file, missedcleavages, etc.
		// this file is pointing to P42681.fasta file
		// change the paths correspondingly
		final File paramFile = new ClassPathResource("blazmass_Q13523.params").getFile();
		// create index
		final String sufix = null;
		final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
		final ProteoformDBIndexInterface proteoformDBIndex = new ProteoformDBIndexInterface(paramFile, sufix, uplr,
				uniprotVersion, maxNumVariationsPerPeptide);

		try {
			System.out.println();
			System.out.println("Looking for the proteins from peptide " + "SPSPDDILERVAADVKEYER");
			Set<IndexedProtein> proteins = proteoformDBIndex.getProteins("SPSPDDILERVAADVKEYER");
			for (final IndexedProtein indexedProtein : proteins) {
				System.out.println(indexedProtein.getAccession());
				System.out.println(indexedProtein.getFastaDefLine());
			}
			System.out.println();
			System.out.println("Looking for the proteins from peptide with natural variance " + "SPSPDDVLERVAADVKEYER");
			proteins = proteoformDBIndex.getProteins("SPSPDDVLERVAADVKEYER");
			for (final IndexedProtein indexedProtein : proteins) {
				System.out.println(indexedProtein.getAccession());
				System.out.println(indexedProtein.getFastaDefLine());
			}
			final double parentMass = IndexUtil.calculateMass("SPSPDDVLERVAADVKEYER");
			System.out.println();
			System.out.println("Looking for the mass " + parentMass + " that is from peptide with natural variance "
					+ "SPSPDDVLERVAADVKEYER");
			final List<IndexedSequence> sequences = proteoformDBIndex.getSequences(parentMass, 0.0001);
			for (final IndexedSequence indexedSequence : sequences) {
				System.out.println(indexedSequence.getSequence() + "\t"
						+ IndexUtil.calculateMass(indexedSequence.getSequence()) + "\t"
						+ indexedSequence.getModSequence() + "\t" + indexedSequence.getMass() + "\t" + parentMass);
			}

			final double phosphorilatedParentMass = parentMass + 79.966331;
			System.out.println();
			System.out.println(
					"Looking for the mass " + phosphorilatedParentMass + " that is from peptide with natural variance "
							+ "SPSPDDVLERVAADVKEYER" + " plus a phosphorilation");
			final List<IndexedSequence> phosphorilatedSequences = proteoformDBIndex
					.getSequences(phosphorilatedParentMass, 0.0001);
			Assert.assertFalse(phosphorilatedSequences.isEmpty());
			for (final IndexedSequence indexedSequence : phosphorilatedSequences) {
				if (indexedSequence instanceof IndexedSequenceWithPTMs) {
					final IndexedSequenceWithPTMs indexedSequenceWithPTM = (IndexedSequenceWithPTMs) indexedSequence;
					System.out.println(indexedSequenceWithPTM.getSequence() + "\t"
							+ IndexUtil.calculateMass(indexedSequenceWithPTM.getSequence()) + "\t"
							+ indexedSequenceWithPTM.getModSequence() + "\t" + indexedSequenceWithPTM.getMass() + "\t"
							+ phosphorilatedParentMass);
				} else {
					System.out.println(
							indexedSequence.getSequence() + IndexUtil.calculateMass(indexedSequence.getSequence())
									+ "\t" + "\t" + indexedSequence.getModSequence() + "\t" + +indexedSequence.getMass()
									+ "\t" + phosphorilatedParentMass);
				}
			}

			final double doublyPhosphorilatedParentMass = parentMass + 79.966331 + 79.966331;
			System.out.println();
			System.out.println("Looking for the mass " + doublyPhosphorilatedParentMass
					+ " that is from peptide with natural variance" + "SPSPDDVLERVAADVKEYER"
					+ " plus 2 phosphorilations");
			final List<IndexedSequence> doublePhosphorilatedSequences = proteoformDBIndex
					.getSequences(doublyPhosphorilatedParentMass, 0.0001);
			Assert.assertFalse(doublePhosphorilatedSequences.isEmpty());
			for (final IndexedSequence indexedSequence : doublePhosphorilatedSequences) {
				if (indexedSequence instanceof IndexedSequenceWithPTMs) {
					final IndexedSequenceWithPTMs indexedSequenceWithPTM = (IndexedSequenceWithPTMs) indexedSequence;
					System.out.println(indexedSequenceWithPTM.getSequence() + "\t"
							+ IndexUtil.calculateMass(indexedSequenceWithPTM.getSequence()) + "\t"
							+ indexedSequenceWithPTM.getModSequence() + "\t" + indexedSequenceWithPTM.getMass() + "\t"
							+ doublyPhosphorilatedParentMass);
				} else {
					System.out.println(
							indexedSequence.getSequence() + IndexUtil.calculateMass(indexedSequence.getSequence())
									+ "\t" + "\t" + indexedSequence.getModSequence() + "\t" + +indexedSequence.getMass()
									+ "\t" + phosphorilatedParentMass);
				}
			}

		} catch (final Exception e) {
			e.printStackTrace();
			fail();
		} catch (final DBIndexStoreException e) {
			e.printStackTrace();
			fail();
		}

	}

	@Test
	public void testingProteoformIndex_UniprotMouse_with_contaminants() throws IOException {
		// paramFile with input FASTA file, missedcleavages, etc.
		// this file is pointing to P42681.fasta file
		// change the paths correspondingly
		final File paramFile = new ClassPathResource("blazmass_uniprot_mouse.params").getFile();
		// create index
		final String sufix = null;
		final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
		final ProteoformDBIndexInterface proteoformDBIndex = new ProteoformDBIndexInterface(paramFile, sufix, uplr,
				uniprotVersion, maxNumVariationsPerPeptide);

		try {
			System.out.println();
			System.out.println("Looking for the proteins from peptide " + "MACARPLISVYSEK");
			final Set<IndexedProtein> proteins = proteoformDBIndex.getProteins("MACARPLISVYSEK");
			for (final IndexedProtein indexedProtein : proteins) {
				System.out.println(indexedProtein.getAccession());
				System.out.println(indexedProtein.getFastaDefLine());
			}

		} catch (final Exception e) {
			e.printStackTrace();
			fail();
		} catch (final DBIndexStoreException e) {
			e.printStackTrace();
			fail();
		}

	}

	@Test
	public void testingContaminantFasta() throws IOException {
		// paramFile with input FASTA file, missedcleavages, etc.
		// this file is pointing to P42681.fasta file
		// change the paths correspondingly
		final File paramFile = new ClassPathResource("blazmass_contaminant.params").getFile();
		// create index
		final String sufix = null;
		final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
		final ProteoformDBIndexInterface proteoformDBIndex = new ProteoformDBIndexInterface(paramFile, sufix, uplr,
				uniprotVersion, maxNumVariationsPerPeptide);

		try {
			System.out.println();
			System.out.println("Looking for the proteins from peptide " + "MAAALLLALAFTLLSGQGACAAAGTIQTSVQEVNSK");
			final Set<IndexedProtein> proteins = proteoformDBIndex.getProteins("MAAALLLALAFTLLSGQGACAAAGTIQTSVQEVNSK");
			for (final IndexedProtein indexedProtein : proteins) {
				System.out.println(indexedProtein.getAccession());
				System.out.println(indexedProtein.getFastaDefLine());
			}

		} catch (final Exception e) {
			e.printStackTrace();
			fail();
		} catch (final DBIndexStoreException e) {
			e.printStackTrace();
			fail();
		}

	}

	@Test
	public void testingProteoformIndex_TitusMar2020_4() throws IOException {
		// paramFile with input FASTA file, missedcleavages, etc.
		// this file is pointing to P42681.fasta file
		// change the paths correspondingly
		final File paramFile = new File(
				"C:\\Users\\salvador\\Desktop\\ProteoformIndex\\UniProt_Human_reviewed_contaminant_04-02-2019_reversed.params");
		// create index
		final String sufix = null;
		final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
		final ProteoformDBIndexInterface proteoformDBIndex = new ProteoformDBIndexInterface(paramFile, sufix, uplr,
				uniprotVersion, maxNumVariationsPerPeptide);
		final double parentMass = 872.508;
		System.out.println(parentMass);
		try {
			final List<IndexedSequence> sequences = proteoformDBIndex.getSequences(parentMass, 0.08725);
			for (final IndexedSequence indexedSequence : sequences) {
				final List<Integer> proteinIds = indexedSequence.getProteinIds();
				for (final Integer proteinID : proteinIds) {
					final IndexedProtein protein = proteoformDBIndex.getIndexedProteinById(proteinID);
					System.out.println(protein.getAccession());
				}
				final List<String> proteinDescArray = indexedSequence.getProteinDescArray();
				for (final String proteinDesc : proteinDescArray) {
					System.out.println(proteinDesc);
				}
				System.out.println(indexedSequence.getSequence());
				final IndexedSequenceWithPTMs pep = (IndexedSequenceWithPTMs) indexedSequence;
				final List<PTM> ptms = pep.getPtms();
				for (final PTM ptm : ptms) {
					System.out.println(ptm.getPosInPeptide());

				}
				System.out.println(pep.getModSequence());

			}
		} catch (final DBIndexStoreException e) {
			e.printStackTrace();
			fail();
		}
	}
}
