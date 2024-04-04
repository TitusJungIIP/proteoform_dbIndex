package edu.scripps.yates.utilities.fasta;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import com.compomics.dbtoolkit.io.EnzymeLoader;
import com.compomics.dbtoolkit.io.implementations.FASTADBLoader;
import com.compomics.dbtoolkit.io.implementations.FASTAHeaderFilter;
import com.compomics.util.protein.Enzyme;
import com.compomics.util.protein.Protein;

import edu.scripps.yates.utilities.venndata.VennData;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class FastaComparator {
	private final static Logger log = Logger.getLogger(FastaComparator.class);
	private final int minPeptideLength = 6;
	private final int maxPeptideLength = Integer.MAX_VALUE;
	private final String DECOY_PREFIX = "Reverse";

	private final Map<String, Set<String>> map;
	private final NumberFormat nf = NumberFormat.getPercentInstance();
	private final Enzyme enzyme;
	private final List<File> files;
	private final File imageFileoutput;
	private VennData vennData;

	// public static void main(String[] args) throws IOException {
	// File fastaFile = new File(
	// "C:\\Users\\Salva\\Desktop\\data\\isotopologues\\databases\\NCBI_RefSeq_Melanogaster__07-01-2014_reversed.fasta");
	//
	// String[] species = { "Drosophila melanogaster", "Drosophila virilis" };
	// FastaComparator fastaComparator = new FastaComparator(fastaFile,
	// species);
	// fastaComparator.analyzeIntersections();
	// }

	public static void main(String[] args) throws IOException {
		final File fastaFile3 = new File(
				"Z:\\share\\Salva\\data\\PINT projects\\DroHybrids\\UniProt_D_simulans_and_melanogaster_11-01-2014.fasta");
		final File fastaFile = new File(
				"Z:\\share\\Salva\\data\\isotopologues\\databases\\D_simulans_canonical_and_isoform.fasta");
		final File fastaFile2 = new File(
				"Z:\\share\\Salva\\data\\isotopologues\\databases\\NCBI_RefSeq_Melanogaster__07-01-2014_reversed.fasta");
		final File venndiagramFile = new File("Z:\\share\\Salva\\data\\isotopologues\\databases\\venn.png");
		final List<File> files = new ArrayList<File>();
		// files.add(fastaFile);
		// files.add(fastaFile2);
		files.add(fastaFile3);
		final String[] species = { // "Drosophila virilis",
				"Drosophila melanogaster", "Drosophila simulans" };

		final FastaComparator fastaComparator = new FastaComparator(files, species, venndiagramFile, "Lys-C");
		fastaComparator.analyzeIntersections();
		fastaComparator.analyzeIntersectionsUsingVennData();
		final URL url = fastaComparator.getVennDiagramURL();
		System.out.println(url);
	}

	public FastaComparator(List<File> files, String[] species, File imageFileOutput, String enzymeName)
			throws IOException {
		this.files = files;
		imageFileoutput = imageFileOutput;

		enzyme = EnzymeLoader.loadEnzyme(enzymeName, null);
		enzyme.setMiscleavages(0);

		map = new THashMap<String, Set<String>>();

		for (final String specie : species) {
			final Set<String> set = new THashSet<String>();
			int numProteinsFromSpecie = 0;
			for (final File file : files) {
				final FASTADBLoader loader = new FASTADBLoader();
				loader.load(file.getAbsolutePath());

				final FASTAHeaderFilter filefilter = new FASTAHeaderFilter(specie);
				Protein entry = loader.nextFilteredProtein(filefilter);
				while (entry != null) {
					if (!entry.getHeader().getRawHeader().contains(DECOY_PREFIX)) {
						numProteinsFromSpecie++;
						final Protein[] peptides = enzyme.cleave(entry, minPeptideLength, maxPeptideLength);
						for (final Protein peptide : peptides) {
							final String sequence = peptide.getSequence().getSequence();
							set.add(sequence);
						}
					}
					entry = loader.nextFilteredProtein(filefilter);
				}
			}
			System.out.println(set.size() + " Different peptides for specie: " + specie + " in " + numProteinsFromSpecie
					+ " proteins");
			map.put(specie, set);
		}
	}

	public void analyzeIntersections() {
		if (map.size() == 2) {
			final Set<String> keySet = map.keySet();
			final List<String> species = new ArrayList<String>();
			species.addAll(keySet);
			Collections.sort(species);
			final Set<String> col1 = map.get(species.get(0));
			final Set<String> col2 = map.get(species.get(1));
			final Set<String> justIn1 = new THashSet<String>();
			final Set<String> justIn2 = new THashSet<String>();
			final Set<String> intersection = new THashSet<String>();
			final Set<String> union = new THashSet<String>();

			for (final String seq1 : col1) {
				if (col2.contains(seq1)) {
					intersection.add(seq1);
				} else {
					justIn1.add(seq1);
				}
				union.add(seq1);
			}
			for (final String seq2 : col2) {
				if (col1.contains(seq2)) {
					intersection.add(seq2);
				} else {
					justIn2.add(seq2);
				}
				union.add(seq2);
			}
			System.out.println("-------------------------------");
			System.out.println("Analysis of " + getSpeciesNames(species));
			System.out.println("Using files: " + getFileNames(files));
			System.out.println("-------------------------------");
			System.out.println("Decoys skipped with prefix: " + DECOY_PREFIX);
			System.out.println("-------------------------------");

			System.out.println("Enzyme used: " + enzyme);
			System.out.println("-------------------------------");
			System.out.println("Peptides exclusive from " + species.get(0) + ": " + justIn1.size() + "("
					+ nf.format(Double.valueOf(justIn1.size()) / union.size()) + ")");
			System.out.println("Peptides exclusive from " + species.get(1) + ": " + justIn2.size() + "("
					+ nf.format(Double.valueOf(justIn2.size()) / union.size()) + ")");
			System.out.println("Peptides present in both species: " + intersection.size() + "("
					+ nf.format(Double.valueOf(intersection.size()) / union.size()) + ")");
			System.out.println("Union: " + union.size());
			System.out.println("-------------------------------");

		}

	}

	public void analyzeIntersectionsUsingVennData() throws IOException {
		if (map.size() <= 3) {

			final Set<String> keySet = map.keySet();
			final List<String> species = new ArrayList<String>();
			species.addAll(keySet);
			Collections.sort(species);
			Set<String> col1 = null;
			Set<String> col2 = null;
			Set<String> col3 = null;
			String name1 = null;
			String name2 = null;
			String name3 = null;

			col1 = map.get(species.get(0));
			name1 = species.get(0);
			if (species.size() > 1) {
				col2 = map.get(species.get(1));
				name2 = species.get(1);
			}
			if (species.size() > 2) {
				col3 = map.get(species.get(2));
				name3 = species.get(2);
			}

			vennData = new VennData(getSpeciesNames(species), name1, col1, name2, col2, name3, col3);
			System.out.println("-------------------------------");
			System.out.println("Analysis of " + getSpeciesNames(species));
			System.out.println("Using files: " + getFileNames(files));
			System.out.println("-------------------------------");
			System.out.println("Decoys skipped with prefix: " + DECOY_PREFIX);
			System.out.println("-------------------------------");

			System.out.println("Enzyme used: " + enzyme);
			System.out.println("-------------------------------");
			System.out.println(vennData.getIntersectionsText(getSpeciesNames(species)));
			System.out.println("-------------------------------");

		}

	}

	public void generateVennDiagram() throws IOException {
		vennData.saveVennDiagramToFile(imageFileoutput);
	}

	public URL getVennDiagramURL() throws IOException {
		return vennData.getImageURL();
	}

	private String getFileNames(List<File> files2) {
		final StringBuilder sb = new StringBuilder();
		for (final File file : files2) {
			if (!"".equals(sb.toString()))
				sb.append(", ");
			sb.append(file.getName());
		}
		return sb.toString();
	}

	private String getSpeciesNames(List<String> species) {
		final StringBuilder sb = new StringBuilder();
		for (final String specie : species) {
			if (!"".equals(sb.toString()))
				sb.append(" vs ");
			sb.append(specie);
		}
		return sb.toString();
	}
}
