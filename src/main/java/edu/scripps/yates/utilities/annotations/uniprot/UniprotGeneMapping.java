package edu.scripps.yates.utilities.annotations.uniprot;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPFile;
import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.files.FileUtils;
import edu.scripps.yates.utilities.files.ZipManager;
import edu.scripps.yates.utilities.ftp.FTPUtils;
import edu.scripps.yates.utilities.strings.StringUtils;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class UniprotGeneMapping {

	private static final String SEPARATOR = "\t";
	/**
	 * Gene name to Uniprot accession map.
	 */
	private final Map<String, Set<String>> geneNameToAccession = new THashMap<String, Set<String>>();
	private final Map<Gene, Set<String>> geneToAccession = new THashMap<Gene, Set<String>>();
	private final Map<String, Set<Gene>> accessionToGenes = new THashMap<String, Set<Gene>>();
	private boolean loaded = false;
	private final File uniprotPath;
	private final List<String> taxonomies = new ArrayList<String>();
	private static final String uniprotHostName = "ftp.uniprot.org";
	private static Map<String, UniprotGeneMapping> instances = new THashMap<String, UniprotGeneMapping>();
	private static String uniprotFTPRemoteFolder = "/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/";
	private static Logger log = Logger.getLogger(UniprotGeneMapping.class);
	public static final String GENE_SYNONM = "Gene_Synonym";
	public static final String GENE_NAME = "Gene_Name";
	public static final String ENSEMBL = "Ensembl";
	private final boolean mapToENSEMBL;
	private final boolean mapToGENENAME;
	private final boolean mapToGENESYNONIM;
	private final static Set<String> notFoundTaxonomies = new THashSet<String>();

	public UniprotGeneMapping(File uniprotPath, String[] taxonomies, boolean mapToENSEMBL, boolean mapToGENENAME,
			boolean mapToGENESYNONIM) {
		this.uniprotPath = uniprotPath;
		for (final String taxonomy : taxonomies) {
			this.taxonomies.add(taxonomy);
		}
		this.mapToENSEMBL = mapToENSEMBL;
		this.mapToGENENAME = mapToGENENAME;
		this.mapToGENESYNONIM = mapToGENESYNONIM;
	}

	public UniprotGeneMapping(File uniprotPath, List<String> taxonomies, boolean mapToENSEMBL, boolean mapToGENENAME,
			boolean mapToGENESYNONIM) {
		this.uniprotPath = uniprotPath;
		this.taxonomies.addAll(taxonomies);
		this.mapToENSEMBL = mapToENSEMBL;
		this.mapToGENENAME = mapToGENENAME;
		this.mapToGENESYNONIM = mapToGENESYNONIM;
	}

	public UniprotGeneMapping(File uniprotPath, String taxonomy, boolean mapToENSEMBL, boolean mapToGENENAME,
			boolean mapToGENESYNONIM) {
		this.uniprotPath = uniprotPath;
		taxonomies.add(taxonomy);
		this.mapToENSEMBL = mapToENSEMBL;
		this.mapToGENENAME = mapToGENENAME;
		this.mapToGENESYNONIM = mapToGENESYNONIM;
	}

	public static UniprotGeneMapping getInstance(File uniprotPath, String taxonomy) {
		return getInstance(uniprotPath, taxonomy, true, true, true);
	}

	public static UniprotGeneMapping getInstance(File uniprotPath, String[] taxonomies) {
		return getInstance(uniprotPath, taxonomies, true, true, true);
	}

	public static UniprotGeneMapping getInstance(File uniprotPath, String[] taxonomies, boolean mapToENSEMBL,
			boolean mapToGENENAME, boolean mapToGENESYNONIM) {
		final String key = StringUtils.getSortedSeparatedValueStringFromChars(taxonomies, "-") + "ENSEMBL="
				+ mapToENSEMBL + " GENENAME=" + mapToGENENAME + " GENESYNONIM=" + mapToGENESYNONIM;
		if (!instances.containsKey(key)) {
			instances.put(key,
					new UniprotGeneMapping(uniprotPath, taxonomies, mapToENSEMBL, mapToGENENAME, mapToGENESYNONIM));
		}
		return instances.get(key);
	}

	public static UniprotGeneMapping getInstance(File uniprotPath, List<String> taxonomies, boolean mapToENSEMBL,
			boolean mapToGENENAME, boolean mapToGENESYNONIM) {
		final String key = StringUtils.getSortedSeparatedValueStringFromChars(taxonomies, "-") + "ENSEMBL="
				+ mapToENSEMBL + " GENENAME=" + mapToGENENAME + " GENESYNONIM=" + mapToGENESYNONIM;
		if (!instances.containsKey(key)) {
			instances.put(key,
					new UniprotGeneMapping(uniprotPath, taxonomies, mapToENSEMBL, mapToGENENAME, mapToGENESYNONIM));
		}
		return instances.get(key);
	}

	public static UniprotGeneMapping getInstance(File uniprotPath, String taxonomy, boolean mapToENSEMBL,
			boolean mapToGENENAME, boolean mapToGENESYNONIM) {
		final String key = taxonomy + "ENSEMBL=" + mapToENSEMBL + " GENENAME=" + mapToGENENAME + " GENESYNONIM="
				+ mapToGENESYNONIM;
		if (!instances.containsKey(key)) {
			instances.put(key,
					new UniprotGeneMapping(uniprotPath, taxonomy, mapToENSEMBL, mapToGENENAME, mapToGENESYNONIM));
		}
		return instances.get(key);
	}

	public Set<String> mapGeneToUniprotACC(String geneNameOrID) throws IOException {
		if (geneNameOrID == null) {
			return Collections.emptySet();
		}
		if (!loaded) {
			if (!taxonomies.isEmpty()) {
				for (final String taxonomy : taxonomies) {
					importFromFile(taxonomy, getMappingFile(taxonomy));
				}

			}
		}

		if (geneNameToAccession.containsKey(geneNameOrID.trim().toUpperCase())) {
			return geneNameToAccession.get(geneNameOrID.trim().toUpperCase());
		}
		return Collections.emptySet();
	}

	public Set<String> mapGeneToUniprotACC(Gene gene) throws IOException {
		if (gene == null) {
			return Collections.emptySet();
		}
		if (!loaded) {
			if (!taxonomies.isEmpty()) {
				for (final String taxonomy : taxonomies) {
					importFromFile(taxonomy, getMappingFile(taxonomy));
				}
			}
		}

		if (geneToAccession.containsKey(gene)) {
			return geneToAccession.get(gene);
		}
		return Collections.emptySet();
	}

	/**
	 * Gets a Set of all gene names mapped to this uniprot acc
	 * 
	 * @param uniprotACC
	 * @return
	 * @throws IOException
	 */
	public Set<String> mapUniprotACCToGene(String uniprotACC) throws IOException {
		if (uniprotACC == null) {
			return Collections.emptySet();
		}
		if (!loaded) {
			if (!taxonomies.isEmpty()) {
				for (final String taxonomy : taxonomies) {
					importFromFile(taxonomy, getMappingFile(taxonomy));
				}
			}
		}
		if (accessionToGenes.containsKey(uniprotACC)) {
			final Set<Gene> genes = accessionToGenes.get(uniprotACC);
			return genes.stream().map(gene -> gene.getName()).collect(Collectors.toSet());

		}
		return Collections.emptySet();
	}

	/**
	 * Gets a map of all gene names mapped to this uniprot acc mapped by its gene
	 * type
	 * 
	 * @param uniprotACC
	 * @return
	 * @throws IOException
	 */
	public Map<String, Set<String>> mapUniprotACCToGeneByType(String uniprotACC) throws IOException {
		if (uniprotACC == null) {
			return Collections.emptyMap();
		}
		if (!loaded) {
			if (!taxonomies.isEmpty()) {
				for (final String taxonomy : taxonomies) {
					importFromFile(taxonomy, getMappingFile(taxonomy));
				}
			}
		}
		if (accessionToGenes.containsKey(uniprotACC)) {
			final Set<Gene> genes = accessionToGenes.get(uniprotACC);

			final Map<String, Set<String>> map = new THashMap<String, Set<String>>();
			for (final Gene gene : genes) {
				if (!map.containsKey(gene.getType())) {
					map.put(gene.getType(), new THashSet<String>());
				}
				map.get(gene.getType()).add(gene.getName());
			}

			return map;

		}
		return Collections.emptyMap();
	}

	private File getMappingFile(String taxonomy) throws IOException {
		// if taxonomy is null, this means that we cannot map genes
		if (taxonomy == null || notFoundTaxonomies.contains(taxonomy)) {
			return null;
		}
		final File mappingFolder = new File(uniprotPath.getAbsolutePath() + File.separator + "genemap");
		final File[] listFiles = mappingFolder.listFiles();
		if (listFiles != null) {
			for (final File file : listFiles) {
				if (file.length() > 0l) {
					if (FilenameUtils.getExtension(file.getAbsolutePath()).equals("dat")) {
						if (FilenameUtils.getBaseName(file.getAbsolutePath().toLowerCase())
								.contains(taxonomy.toLowerCase())) {
							return file;
						}
					}
				}
			}
		}
		log.info("mapping file for taxonomy '" + taxonomy + "' not found. Going to uniprot at "
				+ uniprotFTPRemoteFolder);
		final FTPClient ftpClient = FTPUtils.loginFTPClient(uniprotHostName, "anonymous", "");
		// if not found, go to remote location
		final List<FTPFile> filesInFolderByExtension = FTPUtils.getFilesInFolderByExtension(ftpClient,
				uniprotFTPRemoteFolder, "gz");
		for (final FTPFile ftpFile : filesInFolderByExtension) {
			final String name = ftpFile.getName();
			if (name.toLowerCase().startsWith(taxonomy.toLowerCase())) {
				log.info("File for taxonomy '" + taxonomy + "' found in server. Now downloading.");
				final String remoteFilePath = uniprotFTPRemoteFolder + ftpFile.getName();
				final File mappinfFileCompressed = new File(mappingFolder.getAbsolutePath() + File.separator + name);
				if (!mappinfFileCompressed.getParentFile().exists()) {
					mappinfFileCompressed.getParentFile().mkdirs();
				}
				mappinfFileCompressed.createNewFile();
				final OutputStream outputStream = new FileOutputStream(mappinfFileCompressed, false);
				FTPUtils.downloadFile(uniprotHostName, "anonymous", "", remoteFilePath, outputStream);
				outputStream.close();
				log.info("File of " + FileUtils.getDescriptiveSizeFromBytes(mappinfFileCompressed.length())
						+ " downloaded to  " + mappinfFileCompressed.getAbsolutePath());
				log.info("Decompressing file...");
				final File file = ZipManager.decompressGZipFile(mappinfFileCompressed);
				log.info("File decompressed to: " + file.getAbsolutePath() + " - "
						+ FileUtils.getDescriptiveSizeFromBytes(file.length()));
				if (file != null && file.exists()) {
					log.info("Deleting file " + mappinfFileCompressed.getAbsolutePath());
					final boolean delete = mappinfFileCompressed.delete();
					if (delete) {
						log.info("File deleted");
					} else {
						log.warn("File cannot be deleted!");
					}
					return file;
				}

			}
		}
		notFoundTaxonomies.add(taxonomy);
		return null;
	}

	private void importFromFile(String taxonomy, File file) throws IOException {
		if (file == null) {
			throw new IllegalArgumentException("Mapping file for taxonomy " + taxonomy + " not available");
		}
		log.info("Loading mapping file: " + file.getAbsolutePath() + " ("
				+ FileUtils.getDescriptiveSizeFromBytes(file.length()) + ")");
		// read the species list
		final FileReader r = new FileReader(file);
		try {
			final BufferedReader br = new BufferedReader(r);
			try {

				String line;

				while ((line = br.readLine()) != null) {

					final String[] splittedLine = line.split(SEPARATOR);

					if (splittedLine.length == 3 && !splittedLine[0].equals("") && !splittedLine[1].equals("")) {

						final String accession = splittedLine[0];
						final String type = splittedLine[1];
						String geneName = null;
						if (mapToGENENAME && type.equals(GENE_NAME)) {
							geneName = splittedLine[2].trim().toUpperCase();
							if (geneName != null) {
								if (!geneNameToAccession.containsKey(geneName)) {
									geneNameToAccession.put(geneName, new THashSet<String>());
								}
								geneNameToAccession.get(geneName).add(accession);
								final Gene gene = new Gene(geneName, type);
								if (!geneToAccession.containsKey(gene)) {
									geneToAccession.put(gene, new THashSet<String>());
								}
								geneToAccession.get(gene).add(accession);
							}
						}

						if (mapToGENESYNONIM && type.equals(GENE_SYNONM)) {
							geneName = splittedLine[2].trim().toUpperCase();
							if (geneName != null) {
								if (!geneNameToAccession.containsKey(geneName)) {
									geneNameToAccession.put(geneName, new THashSet<String>());
								}
								geneNameToAccession.get(geneName).add(accession);
								final Gene gene = new Gene(geneName, type);
								if (!geneToAccession.containsKey(gene)) {
									geneToAccession.put(gene, new THashSet<String>());
								}
								geneToAccession.get(gene).add(accession);
							}
						}

						if (mapToENSEMBL && type.equals(ENSEMBL)) {
							geneName = splittedLine[2].trim().toUpperCase();
							if (geneName != null) {
								if (!geneNameToAccession.containsKey(geneName)) {
									geneNameToAccession.put(geneName, new THashSet<String>());
								}
								geneNameToAccession.get(geneName).add(accession);
								final Gene gene = new Gene(geneName, type);
								if (!geneToAccession.containsKey(gene)) {
									geneToAccession.put(gene, new THashSet<String>());
								}
								geneToAccession.get(gene).add(accession);
							}
						}
					}

				}
			} finally {
				br.close();
				loaded = true;
				log.info(geneNameToAccession.size() + " genes mapped to uniprot accessions");
				getReverseMapping();
				log.info(accessionToGenes.size() + " uniprot accessions mapped to genes");
			}

		} finally {
			r.close();
		}
	}

	private void getReverseMapping() {
		for (final Gene gene : geneToAccession.keySet()) {
			final Set<String> accs = geneToAccession.get(gene);
			for (final String acc : accs) {
				if (!accessionToGenes.containsKey(acc)) {
					accessionToGenes.put(acc, new THashSet<Gene>());
				}
				accessionToGenes.get(acc).add(gene);
			}
		}

	}
}
