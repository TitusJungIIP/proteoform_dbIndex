package edu.scripps.yates.utilities.annotations.uniprot;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.utilities.annotations.uniprot.xml.CommentType;
import edu.scripps.yates.utilities.annotations.uniprot.xml.DbReferenceType;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.annotations.uniprot.xml.EvidencedStringType;
import edu.scripps.yates.utilities.annotations.uniprot.xml.FeatureType;
import edu.scripps.yates.utilities.annotations.uniprot.xml.GeneNameType;
import edu.scripps.yates.utilities.annotations.uniprot.xml.GeneType;
import edu.scripps.yates.utilities.annotations.uniprot.xml.IsoformType;
import edu.scripps.yates.utilities.annotations.uniprot.xml.KeywordType;
import edu.scripps.yates.utilities.annotations.uniprot.xml.PropertyType;
import edu.scripps.yates.utilities.annotations.uniprot.xml.ProteinType.AlternativeName;
import edu.scripps.yates.utilities.annotations.uniprot.xml.ProteinType.SubmittedName;
import edu.scripps.yates.utilities.annotations.uniprot.xml.SubcellularLocationType;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Uniprot;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class UniprotEntryUtil {
	public static void removeNonUsedElements(Uniprot uniprot, boolean removeReferences, boolean removeDBReferences) {
		if (uniprot != null && uniprot.getEntry() != null) {
			for (final Entry entry : uniprot.getEntry()) {
				removeNonUsedElements(entry, removeReferences, removeDBReferences);
			}
		}
	}

	public static void removeNonUsedElements(Entry entry, boolean removeReferences, boolean removeDBReferences) {
		if (entry != null) {
			if (removeReferences && entry.getReference() != null) {
				entry.getReference().clear();
			}
			if (removeDBReferences && entry.getDbReference() != null) {
				entry.getDbReference().clear();
			}
		}
	}

	private static final String ENSEMBL = "Ensembl";
	private static final String GENE_ID = "gene ID";

	/**
	 * Gets a map, of ENSGID identifiers from an uniprot entry. They key of the
	 * entry is the uniprot Accession. An entry can have multiple uniprot ids
	 * (including isoforms), so sometimes, for each isoform, we have a different
	 * ensgid
	 * 
	 * @param entry
	 * @return
	 */
	public static Map<String, Set<String>> getENSGIDs(Entry entry) {
		final Map<String, Set<String>> ret = new THashMap<String, Set<String>>();
		if (entry != null) {
			if (entry.getDbReference() != null) {
				for (final DbReferenceType dbReference : entry.getDbReference()) {
					if (ENSEMBL.equals(dbReference.getType())) {
						if (dbReference.getProperty() != null) {
							String acc = UniprotEntryUtil.getPrimaryAccession(entry);
							if (dbReference.getMolecule() != null) {
								acc = dbReference.getMolecule().getId();
							}
							for (final PropertyType property : dbReference.getProperty()) {
								if (GENE_ID.equals(property.getType())) {
									if (!ret.containsKey(acc)) {
										ret.put(acc, new THashSet<String>());
									}
									ret.get(acc).add(property.getValue());
									if ("1".equals(FastaParser.getIsoformVersion(acc))) {
										final String noIsoformACC = FastaParser.getNoIsoformAccession(acc);
										if (!ret.containsKey(noIsoformACC)) {
											ret.put(noIsoformACC, new THashSet<String>());
										}
										ret.get(noIsoformACC).add(property.getValue());
									}
								}
							}
						}
					}
				}
			}
		}
		return ret;
	}

	public static String getProteinSequence(Entry entry) {
		if (entry != null) {
			if (entry.getSequence() != null) {
				if (entry.getSequence().getValue() != null) {
					final String seq = entry.getSequence().getValue().replace("\n", "");
					return seq;
				}
			}
		}
		return null;
	}

	public static Double getMolecularWeightInDalton(Entry entry) {
		if (entry != null) {
			if (entry.getSequence() != null) {
				if (entry.getSequence().getMass() > 0) {
					return Double.valueOf(entry.getSequence().getMass());
				}
			}
		}
		return null;
	}

	/**
	 * 
	 * @param entry
	 * @param justPrimary
	 * @param secondaryIfPrimaryIsNull
	 * @return a list of pairs: first element is the name of the gene and second is
	 *         the type of the gene
	 */
	public static List<Pair<String, String>> getGeneName(Entry entry, boolean justPrimary,
			boolean secondaryIfPrimaryIsNull) {
		final List<Pair<String, String>> ret = new ArrayList<Pair<String, String>>();
		if (entry != null) {
			final List<GeneType> gene = entry.getGene();
			if (gene != null) {
				for (final GeneType geneType : gene) {
					for (final GeneNameType geneNameType : geneType.getName()) {
						boolean isPrimary = false;
						if ("primary".equals(geneNameType.getType())) {
							isPrimary = true;
						}
						if (!justPrimary || (justPrimary && isPrimary)) {
							ret.add(new Pair<String, String>(geneNameType.getValue(), geneNameType.getType()));
							if (justPrimary) {
								return ret;
							}
						}
					}
				}
				if (justPrimary && secondaryIfPrimaryIsNull && ret.isEmpty()) {
					final List<Pair<String, String>> geneNames2 = getGeneName(entry, false, false);
					if (geneNames2.isEmpty()) {
					} else {
						ret.add(geneNames2.get(0));
					}
				}
			}
		}
		return ret;
	}

	/**
	 * <dbReference type="Proteomes" id="UP000005640">
	 * <property type= "component" value="Chromosome 1"/> </dbReference>
	 * 
	 * @param entry
	 * @return
	 */
	public static String getChromosomeName(Entry entry) {
		if (entry != null) {
			if (entry.getDbReference() != null) {
				for (final DbReferenceType dbReference : entry.getDbReference()) {
					if ("Proteomes".equals(dbReference.getType())) {
						if (dbReference.getProperty() != null) {
							for (final PropertyType property : dbReference.getProperty()) {
								if (property.getValue().contains("Chromosome")) {
									return property.getValue();
								}
							}
						}
					}
				}
			}
		}
		return null;
	}

	public static String getUniprotEvidence(Entry entry) {
		if (entry != null) {
			if (entry.getProteinExistence() != null) {
				return entry.getProteinExistence().getType();
			}
		}
		return null;
	}

	public static boolean isSwissProt(Entry entry) {
		if (entry != null) {
			if (entry.getDataset().equals("Swiss-Prot")) {
				return true;
			}
		}
		return false;
	}

	public static String getProteinDescription(Entry entry) {
		if (entry != null) {
			if (entry.getProtein() != null) {
				if (entry.getProtein().getRecommendedName() != null) {
					if (entry.getProtein().getRecommendedName().getFullName() != null) {
						return entry.getProtein().getRecommendedName().getFullName().getValue();
					}
				} else {
					if (entry.getProtein().getSubmittedName() != null) {
						for (final SubmittedName submittedName : entry.getProtein().getSubmittedName()) {
							final EvidencedStringType fullName = submittedName.getFullName();
							if (fullName != null) {
								final String fullNameValue = fullName.getValue();
								if (fullNameValue != null && !"".equals(fullNameValue)) {
									return fullNameValue;
								}
							}
						}
					}
				}
			}
		}
		return null;
	}

	public static List<String> getAlternativeNames(Entry entry) {
		final List<String> ret = new ArrayList<String>();
		if (entry != null) {
			if (entry.getProtein() != null) {
				if (entry.getProtein().getAlternativeName() != null) {

					for (final AlternativeName alternativeName : entry.getProtein().getAlternativeName()) {
						final EvidencedStringType fullName = alternativeName.getFullName();
						if (fullName != null) {
							final String fullNameValue = fullName.getValue();
							if (fullNameValue != null && !"".equals(fullNameValue)) {
								ret.add(fullNameValue);
							}
						}
						if (alternativeName.getShortName() != null) {
							for (final EvidencedStringType shortNameType : alternativeName.getShortName()) {
								if (shortNameType.getValue() != null) {
									ret.add(shortNameType.getValue());
								}
							}
						}
					}
				}
			}
		}
		return ret;

	}

	public static String getPrimaryAccession(Entry entry) {
		if (entry != null) {
			if (entry.getAccession() != null) {
				if (!entry.getAccession().isEmpty()) {
					return entry.getAccession().get(0);
				}
			}
		}
		return null;
	}

	public static String getTaxonomyName(Entry entry) {
		if (entry != null) {
			if (entry.getOrganism() != null && entry.getOrganism().getName() != null) {
				if (!entry.getOrganism().getName().isEmpty()) {
					return entry.getOrganism().getName().get(0).getValue();
				}
			}
		}
		return null;
	}

	public static String getTaxonomyNCBIID(Entry entry) {
		if (entry != null) {
			if (entry.getOrganism() != null && entry.getOrganism().getDbReference() != null) {
				for (final DbReferenceType dbRef : entry.getOrganism().getDbReference()) {
					if (dbRef.getType() != null && "NCBI Taxonomy".equals(dbRef.getType())) {
						return dbRef.getId();
					}
				}
			}
		}
		return null;
	}

	public static List<String> getNames(Entry entry) {
		if (entry != null) {
			if (entry.getName() != null) {
				return entry.getName();
			}
		}
		return null;
	}

	public static List<String> getCellularLocations(Entry entry) {
		final List<String> ret = new ArrayList<String>();
		if (entry != null) {
			if (entry.getComment() != null) {
				for (final CommentType comment : entry.getComment()) {
					if (comment.getType() != null && "subcellular location".equals(comment.getType())) {
						final List<SubcellularLocationType> subcellularLocations = comment.getSubcellularLocation();
						if (subcellularLocations != null) {
							for (final SubcellularLocationType subcellularLocation : subcellularLocations) {
								final List<EvidencedStringType> locations = subcellularLocation.getLocation();
								if (locations != null) {
									for (final EvidencedStringType evidencedStringType : locations) {
										if (!ret.contains(evidencedStringType.getValue())) {
											ret.add(evidencedStringType.getValue());
										}
									}
								}
							}
						}
					}
				}
			}
		}
		return ret;
	}

	/**
	 * Gets a list of {@link CommentType} of comments that have a certain
	 * commentType. commentType can be one of:<br>
	 * Function<br>
	 * Catalytic activity<br>
	 * Cofactor<br>
	 * Activity regulation<br>
	 * Biophysicochemical properties<br>
	 * Subunit structure<br>
	 * Pathway<br>
	 * Subcellular location<br>
	 * Tissue specificity<br>
	 * Developmental stage<br>
	 * Induction<br>
	 * Domain<br>
	 * Post-translational modification<br>
	 * RNA editing<br>
	 * Mass spectrometry<br>
	 * Polymorphism<br>
	 * Involvement in disease<br>
	 * Disruption phenotype<br>
	 * Allergenic properties<br>
	 * Toxic dose<br>
	 * Biotechnological use<br>
	 * Pharmaceutical use<br>
	 * Miscellaneous<br>
	 * Sequence similarities<br>
	 * Caution<br>
	 * Sequence caution<br>
	 * 
	 * @param entry
	 * @param commentType
	 * @return
	 */
	public static List<CommentType> getCommentsByType(Entry entry, String commentType) {
		final List<CommentType> ret = new ArrayList<CommentType>();
		if (entry != null) {
			if (entry.getComment() != null) {
				for (final CommentType comment : entry.getComment()) {
					if (comment.getType() != null && comment.getType().equalsIgnoreCase(commentType)) {
						ret.add(comment);
					}
				}
			}
		}
		return ret;
	}

	/**
	 * Gets a list of {@link FeatureType} of a type provided. featureType can
	 * be:<br>
	 * Initiator methionine<br>
	 * Signal<br>
	 * Transit peptide<br>
	 * Propeptide<br>
	 * Chain<br>
	 * Peptide<br>
	 * <br>
	 * Topological domain<br>
	 * Transmembrane<br>
	 * Intramembrane<br>
	 * Domain<br>
	 * Repeat<br>
	 * Calcium binding<br>
	 * Zinc finger<br>
	 * DNA binding<br>
	 * Nucleotide binding<br>
	 * Region<br>
	 * Coiled coil<br>
	 * Motif<br>
	 * Compositional bias<br>
	 * <br>
	 * Active site<br>
	 * Metal binding<br>
	 * Binding site<br>
	 * Site<br>
	 * <br>
	 * Non-standard residue<br>
	 * Modified residue<br>
	 * Lipidation<br>
	 * Glycosylation<br>
	 * Disulfide bond<br>
	 * Cross-link<br>
	 * <br>
	 * Alternative sequence<br>
	 * Natural variant<br>
	 * <br>
	 * Mutagenesis<br>
	 * Sequence uncertainty<br>
	 * Sequence conflict<br>
	 * Non-adjacent residues<br>
	 * Non-terminal residue<br>
	 * <br>
	 * Helix<br>
	 * Turn<br>
	 * Beta strand<br>
	 * 
	 * @param entry
	 * @param featureType
	 * @return
	 */
	public static List<FeatureType> getFeaturesByType(Entry entry, String featureType) {
		final List<FeatureType> ret = new ArrayList<FeatureType>();
		if (entry != null) {
			if (entry.getFeature() != null) {
				for (final FeatureType feature : entry.getFeature()) {
					if (feature.getType() != null && feature.getType().equals(featureType)) {
						ret.add(feature);
					}
				}
			}
		}
		return ret;
	}

	/**
	 * Get the available isoforms accessions from an {@link Entry}
	 * 
	 * @param entry
	 * @return
	 */
	public static List<String> getIsoforms(Entry entry) {
		final List<String> ret = new ArrayList<String>();
		if (entry.getComment() != null) {
			for (final CommentType comment : entry.getComment()) {
				if (comment.getIsoform() != null) {
					for (final IsoformType isoform : comment.getIsoform()) {
						final List<String> ids = isoform.getId();
						for (final String acc : ids) {
							final String isoformVersion = FastaParser.getIsoformVersion(acc);
							if (isoformVersion != null) {
								if (Integer.valueOf(isoformVersion) > 1) {
									ret.add(acc);
								}
							}
						}

					}
				}
			}
		}
		return ret;
	}

	public static Integer getSequenceVersion(Entry entry) {
		if (entry != null) {
			if (entry.getSequence() != null) {
				if (entry.getSequence().getVersion() > 0) {

					return entry.getSequence().getVersion();
				}
			}
		}
		return null;
	}

	public static ProteinExistence getProteinExistence(Entry entry) {
		if (entry.getProteinExistence() != null) {
			if (entry.getProteinExistence().getType() != null) {
				for (final ProteinExistence pe : ProteinExistence.values()) {
					if (pe.getDescription().equalsIgnoreCase(entry.getProteinExistence().getType())) {
						return pe;
					}
				}
			}
		}
		return null;
	}

	public static String getFullFastaHeader(Entry entry) {

		final String description = UniprotEntryUtil.getProteinDescription(entry);
		final String primaryAcc = UniprotEntryUtil.getPrimaryAccession(entry);

		final List<String> names = UniprotEntryUtil.getNames(entry);
		String name = null;
		if (names == null || names.isEmpty()) {
			name = description;
		} else {
			name = names.get(0);
		}
		String sp = "sp";
		if (!UniprotEntryUtil.isSwissProt(entry)) {
			sp = "tr";
		}
		String defline = ">" + sp + "|" + primaryAcc + "|" + name + " " + description;
		if (FastaParser.isReverse(primaryAcc)) {
			defline = ">" + primaryAcc + " " + description;
		}
		final List<Pair<String, String>> geneNames = UniprotEntryUtil.getGeneName(entry, true, true);
		if (geneNames != null && !geneNames.isEmpty()) {
			defline += " GN=" + geneNames.get(0).getFirstelement();
		}
		final String taxonomy = UniprotEntryUtil.getTaxonomyName(entry);
		if (taxonomy != null) {
			defline += " OS=" + taxonomy;
		}
		final String taxID = UniprotEntryUtil.getTaxonomyNCBIID(entry);
		if (taxonomy != null) {
			defline += " OX=" + taxID;
		}
		final Integer sequenceVersion = UniprotEntryUtil.getSequenceVersion(entry);
		if (sequenceVersion != null) {
			defline += " SV=" + sequenceVersion;
		}
		final ProteinExistence pe = UniprotEntryUtil.getProteinExistence(entry);
		if (pe != null) {
			defline += " PE=" + pe.getNum();
		}
		return defline;
	}

	/**
	 * Gets keywords into pairs key-value in which key is the id and the value is
	 * the value
	 * 
	 * @param entry
	 * @return
	 */
	public static List<Pair<String, String>> getKeyWords(Entry entry) {
		final List<Pair<String, String>> ret = new ArrayList<Pair<String, String>>();
		if (entry.getKeyword() != null) {
			for (final KeywordType key : entry.getKeyword()) {
				final String id = key.getId();
				final String value = key.getValue();
				ret.add(new Pair<String, String>(id, value));
			}
		}
		return ret;
	}

	/**
	 * Looks into the dbReferences of type GO and looks for the ones that have in
	 * its value of the term name the provided name. If so, a string with the GO ID
	 * and the term name separated by a "=" is returned in a list
	 * 
	 * @param entry
	 * @param name
	 * @return
	 */
	public static List<String> getGeneOntologyContainingName(Entry entry, String name) {
		final List<String> ret = new ArrayList<String>();
		if (entry.getDbReference() != null) {
			for (final DbReferenceType dbReference : entry.getDbReference()) {
				if (dbReference.getType().equals("GO") && dbReference.getProperty() != null) {
					final List<PropertyType> properties = dbReference.getProperty();
					for (final PropertyType property : properties) {
						if ("term".equals(property.getType())) {
							if (property.getValue().toLowerCase().contains(name.toLowerCase())) {
								ret.add(dbReference.getId() + "=" + property.getValue());
							}
						}
					}
				}
			}
		}
		return ret;
	}

	public static Set<String> getDBReferenceIDsByType(Entry entry, String dbReferenceType) {
		final Set<String> ret = new THashSet<String>();
		if (entry.getDbReference() != null) {
			for (final DbReferenceType dbReference : entry.getDbReference()) {
				if (dbReference.getType().equals(dbReferenceType) && dbReference.getId() != null) {

					ret.add(dbReference.getId());

				}
			}
		}
		return ret;
	}

	public static List<String> getGeneOntologyMolecularFunction(Entry entry) {
		final List<String> ret = new ArrayList<String>();
		if (entry.getDbReference() != null) {
			for (final DbReferenceType dbReference : entry.getDbReference()) {
				if (dbReference.getType().equals("GO") && dbReference.getProperty() != null) {
					final List<PropertyType> properties = dbReference.getProperty();
					for (final PropertyType property : properties) {
						if ("term".equals(property.getType())) {
							if (property.getValue().startsWith("F:")) {
								ret.add(dbReference.getId() + "=" + property.getValue());
							}
						}
					}
				}
			}
		}
		return ret;
	}

	public static List<String> getGeneOntologyCelularComponent(Entry entry) {
		final List<String> ret = new ArrayList<String>();
		if (entry.getDbReference() != null) {
			for (final DbReferenceType dbReference : entry.getDbReference()) {
				if (dbReference.getType().equals("GO") && dbReference.getProperty() != null) {
					final List<PropertyType> properties = dbReference.getProperty();
					for (final PropertyType property : properties) {
						if ("term".equals(property.getType())) {
							if (property.getValue().startsWith("C:")) {
								ret.add(dbReference.getId() + "=" + property.getValue());
							}
						}
					}
				}
			}
		}
		return ret;
	}

	public static List<String> getGeneOntologyBiologicalProcess(Entry entry) {
		final List<String> ret = new ArrayList<String>();
		if (entry.getDbReference() != null) {
			for (final DbReferenceType dbReference : entry.getDbReference()) {
				if (dbReference.getType().equals("GO") && dbReference.getProperty() != null) {
					final List<PropertyType> properties = dbReference.getProperty();
					for (final PropertyType property : properties) {
						if ("term".equals(property.getType())) {
							if (property.getValue().startsWith("P:")) {
								ret.add(dbReference.getId() + "=" + property.getValue());
							}
						}
					}
				}
			}
		}
		return ret;
	}

}
