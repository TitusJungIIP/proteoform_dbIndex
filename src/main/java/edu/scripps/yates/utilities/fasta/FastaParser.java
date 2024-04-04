package edu.scripps.yates.utilities.fasta;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

import org.apache.commons.lang.math.NumberUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.masses.AssignMass;
import edu.scripps.yates.utilities.proteomicsmodel.Accession;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AccessionType;
import edu.scripps.yates.utilities.proteomicsmodel.factories.AccessionEx;
import edu.scripps.yates.utilities.strings.StringUtils;
import edu.scripps.yates.utilities.taxonomy.UniprotOrganism;
import edu.scripps.yates.utilities.taxonomy.UniprotSpeciesCodeMap;
import edu.scripps.yates.utilities.util.StringPosition;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntDoubleHashMap;
import gnu.trove.set.hash.THashSet;

public class FastaParser {
	private static Set<String> notRecognizedFastas = new THashSet<String>();

	public static enum UNIPROT_FASTA_KEYWORD {
		OS, GN, PE, SV, OX
	};

	public static final String UNKNOWN = "UNKNOWN";
	public static final String CONTAMINANT = "contaminant";
	public static final String UNIPROT = "UNIPROT";
	public static final String NCBI = "NCBI";

	public static final String REVERSE = "Reverse";

	// taken from http://www.uniprot.org/help/accession_numbers
	public static final Pattern UNIPROT_ACC = Pattern
			.compile("[\\W^]([OPQ][0-9][A-Z0-9]{3}[0-9](\\-[0-9]{1,2})?|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})");

	public static final Pattern UNIPROT_ACC_TYPE1 = Pattern
			.compile("([OPQ][0-9][A-Z0-9][A-Z0-9][A-Z0-9][0-9](\\-[0-9]{1,2})?)");
	public static final Pattern UNIPROT_ACC_TYPE2 = Pattern.compile("([A-NR-Z][0-9][A-Z][A-Z0-9][A-Z0-9][0-9])");
	public static final Pattern UNIPROT_ACC_TYPE3 = Pattern
			.compile("([A-NR-Z][0-9][A-Z][A-Z0-9][A-Z0-9][0-9][A-Z][A-Z0-9][A-Z0-9][0-9])");

	public static final Pattern ONLY_UNIPROT_ACC = Pattern
			.compile("([OPQ][0-9][A-Z0-9]{3}[0-9](\\-[0-9]{1,2})?|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})");
	public static final Pattern NCBI_ACC = Pattern.compile(".*(gi\\|[^|]+)\\|.*");

	private static Pattern IPI_ACC = Pattern.compile(".*(IPI[0-9]{8}\\.*[0-9]*).*");

	public static final Pattern OS_FASTA_HEADER = Pattern.compile(".*OS=(.*)(GN)?PE=.*");
	public static final Pattern BRACKETS_FASTA_HEADER = Pattern.compile(".*\\[(.*)\\].*");
	public static final Pattern GENE_FASTA_HEADER = Pattern.compile(".*GN=(\\S+)(PE)?");
	public static final Pattern GENE_SYMBOL_FASTA_HEADER = Pattern.compile(".*Gene_Symbol=(\\S+).*");
	public static final Pattern PE_FASTA_HEADER = Pattern.compile(".*PE=(\\d+).*");
	public static final Pattern SV_FASTA_HEADER = Pattern.compile(".*SV=(\\d+).*");
	public static final Pattern TAXID_FASTA_HEADER = Pattern.compile(".*Tax_Id=(\\d+).*");
	/**
	 * regular expressions to capture data that will be removed since these
	 * annotations can be located in the middle or at the beginning of the fasta
	 * header
	 */
	public static final String[] FASTA_ANNOTATIONS = { "Tax_Id=\\d+", "Gene_Symbol=\\S+" };
	public static Pattern[] FASTA_ANNOTATION_PATTERNS;
	public static final String CONTAMINANT_PREFIX = "contaminant";

	private static Logger log = Logger.getLogger(FastaParser.class);
//	private static final Map<String, String> uniprotAccCache = new THashMap<String, String>();
	private static final Pattern UNIPROT_SP_ACC_TMP = Pattern.compile(".*s.\\|(\\S+)\\|\\S*\\s*.*");
	private static final Pattern UNIPROT_TR_ACC_TMP = Pattern.compile(".*tr\\|(\\S+)\\|\\S*\\s*.*");
	private static final Pattern UNIPROT_SP_ACC_TMP2 = Pattern.compile(".*sp\\|\\S+\\|\\S*\\s(.*)");
	private static final Pattern UNIPROT_TR_ACC_TMP2 = Pattern.compile(".*tr\\|\\S+\\|\\S*\\s(.*)");
	private static final Pattern dashAndNumber = Pattern.compile("^(-[0-9]+).*");
	private static final Pattern startingByWord = Pattern.compile("^\\w+");
	private static final Pattern endingByWord = Pattern.compile("\\w+$");
	private static final Pattern untilSpace = Pattern.compile("^(\\S+)");

	public static final String conflict = "_conflict_";
	public static final String mutated = "_mutated_";
	public static final String variant = "_variant_";
	static {
		FASTA_ANNOTATION_PATTERNS = new Pattern[FASTA_ANNOTATIONS.length];
		int index = 0;
		for (final String fasta_annotation : FASTA_ANNOTATIONS) {
			FASTA_ANNOTATION_PATTERNS[index] = Pattern.compile(fasta_annotation);
			index++;
		}
	}

	/**
	 * From a header like >sp|Q96PG8|BBC3B_HUMAN Bcl-2-binding component 3 OS=Homo
	 * sapiens OX=9606 GN=BBC3 PE=2 SV=2 it returns BBC3B_HUMAN, which is the string
	 * after the second '|' and before the first space<br>
	 * From a header like >contaminant_BASI_MOUSE owl|P18572| BASIGIN PRECURSOR
	 * (BASIC IMMUNOGLOBULIN SUPERFAMILY) (MEMBRANE... take everything after first
	 * space if exists and if not, take everything
	 * 
	 * @param fastaHeader
	 * @return
	 */
	public static String getUniProtProteinName(String fastaHeader) {
		final TIntArrayList allPositionsOfPipe = StringUtils.allPositionsOf(fastaHeader, '|');
		final TIntArrayList allPositionsOfSpace = StringUtils.allPositionsOf(fastaHeader, ' ');

		if (!allPositionsOfPipe.isEmpty()) {
			if (!allPositionsOfSpace.isEmpty()) {
				final int secondPipePosition = allPositionsOfPipe.get(1);
				final int firstPositionOfSpace = allPositionsOfSpace.get(0);
				try {
					final String name = fastaHeader.substring(secondPipePosition, firstPositionOfSpace);
					return name;
				} catch (final Exception e) {
					// it is the case described as second in the comments of the function
					// so take everthing after the space
					if (fastaHeader.indexOf(" ") >= 0) {
						return fastaHeader.substring(fastaHeader.indexOf(" "));
					} else {
						return fastaHeader;
					}
				}
			}
		}
		return null;
	}

	/**
	 * Gets whatever is in between of sp| | or tr| | not checking whether that is a
	 * uniprot format or not.<br>
	 * This has been used because IP2 parses these accessions getting what is
	 * between the two '|', and using this, it will be the same.
	 * 
	 * @param id
	 * @return
	 */
	public static String getSPorTRAccession(String id) {
		if (id != null && !"".equals(id)) {
			if (id.length() > 600) {
//				uniprotAccCache.put(id, null);
				return null;
			}

			final Matcher matcherSP = UNIPROT_SP_ACC_TMP.matcher(id);
			final Matcher matcherTR = UNIPROT_TR_ACC_TMP.matcher(id);
			if (id.contains("|")) {
				if (matcherSP.find()) {
					final String group = matcherSP.group(1);
					if (!"".equals(group))
						id = group.trim();
				} else if (matcherTR.find()) {
					final String group = matcherTR.group(1);
					if (!"".equals(group))
						id = group.trim();
				}
				return id;
			}

		}
		return null;
	}

	/**
	 * This function parses the uniprot accession, getting the resulting string
	 * after applying this regular expression:
	 * ".*([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}).*"
	 * If the regular expresion doens't fit, then the original id is returned.
	 *
	 * @param id
	 * @return nulls if that is not a uniprot entry
	 */
	public static String getUniProtACC(String id) {
//		if (uniprotAccCache.containsKey(id)) {
//			return uniprotAccCache.get(id);
//		}

		if (id != null && !"".equals(id)) {
			if (id.length() > 600) {
//				uniprotAccCache.put(id, null);
				return null;
			}

			final Matcher matcherSP = UNIPROT_SP_ACC_TMP.matcher(id);
			final Matcher matcherTR = UNIPROT_TR_ACC_TMP.matcher(id);
			if (id.contains("|")) {
				if (matcherSP.find()) {
					final String group = matcherSP.group(1);
					if (!"".equals(group))
						id = group.trim();
				} else if (matcherTR.find()) {
					final String group = matcherTR.group(1);
					if (!"".equals(group))
						id = group.trim();
				}
			}
			String trim = null;
			final Matcher matcher2 = UNIPROT_ACC_TYPE3.matcher(id);
			if (matcher2.find()) {
				trim = matcher2.group(1).trim();
			} else {
				final Matcher matcher3 = UNIPROT_ACC_TYPE2.matcher(id);
				if (matcher3.find()) {
					trim = matcher3.group(1).trim();
				} else {
					final Matcher matcher4 = UNIPROT_ACC_TYPE1.matcher(id);
					if (matcher4.find()) {
						trim = matcher4.group(1).trim();
					}
				}
			}

			if (trim != null) {

				// get the string after the trim in the id
				String tmp = id.substring(id.indexOf(trim) + trim.length());
				if (tmp.length() > 0) {
					final Matcher matcher6 = dashAndNumber.matcher(tmp);
					if (matcher6.find()) {
						final String string = trim + matcher6.group(1);
//						uniprotAccCache.put(id, string);
						return string;
					}
					// check that there is no more words before the id
					// ??
					// if (startingByWord.matcher(tmp).find()) {
					// return id;
					// }

				}
				tmp = id.substring(0, id.indexOf(trim));
				if (tmp.length() > 0) {

					if (endingByWord.matcher(tmp).find()) {
//						uniprotAccCache.put(id, null);
						return null;
					}
				}
//				uniprotAccCache.put(id, trim);
				return trim;
			}
		}
//		uniprotAccCache.put(id, null);
		return null;
	}

	public static boolean isUniProtACC(String id) {
		return getUniProtACC(id) != null;

		// if (id != null && !"".equals(id)) {
		// final Matcher matcher = ONLY_UNIPROT_ACC.matcher(id);
		// if (matcher.find()) {
		// return true;
		// }
		// }
		// return false;
	}

	public static String getNCBIACC(String id) {
		if (id != null && !"".equals(id)) {
			if (id.length() > 600) {
				return null;
			}
			if (!id.toLowerCase().contains("gi")) {
				return null;
			}
			final Matcher matcher = NCBI_ACC.matcher(id);
			if (matcher.find()) {
				return matcher.group(1).trim();
			}
		}
		return null;
	}

	public static String getIPIACC(String id) {
		if (id != null && !"".equals(id)) {
			if (id.length() > 600) {
				return null;
			}
			final Matcher matcher = IPI_ACC.matcher(id);
			if (matcher.find()) {
				return matcher.group(1).trim();
			}
		}
		return null;
	}

	/**
	 * Gets a pair where the first element is the parsed protein accession and the
	 * second is the accession type
	 *
	 * @param id
	 * @return never null.
	 */
	public static Accession getACC(String id) {

		if (id == null) {
			return new AccessionEx("", AccessionType.UNKNOWN);
		}

		// if starts by >, remove it
		if (id.startsWith(">"))
			id = id.substring(1);
		// if starts we find Reverse in the string before the first space (or everything
		// if there is no space), take that
		// word (everything until the space) as the id
		final String reverse = REVERSE;
		final boolean reverseFitsInString = id.length() >= reverse.length();
		final boolean idContainsSpace = id.contains(" ");
		if ((reverseFitsInString && idContainsSpace && id.substring(0, id.indexOf(" ")).contains(reverse))
				|| (reverseFitsInString && !idContainsSpace && id.contains(reverse))) {
			final Matcher matcher = untilSpace.matcher(id);
			if (matcher.find()) {
				return new AccessionEx(matcher.group(0).trim(), AccessionType.UNKNOWN);
			} else {
				return new AccessionEx(id.trim(), AccessionType.UNKNOWN);
			}
		}
		if (id.length() >= CONTAMINANT.length()) {
			if (id.toLowerCase().startsWith(CONTAMINANT.toLowerCase())
					|| id.toLowerCase().startsWith("tr|" + CONTAMINANT.toLowerCase())
					|| id.toLowerCase().startsWith("sp|" + CONTAMINANT.toLowerCase())) {

				final Matcher matcher = untilSpace.matcher(id);
				if (matcher.find()) {
					return new AccessionEx(matcher.group(0).trim(), AccessionType.UNKNOWN);
				} else {
					return new AccessionEx(id.trim(), AccessionType.UNKNOWN);
				}
			}
		}
		final String uniProtACC = getUniProtACC(id);
		if (uniProtACC != null) {
			return new AccessionEx(uniProtACC, AccessionType.UNIPROT);
		}
		final String ncbiacc = getNCBIACC(id);
		if (ncbiacc != null) {
			return new AccessionEx(ncbiacc, AccessionType.NCBI);
		}
		final String ipiacc = getIPIACC(id);
		if (ipiacc != null) {
			return new AccessionEx(ipiacc, AccessionType.IPI);
		}

		// if contains an space, take the string before the space
		if (idContainsSpace) {
			return new AccessionEx(id.substring(0, id.indexOf(" ")), AccessionType.UNKNOWN);
		}
		return new AccessionEx(id.trim(), AccessionType.UNKNOWN);
	}

	/**
	 * Gets the OX=1960 name from the fasta header.
	 *
	 * @param fastaheader
	 * @return
	 */
	public static String getOrganismNCBIIDFromFastaHeader(String fastaheader) {
		if (fastaheader != null && notRecognizedFastas.contains(fastaheader)) {
			return null;
		}
		if (fastaheader != null && fastaheader.length() > 600) {
			return null;
		}
		if (fastaheader != null) {
			// new way 4 Feb 2015
			final String[] split = fastaheader.split(" ");
			for (int i = 0; i < split.length; i++) {
				final String string = split[i];
				if (startsWithUniprotKeyword(string)) {
					if (string.startsWith(UNIPROT_FASTA_KEYWORD.OX.name())) {
						final String tmp = getStringAfterEqual(string);
						String ret = "";
						// if (tmp != null) {
						// ret = tmp;
						// }
						if (tmp == null) {
							continue;
						}
						ret = tmp;
						for (int j = i + 1; j < split.length; j++) {
							// get everything until another keyword
							if (!startsWithUniprotKeyword(split[j])) {
								ret += " " + split[j];
							} else {
								break;
							}
						}
						return ret.trim();
					}
				}
			}
		}
		return null;
	}

	/**
	 * Gets the OS=taxonomy name from the fasta header. Note that it is necessary to
	 * be before the GN= annotation. It also takes the taxonomy if it is between
	 * brackets if it is not found before.<br>
	 * if it is not found yet, it uses the fullaccession provided to get the Uniprot
	 * code if available and use the {@link UniprotSpeciesCodeMap} class to get the
	 * species name.
	 *
	 * @param fastaheader
	 * @return
	 */
	public static String getOrganismNameFromFastaHeader(String fastaheader, String fullaccession) {
		if (fastaheader != null && notRecognizedFastas.contains(fastaheader) && fullaccession != null
				&& !notRecognizedFastas.contains(fullaccession)) {
			return null;
		}
		if (fastaheader != null && fastaheader.length() > 600 && fullaccession != null
				&& fullaccession.length() > 600) {
			return null;
		}
		if (fastaheader != null) {
			// new way 4 Feb 2015
			final String[] split = fastaheader.split(" ");
			for (int i = 0; i < split.length; i++) {
				final String string = split[i];
				if (startsWithUniprotKeyword(string)) {
					if (string.startsWith(UNIPROT_FASTA_KEYWORD.OS.name())) {
						final String tmp = getStringAfterEqual(string);
						String ret = "";
						// if (tmp != null) {
						// ret = tmp;
						// }
						if (tmp == null) {
							continue;
						}
						ret = tmp;
						for (int j = i + 1; j < split.length; j++) {
							// get everything until another keyword
							if (!startsWithUniprotKeyword(split[j])) {
								ret += " " + split[j];
							} else {
								break;
							}
						}
						return ret.trim();
					}
				}
			}
		}
		if (fastaheader != null && !"".equals(fastaheader)) {
			final Matcher matcher = OS_FASTA_HEADER.matcher(fastaheader);
			if (matcher.find()) {
				final String osString = matcher.group(1).trim();
				final String taxonomyFromString = getTaxonomyFromString(osString);
				if (taxonomyFromString != null)
					return taxonomyFromString.trim();
			}

			final Matcher matcher3 = TAXID_FASTA_HEADER.matcher(fastaheader);
			if (matcher3.find()) {
				final String taxID = matcher3.group(1).trim();
				final String taxonomyFromString = getTaxonomyFromString(taxID);
				if (taxonomyFromString != null)
					return taxonomyFromString.trim();
			}
			final Matcher matcher2 = BRACKETS_FASTA_HEADER.matcher(fastaheader);
			if (matcher2.find()) {
				final String insideBrackets = matcher2.group(1).trim();
				final String taxonomyFromString = getTaxonomyFromString(insideBrackets);
				if (taxonomyFromString != null)
					return taxonomyFromString.trim();
			}
		}
		if (fullaccession != null && !"".equals(fullaccession)) {
			final Matcher matcher4 = BRACKETS_FASTA_HEADER.matcher(fullaccession);
			if (matcher4.find()) {
				final String insideBrackets = matcher4.group(1).trim();
				final String taxonomyFromString = getTaxonomyFromString(insideBrackets);
				if (taxonomyFromString != null)
					return taxonomyFromString.trim();
			}
			// try to get from the code of the uniprot like _DROVI

			final String code = getUniprotTaxonomyCode(fullaccession);
			if (code != null) {
				final UniprotOrganism organism = UniprotSpeciesCodeMap.getInstance().get(code);
				if (organism != null) {
					return organism.getScientificName().trim();
				}
			}
		}
		if (fastaheader != null && !"".equals(fastaheader)) {

			final String code = getUniprotTaxonomyCode(fastaheader);
			if (code != null) {
				final UniprotOrganism organism = UniprotSpeciesCodeMap.getInstance().get(code);
				if (organism != null) {
					return organism.getScientificName().trim();
				}
			}
		}

		if (!notRecognizedFastas.contains(fastaheader)) {
			notRecognizedFastas.add(fastaheader);
		}
		if (!notRecognizedFastas.contains(fullaccession)) {
			notRecognizedFastas.add(fullaccession);
		}

		return null;
	}

	/**
	 * Gets the string just after an "="
	 *
	 * @param string
	 * @return
	 */
	private static String getStringAfterEqual(String string) {
		if (string.contains("=")) {
			return string.substring(string.indexOf("=") + 1).trim();
		}
		return null;
	}

	/**
	 * True if the string starts with any UNIPROT KEYWORD
	 *
	 * @param string
	 * @return
	 */
	private static boolean startsWithUniprotKeyword(String string) {
		for (final UNIPROT_FASTA_KEYWORD keyw : UNIPROT_FASTA_KEYWORD.values()) {
			if (string.startsWith(keyw.name())) {
				return true;
			}
		}
		return false;
	}

	private static String getTaxonomyFromString(String taxonomyString) {
		final UniprotOrganism uniprotOrganism = UniprotSpeciesCodeMap.getInstance().get(taxonomyString);
		if (uniprotOrganism != null) {
			return uniprotOrganism.getScientificName();
		}
		return null;
	}

	/**
	 * If the uniprot accession is like >Reverse_gi|75022131|sp|Q9GRW7.1|NONA_DROVI,
	 * it gets "_DROVI"
	 *
	 * @param fastaheader
	 * @return
	 */
	private static String getUniprotTaxonomyCode(String string) {
		if (string.lastIndexOf("_") > -1) {
			String trim = string.substring(string.lastIndexOf("_") + 1).trim();
			if (trim.contains(" "))
				trim = trim.substring(0, trim.indexOf(" "));
			return trim;
		}
		return null;
	}

	/**
	 * Gets the GN=gene name from the fasta header. Note that it is necessary to be
	 * before the annotation PE
	 *
	 * @param fastaheader
	 * @return
	 */
	public static String getGeneFromFastaHeader(String fastaheader) {
		if (fastaheader != null && !"".equals(fastaheader)) {
			final Matcher matcher = GENE_FASTA_HEADER.matcher(fastaheader);
			if (matcher.find()) {
				return getFirstWordIfMany(matcher.group(1).trim());
			}
			final Matcher matcher2 = GENE_SYMBOL_FASTA_HEADER.matcher(fastaheader);
			if (matcher2.find()) {
				return getFirstWordIfMany(matcher2.group(1).trim());
			}
		}
		return null;
	}

	/**
	 * If the string contains separators like "," or ";", take just the first word
	 *
	 * @param trim
	 * @return
	 */
	private static String getFirstWordIfMany(String string) {
		if (string.contains(","))
			return string.split(",")[0];
		if (string.contains(";"))
			return string.split(";")[0];
		return string;
	}

	/**
	 * Gets the protein existence (PE=xx) annotation from the fasta header
	 *
	 * @param fastaheader
	 * @return
	 */
	public static String getProteinExistenceFromFastaHeader(String fastaheader) {
		if (fastaheader != null && !"".equals(fastaheader)) {
			final Matcher matcher = PE_FASTA_HEADER.matcher(fastaheader);
			if (matcher.find()) {
				return matcher.group(1).trim();
			}
		}
		return null;

	}

	/**
	 * Gets the sequence version (SV=xx) annotation from the fasta header
	 *
	 * @param fastaheader
	 * @return
	 */
	public static String getSequenceVersionFromFastaHeader(String fastaheader) {
		if (fastaheader != null && !"".equals(fastaheader)) {
			final Matcher matcher = SV_FASTA_HEADER.matcher(fastaheader);
			if (matcher.find()) {
				return matcher.group(1).trim();
			}
		}
		return null;

	}

	/**
	 * This method parse the description string removing annotations from the header
	 * like:<br>
	 * "Tax_id=1234" or "OS=Drosophila virilis" or "GN=Dvir\GJ22546"...<br>
	 * It also take out the sp|ASDFF|ASDF if exists.
	 *
	 * @param description
	 * @return
	 */
	public static String getDescription(String description) {
		if (description == null)
			return null;
		if (description.startsWith(">"))
			description = description.substring(1);
		// log.debug("Trying to get description from '" + description + "'");
		final String ncbiacc = getNCBIACC(description);
		if (ncbiacc != null && !description.equals(ncbiacc)) {
			if (description.indexOf(" ") >= 0)
				description = description.substring(description.indexOf(" ")).trim();
		}
		// take out the sp|ASDF|ASDF or tr|ASDF|ASDF

		final Matcher matcherSP = UNIPROT_SP_ACC_TMP2.matcher(description);
		final Matcher matcherTR = UNIPROT_TR_ACC_TMP2.matcher(description);
		if (matcherSP.find()) {
			final String group = matcherSP.group(1);
			if (!"".equals(group))
				description = group.trim();
		} else if (matcherTR.find()) {
			final String group = matcherTR.group(1);
			if (!"".equals(group))
				description = group.trim();
		}

		for (final Pattern regex : FASTA_ANNOTATION_PATTERNS) {
			final Matcher regexMatcher = regex.matcher(description);
			description = regexMatcher.replaceFirst("").trim();
		}

		// remove all data from the first annotation
		int firstAnnotationIndex = Integer.MAX_VALUE;

		// Annotations that are always after the description
		final String[] annotations = new String[UNIPROT_FASTA_KEYWORD.values().length];
		int i = 0;
		for (final UNIPROT_FASTA_KEYWORD keyword : UNIPROT_FASTA_KEYWORD.values()) {
			annotations[i++] = keyword.name() + "=";
		}

		for (final String annotation : annotations) {
			final int os = description.indexOf(annotation);
			if (os != -1 && os < firstAnnotationIndex) {
				firstAnnotationIndex = os;
			}
		}

		if (firstAnnotationIndex != Integer.MAX_VALUE && firstAnnotationIndex > -1) {
			final String ret = description.substring(0, firstAnnotationIndex).trim();
			if (!"".equals(ret))
				return ret;
		}

		// if there is something between brackets and it matches with a organism
		// name, remove it
		try {
			// log.info(description);

			final Matcher matcher2 = BRACKETS_FASTA_HEADER.matcher(description);
			if (matcher2.find()) {

				final String match = matcher2.group(1);
				final String taxonomy = getTaxonomyFromString(match);
				if (taxonomy != null) {
					final String quote = Pattern.quote("[" + match + "]");
					description = description.replaceFirst(quote, "").trim();
				}
			}
		} catch (final PatternSyntaxException e) {
			e.printStackTrace();
		}
		return description.trim();
	}

	private static boolean tolerantToFakePeptides = false;

	/**
	 * Sets the value of this static variable that if it is true, not exceptions
	 * will be thrown by this class regarding finding weird strings as peptides.
	 * 
	 * @param b
	 */
	public static void setTolerantToFakePeptides(boolean b) {
		tolerantToFakePeptides = b;
	}

	/**
	 * This function allow to get the peptide sequence as <br>
	 * <ul>
	 * <li>K.VDLSFSPSQSLPASHAHLR.V -> VDLSFSPSQSLPASHAHLR</li>
	 * <li>R.LLLQQVSLPELPGEYSMK.V + Oxidation (M) -> LLLQQVSLPELPGEYSMK</li>
	 * <li>(-)TVAAPSVFIFPPSDEQLK(S) -> TVAAPSVFIFPPSDEQLK</li>
	 * <li>K.EKS[167.00]KESAIASTEVK.L -> EKSKESAIASTEVK</li>
	 * </ul>
	 * getting just the sequence without modifications and between the pre and post
	 * AA if available
	 *
	 * @param seq
	 * @return
	 */
	public static String cleanSequence(String seq) {
		if (seq == null)
			return null;
		final String seqTmp = seq.trim();

		if (somethingExtrangeInSequence(seqTmp)) {

			// parenthesis or brackets
			final List<String> outside = getOutside(seqTmp);
			if (!outside.isEmpty()) {
				final String tmp = appendList(outside);
				final String removeBeforeAfterAAs = removeBeforeAfterAAs(tmp);
				if (!removeBeforeAfterAAs.equals(seq)) {
					return cleanSequence(removeBeforeAfterAAs);
				}
			}

		}
		final String errorMessage = "Peptide sequence '" + seq
				+ "' is not supported. Either having not recognizable characteres or in lower case? Has it a non standard PTM enconded on it? PTMs can be encoded as in PEPTID[+45.92]E";
		AssignMass.getInstance(true);
		for (int index = 0; index < seq.length(); index++) {
			final char aa = seq.charAt(index);
			if (!tolerantToFakePeptides && !AssignMass.containsMass(aa)) {
				throw new IllegalArgumentException("'" + aa + "' not recognized. " + errorMessage);
			}
		}
		if (!seq.toUpperCase().equals(seq)) {
			// it has something in lower case

			throw new IllegalArgumentException(errorMessage);

		}
		return seq.toUpperCase();

	}

	/**
	 * This function works the same as cleanSequence, but additionally, if finds
	 * some annotation of a sequence variance like [ADF->ITV], applies the
	 * substitution and returns the resulting sequence with the variance.
	 *
	 * @param seq
	 * @return
	 */
	public static String cleanSequenceAndApplySequenceVariances(String seq) {
		if (seq == null)
			return null;
		final String seqTmp = seq.trim();

		if (somethingExtrangeInSequence(seqTmp)) {

			// parenthesis or brackets
			final List<String> outside = getOutsideKeepingVariants(seqTmp);
			if (!outside.isEmpty()) {

				final String tmp = appendListApplyingVariants(outside);
				final String removeBeforeAfterAAs = removeBeforeAfterAAs(tmp);
				if (!removeBeforeAfterAAs.equals(seq)) {
					return cleanSequence(removeBeforeAfterAAs);
				}
			}

		}
		final String errorMessage = "Peptide sequence '" + seq
				+ "' is not supported. Either having not recognizable characteres or in lower case? Has it a non standard PTM enconded on it? PTMs can be encoded as in PEPTID[+45.92]E";
		AssignMass.getInstance(true);
		for (int index = 0; index < seq.length(); index++) {
			final char aa = seq.charAt(index);
			if (!tolerantToFakePeptides && !AssignMass.containsMass(aa)) {
				throw new IllegalArgumentException("'" + aa + "' not recognized. " + errorMessage);
			}
		}
		if (!seq.toUpperCase().equals(seq)) {
			// it has something in lower case

			throw new IllegalArgumentException(errorMessage);

		}
		return seq.toUpperCase();

	}

	/**
	 * This function works the same as cleanSequence, but additionally, if finds
	 * some annotation of a sequence variance like [ADF->ITV], DO NOT applies the
	 * substitution and returns the original sequence .
	 *
	 * @param seq
	 * @return
	 */
	public static String cleanSequenceAndNotApplySequenceVariances(String seq) {
		if (seq == null)
			return null;
		final String seqTmp = seq.trim();

		if (somethingExtrangeInSequence(seqTmp)) {

			// parenthesis or brackets
			final List<String> outside = getOutsideKeepingVariants(seqTmp);
			if (!outside.isEmpty()) {

				final String tmp = appendListIgnoringVariants(outside);
				final String removeBeforeAfterAAs = removeBeforeAfterAAs(tmp);
				if (!removeBeforeAfterAAs.equals(seq)) {
					return cleanSequence(removeBeforeAfterAAs);
				}
			}

		}
		final String errorMessage = "Peptide sequence '" + seq
				+ "' is not supported. Either having not recognizable characteres or in lower case? Has it a non standard PTM enconded on it? PTMs can be encoded as in PEPTID[+45.92]E";
		AssignMass.getInstance(true);
		for (int index = 0; index < seq.length(); index++) {
			final char aa = seq.charAt(index);
			if (!tolerantToFakePeptides && !AssignMass.containsMass(aa)) {
				throw new IllegalArgumentException("'" + aa + "' not recognized. " + errorMessage);
			}
		}
		if (!seq.toUpperCase().equals(seq)) {
			// it has something in lower case

			throw new IllegalArgumentException(errorMessage);

		}
		return seq.toUpperCase();

	}

	public static List<String> getOutside(String seq) {
		int numOpen = 0;
		final List<String> ret = new ArrayList<String>();
		StringBuffer outside = new StringBuffer();
		for (int i = 0; i < seq.length(); i++) {
			final char charAt = seq.charAt(i);
			if (charAt == '[' || charAt == '(') {
				if (!"".equals(outside.toString())) {
					ret.add(outside.toString());
				}
				numOpen++;
			} else if (charAt == ')' || charAt == ']') {
				numOpen--;
				outside = new StringBuffer();
			} else {
				if (numOpen == 0) {
					outside.append(charAt);
				}
			}
		}
		if (!"".equals(outside.toString())) {
			ret.add(outside.toString());
		}
		return ret;
	}

	/**
	 * Makes the same as getOutside but keeps the resulting sequence of the found
	 * variants.<br>
	 * ABDSRE[+89.45]ASDF[AS->PRT]ABC --> returned list: [ABDSRE, ASDF, AS->PRT,
	 * ABC]
	 * 
	 * @param seq
	 * @return
	 */
	public static List<String> getOutsideKeepingVariants(String seq) {
		int numOpen = 0;
		final List<String> ret = new ArrayList<String>();
		StringBuffer outside = new StringBuffer();
		boolean isVariant = false;
		for (int i = 0; i < seq.length(); i++) {
			final char charAt = seq.charAt(i);
			if (charAt == '[' || charAt == '(') {
				if (!"".equals(outside.toString())) {
					ret.add(outside.toString());
				}
				if (charAt == '[' && isVariant(seq, i + 1)) {
					outside = new StringBuffer();
					isVariant = true;
				} else {
					numOpen++;
				}
			} else if (charAt == ')' || charAt == ']') {
				if (isVariant) {
					isVariant = false;
					ret.add(outside.toString());
				} else {
					numOpen--;
				}
				outside = new StringBuffer();
			} else {
				if (numOpen == 0) {
					outside.append(charAt);
				}
			}
		}
		if (!"".equals(outside.toString())) {
			ret.add(outside.toString());
		}
		return ret;
	}

	/**
	 * Returns true if looking to the substring starting from index to the next ']'
	 * it contains '->' as this example: [ASDF->WER]
	 * 
	 * @param seq
	 * @param index
	 * @return
	 */
	private static boolean isVariant(String seq, int index) {
		if (seq.length() <= index) {
			return false;
		}
		final String substring = seq.substring(index);
		final int indexOf = substring.indexOf(']');
		if (indexOf > -1) {
			final String substringVariant = substring.substring(0, indexOf);
			if (substringVariant.contains("->")) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Gets a list of {@link StringPosition} objects with the information inside of
	 * parenthesis or brackets. The information is the text and the position in the
	 * text, not counting the text inside the parenthesis or brackets by itself.<br>
	 * The position is based on 1, that is, starting from 1 in the first character.
	 *
	 * @param seq
	 * @return
	 */
	public static List<StringPosition> getInside(String seq) {
		int numOpen = 0;
		final List<StringPosition> ret = new ArrayList<StringPosition>();
		StringBuffer inside = new StringBuffer();
		int lastNormal = 0;
		int lenthInsides = 0;
		for (int i = 0; i < seq.length(); i++) {
			final char charAt = seq.charAt(i);
			if (charAt == '[' || charAt == '(') {

				if (numOpen > 0) {
					inside.append(charAt);
				} else {
					lastNormal = i - 1 - lenthInsides;
					inside = new StringBuffer();
				}
				numOpen++;
				lenthInsides++;
			} else if (charAt == ']' || charAt == ')') {
				numOpen--;
				if (numOpen > 0) {
					inside.append(charAt);
				} else {
					if (!"".equals(inside.toString())) {
						ret.add(new StringPosition(inside.toString(), lastNormal + 1));
					}
				}
				lenthInsides++;
			} else {
				if (numOpen > 0) {
					inside.append(charAt);
					lenthInsides++;
				}
			}
		}
		return ret;
	}

	public static boolean somethingExtrangeInSequence(String seq) {
		for (int i = 0; i < seq.length(); i++) {
			final char charAt = seq.charAt(i);
			if (charAt == '[' || charAt == ']')
				return true;
			if (charAt == '(' || charAt == ')')
				return true;
			if (charAt == '.')
				return true;
		}
		return false;
	}

	/**
	 * // LLLQQVSL(+80)PELPGEYSMK --> Map<8,+80>
	 *
	 * @param seq
	 * @return a map with the positions and modifications.<br>
	 *         Note that positions start by 1 in the sequence.
	 */
	public static TIntDoubleHashMap getPTMsFromSequence(String rawSeq) {
		return getPTMsFromSequence(rawSeq, null);
	}

	/**
	 * // LLLQQVSL(+80)PELPGEYSMK --> Map<8,+80>
	 *
	 * @param seq
	 * @param isDeltaMass
	 * @return a map with the positions and modifications.<br>
	 *         Note that positions start by 1 in the sequence.
	 */
	public static TIntDoubleHashMap getPTMsFromSequence(String rawSeq, Boolean isDeltaMass) {
		// get the peptide inside '.'
		final String seq = removeBeforeAfterAAs(rawSeq);

		final TIntDoubleHashMap ret = new TIntDoubleHashMap();
		boolean isPTM = false;
		String ptmString = "";
		int position = 1;
		char currentAA = 0;
		for (int i = 0; i < seq.length(); i++) {
			final char charAt = seq.charAt(i);
			if (charAt == '[' || charAt == '(') {
				isPTM = true;
			} else if (charAt == ')' || charAt == ']') {
				isPTM = false;

				try {
					// if the ptm is like [+18.00] is a difference
					// if the ptm is like [180.00] is the modified aminoacid
					boolean modifiedAAMass = false;
					if (isDeltaMass != null) {
						modifiedAAMass = !isDeltaMass;
					} else if (!ptmString.contains("-") && !ptmString.contains("+")) {
						modifiedAAMass = true;
					}
					Double ptm = Double.valueOf(ptmString.replace("+", "").replace("-", ""));

					if (ptmString.contains("-")) {
						ptm = -ptm;
					}
					if (modifiedAAMass) {
						AssignMass.getInstance(true);
						final double AAMass = AssignMass.getMass(currentAA);
						ptm = ptm - AAMass;
					}
					ret.put(position - 1, ptm);
					ptmString = "";
				} catch (final NumberFormatException e) {
					log.warn("Error parsing PTM string '" + ptmString + "' in peptide '" + seq + "'. Ignoring it...");
				}
			} else if (isPTM) {
				ptmString += charAt;
			} else {
				position++;
				currentAA = charAt;
			}
		}

		return ret;
	}

	private static String appendList(List<String> list) {
		final StringBuffer sb = new StringBuffer();
		for (final String string : list) {
			sb.append(string);
		}
		return sb.toString();
	}

	/**
	 * Makes the same as appendList but if one of the elements contains a "->", then
	 * it appends only the part after the "->"
	 * 
	 * @param list
	 * @return
	 */
	private static String appendListApplyingVariants(List<String> list) {
		final StringBuffer sb = new StringBuffer();
		for (final String string : list) {
			if (string.contains("->")) {
				final String[] split = string.split("->");
				if (split.length == 2) {
					sb.append(split[1]);
				}
			} else {
				sb.append(string);
			}
		}
		return sb.toString();
	}

	/**
	 * Makes the same as appendList but if one of the elements contains a "->", then
	 * it appends only the part before the "->"
	 * 
	 * @param list
	 * @return
	 */
	private static String appendListIgnoringVariants(List<String> list) {
		final StringBuffer sb = new StringBuffer();
		for (final String string : list) {
			if (string.contains("->")) {
				final String[] split = string.split("->");
				if (split.length == 2) {
					sb.append(split[0]);
				}
			} else {
				sb.append(string);
			}
		}
		return sb.toString();
	}

	/**
	 * // R.LLLQQVSLPELPGEYSMK.V --> LLLQQVSLPELPGEYSMK
	 *
	 * @param seq
	 * @return
	 */
	private static String removeBeforeAfterAAs(String seq) {
		final String point = ".";
		if (seq.contains(point)) {
			final int firstPoint = seq.indexOf(point);
			final int lastPoint = seq.lastIndexOf(point);

			if (firstPoint != lastPoint) {
				// check that there are no numbers before or afer the points,
				// which would indicate a modification
				if ((firstPoint > 0 && NumberUtils.isNumber(String.valueOf(seq.charAt(firstPoint - 1))))
						|| (seq.length() < firstPoint + 1
								&& !NumberUtils.isNumber(String.valueOf(seq.charAt(firstPoint + 1))))) {
					return seq;
				}
				if ((lastPoint > 0 && NumberUtils.isNumber(String.valueOf(seq.charAt(lastPoint - 1))))
						|| (seq.length() < lastPoint + 1
								&& !NumberUtils.isNumber(String.valueOf(seq.charAt(lastPoint + 1))))) {
					return seq;
				}
				return seq.substring(firstPoint + 1, lastPoint);
			} else if (seq.startsWith(".")) {
				return seq.substring(1);
			} else if (seq.endsWith(".")) {
				return seq.substring(0, seq.length() - 1);
			}
		}
		return seq;
	}

	/**
	 * // R.LLLQQVSLPELPGEYSMK.V --> R
	 *
	 * @param seq
	 * @return
	 */
	public static String getBeforeSeq(String seq) {
		final String point = ".";
		if (seq.contains(point)) {
			final Integer firstPoint = getBeforeSeqPointIndex(seq);
			final Integer lastPoint = getAfterSeqPointIndex(seq);

			if (firstPoint != null && lastPoint != null && firstPoint != lastPoint) {
				final String substring = seq.substring(0, firstPoint);
				return substring;
			}
		}
		return null;
	}

	/**
	 * // R.LLLQQVSL(+80.123)PELPGEYSMK.V --> LLLQQVSL(+80.32)PELPGEYSMK
	 *
	 * @param seq
	 * @return
	 */
	public static String getSequenceInBetween(String seq) {
		final String point = ".";
		if (seq.contains(point)) {
			final int firstPoint = seq.indexOf(point);
			final int lastPoint = seq.lastIndexOf(point);

			if (firstPoint != lastPoint) {
				boolean cutSequence = true;
				// only if the previous and following character of the point are
				// not numbers, otherwise, the point is from a PTM mass
				int previous = firstPoint - 1;
				int following = firstPoint + 1;
				if (previous >= 0) {
					if (NumberUtils.isDigits(String.valueOf(seq.charAt(previous)))) {
						cutSequence = false;
					}
				}
				if (following < seq.length()) {
					if (NumberUtils.isDigits(String.valueOf(seq.charAt(following)))) {
						cutSequence = false;
					}
				}
				previous = lastPoint - 1;
				following = lastPoint + 1;
				if (previous >= 0) {
					if (NumberUtils.isDigits(String.valueOf(seq.charAt(previous)))) {
						cutSequence = false;
					}
				}
				if (following < seq.length()) {
					if (NumberUtils.isDigits(String.valueOf(seq.charAt(following)))) {
						cutSequence = false;
					}
				}
				if (cutSequence) {
					final String substring = seq.substring(firstPoint + 1, lastPoint);
					return substring;
				} else {
					return seq;
				}
			}
		}
		return seq;
	}

	/**
	 * // R.LLLQQVSLPELPGEYSMK.V --> V
	 *
	 * @param seq
	 * @return
	 */
	public static String getAfterSeq(String seq) {
		final String point = ".";
		if (seq.contains(point)) {
			final Integer firstPoint = getBeforeSeqPointIndex(seq);
			final Integer lastPoint = getAfterSeqPointIndex(seq);

			if (firstPoint != null && lastPoint != null && firstPoint != lastPoint) {
				final String substring = seq.substring(lastPoint + 1, seq.length());
				return substring;
			}
		}
		return null;
	}

	/**
	 * R.LLLQQVSLPELPGEYSMK.V --> 1 <br>
	 * LLLQQVSLPELPGEYSMK --> null
	 *
	 * @param sequence
	 * @return
	 */
	private static Integer getBeforeSeqPointIndex(String sequence) {
		final TIntArrayList allPositionsOf = StringUtils.allPositionsOf(sequence, '.');
		for (int i = 0; i < allPositionsOf.size(); i++) {
			final Integer index = allPositionsOf.get(i) - 1;
			if (index < sequence.length() - 1) {
				// is followed by a number?
				final char charAt = sequence.charAt(index + 1);
				if (!Character.isDigit(charAt)) {
					return index;
				}

			}
		}
		return null;
	}

	/**
	 * R.LLLQQVSLPELPGEYSMK.V --> 1 <br>
	 * LLLQQVSLPELPGEYSMK --> null
	 *
	 * @param sequence
	 * @return
	 */
	private static Integer getAfterSeqPointIndex(String sequence) {
		final TIntArrayList allPositionsOf = StringUtils.allPositionsOf(sequence, '.');
		for (int i = allPositionsOf.size() - 1; i >= 0; i--) {
			final Integer index = allPositionsOf.get(i) - 1;
			if (index < sequence.length() - 1) {
				// is followed by a number? then is a modification like (80.009)
				final char charAt = sequence.charAt(index + 1);
				if (!Character.isDigit(charAt)) {
					return index;
				}

			}
		}
		return null;
	}

	public static void main(String[] args) {
		System.out.println(FastaParser.getACC("sp|P54399-2|PDI_DROME"));
		System.out
				.println(FastaParser.getACC("asdfaksjdfhasjkdf_ASDFasdf_ASDfsp|Pasdf22532|SPR2D_sdfsdfHUMAN asdfasdf"));

		System.out.println(FastaParser.getNCBIACC(">gi|194150008|gb|EDW65699.1| GJ18760 [Drosophila virilis"));
		System.out.println(FastaParser.getDescription(">gi|194150008|gb|EDW65699.1| GJ18760 [Drosophila virilis]"));

		System.out.println(FastaParser.getACC(">gi|194150008|gb|EDW65699.1| GJ18760 [Drosophila virilis]"));
		System.out.println(FastaParser.getDescription(
				">sp|P22532|SPR2D_HUMAN Small proline-rich protein 2D OS=Homo sapiens GN=SPRR2D PE=2 SV=2"));
		System.out.println(FastaParser
				.getACC(">sp|P22532|SPR2D_HUMAN Small proline-rich protein 2D OS=Homo sapiens GN=SPRR2D PE=2 SV=2"));
		System.out.println(FastaParser.getACC("IPI:IPI00930860.1|REFSEQ:NP_001155007"));

		System.out.println(FastaParser.getACC("IPI:IPI00930860|REFSEQ:NP_001155007"));

		System.out.println(FastaParser
				.getOrganismNameFromFastaHeader("CG7415 (Fragment) OS=Drosophila simulans PE=2 SV=1 ", "CG7415"));

		System.out.println(FastaParser
				.getOrganismNameFromFastaHeader("CG7415 (Fragment) OS=Drosophila simulans PE=2 SV=1 ", "CG7415"));

		System.out.println(FastaParser.getOrganismNameFromFastaHeader(
				">tr|A3SA23|A3SA23_9RHOB TonB dependent, hydroxamate-type ferrisiderophore, outer membrane receptor OS=Sulfitobacter sp. EE-36 GN=EE36_08023 PE=3 SV=1",
				"asdf"));

		System.out.println(FastaParser
				.getACC("A6NKZ8\tA6NKZ8;  RecName: Full=Putative tubulin beta chain-like protein ENSP00000290377;"));

		final List<String> ret = FastaParser.getOutsideKeepingVariants("ABDSRE[+89.45]ASDF[AS->PRT]ABC");
		for (final String string : ret) {
			System.out.println(string);
		}
		System.out.println(FastaParser.cleanSequenceAndApplySequenceVariances("ABDSRE[+89.45]ASDF[AS->PRT]ABC"));
		System.out.println(FastaParser.cleanSequenceAndApplySequenceVariances("ABDSRE[+89.45]ASDF[AS->]ABC"));
	}

	public static boolean isContaminant(String accession) {

		if (getACC(accession).getAccession().startsWith(CONTAMINANT_PREFIX)) {
			return true;
		}
		return false;
	}

	public static boolean isReverse(String id) {

		if (getACC(id).getAccession().contains("Reverse")) {
			return true;
		}
		return false;
	}

	/**
	 * Gets the scan number from a string like 'rawfile.scan.scan.charge'
	 *
	 * @param psmId
	 * @return
	 */
	public static String getScanFromPSMIdentifier(String psmId) {
		if (psmId == null) {
			return null;
		}
		if (psmId.contains(".")) {
			final String[] split = psmId.split("\\.");
			return split[split.length - 3];
		}
		return null;
	}

	/**
	 * Gets the scan2 number from a string like 'rawfile.scan.scan2.charge'
	 *
	 * @param psmId
	 * @return
	 */
	public static String getSecondScanFromPSMIdentifier(String psmId) {
		if (psmId == null) {
			return null;
		}
		if (psmId.contains(".")) {
			final String[] split = psmId.split("\\.");
			return split[split.length - 2];
		}
		return null;
	}

	/**
	 * Gets the file name from a string like 'rawfile.scan.scan.charge'
	 *
	 * @param psmId
	 * @return
	 */
	public static String getFileNameFromPSMIdentifier(String psmId) {
		if (psmId == null) {
			return null;
		}
		if (psmId.contains(".")) {
			final int indexOf = psmId.indexOf(".");
			final String fileName = psmId.substring(0, indexOf);
			return fileName;
		} else {
			if (psmId.contains("-")) {
				final int index = psmId.lastIndexOf("-");
				return psmId.substring(0, index);
			}
		}
		return null;
	}

	public static Integer getChargeStateFromPSMIdentifier(String psmId) {
		if (psmId == null) {
			return null;
		}
		if (psmId.contains(".")) {
			final int indexOf = psmId.lastIndexOf(".");
			final String chargeStateString = psmId.substring(indexOf + 1);
			try {
				return Integer.valueOf(chargeStateString);
			} catch (final NumberFormatException e) {

			}
		}
		return null;
	}

	private final static Pattern isoformPattern = Pattern.compile("\\w+\\-\\w*");

	/**
	 * If an accession like P12334-1 is found, P12334 is returned.<br>
	 * If an accession like P12334-4 is found, P12334 is returned.<br>
	 * If an accession like P12345 is found, the same P12345 is returned
	 *
	 * @param uniprotAcc
	 * @return
	 */
	public static String getNoIsoformAccession(String uniprotAcc) {

		if (uniprotAcc.indexOf("-") >= 0 && isoformPattern.matcher(uniprotAcc).matches()) {
			return uniprotAcc.substring(0, uniprotAcc.indexOf("-"));
		}
		return uniprotAcc;
	}

	/**
	 * If an accession like P12334-1 is found, 1 is returned.<br>
	 * If an accession like P12334-4 is found, 4 is returned.<br>
	 * If an accession like P12345 is found, null
	 *
	 * @param uniprotAcc
	 * @return
	 */
	public static String getIsoformVersion(String uniprotAcc) {

		if (uniprotAcc.indexOf("-") >= 0 && isoformPattern.matcher(uniprotAcc).matches()) {
			return uniprotAcc.substring(uniprotAcc.indexOf("-") + 1);
		}
		return null;
	}

	/**
	 * Returns true if the acc contains _mutated_ _conflict_ or _variant_
	 * 
	 * @param acc
	 * @return
	 */
	public static boolean isProteoform(String acc) {
		if (acc.contains(mutated) || acc.contains(conflict) || acc.contains(variant)) {
			return true;
		}
		return false;
	}

	public static boolean isSwissProt(String fastaHeader) {
		final String acc = getUniProtACC(fastaHeader);
		if (acc == null) {
			return false;
		}
		final Matcher matcherSP = UNIPROT_SP_ACC_TMP.matcher(fastaHeader);
		final Matcher matcherTR = UNIPROT_TR_ACC_TMP.matcher(fastaHeader);
		if (fastaHeader.contains("|")) {
			if (matcherSP.find()) {
				final String group = matcherSP.group(1);
				if (!"".equals(group))
					return true;
			} else if (matcherTR.find()) {
				final String group = matcherTR.group(1);
				if (!"".equals(group))
					return false;
			}
		}
		return false;
	}

	public static boolean isFragment(String fullHeaderWithAddenda) {
		if (fullHeaderWithAddenda.toLowerCase().contains("(fragment)")) {
			return true;
		}
		return false;
	}
}
