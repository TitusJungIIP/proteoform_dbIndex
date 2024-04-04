package edu.scripps.yates.utilities.annotations;

import java.util.Collection;
import java.util.Map;

import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;

public interface UniprotProteinLocalRetrieverInterface {

	Map<String, Entry> getAnnotatedProtein(String uniprotVersion, String accession);

	Map<String, Entry> getAnnotatedProteins(String uniprotVersion, Collection<String> accessions);

	Map<String, Entry> getAnnotatedProteins(String uniprotVersion, Collection<String> accessions,
			boolean retrieveFastaIsoforms, boolean lookForIsoformsFromMainForms);

	boolean isCacheEnabled();

	void setCacheEnabled(boolean cacheEnabled);

	void setRetrieveFastaIsoforms(boolean b);

}
