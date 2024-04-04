package edu.scripps.yates.utilities.ontologies;

import java.util.List;

import org.apache.log4j.Logger;

import uk.ac.ebi.pride.utilities.ols.web.service.client.OLSClient;
import uk.ac.ebi.pride.utilities.ols.web.service.config.OLSWsConfig;
import uk.ac.ebi.pride.utilities.ols.web.service.model.Identifier;
import uk.ac.ebi.pride.utilities.ols.web.service.model.Term;

public class OntologyQueryInterface {
	private final static Logger log = Logger.getLogger(OntologyQueryInterface.class);
	private static final OLSClient olsClient = new OLSClient(new OLSWsConfig());

	public static boolean containsAsParent(String termID, String parentID, boolean retryiffails, int distance) {
		try {
			if (!termID.contains(":")) {
				throw new IllegalArgumentException(
						"TermID should have the ontology ID followed by ':' and the term ID number");
			}
			final String ontologyID = termID.split(":")[0];
			final List<Term> parent = olsClient.getTermParents(new Identifier(termID, Identifier.IdentifierType.OBO),
					ontologyID, distance);
			if (parent.isEmpty()) {
				throw new IllegalArgumentException("term " + termID + " has no parents?? Maybe the id is wrong");
			}
			for (final Term term : parent) {
				if (term.getOboId().getIdentifier().equals(parentID)) {
					return true;
				}
			}

		} catch (final Exception e) {
			e.printStackTrace();
			log.error("Error getting term " + termID);
			if (retryiffails) {
				return containsAsParent(termID, parentID, false, distance);
			}
		}
		return false;
	}

}
