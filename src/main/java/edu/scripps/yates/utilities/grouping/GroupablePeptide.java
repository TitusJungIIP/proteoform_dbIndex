package edu.scripps.yates.utilities.grouping;

import java.util.List;

public interface GroupablePeptide {

	String getSequence();

	String getIdentifier();

	void setRelation(PeptideRelation relation);

	PeptideRelation getRelation();

	List<GroupableProtein> getGroupableProteins();

}
