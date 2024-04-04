package edu.scripps.yates.utilities.grouping;

import java.util.List;

public interface GroupableProtein {

	public List<GroupablePeptide> getGroupablePeptides();

	public ProteinGroup getProteinGroup();

	public String getUniqueID();

	public String getAccession();

	public void setEvidence(ProteinEvidence evidence);

	public ProteinEvidence getEvidence();

	public void setProteinGroup(ProteinGroup proteinGroup);

}
