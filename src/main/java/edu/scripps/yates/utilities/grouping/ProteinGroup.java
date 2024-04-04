package edu.scripps.yates.utilities.grouping;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import gnu.trove.set.hash.THashSet;

public class ProteinGroup extends ArrayList<GroupableProtein> {
	/**
	 *
	 */
	private static final long serialVersionUID = 1516424786690373161L;
	private static final Logger log = Logger.getLogger(ProteinGroup.class);
	private ProteinEvidence evidence;
	private List<String> accessions;
	private String key;

	public ProteinGroup() {
		super();
	}

	public ProteinGroup(ProteinEvidence e) {
		super();
		evidence = e;
	}

	public ProteinGroup(ProteinGroupInference iProteinGroup) {
		if (iProteinGroup == null)
			throw new IllegalArgumentException("group is null");

		for (final InferenceProtein inferenceProtein : iProteinGroup) {
			final List<GroupableProtein> proteinsMerged = inferenceProtein.getProteinsMerged();
			for (final GroupableProtein protein : proteinsMerged) {
				protein.setProteinGroup(this);
				protein.setEvidence(inferenceProtein.getEvidence());
			}
			this.addAll(proteinsMerged);
		}
		evidence = iProteinGroup.getEvidence();

	}

	public String getKey() {
		if (key != null)
			return key;
		String ret = "";

		final List<String> accessions2 = getAccessions();

		for (final String accession : accessions2) {
			if (!"".equals(ret)) {
				ret = ret + ",";
			}
			ret = ret + accession;
		}
		key = ret;
		if (evidence != null) {
			key += "[" + evidence.toString() + "]";
		}
		return key;
	}

	// public ProteinGroup() {
	// this(ProteinEvidence.AMBIGUOUSGROUP);
	// }

	@Override
	public String toString() {
		final StringBuilder sb = new StringBuilder();
		String shareOrContain = "sharing";
		if (size() == 1)
			shareOrContain = "containing";
		final Set<String> seqs = new THashSet<String>();
		final List<GroupablePeptide> psMs = getPSMs();
		for (final GroupablePeptide groupablePSM : psMs) {
			seqs.add(groupablePSM.getSequence());
		}
		String evidenceString = "";
		if (getEvidence() != null) {
			evidenceString = getEvidence().name() + " - ";
		}
		sb.append("Group: " + evidenceString + shareOrContain + " " + psMs.size() + " PSMs (" + seqs.size()
				+ " different sequences)\n");
		int memberGroup = 1;
		for (final GroupableProtein protein : this) {
			final StringBuilder sb2 = new StringBuilder();
			for (final GroupablePeptide groupablePSM : protein.getGroupablePeptides()) {
				sb2.append(groupablePSM.getSequence() + " ");
			}
			sb.append("PRT\t\t" + memberGroup++ + " " + protein.getAccession() + "\t" + sb2.toString() + "\t"
					+ protein.getEvidence() + "\n");
		}
		return sb.toString();
	}

	@Override
	public boolean equals(Object object) {
		// return shareOneProtein(object);
		return shareAllProteins(object);
		// return theyformJustOneGroup(object);
	}

	private boolean theyformJustOneGroup(Object object) {
		if (object instanceof ProteinGroup) {
			if (shareOneProtein(object)) {

				final List<GroupableProtein> proteins = new ArrayList<GroupableProtein>();
				proteins.addAll(this);
				proteins.addAll((ProteinGroup) object);
				final PAnalyzer pa = new PAnalyzer(false);
				final List<ProteinGroup> groups = pa.run(proteins);
				if (groups.size() == 1)
					return true;
			}
			return false;
		}
		return super.equals(object);
	}

	public boolean shareAllProteins(Object object) {
		if (object instanceof ProteinGroup) {
			final ProteinGroup pg2 = (ProteinGroup) object;

			if (getKey().equals(pg2.getKey())) {
				// if (this.evidence == pg2.evidence)
				return true;
			}
			return false;
		}
		// else if (object instanceof ProteinGroupOccurrence) {
		// ProteinGroupOccurrence pgo2 = (ProteinGroupOccurrence) object;
		// if (equals(pgo2.getFirstOccurrence()))
		// return true;
		// return false;
		// }
		return super.equals(object);
	}

	/**
	 * This method will determine how comparisons are made between proteinGroups! In
	 * this case, two groups are equals if share at least one protein.
	 */
	public boolean shareOneProtein(Object object) {
		if (object instanceof ProteinGroup) {
			final ProteinGroup pg2 = (ProteinGroup) object;

			// At least share one protein
			for (final String acc : getAccessions()) {
				if (pg2.getAccessions().contains(acc))
					return true;
			}
			return false;
		}
		// else if (object instanceof ProteinGroupOccurrence) {
		// ProteinGroupOccurrence pgo2 = (ProteinGroupOccurrence) object;
		// if (equals(pgo2.getFirstOccurrence()))
		// return true;
		// return false;
		// }
		return super.equals(object);
	}

	/*
	 * public int updateMinimum() { } List<ProteinGroup> getRecursive( ProteinGroup
	 * group, Iterator<Protein> it ) { List<ProteinGroup> res = new
	 * ArrayList<ProteinGroup>(); if( group == null ) group = new ProteinGroup();
	 * res.add(group); if( !it.hasNext() ) return res; ProteinGroup group2 =
	 * (ProteinGroup)group.clone(); group2.add(it.next());
	 * res.addAll(getRecursive(group, it)); res.addAll(getRecursive(group2, it));
	 * return res; }
	 */

	public List<String> getAccessions() {

		if (accessions != null) {
			return accessions;
		}
		accessions = new ArrayList<String>();
		for (final GroupableProtein protein : this) {
			if (!accessions.contains(protein.getAccession())) {
				accessions.add(protein.getAccession());
			}
		}
		Collections.sort(accessions);
		return accessions;

	}

	public String getAccessionString(String separator) {
		final List<String> accs = getAccessions();
		final StringBuilder sb = new StringBuilder();
		for (final String acc : accs) {
			if (!"".equals(sb.toString())) {
				sb.append(separator);
			}
			sb.append(acc);
		}
		return sb.toString();
	}

	@Override
	public void add(int index, GroupableProtein element) {
		accessions = null;
		super.add(index, element);
		key = null;
	}

	public ProteinEvidence getEvidence() {
		return evidence;
	}

	public void setEvidence(ProteinEvidence evidence) {
		this.evidence = evidence;

	}

	/**
	 * Gets all peptides from the proteins of the group
	 *
	 * @return
	 */
	public List<GroupablePeptide> getPSMs() {
		// if (this.peptides == null || this.peptides.isEmpty()) {
		final List<GroupablePeptide> ret = new ArrayList<GroupablePeptide>();
		final Set<String> peptideIds = new THashSet<String>();
		for (final GroupableProtein protein : this) {
			final List<GroupablePeptide> psms = protein.getGroupablePeptides();
			if (psms != null)
				for (final GroupablePeptide psm : psms) {
					if (!peptideIds.contains(psm.getIdentifier())) {
						peptideIds.add(psm.getIdentifier());
						ret.add(psm);
					}
				}

		}
		// }

		return ret;
	}

	@Override
	public boolean add(GroupableProtein e) {
		key = null;
		return super.add(e);
	}

	public String getAccessionStringByEvidence(Map<String, Entry> uniprotEntries, String separator) throws IOException {

		final List<Boolean> validArray = filterAccessionsByEvidence(uniprotEntries);
		int index = 0;
		final StringBuilder sb = new StringBuilder();
		for (final String acc : getAccessions()) {
			if (!validArray.get(index++)) {
				continue;
			}
			if (!"".equals(sb.toString())) {
				sb.append(separator);
			}
			sb.append(acc);

		}
		return sb.toString();
	}

	public List<String> getAccessionsFilteredByEvidence(Map<String, Entry> uniprotEntries, String separator) {

		final List<Boolean> validArray = filterAccessionsByEvidence(uniprotEntries);
		int index = 0;
		final List<String> ret = new ArrayList<String>();
		for (final String acc : getAccessions()) {
			if (!validArray.get(index++)) {
				continue;
			}
			ret.add(acc);

		}
		return ret;
	}

	private List<Boolean> filterAccessionsByEvidence(Map<String, Entry> uniprotEntries) {

		final List<Boolean> ret = new ArrayList<>();
		final List<Boolean> groupEvidenceArray = new ArrayList<>();
		final List<Boolean> uniprotEvidenceArray = new ArrayList<>();
		// only swissprot is valid
		boolean thereIsASwissProt = false;
		boolean thereIsAConclusiveProt = false;
		for (final String acc : getAccessions()) {
			boolean valid = false;
			try {

				final ProteinEvidence evidence = getEvidence(acc);
				if (evidence == ProteinEvidence.CONCLUSIVE) {
					groupEvidenceArray.add(true);
					thereIsAConclusiveProt = true;

					valid = true;
				} else if (evidence == ProteinEvidence.NONCONCLUSIVE) {
					groupEvidenceArray.add(false);
					valid = false;
				} else {
					groupEvidenceArray.add(false);
					if (uniprotEntries.containsKey(acc)) {
						final Entry protein2 = uniprotEntries.get(acc);
						if (protein2 != null) {
							final String dataset = protein2.getDataset();
							if (dataset != null) {
								if (dataset.toLowerCase().equals("swiss-prot")) {
									thereIsASwissProt = true;
									uniprotEvidenceArray.add(true);
									valid = true;
								} else {
									uniprotEvidenceArray.add(false);
								}
							}
						} else {
							uniprotEvidenceArray.add(false);
						}
					} else {
						uniprotEvidenceArray.add(false);
					}
				}
			} finally {
				ret.add(valid);
			}
		}
		boolean allSwissprot = true;
		for (final Boolean uniprotEvidence : uniprotEvidenceArray) {
			if (!uniprotEvidence) {
				allSwissprot = false;
			}
		}

		if (!thereIsASwissProt && !thereIsAConclusiveProt) {
			ret.set(0, true);
		}
		if (allSwissprot) {
			if (thereIsAConclusiveProt) {
				// only report the conclusive ones
				return groupEvidenceArray;
			} else {
				// do not report nonconclusive

				if (uniprotEvidenceArray.size() == ret.size()) {
					return uniprotEvidenceArray;
				}
			}
		}

		return ret;
	}

	private ProteinEvidence getEvidence(String acc) {

		for (final GroupableProtein groupableProtein : this) {
			if (groupableProtein.getAccession().equals(acc)) {
				return groupableProtein.getEvidence();
			}
		}

		return null;
	}
}
