package edu.scripps.yates.utilities.proteomicsmodel.factories;

import java.io.Serializable;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.utilities.proteomicsmodel.PTM;
import edu.scripps.yates.utilities.proteomicsmodel.PTMPosition;
import edu.scripps.yates.utilities.proteomicsmodel.PTMSite;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;
import uk.ac.ebi.pride.utilities.pridemod.ModReader;
import uk.ac.ebi.pride.utilities.pridemod.model.Specificity;

public class PTMEx implements PTM, Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = -4442800986087569468L;
	private String name;
	private final double massShift;
	private String cvId;
	private List<PTMSite> ptmSites;
	private String residues;
	private static final ModReader modReader = ModReader.getInstance();
	public static final Double PRECISION = 0.001;
	public static final String UNKNOWN = "Unknown";
	public static final DecimalFormat formatter = new DecimalFormat("#.###");
	public static Map<String, List<uk.ac.ebi.pride.utilities.pridemod.model.PTM>> ptmsByMonoDeltaMass = new THashMap<String, List<uk.ac.ebi.pride.utilities.pridemod.model.PTM>>();

	public PTMEx(String name, double massShift) {
		this.name = name;
		this.massShift = massShift;
	}

	public PTMEx(double massShift, String aa, int position) {
		this(massShift, aa, position, position == 0 ? PTMPosition.NTERM : PTMPosition.NONE);
	}

	public PTMEx(double massShift, String aa, int position, PTMPosition ptmPosition) {
		List<uk.ac.ebi.pride.utilities.pridemod.model.PTM> ptms = null;
		final String deltaMassString = formatter.format(massShift);
		if (ptmsByMonoDeltaMass.containsKey(deltaMassString)) {
			ptms = ptmsByMonoDeltaMass.get(deltaMassString);
		} else {
			ptms = modReader.getPTMListByMonoDeltaMass(massShift, PRECISION / 10.0);
			if (ptms == null || ptms.isEmpty()) {
				ptms = modReader.getPTMListByMonoDeltaMass(massShift, PRECISION);
			}
		}
		if (ptms != null && !ptms.isEmpty()) {
			ptmsByMonoDeltaMass.put(deltaMassString, ptms);
			for (final uk.ac.ebi.pride.utilities.pridemod.model.PTM ptm : ptms) {

				if (ptm.getSpecificityCollection() != null && !ptm.getSpecificityCollection().isEmpty()) {
					boolean specificityOK = false;
					for (final Specificity specificity : ptm.getSpecificityCollection()) {
						if (specificity.getName().name().equalsIgnoreCase(aa)) {
							specificityOK = true;
							break;
						}
					}
					if (specificityOK) {
						name = ptm.getName();
						cvId = ptm.getAccession();
						break;
					}
				}
			}
		} else {
			if (Double.compare(0.0, massShift) != 0) {
				name = deltaMassString;
			}
		}
		if (name == null) {
			name = UNKNOWN;
		}

		addPtmSite(new PTMSiteEx(aa, position, ptmPosition));
		this.massShift = massShift;
	}

	@Override
	public Double getMassShift() {
		return massShift;
	}

	@Override
	public String getName() {
		return name;
	}

	@Override
	public String getCVId() {
		return cvId;
	}

	@Override
	public List<PTMSite> getPTMSites() {
		return ptmSites;
	}

	/**
	 * @param cvId the cvId to set
	 */
	public void setCvId(String cvId) {
		this.cvId = cvId;
	}

	/**
	 * @param ptmSites the ptmSites to set
	 */
	public void setPtmSites(List<PTMSite> ptmSites) {
		this.ptmSites = ptmSites;
	}

	public void addPtmSite(PTMSite ptmSite) {
		if (ptmSites == null)
			ptmSites = new ArrayList<PTMSite>();
		ptmSites.add(ptmSite);
	}

	@Override
	public String getResidues() {
		if (residues == null) {
			final Set<String> aas = new THashSet<String>();
			final List<PTMSite> ptmSites2 = getPTMSites();
			if (ptmSites2.size() == 1) {
				residues = ptmSites2.iterator().next().getAA();
				return residues;
			}
			for (final PTMSite ptmSite : ptmSites2) {
				if (!aas.contains(ptmSite.getAA())) {
					aas.add(ptmSite.getAA());
				}
			}
			final StringBuilder sb = new StringBuilder();
			aas.stream().forEach(aa -> sb.append(aa));
			residues = sb.toString();
		}
		return residues;
	}

}
