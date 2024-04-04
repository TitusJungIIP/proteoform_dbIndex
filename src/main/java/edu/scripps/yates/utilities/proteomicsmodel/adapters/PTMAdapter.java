package edu.scripps.yates.utilities.proteomicsmodel.adapters;

import java.util.List;

import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.proteomicsmodel.PTM;
import edu.scripps.yates.utilities.proteomicsmodel.PTMPosition;
import edu.scripps.yates.utilities.proteomicsmodel.Score;
import edu.scripps.yates.utilities.proteomicsmodel.factories.PTMEx;
import gnu.trove.map.TDoubleObjectMap;
import gnu.trove.map.hash.TDoubleObjectHashMap;
import uk.ac.ebi.pride.utilities.pridemod.ModReader;

public class PTMAdapter implements edu.scripps.yates.utilities.pattern.Adapter<PTM> {
	private final Double massShift;
	private final int position;
	private final String aa;
	private final Score score;
	private final static Logger log = Logger.getLogger(PTMAdapter.class);
	private final uk.ac.ebi.pride.utilities.pridemod.model.PTM prideModPTM;
	private static final String MOD0 = "MOD:00000";
	private final PTMPosition ptmPosition;
	private static TDoubleObjectMap<uk.ac.ebi.pride.utilities.pridemod.model.PTM> staticPtmsByMonoDelta = new TDoubleObjectHashMap<uk.ac.ebi.pride.utilities.pridemod.model.PTM>();

	public PTMAdapter(double massShift, String aa, int position, PTMPosition ptmPosition) {
		this(massShift, aa, position, ptmPosition, null);
	}

	public PTMAdapter(double massShift, char aa, int position, PTMPosition ptmPosition) {
		this(massShift, String.valueOf(aa), position, ptmPosition);
	}

	public PTMAdapter(String name, char aa, int position, PTMPosition ptmPosition) {
		this(name, String.valueOf(aa), position, ptmPosition);
	}

	public PTMAdapter(String name, String aa, int position, PTMPosition ptmPosition) {
		this(name, aa, position, ptmPosition, null);
	}

	public PTMAdapter(String name, char aa, int position, PTMPosition ptmPosition, Score score) {
		this(name, String.valueOf(aa), position, ptmPosition, score);
	}

	public PTMAdapter(String name, String aa, int position, PTMPosition ptmPosition, Score score) {
		this.position = position;
		this.ptmPosition = ptmPosition;
		this.aa = aa;
		this.score = score;
		final ModReader modReader = ModReader.getInstance();
		final List<uk.ac.ebi.pride.utilities.pridemod.model.PTM> ptmListByMonoDeltaMass = modReader
				.getPTMListByPatternName(name);
		if (ptmListByMonoDeltaMass != null && !ptmListByMonoDeltaMass.isEmpty()) {
			prideModPTM = ptmListByMonoDeltaMass.get(0);
			massShift = prideModPTM.getMonoDeltaMass();
		} else {
			log.warn("PTM with name '" + name + "' is not recognized");
			prideModPTM = null;
			massShift = null;
		}
	}

	public PTMAdapter(double massShift, char aa, int position, PTMPosition ptmPosition, Score score) {
		this(massShift, String.valueOf(aa), position, ptmPosition, score);
	}

	public PTMAdapter(double massShift, String aa, int position, PTMPosition ptmPosition, Score score) {
		this.massShift = massShift;
		this.position = position;
		this.ptmPosition = ptmPosition;
		this.aa = aa;
		this.score = score;
		final ModReader modReader = ModReader.getInstance();

		if (staticPtmsByMonoDelta.containsKey(massShift)) {
			prideModPTM = staticPtmsByMonoDelta.get(massShift);
		} else {
			final List<uk.ac.ebi.pride.utilities.pridemod.model.PTM> ptmListByMonoDeltaMass = modReader
					.getPTMListByMonoDeltaMass(massShift, PTMEx.PRECISION);
			if (ptmListByMonoDeltaMass != null && !ptmListByMonoDeltaMass.isEmpty()) {

				prideModPTM = ptmListByMonoDeltaMass.get(0);
			} else {
				prideModPTM = null;
			}
			staticPtmsByMonoDelta.put(massShift, prideModPTM);
		}
	}

	@Override
	public PTM adapt() {
		String name = null;
		if (massShift != null && Double.compare(0.0, massShift) != 0) {
			name = PTMEx.formatter.format(massShift);
		}
		if (name == null) {
			name = PTMEx.UNKNOWN;
		}
		if (prideModPTM != null) {
			name = prideModPTM.getName();
		}
		final PTMEx ptm = new PTMEx(name, massShift);
		if (prideModPTM != null) {
			ptm.setCvId(prideModPTM.getAccession());
		}

		ptm.addPtmSite(new PTMSiteAdapter(aa, position, ptmPosition, score).adapt());
		return ptm;
	}

}
