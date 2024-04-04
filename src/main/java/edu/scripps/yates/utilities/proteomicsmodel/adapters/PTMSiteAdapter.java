package edu.scripps.yates.utilities.proteomicsmodel.adapters;

import edu.scripps.yates.utilities.proteomicsmodel.PTMPosition;
import edu.scripps.yates.utilities.proteomicsmodel.PTMSite;
import edu.scripps.yates.utilities.proteomicsmodel.Score;
import edu.scripps.yates.utilities.proteomicsmodel.factories.PTMSiteEx;

public class PTMSiteAdapter implements edu.scripps.yates.utilities.pattern.Adapter<PTMSite> {
	private final String aa;
	private final int position;
	private final Score score;
	private final PTMPosition ptmPosition;

	public PTMSiteAdapter(String aa, int position, PTMPosition ptmPosition, Score score) {
		this.aa = aa;
		this.position = position;
		this.score = score;
		this.ptmPosition = ptmPosition;
	}

	public PTMSiteAdapter(String aa, int position, PTMPosition ptmPosition) {
		this(aa, position, ptmPosition, null);
	}

	@Override
	public PTMSite adapt() {
		final PTMSiteEx ptmSite = new PTMSiteEx(aa, position, ptmPosition);
		if (score != null)
			ptmSite.setScore(score);
		return ptmSite;
	}

}
