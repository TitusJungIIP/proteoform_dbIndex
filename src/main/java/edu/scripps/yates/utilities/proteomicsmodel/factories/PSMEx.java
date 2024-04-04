package edu.scripps.yates.utilities.proteomicsmodel.factories;

import java.io.Serializable;

import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.proteomicsmodel.AbstractPSM;

public class PSMEx extends AbstractPSM implements Serializable {
	/**
	 *
	 */
	private static final long serialVersionUID = -1713998292042049510L;

	private final Logger log = Logger.getLogger(PSMEx.class);

	public PSMEx(String psmID, String sequence, String fullSequence, boolean distinguishModifiedSequence,
			boolean chargeStateSensible) {
		super(distinguishModifiedSequence, chargeStateSensible);
		setKey(psmID);
		setIdentifier(psmID);
		if (sequence != null) {
			setSequence(sequence);
		}
		if (fullSequence != null) {
			setFullSequence(fullSequence);
		}
	}

}
