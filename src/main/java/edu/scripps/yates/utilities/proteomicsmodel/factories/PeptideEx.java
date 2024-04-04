package edu.scripps.yates.utilities.proteomicsmodel.factories;

import java.io.Serializable;

import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.proteomicsmodel.AbstractPeptide;

public class PeptideEx extends AbstractPeptide implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = -6127166815987338850L;

	private static Logger log = Logger.getLogger(PeptideEx.class);

	/**
	 * 
	 * @param fullSequence
	 * @param key          the key used to store the peptide (can be charge-state
	 *                     sensible and/or PTM sensible)
	 */
	public PeptideEx(String fullSequence, String key) {
		super(key);
		setFullSequence(fullSequence);

	}

}
