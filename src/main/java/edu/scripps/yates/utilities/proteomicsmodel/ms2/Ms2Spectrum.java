package edu.scripps.yates.utilities.proteomicsmodel.ms2;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import edu.scripps.yates.utilities.proteomicsmodel.spectra.Spectrum;
import gnu.trove.map.hash.THashMap;

public class Ms2Spectrum extends Spectrum {
	private final Map<String, String> iProperties = new THashMap<String, String>();
	private final List<String> iPropertyList = new ArrayList<String>();

	public Ms2Spectrum() {
		super();
	}

	public Ms2Spectrum(String scanNumber, double rt, double precursorMZ) {
		super(scanNumber, rt, precursorMZ);
	}

	public void addIProperty(String propertyName, String propertyValue) {
		this.iPropertyList.add(propertyName);
		this.iProperties.put(propertyName, propertyValue);
	}

	public Map<String, String> getIProperties() {
		return iProperties;
	}

	public List<String> getIPropertyList() {
		return iPropertyList;
	}
}
