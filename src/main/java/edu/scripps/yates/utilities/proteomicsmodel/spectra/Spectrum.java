package edu.scripps.yates.utilities.proteomicsmodel.spectra;

import edu.scripps.yates.utilities.masses.AssignMass;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.TDoubleDoubleMap;
import gnu.trove.map.hash.TDoubleDoubleHashMap;

public class Spectrum {
	private String scanNumber;
	private double rt;
	private double precursorMZ;
	private double precursorIntensity;
	private int z;
	private final TDoubleDoubleMap peaks = new TDoubleDoubleHashMap();
	private final TDoubleList masses = new TDoubleArrayList();

	public Spectrum() {

	}

	public void setScanNumber(String scanNumber) {
		this.scanNumber = scanNumber;
	}

	public void setRt(double rt) {
		this.rt = rt;
	}

	/**
	 * Sets the m/z mass which is m/z = (MW + nH+)/n
	 * 
	 * @param precursorMZ
	 */
	public void setPrecursorMZ(double precursorMZ) {
		this.precursorMZ = precursorMZ;
	}

	public Spectrum(String scanNumber, double rt, double precursorMZ) {
		this.scanNumber = scanNumber;
		this.rt = rt;
		this.precursorMZ = precursorMZ;
	}

	public void addPeak(double mass, double peak) {
		this.peaks.put(mass, peak);
		this.masses.add(mass);
	}

	public String getScanNumber() {
		return scanNumber;
	}

	public double getRt() {
		return rt;
	}

	/**
	 * Gets the precursor M/z mass which is m/z = (MW + nH+)/n
	 * 
	 * @return
	 */
	public double getPrecursorMZ() {
		return precursorMZ;
	}

	/**
	 * Gets the precursor single charged mass (M + H+)
	 * 
	 * @return
	 */
	public double getPrecursorMH() {
		final double m_zh = precursorMZ * getZ();
		final double m = m_zh - (getZ() * AssignMass.H);
		final double mh = m + AssignMass.H;
		return mh;
	}

	public TDoubleDoubleMap getPeaks() {
		return peaks;
	}

	public TDoubleList getMasses() {
		return masses;
	}

	public double getMass(int index) {
		return masses.get(index);
	}

	public double getIntensity(int index) {
		return peaks.get(getMass(index));
	}

	public int getZ() {
		return z;
	}

	public void setZ(int z) {
		this.z = z;
	}

	public double getPrecursorIntensity() {
		return precursorIntensity;
	}

	public void setPrecursorIntensity(double precursorIntensity) {
		this.precursorIntensity = precursorIntensity;
	}

}
