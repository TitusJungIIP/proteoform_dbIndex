package edu.scripps.yates.utilities.proteomicsmodel.ms2;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.masses.AssignMass;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;

public class MS2Reader {
	private final static Logger log = Logger.getLogger(MS2Reader.class);
	private final File ms2File;
	private Map<String, TIntList> spectrumIndexesByScanNumberMap;
	private TIntObjectHashMap<String> scanNumberBySpectrumIndexMap;
	private TObjectDoubleHashMap<String> rtByScanNumberMap;
	private TObjectDoubleHashMap<String> precursorIntensityByScanNumberMap;
	private TObjectDoubleHashMap<String> precursorMassByScanNumberMap;
	private final List<Ms2Spectrum> spectra = new ArrayList<Ms2Spectrum>();
	private final Map<String, Ms2Spectrum> spectrumByScanNumber = new THashMap<String, Ms2Spectrum>();
	private boolean grabPrecursorIntensitiesAndMass;
	private final boolean readPeaks;
	private TIntList spectrumIndexes;
	private Map<String, String> hPropertiesMap = new THashMap<String, String>();
	private final static DecimalFormat df = new DecimalFormat("000000");

	/**
	 * Reads ms2 file without reading the peaks
	 * 
	 * @param ms2File
	 * @param readPeaks whether to read peaks or not
	 */
	public MS2Reader(File ms2File) {
		this(ms2File, false);
	}

	/**
	 * 
	 * @param ms2File
	 * @param readPeaks whether to read peaks or not
	 */
	public MS2Reader(File ms2File, boolean readPeaks) {
		this.readPeaks = readPeaks;
		if (readPeaks) {
			setGrabPrecursorIntensities(true);
		}
		this.ms2File = ms2File;
		if (!ms2File.exists()) {
			throw new IllegalArgumentException("File not found at '" + ms2File.getAbsolutePath() + "'");
		}
	}

	public TIntList getSpectrumIndexByScan(String scanNumber) {
		process();
		if (spectrumIndexesByScanNumberMap.containsKey(scanNumber)) {
			return spectrumIndexesByScanNumberMap.get(scanNumber);
		}
		return null;
	}

	private void process() {
		if (spectrumIndexesByScanNumberMap == null || rtByScanNumberMap == null) {
			log.info("Reading MS2 file: " + ms2File.getAbsolutePath());
			spectrumIndexesByScanNumberMap = new THashMap<String, TIntList>();
			scanNumberBySpectrumIndexMap = new TIntObjectHashMap<String>();
			precursorIntensityByScanNumberMap = new TObjectDoubleHashMap<String>();
			precursorMassByScanNumberMap = new TObjectDoubleHashMap<String>();
			rtByScanNumberMap = new TObjectDoubleHashMap<String>();
			BufferedReader br = null;
			int charge = -1;
			int scan1Num = -1;
			int scan2Num = -1;
			int numScan = -1;
			Double rt = null;
			Double intensity = null;
			Double mz = null;
			final Map<String, String> iPropertiesMap = new THashMap<String, String>();
			final List<String> iProperties = new ArrayList<String>();
			Ms2Spectrum spectrum = null;
			try {
				br = new BufferedReader(new FileReader(ms2File));
				String line;
				while ((line = br.readLine()) != null) {
					if (line.startsWith("S")) {

						numScan++;
						final String[] split = line.split("\t");
						final String scan1 = split[1];
						scan1Num = Integer.valueOf(scan1);
						final String scan2 = split[2];
						scan2Num = Integer.valueOf(scan2);
						charge = -1;

						if (grabPrecursorIntensitiesAndMass && split.length > 3) {
							// this mass is m/z
							mz = Double.valueOf(split[3]);

						}
					} else if (line.startsWith("Z")) {

						final String[] split = line.split("\t");
						charge = Integer.valueOf(split[1]);

						final String key = scan1Num + "." + scan2Num + "." + charge;
						if (readPeaks) {

							spectrum = new Ms2Spectrum();
							spectra.add(spectrum);
							spectrum.setZ(charge);
							spectrum.setScanNumber(df.format(scan1Num) + "." + df.format(scan2Num));
							spectrumByScanNumber.put(spectrum.getScanNumber(), spectrum);

						}
						if (!spectrumIndexesByScanNumberMap.containsKey(key)) {
							spectrumIndexesByScanNumberMap.put(key, new TIntArrayList());
						}
						spectrumIndexesByScanNumberMap.get(key).add(numScan);
						scanNumberBySpectrumIndexMap.put(numScan, key);
						if (rt != null) {
							rtByScanNumberMap.put(key, rt);
							if (readPeaks) {
								spectrum.setRt(rt);

							}
						}
						if (intensity != null) {
							precursorIntensityByScanNumberMap.put(key, intensity);
							if (readPeaks) {
								spectrum.setPrecursorIntensity(intensity);
							}

						}
						if (mz != null) {
							if (readPeaks) {
								spectrum.setPrecursorMZ(mz);
							}
							precursorMassByScanNumberMap.put(key, mz);

						} else {
							if (split.length > 2) {
								try {
									// this mass is M+H
									// we follow this equation m/z = (MW + nH+)/n
									// in order to get m/z
									final double mh = Double.valueOf(split[2]);
									final double m = mh - AssignMass.H;
									final double m_nh = m + (charge * AssignMass.H);
									mz = m_nh / charge;
									if (readPeaks) {
										spectrum.setPrecursorMZ(mz);
									}
								} catch (final NumberFormatException e) {
									log.warn(e);
								}
							}
						}
						rt = null;
						intensity = null;
						mz = null;

						iProperties.clear();

						iPropertiesMap.clear();
					} else if (line.startsWith("I")) {

						// separate by TAB here, not space, since some properties have spaces
						final String[] split = line.split("\t");
						if (split[1].equalsIgnoreCase("RetTime") || split[1].equalsIgnoreCase("RTime")) {
							try {
								rt = Double.valueOf(split[2]);

							} catch (final NumberFormatException e) {
								log.warn(e);
								rt = null;
							}
						} else if (grabPrecursorIntensitiesAndMass) {
							if (split[1].equalsIgnoreCase("PrecursorInt")) {
								try {
									intensity = Double.valueOf(split[2]);

								} catch (final NumberFormatException e) {
									log.warn(e);
									intensity = null;
								}
							}
						}
						if (readPeaks && split.length > 1) {
							final String iProperty = split[1];
							String iPropertyValue = null;
							if (split.length > 2) {
								iPropertyValue = split[2];
							}
							iPropertiesMap.put(iProperty, iPropertyValue);
						}
					} else if (line.startsWith("H")) {

						// separate by TAB here, not space, since some properties have spaces
						final String[] split = line.split("\t");
						if (readPeaks && split.length > 1) {
							final String hProperty = split[1];
							String hPropertyValue = null;
							if (split.length > 2) {
								hPropertyValue = split[2];
							}
							hPropertiesMap.put(hProperty, hPropertyValue);
						}
					} else {

						if (readPeaks) {
							// if different than empty string it is a peak
							if (!"".equals(line)) {
								final String[] split = line.split("\\s");
								try {
									final double ms2mass = Double.valueOf(split[0]);
									final double ms2Intensity = Double.valueOf(split[1]);
									spectrum.addPeak(ms2mass, ms2Intensity);
								} catch (final NumberFormatException e) {

								}
							}
						}
					}
				}
				log.info(spectrumIndexesByScanNumberMap.size() + " spectra read in MS2 file. "
						+ ms2File.getAbsolutePath());
				log.info("numScan counter= " + numScan + "\tNum spectra=" + spectra.size());
			} catch (final FileNotFoundException e) {
				e.printStackTrace();
			} catch (final IOException e) {
				e.printStackTrace();
			} finally {
				if (br != null) {
					try {
						br.close();
					} catch (final IOException e) {
						e.printStackTrace();
					}
				}
			}
		}

	}

	/**
	 *
	 * @param scanNumber in the form of "scan.scan.charge"
	 * @return
	 */
	public Double getSpectrumRTByScan(String scanNumber) {
		process();
		if (rtByScanNumberMap.containsKey(scanNumber)) {
			return rtByScanNumberMap.get(scanNumber);
		}
		return null;
	}

	/**
	 *
	 * @param scanNumber in the form of "scan.scan.charge"
	 * @return
	 */
	public Double getPrecursorIntensityByScan(String scanNumber) {
		process();
		if (precursorIntensityByScanNumberMap.containsKey(scanNumber)) {
			return precursorIntensityByScanNumberMap.get(scanNumber);
		}
		return null;
	}

	/**
	 *
	 * @param scanNumber in the form of "scan.scan.charge"
	 * @return
	 */
	public Double getPrecursorMassByScan(String scanNumber) {
		process();
		if (precursorMassByScanNumberMap.containsKey(scanNumber)) {
			return precursorMassByScanNumberMap.get(scanNumber);
		}
		return null;
	}

	public String getFileName() {
		return FilenameUtils.getBaseName(ms2File.getAbsolutePath());
	}

	public Map<String, TIntList> getSpectrumIndexByScanNumberMap() {
		if (spectrumIndexesByScanNumberMap == null) {
			process();
		}
		return spectrumIndexesByScanNumberMap;
	}

	public TObjectDoubleHashMap<String> getRtByScanNumberMap() {
		if (rtByScanNumberMap == null) {
			process();
		}
		return rtByScanNumberMap;
	}

	public boolean isGrabPrecursorIntensities() {
		return grabPrecursorIntensitiesAndMass;
	}

	public void setGrabPrecursorIntensities(boolean grabPrecursorIntensities) {
		this.grabPrecursorIntensitiesAndMass = grabPrecursorIntensities;
	}

	public List<Ms2Spectrum> getSpectra() {
		if (spectra.isEmpty()) {
			if (!readPeaks) {
				throw new IllegalArgumentException("Initialize the reader with readPeaks=true");
			}
			process();
		}
		return spectra;
	}

	public Ms2Spectrum getSpectrumByScanNumber(String scanNumber) {
		if (spectrumByScanNumber == null) {
			process();
		}
		return spectrumByScanNumber.get(scanNumber);
	}

	public TIntList getSpectrumIndexes() {
		if (spectrumIndexes == null) {
			if (scanNumberBySpectrumIndexMap == null) {
				process();
			}
			spectrumIndexes = new TIntArrayList();

			for (final int index : scanNumberBySpectrumIndexMap.keys()) {
				spectrumIndexes.add(index);
			}
			spectrumIndexes.sort();
		}
		return spectrumIndexes;
	}

	public Map<String, String> gethPropertiesMap() {
		return hPropertiesMap;
	}

	public void sethPropertiesMap(Map<String, String> hPropertiesMap) {
		this.hPropertiesMap = hPropertiesMap;
	}
}
