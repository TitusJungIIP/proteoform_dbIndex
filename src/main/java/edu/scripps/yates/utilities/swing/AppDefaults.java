package edu.scripps.yates.utilities.swing;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.Properties;

import edu.scripps.yates.utilities.properties.PropertiesUtil;
import gnu.trove.map.hash.THashMap;

class AppDefaults {
	private static AppDefaults instance;

	private static final String PROPERTIES_FILE = "default.properties";

	private Properties properties;
	private final Map<String, String> propertyValues = new THashMap<String, String>();

	private File propertiesFile;
	private final static String comments = "# properties file from HIVPTMAnalyzer";

	public AppDefaults(Collection<String> propertiesNames) {
		try {
			propertiesFile = new File(System.getProperty("user.dir") + File.separator + PROPERTIES_FILE);
			if (propertiesFile.exists()) {
				properties = PropertiesUtil.getProperties(propertiesFile);
				for (final String propertyName : propertiesNames) {
					if (properties.containsKey(propertyName)) {
						this.propertyValues.put(propertyName, properties.getProperty(propertyName));
					}
				}
			} else {
				// create file
				final FileWriter fw = new FileWriter(propertiesFile);
				fw.write(comments);
				fw.close();
				properties = PropertiesUtil.getProperties(propertiesFile);
			}
		} catch (final Exception e) {
			e.printStackTrace();
		}
	}

	public String getPropertyValue(String propertyName) {
		return this.propertyValues.get(propertyName);
	}

	public void setProperty(String propertyName, String propertyValue) {
		propertyValues.put(propertyName, propertyValue);
		properties.put(propertyName, propertyValue);
		try {
			savePropetiesFile();
		} catch (final IOException e) {
			e.printStackTrace();
		}
	}

	private synchronized void savePropetiesFile() throws IOException {
		properties.store(new FileWriter(this.propertiesFile), comments);
	}

}
