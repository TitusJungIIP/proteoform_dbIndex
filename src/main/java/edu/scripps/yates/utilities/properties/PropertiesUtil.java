package edu.scripps.yates.utilities.properties;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class PropertiesUtil {
	private static final String PROPERTIES_FILE = "proteindb.properties";

	public static Properties getProperties() throws Exception {
		final ClassLoader cl = PropertiesUtil.class.getClassLoader();
		InputStream is;

		is = cl.getResourceAsStream(PROPERTIES_FILE);
		if (is == null)
			throw new Exception(PROPERTIES_FILE + " file not found");

		return getProperties(is);
	}

	public static Properties getProperties(String classPathFileName) {
		final ClassLoader cl = PropertiesUtil.class.getClassLoader();
		InputStream is;

		is = cl.getResourceAsStream(classPathFileName);
		if (is == null)
			throw new IllegalArgumentException(classPathFileName + " file not found");

		final Properties prop = new Properties();
		try {
			prop.load(is);
		} catch (final IOException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(e);
		}
		return prop;
	}

	public static Properties getProperties(InputStream is) throws Exception {

		if (is == null) {
			throw new Exception("Input stream is null");
		}
		final Properties prop = new Properties();
		try {
			prop.load(is);
			is.close();
		} catch (final IOException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(e);
		}
		return prop;
	}

	public static Properties getProperties(File propertiesFile) throws Exception {
		InputStream is;

		is = new FileInputStream(propertiesFile);

		return getProperties(is);
	}

	public static String getPropertyValue(String propertyName) {
		try {
			return getProperties().getProperty(propertyName);
		} catch (final Exception e) {
			e.printStackTrace();
		}
		return null;
	}
}
