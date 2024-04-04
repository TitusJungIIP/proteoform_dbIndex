package edu.scripps.yates.utilities.env;

import java.util.Map;

import org.apache.log4j.Logger;

public class EnvironmentalVariables {
	private final static Logger log = Logger.getLogger(EnvironmentalVariables.class);

	public static String getEnvironmentalVariableValue(String variableName) {
		Map<String, String> env = System.getenv();

		final String ret = env.get(variableName);
		log.info("Getting environmental variable value of  " + variableName + " = " + ret);
		return ret;
	}
}
