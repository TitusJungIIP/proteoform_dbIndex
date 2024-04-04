package edu.scripps.yates.utilities.parsers.idparser;

import java.util.Map;

public interface CommandLineParameters {

	Map<String, String> getParametersMap();

	String getParameterValue(String parameterName);

}
