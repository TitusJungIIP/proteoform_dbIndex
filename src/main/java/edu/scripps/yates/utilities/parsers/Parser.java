package edu.scripps.yates.utilities.parsers;

import java.io.File;

public interface Parser {
	/**
	 * This function will try to guess whether the parser is able to read this
	 * file or not using a quick reading of the file
	 * 
	 * @param file
	 * @return
	 */
	public boolean canRead(File file);
}
