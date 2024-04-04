package edu.scripps.yates.utilities.fasta;

public interface Fasta extends Comparable<Fasta> {

	public void setMPlusH(double mh);

	public double getMPlusH();

	public String getAccessionWithNoVersion();

	public String getSequence();

	public byte[] getSequenceAsBytes();

	public String getOriginalDefline();

	public String getDefline();

	public byte byteAt(int index);

	public int getLength();

	// get accession without version
	public String getAccession();

	public String getSequestLikeAccession();

	public boolean isReversed();

	public boolean isProteoform();

	public static String getSequestLikeAccession(String acc) {
		final String[] arr = acc.split("\t");
		final String[] arr1 = arr[0].split(" ");
		String newacc = arr1[0];
		if (newacc != null && newacc.length() > 40) {
			newacc = newacc.substring(0, 40);
		}
		return newacc;

	}

	public static String getAccessionWithNoVersion(String ac) {
		final int index = ac.indexOf(".");
		if (index == -1) {
			return ac;
		} else {
			return ac.substring(0, index);
		}
	}

	public static String getAccession(String accession) {
		// NCBI, IPI, or others such as UNIT_PROT, SGD, NCI
		// accession = getDefline().substring( getDefline().indexOf('>')+1 );
		// accession = getDefline();

		// There are many corruptted sqt file. Ignore it.
		try {
			if (accession.startsWith("gi") && accession.contains("|")) // NCBI
			{
				String[] arr = accession.split("\\|");

				if (arr.length >= 4
						&& ("gb".equals(arr[2]) || "ref".equals(arr[2]) || "emb".equals(arr[2]) || "dbj".equals(arr[2])
								|| "prf".equals(arr[2]) || "sp".equals(arr[2]))
						|| "tpd".equals(arr[2]) || "tpg".equals(arr[2]) || "tpe".equals(arr[2]))
					accession = arr[3];
				else {
					arr = accession.split(" ");
					accession = arr[1];
				}

				// Accession # should end with digit. If accession # does not
				// end with digit,
				// grap next string (We assume this next one ends with digit.)
				/*
				 * if( pattern.matcher(arr[3]).matches() ) accession = arr[3];
				 * else accession = arr[4].substring(0, arr[4].indexOf(" "));
				 */

			} else if (accession.startsWith("IPI")) // IPI
			{
				final String arr[] = accession.split("\\|");
				final String subArr[] = arr[0].split(":");

				if (subArr.length > 1)
					accession = subArr[1];
				else
					accession = subArr[0];
			} else if (accession.startsWith("Re") || accession.startsWith("contam") || accession.startsWith("Contam")) // Reverse
																														// database
			{
				int space = accession.indexOf(" ");
				int tab = accession.indexOf("\t");

				if (space < 0)
					space = 40;
				if (tab < 0)
					tab = 40;

				final int index = (tab > space) ? space : tab;

				int end;

				if (index <= 0 || index >= 40) // no space
				{
					final int length = accession.length();
					end = (length > 40) ? 40 : length;
				} else
					// cut by the first space
					end = index;

				accession = accession.substring(0, end);
			} else // UNIT_PROT, NCI or SGD

			{
				final int spaceIndex = accession.indexOf(" ");
				int tabIndex;

				if (spaceIndex > 0) {
					tabIndex = accession.indexOf("\t");

					if (tabIndex > 0 && spaceIndex > tabIndex)
						accession = accession.substring(0, tabIndex);
					else
						accession = accession.substring(0, spaceIndex);
				}
			}
		} catch (final Exception e) {
			// System.out.println("No Correct Accession found, but this will be
			// handled by MSP system."
			// + accession + " " + e);

			final int i = accession.indexOf(" ");
			if (i < 0)
				return accession;
			else
				return accession.substring(0, i);

		}

		return accession;
	}

}
