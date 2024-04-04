package edu.scripps.yates.utilities.fasta;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

public class GeneBankToFasta {

	public static void main(String[] args) {
		try {
			final File folder = new File(args[0]);
			final FileWriter fw = new FileWriter(new File(folder.getAbsolutePath() + File.separator + "merged.fasta"));
			boolean inTranslation = false;
			String locus = null;
			String description = null;
			int total = 0;
			String proteinID = null;
			for (final File file : folder.listFiles()) {
				int numEntries = 0;
				final BufferedReader br = Files.newBufferedReader(Paths.get(file.toURI()));
				String line = null;
				while ((line = br.readLine()) != null) {
					line = line.trim();
					if (line.startsWith("BASE COUNT")) {
						inTranslation = false;
						if (proteinID != null) {
							fw.write("\n");
						}
						proteinID = null;
					}
					if (inTranslation) {
						fw.write(line.replace("\"", ""));
					}
					if (line.startsWith("LOCUS")) {
						final String[] split = line.split("\\s+");
						locus = split[1];

					}
					if (line.startsWith("/protein_id=\"")) {
						proteinID = line.split("\"")[1].replace("\"", "");
						fw.write(">" + proteinID + " locus=" + locus + " " + description + "\n");
						numEntries++;
						total++;
					}
					if (line.startsWith("/product=\"")) {
						description = line.split("\"")[1];

					}
					if (line.startsWith("/translation=\"")) {
						inTranslation = true;
						final String translation = line.split("\"")[1];
						fw.write(translation);
					}
				}
				System.out.println(numEntries + "/" + total + " proteins in file " + file.getAbsolutePath());
				br.close();
			}
			fw.close();
		} catch (final IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
}
