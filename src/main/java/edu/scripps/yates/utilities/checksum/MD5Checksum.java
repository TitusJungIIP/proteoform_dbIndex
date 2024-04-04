package edu.scripps.yates.utilities.checksum;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;

public class MD5Checksum {

	public static byte[] createChecksumFromFile(String filename) throws Exception {
		InputStream fis = new FileInputStream(filename);

		return createChecksum(fis);
	}

	public static byte[] createChecksumFromFile(File file) throws Exception {
		InputStream fis = new FileInputStream(file);

		return createChecksum(fis);
	}

	public static byte[] createChecksumFromString(String string) throws Exception {
		InputStream is = new ByteArrayInputStream(string.getBytes(StandardCharsets.UTF_8));
		return createChecksum(is);
	}

	public static byte[] createChecksum(InputStream fis) throws Exception {

		byte[] buffer = new byte[1024];
		MessageDigest complete = MessageDigest.getInstance("MD5");
		int numRead;

		do {
			numRead = fis.read(buffer);
			if (numRead > 0) {
				complete.update(buffer, 0, numRead);
			}
		} while (numRead != -1);

		fis.close();
		return complete.digest();
	}

	// see this How-to for a faster way to convert
	// a byte array to a HEX string
	public static String getMD5ChecksumFromFileName(String filename) throws Exception {
		byte[] b = createChecksumFromFile(filename);
		return getStringFromBytes(b);
	}

	public static String getMD5ChecksumFromString(String string) throws Exception {
		byte[] b = createChecksumFromString(string);
		return getStringFromBytes(b);
	}

	public static String getMD5ChecksumFromFileName(File file) throws Exception {
		byte[] b = createChecksumFromFile(file);
		return getStringFromBytes(b);

	}

	public static String getStringFromBytes(byte[] b) {
		StringBuilder result = new StringBuilder();

		for (int i = 0; i < b.length; i++) {
			result.append(Integer.toString((b[i] & 0xff) + 0x100, 16).substring(1));
		}
		return result.toString();
	}

}
