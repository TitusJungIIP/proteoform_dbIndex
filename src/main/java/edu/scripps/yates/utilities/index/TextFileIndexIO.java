package edu.scripps.yates.utilities.index;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel.MapMode;
import java.nio.channels.FileLock;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.files.FileUtils;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

/**
 * Indexer of text files. This class will read a text file, usually a large one,
 * and it will get the start and end offsets for the items found in the text
 * files. The items are determined by a start and an end string tokens. The
 * offsets are then retrieved in a {@link Map} while the key of the offsets is,
 * by default, the number of the item.<br>
 * In order to index by other keys, override the getKeys() method.<br>
 * This class doesn't keep the indexes in memory.
 * 
 * @author Salva
 * 
 */
public class TextFileIndexIO {
	private static final Logger log = Logger.getLogger(TextFileIndexIO.class);
	protected int numEntries;
	protected final String beginToken;
	protected final String endToken;
	protected final File fileToIndex;

	/**
	 * Constructor of the indexer of text files. It takes, the path to the file, and
	 * two strings specifying the start and end string tokens from which the entries
	 * in the index are going to be taken.
	 * 
	 * @param path       path to the file to index
	 * @param beginToken
	 * @param endToken   the string from which an indexed entry of the text file is
	 *                   starting
	 * @throws IOException the string from which an indexed entry of the text file
	 *                     is ending
	 */
	public TextFileIndexIO(String path, String beginToken, String endToken) throws IOException {
		this(new File(path), beginToken, endToken);

	}

	/**
	 * Constructor of the indexer of text files. It takes, the path to the file, and
	 * two strings specifying the start and end string tokens from which the entries
	 * in the index are going to be taken.
	 * 
	 * @param path       path to the file to index
	 * @param beginToken
	 * @param endToken   the string from which an indexed entry of the text file is
	 *                   starting
	 * @throws IOException the string from which an indexed entry of the text file
	 *                     is ending
	 */
	public TextFileIndexIO(File file, String beginToken, String endToken) throws IOException {
		this.beginToken = beginToken;
		this.endToken = endToken;
		numEntries = 0;
		fileToIndex = file;
	}

	/**
	 * Gets the keys to use in the index for each entry.<br>
	 * In this case, the index will not have nothing to do with the entry string,
	 * and it will be the number of entry, counted from 1.<br>
	 * Override this function in order to index the entry by another custom keys.
	 * 
	 * @param string
	 * @return
	 */
	protected Set<String> getKeys(String string) {
		final Set<String> set = new THashSet<String>();
		set.add(String.valueOf(numEntries));
		return set;
	}

	/**
	 * @return the indexMap
	 * @throws IOException
	 */
	public Map<String, Pair<Long, Long>> getIndexMap() throws IOException {
		final Map<String, Pair<Long, Long>> ret = new THashMap<String, Pair<Long, Long>>();
		if (!fileToIndex.exists())
			return ret;
		final long totalLength = fileToIndex.length();
		final RandomAccessFile raf = new RandomAccessFile(fileToIndex, "r");
		FileLock lock = raf.getChannel().tryLock(0, Long.MAX_VALUE, true);
		while (lock == null) {
			lock = raf.getChannel().tryLock(0, Long.MAX_VALUE, true);
			try {
				Thread.sleep(1000);
			} catch (final InterruptedException e) {
			}
			log.info("Waiting for writting access to file " + fileToIndex.getAbsolutePath());
		}
		String line;
		try {
			long offset = 0;
			long init = 0;
			long end = 0;
			final ProgressCounter counter = new ProgressCounter(totalLength, ProgressPrintingType.PERCENTAGE_STEPS, 0);
			StringBuilder sb = new StringBuilder();
			while ((line = raf.readLine()) != null) {

				line = line.trim();
				sb.append(line + "\n");
				if (line.startsWith(beginToken)) {
					sb = new StringBuilder();
					sb.append(line + "\n");
					numEntries++;
					// log.debug(offset + " " + line);
					init = offset;
				}
				offset = raf.getFilePointer();
				counter.setProgress(offset);
				final String printIfNecessary = counter.printIfNecessary();
				if (!"".equals(printIfNecessary)) {
					log.debug(FileUtils.getDescriptiveSizeFromBytes(offset) + "/"
							+ FileUtils.getDescriptiveSizeFromBytes(totalLength) + " (" + printIfNecessary
							+ ") readed...");
				}
				if ((line.endsWith(endToken) && !"".equals(endToken)) || ("".equals(endToken) && "".equals(line))) {
					end = offset;
					final Pair<Long, Long> pair = new Pair<Long, Long>(init, end);
					final Set<String> keys = getKeys(sb.toString());
					for (final String key : keys) {
						ret.put(key, pair);
					}
				}
			}

		} finally {
			if (lock != null) {
				lock.release();
			}
			raf.close();
		}
		return ret;
	}

	/**
	 * Adds a new item in the file to index and return the keys and the position in
	 * which the item has been written
	 * 
	 * @param item
	 * @return a map containing the keys and position in where the item was stored
	 * @throws IOException
	 */
	public Map<String, Pair<Long, Long>> addNewItem(String item, Set<String> keys) throws IOException {

		if (!item.startsWith(beginToken)) {
			if (item.contains(beginToken)) {
				item = item.substring(item.indexOf(beginToken));
			} else {
				throw new IllegalArgumentException(
						"The provided item '" + item + "' is not starting with the begin Token '" + beginToken + "'");
			}
		}
		if (!item.endsWith(endToken))
			throw new IllegalArgumentException(
					"The provided item '" + item + "' is not ending with the end Token '" + endToken + "'");

		item = "\n" + item;
		final byte[] bytes = item.getBytes();

		MappedByteBuffer buffer = null;
		final RandomAccessFile raf = new RandomAccessFile(fileToIndex, "rws");
		FileLock lock = raf.getChannel().tryLock();
		while (lock == null) {
			lock = raf.getChannel().tryLock();
			try {
				Thread.sleep(1000);
			} catch (final InterruptedException e) {
			}
			log.info("Waiting for writting access to file " + fileToIndex.getAbsolutePath());
		}
		try {

			buffer = raf.getChannel().map(MapMode.READ_WRITE, raf.length(), bytes.length);

			final Map<String, Pair<Long, Long>> ret = new THashMap<String, Pair<Long, Long>>();
			// char[] itemInChars = item.toCharArray();

			// go to the end
			// raf.seek(raf.length());
			// raf.writeBytes("\n");
			// long init = raf.length();
			// long end = raf.length() + itemInChars.length;
			// raf.writeBytes(item);

			// go to the end
			final long init = raf.length();
			final long end = raf.length() + bytes.length;
			buffer.put(bytes);

			// it is important to increase this variable before getKeys() is
			// called
			numEntries++;
			// if everything is fine, store in the map
			final Pair<Long, Long> pair = new Pair<Long, Long>(init, end);
			if (keys == null || keys.isEmpty()) {
				keys = getKeys(item);
			}
			for (final String key : keys) {
				ret.put(key, pair);
			}
			return ret;
		} finally {
			if (lock != null) {
				lock.release();
			}
			raf.close();

		}
	}

	public File getFileToIndex() {
		return fileToIndex;
	}

	public String getItem(Pair<Long, Long> positions) throws IOException {
		return getItem(positions.getFirstelement(), positions.getSecondElement());
	}

	/**
	 * 
	 * @param firstelement
	 * @param secondElement
	 * @throws IOException
	 */
	public String getItem(Long start, Long end) throws IOException {

		final RandomAccessFile raf = new RandomAccessFile(fileToIndex, "r");
		try {
			// go to the start
			raf.seek(start);
			final int lenthToRead = new Long(end - start).intValue();
			// array to store the readed item
			final byte[] bytesToRead = new byte[lenthToRead];
			raf.read(bytesToRead);
			final String readed = new String(bytesToRead);
			return readed;
		} finally {
			raf.close();
		}
	}
}
