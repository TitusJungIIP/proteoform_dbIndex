package edu.scripps.yates.utilities.index;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel.MapMode;
import java.nio.channels.FileLock;
import java.nio.channels.OverlappingFileLockException;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.bytes.DynByteBuffer;
import edu.scripps.yates.utilities.files.FileUtils;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressNumberFormatter;
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
public class TextFileIndexMultiThreadSafeIO {
	private static final Logger log = Logger.getLogger(TextFileIndexMultiThreadSafeIO.class);
	private static Map<String, FileRecordReservation> fileRecordReservation = new THashMap<String, FileRecordReservation>();
	private int numEntries;
	private final String beginToken;
	private final String endToken;
	private final File fileToIndex;

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
	public TextFileIndexMultiThreadSafeIO(String path, String beginToken, String endToken) throws IOException {
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
	public TextFileIndexMultiThreadSafeIO(File file, String beginToken, String endToken) throws IOException {
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
		log.info("Reading file of " + FileUtils.getDescriptiveSizeFromBytes(totalLength));
		final RandomAccessFile raf = new RandomAccessFile(fileToIndex, "r");
		FileLock lock = raf.getChannel().tryLock(0, Long.MAX_VALUE, true);
		while (lock == null) {
			lock = raf.getChannel().tryLock(0, Long.MAX_VALUE, true);
			try {
				Thread.sleep(1000);
			} catch (final InterruptedException e) {
			}
			log.info("Waiting for reading access to file " + fileToIndex.getAbsolutePath());
		}
		String line;
		try {
			final ProgressCounter counter = new ProgressCounter(totalLength, ProgressPrintingType.PERCENTAGE_STEPS, 0);
			counter.setProgressNumberFormatter(new ProgressNumberFormatter() {

				@Override
				public String format(long number) {
					return FileUtils.getDescriptiveSizeFromBytes(number);
				}
			});
			long offset = 0;
			long init = 0;
			long end = 0;

			StringBuilder sb = new StringBuilder();
			while ((line = raf.readLine()) != null) {
				counter.setProgress(offset);
				final String printIfNecessary = counter.printIfNecessary();
				if (!"".equals(printIfNecessary)) {
					log.info("File index progress: " + printIfNecessary);
				}

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
				if (line.endsWith(endToken)) {
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
	public synchronized Map<String, Pair<Long, Long>> addNewItem(String item, Set<String> keys) throws IOException {

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
		long init = 0l;
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
			// request the record reservation system
			final FileRecordReservation fileRecordReservation = getFileRecordReservation(fileToIndex);
			// book the position in the file (thread safe)
			log.debug("Requesting reservation for writting in position " + fileRecordReservation.getCurrentposition()
					+ " writting " + bytes.length + " bytes in thread " + Thread.currentThread().getId());
			init = fileRecordReservation.reserveRecord(bytes);
			log.debug("Reserved from " + init + " to " + (init + bytes.length) + " in thread "
					+ Thread.currentThread().getId());
			// map the from position to position+bytes.length in the buffer
			buffer = raf.getChannel().map(MapMode.READ_WRITE, init, bytes.length);

			// char[] itemInChars = item.toCharArray();

			buffer.put(bytes);
		} catch (final IOException e) {
			e.printStackTrace();
			log.error("Error from thread " + Thread.currentThread().getId() + ": " + e.getMessage());
		} finally {
			if (lock != null) {
				lock.release();
			}
			raf.close();
			log.debug("Closing file access from thread " + Thread.currentThread().getId());

		}
		final Map<String, Pair<Long, Long>> ret = new THashMap<String, Pair<Long, Long>>();

		// it is important to increase this variable before getKeys() is
		// called
		numEntries++;
		// if everything is fine, store in the map
		final long end = init + bytes.length;
		final Pair<Long, Long> pair = new Pair<Long, Long>(init, end);
		if (keys == null || keys.isEmpty()) {
			keys = getKeys(item);
		}
		for (final String key : keys) {
			ret.put(key, pair);
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
	public synchronized Map<String, Pair<Long, Long>> addNewItems(Map<String, Set<String>> keysByItems)
			throws IOException {
		final Map<String, Pair<Long, Long>> ret = new THashMap<String, Pair<Long, Long>>();

		long firstInit = -1l;
		final DynByteBuffer byteBuffer = new DynByteBuffer();

		for (String item : keysByItems.keySet()) {
			Set<String> keys = keysByItems.get(item);

			// prepare item
			if (!item.startsWith(beginToken)) {
				if (item.contains(beginToken)) {
					item = item.substring(item.indexOf(beginToken));
				} else {
					throw new IllegalArgumentException("The provided item '" + item
							+ "' is not starting with the begin Token '" + beginToken + "'");
				}
			}
			if (!item.endsWith(endToken)) {
				throw new IllegalArgumentException(
						"The provided item '" + item + "' is not ending with the end Token '" + endToken + "'");
			}
			item = "\n" + item;
			final byte[] bytes = item.getBytes();
			byteBuffer.add(bytes);
			// item is ready

			// request the record reservation system
			final FileRecordReservation fileRecordReservation = getFileRecordReservation(fileToIndex);
			// book the position in the file (thread safe)
			log.debug(" Requesting reservation for writting in position " + fileRecordReservation.getCurrentposition()
					+ " writting " + bytes.length + " bytes in thread " + Thread.currentThread().getId());
			final long init = fileRecordReservation.reserveRecord(bytes);
			if (firstInit == -1l) {
				firstInit = init;
			}
			final long end = init + bytes.length;
			log.debug("Reserved from " + init + " to " + end + " in thread " + Thread.currentThread().getId());
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
		}
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

			buffer = raf.getChannel().map(MapMode.READ_WRITE, firstInit, byteBuffer.getSize());
			buffer.put(byteBuffer.getData());
		} catch (final IOException e) {
			e.printStackTrace();
			log.warn("Error from thread " + Thread.currentThread().getId());
		} finally {
			if (lock != null) {
				lock.release();
			}
			raf.close();
			log.debug("Closing file access from thread " + Thread.currentThread().getId());

		}

		return ret;

	}

	private synchronized static FileRecordReservation getFileRecordReservation(File file) {
		if (!fileRecordReservation.containsKey(file.getAbsolutePath())) {
			fileRecordReservation.put(file.getAbsolutePath(), new FileRecordReservation(file));
		}
		return fileRecordReservation.get(file.getAbsolutePath());
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
		FileLock lock = null;
		try {
			lock = raf.getChannel().tryLock(start, end - start, true);
		} catch (final OverlappingFileLockException e) {
			log.info(e);
		}
		while (lock == null) {
			try {
				lock = raf.getChannel().tryLock(start, end - start, true);
			} catch (final OverlappingFileLockException e) {
				log.info(e);
			}
			try {
				Thread.sleep(1000);
			} catch (final InterruptedException e) {
			}
			log.info("Waiting for reading access to file " + fileToIndex.getAbsolutePath());
		}
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
			if (lock != null) {
				lock.release();
			}
			raf.close();

		}
	}
}
