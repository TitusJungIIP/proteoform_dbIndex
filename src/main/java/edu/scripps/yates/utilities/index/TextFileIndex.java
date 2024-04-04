package edu.scripps.yates.utilities.index;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.nio.channels.FileLock;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.files.FileUtils;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.map.hash.THashMap;

public class TextFileIndex implements FileIndex<String> {
	private final Logger log = Logger.getLogger(TextFileIndex.class);
	private final File fileToIndex;
	private final File indexFile;
	private final Map<String, Pair<Long, Long>> indexMap = new THashMap<String, Pair<Long, Long>>();
	private final static String INDEX_EXT = ".idx";
	private static final String TAB = "\t";
	private static final String NEWLINE = "\n";
	private final TextFileIndexIO textFileIndexIO;

	private enum Status {
		READY, NOT_READY
	};

	private Status status = Status.NOT_READY;

	public TextFileIndex(File file, TextFileIndexIO textFileIndexIO) throws IOException {
		fileToIndex = file;
		// create the index file .idx
		indexFile = new File(getIndexPathName(file));
		this.textFileIndexIO = textFileIndexIO;
	}

	public TextFileIndex(File file, String beginToken, String endToken) throws IOException {
		fileToIndex = file;
		// create the index file
		indexFile = new File(getIndexPathName(file));
		textFileIndexIO = new TextFileIndexIO(fileToIndex, beginToken, endToken);
	}

	private String getIndexPathName(File file) {
		final String pathName = file.getParent() + File.separator + FilenameUtils.getBaseName(file.getAbsolutePath())
				+ INDEX_EXT;
		return pathName;
	}

	public TextFileIndex(String path, String beginToken, String endToken) throws IOException {
		this(new File(path), beginToken, endToken);
	}

	public void indexFile() throws IOException {
		log.info("Indexing file " + FilenameUtils.getName(fileToIndex.getAbsolutePath()) + "...");
		// read the index, getting the positions of the items
		final Map<String, Pair<Long, Long>> indexMap = textFileIndexIO.getIndexMap();
		// add to the map
		this.indexMap.putAll(indexMap);
		// write the index file without appending
		writePositionsInIndex(indexMap, false);

	}

	private void writePositionsInIndex(Map<String, Pair<Long, Long>> itemPositions, boolean appendOnIndexFile)
			throws IOException {
		final FileOutputStream fos = new FileOutputStream(indexFile, appendOnIndexFile);
		FileLock lock = fos.getChannel().tryLock();
		while (lock == null) {
			lock = fos.getChannel().tryLock();
			try {
				Thread.sleep(1000);
			} catch (final InterruptedException e) {
			}
			log.info("Waiting for writting access to file " + indexFile.getAbsolutePath());
		}
		final BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));
		try {
			for (final String key : itemPositions.keySet()) {
				final Pair<Long, Long> pair = itemPositions.get(key);
				bw.write(key + TAB + pair.getFirstelement() + TAB + pair.getSecondElement() + NEWLINE);
			}
		} finally {
			if (lock != null) {
				lock.release();
			}
			bw.close();
			status = Status.READY;
			log.info("Indexing done. File of index: " + FileUtils.getDescriptiveSizeFromBytes(indexFile.length()));

		}

	}

	@Override
	public String getItem(String key) {
		try {
			// load index file
			final Map<String, Pair<Long, Long>> indexMap = loadIndexFile();
			// look for the provided key
			if (indexMap.containsKey(key)) {
				final Pair<Long, Long> pair = indexMap.get(key);
				final String item = textFileIndexIO.getItem(pair.getFirstelement(), pair.getSecondElement());
				return item;
			}
		} catch (final IOException e) {
			e.printStackTrace();
			log.error(e.getMessage());
		}
		log.info("Item with key '" + key + "' is not in the index");
		return null;
	}

	private Map<String, Pair<Long, Long>> loadIndexFile() throws IOException {
		// if not ready, means that the index file has to be updated
		if (status == Status.NOT_READY) {
			indexMap.clear();
			if (indexFile == null || !indexFile.exists() || indexFile.length() <= 0) {
				indexFile();
			}
		}
		// if index Map is empty, read the index file
		if (indexMap.isEmpty()) {
			final FileInputStream fis = new FileInputStream(indexFile);
			FileLock lock = fis.getChannel().tryLock(0, Long.MAX_VALUE, true);
			while (lock == null) {
				lock = fis.getChannel().tryLock(0, Long.MAX_VALUE, true);
				try {
					Thread.sleep(1000);
				} catch (final InterruptedException e) {
				}
				log.info("Waiting for reading access to file " + indexFile.getAbsolutePath());
			}
			final BufferedReader fr = new BufferedReader(new InputStreamReader(fis));
			try {
				String line;
				while ((line = fr.readLine()) != null) {
					final String[] split = line.split(TAB);
					final String key = split[0];
					final long start = new Long(split[1]);
					final long end = new Long(split[2]);
					final Pair<Long, Long> pair = new Pair<Long, Long>(start, end);
					indexMap.put(key, pair);
				}
			} finally {
				if (lock != null) {
					lock.release();
				}
				fr.close();

			}
		}
		return indexMap;
	}

	/**
	 * Adds an item to the index. It will be written in the indexed file, and the
	 * index file will be updated.
	 * 
	 * @param item
	 * @return
	 * @throws IOException
	 */
	@Override
	public Map<String, Pair<Long, Long>> addItem(String item, Set<String> keys) {
		// load index file
		try {
			loadIndexFile();

			// add into the file to index

			final Map<String, Pair<Long, Long>> itemPositions = textFileIndexIO.addNewItem(item, keys);

			// add to the map
			indexMap.putAll(itemPositions);

			// write the index file appending
			writePositionsInIndex(itemPositions, true);

			// return the positions
			return itemPositions;
		} catch (final IOException e) {
			e.printStackTrace();
			log.error(e.getMessage());
		}
		return null;
	}

	@Override
	public boolean isEmpty() {
		if (status == Status.READY && indexMap.isEmpty()) {
			return true;
		}
		return false;
	}

	@Override
	public List<String> getItems(Collection<String> keys) {
		final List<String> ret = new ArrayList<String>();
		for (final String key : keys) {
			final String item = getItem(key);
			if (item != null) {
				ret.add(item);
			}
		}
		return ret;
	}

	public List<String> getAllItemKeys() throws IOException {

		final List<String> ret = new ArrayList<String>();
		ret.addAll(loadIndexFile().keySet());
		return ret;
	}

}
