package edu.scripps.yates.utilities.matrix.sparse;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipException;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.dates.DatesUtil;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;

public abstract class SparseMatrix<T> {
	private final static Logger log = Logger.getLogger(SparseMatrix.class);
	// Maps with the indexes of the rows and the columns
	private final TObjectIntHashMap<String> rowIndexes = new TObjectIntHashMap<String>();
	private final TObjectIntHashMap<String> colIndexes = new TObjectIntHashMap<String>();

	// data by rows, then by cols
	private final TIntObjectHashMap<TIntObjectHashMap<T>> data = new TIntObjectHashMap<TIntObjectHashMap<T>>();

	// indexes where to store new cols and rows
	private int maxRowIndex = -1;
	private int maxColIndex = -1;

	// number of non-zero values
	private int numValues;
	private File matrixFile;
	private boolean overrideAllowed = true;
	private boolean synchronizedWithFile;
	private File rowsFile;
	private File colsFile;
	private BufferedWriter matrixFileWriter;
	private BufferedWriter columnsWriter;
	private BufferedWriter rowsWriter;
	private boolean compressFiles = true;

	/**
	 * Constructor that creates a matrix from scratch.
	 */
	public SparseMatrix() {
	}

	/**
	 * Constructor in which a matrixFile is loaded (if exists) and a boolean
	 * parameter determines if the file will be synchronized while adding the data.
	 * A flag accessible by <code>isSynchronizedFile()</code> will return whether
	 * the file should be written again or not.
	 * 
	 * @param matrixFile          matrix file
	 * @param synchronizeWithFile
	 * @param overrideAllowed
	 * @param compressFiles
	 * @throws IOException
	 */
	public SparseMatrix(File matrixFile, boolean synchronizeWithFile, boolean allowOverrideData, boolean compressFiles)
			throws IOException {
		this.matrixFile = matrixFile;
		this.synchronizedWithFile = synchronizeWithFile;
		this.overrideAllowed = allowOverrideData;
		loadDataMatrixFromFile();
	}

	/**
	 * A flag telling whether the information in the matrix is synchronized with the
	 * file.
	 * 
	 * @return
	 */
	public boolean isSynchronizedFile() {
		if (matrixFile == null) {
			return false;
		}
		return synchronizedWithFile;
	}

	public int getRowIndex(String rowName) {
		if (rowIndexes.containsKey(rowName)) {
			return rowIndexes.get(rowName);
		}
		return -1;
	}

	public int getColIndex(String colName) {
		if (colIndexes.containsKey(colName)) {
			return colIndexes.get(colName);
		}
		return -1;
	}

	public void addValue(String rowName, String colName, T datapoint) throws IOException {
		int row = getRowIndex(rowName);
		if (row == -1) {
			row = ++maxRowIndex;
			rowIndexes.put(rowName, row);
			if (synchronizedWithFile) {
				writeRow(rowName, true);
			}
		}
		int col = getColIndex(colName);
		if (col == -1) {
			col = ++maxColIndex;
			colIndexes.put(colName, col);
			if (synchronizedWithFile) {
				writeColumn(colName, true);
			}
		}
		addValue(row, col, datapoint, true);

	}

	private void writeColumn(String colName, boolean append) throws IOException {
		final BufferedWriter colsFileWriter = getColsFileWriter(append);
		colsFileWriter.write(colName + "\n");
	}

	private void writeRow(String rowName, boolean append) throws IOException {
		final BufferedWriter rowsFileWriter = getRowsFileWriter(append);
		rowsFileWriter.write(rowName + "\n");
	}

	public void addValue(int rowIndex, int colIndex, T value) throws IOException {
		addValue(rowIndex, colIndex, value, true);
	}

	private void addValue(int rowIndex, int colIndex, T value, boolean writeToFile) throws IOException {
		if (data.containsKey(rowIndex)) {
			final T previousValue = data.get(rowIndex).put(colIndex, value);
			if (previousValue == null) {
				numValues++;
				if (writeToFile && synchronizedWithFile) {
					writeValue(rowIndex, colIndex, value, true);
				}
			} else if (!overrideAllowed) {
				data.get(rowIndex).put(colIndex, previousValue);
			} else {
				if (writeToFile && synchronizedWithFile) {
					writeValue(rowIndex, colIndex, value, true);
				}
			}
		} else {
			final TIntObjectHashMap<T> valuesPerColumn = new TIntObjectHashMap<T>();
			valuesPerColumn.put(colIndex, value);
			data.put(rowIndex, valuesPerColumn);
			numValues++;
			if (writeToFile && synchronizedWithFile) {
				writeValue(rowIndex, colIndex, value, true);
			}
		}
		if (numValues % 1000000 == 0) {
			System.out.println(numValues + " values in matrix with " + maxRowIndex + " rows and " + maxColIndex
					+ " cols (" + percentageFormat.format(getSparsity()) + " sparse values)");
		}
	}

	private double getSparsity() {
		final double denseMatrixSize = 1.0 * maxColIndex * maxRowIndex;
		final double numZeros = denseMatrixSize - numValues;
		return numZeros / denseMatrixSize;
	}

	private final static DecimalFormat percentageFormat = new DecimalFormat("#.#%");

	public static void main(String[] args) {
		final File file = new File(
				"C:\\Users\\salvador\\Desktop\\casimir\\SARS_Cov2\\singlecell_human\\singlecells_human_ALL_fullTable_test.txt.gz");
		SparseMatrix<Short> matr = null;
		try {
			matr = new SparseMatrix<Short>(file, true, true, true) {

				@Override
				public Short readFromString(String stringRepresentationOfTheData) {
					return Short.valueOf(stringRepresentationOfTheData);
				}
			};

			final int rowIndex = 6101;
			final int colIndex = 13921;
			final Short value = matr.getValue(rowIndex, colIndex);
			final String cell = matr.getReverseRowColIndexes(matr.rowIndexes).get(rowIndex);
			final String gene = matr.getReverseRowColIndexes(matr.colIndexes).get(colIndex);
			System.out.println(gene + "\t" + cell + "\t" + value);

			final Map<String, Short> cellExpressions = matr.getRow(cell);
			for (final String gene2 : cellExpressions.keySet()) {
				if (Short.compare(matr.getValue(cell, gene2), cellExpressions.get(gene2)) == 0) {
					System.out.println(cell + "\t" + gene2 + "\t" + cellExpressions.get(gene2));
				} else {
					System.err.println("ERROR");
				}
			}
		} catch (final IOException e) {
			e.printStackTrace();
		}
		final int colIndex = matr.getColIndex("asdf");
		System.out.println(colIndex);
	}

	public T deleteValue(String rowName, String colName) {
		return deleteValue(getRowIndex(rowName), getColIndex(colName));
	}

	public T deleteValue(int rowIndex, int colIndex) {
		if (data.containsKey(rowIndex)) {
			if (data.get(rowIndex).containsKey(colIndex)) {
				final T ret = data.get(rowIndex).remove(colIndex);
				if (ret != null) {
					numValues--;
				}
				if (data.get(rowIndex).isEmpty()) {
					data.remove(rowIndex);
				}
				return ret;
			}
		}
		return null;
	}

	private Map<String, T> getRow(String rowName) {
		final TIntObjectHashMap<String> reverseColIndexes = getReverseRowColIndexes(colIndexes);
		final Map<String, T> ret = new THashMap<String, T>();
		final int rowIndex = getRowIndex(rowName);
		if (rowIndex >= 0) {
			final TIntObjectHashMap<T> row = data.get(rowIndex);
			for (final int colIndex : row.keys()) {
				final String colName = reverseColIndexes.get(colIndex);
				final T value = row.get(colIndex);
				ret.put(colName, value);
			}
		}
		return ret;
	}

	public T getValue(String rowName, String colName) {
		return getValue(getRowIndex(rowName), getColIndex(colName));

	}

	public T getValue(int rowIndex, int colIndex) {
		if (rowIndex < 0 || colIndex < 0) {
			return null;
		}
		if (data.containsKey(rowIndex)) {
			return data.get(rowIndex).get(colIndex);
		}
		return null;
	}

	private BufferedWriter getMatrixFileWriter(boolean append) throws IOException {
		if (matrixFileWriter == null) {
			if (matrixFile != null) {
				// third, the data in the matrix file
				if (this.compressFiles) {
					final GZIPOutputStream zip = new GZIPOutputStream(new FileOutputStream(this.matrixFile, append));
					matrixFileWriter = new BufferedWriter(new OutputStreamWriter(zip, "UTF-8"));
				} else {
					matrixFileWriter = new BufferedWriter(
							new OutputStreamWriter(new FileOutputStream(matrixFile, append), "UTF-8"));
				}
			}
		}
		return matrixFileWriter;
	}

	private File getRowsFile() {

		if (rowsFile == null && matrixFile != null) {
			String baseName = FilenameUtils.getBaseName(matrixFile.getAbsolutePath());
			String baseName2 = "";
			while (true) {
				baseName2 = FilenameUtils.getBaseName(baseName);
				if (!baseName.equals(baseName2)) {
					baseName = baseName2;

				} else {
					break;
				}
			}
			rowsFile = new File(matrixFile.getParent() + File.separator + baseName + "_rows.txt.gz");
		}

		return rowsFile;
	}

	private BufferedWriter getRowsFileWriter(boolean append) throws IOException {
		if (matrixFile != null && rowsWriter == null) {
			getRowsFile();
			if (this.compressFiles) {
				final GZIPOutputStream zip = new GZIPOutputStream(new FileOutputStream(this.rowsFile, append));
				rowsWriter = new BufferedWriter(new OutputStreamWriter(zip, "UTF-8"));
			} else {
				rowsWriter = new BufferedWriter(
						new OutputStreamWriter(new FileOutputStream(rowsFile, append), "UTF-8"));
			}
		}
		return rowsWriter;
	}

	private File getColsFile() {
		if (colsFile == null && matrixFile != null) {
			String baseName = FilenameUtils.getBaseName(matrixFile.getAbsolutePath());
			String baseName2 = "";
			while (true) {
				baseName2 = FilenameUtils.getBaseName(baseName);
				if (!baseName.equals(baseName2)) {
					baseName = baseName2;

				} else {
					break;
				}
			}
			colsFile = new File(matrixFile.getParent() + File.separator + baseName + "_cols.txt.gz");
		}
		return colsFile;
	}

	private BufferedWriter getColsFileWriter(boolean append) throws IOException {
		if (matrixFile != null && columnsWriter == null) {
			getColsFile();
			if (this.compressFiles) {
				final GZIPOutputStream zip = new GZIPOutputStream(new FileOutputStream(this.colsFile, append));
				columnsWriter = new BufferedWriter(new OutputStreamWriter(zip, "UTF-8"));
			} else {
				columnsWriter = new BufferedWriter(
						new OutputStreamWriter(new FileOutputStream(colsFile, append), "UTF-8"));
			}
		}
		return columnsWriter;
	}

	public void printToFile(File outputFile, boolean gzip, Boolean append) throws IOException {
		if (append == null) {
			append = false;
			if (matrixFile != null) {
				append = true;
			}
		}
		if (matrixFile.equals(outputFile) && isSynchronizedFile()) {
			return;
		}
		this.matrixFile = outputFile;
		this.rowsFile = new File(
				outputFile.getParent() + File.separator + FilenameUtils.getBaseName(outputFile.getAbsolutePath())
						+ "_rows." + FilenameUtils.getExtension(outputFile.getAbsolutePath()));
		this.colsFile = new File(
				outputFile.getParent() + File.separator + FilenameUtils.getBaseName(outputFile.getAbsolutePath())
						+ "_rows." + FilenameUtils.getExtension(outputFile.getAbsolutePath()));

	}

	private void writeValue(int rowIndex, int colIndex, T value, boolean append) throws IOException {

		final BufferedWriter writer = getMatrixFileWriter(append);

		writer.write(rowIndex + "\t" + colIndex + "\t" + value + "\n");

	}

	public void writeDataToFile(File matrixFile) throws IOException {
		this.matrixFile = matrixFile;
		try {

			// third we print all the triplets: row, col, value, one per line

			final TIntList rows = new TIntArrayList(data.keys());
			rows.sort();
			final ProgressCounter counter = new ProgressCounter(numValues, ProgressPrintingType.PERCENTAGE_STEPS, 0);
			counter.setShowRemainingTime(true);
			counter.setSuffix("writting matrix");
			for (final int row : rows.toArray()) {
				for (final int col : data.get(row).keys()) {
					final T value = data.get(row).get(col);
					writeValue(row, col, value, false);
					counter.increment();
					final String printIfNecessary = counter.printIfNecessary();
					if (!"".equals(printIfNecessary)) {
						System.out.println(printIfNecessary);
					}
				}
			}

		} finally {
			if (matrixFileWriter != null) {
				matrixFileWriter.close();
				this.synchronizedWithFile = true;
				log.info("Matrix with " + numValues + " values, " + maxRowIndex + " rows and " + maxColIndex
						+ " columns and " + percentageFormat.format(getSparsity()) + " sparse values");
				log.info("File writen at: " + matrixFile.getAbsolutePath() + " ("
						+ edu.scripps.yates.utilities.files.FileUtils.getDescriptiveSizeFromBytes(matrixFile.length())
						+ ")");
			}
		}
	}

	private void loadDataMatrixFromFile() throws IOException {
		if (this.matrixFile == null || !matrixFile.exists()) {
			return;
		}
		loadColsFromFile();
		loadRowsFromFile();
		final long t0 = System.currentTimeMillis();
		log.info("Loading matrix from file " + matrixFile.getAbsolutePath());
		BufferedReader br = null;
		try {
			InputStreamReader isr = null;
			try {

				final GZIPInputStream gzipIS = new GZIPInputStream(new FileInputStream(matrixFile));
				isr = new InputStreamReader(gzipIS);
			} catch (final ZipException e) {
				// try with regular reader
				isr = new InputStreamReader(new FileInputStream(matrixFile));
			}
			br = new BufferedReader(isr);
			String line = null;
			while ((line = br.readLine()) != null) {
				line = line.trim();
				if ("".equals(line) || line.startsWith("#")) {
					continue;
				}
				final String[] split = line.split("\t");

				// then the data starts in triplets
				final int row = Integer.valueOf(split[0]);
				final int col = Integer.valueOf(split[1]);
				final T datapoint = readFromString(split[2]);
				addValue(row, col, datapoint, false);
			}
		} finally {
			if (br != null) {
				br.close();
				this.synchronizedWithFile = true;
				final long t1 = System.currentTimeMillis() - t0;
				log.info("Matrix read with " + numValues + " values, " + maxRowIndex + " rows and " + maxColIndex
						+ " columns and " + percentageFormat.format(getSparsity()) + " sparse values in "
						+ DatesUtil.getDescriptiveTimeFromMillisecs(t1));

			}
		}
	}

	private void loadRowsFromFile() throws IOException {
		if (this.matrixFile == null || !matrixFile.exists()) {
			return;
		}
		final File file = getRowsFile();
		log.info("Loading rows from file " + file.getAbsolutePath());
		BufferedReader br = null;
		try {
			InputStreamReader isr = null;
			try {
				final GZIPInputStream gzipIS = new GZIPInputStream(new FileInputStream(file));
				isr = new InputStreamReader(gzipIS);
			} catch (final ZipException e) {
				// try with regular reader
				isr = new InputStreamReader(new FileInputStream(file));
			}
			br = new BufferedReader(isr);
			String line = null;
			int row = 0;
			while ((line = br.readLine()) != null) {
				line = line.trim();
				if ("".equals(line) || line.startsWith("#")) {
					continue;
				}
				final String rowName = line;
				rowIndexes.put(rowName, row);
				maxColIndex = row;
				row++;
			}
		} finally {
			if (br != null) {
				br.close();
				log.info(maxColIndex + " rows read");
			}
		}
	}

	private void loadColsFromFile() throws IOException {
		if (this.matrixFile == null || !matrixFile.exists()) {
			return;
		}
		final File file = getColsFile();
		log.info("Loading columns from file " + file.getAbsolutePath());
		BufferedReader br = null;
		try {
			InputStreamReader isr = null;
			try {
				final GZIPInputStream gzipIS = new GZIPInputStream(new FileInputStream(file));
				isr = new InputStreamReader(gzipIS);
			} catch (final ZipException e) {
				// try with regular reader
				isr = new InputStreamReader(new FileInputStream(file));
			}
			br = new BufferedReader(isr);
			String line = null;
			int col = 0;
			while ((line = br.readLine()) != null) {
				line = line.trim();
				if ("".equals(line) || line.startsWith("#")) {
					continue;
				}
				final String colName = line;
				colIndexes.put(colName, col);
				maxRowIndex = col;

				col++;
			}
		} finally {
			if (br != null) {
				br.close();
				log.info(maxColIndex + " rows read");
			}
		}
	}

	private TIntObjectHashMap<String> getReverseRowColIndexes(TObjectIntHashMap<String> map) {
		final TIntObjectHashMap<String> reverseIndexes = new TIntObjectHashMap<String>();
		final Object[] rowNames = map.keys();
		for (final Object rowName : rowNames) {
			final int index = map.get(rowName);
			reverseIndexes.put(index, (String) rowName);
		}
		return reverseIndexes;
	}

	/**
	 * This method needs to be defined so that when reading the file the object is
	 * read properly
	 * 
	 * @param stringRepresentationOfTheData
	 * @return
	 */
	public abstract T readFromString(String stringRepresentationOfTheData);

	public File getMatrixFile() {
		return this.matrixFile;
	}

	public boolean isCompressMatrixFile() {
		return compressFiles;
	}

	public void setCompressMatrixFile(boolean compressMatrixFile) {
		this.compressFiles = compressMatrixFile;
	}

	public void finishSynchronization() throws IOException {
		if (this.matrixFileWriter != null) {
			matrixFileWriter.close();
		}
		if (this.rowsWriter != null) {
			rowsWriter.close();
		}
		if (this.columnsWriter != null) {
			columnsWriter.close();
		}

	}
}
