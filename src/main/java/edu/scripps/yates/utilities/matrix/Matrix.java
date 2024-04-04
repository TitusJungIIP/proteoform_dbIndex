package edu.scripps.yates.utilities.matrix;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.files.FileUtils;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public abstract class Matrix<K, T> {

	private final Map<K, Map<K, T>> matrix = new THashMap<K, Map<K, T>>();
	private final List<K> rowKeys = new ArrayList<K>();
	private final List<K> colKeys = new ArrayList<K>();
	private final Set<K> colKeysSet = new THashSet<K>();
	private Boolean isSymetrix;
	private Object emptyValue;
	private static final Logger log = Logger.getLogger(Matrix.class);
	private String newlineString = "\n";

	/**
	 * 
	 * @param emptyValue value that will be printed on the file when there is no
	 *                   value in the cell
	 */
	public Matrix(Object emptyValue) {
		this.emptyValue = emptyValue;
	}

	public void setEmptyValue(Object emptyValue) {
		this.emptyValue = emptyValue;
	}

	public Object getEmptyValue() {
		return this.emptyValue;
	}

	public T getValue(int rowIndex, int columnIndex) {
		try {
			return matrix.get(rowKeys.get(rowIndex)).get(colKeys.get(columnIndex));
		} catch (final Exception e) {
			return null;
		}
	};

	public T getValue(K rowKey, K columnKey) {
		try {
			return matrix.get(rowKey).get(columnKey);
		} catch (final Exception e) {
			return null;
		}
	};

	public T add(K rowKey, K columnKey, T object) {
		return add(rowKey, columnKey, object, false);
	}

	public T add(K rowKey, K columnKey, T object, boolean addSymetrically) {
		if (!matrix.containsKey(rowKey)) {
			rowKeys.add(rowKey);
			matrix.put(rowKey, new THashMap<K, T>());
		}
		final T previousValue = matrix.get(rowKey).put(columnKey, object);
		if (previousValue == null) {
			if (!colKeysSet.contains(columnKey)) {
				colKeysSet.add(columnKey);
				colKeys.add(columnKey);
			}
		}
		isSymetrix = null;
		if (addSymetrically) {
			add(columnKey, rowKey, object);
		}
		return previousValue;
	};

	public int getNCols() {
		return colKeys.size();
	}

	public int getNRows() {
		return rowKeys.size();
	}

	public boolean isSymetric() {
		if (isSymetrix == null) {
			if (rowKeys.size() == colKeys.size()) {
				for (int i = 0; i < rowKeys.size(); i++) {
					if (!rowKeys.get(i).equals(colKeys.get(i))) {
						isSymetrix = false;
						return isSymetrix;
					}
				}
				isSymetrix = true;
			} else {
				isSymetrix = false;
			}
		}
		return isSymetrix;
	}

	/**
	 * Function to override that provides the String object that will be printed
	 * from the T objects in the output file when calling saveToFile(File)
	 * 
	 * @param t
	 * @return
	 */
	public abstract String printValue(T t);

	public void saveToFile(File file) throws IOException {
		if (file == null) {
			throw new IllegalArgumentException("file is null");
		}
		if (!file.getParentFile().exists()) {
			log.info("Creating folder " + file.getParentFile().getAbsolutePath());
			file.getParentFile().mkdirs();
		}
		FileWriter fw = null;
		try {
			fw = new FileWriter(file);
			// header
			for (final K col : colKeys) {
				fw.write("\t" + col.toString());
			}
			fw.write(newlineString);
			// row by row
			for (final K row : rowKeys) {
				fw.write(row.toString() + "\t");
				for (final K col : colKeys) {
					final T t = getValue(row, col);
					if (t != null) {
						fw.write(printValue(t));
					} else {
						fw.write(emptyValue.toString());
					}
					fw.write("\t");
				}
				fw.write(newlineString);
			}
		} finally {
			fw.close();
			log.info("Matrix printed at: " + file.getAbsolutePath() + " ("
					+ FileUtils.getDescriptiveSizeFromBytes(file.length()) + ")");
		}
	}

	public void removeRow(int rowIndex) {
		final K rowKey = rowKeys.get(rowIndex);
		matrix.remove(rowKey);
		rowKeys.remove(rowIndex);
		isSymetrix = null;
	}

	public void removeColumn(int colIndex) {
		final K colKey = colKeys.get(colIndex);
		for (final K rowKey : rowKeys) {
			matrix.get(rowKey).remove(colKey);
		}

		colKeys.remove(colIndex);
		colKeysSet.remove(colKey);
		isSymetrix = null;
	}

	public void sortColumns(Comparator<K> comparator) {
		Collections.sort(colKeys, comparator);
	}

	public void sortRows(Comparator<K> comparator) {
		Collections.sort(rowKeys, comparator);
	}

	public String getNewlineString() {
		return newlineString;
	}

	public void setNewlineString(String newlineString) {
		this.newlineString = newlineString;
	}

	public List<K> getColNames() {
		return colKeys;
	}

	public List<K> getRowNames() {
		return rowKeys;
	}
}
