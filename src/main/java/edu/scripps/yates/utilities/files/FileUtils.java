package edu.scripps.yates.utilities.files;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.nio.charset.MalformedInputException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;
import org.apache.poi.openxml4j.exceptions.InvalidFormatException;
import org.apache.poi.openxml4j.opc.OPCPackage;
import org.apache.poi.openxml4j.opc.PackageAccess;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellType;
import org.apache.poi.ss.usermodel.DateUtil;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.xssf.usermodel.XSSFCell;
import org.apache.poi.xssf.usermodel.XSSFRow;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.apache.xmlbeans.ThreadLocalUtil;

import edu.scripps.yates.utilities.dates.DatesUtil;

public class FileUtils {
	private final static Logger log = Logger.getLogger(FileUtils.class);
	private static DecimalFormat df;

	public static void mergeFiles(Collection<File> files, File mergedFile, boolean skipHeaderOfNotFirstFiles) {
		final File[] fileArray = new File[files.size()];
		int i = 0;
		for (final File file : files) {
			fileArray[i] = file;
			i++;
		}
		mergeFiles(fileArray, mergedFile, skipHeaderOfNotFirstFiles);
	}

	private static DecimalFormat getDecimalFormat() {
		if (df == null) {
			df = new DecimalFormat("#.#");
		}
		return df;
	}

	public static void mergeFiles(File[] files, File mergedFile, boolean skipHeaderOfNotFirstFiles) {
		// if the file already exists, delete it
		if (mergedFile.exists()) {
			mergedFile.delete();
		}
		FileWriter fstream = null;
		BufferedWriter out = null;
		try {
			try {
				fstream = new FileWriter(mergedFile, true);
				out = new BufferedWriter(fstream);
			} catch (final IOException e1) {
				e1.printStackTrace();
			}

			boolean firstFile = true;
			for (final File f : files) {
				log.debug("merging: " + f.getName());
				FileInputStream fis;
				BufferedReader in = null;
				try {
					fis = new FileInputStream(f);
					in = new BufferedReader(new InputStreamReader(fis));

					String aLine;
					boolean firstLine = true;
					while ((aLine = in.readLine()) != null) {
						if (skipHeaderOfNotFirstFiles && (!firstFile && firstLine)) {
							firstLine = false;
							continue;
						}

						out.write(aLine);
						out.newLine();

					}
					firstLine = false;

				} catch (final IOException e) {
					e.printStackTrace();
				} finally {
					if (in != null) {
						try {
							in.close();
						} catch (final IOException e) {
							e.printStackTrace();
							log.error(e);
						}
					}

				}
				firstFile = false;
			}
		} finally {
			try {
				if (out != null) {
					out.close();
					log.debug("File merged at: " + mergedFile.getAbsolutePath());
				}
			} catch (final IOException e) {
				e.printStackTrace();
				log.error(e);
			}
		}

	}

	/**
	 * Delete a folder if it is empty. if recursivelly is true, also delete the
	 * parent folders if they are empty.
	 *
	 * @param folder
	 * @param recursivelly
	 */
	public static void removeFolderIfEmtpy(File folder, boolean recursivelly) {
		if (folder == null || !folder.isDirectory())
			return;
		if (folder.exists() && folder.list().length == 0) {
			folder.delete();

			if (recursivelly) {
				removeFolderIfEmtpy(folder.getParentFile(), recursivelly);
			}
		}
	}

	/**
	 * By default File#delete fails for non-empty directories, it works like "rm".
	 * We need something a little more brutual - this does the equivalent of "rm -r"
	 *
	 * @param path Root File Path
	 * @return true iff the file and all sub files/directories have been removed
	 * @throws FileNotFoundException
	 */
	public static boolean deleteFolderRecursive(File path) throws FileNotFoundException {
		if (!path.exists())
			throw new FileNotFoundException(path.getAbsolutePath());
		boolean ret = true;
		if (path.isDirectory()) {
			for (final File f : path.listFiles()) {
				ret = ret && FileUtils.deleteFolderRecursive(f);
			}
		}
		return ret && path.delete();
	}

	public static String getDescriptiveSizeFromBytes(long sizeInBytes) {
		log.debug("Getting file descriptive string from a number of bytes: " + sizeInBytes);
		if (sizeInBytes < 1024) {
			return getDecimalFormat().format(sizeInBytes) + " bytes";
		}
		final double sizeInKBytes = sizeInBytes / 1024;
		if (sizeInKBytes < 1024) {
			return getDecimalFormat().format(sizeInKBytes) + " Kb";
		}
		final double sizeInMBytes = sizeInKBytes / 1024;
		if (sizeInMBytes < 1024) {
			return getDecimalFormat().format(sizeInMBytes) + " Mb";
		}
		final double sizeInGBytes = sizeInMBytes / 1024;
		if (sizeInGBytes < 1024) {
			return getDecimalFormat().format(sizeInGBytes) + " Gb";
		}
		final double sizeInTBytes = sizeInGBytes / 1024;
		if (sizeInTBytes < 1024) {
			return getDecimalFormat().format(sizeInTBytes) + " Tb";
		}
		final double sizeInPBytes = sizeInTBytes / 1024;
		// if (sizeInPBytes < 1024) {
		return getDecimalFormat().format(sizeInPBytes) + " Pb";
		// }
	}

	/**
	 * returns a {@link File} object where having an input file like
	 * /path/myfile.extension, it returns /path/prefixmyfilesuffix.extension
	 *
	 * @param file
	 * @param prefix
	 * @param sufix
	 * @return
	 */
	public static File appendToFileName(File file, String prefix, String suffix) {
		if (file != null) {
			final String absolutePath = file.getAbsolutePath();
			final String baseName = FilenameUtils.getBaseName(absolutePath);
			final String extension = FilenameUtils.getExtension(absolutePath);
			final String prefix2 = prefix != null ? prefix : "";
			final String suffix2 = suffix != null ? suffix : "";
			return new File(file.getParentFile().getAbsolutePath() + File.separator + prefix2 + baseName + suffix2 + "."
					+ extension);
		}
		return null;
	}

	public static List<String> readColumnFromTextFile(File inputFile, final String separator, int columnIndex,
			String startAfterLineBegginingBy) throws IOException {
		int index = 0;

		Stream<String> stream = null;
		try {
			stream = Files.lines(Paths.get(inputFile.getAbsolutePath()), Charset.defaultCharset());
			final Iterator<String> iterator = stream.iterator();
			while (iterator.hasNext()) {
				index++;
				final String line = iterator.next();
				if (line.startsWith(startAfterLineBegginingBy)) {
					break;
				}
			}
		} finally {
			if (stream != null) {
				stream.close();
			}
		}

		return readColumnFromTextFile(inputFile, separator, columnIndex, index);
	}

	public static List<String> readColumnFromTextFile(File inputFile, final String separator, int columnIndex,
			int skipFirstLines) throws IOException {
		List<String> lines = null;
		try {
			lines = Files.readAllLines(Paths.get(inputFile.getAbsolutePath()), Charset.defaultCharset());
			final List<String> ret = new ArrayList<String>();
			int numLine = 0;
			for (final String line : lines) {
				numLine++;
				if (numLine <= skipFirstLines) {
					continue;
				}
				final String[] split = line.split(separator);
				if (split.length > columnIndex) {
					ret.add(split[columnIndex]);
				} else {
					ret.add(null);
				}
			}

			return ret;
		} catch (final IndexOutOfBoundsException e) {
			log.error(e);
			return null;
		}
	}

	public static List<String> readColumnFromTextFile(File inputFile, final String separator, int columnIndex,
			boolean skipHeader) throws IOException {
		Stream<String> lines = null;
		try {
			lines = Files.lines(Paths.get(inputFile.getAbsolutePath()), Charset.defaultCharset());
			if (skipHeader) {
				lines = lines.skip(1);
			}
			final List<String> collect = lines.map(line -> line.split(separator)[columnIndex])
					.collect(Collectors.toList());
			return collect;
		} catch (final IndexOutOfBoundsException e) {
			log.error(e);
			return null;
		} finally {
			if (lines != null) {
				lines.close();
			}
		}
	}

	public static File getFileFromInputStream(InputStream is) throws IOException {
		final File tempFile = File.createTempFile("temp", "tmp");
		return getFileFromInputStream(is, tempFile);
	}

	public static File getFileFromInputStream(InputStream is, File targetFile) throws IOException {

		Files.copy(is, targetFile.toPath(), StandardCopyOption.REPLACE_EXISTING);

		IOUtils.closeQuietly(is);
		return targetFile;
	}

	public static List<String> readFirstLines(InputStream inputStream, long maxLines) {
		final long t1 = System.currentTimeMillis();
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new InputStreamReader(inputStream));
			String line = null;
			final List<String> ret = new ArrayList<String>();
			while ((line = reader.readLine()) != null) {
				ret.add(line);
				if (ret.size() == maxLines) {
					break;
				}
			}
			return ret;
		} catch (final IOException e) {
			e.printStackTrace();
			return Collections.emptyList();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (final IOException e) {
					e.printStackTrace();
				}
			}
			log.info("File readed in " + DatesUtil.getDescriptiveTimeFromMillisecs((System.currentTimeMillis() - t1)));
		}
	}

	/**
	 * Read the first (maxLines) lines from the file.<br>
	 * It supports GZIP and ZIP files
	 * 
	 * @param file
	 * @param maxLines
	 * @return
	 */
	public static List<String> readFirstLines(File file, long maxLines) {

		InputStream inputStream;
		try {
			inputStream = getInputStream(file);
			if (inputStream != null) {
				return readFirstLines(inputStream, maxLines);
			}
		} catch (final FileNotFoundException e) {

		}

		return Collections.emptyList();
	}

	/**
	 * Gets the appropriate {@link InputStream} from a {@link File}, taking into
	 * account GZIP or ZIP files
	 * 
	 * @param file
	 * @return
	 * @throws FileNotFoundException
	 */
	public static InputStream getInputStream(File file) throws FileNotFoundException {

		InputStream in = null;
		try {
			in = new FileInputStream(file);
			final GZIPInputStream gzip = new GZIPInputStream(in);
			return gzip;
		} catch (final IOException e) {
			try {
				in.close();
			} catch (final IOException ee) {
			}
			ZipFile zip = null;
			try {
				zip = new ZipFile(file);
				final Enumeration<? extends ZipEntry> entries = zip.entries();
				while (entries.hasMoreElements()) {
					final ZipEntry zipEntry = entries.nextElement();
					if (!zipEntry.isDirectory()) {
						in = zip.getInputStream(zipEntry);
						return in;
					}
				}
			} catch (final IOException e2) {

				if (zip != null) {
					try {
						zip.close();
					} catch (final IOException e1) {
					}
				}

			} finally {
				// do not close zip, because that would close the stream
			}
			in = new FileInputStream(file);
			return in;
		}
	}

	/**
	 * Convert a TSV or CSV to Excel XLSX format
	 * 
	 * @param csvFilePath
	 * @param outputXlsFilePath
	 * @param separator
	 * @param sheetName         if null, a new sheet will be created with standard
	 *                          name "sheet1", "sheet2"...
	 * @throws IOException
	 */
	public static XSSFWorkbook separatedValuesToXLSX(String csvFilePath, String outputXlsFilePath, String separator,
			String sheetName) throws IOException {
		log.info("Converting file '" + csvFilePath + "' to Excel file at '" + outputXlsFilePath + "'");
		final boolean createNewExcelFile = !(new File(outputXlsFilePath).exists());
		if (createNewExcelFile) {
			log.info("File " + outputXlsFilePath + " already exists. The new sheet will be added");
		}
		XSSFWorkbook workBook = null;

		// if the excel file already exist, then create the XSSFWorkbook
		// with the constructor using the path,
		if (!createNewExcelFile) {
			workBook = new XSSFWorkbook(new FileInputStream(outputXlsFilePath));
		} else {
			workBook = new XSSFWorkbook();
		}

		BufferedReader br = null;
		FileOutputStream fileOutputStream = null;
		try {
			XSSFSheet sheet = null;
			if (sheetName == null) {
				// figure out what is the next sheet name available
				for (int i = 0; i < workBook.getNumberOfSheets(); i++) {
					sheetName = "sheet" + (i + 1);
					if (workBook.getSheet(sheetName) == null) {
						break;
					}
				}
			}
			sheet = workBook.createSheet(sheetName);

			String currentLine = null;
			int RowNum = -1;
			br = new BufferedReader(new FileReader(csvFilePath));
			while ((currentLine = br.readLine()) != null) {
				final String str[] = currentLine.split(separator);
				RowNum++;
				final XSSFRow currentRow = sheet.createRow(RowNum);
				for (int i = 0; i < str.length; i++) {
					try {
						// try as number first
						final double d = Double.parseDouble(str[i]);
						currentRow.createCell(i).setCellValueImpl(d);
					} catch (final NumberFormatException e) {
						final XSSFCell newcell = currentRow.createCell(i);
						newcell.setCellValue(str[i]);
					}

				}
			}

			fileOutputStream = new FileOutputStream(outputXlsFilePath);
			workBook.write(fileOutputStream);
			fileOutputStream.close();
			workBook.close();

			ThreadLocalUtil.clearAllThreadLocals();
			log.info("Excel file written successfully at: " + outputXlsFilePath);

		} catch (

		final IOException ex) {
			log.error(ex);
			throw ex;
		} finally {
			if (br != null) {
				try {
					br.close();
				} catch (final IOException e) {
					e.printStackTrace();
				}
			}

		}
		return workBook;
	}

	public static void separatedValuesToXLSX(String csvFilePath, String outputXlsFilePath, String separator)
			throws IOException {
		separatedValuesToXLSX(csvFilePath, outputXlsFilePath, separator, "sheet1");
	}

	public static List<String> readLinesFromXLSX(String excelFile, String cellSeparator, Integer sheetIndex)
			throws FileNotFoundException, IOException, InvalidFormatException {
		return readLinesFromXLSX(new File(excelFile), cellSeparator, sheetIndex);
	}

	public static List<String> readLinesFromXLSX(String excelFile, String cellSeparator)
			throws FileNotFoundException, IOException, InvalidFormatException {
		return readLinesFromXLSX(excelFile, cellSeparator, 0);
	}

	/**
	 * Tells whether a file is an Excel file by trying to open it with <br>
	 * <code>OPCPackage pkg = OPCPackage.open(file,PackageAccess.READ)</code>
	 * 
	 * @param potentialExcelFile
	 * @return
	 */
	public static boolean isExcelFile(File potentialExcelFile) {
		try {
			final OPCPackage pkg = OPCPackage.open(potentialExcelFile.getAbsolutePath(), PackageAccess.READ);
			pkg.revert();
			return true;
		} catch (final Exception e) {
			return false;
		}
	}

	/**
	 * Returns a list of string in which each string represents the values of the
	 * cells of a row, separated by a cellSeparator
	 * 
	 * @param excelFile
	 * @param cellSeparator
	 * @param sheetIndex
	 * @return null if it is not a excel file. Empty list if there is no data in the
	 *         excel file. a list if there is data
	 * @throws FileNotFoundException
	 * @throws IOException
	 * @throws InvalidFormatException
	 */
	public static List<String> readLinesFromXLSX(File excelFile, String cellSeparator, Integer sheetIndex)
			throws FileNotFoundException, IOException {
		final List<String> lines = new ArrayList<String>();

		if (!isExcelFile(excelFile)) {
			return null;
		}
		OPCPackage pkg = null;
		try {

			pkg = OPCPackage.open(excelFile.getAbsolutePath(), PackageAccess.READ);

			final XSSFWorkbook workBook = new XSSFWorkbook(pkg);
			int index = 0;
			if (sheetIndex != null && sheetIndex >= 0) {
				index = sheetIndex;
			}
			int rowIndex = 0;
			String rowValues = getStringValuesFromRow(workBook, index, rowIndex, cellSeparator);
			while (!rowValues.isEmpty()) {
				lines.add(rowValues);
				rowIndex++;
				rowValues = getStringValuesFromRow(workBook, index, rowIndex, cellSeparator);
			}
			return lines;
		} catch (final InvalidFormatException e) {
			throw new IOException(e);
		} finally {
			if (pkg != null) {
				// close without saving changes
				pkg.revert();
			}
		}
	}

	/**
	 * Gets a string with the values of the row, separated by a provided separator
	 * 
	 * @param wb
	 * @param sheetIndex
	 * @param rowIndex
	 * @param cellSeparator
	 * @return
	 * @throws IOException
	 */
	public static String getStringValuesFromRow(XSSFWorkbook wb, int sheetIndex, int rowIndex, String cellSeparator)
			throws IOException {
		final StringBuilder sb = new StringBuilder();
		final Sheet sheet = wb.getSheetAt(sheetIndex);
		final Row row = sheet.getRow(rowIndex);
		if (row != null) {
			final int lastCellNumber = row.getLastCellNum();
			for (int i = 0; i < lastCellNumber; i++) {
				if (!"".equals(sb.toString())) {
					sb.append(cellSeparator);
				}
				try {
					final Cell cell = row.getCell(i, Row.MissingCellPolicy.RETURN_NULL_AND_BLANK);
					if (cell != null) {
						final CellType cellType = cell.getCellType();
						switch (cellType) {

						case STRING:
						case BLANK:
							sb.append(cell.getStringCellValue());
							break;
						case NUMERIC:
							if (DateUtil.isCellDateFormatted(cell)) {
								sb.append(cell.getDateCellValue().toString());
							} else {
								sb.append(String.valueOf(cell.getNumericCellValue()));
							}
							break;
						case FORMULA:
							sb.append("" + cell.getNumericCellValue());
							break;
						case BOOLEAN:
							sb.append(String.valueOf(cell.getBooleanCellValue()));
							break;
						case ERROR:
							sb.append("");
						case _NONE:
							sb.append("");
						default:
							sb.append("");
						}
					} else {
						sb.append("");
					}
				} catch (final Exception e) {
					sb.append("");
				}

			}

		}

		return sb.toString();
	}

	/**
	 * Transposes a table and writes it in a file
	 * 
	 * @param fileToTranspose
	 * @param fileToWrite
	 * @throws IllegalArgumentException if the table has different number of columns
	 *                                  in any row compared to the first one.
	 */
	public static void transposeTable(File fileToTranspose, File fileToWrite, String separator) throws IOException {
		final List<List<String>> matrix = new ArrayList<List<String>>();
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(fileToTranspose));

			String line;
			int size = -1;
			int numLine = 0;
			while ((line = reader.readLine()) != null) {
				numLine++;
				final List<String> list = Arrays.asList(line.split(separator));
				matrix.add(list);
				if (size == -1) {
					size = list.size();
				} else {
					if (size != list.size()) {
						throw new IllegalArgumentException("Different number of columns in line number " + numLine + "("
								+ list.size() + ") compared to columns in row 0 (" + size + ")");
					}
				}
			}
			log.debug("File readed with " + matrix.size() + " and lines of " + matrix.get(0).size() + " columns");
		} finally {
			if (reader != null) {
				reader.close();
			}
		}

		// write now
		log.debug("Writting transposed table in file " + fileToWrite.getAbsolutePath());
		final FileWriter fw = new FileWriter(fileToWrite);
		for (int j = 0; j < matrix.get(0).size(); j++) {
			for (int i = 0; i < matrix.size(); i++) {
				fw.write(matrix.get(i).get(j));
				if (i == matrix.size() - 1) {
					fw.write("\n");
				} else {
					fw.write(separator);
				}
			}
		}
		fw.close();
	}

	/**
	 * Read lines of file (Path) trying all the standard charsets. So this is a
	 * error-less (related to charset) read lines utility method.
	 * {@link StandardCharsets}
	 * 
	 * @param path
	 * @return a list of string, each of them being a line from the file.
	 * @throws IOException
	 */
	public static List<String> readAllLines(Path path) throws IOException {
		final List<Charset> charsets = new ArrayList<Charset>();
		charsets.add(StandardCharsets.ISO_8859_1);
		charsets.add(StandardCharsets.US_ASCII);
		charsets.add(StandardCharsets.UTF_16);
		charsets.add(StandardCharsets.UTF_16BE);
		charsets.add(StandardCharsets.UTF_16LE);
		charsets.add(StandardCharsets.UTF_8);
		int i = 0;
		IOException exception = null;
		while (i < charsets.size()) {
			try {
				final List<String> lines = Files.readAllLines(path, charsets.get(i));
				return lines;
			} catch (final MalformedInputException e) {
				exception = e;
			} catch (final IOException e) {
				e.printStackTrace();
				exception = e;

			} finally {
				i++;
			}
		}
		throw exception;
	}

	private static char[] invalidFileNameCharacters = { 34, 60, 62, 124, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
			14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 58, 42, 63, 92, 47 };
	private static char[] invalidFolderNameCharacters = { 34, 60, 62, 124, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
			14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31 };

	public static char[] getInvalidFileNameCharacters() {
		return invalidFileNameCharacters;
	}

	public static char[] getInvalidFolderNameCharacters() {
		return invalidFolderNameCharacters;
	}

	public static String checkInvalidCharacterNameForFileName(String fileName) {
		String ret = fileName;
		for (final char c : invalidFileNameCharacters) {
			ret = ret.replace(c, '_');
		}
		return ret;
	}

	public static String checkInvalidCharacterNameForFolderName(String fileName) {
		String ret = fileName;
		for (final char c : invalidFolderNameCharacters) {
			ret = ret.replace(c, '_');
		}
		return ret;
	}

}
