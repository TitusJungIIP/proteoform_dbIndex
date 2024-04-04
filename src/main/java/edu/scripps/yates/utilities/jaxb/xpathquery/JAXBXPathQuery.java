package edu.scripps.yates.utilities.jaxb.xpathquery;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import org.apache.log4j.Logger;

import gnu.trove.set.hash.THashSet;

/**
 * Class that perform a query over an object created with JAXB from an XML.<br>
 * It will return a list of values specified by two XPaths, one that indicates
 * the element (and optionally perform a value check) and the other which
 * specifies the final element to return
 * 
 * @author Salva
 *
 */
public class JAXBXPathQuery {
	private final static Logger log = Logger.getLogger(JAXBXPathQuery.class);

	/**
	 * Example, using the Entry object from uniprot xml schema: String xpath =
	 * "$id"; String filterXPath = "dbReference$type=GO";
	 * 
	 * @param jaxbObject
	 *            the jaxb object
	 * @param xpath
	 *            the main XPath string. It can contain an attribute
	 *            (dbReference$type). If it contains an attribute (separated by
	 *            $), a condition in the value of the condition can be written
	 *            as "dbReference$type=GO"
	 * @param subXpath
	 *            an XPath applied to the deepest element in the queryXPath.
	 *            This is the value that will be returned.
	 * @return
	 */
	public static List<String> query(Object jaxbObject, String xpath, String subXpath) {
		final XPathParser xpathParserFilter = new XPathParser(xpath);
		final XPathParser xpathParser = new XPathParser(subXpath);

		final Set<Object> set = new THashSet<Object>();
		set.add(jaxbObject);
		try {
			final List<Object> explore = explore(xpathParserFilter, xpathParser, set, set);
			final List<String> ret = new ArrayList<String>();
			for (final Object object : explore) {
				if (!ret.contains(object.toString())) {
					ret.add(object.toString());
				}
			}
			return ret;
		} catch (NoSuchMethodException | SecurityException | IllegalAccessException | IllegalArgumentException
				| InvocationTargetException | CloneNotSupportedException e) {
			log.error(e);
		}
		return Collections.emptyList();
	}

	private static List<Object> explore(XPathParser xPathParserFilter, XPathParser xPathParser,
			Collection<Object> originalObjects, Collection<Object> objs)
			throws NoSuchMethodException, SecurityException, IllegalAccessException, IllegalArgumentException,
			InvocationTargetException, CloneNotSupportedException {
		final List<Object> list = new ArrayList<Object>();
		final boolean hasNext = xPathParserFilter.hasNext();
		XPathElement element = null;
		if (hasNext) {
			element = xPathParserFilter.next();
		}
		final String filterValue = xPathParserFilter.getFilterValue();
		for (final Object object : objs) {
			if (hasNext) {
				final String methodToCall = "get" + capitalizeFirstLetter(element.getName());
				final Method method = object.getClass().getMethod(methodToCall, null);
				final Object returnedObj = method.invoke(object, null);
				if (returnedObj != null && returnedObj instanceof Collection) {
					final Collection<Object> collection = (Collection<Object>) returnedObj;
					if (filterValue != null) {
						for (final Object object2 : collection) {
							final String methodToCall2 = "get" + capitalizeFirstLetter(filterValue);
							final Method method2 = object2.getClass().getMethod(methodToCall2, null);
							final Object returnedObj2 = method2.invoke(object2, null);
							if (returnedObj2 != null) {
								if (filterValue != null) {
									if (filterValue.equals(returnedObj2.toString())) {
										final Set<Object> set = new THashSet<Object>();
										set.add(object);
										list.addAll(explore((XPathParser) xPathParser.clone(), null, null, set));
									}
								} else {
									list.add(returnedObj2);
								}
							}
						}
					}
					list.addAll(explore(xPathParserFilter, xPathParser, originalObjects, collection));

				} else {
					if (filterValue != null) {
						if (filterValue.equals(returnedObj)) {
							final Set<Object> set = new THashSet<Object>();
							set.add(object);
							list.addAll(explore(xPathParserFilter, xPathParser, originalObjects, set));
						}
					} else {
						final Set<Object> set = new THashSet<Object>();
						set.add(returnedObj);
						list.addAll(explore(xPathParserFilter, xPathParser, originalObjects, set));
					}
				}

			} else {
				// getValue
				final Method method = object.getClass().getMethod("getValue", null);
				final Object returnedObj = method.invoke(object, null);
				if (returnedObj != null) {
					if (filterValue != null) {
						if (filterValue.equals(returnedObj)) {
							list.addAll(explore(xPathParser, null, null, originalObjects));
						}
					} else {
						list.add(returnedObj);
					}
				}
			}
		}
		return list;

	}

	private static String capitalizeFirstLetter(String string) {
		final String string2 = Character.toUpperCase(string.charAt(0)) + string.substring(1);
		return string2;
	}
}
