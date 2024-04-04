package edu.scripps.yates.annotations.uniprot;

import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.annotations.uniprot.xml.ObjectFactory;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Uniprot;
import edu.scripps.yates.utilities.dates.DatesUtil;
import org.apache.commons.io.IOUtils;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.Buffer;
import java.util.Collections;
import java.util.List;

public class TestReadEntries {

    public static synchronized List<Entry> parseResponse(String is) throws JAXBException {


        if (is.startsWith("ERROR")) {
            return Collections.emptyList();
        }
        // OutputStream outputStream = null;
        JAXBContext jaxbContext = JAXBContext.newInstance(ObjectFactory.class);

        Unmarshaller unmarshaller = jaxbContext.createUnmarshaller();
        try {

            final long t2 = System.currentTimeMillis();
            // final File targetFile = File.createTempFile("uniprot", ".tmp");

            // java.nio.file.Files.copy(is, targetFile.toPath(),
            // StandardCopyOption.REPLACE_EXISTING);

            // final XMLStreamReader xmlStreamReader =
            // XMLInputFactory.newInstance().createXMLStreamReader(is);
            final Uniprot uniprot = (Uniprot) unmarshaller.unmarshal(IOUtils.toInputStream(is));


            return uniprot.getEntry();
            // } catch (final JAXBException e) {
            // e.printStackTrace();
            // log.warn("Error sending " + url);
            // log.warn(e.getMessage() + "\t" +
            // e.getLinkedException().getMessage());
        } catch (final JAXBException e) {
            throw e;
        } finally {
            // if (is != null) {
            // IOUtils.closeQuietly(is);
            // }
        }
    }

    public static void main(String []args) throws IOException, JAXBException {
        String path = "/home/yateslab/projectData/prolucid/0812ProteoformRegularCompare/test.xml";
        StringBuilder xmlFileContents = new StringBuilder();
        String line;
        BufferedReader br = new BufferedReader(new FileReader(path));
        while((line = br.readLine())!=null)
        {
            xmlFileContents.append(line).append("\n");
        }
        br.close();
        String testStr = UniprotEntryRetrieverThread.replaceDuplicateXMLHeader(xmlFileContents.toString());
        List<Entry> entries = parseResponse(testStr);
        for(Entry entry: entries)
        {
            System.out.println(entry.getAccession());
        }
    }
}

