package edu.scripps.yates.proteoform_dbindex;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.annotations.uniprot.proteoform.Proteoform;
import edu.scripps.yates.annotations.uniprot.proteoform.ProteoformType;
import edu.scripps.yates.annotations.uniprot.proteoform.ProteoformUtil;
import edu.scripps.yates.annotations.uniprot.proteoform.xml.UniprotProteoformRetrieverFromXML;
import edu.scripps.yates.proteoform_dbindex.util.ProteoformDBIndexUtil;
import gnu.trove.map.hash.TIntObjectHashMap;
import umich.ms.fileio.filetypes.agilent.cef.jaxb.P;

import java.io.File;
import java.util.List;
import java.util.Map;

public class Scratch {

    public static void main(String []args)
    {
        File dir = new File("/home/yateslab/projectData/prolucid/2601LatestProteoform");
        final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(
                dir, true);
        final String uniprotVersion = "2024_01";
        UniprotProteoformRetrieverFromXML proteoformRetriever = new UniprotProteoformRetrieverFromXML(uplr, uniprotVersion);
        String protAccession = "Q9BRY0";
        Map<String, List<Proteoform>>  proteoformMap = proteoformRetriever.getProteoforms(protAccession);
        final List<Proteoform> proteoforms = proteoformMap.get(protAccession);

        final List<Proteoform> isoformProteoforms = ProteoformUtil.getProteoformsAs(proteoforms,
                ProteoformType.ISOFORM);
        // others here
        final List<Proteoform> nonIsoformProteoforms = ProteoformUtil.getProteoformsDifferentThan(proteoforms,
                ProteoformType.ISOFORM);

        final TIntObjectHashMap<List<Proteoform>> nonIsoformsProteoformsByPositionInMainProtein = ProteoformDBIndexUtil
                .getInstance().getProteoformsByPositionInProtein(nonIsoformProteoforms);
        for(Map.Entry<String,List<Proteoform>> entry: proteoformMap.entrySet())
        {
            System.out.println(">> "+entry.getKey());
            for(Proteoform proteoform: entry.getValue())
            {
                System.out.println(">> "+entry.getKey() +"\t"+proteoform);

            }
        }


    }
}
