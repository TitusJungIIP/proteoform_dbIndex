package edu.scripps.yates.proteoform_dbindex;


import org.apache.commons.lang.math.NumberUtils;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

public class SearchParams {
    private String databaseName;
    private int primaryScoreType = 0;
    private int secondaryScoreType = 1;
    private int minMatch;
    private String precursorIsotope;
    private String fragmentIsotope;
    private int numIsotopicPeaks;
    private int preprocess;   // preprocess mode
    private double highPrecursorTolerance;
    private double lowPrecursorTolerance;
    private double fragmentTolerance;
    private String useLowFragmentIons = "no";
    //high_resolution_data
    private String highResolutionData = "";
    private double fragmentToleranceMillionth;
    private double precursorTolerance;
    private double maxPrecursorMass = 10000;
    private double minPrecursorMass = 600;
    private int minimumPeptideLength = 5;
    private int maxNumSpectra = 5000;
    private int minNumSpectra = 25;
    private int maxAlter;
    private int maxPrecursorCharge = 10000;
    private int minPrecursorCharge = 1;
    private boolean isDatabaseIndexed;
    private int peakRankThreshold;
    private int candidatePeptideThreshold;
    private int numOutput = 5;
    private int locusType = 0;
    private int displayMod = 0;
    private String atomicEnrichment = "0";
    //    private int incompletelabeling = 0;
    private int leftshift1 = 1000000;
    private int leftshift2 = 1000000;
    private int leftshift3 = 1000000;
    //    private Modifications mods = new Modifications();
    private boolean isDeCharged = false;
    private boolean isEtdSearch = false;
    private int multistageActivationMod = 0;
    private String paramFilePath;
    private String paramFileName;
    private boolean chargeDisambiguation = false;
    private int enzymeSpecificity;
    private int maxInternalMisCleavage = -1; // -1 for unlimited
    private String residues = null;
    public static final String AVGISOTOPE = "avg";
    public static final String MONOISOTOPE = "mono";
    private boolean lowFragIonMode = false;

    private double ntermStaticMod;
    private int offset =0;
    private double ctermStaticMod;

    private double maxMassShift = -10000;
    public static boolean isHighResolution = false;
    private String ptmDBName = null;
    private boolean useProteoformDB = false;
    private String uniprotVersion = null;

    public SearchParams(String path, String paramFile) throws IOException, JDOMException {
        paramFilePath = path;
        paramFileName = paramFile;
        readParams();
    }

    public SearchParams(String paramFile) throws IOException, JDOMException {
        paramFileName = paramFile;
        //paramFilePath = ".";
        readParams();
    }
    public SearchParams() {

    }
    public boolean isEtdSearch() {
        return isEtdSearch;
    }

    public void setMaxInternalMisCleavage(int msc) {
        maxInternalMisCleavage = msc;
    }

    public int getMaxInternalMisCleavage() {
        if(enzymeSpecificity == 0) return -1; // -1 for unlimited
        return maxInternalMisCleavage;
    }

    public boolean isDeCharged() {
        return isDeCharged;
    }
    public int getMultistageActivationMod() {
        return multistageActivationMod;
    }


    public double getStaticCTermMod() {
        return ctermStaticMod;
    }
    public double getStaticNTermMod() {
        return ntermStaticMod;
    }

    public int getNumOutput() {
        return numOutput;
    }



    private void readParams() throws IOException, JDOMException {
        String file = paramFilePath == null? paramFileName : paramFilePath + "/" + paramFileName;
        //File paraFile = new File(paramFilePath + "/" + paramFileName);
        File paraFile = new File(file);
        FileInputStream paramInput = new FileInputStream(paraFile);
        Document doc = new SAXBuilder().build(paramInput);

        Element root = doc.getRootElement();
        readDatabase(root.getChild("database"));
        readSearchMode(root.getChild("search_mode"));
        readIsotopes(root.getChild("isotopes"));
        readTolerance(root.getChild("tolerance"));
        readPrecursorMassLimits(root.getChild("precursor_mass_limits"));
        readPrecursorChargeLimits(root.getChild("precursor_charge_limits"));
        readPeptideLengthLimits(root.getChild("peptide_length_limits"));

        Element numpeaklimit = root.getChild("num_peak_limits");
        numpeaklimit = numpeaklimit == null? root.getChild("num_spectra_limits") : numpeaklimit;
        if(numpeaklimit != null) {
            readNumSpectraLimits(numpeaklimit);
        }

        String maxnumdiffmods = root.getChildTextTrim("max_num_diffmod");
        maxnumdiffmods = maxnumdiffmods == null? root.getChildTextTrim("max_alter") :  maxnumdiffmods;

        maxAlter = (maxnumdiffmods == null || "".equals(maxnumdiffmods))? 0 : Integer.parseInt(maxnumdiffmods);

        readEnzymeInfo(root.getChild("enzyme_info"));

        paramInput.close();
    }
    private void readEnzymeInfo(Element e) {
        enzymeSpecificity = Integer.parseInt(e.getChildTextTrim("specificity"));
        String maxmiscleavage =  e.getChildTextTrim("max_num_internal_mis_cleavage");

        System.out.println("Maximum number of internal mis : " + maxmiscleavage);
        if(maxmiscleavage != null && !"".equals(maxmiscleavage.trim())) {
            maxInternalMisCleavage = Integer.parseInt(maxmiscleavage);
        }
        System.out.println("maxInternalMisCleavage value is now: " + maxInternalMisCleavage);
        String enzymeName = e.getChildTextTrim("name");
        this.offset =  Boolean.parseBoolean(e.getChildTextTrim("type"))? 1:0;
        Element re = e.getChild("residues");
        StringBuilder sb = new StringBuilder();
        List residues = re.getChildren();
        List<Character> characters = new ArrayList<>();
        for(Iterator i = residues.iterator(); i.hasNext();) {
            Element r = (Element) i.next();
            String s = r.getTextTrim();
            characters.add(s.charAt(0));
        }
        characters.sort(Comparator.naturalOrder());
        for(char c: characters)
        {
            sb.append(c);
        }
        this.residues = sb.toString();
    }
    private void readDatabase(Element e) {
        databaseName = e.getChildTextTrim("database_name");
        String s = e.getChildTextTrim("is_indexed");
        if(s != null) {
            char c = s.charAt(0);
            if(c == 'T' || c == 't' || c == 'Y' || c == 'y') {
                isDatabaseIndexed = true;
            } else if(c == 'N' || c == 'n' || c == 'F' || c == 'f') {
                isDatabaseIndexed = false;
            }
        } else {
            isDatabaseIndexed = false; // by default, database is not indexed
        }
        String temp = e.getChildTextTrim("ptm_db_name");
        if(temp!=null)
            ptmDBName = temp.trim();
        String proteoformElement = e.getChildTextTrim("use_proteoform_database");
        if(proteoformElement!=null)
        {
            useProteoformDB = Boolean.parseBoolean(proteoformElement);
        }
        uniprotVersion = e.getChildTextTrim("uniprot_version");
    }

    private void readNumSpectraLimits(Element ne) {
        if(ne == null) return;
        maxNumSpectra =  Integer.parseInt(ne.getChildTextTrim("maximum"));
        minNumSpectra =  Integer.parseInt(ne.getChildTextTrim("minimum"));
    }
    private void readPrecursorMassLimits(Element pe) {
        if(pe == null) return;
//        maxPrecursorMass = Double.parseDouble(pe.getChildTextTrim("maximum"));
        maxPrecursorMass = Double.parseDouble(pe.getChildTextTrim("maximum"));
//        minPrecursorMass = Double.parseDouble(pe.getChildTextTrim("minimum"));
        minPrecursorMass = Double.parseDouble(pe.getChildTextTrim("minimum"));
    }
    private void readPrecursorChargeLimits(Element pe) {
//        maxPrecursorMass = Double.parseDouble(pe.getChildTextTrim("maximum"));
        if (pe == null) {
            return;
        }
        try {
            String precursorChargeMaxValue = pe.getChildTextTrim("maximum");
            String precursorChargeMinValue = pe.getChildTextTrim("minimum");
            if(NumberUtils.isNumber(precursorChargeMaxValue) && NumberUtils.isNumber(precursorChargeMinValue) ){
                maxPrecursorCharge = NumberUtils.createNumber(precursorChargeMaxValue).intValue();
                int cs = NumberUtils.createNumber(precursorChargeMinValue).intValue();
                minPrecursorCharge = cs > 0 ? cs : 1;
            }
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    private void readPeptideLengthLimits(Element pe) {
        if(pe == null) return;
        minimumPeptideLength = Integer.parseInt(pe.getChildTextTrim("minimum"));
    }

    private void readTolerance(Element te) {
//        precursorTolerance = Double.parseDouble(te.getChildTextTrim("precursor"));
        //precursorTolerance = Double.parseDouble(te.getChildTextTrim("precursor"));
        highPrecursorTolerance = Double.parseDouble(te.getChildTextTrim("precursor_high"));
        lowPrecursorTolerance = Double.parseDouble(te.getChildTextTrim("precursor_low"));
//        fragmentTolerance = Double.parseDouble(te.getChildTextTrim("fragment"));
        precursorTolerance = Double.parseDouble(te.getChildTextTrim("precursor"));
        fragmentTolerance = Double.parseDouble(te.getChildTextTrim("fragment"));

        useLowFragmentIons = te.getChildTextTrim("low_frag_ion");
        highResolutionData = te.getChildTextTrim("high_resolution_data");

        fragmentToleranceMillionth = fragmentTolerance/1_000_000.0;

        String lowFragIonModeString = te.getChildTextTrim("low_frag_ion");
        lowFragIonMode  = lowFragIonModeString == null ? false : Boolean.parseBoolean(lowFragIonModeString);

        String highResStr = te.getChildTextTrim("high_resolution_data");
        isHighResolution = highResStr== null ? false : Boolean.parseBoolean(highResStr);


    }
    private void readIsotopes(Element ie)
    {
        String mono = MONOISOTOPE;
        String avg = AVGISOTOPE;
        String isotope =  ie.getChildTextTrim("precursor");
        numIsotopicPeaks = Integer.parseInt(ie.getChildTextTrim("num_peaks"));
        // set the precursor isotope type
        if (mono.equals(isotope) || avg.equals(isotope)) {
            precursorIsotope = isotope;
        }

        // set fgragment isotope type
        isotope = ie.getChildTextTrim("fragment");
        if (mono.equals(isotope) || avg.equals(isotope)) {
            fragmentIsotope =  isotope;
        }
    }
    private void readSearchMode(Element me) {
        String primaryscoretype = me.getChildTextTrim("primary_score_type");
        String secondaryscoretype = me.getChildTextTrim("secondary_score_type");
        primaryScoreType = primaryscoretype == null? 0 : Integer.parseInt(primaryscoretype);
        if(secondaryscoretype != null) {
            secondaryScoreType = Integer.parseInt(secondaryscoretype);
        } else {
            secondaryScoreType = primaryScoreType == 0? 1 : 0;
        }
        String peakrankt = me.getChildTextTrim("peak_rank_threshold");
        peakRankThreshold =  peakrankt == null? 200: Integer.parseInt(peakrankt);

        String candidatet = me.getChildTextTrim("candidate_peptide_threshold");
        candidatePeptideThreshold = candidatet == null? 500: Integer.parseInt(candidatet);

        String numoutput = me.getChildTextTrim("num_output");
        numOutput = numoutput == null? 5 : Integer.parseInt(numoutput);
        minMatch =  Integer.parseInt(me.getChildTextTrim("min_match"));

        String preproc = me.getChildTextTrim("preprocess");
        preprocess  =  preproc == null? 1 : Integer.parseInt(preproc);
        String fragMethod = me.getChildTextTrim("fragmentation_method");
        isEtdSearch = fragMethod != null && (fragMethod.startsWith("ETD") || fragMethod.startsWith("etd"));

        String locus = me.getChildTextTrim("locus_type");
        locusType = locus == null? 0 : Integer.parseInt(locus);

        atomicEnrichment = me.getChildTextTrim("atomic_enrichement");
        String leftshift = atomicEnrichment;
        if(leftshift != null) {
            //System.out.println("Atomic enrichement: " + leftshift);

            double atomicEnrichment = Double.parseDouble(leftshift);
            int ae = 0;
            if(atomicEnrichment != 0) {
                ae = (int)(atomicEnrichment);
            }

            switch(ae) {
                case 90 : leftshift1 = 2; leftshift2 = 4; leftshift3 = 8; break;
                case 91 : leftshift1 = 2; leftshift2 = 5; leftshift3 = 9; break;
                case 92 : leftshift1 = 2; leftshift2 = 6; leftshift3 = 10; break;
                case 93 : leftshift1 = 2; leftshift2 = 6; leftshift3 = 11; break;
                case 94 : leftshift1 = 3; leftshift2 = 8; leftshift3 = 14; break;
                case 95 : leftshift1 = 3; leftshift2 = 9; leftshift3 = 16; break;
                case 96 : leftshift1 = 4; leftshift2 = 12; leftshift3 = 21; break;
                case 97 : leftshift1 = 5; leftshift2 = 16; leftshift3 = 29; break;
                case 98 : leftshift1 = 7; leftshift2 = 24; leftshift3 = 36; break;
                case 99 : leftshift1 = 10; leftshift2 = 36; leftshift3 = 1000; break;

            }
        }
        //System.out.println("leftshift1: " + leftshift1 + "\tleftshift2: " + leftshift2 + "\tleftshift3: " + leftshift3);
        String chargedisamb = me.getChildTextTrim("charge_disambiguation");
        if(chargedisamb != null && Integer.parseInt(chargedisamb) != 0) {
            chargeDisambiguation = true;
        }

        String decharged = me.getChildTextTrim("is_decharged");
        int dechargestate = decharged == null? 0 : Integer.parseInt(decharged);
        isDeCharged = dechargestate == 1? true : false;

        String multistageActivation = me.getChildTextTrim("multistage_activation_mode");
        multistageActivationMod = multistageActivation == null? 0 : Integer.parseInt(multistageActivation);

        String lowFragIonModeString = me.getChildTextTrim("low_frag_ion");
        lowFragIonMode  = lowFragIonModeString == null ? false : Boolean.parseBoolean(lowFragIonModeString);

        String highResStr = me.getChildTextTrim("high_resolution_data");
        isHighResolution = highResStr== null ? false : Boolean.parseBoolean(highResStr);

    }
    public boolean getChargeDisambiguation() {

        return chargeDisambiguation;
    }
    public int getLeftShift1() {
        return leftshift1;
    }
    public int getLeftShift2() {
        return leftshift2;
    }
    public int getLeftShift3() {
        return leftshift3;
    }
    public void setDbName(String dbName) {
        databaseName = dbName;
    }
    public void setPrimaryScoreType(int primaryScoreType) {
        this.primaryScoreType = primaryScoreType;
    }

    public void setsecondaryScoreType(int secondaryScoreType) {
        this.secondaryScoreType = secondaryScoreType;
    }
    public void setMinMatch(int minMatch) {
        this.minMatch = minMatch;
    }



    public void setPreprocess(int mode) {
        preprocess = mode;
    }
    public void setFragmentTolerance(double tolerance) {
        fragmentTolerance = tolerance;
    }
    public void setMaxPrecursorMass(double mass) {
        maxPrecursorMass = mass;
    }
    public void setMinPrecursorMass(double mass) {
        minPrecursorMass = mass;
    }
    public void setMaxPrecursorCharge(int charge) {
        maxPrecursorCharge = charge;
    }
    public void setMinPrecursorCharge(int charge) {
        minPrecursorCharge = charge;
    }
    public void setMinimumPeptideLength(int len) {
        minimumPeptideLength = len;
    }
    public void setMaxNumSpectra(int numSpectra) {
        maxNumSpectra = numSpectra;
    }
    public void setMinNumSpectra(int numSpectra) {
        minNumSpectra = numSpectra;
    }
    public void setMaxAlter(int maxAlter) {
        this.maxAlter = maxAlter;
    }
    public int getLocusType() {
        return locusType;
    }
    public int getPeakRankThreshold() {
        return peakRankThreshold;
    }
    public int getCandidatePeptideThreshold() {
        return candidatePeptideThreshold;
    }

    public String getDbName() {
        return databaseName;
    }
    public boolean isDatabaseIndexed() {
        return isDatabaseIndexed;
    }

    public int getPrimaryScoreType() {
        return primaryScoreType;
    }

    public int getSecondaryScoreType() {
        return secondaryScoreType;
    }

    public int getMinMatch() {
        return minMatch;
    }
    public int getNumIsotopicPeaks() {
        return numIsotopicPeaks;
    }
    public String getPrecursorIsotope() {
        return precursorIsotope;
    }

    public String getFragmentIsotope() {
        return fragmentIsotope;
    }
    public int getPreprocess() {
        return preprocess;
    }

    public double getHighPrecursorTolerance() {
        return highPrecursorTolerance;
    }
    public double getLowPrecursorTolerance() {
        return lowPrecursorTolerance;
    }
    public double getPrecursorTolerance() {
        return precursorTolerance;
    }
    public double getFragmentTolerance() {
        return fragmentTolerance;
    }

    public int getMinimumPeptideLength() {
        return minimumPeptideLength;
    }
    public double getMaxPrecursorMass() {
        return maxPrecursorMass;
    }
    public double getMinPrecursorMass() {
        return minPrecursorMass;
    }
    public int getMaxPrecursorCharge() {
        return maxPrecursorCharge;
    }
    public int getMinPrecursorCharge() {
        return minPrecursorCharge;
    }

    public int getMaxNumSpectra() {
        return maxNumSpectra;
    }
    public int getMinNumSpectra() {
        return minNumSpectra;
    }
    public int getMaxAlter() {
        return maxAlter;
    }

    public int getEnzymeSpecificity() {
        return enzymeSpecificity;
    }


    public double getStaticTerminalMods() {
        //return mods.getStaticTerminalMods();
        return ctermStaticMod + ntermStaticMod;
    }
    public String getDatabaseName() {
        return databaseName;
    }
    public int getDisplayMod () {
        return displayMod ;
    }
    public String getAtomicEnrichment () {
        return atomicEnrichment;
    }
    public void setDatabaseName(String dbname) {
        databaseName = dbname;
    }
    public String output() {
        StringBuffer sb = new StringBuffer(10000);
        sb.append("<?xml version=\"1.0\" encoding=\"ISO-8859-1\" ?>\n");
        sb.append("<!--  Parameters for ProLuCID Seach  -->\n");
        sb.append("<parameters>\n");
        sb.append("\t<database>\n");
        sb.append("\t\t<database_name>" + getDatabaseName() + "</database_name>\n");
        sb.append("\t</database>\n");

        sb.append("</parameters>\n");
        return sb.toString();
    }
    public static void main(String args[]) throws Exception {
        /*
        SearchParams sp = new SearchParams(args[0], args[1]);
        System.out.println("Database is: " + sp.getDbName());
        System.out.println("maxAlter is: " + sp.getMaxAlter());
        System.out.println("highPrecursorTolerance is: " + sp.getHighPrecursorTolerance());
        System.out.println("Fragment isotope is: " + sp.getFragmentIsotope());
        Iterator <Modification> modifications = sp.getStaticMods();
        while(modifications.hasNext()) {
            Modification m = modifications.next();
            System.out.println("Residue " + (char)m.getResidue() + "\t" + m.getMassShift() + "\n");
        }

        String a = "java";
        String b = new StringBuffer(a).toString();
        a = b;
        if(a.equals(b)) {
            System.out.println("a equal b");

        }
        if(a == b) {

            System.out.println("a == b");
        }else {

            System.out.println("a != b");
        }
        */
        SearchParams sp = new SearchParams("search.xml");
        System.out.println(sp.output());
        System.out.println(sp.getDbName());
    }

    public boolean isHighResolution() {
        return isHighResolution;
    }

    public void setHighResolution(boolean highResolution) {
        isHighResolution = highResolution;
    }

    public double getFragmentToleranceMillionth() {
        return fragmentToleranceMillionth;
    }


    public boolean isLowFragIonMode() {
        return lowFragIonMode;
    }

    public String getUseLowFragmentIons() {
        return useLowFragmentIons;
    }

    public void setUseLowFragmentIons(String useLowFragmentIons) {
        this.useLowFragmentIons = useLowFragmentIons;
    }

    public String getHighResolutionData() {
        return highResolutionData;
    }

    public void setHighResolutionData(String highResolutionData) {
        this.highResolutionData = highResolutionData;
    }

    public String getPtmDBName() {
        return ptmDBName;
    }

    public void setPtmDBName(String ptmDBName) {
        this.ptmDBName = ptmDBName;
    }

    public boolean isUseProteoformDB() {
        return useProteoformDB;
    }

    public String getUniprotVersion() {
        return uniprotVersion;
    }

    public String getResidues() {
        return residues;
    }

    public int getOffset() {
        return offset;
    }
}




