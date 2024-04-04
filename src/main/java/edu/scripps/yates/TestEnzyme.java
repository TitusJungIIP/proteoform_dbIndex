package edu.scripps.yates;

import edu.scripps.yates.utilities.sequence.MyEnzyme;

import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class TestEnzyme {

        public static final String protein = "MREIVHIQAGQCGNQIGAKFWEVISDEHGIDPTGSYHGDSDLQLERINVYYNEAAGNKYVPRAILVDLEPGTMDSVRSGPFGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELVDSVLDVVRKESESCDCLQGFQLTHSLGGGTGSGMGTLLISKIREEYPDRIMNTFSVMPSPKVSDTVVEPYNATLSVHQLVENTDETYSIDNEALYDICFRTLKLTTPTYGDLNHLVSATMSGVTTCLRFPGQLNADLRKLAVNMVPFPRLHFFMPGFAPLTSRGSQQYRALTVPELTQQMFDSKNMMAACDPRHGRYLTVAAIFRGRMSMKEVDEQMLNVQNKNSSYFVEWIPNNVKTAVCDIPPRGLKMSATFIGNSTAIQELFKRISEQFTAMFRRKAFLHWYTGEGMDEMEFTEAESNMNDLVSEYQQYQDATADEQGEFEEEEGEDEA";
        public static void main(String []args)
        {
            MyEnzyme  enzyme = new MyEnzyme(null, "KR",  null, "Cterm",
                    1);
            Collection<String> petpideList =  enzyme.cutSeq(protein, 6,100, null, true);
            int length = petpideList.size();
            Set<String> peptideSet = new HashSet<>();
            for(String peptide: petpideList)
            {
                peptideSet.add(peptide);
            //    System.out.println(peptide);
            }
            System.out.println(length);

            System.out.println(peptideSet.size());
        }


}


