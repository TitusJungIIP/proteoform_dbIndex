package edu.scripps.yates.proteoform_dbindex.model;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.junit.Test;

import junit.framework.Assert;

public class PTMTest {

	@Test
	public void extractPTMsFromSequence() throws IOException {
		final File ptmCodeFiles = new File(
				"D:\\titus\\UniProt_Human_sprot_11-08-2010_reversed.fasta_d0419382f1d8fa4cd8d1703f64597c4\\ptmCodes.txt");
		final ExtendedAssignMass extendedAssignMass = ExtendedAssignMass.getInstance(true, ptmCodeFiles);
		String modifiedSequence = "ABC[D->]EF[+30]";
		List<PTM> ptms = PTM.extractPTMsFromSequenceBasedOnResultingSequence(modifiedSequence, extendedAssignMass);
		Assert.assertEquals(2, ptms.size());
		PTM ptm1 = ptms.get(0);
		Assert.assertEquals(3, ptm1.getPosInPeptide());
		PTM ptm2 = ptms.get(1);
		Assert.assertEquals(5, ptm2.getPosInPeptide());

		//
		modifiedSequence = "ABC[D->KK]EF[+30]";
		ptms = PTM.extractPTMsFromSequenceBasedOnResultingSequence(modifiedSequence, extendedAssignMass);
		Assert.assertEquals(2, ptms.size());
		ptm1 = ptms.get(0);
		Assert.assertEquals(3, ptm1.getPosInPeptide());
		ptm2 = ptms.get(1);
		Assert.assertEquals(7, ptm2.getPosInPeptide());

		//
		modifiedSequence = "ABC[D->KK]EF[+30][R->][TT->EEE]";
		ptms = PTM.extractPTMsFromSequenceBasedOnResultingSequence(modifiedSequence, extendedAssignMass);
		Assert.assertEquals(4, ptms.size());
		ptm1 = ptms.get(0);
		Assert.assertEquals(3, ptm1.getPosInPeptide());
		ptm2 = ptms.get(1);
		Assert.assertEquals(7, ptm2.getPosInPeptide());
		final PTM ptm3 = ptms.get(2);
		Assert.assertEquals(7, ptm3.getPosInPeptide());
		final PTM ptm4 = ptms.get(3);
		Assert.assertEquals(7, ptm4.getPosInPeptide());
	}
}
