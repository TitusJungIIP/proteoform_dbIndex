package edu.scripps.yates.proteoform_dbindex.util;

import org.junit.Assert;
import org.junit.Test;

public class ProteoformUtilTest {

	@Test
	public void test() {
		final String sequenceWithModification = "VLATVTK[+100.016]PVGGD[K->Q]";
		final String newPeptideSequence = "VLATVTKPVGGDQNGGTRVVK";
		final String originalPeptideSequence = ProteoformDBIndexUtil
				.getOriginalPeptideSequence(sequenceWithModification, newPeptideSequence);
		Assert.assertEquals("VLATVTKPVGGDKNGGTRVVK", originalPeptideSequence);
	}

	@Test
	public void test2() {
		final String sequenceWithModification = "VLATVTK[+100.016]PVGGDK";
		final String newPeptideSequence = "VLATVTKPVGGDKNGGTRVVK";
		final String originalPeptideSequence = ProteoformDBIndexUtil
				.getOriginalPeptideSequence(sequenceWithModification, newPeptideSequence);
		Assert.assertEquals("VLATVTKPVGGDKNGGTRVVK", originalPeptideSequence);
	}

	@Test
	public void test3() {
		final String sequenceWithModification = "VLATVTKPVGGD[K->Q]N";
		final String newPeptideSequence = "VLATVTKPVGGDQNGGTRVVK";
		final String originalPeptideSequence = ProteoformDBIndexUtil
				.getOriginalPeptideSequence(sequenceWithModification, newPeptideSequence);
		Assert.assertEquals("VLATVTKPVGGDKNGGTRVVK", originalPeptideSequence);
	}

	@Test
	public void test4() {
		final String sequenceWithModification = "VLATVTK[+100.016]PVGGD[K->Q]N[+123.23]";
		final String newPeptideSequence = "VLATVTKPVGGDQNGGTRVVK";
		final String originalPeptideSequence = ProteoformDBIndexUtil
				.getOriginalPeptideSequence(sequenceWithModification, newPeptideSequence);
		Assert.assertEquals("VLATVTKPVGGDKNGGTRVVK", originalPeptideSequence);
	}

	@Test
	public void test5() {
		final String sequenceWithModification = "QLASGLLLVT[GPLVLNR->DLWSSIE]";
		final String newPeptideSequence = "QLASGLLLVTDLWSSIEVPLRR";
		final String originalPeptideSequence = ProteoformDBIndexUtil
				.getOriginalPeptideSequence(sequenceWithModification, newPeptideSequence);
		Assert.assertEquals("QLASGLLLVTGPLVLNRVPLRR", originalPeptideSequence);
	}

	@Test
	public void test6() {
		final String sequenceWithModification = "TCQSWSSMTPHWH[R->Q]";
		final String newPeptideSequence = "TCQSWSSMTPHWHQRTTEYYPNGGLTR";
		final String originalPeptideSequence = ProteoformDBIndexUtil
				.getOriginalPeptideSequence(sequenceWithModification, newPeptideSequence);
		Assert.assertEquals("TCQSWSSMTPHWHRRTTEYYPNGGLTR", originalPeptideSequence);
	}

}
