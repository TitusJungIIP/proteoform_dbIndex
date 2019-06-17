# proteoform_dbIndex
This is an extension of dbIndex module to index proteoforms, getting both sequence variances and PTMs from UniprotKB annotations.
Given a certain FASTA file with Uniprot accessions, it automatically downloads the UniprotKB annotations and applies them to the proteins, generating all posible combinations of annotations for the in-silico digested peptides and indexing them.  
Therefore, this index code could be used by any algorithm or search engine that gets experimental parent masses and wants to get the peptide candidates in the proteome, taking into account that it could contain a PTM or a sequence variance that is not included in the stardard FASTA file.

## How to get it
Using maven, add this to your *pom.xml* file:
 * Artifact (version may change over time, so we recommend you to look into the maven artifact URLs to get the latest version):
```
<dependency>
  <groupId>edu.scripps.yates</groupId>
  <artifactId>proteoform_dbindex</artifactId>
  <version>0.0.2-SNAPSHOT</version>
</dependency>
```
 * Repositories:
```
<repositories>
  <repository>
    <id>internal</id>
    <name>John Yates's lab releases maven repository</name>
    <url>http://sealion.scripps.edu/archiva/repository/internal/</url>
  </repository>
  <repository>
    <id>snapshots</id>
    <name>John Yates's lab snapshots maven repository</name>
    <url>http://sealion.scripps.edu/archiva/repository/snapshots/</url>
  </repository>
</repositories>
```

# How it works:

This module contains some test classes with some examples using test data. You can find it on the */src/test/java* folder, at the [*ProteoformDBIndexTest.java*](https://github.com/proteomicsyates/proteoform_dbIndex/blob/master/src/test/java/proteoform_dbindex/ProteoformDBIndexTest.java) java file.  

Firstly, you may know the steps this indexing is doing for every FASTA file:
 * **Reading input parameters**: The input parameters file is one of the parameters file used by BlazzMass search engine. However, many of the parameters included on it, are ignored, and so the only ones used in the indexing are (`with some example values`):
   - *database_name*: Full path to FASTA file. `D:\\Salva\\git_projects\\proteoform_dbindex\\src\\test\\resources\\Q8WZ42.fasta`
   - *enzyme_residues*: Amino acids in which the enzyme should cut. `KR`.
   - *max_num_internal_cleavage_sites*: Number of allowed missed-cleavages. `3`.
   - *miscleavage*: This parameter name is a little bit confusing. It means 'semi-tryptic' or 'semi-cleavage'. If it is 'false', semi-cleavages are not allowed, and so in case of trypsin, it will allow only fully-tryptic peptides. In case of being 'true', semi-cleavages are allowed.
   - *enzyme_nocut_residues*: amino acid that if present before the cleavage amino acid, will make the cleavage to not ocurr.
 * **Reading FASTA file**: It reads fasta file and extract UniprotKB accession numbers. Note that if the fasta file doesn't have accession numbers from UniprotKB, it will not be able to get the annotations of the sequence variants or PTMs.
 * **Getting UniprotKB annotations**: It retrieves the UniprotKB entries of the proteins in the FASTA file. It will store them in the folder you specify. By default, it will check for the latest release of the database. However, you can specify to use the same release even though a new version has been released.
 * **Parsing protein sequences**: It reads the protein sequences and performs an in-silico digestion of the proteins, taking into account the parameters specified in the parameters file. While doing it, it gets the protein sequence variances and PTMs described in UniprotKB for each protein, and maps them to the sequence so that the generated peptides may contain these annotations.
 * **Peptide sequences are indexed**: Peptide sequences and their possible annotations are indexed into a SQLite database so that then, they can be retrieved quickly searching for a parent mass.
 
 # How to use it:
 
 As we have mentioned, you have to use a parameters file. You can find one [here](https://github.com/proteomicsyates/proteoform_dbIndex/blob/master/src/test/resources/blazmass_Q13523.params).  
   
 Then, you can see this code snippet where an instance of ***ProteoformDBIndexInterface***: 
 ```
 // input parameters file
 File paramFile = new ClassPathResource("blazmass_Q13523.params").getFile();
 
 // maximum number of variations (sequence variations and PTMs) per peptide
 int maxNumVariationsPerPeptide = 4;
 
 // Uniprot repository version release
 // null for latest version or "2019_05" for May 2019 version, for example
 String uniprotVersion = "2019_05";
	
 // Uniprot annotations retriever. It will retrieve the annotations to folder uniprotReleasesFolder
 UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, true);
 
 // Create ProteoformDBIndexInterface instance
 ProteoformDBIndexInterface proteoformDBIndex = new ProteoformDBIndexInterface(paramFile, uplr, uniprotVersion, maxNumVariationsPerPeptide);
 ```  
   
 Once having the *proteoformDBIndex* we can ask for the proteins containing a certain peptide sequence:  
   
 ```
 System.out.println("Looking for the proteins from peptide " + "SPSPDDILERVAADVKEYER");
 final Set<IndexedProtein> proteins = proteoformDBIndex.getProteins("SPSPDDILERVAADVKEYER");
 for (final IndexedProtein indexedProtein : proteins) {
    System.out.println(indexedProtein.getAccession());
    System.out.println(indexedProtein.getFastaDefLine());
 }
 ```
   
 Output:  
   
 ```
 Looking for the proteins from peptide SPSPDDILERVAADVKEYER
 Q13523
 >sp|Q13523|PRP4B_HUMAN Serine/threonine-protein kinase PRP4 homolog OS=Homo sapiens OX=9606 GN=PRPF4B PE=1 SV=3
 ```
   
Protein *[Q13523](https://www.uniprot.org/uniprot/Q13523)* contains peptide *SPSPDDILERVAADVKEYER* starting at position 578.  
  
Now, as we can check in UniprotKB, protein (Q13523)[https://www.uniprot.org/uniprot/Q13523] contains a **polymorfism at position [584](https://www.uniprot.org/blast/?about=Q13523[584]&key=Natural%20variant&id=VAR_047798) ([VAR_047798](https://web.expasy.org/variant_pages/VAR_047798.html))**, so we can ask for that modified peptide SPSPDD***V***LERVAADVKEYER containing a ***'V'*** instead of a ***'I'*** (***I → V***):
  
```
 System.out.println("Looking for the proteins from peptide with natural variance " + "SPSPDDVLERVAADVKEYER");
 final Set<IndexedProtein> proteins = proteoformDBIndex.getProteins("SPSPDDVLERVAADVKEYER");
 for (final IndexedProtein indexedProtein : proteins) {
    System.out.println(indexedProtein.getAccession());
    System.out.println(indexedProtein.getFastaDefLine());
 }
 ```  
   
 Output:  
   
 ```
 Looking for the proteins from peptide with natural variance SPSPDDVLERVAADVKEYER
 Q13523
 >sp|Q13523|PRP4B_HUMAN Serine/threonine-protein kinase PRP4 homolog OS=Homo sapiens OX=9606 GN=PRPF4B PE=1 SV=3
 ```
As you can see, we found that the peptide containing that variant is also found for that protein.    
  
We can also use the index with masses instead of sequences. Therefore, we can ask for the mass of that peptide with the natural variance and get an *IndexedSequence* object:  
  
```
double parentMass = IndexUtil.calculateMass("SPSPDDVLERVAADVKEYER");
System.out.println("Looking for the mass " + parentMass + " that is from peptide with natural variance " + "SPSPDDVLERVAADVKEYER");
List<IndexedSequence> sequences = proteoformDBIndex.getSequences(parentMass, 0.0001);
for (final IndexedSequence indexedSequence : sequences) {
  System.out.println(indexedSequence.getSequence() + "\t"
	+ IndexUtil.calculateMass(indexedSequence.getSequence()) + "\t"
	+ indexedSequence.getModSequence() + "\t" + indexedSequence.getMass());
}
```
  
Output:
  
```
Looking for the mass 2275.1242204659998 that is from peptide SPSPDDVLERVAADVKEYER
SPSPDDVLERVAADVKEYER	2275.1242204659998	SPSPDD[I->V]LERVAADVKEYER	2275.1242204659993
```  
  
As you can see, the returned *IndexedSequence* contains the annotation of the sequence variation (**SPSPDD\[I->V\]LERVAADVKEYER**) when calling to *.getModSequence()* method.  
  
Then, we also know from UniprotKB that this protein has been annotated as having a phosphorilation at positions [578](https://www.uniprot.org/blast/?about=Q13523[578]&key=Modified%20residue) and [580](https://www.uniprot.org/blast/?about=Q13523[580]&key=Modified%20residue) and so, we can ask for the same parent mass plus a phosphorilation (+79.966):
  
```
final double phosphorilatedParentMass = parentMass + 79.966331;
System.out.println("Looking for the mass " + phosphorilatedParentMass + " that is from peptide with natural variance "
		+ "SPSPDDVLERVAADVKEYER" + " plus a phosphorilation");
final List<IndexedSequence> phosphorilatedSequences = proteoformDBIndex.getSequences(phosphorilatedParentMass, 0.0001);
for (final IndexedSequence indexedSequence : phosphorilatedSequences) {
	if (indexedSequence instanceof IndexedSequenceWithPTMs) {
		final IndexedSequenceWithPTMs indexedSequenceWithPTM = (IndexedSequenceWithPTMs) indexedSequence;
		System.out.println(indexedSequenceWithPTM.getSequence() + "\t"
			+ IndexUtil.calculateMass(indexedSequenceWithPTM.getSequence()) + "\t"
			+ indexedSequenceWithPTM.getModSequence() + "\t" + indexedSequenceWithPTM.getMass() + "\t"
			+ phosphorilatedParentMass);
	} else {
		System.out.println(indexedSequence.getSequence() + IndexUtil.calculateMass(indexedSequence.getSequence())
			+ "\t" + indexedSequence.getModSequence() + "\t" + +indexedSequence.getMass()
			+ "\t" + phosphorilatedParentMass);
	}
}
```  

This will output:  
  
```
Looking for the mass 2355.090551466 that is from peptide SPSPDDVLERVAADVKEYER plus a phosphorilation
SPSPDDVLERVAADVKEYER	2275.1242204659998	SPS[+79.9663]PDD[I->V]LERVAADVKEYER	2355.0905204659994	2355.090551466
SPSPDDVLERVAADVKEYER	2275.1242204659998	S[+79.9663]PSPDD[I->V]LERVAADVKEYER	2355.0905204659994	2355.090551466
```  

As you can see, the index returns two different *IndexedSequence* objects that correspond to the protein [Q13523](https://www.uniprot.org/uniprot/Q13523) with a phosphorilation at position [578](https://www.uniprot.org/blast/?about=Q13523[578]&key=Modified%20residue) and [580](https://www.uniprot.org/blast/?about=Q13523[580]&key=Modified%20residue).  
  
Then, if we ask for the mass corresponding to the doubly phosphorilated mass:
  
```
final double doublyPhosphorilatedParentMass = parentMass + 79.966331 + 79.966331;
System.out.println("Looking for the mass " + doublyPhosphorilatedParentMass
		+ " that is from peptide with natural variance" + "SPSPDDVLERVAADVKEYER"
		+ " plus 2 phosphorilations");
List<IndexedSequence> doublePhosphorilatedSequences = proteoformDBIndex
		.getSequences(doublyPhosphorilatedParentMass, 0.0001);
for (final IndexedSequence indexedSequence : doublePhosphorilatedSequences) {
	if (indexedSequence instanceof IndexedSequenceWithPTMs) {
		final IndexedSequenceWithPTMs indexedSequenceWithPTM = (IndexedSequenceWithPTMs) indexedSequence;
		System.out.println(indexedSequenceWithPTM.getSequence() + "\t"
			+ IndexUtil.calculateMass(indexedSequenceWithPTM.getSequence()) + "\t"
			+ indexedSequenceWithPTM.getModSequence() + "\t" + indexedSequenceWithPTM.getMass() + "\t"
			+ doublyPhosphorilatedParentMass);
	} else {
		System.out.println(indexedSequence.getSequence() + IndexUtil.calculateMass(indexedSequence.getSequence())
			+ "\t" + "\t" + indexedSequence.getModSequence() + "\t" + +indexedSequence.getMass()							+ "\t" + phosphorilatedParentMass);
	}
}
```
  
The output will be:  
  
```
Looking for the mass 2435.056882466 that is from peptide with natural varianceSPSPDDVLERVAADVKEYER plus 2 phosphorilations
SPSPDDVLERVAADVKEYER	2275.1242204659998	S[+79.9663]PS[+79.9663]PDD[I->V]LERVAADVKEYER	2435.0568204659994	2435.056882466
```  
  
As you can see the doubly phosphorilated peptide is also found in the index and is correctly annotated.



