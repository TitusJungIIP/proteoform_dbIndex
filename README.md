# proteoform_dbIndex
This is an extension of dbIndex module to index proteoforms, getting both sequence variances and PTMs from UniprotKB annotations.
Given a certain FASTA file with Uniprot accessions, it automatically downloads the UniprotKB annotations and applies them to the proteins, generating all posible combinations of annotations for the in-silico digested peptides and indexing them.  
Therefore, this index code could be used by any algorithm or search engine that gets experimental parent masses and wants to get the peptide candidates in the proteome, taking into account that it could contain a PTM or a sequence variance that is not included in the stardard FASTA file.

## How to get it
Using maven, add this to your prom.xml file:
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

This module contains some test classes with some examples using test data. You can find it on the */src/test/java* folder, at the *ProteoformDBIndexTest.java* java file.  

Firstly, you may know the steps this indexing is doing for every FASTA file:
 * **Reading input parameters**: The input parameters file is one of the parameters file used by BlazzMass search engine. However, many of the parameters included on it, are ignored, and so the only ones used in the indexing are (`with some example values`):
   - *database_name*: Full path to FASTA file. `D:\\Salva\\git_projects\\proteoform_dbindex\\src\\test\\resources\\Q8WZ42.fasta`
   - *enzyme_residues*: Amino acids in which the enzyme should cut. `KR`.
   - *miscleavage*: Number of allowed missed cleavages. `3`.
   - *enzyme_nocut_residues*: amino acid that if present before the cleavage amino acid, will make the cleavage to not ocurr.
 * **Reading FASTA file**: It reads fasta file and extract UniprotKB accession numbers. Note that if the fasta file doesn't have accession numbers from UniprotKB, it will not be able to get the annotations of the sequence variants or PTMs.
 * **Getting UniprotKB annotations**: It retrieves the UniprotKB entries of the proteins in the FASTA file. It will store them in the folder you specify. By default, it will check for the latest release of the database. However, you can specify to use the same release even though a new version has been released.
 * **Parsing protein sequences**: It reads the protein sequences and performs an in-silico digestion of the proteins, taking into account the parameters specified in the parameters file. While doing it, it gets the protein sequence variances and PTMs described in UniprotKB for each protein, and maps them to the sequence so that the generated peptides may contain these annotations.
 * **Peptide sequences are indexed**: Peptide sequences and their possible annotations are indexed into a SQLite database so that then, they can be retrieved quickly searching for a parent mass.
 
 
   
