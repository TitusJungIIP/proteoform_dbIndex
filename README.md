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

# How to use it

This module contains some test classes with some examples using test data. You can find it on the */src/test/java* folder, at the *ProteoformDBIndexTest.java* java file.  
