//
// This file was generated by the JavaTM Architecture for XML Binding(JAXB) Reference Implementation, vhudson-jaxb-ri-2.2-7 
// See <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Any modifications to this file will be lost upon recompilation of the source schema. 
// Generated on: 2019.02.27 at 02:24:29 PM PST 
//


package edu.scripps.yates.utilities.annotations.uniprot.xml;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlType;


/**
 * Describes the source of the evidence, when it is not assigned by UniProt, but imported from an external database.
 * 
 * <p>Java class for importedFromType complex type.
 * 
 * <p>The following schema fragment specifies the expected content contained within this class.
 * 
 * <pre>
 * &lt;complexType name="importedFromType">
 *   &lt;complexContent>
 *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       &lt;sequence>
 *         &lt;element name="dbReference" type="{http://uniprot.org/uniprot}dbReferenceType"/>
 *       &lt;/sequence>
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "importedFromType", propOrder = {
    "dbReference"
})
public class ImportedFromType {

    @XmlElement(required = true)
    protected DbReferenceType dbReference;

    /**
     * Gets the value of the dbReference property.
     * 
     * @return
     *     possible object is
     *     {@link DbReferenceType }
     *     
     */
    public DbReferenceType getDbReference() {
        return dbReference;
    }

    /**
     * Sets the value of the dbReference property.
     * 
     * @param value
     *     allowed object is
     *     {@link DbReferenceType }
     *     
     */
    public void setDbReference(DbReferenceType value) {
        this.dbReference = value;
    }

}