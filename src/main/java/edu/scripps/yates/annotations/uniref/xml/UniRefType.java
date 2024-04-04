//
// This file was generated by the JavaTM Architecture for XML Binding(JAXB) Reference Implementation, vhudson-jaxb-ri-2.2-7 
// See <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Any modifications to this file will be lost upon recompilation of the source schema. 
// Generated on: 2015.07.08 at 08:56:44 AM PDT 
//


package edu.scripps.yates.annotations.uniref.xml;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java class for UniRefType complex type.
 * 
 * <p>The following schema fragment specifies the expected content contained within this class.
 * 
 * <pre>
 * &lt;complexType name="UniRefType">
 *   &lt;complexContent>
 *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       &lt;choice>
 *         &lt;element ref="{http://uniprot.org/uniref}UniRef100"/>
 *         &lt;element ref="{http://uniprot.org/uniref}UniRef90"/>
 *         &lt;element ref="{http://uniprot.org/uniref}UniRef50"/>
 *         &lt;element ref="{http://uniprot.org/uniref}UniRef"/>
 *       &lt;/choice>
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "UniRefType", propOrder = {
    "uniRef100",
    "uniRef90",
    "uniRef50",
    "uniRef"
})
public class UniRefType {

    @XmlElement(name = "UniRef100")
    protected UniRef100 uniRef100;
    @XmlElement(name = "UniRef90")
    protected UniRef90 uniRef90;
    @XmlElement(name = "UniRef50")
    protected UniRef50 uniRef50;
    @XmlElement(name = "UniRef")
    protected UniRef uniRef;

    /**
     * Gets the value of the uniRef100 property.
     * 
     * @return
     *     possible object is
     *     {@link UniRef100 }
     *     
     */
    public UniRef100 getUniRef100() {
        return uniRef100;
    }

    /**
     * Sets the value of the uniRef100 property.
     * 
     * @param value
     *     allowed object is
     *     {@link UniRef100 }
     *     
     */
    public void setUniRef100(UniRef100 value) {
        this.uniRef100 = value;
    }

    /**
     * Gets the value of the uniRef90 property.
     * 
     * @return
     *     possible object is
     *     {@link UniRef90 }
     *     
     */
    public UniRef90 getUniRef90() {
        return uniRef90;
    }

    /**
     * Sets the value of the uniRef90 property.
     * 
     * @param value
     *     allowed object is
     *     {@link UniRef90 }
     *     
     */
    public void setUniRef90(UniRef90 value) {
        this.uniRef90 = value;
    }

    /**
     * Gets the value of the uniRef50 property.
     * 
     * @return
     *     possible object is
     *     {@link UniRef50 }
     *     
     */
    public UniRef50 getUniRef50() {
        return uniRef50;
    }

    /**
     * Sets the value of the uniRef50 property.
     * 
     * @param value
     *     allowed object is
     *     {@link UniRef50 }
     *     
     */
    public void setUniRef50(UniRef50 value) {
        this.uniRef50 = value;
    }

    /**
     * Gets the value of the uniRef property.
     * 
     * @return
     *     possible object is
     *     {@link UniRef }
     *     
     */
    public UniRef getUniRef() {
        return uniRef;
    }

    /**
     * Sets the value of the uniRef property.
     * 
     * @param value
     *     allowed object is
     *     {@link UniRef }
     *     
     */
    public void setUniRef(UniRef value) {
        this.uniRef = value;
    }

}
