package edu.scripps.yates.utilities.sequence;

import java.io.IOException;
import java.util.*;
import java.util.logging.Level;

import com.compomics.util.protein.Enzyme;

import gnu.trove.map.hash.THashMap;

import static edu.scripps.yates.utilities.masses.AssignMass.H2O_PROTON;

public class MyEnzyme extends Enzyme {
	private final Map<String, Collection<String>> cleavagesBySequence = new THashMap<String, Collection<String>>();
	private boolean cacheEnabled = false;

	public MyEnzyme(String aTitle, String aCleavage, String aRestrict, String aPosition, int aMiscleavages) {
		this(aTitle, aCleavage, aRestrict, aPosition, aMiscleavages, false);
	}

	public MyEnzyme(String aTitle, String aCleavage, String aRestrict, String aPosition, int aMiscleavages,
			boolean cacheEnabled) {
		super(aTitle, aCleavage, aRestrict, aPosition, aMiscleavages);
		this.cacheEnabled = cacheEnabled;
	}

	public MyEnzyme(com.compomics.util.experiment.biology.Enzyme enzyme, int maxMissedCleavages) {
		this(enzyme, maxMissedCleavages, false);
	}

	public MyEnzyme(com.compomics.util.experiment.biology.Enzyme enzyme, int maxMissedCleavages, boolean cacheEnabled) {
		super(enzyme, maxMissedCleavages);
		this.cacheEnabled = cacheEnabled;
	}

	public MyEnzyme(String aTitle, String aCleavage, String aRestrict, String aPosition) {
		this(aTitle, aCleavage, aRestrict, aPosition, false);
	}

	public MyEnzyme(String aTitle, String aCleavage, String aRestrict, String aPosition, boolean cacheEnabled) {
		super(aTitle, aCleavage, aRestrict, aPosition);
		this.cacheEnabled = cacheEnabled;
	}

	public boolean isCacheEnabled() {
		return cacheEnabled;
	}

	public void setCacheEnabled(boolean cacheEnabled) {
		this.cacheEnabled = cacheEnabled;
	}

	public void clearCache() {
		cleavagesBySequence.clear();
	}

	/**
	 * this method is the focus of the enzyme instance. it can perform an
	 * <i>in-silico</i> digest of a protein sequence according to the specifications
	 * detailed in the construction or via the setters. only returns peptides
	 * between the minimum and maximum peptide lengths.
	 *
	 * @param aprotein         protein instance to cleave.
	 * @param minpeptidelength the minimum peptide length to consider
	 * @param maxpeptidelength the maximum peptide length to consider
	 * @return protein[] with the resultant peptides.
	 */

	public Collection<String>  cutSeq(String aProtein, int minPeptideLength, int maxPeptideLength,
						Collection<String> result, boolean semicleavage)  {

		//Enzyme enz = sparam.getEnzyme();
//
		//AssignMass aMass = AssignMass.getInstance(true);
		if (cacheEnabled && cleavagesBySequence.containsKey(aProtein)) {
			return cleavagesBySequence.get(aProtein);
		}
		if (result == null) {
			result = new ArrayList<String>();
		}
		final int length = aProtein.length();

		final char[] pepSeq = new char[maxPeptideLength]; //max seq length
		int curSeqI = 0;

		final int maxMissedCleavages = this.getMiscleavages();
		int maxIntCleavage = getMiscleavages();


//	    System.out.println(fasta.getSequestLikeAccession());
			//System.out.println(fasta.getDefline());
			int count =1;
			for (int start = 0; start < length; ++start) {
				int end = start;

				//reset the preallocated seq byte array
				//Arrays.fill(seq, 0, curSeqI > 0?curSeqI-1:0, (byte) 0); //no need, we copy up to curSeqI nowu
				curSeqI = 0;

				//float precMass = Constants.H2O_PROTON_SCALED_DOWN;



				// System.out.println("===>>" + precMass + "\t" + Constants.MAX_PRECURSOR);
				//System.out.println("==" + j + " " + length + " " + (j < length));

				//int testC=0;
				int pepSize = 0;

				int intMisCleavageCount = 0;



				//while (precMass <= Constants.MAX_PRECURSOR_MASS && end < length) {
				while ( end < length && curSeqI<maxPeptideLength)
				{
					pepSize++;

					final char curIon = aProtein.charAt(end);
					pepSeq[curSeqI++] = curIon;


					//if (precMass > Constants.MAX_PRECURSOR_MASS) {

					// System.out.println(""+precMass);
					if (iCleavables.get(aProtein.charAt(end))!= null) {
						intMisCleavageCount++;
					}

					//       System.out.println("--->>" + String.valueOf(Arrays.copyOf(pepSeq, curSeqI)) + "\t" + intMisCleavageCount + "\t" + maxIntCleavage);
					if ((intMisCleavageCount > maxIntCleavage && maxIntCleavage !=-1) && (iCleavables.get(aProtein.charAt(end))== null)) {
						break;
					}

					if (pepSize >= minPeptideLength  && pepSize<=maxPeptideLength) {
						//Constants.MIN_PRECURSOR ) {


						final boolean cleavageStatus = checkCleavage(aProtein, start, end, semicleavage);

						if (cleavageStatus ) {
							//if (cleavageStatus == 2) {
							//qualifies based on params

							final String peptideSeqString = String.valueOf(Arrays.copyOf(pepSeq, curSeqI));
							result.add(peptideSeqString);
						}
					}


					if ((intMisCleavageCount > maxIntCleavage && maxIntCleavage !=-1) ) {
						break;
					}

					++end;

				}
//
		}
		return result;

	}



	public Collection<String> cleave(String aProtein, int minPeptideLength, int maxPeptideLength,
			Collection<String> result) {
		if (cacheEnabled && cleavagesBySequence.containsKey(aProtein)) {
			return cleavagesBySequence.get(aProtein);
		}
		if (result == null) {
			result = new ArrayList<String>();
		}
		// We'll need a lot of stuff here.
		// - a Vector for all the startindices
		// - a Vector for the stopindices
		// - a Vector of intermediate results.
		final Vector startIndices = new Vector(20, 10);
		final Vector endIndices = new Vector(20, 10);
		final Vector interMed = new Vector(20, 10);

		// We will also feed the current Protein sequence into a
		// char[] for easy iteration.
		final char[] sequence = aProtein.toCharArray();

		// Check for a header that contains locations.
		final int headerStart = 0;

		// Okay, I guess we've set the stage now.
		// Let's start cleaving!
		int walkingIndex = 0;

		for (int i = 0; i < sequence.length; i++) {
			// Transform the current char into the corresponding wrapper.
			final Character current = Character.valueOf(sequence[i]);

			// See whether it is a cleavable residu!
			if (iCleavables.get(current) != null) {
				// Okay, this should be cleavable.
				// First of all however, we need to check
				// for the possible presence of a restrictor!
				// (And, of course, first check to see whether there is a
				// next character at all!)
				if ((i + 1) < sequence.length) {
					final Character next = Character.valueOf(sequence[i + 1]);
					if (iRestrictors.get(next) != null) {
						// It is a restrictor!
						// Just let the loop continue!
						continue;
					}
				}

				// Since we've gotten to here, we need to cleave here!
				// So do it!
				// Oh yeah, and mind the position of cleaving!
				String temp = null;
				int start = -1;
				int end = -1;

				if (getPosition() == Enzyme.CTERM) {
					// Take the part, starting from walkingIndex up to the
					// current
					// as a new peptide and store it in the interMed Vector.
					temp = new String(sequence, walkingIndex, ((i - walkingIndex) + 1));
					// Start index is human-readable (starting from 1),
					// hence the '+1'.
					start = headerStart + walkingIndex + 1;
					end = headerStart + i + 1;
					// Start the next peptide after the current one.
					// An index that so happens to
					walkingIndex = i + 1;
				} else if (getPosition() == Enzyme.NTERM) {
					temp = new String(sequence, walkingIndex, (i - walkingIndex));
					// Start index is human readable: starting from 1.
					start = headerStart + walkingIndex + 1;
					end = headerStart + i;
					walkingIndex = i;
				}

				// Add each retrieved value to the correct
				// Vector.
				interMed.add(temp);
				startIndices.add(Integer.valueOf(start));
				endIndices.add(Integer.valueOf(end));
			}
		}

		// Add this point, we should check whether we have
		// the entire sequence.
		// We probably don't, because the last cleavable residu will
		// probably not have been the last residu in the sequence.
		// That's why we should append the 'remainder' of our cleavage
		// as well (and the corresponding indices as well, of course).
		if (walkingIndex < sequence.length) {
			interMed.add(new String(sequence, walkingIndex, (sequence.length - walkingIndex)));
			startIndices.add(Integer.valueOf(headerStart + walkingIndex + 1));
			endIndices.add(Integer.valueOf(headerStart + sequence.length));
		}

		// Allright, now we should have all the individual peptides.
		// Now we should take into account the specified number of miscleavages.

		// Get all the sequences up to now.
		final String[] imSequences = (String[]) interMed.toArray(new String[interMed.size()]);

		// Cycle the current sequences.
		for (int j = 0; j < imSequences.length; j++) {

			String temp = imSequences[j];

			// Apply the number of allowed missed cleavages sequentially from
			// this sequence.
			for (int k = 0; k < getMiscleavages(); k++) {

				// If we fall outside of the range of current sequences
				// (for instance if we try to apply a second allowed missed
				// cleavage to the penultimate peptide, we fall outside of
				// the available peptides!)
				// we break the loop.
				if ((j + k + 1) >= imSequences.length) {
					break;
				}

				// Add our constructed sequence.
				temp += imSequences[j + k + 1];
				interMed.add(temp);
				startIndices.add(startIndices.get(j));
				endIndices.add(endIndices.get(j + k + 1));
			}
		}

		// Cycle all to check for

		// We've got all sequences.
		// Let's construct the Protein instances for them and
		// then return them!
		final int liSize = interMed.size();

		// Create the Proteins and store them.
		for (int i = 0; i < liSize; i++) {

			// If the sequence comes from a translation, it will contain an '_'
			// if a stopcodon is present.
			// Omit all sequences containing these.
			final String pepSequence = (String) interMed.get(i);
			if (pepSequence.indexOf("_") < 0) {

				// only include peptides within the min and max peptide lengths
				if (pepSequence.length() >= minPeptideLength && pepSequence.length() <= maxPeptideLength) {

					result.add(pepSequence);
				}
			}
		}
		if (cacheEnabled) {
			cleavagesBySequence.put(aProtein, result);
		}
		return result;
	}
	private  boolean checkCleavage(String proteinSequence, int start, int end, boolean semiCleave) {

		int numCleavages  = 0;
		if (start > 0) {
			// that is the case when is fully triptic but the preAA is not a
			// cleave site
			final char preAA = proteinSequence.charAt(start - 1);
			if (iCleavables.get(preAA)!=null) {
				numCleavages++;
			}
		} else {
			numCleavages++;
			// if the peptide starts at the N-terminal of the protein
		}

		// if the last AA is not a cleavage site
		final char lastAA = proteinSequence.charAt(end);
		if (iCleavables.get(lastAA)!=null || proteinSequence.length() == end + 1) {
			// if the lastAA of the peptide is not the Cterminal of the protein
			// and is fully triptic
			numCleavages++;
		}
		int minCleavage = semiCleave ? 1:2;
		if(numCleavages<minCleavage)
		{
			return false;
		}



		// check the postAA to see if is in the enzymeNoCutAA
		if ( proteinSequence.length() > end + 1) {
			final Character postAA = proteinSequence.substring(end + 1, end + 2).charAt(0);
			if (iRestrictors.get(postAA)!=null) {
				return false;
			}
		}
		int numMissedCleavages = 0;
		final String peptideSequence = proteinSequence.substring(start, end + 1);
		for (int index = 0; index < peptideSequence.length() - 1; index++) {
			char aa = peptideSequence.charAt(index);
			if (iCleavables.get(aa)!=null) {
				numMissedCleavages++;
			}
		}
		if (numMissedCleavages > this.getMiscleavages() ) {
			return false;
		}


		return true;
	}


}
