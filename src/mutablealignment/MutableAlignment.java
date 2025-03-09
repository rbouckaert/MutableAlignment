package mutablealignment;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import beast.base.core.Description;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;

@Description("Alignment that can be sampled by MCMC")
public class MutableAlignment extends Alignment {
	
	/**
	 * Set 1 character in the alignment at a specific taxon and site
	 */
	public void setSiteValue(int taxonNr, int siteNr, int newValue) {
		startEditing(null);
		editList.add(new Edit(siteNr, taxonNr, sitePatterns[siteNr][taxonNr], newValue));
		sitePatterns[siteNr][taxonNr] = newValue;
	}

	public int getSiteValue(int taxonNr, int siteNr) {
		return sitePatterns[siteNr][taxonNr];
	}


	/**
	 * Set characters in the alignment at a specific taxon for all sites
	 */
	public void setSiteValuesByTaxon(int taxonNr, int [] newValues) {
		startEditing(null);
		int [] oldValues = new int[newValues.length];
		for (int i = 0; i < newValues.length; i++) {
			oldValues[i] = sitePatterns[i][taxonNr]; 
		}
		editList.add(new Edit(oldValues, taxonNr, newValues));
		for (int i = 0; i < newValues.length; i++) {
			sitePatterns[i][taxonNr] = newValues[i];
		}
	}

	public int [] getSiteValuesByTaxon(int taxonNr) {
		int [] seq = new int[sitePatterns.length];
		for (int i = 0; i < seq.length; i++) {
			seq[i] = sitePatterns[i][taxonNr];
		}
		return seq;
	}

	
	/**
	 * Set characters in the alignment at a specific site for all taxa
	 */
	public void setSiteValuesBySite(int siteNr, int [] newValues) {
		startEditing(null);
		editList.add(new Edit(siteNr, sitePatterns[siteNr].clone(), newValues));
		System.arraycopy(newValues, 0, sitePatterns[siteNr], 0, sitePatterns[siteNr].length);
	}

	public int [] getSiteValuesBySite(int siteNr) {
		return getPattern(siteNr);
	}


	/**
	 * Set all characters in the alignment 
	 */
	public void setSiteValues(int [][] newValues) {
		startEditing(null);
		int n = newValues[0].length;
		int [][] oldValues = new int[newValues.length][n];
		for (int i = 0; i < oldValues.length; i++) {
			System.arraycopy(sitePatterns[i], 0, oldValues[i], 0, n);
		}
		editList.add(new Edit(oldValues, newValues));
		for (int i = 0; i < newValues.length; i++) {
			System.arraycopy(newValues[i], 0, sitePatterns[i], 0, sitePatterns[i].length);
		}
	}
	
	
	
	/** reset site patterns to previous values **/
	public void resetSitePatterns(int siteNr, int taxonNr, int oldValue) {
		sitePatterns[siteNr][taxonNr] = oldValue;
	}

	public void resetSitePatterns(int siteNr, int[] oldValues) {
		System.arraycopy(oldValues, 0, sitePatterns[siteNr], 0, sitePatterns[siteNr].length);
	}

	public void resetSitePatterns(int[][] oldValues) {
		for (int i = 0; i < oldValues.length; i++) {
			System.arraycopy(oldValues[i], 0, sitePatterns[i], 0, sitePatterns[i].length);
		}
	}
	
	
	
	/** Loggable Interface **/
	@Override
	public void init(PrintStream out) {
		for (int i = 0; i < taxaNames.size(); i++) {
			out.print(taxaNames.get(i) + "\t");
		}
	}

	@Override
	public void log(long sample, PrintStream out) {
		for (int i = 0; i < taxaNames.size(); i++) {
			for (int j = 0; j < sitePatterns.length; j++) {
				out.print(getDataType().getCharacter(sitePatterns[j][i]));
			}
			out.print("\t");
		}
	}
	
	@Override
	public void close(PrintStream out) {
	}
	

	/** StateNode stuff **/
	
	enum EditType {singleSite, allTaxa, allSites, all}
	
	protected List<Edit> editList = new ArrayList<>();
	
	@Override
	protected void store() {
		editList.clear();
	}
	
	@Override
	public void restore() {
		for (int i = editList.size()-1; i>=0; i--) {
			editList.get(i).undo(this);
		}
		editList.clear();
	}
	
	/**
	 * calculate patterns from sequence data *
	 */
	protected void calcPatterns(boolean log) {
		int taxonCount = counts.size();
		int siteCount = counts.get(0).size();

		// convert data to transposed int array
		sitePatterns = new int[siteCount][taxonCount];
		for (int i = 0; i < taxonCount; i++) {
			List<Integer> sites = counts.get(i);
			for (int j = 0; j < siteCount; j++) {
				sitePatterns[j][i] = sites.get(j);
			}
		}
		int patterns = siteCount;
		
		// reserve memory for patterns
		patternWeight = new int[patterns];
		Arrays.fill(patternWeight, 1);

		// find patterns for the sites
		patternIndex = new int[siteCount];
		for (int i = 0; i < siteCount; i++) {
			patternIndex[i] = i;
		}

		if (siteWeights != null) {
			Arrays.fill(patternWeight, 0);
			for (int i = 0; i < siteCount; i++) {
				patternWeight[patternIndex[i]] += siteWeights[i];
			}
		}

		// determine maximum state count
		// Usually, the state count is equal for all sites,
		// though for SnAP analysis, this is typically not the case.
		maxStateCount = 0;
		for (int m_nStateCount1 : stateCounts) {
			maxStateCount = Math.max(maxStateCount, m_nStateCount1);
		}
		// report some statistics
		if (log && taxaNames.size() < 30) {
			for (int i = 0; i < taxaNames.size(); i++) {
				Log.info.println(taxaNames.get(i) + ": " + counts.get(i).size() + " " + stateCounts.get(i));
			}
		}

		if (stripInvariantSitesInput.get()) {
			// don't add patterns that are invariant, e.g. all gaps
			if (log)
				Log.info.println("Stripping invariant sites");

			int removedSites = 0;
			for (int i = 0; i < patterns; i++) {
				int[] pattern = sitePatterns[i];
				int value = pattern[0];
				boolean isInvariant = true;
				for (int k = 1; k < pattern.length; k++) {
					if (pattern[k] != value) {
						isInvariant = false;
						break;
					}
				}
				if (isInvariant) {
					removedSites += patternWeight[i];
					patternWeight[i] = 0;

					if (log)
						Log.info.print(" <" + value + "> ");
				}
			}
			if (log)
				Log.info.println(" removed " + removedSites + " sites ");
		}
	} // calcPatterns

	public Integer [] getDirtySequenceIndices() {
		Set<Integer> dirtySequences = new HashSet<>();
		for (Edit edit : editList) {
			switch (edit.type) {
			case all:
			case allSites:
				for (int i = 0; i < getTaxonCount(); i++) {
					dirtySequences.add(i);
				}
				break;
			case allTaxa:
			case singleSite:
				dirtySequences.add(edit.taxonNr);
				break;
			}
		}
		return dirtySequences.toArray(new Integer[] {});
	}

}
