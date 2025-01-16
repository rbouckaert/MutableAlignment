package mutablealignment;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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
		for (int i = 0; i < newValues.length; i++) {
			sitePatterns[i][taxonNr] = newValues[i];
		}
	}

	public int [] getSiteValueByTaxon(int taxonNr) {
		return getPattern(taxonNr);
	}

	
	/**
	 * Set characters in the alignment at a specific site for all taxa
	 */
	public void setSiteValuesBySite(int siteNr, int [] newValues) {
		startEditing(null);
		System.arraycopy(newValues, 0, sitePatterns[siteNr], 0, sitePatterns[siteNr].length);
	}

	public int [] getSiteValuesBySite(int siteNr) {
		return getPattern(siteNr);
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
	
	enum EditType {singleSite, allTaxa, allSites}
	
	/** 
	 * class for tracking edits to the alignment 
	 * used to restore an alignment if necessary
	 **/
	class Edit {
		EditType type;
		int siteNr;
		int taxonNr;
		Object oldValue;
		
		Edit(int siteNr, int taxonNr, int oldValue) {
			this.type = EditType.singleSite;
			this.siteNr = siteNr;
			this.taxonNr = taxonNr;
			this.oldValue = oldValue;
		}

		Edit(int siteNr, int [] oldValue) {
			this.type = EditType.allSites;
			this.siteNr = siteNr;
			this.oldValue = oldValue;
		}

		Edit(int [] oldValue, int taxonNr) {
			this.type = EditType.allTaxa;
			this.taxonNr = taxonNr;
			this.oldValue = oldValue;
		}
		
		void undo() {
			switch(type) {
			case singleSite:
				sitePatterns[siteNr][taxonNr] = (int) oldValue;
				break;
			case allSites:;
			int [] oldValues = (int[]) oldValue;
			for (int i = 0; i < oldValues.length; i++) {
				sitePatterns[i][taxonNr] = oldValues[i];
			}
			break;
			case allTaxa:;
			System.arraycopy((int[])oldValue, 0, sitePatterns[siteNr], 0, sitePatterns[siteNr].length);
			break;
			}
		}
	}
	
	protected List<Edit> editList = new ArrayList<>();
	
	@Override
	protected void store() {
		editList.clear();
	}
	
	@Override
	public void restore() {
		for (int i = editList.size()-1; i>=0; i--) {
			editList.get(i).undo();
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

}
