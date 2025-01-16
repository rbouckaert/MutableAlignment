package mutablealignment;

import mutablealignment.MutableAlignment.EditType;

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
	
	Edit(int [][] oldValue) {
		this.type = EditType.all;
		this.oldValue = oldValue;
	}

	void undo(MutableAlignment mutableAlignment) {
		switch(type) {
		case singleSite:
			mutableAlignment.resetSitePatterns(siteNr, taxonNr, (int) oldValue);
			break;
		case allSites:;
			int [] oldValues = (int[]) oldValue;
			for (int i = 0; i < oldValues.length; i++) {
				mutableAlignment.resetSitePatterns(i, taxonNr, oldValues[i]);
			}
			break;
		case allTaxa:;
			mutableAlignment.resetSitePatterns(siteNr, (int[])oldValue);
			break;
		case all:;
			mutableAlignment.resetSitePatterns((int[][])oldValue);
			break;
		}
	}
}