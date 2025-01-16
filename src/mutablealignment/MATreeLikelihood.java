package mutablealignment;

import beast.base.core.Description;
import beast.base.evolution.likelihood.TreeLikelihood;

@Description("Tree likelihood that can efficiently recalculate changes in a mutable alignment")
public class MATreeLikelihood extends TreeLikelihood {
	
	private MutableAlignment alignment;
	private boolean alignmentNeedsUpdate;
	
	@Override
	public void initAndValidate() {
		if (!(dataInput.get() instanceof MutableAlignment)) {
			throw new IllegalArgumentException("Expected MutableAlignment as data, not " + dataInput.get().getClass().getName());
		}
		alignment = (MutableAlignment) dataInput.get();
		
		boolean useJava = System.getProperty("java.only") == null ? false : Boolean.valueOf(System.getProperty("java.only"));
		System.setProperty("java.only", "true");
		super.initAndValidate();
		System.setProperty("java.only", useJava + "");
		
		alignmentNeedsUpdate = false;
	}
	
	
	@Override
	protected boolean requiresRecalculation() {
		boolean isDirty =  super.requiresRecalculation();
		if (alignment.somethingIsDirty()) {
			alignmentNeedsUpdate = true;
			isDirty = true;
		}
		return isDirty;
	}
	
}
