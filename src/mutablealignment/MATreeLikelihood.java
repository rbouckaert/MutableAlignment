package mutablealignment;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;

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
	public double calculateLogP() {
		if (alignmentNeedsUpdate) {
			updateAlignment();
			alignmentNeedsUpdate = false;
		}
		return super.calculateLogP();
	}
	
	
	private List<Integer> dirtySequences = new ArrayList<>();
	
	private void updateAlignment() {
        dirtySequences.clear();
        for (Integer i : alignment.getDirtySequenceIndices()) {
        	dirtySequences.add(i);
        }
        updateTipData();
        
	}
    private void updateTipData() {
        TreeInterface tree = treeInput.get();
    	int patternCount = alignment.getPatternCount();
        int stateCount = alignment.getDataType().getStateCount();
    	for (int nodeNr: dirtySequences) {
    		Node node = tree.getNode(nodeNr);
            int taxonIndex = alignment.getTaxonIndex(node.getID());

            if (m_useAmbiguities.get()) {
	            likelihoodCore.setNodePartialsForUpdate(nodeNr);
	            double[] partials = new double[patternCount * stateCount];
	            
	            int k = 0;
	            for (int patternIndex_ = 0; patternIndex_ < patternCount; patternIndex_++) {                
	                double[] tipLikelihoods = alignment.getTipLikelihoods(taxonIndex, patternIndex_);
	                if (tipLikelihoods != null) {
	                	for (int state = 0; state < stateCount; state++) {
	                		partials[k++] = tipLikelihoods[state];
	                	}
	                } else {
		                boolean[] stateSet = alignment.getStateSet(stateCount);
		                for (int state = 0; state < stateCount; state++) {
		                	 partials[k++] = (stateSet[state] ? 1.0 : 0.0);                
		                }
	                }
	            }
	            likelihoodCore.setNodePartials(nodeNr, partials);
	            
            } else {
            	
                int[] states = new int[patternCount];
                DataType dataType = alignment.getDataType();
                for (int i = 0; i < patternCount; i++) {
                    int code = alignment.getPattern(taxonIndex, i);
                    int[] statesForCode = dataType.getStatesForCode(code);
                    if (statesForCode.length==1)
                        states[i] = statesForCode[0];
                    else
                        states[i] = code; // Causes ambiguous states to be ignored.
                }
                likelihoodCore.setNodeStates(nodeNr, states);
                node.makeDirty(Tree.IS_DIRTY);
            }
        }
    }

	
	
	
	@Override
	public void restore() {
		super.restore();
		
        if (!m_useAmbiguities.get()) {
        	updateTipData();
        	dirtySequences.clear();
        }
	}
	

	@Override
	protected boolean requiresRecalculation() {
		boolean isDirty =  super.requiresRecalculation();
		if (alignment.somethingIsDirty()) {
			alignmentNeedsUpdate = true;
            hasDirt = Tree.IS_DIRTY;
			isDirty = true;
		}
		return isDirty;
	}
	
}
