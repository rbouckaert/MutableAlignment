# MutableAlignment

This is a package for [BEAST 2](http://beast2.org) for for alignments that can change during MCMC.


## Installing

TODO

## Build from code

An alternative way to install is to build from the source code. 
Frist, get code for beast2, BeastFX and MutableAlignment. Then run

```
ant install
```

to install the package.



## API

To use it, replace the standard `Alignment` with `mutablealignment.MutableAlignment` and `TreeLikelihood` with `mutablealignment.MATreeLikelihood`.

The `MutableAlignment` can be changed by operators using:

* `setSiteValue()` to set a single character in the alignment,
* `setSiteValuesByTaxon()` to set a sequence for a taxon,
* `setSiteValuesBySite()` to set site values for all taxa (but only a single site), and
* `setSiteValues()` to set the whole alignment.

There are some test classes that suggest the `MATreeLikelihood` should update itself correctly.

Note that `setSiteValue()` and `setSiteValuesByTaxon()` are inefficiently implemented: partials will be recalculated for all sites. A future BEAGLE aware implementation will probably suffer from the same limitation, just because of the BEAGLE API not allowing efficient updates of single sites.

## Paper

TODO

## Data

Data used in the paper:

TODO