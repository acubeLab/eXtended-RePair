# eXtended RePair

More than twenty years after its introduction, the RePair algorithm [1] is still widely used especially for grammar comression. This repository contains the C implementation by Gonzalo Navarro with some extra functionalities provided by the following command line options:

* `-x X` never form a pair involving a symbol $\leq X$
* `-r R` do not generate more than R rules  
* `-m M` use at most `M` Megabytes of RAM

The `M` paramater was introduced by Gonzalo and it works as follows: for an input of size $n$ the fast (linear) update strategy needs $24n$ bytes of RAM; if `M` is smaller, a quadratic update algorithm is used until the sequence is sufficiently short so that the machinery for the linear algorithm fits in `M` megabytes. 


Additional functionalities and algorithmic improvements will be added to this code base trying to maintain backward compatibility.



## Installation 

Clone/download the repository then `make`. The executables `xrepair.x` (compression) and `despair.x` (decompression) will be created in the same directory. Run without arguments to get a minimal help.



## Output format

The file `name.R` contains the size of the orginal alphabet stored in a `uint32_t` followed by pairs of `uint_32t`'s representing the rules in the format (left symbol, right symbol). The new symbol is assigned to the next available codeword so it does not neeed to be specified. 


The file `name.C` contains the sequence of terminal/non-terminal symbols when the grammar construction stops. Note that grammar construction will stop when either:

 1. no pair of symbols appears more than once (excluding pairs contaning symbols indicated by the `-x` options)

 2. the option `-r R` was used and we have already generated `R` rules

 3. the number of distinc symbols (terminal and non-terminal) reaches $2^{32}$ 


## References

[1] N.J. Larsson; A. Moffat *Offline dictionary-based compression*. IEEE Data compression conference 1999,  DOI: [10.1109/DCC.1999.755679](https://doi.org/10.1109/DCC.1999.755679) and Proceedings IEEE, Vol. 88, Issue 11, 2000, DOI: [10.1109/5.892708](https://doi.org/10.1109/5.892708)

