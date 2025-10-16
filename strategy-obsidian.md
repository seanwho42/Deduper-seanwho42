## Overview

We need to be able to read through *sorted* sam files and remove the PCR duplicates, without undue burden on the memory on the machines we are running this on. Example input and output files are included in `unit-test`.

In broad terms, we need to read through the alignments and compare with other reads at the same position and direction, and check the UMIs to see if they are PCR duplicates (or if the UMIs are erroneous).

Some context:
- Samtools sorts by reference, position in the reference, then by reverse flag per [documentation](https://www.htslib.org/doc/1.14/samtools-sort.html).
- Reads in sam files are tab separated, so we can split off of our tabs and then pull the relevant fields to get the information we need.


## The algorithm

Read through the sorted sam file line by line. If the line starts with an "@" sign, just toss it straight into our output file since none of that is changing.

Let's initialize a position variable, `left_pos`, and a read direction variable/if it is reverse complimented, `rev`, that we can use for our comparisons here.

Does our position (column 4) AND read direction `rev` (column 2 -- bitwise flag: extract bit 16) match the ones we've got saved from the last read?
    No?
        Okay, we're treading new territory here, time to set some stuff up for comparisons.
        Let's set `left_pos` to our current position (whenever we change directions on the same position this will rewrite to the same value as before but that's okay).
        Let's set our read direction to our current read direction.
        Let's set up/reset an empty set of UMIs, `current_umis`, that we have found at this position and direction.
        Now we're good to move on to checking our UMI info.
Now we can check against UMIs
    Cool. We need to see what the UMI is in our read. That's not too bad, but it has a few complications, so we'll need to make a function `extract_umi` for this (described in the functions section below).
    Is that umi in the set of known umis?
        Yes?
            Cool. But is it in the set of umis we've already had at this position and read direction?
                Yes?
                    Looks like a PCR duplicate -- anything in the set of UMIs also shares the same position and direction as this read. Don't write this to the output file.
                Nope.
                    Cool cool. Let's write it to the output file and add this to our `current_umis`.

## Functions

`extract_umi`

```python
def extract_umi(rev: bool, cigar: str, seq: str) -> str:
    '''
    Takes whether the sequence is reverse complemented,
    using the CIGAR string, finds the start of the 
    '''
    return umi
```
