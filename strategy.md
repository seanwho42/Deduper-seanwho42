## Overview

We need to be able to read through *sorted* sam files and remove the PCR duplicates, without undue burden on the memory on the machines we are running this on. Example input and output files are included in `unit-test`. More details on those in the test files section.

In broad terms, we need to read through the alignments and compare with other reads at the same position and direction, and check the UMIs to see if they are PCR duplicates (or if the UMIs are erroneous).

Some context:
- Samtools sorts by reference, position in the reference, then by reverse flag per [documentation](https://www.htslib.org/doc/1.14/samtools-sort.html).
- Reads in sam files are tab separated, so we can split off of our tabs and then pull the relevant fields to get the information we need.
- This is banking on the (perhaps incorrect assumption due to sequencing error etc.) that PCR duplicates are aligning to the same position, so we check all reads at the same position which have the same soft clipping. Potential fixes down the line if this is a problem might look something like a sliding window situation instead of just looking at each position where the sliding window is as large as the largest length of softclipping?


## The algorithm
Read through `STL96.txt` to get the UMIs. Throw them into a set. (*Might make this into a function for tidyness reasons if it ends up making sense when the code is written*)

Read through the sorted sam file line by line. If the line starts with an "@" sign, just toss it straight into our output file since none of that is changing.

Let's initialize a chromosome variable, `chr`, position variable, `left_pos`, and a read direction variable/if it is reverse complimented, `rev`, that we can use for our comparisons here -- all as None so they won't match in the comparison deal

Does our chromosome (column 3), position (column 4) AND read direction `rev` (column 2 -- bitwise flag: extract bit 16) match the ones we've got saved from the last read?
- No?
    - Okay, we're treading new territory here, time to set some stuff up for comparisons.
    - Let's set the chromosome to our current chromosome
    - Let's set `left_pos` to our current position (whenever we change directions on the same position this will rewrite to the same value as before but that's okay).
    - Let's set our read direction to our current read direction.
    - Let's set up/reset an empty dictionary of UMIs, `current_umis`, that we have found at this position and direction. Keys will be the soft clipping number, values will be a set of UMIs

Let's get out soft clipping number. This is done with our function `get_soft_clipping` (described in functions section below)

We need to see what the UMI is in our read. To do this we can just take the last 8 characters off of the QNAME (column 1)

Now we can check against UMIs (*Might make this into a function for tidyness reasons if it ends up making sense when the code is written*)
- Is this umi in the set of known umis that should exist in this data?
    - Yes?
        - Is this soft clipping number a key in the dictionary of umis?
            - Yes?
                - Cool. But is the umi in the set of umis at this soft clipping number (**need to check our data to see if we have to search for barcodes as reverse compliments if direction is reverse -- probably do**)?
                    - Yes?
                        - Looks like a PCR duplicate -- anything in the set of UMIs also shares the same position and direction as this read. Don't write this to the output file.
                    - Nope.
                        - Cool cool. Let's write it to the output file and add this to our set of umis for this soft clipping number
            - No!
                - Cool cool. Write this to the output file, make a key value pair in the dictionary. Key is the soft clipping number, value is a set with the UMI in it



## Functions

`get_soft_clipping`

```python
def get_soft_clipping(cigar: str) -> int: # cigar string is column 6
    '''
    Takes CIGAR string, returns soft clipping number at the 5' end. CIGAR strings are 
    reversed for reverse complimented reads, so our 5' end is always the left end of
    the CIGAR string.

    This is done using regex. If there is not soft clipping, the number returned is 0.
    '''
    return soft_clip_number
```


## Test files
Unit files are taken from the provided test.sam file. Just took the first read and made modifications to it to support the different test cases. Lines are as such:
- header lines -- keep
- template line -- soft clipping is 0, not reverse complimented, cigar 0 -- keep
- same everything -- do not keep
- same everything except read direction (bitwise flag) -- keep
- same as 4 except including all things in cigar string which consume reference -- keep
- same as above except also adding soft clipping to the right (5') end -- keep
- same as 4 except including all things in cigar string which don't consume reference (except hard clipping) -- don't keep
- same everything except UMI (is in provided UMIs) -- keep
- same everything except UMI (is NOT in provided UMIs) -- do not keep
- soft clipping is 2 instead of 0 -- keep
- same as above -- do not keep
- same everything except chromosome -- keep
