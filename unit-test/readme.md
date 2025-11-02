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
