## Test files

Here are the unit test files for Deduper.

If there is a missing case which you have encountered which needs to be tested for, please contact Sean Bergan (sbergan@uoregon.edu).

### test.sam

Unit files are taken from the provided test.sam file. Just took the first read and made modifications to it to support the different test cases. Lines are as such:
- header lines -- keep
- template line -- soft clipping is 0, not reverse complimented, cigar 0 -- keep
- same everything -- do not keep
- same everything except read direction (bitwise flag) -- keep
- same as 4 except including all things in cigar string which consume reference -- do not keep
- same as above except also adding soft clipping to the right (5') end -- keep
- same as 4 except including all things in cigar string which don't consume reference (except hard clipping) -- do not keep
- same everything except UMI (is in provided UMIs) -- keep
- same everything except UMI (is NOT in provided UMIs + too big Hamming distance for EC) -- do not keep
- soft clipping is 2 instead of 0 -- keep
- same as above -- do not keep
- same everything except chromosome -- keep
- same everything except UMI (is error by two characters from provided UMIs) -- do not keep
- same everything except UMI (is error by one character from provided UMIs) -- keep if error correction flag present
- same as above -- do not keep (checking if we still recognize which UMI is being used)
- same as above but with no errors in UMI -- (do not keep)

### provided_umis.txt

Properly formatted UMIs file.

### inconsistent_length_umis.txt

Improperly formatted UMIs file (inconsistent UMI length).

### expected-output-ec.sam

What we expect to see output from the provided sam and UMI file (with error correction).

### expected-output-no-ec.sam

What we expect to see output from the provided sam and UMI file (without error correction).
