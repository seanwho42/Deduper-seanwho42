#!/usr/bin/env python
import argparse
import re
# for error correction
from rapidfuzz.process import extract
from rapidfuzz.distance import Hamming


def get_args():
    parser = argparse.ArgumentParser(description="Generates a copy of a sorted SAM file after removing PCR duplications given provided UMIs.")
    # TODO: include the file formatting info for the index file in the help? or reference to the relevant function
    parser.add_argument("-f", "--file", help = "Sorted SAM file to have PCR duplicates removed from", required = True)
    parser.add_argument("-o", "--out", help = "Output path to sorted SAM file without PCR duplicates", required = True)
    parser.add_argument("-u", "--umi", help = "File path for file containing the list of valid UMIs", required=True)
    parser.add_argument("-c", "--correction", action='store_true', 
                        help = "Enable error correction to find the closeset fuzzy match (using Hamming distance <= 1)")
    # parser.add_argument("-w", "--window", help = "Window size (in base pairs) to compare 5' start position within", type = int)
    # TODO: add the help text here
    # parser.add_argument("-h", "--help", help="")
    return parser.parse_args()

args = get_args()


def main():
    '''
    Main function. Generates a copy of a sorted SAM file after removing PCR duplications given provided UMIs.
    '''
    # initialize summary statistics to be reported at the end
    header_lines = 0
    unique_records = 0
    wrong_umis = 0
    duplicates_removed = 0
    # keys chromosome, value is count
    chrom_counts = {}
    chrom_read_count = 0

    errors_ok = args.correction
    print(f"errors_ok = {errors_ok}")

    # debugging
    unique_umis = set()

    # keys will be 5' position, values will be a set of tuples for each non-pcr duplicate in that position
    # positions will be cleared as new 5' positions get further away from the 
    # read in the file with the umis
    with open(args.umi, "r") as umi_file:
        umis = set()
        umi_length = None
        # support different UMI lengths as long as they are consistent
        for line in umi_file:
            umi = line.strip()
            if umi_length is None:
                umi_length = len(umi)
            elif len(umi) != umi_length:
                raise Exception("Provided UMI files have inconsistent lengths. Deduper does not support this.")
            umis.add(umi)

    umis_dict = {}
    if errors_ok:
        # need to be a dictionary or a list for rapidfuzz.. not sure which is faster
        # TODO: check if dictionary is actually faster here or if it doesn't matter
        # for the library since it's converting it all to its own stuff anyways
        for umi in list(umis):
            # key value is what needs to be returned, string is what we are searching
            umis_dict[umi] = umi
    
    # initializing last known chromosome variable to compare against
    last_known_chrom = None
    with open(args.file, "r") as f, open(args.out, "w") as o:
        # setting this up as a dictionary with 5' pos as key for future proofing to make implementing a window later easier
        reads_dict = {}
        for line in f:
            # print(unique_records)
            # print(reads_dict)
            if line[0] == "@":
                o.write(line)
                header_lines += 1
            else:
                umi, read_rev, chrom, left_pos, cigar_str = parse_read(line, umi_length)
                unique_umis.add(umi)
                # wipe out the reads dictionary if we're looking at a new chromosome
                # this is necessary even with the sliding window approach, since the
                # position restarts for every chromosome
                if chrom != last_known_chrom:
                    reads_dict = {}
                    if last_known_chrom is not None:
                        chrom_counts[last_known_chrom] = chrom_read_count
                        chrom_read_count = 0

                # only want to continue our analysis if the UMI is valid

                umi, umi_valid = check_umi(umi, errors_ok, umis, umis_dict)
                
                if umi_valid:
                    # TODO: condense into helper function
                    # now we can get going on our process to check for pcr duplicates
                    five_p_pos = get_five_p_pos(left_pos, cigar_str, read_rev)
                    # metadata of relevance to identifying PCR duplicates that isn't the 5' start position to be saved in the dict
                    read_meta = (umi, read_rev, chrom)
                    if five_p_pos in reads_dict:
                        if read_meta not in reads_dict[five_p_pos]:
                            reads_dict[five_p_pos].add(read_meta)
                            o.write(line)
                            chrom_read_count += 1
                            unique_records += 1
                        else:
                            duplicates_removed += 1
                    else:
                        chrom_read_count += 1
                        unique_records += 1
                        reads_dict[five_p_pos] = set()
                        reads_dict[five_p_pos].add(read_meta)
                        o.write(line)
                else:
                    print(f"umi not found: {umi}")
                    wrong_umis += 1
                last_known_chrom = chrom
    print(f'Header lines: {header_lines}\n' +
          f'Unique records: {unique_records}\n' +
          f'Wrong UMIs: {wrong_umis}\n' +
          f'Duplicates removed: {duplicates_removed}')

    # add info for the last chromosome
    chrom_counts[last_known_chrom] = chrom_read_count

    # for chrom, count in chrom_counts.items().sorted():
    #     print(f'{chrom}\t{count}')

    for key in sorted(chrom_counts.keys()):
        print(f"{key}\t{chrom_counts[key]}")

    unique_umis = list(unique_umis)
    unique_umis.sort()

    with open("unique-umis.txt", "w") as umi_file:
        for umi in unique_umis:
            umi_file.write(f"{umi}\n")

# and parsing some basic info from our read
def parse_read(line: str, umi_length):
    # TODO: write the docstring
    '''
    Parses basic information from a read given that read as a string.
    Returns a the UMI, read direction, chromosome, left position, and cigar string.
    '''
    split_line = line.split("\t")

    qname = split_line[0]
    # now parse the umi from the qname
    umi = qname[-umi_length:]

    bitwise_flag = int(split_line[1])
    # and now to parse bit 16 from this
    read_rev = False
    if (bitwise_flag & 16) == 16:
        read_rev = True

    chrom = split_line[2]
    left_pos = int(split_line[3])
    cigar_str = split_line[5]
    # print((umi, read_rev, chrom, left_pos, cigar_str))

    return umi, read_rev, chrom, left_pos, cigar_str

def get_five_p_pos(left_pos: int, cigar_str: str, read_rev):
    # TODO: write the docstring
    '''
    Given the left position aligned to the reference, the cigar string, and the read direction,
    of a given read determines and returns its 5' start position.
    '''
    # check for read direction in here

    # finds all instances of a pattern and returns the groups -- first being the digit and second being the alignment "type"
    cigar_list = re.findall(r'(\d+)(\D)', cigar_str)
    five_p_pos = left_pos

    if read_rev:
        # reverse read -- cigar string runs left to right, which is 3' to 5'
        # do a loop through each item in the list, add it to the count if it consumes reference per documentation
        # also add to count for soft clipping on the end
        consume_ref_set = set(['M', 'D', 'N', '=', 'X'])
        for i, (cigar_val, cigar_char) in enumerate(cigar_list):
            cigar_val = int(cigar_val)
            if cigar_char == 'S':
                if i != 0:
                    # it is at the end -- so the 5' side, so we add it
                    five_p_pos += cigar_val
            else:
                if cigar_char in consume_ref_set:
                    # consumes reference, so we add the value associated with it
                    five_p_pos += cigar_val
    else:
        # forward read -- cigar string runs left to right, which is 5' to 3'
        five_p_int, five_p_char = cigar_list[0]
        five_p_int = int(five_p_int)
        if five_p_char == 'S':
            five_p_pos -= five_p_int
    return five_p_pos

def check_umi(umi, errors_ok, umis, umis_dict):
    # TODO: add support for custom Hamming distances (with warnings/suggested ranges in the help file?)
    '''
    With error correcting enabled:
    Checks and rewrites/finds UMIs within Hamming distance of one. Returns if UMI is valid/found in the provided
    UMI list with a hamming distance cutoff of one. If multiple UMIs share a hamming distance of one, the UMI
    is found to be invalid.

    Without error correcting enabled:
    Returns tuple of UMI and validity
    '''
    umi_valid = False

    if errors_ok:
        # do the fuzzy match off of the given UMIs
        matched_umis = extract(umi, umis_dict, scorer=Hamming.distance, limit=2, score_cutoff=1)
        # print(matched_umis)
        if len(matched_umis) == 2:
            # check to see if the score is the same for these two...
            if matched_umis[0][1] == matched_umis[1][1]:
                # TODO: log file instead of print/only print with verbose flag
                print(f'UMI "{umi}" shares Hamming distance with multiple known UMIs. Omitting.')
            else:
                umi_valid = True
                umi = matched_umis[0][2]
        elif len(matched_umis) == 1:
            umi_valid = True
            umi = matched_umis[0][2]
    else:
        if umi in umis:
            umi_valid = True
    return umi, umi_valid
    

main()