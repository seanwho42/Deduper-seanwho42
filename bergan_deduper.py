#!/usr/bin/env python
import argparse
import re


def get_args():
    parser = argparse.ArgumentParser(description="Generates a copy of a sorted SAM file after removing PCR duplications given provided UMIs.")
    # TODO: include the file formatting info for the index file in the help? or reference to the relevant function
    parser.add_argument("-f", "--file", help = "Sorted SAM file to have PCR duplicates removed from", required = True)
    parser.add_argument("-o", "--out", help = "Output path to sorted sam file without PCR duplicates", required = True)
    parser.add_argument("-u", "--umi", help = "File path for file containing the list of valid UMIs", required=True)
    # parser.add_argument("-w", "--window", help = "Window size (in base pairs) to compare 5' start position within", type = int)
    # TODO: add the help text here
    # parser.add_argument("-h", "--help", help="")
    return parser.parse_args()

args = get_args()


def main():
    '''
    Main function. Generates a copy of a sorted SAM file after removing PCR duplications given provided UMIs.
    '''
    # keys will be 5' position, values will be a set of tuples for each non-pcr duplicate in that position
    # positions will be cleared as new 5' positions get further away from the 
    # read in the file with the umis
    with open(args.umi, "r") as umi_file:
        umis = set()
        for line in umi_file:
            umis.add(line.strip())
    
    # initializing last known chromosome variable to compare against
    last_known_chrom = None
    with open(args.file, "r") as f, open(args.out, "w") as o:
        for line in f:
            # write the line if it is a header line
            if line[0] == "@":
                o.write(line)
            else:
                umi, read_rev, chrom, left_pos, cigar_str = parse_read(line)

                # wipe out the reads dictionary if we're looking at a new chromosome
                # this is necessary even with the sliding window approach, since the
                # position restarts for every chromosome
                if chrom != last_known_chrom:
                    reads_dict = {}

                # only want to continue our analysis if the UMI is valid
                # TODO: implement error correction potentially?
                if umi in umis:
                    # now we can get going on our process to check for pcr duplicates
                    five_p_pos = get_five_p_pos(left_pos, cigar_str, read_rev)
                    # metadata of relevance to identifying PCR duplicates that isn't the 5' start position to be saved in the dict
                    read_meta = (umi, read_rev, chrom)
                    if five_p_pos in reads_dict:
                        if read_meta not in reads_dict[five_p_pos]:
                            reads_dict[five_p_pos].add(read_meta)
                    else:
                        reads_dict[five_p_pos] = set(read_meta)
                        o.write(line)
                last_known_chrom = chrom

# and parsing some basic info from our read
def parse_read(line: str):
    '''
    Parses basic information from a read given that read as a string.
    Returns a the UMI, read direction, chromosome, left position, and cigar string.
    '''
    split_line = line.split("\t")

    qname = split_line[0]
    # now parse the umi from the qname
    umi = qname[-8:]

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
    '''
    Given the left position aligned to the reference, the cigar string, and the read direction,
    of a given read determines and returns its 5' start position.
    '''
    # check for read direction in here

    # finds all instances of a pattern and returns the groups -- first being the digit and second being the alignment "type"
    cigar_list = re.findall(r'(\d+)(\w)', cigar_str)
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



main()