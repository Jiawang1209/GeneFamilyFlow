# get_motif_info.py
import argparse
import re

def get_motif_info(arg1):
    # output file
    fw1 = open("meme_info.txt", "w")
    fw2 = open("meme_location.txt", "w")

    # read file
    with open(arg1, "rt") as fr:
        for line in fr:
            line_tmp = line.strip()
            # motif information
            p = re.compile('Motif.*Description')
            p2 = re.compile(r'Sobic.\d+\w+\d+\.\d\.p\s+\d+\s+\d+.*')
            if p.match(line_tmp):
                motif_info = line_tmp.split()
                fw1.write(motif_info[2] + "\t" + motif_info[1] + "\t" + str(len(motif_info[1])) + "\n")
                motif_name = motif_info[2]

            # gene motif location
            if p2.match(line_tmp):
                motif_location = line_tmp.split()
                fw2.write(motif_location[0] + "\t" + motif_location[1] + "\t" + str(int(motif_location[1]) + len(motif_location[4]) -1) + "\t" + motif_name + "\n")

    fw1.close()
    fw2.close()




if __name__ == "__main__":
    # create ArgumentParser
    parser = argparse.ArgumentParser(description="This Program will get motif information from meme output")

    # add argument
    parser.add_argument("--arg1", help="Input meme output file")

    # parser argument
    args = parser.parse_args()

    get_motif_info(args.arg1)
