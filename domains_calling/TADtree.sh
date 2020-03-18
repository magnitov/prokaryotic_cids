python <path to TADtree.py> <TADtree_control_file.txt>
awk 1 N*.txt | sort -k1,1 -k2,2n | uniq > all_TADs.txt


