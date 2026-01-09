#
#this is to simulation insertion sites based on cleavage seq frequency
##
##
import pysam
from global_values import *
from x_annotation import *
from x_black_list import *
from intervaltree import *
import random

####
MAX_DIV_RATE=5
MIN_COPY_LEN=250

####
def simulate_sites(sf_ref, sf_motif, sf_rmsk, sf_black_list, i_rounds, i_sites, sf_out):
    m_motif, m_motif_rc = _load_in_cleavage_motifs(sf_motif)
    l_seg1, l_seg2, l_seg3, i_total = gnrt_motif_segmts_list(m_motif)

    m_gnrted_motif={}
    m_gnrted_motif_rc={}
    l_gnrted_motif_all_rounds=[]
    for i in range(i_rounds):#for example, run 100000 rounds
        l_gnrted_motif_round_i=[]
        for j in range(i_sites):#each round, generate the specific number of insertions
            #first generate the motif
            s_motif=_simulate_motif_of_one_site(l_seg1, l_seg2, l_seg3)
            l_gnrted_motif_round_i.append(s_motif)
            #second get the positions of this motif
            #randomly select one position and save it
            s_motif_rc=_gnrt_reverse_complementary(s_motif)
            if s_motif not in m_gnrted_motif:
                m_gnrted_motif[s_motif]=1
            if s_motif_rc not in m_gnrted_motif_rc:
                m_gnrted_motif_rc[s_motif_rc]=s_motif

        l_gnrted_motif_all_rounds.append(l_gnrted_motif_round_i)

    sf_gnrted_motifs = sf_out + ".gnrted_motifs_all_rounds.txt"
    with open(sf_gnrted_motifs,"w") as fout_motifs:
        i_round=1
        for l_round_i in l_gnrted_motif_all_rounds:
            for s_motif in l_round_i:
                fout_motifs.write(s_motif+" "+str(i_round)+"\n")
            i_round+=1

    ####load in blacklist regions
    xbl=XBlackList()
    xbl.load_index_regions(sf_black_list)
    ####
    ####prepare the repeats annotation
    xannotation = _prepare_repeats_annotation(sf_rmsk)
    ####
    sf_motif_sites_table=sf_out+".gnrted_motif_ref_position.txt"
    m_motif_pos = _collect_motif_position_from_ref(sf_ref, m_gnrted_motif, m_gnrted_motif_rc, xbl, xannotation,
                                                   sf_motif_sites_table)

    sf_out_sites=sf_out+".gnrted_sites.txt"
    with open(sf_out_sites,"w") as fout_sites:
        for l_motifs in l_gnrted_motif_all_rounds:
            for s_motif in l_motifs:
                i_slct_chrm, i_slct_pos=random_slct_one_position(m_motif_pos, s_motif)
                s_site="%s:%d"%(i_slct_chrm, i_slct_pos)
                fout_sites.write(s_site+" ")
            fout_sites.write("\n")
    ####

def _prepare_repeats_annotation(sf_rmsk):
    boundary_extnd=0
    i_min_copy_len=MIN_COPY_LEN
    xannotation = XAnnotation(sf_rmsk)
    b_with_chr = True
    xannotation.set_with_chr(b_with_chr)
    xannotation.load_rmsk_annotation_with_extnd_div_with_lenth_cutoff(boundary_extnd, i_min_copy_len)
    xannotation.index_rmsk_annotation_interval_tree()
    return xannotation

####
def check_whether_in_repetitive_regions(xannotation, ins_chrm, ins_pos):
    b_in_rep, i_pos = xannotation.is_within_repeat_region_interval_tree(ins_chrm, int(ins_pos))
    if b_in_rep==False:
        return False
    div_rate, sub_family, family, pos_start, pos_end = xannotation.get_div_subfamily(ins_chrm, i_pos)
    if div_rate<MAX_DIV_RATE:
        return True
    return False
#
####
def random_slct_one_position(m_motif_pos, s_motif):
    m_chrm_pos=m_motif_pos[s_motif]
    l_fixed_chrm=[]
    i_total_pos=0
    for s_tmp_chrm in m_chrm_pos:
        i_total_pos+=len(m_chrm_pos[s_tmp_chrm])
        l_fixed_chrm.append(s_tmp_chrm)
    i_slct_idx=random.randrange(i_total_pos)

    i_cur_pos=0
    i_slct_chrm=""
    i_slct_pos=0
    for s_tmp_chrm in l_fixed_chrm:
        i_len=len(m_chrm_pos[s_tmp_chrm])
        if i_slct_idx>=i_cur_pos and i_slct_idx<(i_cur_pos+i_len):
            i_slct_chrm=s_tmp_chrm
            i_diff=i_slct_idx-i_cur_pos
            i_slct_pos=m_chrm_pos[s_tmp_chrm][i_diff]
            break
        else:
            i_cur_pos+=i_len
    return i_slct_chrm, i_slct_pos #
####

####
def _simulate_motif_of_one_site(l_seg1, l_seg2, l_seg3):
    #first simulate the motif, and get the frequency
    s_seg1=random.choice(l_seg1)
    s_seg2=random.choice(l_seg2)
    s_seg3=random.choice(l_seg3)
    s_motif=s_seg1+s_seg2+s_seg3
    return s_motif


#three groups: first char, 2-5, and 6-7
#generate the lists for random select
def gnrt_motif_segmts_list(m_motif):
    m_seg1={}
    m_seg2={}
    m_seg3={}
    i_total=0
    for s_motif in m_motif:
        i_freq=m_motif[s_motif]
        i_total += i_freq
        s_seg1=s_motif[0]
        s_seg2 = s_motif[1:5]
        s_seg3=s_motif[5:]

        ####
        if s_seg1 not in m_seg1:
            m_seg1[s_seg1]=i_freq
        else:
            m_seg1[s_seg1] += i_freq

        ####
        if s_seg2 not in m_seg2:
            m_seg2[s_seg2] = i_freq
        else:
            m_seg2[s_seg2] += i_freq
        #####
        if s_seg3 not in m_seg3:
            m_seg3[s_seg3] = i_freq
        else:
            m_seg3[s_seg3] += i_freq

    ####
    l_seg1=[]
    l_seg2=[]
    l_seg3=[]
    for s_seg1 in m_seg1:
        i_freq=m_seg1[s_seg1]
        for i in range(i_freq):
            l_seg1.append(s_seg1)
    for s_seg2 in m_seg2:
        i_freq=m_seg2[s_seg2]
        for i in range(i_freq):
            l_seg2.append(s_seg2)
    ####
    for s_seg3 in m_seg3:
        i_freq=m_seg3[s_seg3]
        for i in range(i_freq):
            l_seg3.append(s_seg3)
    ####
    return l_seg1, l_seg2, l_seg3, i_total

####
def _collect_motif_position_from_ref(sf_ref, m_motif, m_motif_rc, xblacklist, xannotation, sf_out):
    m_motif_pos={}
    f_fa = pysam.FastaFile(sf_ref)##
    n_total_sites=0
    for tmp_chrm in f_fa.references:
        l_tmp=tmp_chrm.split("_")
        if len(l_tmp)>1:
            continue
        s_chrm_seq = f_fa.fetch(tmp_chrm)
        i_chrm_len=len(s_chrm_seq)-100
        for i in range(i_chrm_len):
            s_tmp_kmer=s_chrm_seq[i:i+7]

            i_cleavage_pos = -1
            s_motif=s_tmp_kmer
            if s_tmp_kmer in m_motif:
                i_cleavage_pos=i+5
            if s_tmp_kmer in m_motif_rc:
                i_cleavage_pos=i+2
                s_motif=m_motif_rc[s_tmp_kmer]
            if i_cleavage_pos>0:
                #if fall in the blacklist regions, then skip
                b_hit, i_hit_pos = xblacklist.fall_in_region(tmp_chrm, i_cleavage_pos)
                if b_hit==True:
                    continue
                #if fall in low diverged the repetitive regions of the same repeat type, then skip
                b_hit_low_rep=check_whether_in_repetitive_regions(xannotation, tmp_chrm, i_cleavage_pos)
                if b_hit_low_rep==True:
                    continue

                if s_motif not in m_motif_pos:
                    m_motif_pos[s_motif]={}
                if tmp_chrm not in m_motif_pos[s_motif]:
                    m_motif_pos[s_motif][tmp_chrm]=[]
                m_motif_pos[s_motif][tmp_chrm].append(i_cleavage_pos)
                n_total_sites+=1
    f_fa.close()

    print("Total number of potential cleavage sites is %d\n" % n_total_sites)

    with open(sf_out, "w") as fout:
        for s_motif in m_motif_pos:
            s_info=s_motif
            for s_chrm in m_motif_pos[s_motif]:
                s_info+=(" "+s_chrm)
                for i_pos in m_motif_pos[s_motif][s_chrm]:
                    s_info+=(" "+str(i_pos))
            fout.write(s_info+"\n")
    return m_motif_pos
#
####
def _load_in_cleavage_motifs(sf_motif):
    m_motif={}
    m_motif_rc={}
    with open(sf_motif) as fin_motif:
        for line in fin_motif:
            fields=line.rstrip().split()
            s_motif=fields[0]
            i_freq=int(fields[1])
            m_motif[s_motif]=i_freq
            s_motif_rc=_gnrt_reverse_complementary(s_motif)
            m_motif_rc[s_motif_rc]=s_motif
    return m_motif, m_motif_rc

####
def _set_chr_map(chr_map):
    chr_map['A'] = 'T'
    chr_map['a'] = 'T'
    chr_map['T'] = 'A'
    chr_map['t'] = 'A'
    chr_map['C'] = 'G'
    chr_map['c'] = 'G'
    chr_map['G'] = 'C'
    chr_map['g'] = 'C'
    chr_map['N'] = 'N'
    chr_map['n'] = 'N'

def _gnrt_reverse_complementary(s_seq):
    chr_map = {}
    _set_chr_map(chr_map)
    s_rc = ""
    for s in s_seq[::-1]:
        if s not in chr_map:
            s_rc += "N"
        else:
            s_rc += chr_map[s]
    return s_rc


if __name__ == '__main__':
    sf_ref=sys.argv[1]
    sf_motif=sys.argv[2]
    sf_rmsk=sys.argv[3]
    sf_black_list=sys.argv[4]
    i_sites=int(sys.argv[5]) #i_sites_L1=30, #i_sites_SVA=35, #i_sites_Alu=95
    sf_out=sys.argv[6]
    i_rounds=10000

    simulate_sites(sf_ref, sf_motif, sf_rmsk, sf_black_list, i_rounds, i_sites, sf_out)
####