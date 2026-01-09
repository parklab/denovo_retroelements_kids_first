#this is to do gene annotation of the simulated sites, and could by given criteria
from x_gene_annotation import *
import sys

####
def annotate_sites(sf_gene_annotation, sf_ori_input, sf_pLI_genes, sf_wfolder, f_min_pLI=0.9):
    if len(sf_wfolder)>0 and sf_wfolder[-1]!="/":
        sf_wfolder+="/"
    gff = GFF3(sf_gene_annotation)
    iextnd = global_values.UP_DOWN_GENE
    gff.load_gene_annotation_with_extnd(iextnd)
    gff.index_gene_annotation_interval_tree()

    m_high_pLI_genes={}
    with open(sf_pLI_genes) as fin_pLI:
        for line in fin_pLI:
            if "pLI" in line:
                continue
            fields=line.rstrip().split()
            s_gene=fields[1]
            f_pLI=float(fields[-4])
            if f_pLI>=f_min_pLI:
                m_high_pLI_genes[s_gene]=f_pLI

    l_gene_count=[]
    sf_out_cnt=sf_wfolder+"summarized_cnt_all_rounds.txt"
    with open(sf_ori_input) as fin_ori, open(sf_out_cnt,"w") as fout_all_rounds:
        i_round=1
        for line in fin_ori:
            fields=line.rstrip().split()
            sf_sites_round_i=sf_wfolder + "round_%d_sites.txt"%i_round
            with open(sf_sites_round_i,"w") as fout_round_i:
                for s_term in fields:
                    l_term=s_term.split(":")
                    fout_round_i.write(l_term[0]+" "+l_term[1]+"\n")
            ####
            sf_out_round_i=sf_wfolder + "round_%d_gene_annotation.txt"%i_round
            gff.annotate_results(sf_sites_round_i, sf_out_round_i)
            i_exon_cnt, i_high_pLI_cnt=cnt_annotated_sites(sf_out_round_i, m_high_pLI_genes)
            l_gene_count.append((i_exon_cnt, i_high_pLI_cnt))
            fout_all_rounds.write("%d %d\n" % (i_exon_cnt, i_high_pLI_cnt))
            i_round+=1
    return l_gene_count

####
def cnt_annotated_sites(sf_input, m_pLI_gene):
    i_exon_cnt=0
    i_high_pLI_cnt=0
    with open(sf_input) as fin_sites:
        for line in fin_sites:
            fields=line.rstrip().split()
            s_gene=fields[-1]#intron:ENSG00000198947.15:DMD
            l_gene_fields=s_gene.split(":")
            if len(l_gene_fields)<3:
                continue
            s_region=l_gene_fields[0]
            if (s_region=="exon") or ("UTR" in s_region):
                i_exon_cnt+=1
            s_gene=l_gene_fields[-1]
            if s_gene in m_pLI_gene:
                if (s_region=="exon") or (s_region=="intron") or ("UTR" in s_region):
                    i_high_pLI_cnt+=1
    return i_exon_cnt, i_high_pLI_cnt
####

####
if __name__ == '__main__':
    sf_gene_annotation=sys.argv[1]
    sf_ori_input=sys.argv[2]
    sf_pLI_genes=sys.argv[3]
    sf_wfolder = sys.argv[4]
    annotate_sites(sf_gene_annotation, sf_ori_input, sf_pLI_genes, sf_wfolder)