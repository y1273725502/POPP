import os

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from matplotlib import pyplot as plt


def nj_Tree(time):
    # ��ȡ fasta �ļ�
    records = list(SeqIO.parse(f"{time}\\temp.fasta", "fasta"))

    # ���������бȶԶ���
    alignment = MultipleSeqAlignment(records)

    # ����������
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)

    # ʹ�� Neighbor-Joining ������
    nj_constructor = DistanceTreeConstructor()
    nj_tree = nj_constructor.nj(dm)


    # ���� Neighbor-Joining ��Ϊ Newick �ļ�
    nj_tree_file = f"{time}/temp.fasta.treefile"
    Phylo.write(nj_tree, nj_tree_file, "newick")
    draw_tree(nj_tree_file,time)

def draw_tree(tree_file, time):
    tree = Phylo.read(tree_file, "newick")
    Phylo.draw(tree,branch_labels=lambda c: None, do_show=False)
    # ɾ��ԭ����png�ļ�
    if os.path.exists(f"{time}/tree.png"):
        os.remove(f"{time}/tree.png")
    plt.savefig(f"{time}/tree.png")