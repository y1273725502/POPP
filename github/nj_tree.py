import os

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from matplotlib import pyplot as plt


def nj_Tree(time):
    # 读取 fasta 文件
    records = list(SeqIO.parse(f"{time}\\temp.fasta", "fasta"))

    # 构建多序列比对对象
    alignment = MultipleSeqAlignment(records)

    # 计算距离矩阵
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)

    # 使用 Neighbor-Joining 构建树
    nj_constructor = DistanceTreeConstructor()
    nj_tree = nj_constructor.nj(dm)


    # 保存 Neighbor-Joining 树为 Newick 文件
    nj_tree_file = f"{time}/temp.fasta.treefile"
    Phylo.write(nj_tree, nj_tree_file, "newick")
    draw_tree(nj_tree_file,time)

def draw_tree(tree_file, time):
    tree = Phylo.read(tree_file, "newick")
    Phylo.draw(tree,branch_labels=lambda c: None, do_show=False)
    # 删除原来的png文件
    if os.path.exists(f"{time}/tree.png"):
        os.remove(f"{time}/tree.png")
    plt.savefig(f"{time}/tree.png")