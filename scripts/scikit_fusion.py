import numpy as np
import os
import sys

from skfusion import fusion


if len(sys.argv) < 5:
    print("Arguments: <input folder> <output folder> <nb factors> <separator> <datasets...>")
    sys.exit()

source_folder = sys.argv[1]
output_folder = sys.argv[2]
nfactors = sys.argv[3]
sep = sys.argv[4]
datasets = sys.argv[5:]


def strip_first_col(fname, delimiter=sep):
    with open(fname, 'r') as fin:
        for line in fin:
            try:
               yield line.split(delimiter, 1)[1]
            except IndexError:
               continue



t1 = fusion.ObjectType('Type 1', nfactors)
tdata = [ fusion.ObjectType(dataset, nfactors) for dataset in datasets ]
relations = []
for i in range(len(datasets)):
    relations.append(
        fusion.Relation(
            np.transpose(
                np.loadtxt( strip_first_col(os.path.join(source_folder, datasets[i])), delimiter=sep, skiprows=1)
            ), 
            t1, tdata[i]
        )
    )

fusion_graph = fusion.FusionGraph()
fusion_graph.add_relations_from(relations)
print(fusion_graph)

fuser = fusion.Dfmf()
fuser.fuse(fusion_graph)


np.savetxt(os.path.join(output_folder,"signals.txt"), fuser.factor(t1), delimiter='\t')
for i in range(len(datasets)):
    np.savetxt(os.path.join(output_folder, "proj%s" %datasets[i]), fuser.factor(tdata[i]), delimiter='\t')

