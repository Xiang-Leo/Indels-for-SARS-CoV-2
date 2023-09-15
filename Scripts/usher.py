# Copyright Contributors to the Pyro-Cov project.
# SPDX-License-Identifier: Apache-2.0
import sys
sys.path.append('../')
from params import *
import gzip
import heapq
import logging
import shutil
import warnings
from collections import defaultdict, namedtuple
from typing import Dict, FrozenSet, Optional, Set, Tuple
import os,pickle
import tqdm,copy
from Bio.Phylo.NewickIO import Parser, Writer

from . import pangolin
from .external.usher import parsimony_pb2


logger = logging.getLogger(__name__)

Mutation = namedtuple("Mutation", ["position", "ref", "mut"])

NUCLEOTIDE = "ACGTN"

MUTATION_KEEP = True


def load_proto(filename):
    open_ = gzip.open if filename.endswith(".gz") else open
    with open_(filename, "rb") as f:
        proto = parsimony_pb2.data.FromString(f.read())  # type: ignore
    newick = proto.newick.replace(";", "")  # work around unescaped node names
    tree = next(Parser.from_string(newick).parse())
    return proto, tree


def load_usher_clades(filename: str) -> Dict[str, Tuple[str, str]]:
    """
    Loads usher's output clades.txt and extracts the best lineage and a list of
    possible lineages, for each sequence.
    """
    clades: Dict[str, Tuple[str, str]] = {}
    with open(filename) as f:
        for line in f:
            name, lineages = line.strip().split("\t")
            # Split histograms like B.1.1.161*|B.1.1(2/3),B.1.1.161(1/3) into points
            # like B.1.1.161 and lists like B.1.1,B.1.1.161.
            if "*|" in lineages:
                lineage, lineages = lineages.split("*|")
                lineages = ",".join(part.split("(")[0] for part in lineages.split(","))
            else:
                assert "*" not in lineages
                assert "|" not in lineages
                lineage = lineages
            clades[name] = lineage, lineages
    return clades


def load_mutation_tree(
    filename: str,
) -> Tuple[Dict[str, FrozenSet[Mutation]], object, object]:
    """
    Loads an usher lineageTree.pb or lineageTree.pb.gz annotated with mutations
    and pango lineages, and creates a mapping from lineages to their set of
    mutations.
    """
    logger.info(f"Loading tree from {filename}")
    proto, tree = load_proto(filename)

    # Extract phylogenetic tree.
    clades = list(tree.find_clades())
    assert len(proto.metadata) == len(clades)
    assert len(proto.node_mutations) == len(clades)

    # Map lineages to clades.
    lineage_to_clade = {
        str(meta.clade): clade
        for clade, meta in zip(clades, proto.metadata)
        if meta and meta.clade
    }
    # proto_clade_list = [meta.clade for meta in proto.metadata if meta and meta.clade]
    # if DEBUG:
    #     logger.info(f"proto clade list:{len(proto_clade_list)}, set: {len(set(proto_clade_list))}")
        # logger.info(f"{list(lineage_to_clade.items())[:10]}")
    # Accumulate mutations in each clade, which are overwritten at each position.
    logger.info(f"Accumulating mutations on {len(clades)} nodes")
    clade_to_muts: Dict[object, Dict[int, Mutation]] = defaultdict(dict)
    # if DEBUG:
    #     # logger.info(f"type mut {type(proto.node_mutations)}")
    #     clade_count = 0
    #     for clade in clades:
    #         clade_count += len(list(clade.clades))
    #     # logger.info(f"len {clade_count} {len(list(clades[0].clades))} clade {list(clades[0].clades)}")
    #     logger.info(f"clade count{clade_count} ")
    for i_clade,(clade, muts) in enumerate(zip(tqdm.tqdm(clades), proto.node_mutations)):
        for mut in muts.mutation:
            # if DEBUG:
                # print(f"pos {mut.position} ref {NUCLEOTIDE[mut.ref_nuc]} mut {''.join(NUCLEOTIDE[n] for n in mut.mut_nuc)}")
            if 5 in mut.mut_nuc: #for del_kind==1
                # clade_to_muts[clade][mut.position] = Mutation(mut.position,NUCLEOTIDE[mut.ref_nuc],"".join(["N",str(mut.mut_nuc)]))
                clade_to_muts[clade][mut.position] = Mutation(mut.position,NUCLEOTIDE[mut.ref_nuc],"N"*mut.mut_nuc[-1]) #mut.mut_nuc = [4,5,len]
            else:
                clade_to_muts[clade][mut.position] = Mutation(
                    mut.position,
                    NUCLEOTIDE[mut.ref_nuc],
                    "".join(NUCLEOTIDE[n] for n in mut.mut_nuc),
                )

        for c in clade.clades:
            clade_to_muts[c].update(clade_to_muts[clade])
        # if DEBUG:
        #     logger.info(f"{i_clade}th clade len clade {len(list(clade.clades))} mutation {len(list(muts.mutation))}")

    mutations_by_lineage = {
        k: frozenset(clade_to_muts[v].values()) for k, v in lineage_to_clade.items()
    }
    if DEBUG:
        logger.info(f"mutations_by_lineage: {len(mutations_by_lineage)} clade: {len(clades)} metadata: {len(proto.metadata)}")
        logger.info(f"clade set: {len(set(clades))} lineage2clade: {len(lineage_to_clade)}")
    return mutations_by_lineage, proto, tree

def load_reference_sequence():
    NEXTCLADE_DATA = os.path.expanduser("~/github/nextstrain/nextclade/data/sars-cov-2")
    print("reference path:",os.path.join(NEXTCLADE_DATA, "reference.fasta"))
    with open(os.path.join(NEXTCLADE_DATA, "reference.fasta")) as f:
        ref = "".join(line.strip() for line in f if not line.startswith(">"))
    assert len(ref) == 29903, len(ref)
    return ref

ref = load_reference_sequence()

def parse_indel(indel: str, dkind=0):
    """

    :param indel:
    :param dkind: 0 for single pos, 1 for multi pos
    :return:
    """
    assert indel, "no indel"
    indel = indel[:-1]
    gene, start, type, seq = indel.split(':')
    global ref
    indel_list = []
    if type == "ins":
        pos = int(start) # 1-based index
        ref_nuc = NUCLEOTIDE.index(ref[pos-1])
        par_nuc = NUCLEOTIDE.index(ref[pos-1])
        mut_nuc = [NUCLEOTIDE.index(c) for c in seq]
        indel_list.append([pos,ref_nuc,par_nuc,mut_nuc])
    elif type == "dels":
        if dkind == 0:
            pos = int(start)-1 #1-based index
            for del_ in seq:
                pos += 1
                ref_nuc = NUCLEOTIDE.index(ref[pos-1])
                assert NUCLEOTIDE.index(del_) == ref_nuc, f"del nuclead not equal ref \tpos: {pos}, del: {del_} ref: {ref_nuc}"
                par_nuc = NUCLEOTIDE.index(ref[pos-1])
                mut_nuc = [4]
                indel_list.append([pos,ref_nuc,par_nuc,mut_nuc])
        elif dkind == 1:
            pos = int(start) # 1-based index
            ref_nuc = NUCLEOTIDE.index(ref[pos-1])
            par_nuc = NUCLEOTIDE.index(ref[pos-1])
            mut_nuc = [4,5,len(seq)]
            indel_list.append([pos, ref_nuc, par_nuc, mut_nuc])
    else:
        print("indel type not found")
        exit()
    return indel_list

def mut_merge(clade,mutations):
    """merge mutations to parent node"""
    # print("merge muts to parent")
    if len(clade.clades) == 0:
        # pass
        return mutations[clade].mutation
    # mut_list = merge_muts(clade.clades[0], mutations).append(mutations[clade.clades[0]])
    clade_mut = mutations[clade].mutation
    child_muts = {}
    for clade_ in clade.clades:
        muts_ = mut_merge(clade_,mutations)
        child_muts.update({clade_:muts_})

    mut_list_comp = child_muts[clade_]

    # mut_list_kept = copy.deepcopy(mut_list_comp)
    # for child, mut_list_ in child_muts.items():
    #     for mut_ in mut_list_comp:
    #         if (not (mut_ in mut_list_)) and (mut_ in mut_list_kept):
    #             mut_list_kept.remove(mut_)
    #         if len(mut_list_kept) == 0:
    #             break
    #     if len(mut_list_kept) == 0:
    #         break
    #
    # if len(mut_list_kept) > 0:
    #     for mut_kept in mut_list_kept:
    #         clade_mut.append(mut_kept) # add to parent
    #         for clade_ in clade.clades:
    #             mutations[clade_].mutation.remove(mut_kept)

    child_num = len(child_muts)
    reserve_count = 0
    mut_list_kept = []
    for mut in mut_list_comp:
        for child, mut_list_ in child_muts.items():
            if mut in mut_list_:
                reserve_count += 1
        if reserve_count > 0.1*res_thresh*child_num:
            mut_list_kept.append(mut)
        reserve_count = 0
    # for child, mut_list_ in child_muts.items():
    #     for mut_ in mut_list_comp:
    #         if (not (mut_ in mut_list_)) and (mut_ in mut_list_kept):
    #             mut_list_kept.remove(mut_)
    #         if len(mut_list_kept) == 0:
    #             break
    #     if len(mut_list_kept) == 0:
    #         break

    if len(mut_list_kept) > 0:
        for mut_kept in mut_list_kept:
            if (len(mut_kept.mut_nuc) == 1) and MUTATION_KEEP: #snp not move in the tree
                continue
            clade_mut.append(mut_kept) # add to parent
            for clade_ in clade.clades:
                if mut_kept in mutations[clade_].mutation:
                    mutations[clade_].mutation.remove(mut_kept) # remove from child
    return mutations[clade].mutation


def refine_mutation_tree(filename_in: str, filename_out: str, indel_file: str) -> Dict[str, str]:
    """
    Refines a mutation tree clade metadata from pango lineages like B.1.1 to
    full node addresses like fine.0.12.4.1. Among clones, only the basal clade
    with have a .clade attribute, all descendents will have metadata.clade ==
    "". The tree structure remains unchanged.
    """
    proto, tree = load_proto(filename_in)
    logger.info(f"len proto meta {len(proto.metadata)}")
    # logger.info(f"len node meta {len(proto.node_metadata)}")
    # exit()

    # Extract phylogenetic tree.
    clades = list(tree.find_clades())
    logger.info(f"Refining a tree with {len(clades)} nodes")
    # if DEBUG:
    #     logger.info(f"metadata {len(proto.metadata)} clades {len(clades)}")
    #     logger.info(f"clades {prote.metadata}")
    assert len(proto.metadata) == len(clades)
    # logger.info(f"proto metadata {list(proto.metadata.keys())}")
    assert len(proto.node_mutations) == len(clades)
    metadata = dict(zip(clades, proto.metadata))
    mutations = dict(zip(clades, proto.node_mutations))

    if rsnp:
        logger.info("remove snp from clades")
        for clade, mut in tqdm.tqdm(mutations.items()):
            mut_copy = copy.deepcopy(mut.mutation)
            for mutation in mut_copy:
                mut.mutation.remove(mutation)

    if ADD_INDEL:
        logger.info("add indels in clades")
        leaf_count = 0
        inter_leaf_count = 0
        # for clade, mut in tqdm.tqdm(mutations.items()):
        #     if len(clade.clades) == 0:
        #         if not(clade.name in clade_name):
        #             leaf_count += 1
        #             clade_name.append(clade.name)
        #         else:
        #             inter_leaf_count -= 1
        #     else:
        #         for c_ in clade.clades:
        #             if len(c_.clades) == 0 and not (c_.name in clade_name):
        #                 clade_name.append(c_.name)
        #                 leaf_count +=1
        #                 inter_leaf_count += 1
        addindel_debug = "/home/yanhongliang/projects/pyro-cov/debug/addindel_debug.txt"
        with open(indel_file,"rb")as indel_f, open(addindel_debug,"w") as debug_file:
            indel_dict = pickle.load(indel_f)
            debug_parttn = "seq id: {}\n"
            for clade, mut in tqdm.tqdm(mutations.items()):
                # if rsnp:
                #     mut_copy = copy.deepcopy(mut.mutation)
                #     for mutation in mut_copy:
                #         mut.mutation.remove(mutation)
                if len(clade.clades) == 0:
                    leaf_count += 1
                    # add mutations to leaf to do: from leaf to parent
                    seq_id = clade.name
                    raw_mut_num = len(mut.mutation)
                    #if len(mutations[clade].mutation) > 0:
                    mut_tmp = parsimony_pb2.mut(position=1, ref_nuc=2, par_nuc=3, mut_nuc=[4], chromosome="NC_045512.2")
                    if not (seq_id in indel_dict):
                        debug_file.write(debug_parttn.format(seq_id))
                        continue
                    indels = indel_dict[seq_id]
                    for indel in indels:
                        indel_list = parse_indel(indel,dkind)
                        for indel_ in indel_list:
                            mut_tmp.position, mut_tmp.ref_nuc, mut_tmp.par_nuc, mut_tmp.mut_nuc[:] = indel_
                            mut.mutation.append(mut_tmp)#add indels
                    new_mut_num = len(mut.mutation)
                    if DEBUG:
                        logger.info(f"seq: {seq_id} indels {indels} #raw {raw_mut_num} #new {new_mut_num} indel {indel_}")

        # depth = len(clades[0].clades) + 1
        logger.info("merge muts in the tree")
        raw_mut = len(mutations[clades[0]].mutation)
        if MUTATION_KEEP:
            logger.info("keep mutation unchanged")
        root_mut = mut_merge(clades[0],mutations)
        logger.info(f"root: {clades[0]} before merge {raw_mut} after merge {len(root_mut)}")
        logger.info(f"leaf count {leaf_count} inter_leaf_count {inter_leaf_count}")

    # Add refined clades, collapsing clones.
    num_children: Dict[str, int] = defaultdict(int)
    clade_to_fine = {clades[0]: "fine"}
    fine_to_clade = {"fine": clades[0]}
    for parent in clades:
        parent_fine = clade_to_fine[parent]
        for child in parent.clades:
            if mutations[child].mutation:
                # Create a new fine id.
                n = num_children[parent_fine]
                fine = f"{parent_fine}.{n - 1}" if n else parent_fine + "."
                num_children[parent_fine] += 1
                clade_to_fine[child] = fine
                fine_to_clade[fine] = child
            else:
                # Collapse clone into parent.
                clade_to_fine[child] = parent_fine

    # Save basal fine clades and the fine -> coarse mapping.
    fine_to_coarse = {}
    with open("debug/refine_tree", "w") as refine_f:
        for i, (clade, meta) in enumerate(metadata.items()):
            fine = clade_to_fine[clade]
            if meta.clade and pangolin.is_pango_lineage(meta.clade):
                if DEBUG:
                    refine_f.write(f"{i}th clade {clade} fine {fine} meta.clade {meta.clade}\n")
                fine_to_coarse[fine] = pangolin.compress(meta.clade)
            meta.clade = fine if clade is fine_to_clade[fine] else ""
        if DEBUG:
            refine_f.write(f"key in fine_to_coarse {list(fine_to_coarse.keys())}")
    # Propagate basal clade metadata downward.
    for parent in clades:
        parent_coarse = fine_to_coarse[clade_to_fine[parent]]
        for child in parent.clades:
            fine_to_coarse.setdefault(clade_to_fine[child], parent_coarse)

    with open(filename_out.format(max_num_clades,dkind,rsnp,aa2muc,ADD_INDEL,res_thresh), "wb") as f:
        f.write(proto.SerializeToString())

    logger.info(f"Found {len(clades) - len(fine_to_coarse)} clones")
    logger.info(f"Refined {len(set(fine_to_coarse.values()))} -> {len(fine_to_coarse)}")
    return fine_to_coarse


def prune_mutation_tree(
    filename_in: str,
    filename_out: str,
    max_num_nodes: int,
    weights: Optional[Dict[str, int]] = None,
) -> Set[str]:
    """
    Condenses a mutation tree by greedily pruning nodes with least value
    under the error-minimizing objective function::

        value(node) = num_mutations(node) * weights(node)

    Returns a restricted set of clade names.
    """
    proto, tree = load_proto(filename_in)
    num_pruned = len(proto.node_mutations) - max_num_nodes
    logger.info(f"node mutation {len(proto.node_mutations)} max node {max_num_nodes} num_pruned {num_pruned}")
    # assert proto.node_metadata
    # exit()
    if num_pruned < 0:
        shutil.copyfile(filename_in, filename_out)
        return {m.clade for m in proto.node_metadata if m.clade}
        # return {m.clade for m in proto.metadata if m.clade}
    # logger.info("run succeed")
    # exit()
    # Extract phylogenetic tree.

    clades = list(tree.find_clades())
    logger.info(f"Pruning {num_pruned}/{len(clades)} nodes")
    assert len(clades) == len(set(clades))
    clade_to_id = {c: i for i, c in enumerate(clades)}
    assert len(proto.metadata) == len(clades)
    assert len(proto.node_mutations) == len(clades)
    metadata = dict(zip(clades, proto.metadata))
    mutations = dict(zip(clades, proto.node_mutations))
    name_set = {m.clade for m in proto.metadata if m.clade}

    # Initialize weights and topology.
    if weights is None:
        weights = {c: 1 for c in clades}
    else:
        assert set(weights).issubset(name_set)
        old_weights = weights.copy()
        weights = {}
        for c, m in metadata.items():
            weights[c] = old_weights.pop(m.clade, 0) if m.clade else 0
        assert not old_weights, list(old_weights)
    parents = {c: parent for parent in clades for c in parent.clades}
    assert tree.root not in parents

    def get_loss(clade):
        return weights[clade] * len(mutations[clade].mutation)

    # Greedily prune nodes.
    heap = [(get_loss(c), clade_to_id[c]) for c in clades[1:]]  # don't prune the root
    heapq.heapify(heap)
    for step in tqdm.tqdm(range(num_pruned)):
        # Find the clade with lowest loss.
        stale_loss, i = heapq.heappop(heap)
        clade = clades[i]
        loss = get_loss(clade)
        while loss != stale_loss:
            # Reinsert clades whose loss was stale.
            stale_loss, i = heapq.heappushpop(heap, (loss, i))
            clade = clades[i]
            loss = get_loss(clade)

        # Prune this clade.
        parent = parents.pop(clade)
        weights[parent] += weights.pop(clade, 0)  # makes the parent loss stale
        parent.clades.remove(clade)
        parent.clades.extend(clade.clades)
        mutation = list(mutations.pop(clade).mutation)
        for child in clade.clades:
            parents[child] = parent
            m = mutations[child].mutation
            cat = mutation + list(m)  # order so as to be compatible with reversions
            del m[:]
            m.extend(cat)
    clades = list(tree.find_clades())
    assert len(clades) == max_num_nodes

    # Create the pruned proto.
    proto.newick = next(iter(Writer([tree]).to_strings()))
    del proto.metadata[:]
    del proto.node_mutations[:]
    proto.metadata.extend(metadata[clade] for clade in clades)
    proto.node_mutations.extend(mutations[clade] for clade in clades)
    with open(filename_out, "wb") as f:
        f.write(proto.SerializeToString())

    return {metadata[clade].clade for clade in clades if metadata[clade].clade}


def apply_mutations(ref: str, mutations: FrozenSet[Mutation]) -> str:
    """
    Applies a set of mutations to a reference sequence.
    """
    seq = list(ref)
    for m in mutations:
        if m.mut == m.ref:
            continue
        if m.ref != seq[m.position - 1]:
            warnings.warn(f"invalid reference: {m.ref} vs {seq[m.position - 1]}")
        seq[m.position - 1] = m.mut
    return "".join(seq)


class FineToMeso:
    """
    Mapping from fine clade names like ``fine.1...3.`` to ancestors in
    ``meso_set`` like ``fine.1..`` .
    """

    def __init__(self, meso_set):
        self.meso_set = frozenset(meso_set)
        self._cache = {}

    def __call__(self, fine):
        meso = self._cache.get(fine, None)
        if meso is None:
            meso = fine if fine in self.meso_set else self(fine.rsplit(".", 1)[0])
            self._cache[fine] = meso
        return meso

