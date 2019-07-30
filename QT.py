#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
.. note ::

  | **Created on  :** Fri Apr 26 19:03:42 2019
  | **Author      :** Roy Gonzalez Aleman
  | **Contact     :** [roy_gonzalez@fq.uh.cu, roy.gonzalez.aleman@gmail.com]
'''
# import glob
import argparse
import numpy as np
from copy import deepcopy as dc
from bisect import bisect

# import matplotlib as mpl
# import matplotlib.pyplot as plt
# from matplotlib.ticker import MaxNLocator

import mdtraj as md


def parse_arguments():
    '''
    DESCRIPTION
    Parse all user arguments from the command line.

    Return:
        user_inputs (parser.argparse): namespace with user input arguments.
    '''

    # Initializing argparse ---------------------------------------------------
    desc = '\nQT: Quality Threshold Clustering algorithm by Heyer et. al.'
    parser = argparse.ArgumentParser(description=desc,
                                     add_help=True,
                                     epilog='As simple as that ;)')
    # Arguments: loading trajectory -------------------------------------------
    parser.add_argument('-top', dest='topology', action='store',
                        help='path to topology file (psf/pdb)',
                        type=str, required=False)
    parser.add_argument('-traj', dest='trajectory', action='store',
                        help='path to trajectory file',
                        type=str)
    parser.add_argument('-first', dest='first',  action='store',
                        help='first frame to analyze (starting from 0)',
                        type=int, required=False, default=0)
    parser.add_argument('-last', dest='last', action='store',
                        help='last frame to analyze (starting from 0)',
                        type=int, required=False, default=-1)
    parser.add_argument('-stride', dest='stride', action='store',
                        help='stride of frames to analyze',
                        type=int, required=False, default=1)
    parser.add_argument('-sel', dest='selection', action='store',
                        help='atom selection (MDTraj syntax)',
                        type=str, required=False, default='all')
    parser.add_argument('-rmwat', dest='remove_waters', action='store',
                        help='remove waters from trajectory?',
                        type=bool, required=False, default=0,
                        choices=[True, False])
    # Arguments: clustering ---------------------------------------------------
    parser.add_argument('-cutoff', action='store', dest='cutoff',
                        help='RMSD cutoff for pairwise comparisons in A',
                        type=float, required=False, default=1.0)
    parser.add_argument('-minsize', action='store', dest='minsize',
                        help='minimum number of frames inside returned clusters',
                        type=int, required=False, default=2)
    parser.add_argument('-ref', action='store', dest='reference',
                        help='reference frame to align trajectory',
                        type=int, required=False, default=0)
    # Arguments: analysis -----------------------------------------------------
    parser.add_argument('-odir', action='store', dest='outdir',
                        help='output directory to store analysis',
                        type=str, required=False, default='./')
    user_inputs = parser.parse_args()
    return user_inputs


def load_trajectory(args):
    '''
    DESCRIPTION
    Loads trajectory file using MDTraj. If trajectory format is h5, lh5 or
    pdb, topology file is not required. Otherwise, you should specify a
    topology file.

    Arguments:
        args (argparse.Namespace): user input parameters parsed by argparse.
    Return:
        trajectory (mdtraj.Trajectory): trajectory object for further analysis.
    '''

    traj_file = args.trajectory
    traj_ext = traj_file.split('.')[-1]
    # Does trajectory file format need topology ? -----------------------------
    if traj_ext in ['h5', 'lh5', 'pdb']:
        trajectory = md.load(traj_file)
    else:
        trajectory = md.load(traj_file, top=args.topology)
    # Reduce RAM consumption by loading all atoms except water ----------------
    if args.remove_waters:
        nowat_indx = trajectory.topology.select('all != water')
        nowat_traj = trajectory.restrict_atoms(nowat_indx)
        del trajectory
        trajectory = nowat_traj
    # Reduce RAM consumption by loading selected atoms only -------------------
    if args.selection != 'all':
        sel_indx = trajectory.topology.select(args.selection)
        sel_traj = trajectory.restrict_atoms(sel_indx)
        del trajectory
        trajectory = sel_traj[args.first:args.last:args.stride]
    else:
        trajectory = trajectory[args.first:args.last:args.stride]
    # Center coordinates of loaded trajectory ---------------------------------
    trajectory.center_coordinates()
    return trajectory


def find_max_kcore(degrees):
    '''
    DESCRIPTION
    Load a numpay array of degrees of the graph and does a binary search to
    return its maximum kcore, i.e. how many vertex have degree >= k.

    Args:
        degrees: numpy array of graph's degrees.
    Return:
        maximum kcore of the graph.
    '''
    sorted_degrees = dc(degrees)
    sorted_degrees.sort()
    if len(degrees) == 0:
        return 0

    new_degrees = np.asarray(range(sorted_degrees.max(), 0, -1))
    kcores = []
    while len(new_degrees) != 0:
        half = len(new_degrees)//2
        count = len(sorted_degrees) - bisect(sorted_degrees,
                                             new_degrees[half]-1)
        if count >= new_degrees[half]:
            kcores.append(new_degrees[half])
            new_degrees = new_degrees[:half:]
        elif count < new_degrees[half]:
            new_degrees = new_degrees[half+1::]
    try:
        max_kcore = max(kcores)
    except ValueError:
        max_kcore = 0
    return max_kcore


# ---- Load user arguments ----------------------------------------------------
inputs = parse_arguments()
sms = '\n\n ATTENTION !!! No trajectory passed.Run with -h for help.'
assert inputs.trajectory, sms

# ---- Load trajectory --------------------------------------------------------
trajectory = load_trajectory(inputs)
N = trajectory.n_frames

# ---- Calculate float matrix -------------------------------------------------
hard_matrix = np.ndarray((N, N), dtype=np.float16)
for i in range(N):
    rmsd_ = md.rmsd(trajectory, trajectory, i, precentered=True)*10
    hard_matrix[i] = rmsd_

# ---- QT-vectorized algorithm ------------------------------------------------
clusters = []
nclustered = 0
max_cluster = np.array(range(N+1))

# ---- Repeat for each cluster found ------------------------------------------
while nclustered < N:
    nclustered = sum([len(x) for x in clusters])
    # print(nclustered, ' frames have been clustered')
    pre_clusters = []
    soft_matrix = dc(hard_matrix)
    index_xy = np.asarray(range(N))
    max_found = 0

    # ---- Iterate until max cluster is detected for each step ----------------
    while True:
        members = []
        # Get degrees ---------------------------------------------------------
        degrees = np.asarray([(i <= inputs.cutoff).sum() for i in soft_matrix])
        # Get real degrees ----------------------------------------------------
        real_degrees = [find_max_kcore(degrees[np.where(row <= inputs.cutoff)[0]])
                        for row in soft_matrix]
        real_degrees = np.asarray(real_degrees)

        # STOP CONDITION #1:
        # When real_degrees can not form clusters bigger than max_cluster_len
        if find_max_kcore(real_degrees) <= max_found:
            # print('max k-core of real degrees are shorter than max_found')
            break

        # Get biggest node and use it as seed ---------------------------------
        if real_degrees.max() > max_found:
            biggest_node = real_degrees.argmax()
        else:
            # STOP CONDITION #2:
            # If biggest node has degree less than max_cluster_len
            # print('max degree became less than max_found')
            break
        members.append(biggest_node)
        soft_matrix[biggest_node][biggest_node] = np.inf

        # Get potential cluster -----------------------------------------------
        candidates = np.where(soft_matrix[biggest_node] <= inputs.cutoff)[0]
        while len(candidates):
            # Detect next node to enter ---------------------------------------
            distances = np.asarray([soft_matrix[candidate][members].sum()
                                    for candidate in candidates])
            min_distances = distances.argmin()
            next_one = candidates[min_distances]
            members.append(next_one)
            # Discard conflicting nodes with those already inside cluster -----
            soft_matrix[biggest_node][next_one] = np.inf
            candidates = np.delete(candidates, min_distances)
            candidates = np.intersect1d(
                    candidates, np.where(soft_matrix[next_one] <= inputs.cutoff)[0])

        # Track promising clusters only ---------------------------------------
        if len(members) > max_found:
            max_found = len(members)
            max_node = index_xy[biggest_node]
            max_cluster = index_xy[members]

        # Exclude nodes with length less than max_cluster detected ------------
        pre_clusters.append(index_xy[members])
        excluded = np.where(real_degrees <= max_found)[0]
        soft_matrix = np.delete(soft_matrix, excluded, axis=0)
        soft_matrix = np.delete(soft_matrix, excluded, axis=1)
        index_xy = np.delete(index_xy, excluded)

    # STOP CONDITION #3:
    # Length of cluster is less than user specification
    if max_cluster.size < inputs.minsize:
        break

    # ---- At max cluster detection, save it and update hard_matrix -----------
    clusters.append(max_cluster)
    hard_matrix[max_cluster] = np.inf
    print(len(max_cluster))

# ---- After all iterations, organize results ---------------------------------
clusters_arr = np.ndarray(N, dtype=np.int64)
clusters_arr.fill(-1)
for i, cluster in enumerate(clusters):
    clusters_arr[cluster] = i
np.savetxt('QT_Clusters.txt', clusters_arr, fmt='%i')
