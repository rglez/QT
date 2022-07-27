#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
.. note ::

  | **Author      :** Roy Gonzalez Aleman
  | **Contact     :** [roy_gonzalez@fq.uh.cu, roy.gonzalez.aleman@gmail.com]
'''
import sys
import argparse
import numpy as np
import numpy.ma as ma

import mdtraj as md

# =============================================================================
# Useful functions
# =============================================================================


def parse_arguments():
    '''
    DESCRIPTION
    Parse all user arguments from the command line.

    Return:
        user_inputs (parser.argparse): namespace with user input arguments.
    '''

    # Initializing argparse ---------------------------------------------------
    desc = '\nQT: Implementation of the Quality Threshold Clustering algorithm by Heyer et. al.'
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
    # Arguments: clustering ---------------------------------------------------
    parser.add_argument('-cutoff', action='store', dest='cutoff',
                        help='RMSD cutoff for pairwise comparisons in A',
                        type=float, required=False, default=1.0)
    parser.add_argument('-minsize', action='store', dest='minsize',
                        help='minimum number of frames inside returned clusters',
                        type=int, required=False, default=2)
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

    # Reduce RAM consumption by loading selected atoms only -------------------
    if args.selection != 'all':
        try:
            sel_indx = trajectory.topology.select(args.selection)
        except ValueError:
            print('Specified selection is invalid')
            sys.exit()
        if sel_indx.size == 0:
            print('Specified selection in your system corresponds to no atoms')
            sys.exit()
        trajectory = trajectory.atom_slice(sel_indx)[args.first:args.last:args.stride]
    else:
        trajectory = trajectory[args.first:args.last:args.stride]

    # Center coordinates of loaded trajectory ---------------------------------
    trajectory.center_coordinates()
    return trajectory


# ---- Check that trajectory exists inside inputs -----------------------------


# =============================================================================
# Preprocessing
# =============================================================================

# ---- Load trajectory --------------------------------------------------------
inputs = parse_arguments()
sms = '\n\n ATTENTION !!! No trajectory passed.Run with -h for help.'
assert inputs.trajectory, sms
trajectory = load_trajectory(inputs)
N = trajectory.n_frames

# ---- Calculate matrix pairwise distances ------------------------------------
matrix = np.ndarray((N, N), dtype=np.float16)
for i in range(N):
    rmsd_ = md.rmsd(trajectory, trajectory, i, precentered=True)*10
    matrix[i] = rmsd_
print('>>> Calculation of the RMSD matrix completed <<<')

# ---- Delete unuseful values from matrix (diagonal &  x>threshold) -----------
matrix[matrix > inputs.cutoff] = np.inf
matrix[matrix == 0] = np.inf
degrees = (matrix < np.inf).sum(axis=0)

# =============================================================================
# QT algotithm
# =============================================================================

clusters_arr = np.ndarray(N, dtype=np.int64)
clusters_arr.fill(-1)

ncluster = 0
while True:
    # This while executes for every cluster in trajectory ---------------------
    len_precluster = 0
    while True:
        # This while executes for every potential cluster analyzed ------------
        biggest_node = degrees.argmax()
        precluster = []
        precluster.append(biggest_node)
        candidates = np.where(matrix[biggest_node] < np.inf)[0]
        next_ = biggest_node
        distances = matrix[next_][candidates]
        while True:
            # This while executes for every node of a potential cluster -------
            next_ = candidates[distances.argmin()]
            precluster.append(next_)
            post_distances = matrix[next_][candidates]
            mask = post_distances > distances
            distances[mask] = post_distances[mask]
            if (distances == np.inf).all():
                break
        degrees[biggest_node] = 0
        # This section saves the maximum cluster found so far -----------------
        if len(precluster) > len_precluster:
            len_precluster = len(precluster)
            max_precluster = precluster
            max_node = biggest_node
            degrees = ma.masked_less(degrees, len_precluster)
        if not degrees.max():
            break
    # General break if minsize is reached -------------------------------------
    if len(max_precluster) < inputs.minsize:
        break

    # ---- Store cluster frames -----------------------------------------------
    clusters_arr[max_precluster] = ncluster
    ncluster += 1
    print('>>> Cluster # {} found with {} frames at center {} <<<'.format(
          ncluster, len_precluster, max_node))

    # ---- Update matrix & degrees (discard found clusters) -------------------
    matrix[max_precluster, :] = np.inf
    matrix[:, max_precluster] = np.inf

    degrees = (matrix < np.inf).sum(axis=0)
    if (degrees == 0).all():
        break

# ---- Write clustering results to disk ---------------------------------------

# simple format
np.savetxt('QT_Clusters.txt', clusters_arr, fmt='%i')

# NMRcluster format. VMD interface
with open('QT_Visualization.log', 'wt') as clq:
    for numcluster in np.unique(clusters_arr):
        clq.write('{}:\n'.format(numcluster))
        members = ' '.join([str(x + 1)
                            for x in np.where(clusters_arr == numcluster)[0]])
        clq.write('Members: ' + members + '\n\n')
