#!/usr/bin/env python
"""
When running a replica exchange molecular dynamics (REMD) simulation, making each
replica start with the configuratin at the corresponding lambda values would shorten
the correlation time. For example, the initial configuration for the simulation of 
replica 40 should be a configuration having lambda values corresponding to replica 40.
Given the log file from a simulation which samples different intermediate states (EXE
or REMD simulation) and the corresponding tpr, xtc files, this Python code is able
to dump the configurations needed for a REMD simulation.
"""

import argparse
import os
import sys
import copy
import numpy as np


def initialize(args):
    parser = argparse.ArgumentParser(
        description='Extract .gro files at certain time frames from one trajectory.')
    parser.add_argument('-nm',
                        '--mdp',
                        help='The .mdp file to use in the REMD simulation.')
    parser.add_argument('-l',
                        '--log',
                        help='The file name of the log file.')
    parser.add_argument('-f',
                        '--trj',
                        help='The file name of the trajectory file (.trr or .xtc).')
    parser.add_argument('-s',
                        '--tpr',
                        help='The file name of the structure file (.tpr).')
    parser.add_argument('-p',
                        '--prefix',
                        help='The common prefix of the simulation files.')
    parser.add_argument('-m',
                        '--mpi',
                        default=False,
                        action='store_true',
                        help='Whether the MPI-version of GROMACS is used. "-m" \
                            specified means that MPI is used.')

    args_parse = parser.parse_args(args)

    if args_parse.mdp is None:
        for file in os.listdir('.'):
            if file.endswith('.mdp'):
                args_parse.mdp = file
        if args_parse.mdp is None:
            print('Error: No log file provided or found!')
            sys.exit()

    if args_parse.log is None:
        for file in os.listdir('.'):
            if file.endswith('.log'):
                args_parse.log = file
        if args_parse.log is None:
            print('Error: No log file provided or found!')
            sys.exit()

    if args_parse.prefix is None:
        args_parse.prefix = args_parse.log.split('.')[0]

    if args_parse.trj is None:
        if os.path.isfile('%s.trr' % args_parse.prefix) is True:
            args_parse.trj = '%s.trr' % args_parse.prefix
        elif os.path.isfile('%s.xtc' % args_parse.prefix) is True:
            args_parse.trj = '%s.xtc' % args_parse.prefix
        else:
            print('Error: No trajectory file (.xtc or .trr file) provided or found!')
            sys.exit()

    if args_parse.tpr is None:
        if os.path.isfile('%s.tpr' % args_parse.prefix) is True:
            args_parse.tpr = '%s.tpr' % args_parse.prefix
        else:
            print('Error: No structure file (.tpr file) provided or found!')
            sys.exit()

    return args_parse


class ExtractState:
    def __init__(self, mdp, log):
        # self.lambda = [coul[i], vdw[i], restr[i]]
        self.mdp_lambda, self.log_lambda = [], []

        # info in mdp: for the new simulation
        f = open(mdp, 'r')
        lines = f.readlines()
        f.close()

        for l in lines:
            if 'coul-lambdas' in l:
                mdp_coul = [float(i) for i in l.split('=')[1].split()]
            if 'vdw-lambdas' in l:
                mdp_vdw = [float(i) for i in l.split('=')[1].split()]
            if 'restraint-lambdas' in l:
                mdp_restr = [float(i) for i in l.split('=')[1].split()]

        self.mdp_N = len(mdp_coul)

        # combine all the lambda values in self.mdp_lambda

        for i in range(self.mdp_N):
            self.mdp_lambda.append([])
            self.mdp_lambda[i] = [mdp_coul[i], mdp_vdw[i], mdp_restr[i]]

        # info in log: out of the exisiting simulation
        f = open(log, 'r')
        lines = f.readlines()
        f.close()

        for l in lines:   # find the information from the top of the log file
            if 'coul-lambdas = ' in l:
                if len(l.split('=')[1].split()) > 1:  # to prevent catching 'coul-lambdas = TRUE'
                    self.log_coul = [float(i) for i in l.split('=')[1].split()]

            if 'vdw-lambdas = ' in l:
                if len(l.split('=')[1].split()) > 1:  # to prevent catching 'vdw-lambdas = TRUE'
                    self.log_vdw = [float(i) for i in l.split('=')[1].split()]

            if 'restraint-lambdas = ' in l:
                if len(l.split('=')[1].split()) > 1:  # to prevent catching 'restraint-lambdas = TRUE'
                    self.log_restr = [float(i) for i in l.split('=')[1].split()]

            if 'Started mdrun' in l:
                break

        self.log_N = len(self.log_coul)

        for i in range(self.log_N):
            self.log_lambda.append([])
            self.log_lambda[i] = [self.log_coul[i], self.log_vdw[i], self.log_restr[i]]

    def get_extract_state(self, mdp, log):
        extract_state = np.zeros(self.mdp_N)
        # extract_state[i]: the configuration i was extracted from the state i in the exisiting simulation

        for i in range(len(self.mdp_lambda)):
            if self.mdp_lambda[i] in self.log_lambda:
                extract_state[i] = self.log_lambda.index(self.mdp_lambda[i])
            else:
                if self.mdp_lambda[i][0] < 1 and self.mdp_lambda[i][1] == 0:   # coul region
                    for j in range(len(self.log_coul) - 1):     # to prevent index out of range
                        if self.mdp_lambda[i][0] > self.log_coul[j] and self.mdp_lambda[i][0] < self.log_coul[j + 1]:
                            if j > len(self.log_coul) / 2:  # closer to uncoupled state
                                extract_state[i] = j + 1
                            else:                           # closer to coupled state
                                extract_state[i] = j

                if self.mdp_lambda[i][0] == 1 and self.mdp_lambda[i][1] < 1:   # vdw region
                    for j in range(len(self.log_vdw) - 1):     # to prevent index out of range
                        if self.mdp_lambda[i][1] > self.log_vdw[j] and self.mdp_lambda[i][1] < self.log_vdw[j + 1]:
                            if j > len(self.log_vdw) / 2:   # closer to uncoupled state
                                extract_state[i] = j + 1
                            else:                           # closer to coupled state
                                extract_state[i] = j

        extract_state = np.array(extract_state) + 1  # make the state index starts from 1 instead of 0

        return extract_state

    def get_extract_time(self, log, extract_state):
        # extract_time[i]: the time frame from which the confiuration i was extracted
        # Note that the list extract_state might look like:
        # [1, 1, ,2, 2, 3, 4, 5, 5, 5, 6, ...]
        # So we will create a diction like {n: number of times that n repeats}
        # like: {1: 2, 2: 2, 3: 1, 4:1, 5: 3, ...}
        # The function will make sure that all the confiuguration are extracted from different time frame
        # even if the configurations are at the same states.
        state_dict = {}   # key: state to find
        for i in range(len(extract_state)):
            if extract_state[i] not in state_dict:
                state_dict[extract_state[i]] = 1
            else:
                state_dict[extract_state[i]] += 1

        found_dict = copy.copy(state_dict)   # for checking if each configuration has been found
        for i in found_dict:
            found_dict[i] = 0

        f = open(log, 'r')
        lines = f.readlines()
        f.close()

        line_n = 0
        extract_time = np.zeros(len(extract_state))

        for l in lines:
            line_n += 1
            if 'MC-lambda information' in l:
                for i in range(self.log_N):
                    # the lambdas value start from lines[line_n + 2]
                    # 'MC-lambda information' is lines[line_n - 1]
                    # the time frame is in lines[line_n - 3]  (units: ps)
                    if lines[line_n + 2 + i].split()[-1] == '<<':
                        N = int(lines[line_n + 2 + i].split()[0])
                        if N in state_dict and found_dict[N] < state_dict[N]:
                            N_idx = list(extract_state).index(N) + found_dict[N]
                            # for example, for the second time that state N is found
                            # N_idx = extract_state.index(N_idx) + 1
                            extract_time[N_idx] = float(lines[line_n - 3].split()[1])   # (units: ps)
                            found_dict[N] += 1   # should be after extract_time[N_idx]
                    if found_dict == state_dict:
                        # all the configurations have been found
                        break

        return extract_time


def main():
    args = initialize(sys.argv[1:])
    print("EXE simulation output files from which the configurations will be extracted:")
    print("%s, %s, and %s." % (args.log, args.tpr, args.trj) )
    print("The .mdp file for the REMD simulation is : %s\n" % args.mdp)
    print("Initializing the preparation of the .gro files...")
    ES = ExtractState(args.mdp, args.log)
    print("Identifying the time frames from which the configurations will be extracted ...")
    extract_state = ES.get_extract_state(args.mdp, args.log)
    extract_time = ES.get_extract_time(args.log, extract_state)
    print("Start extracting configurations from the trajectory of the provided EXE simulation ...")
    
    for i in range(len(extract_time)):
        # note that we use os.popen intead of os.system to use pip to answer the prompts
        # raised by gromacs
        os.system("mkdir state_%s" % i)
        if args.mpi:
            p = os.popen("gmx_mpi trjconv -f %s -s %s -o state_%s/%s.gro -dump %s"
                         % (args.trj, args.tpr, i, args.prefix, extract_time[i]), 'w')
            p.write("0")

        if not args.mpi:
            p = os.popen("gmx trjconv -f %s -s %s -o state_%s/%s.gro -dump %s"
                         % (args.trj, args.tpr, i, args.prefix, extract_time[i]), 'w')
            p.write("0")


    # Might not be reached here since GROMCAS does not send signal to Python after 
    # the trjconv commands are done. Press Ctrl + C to exit. (The work can be properly done though.)
    print("Complete!")  
