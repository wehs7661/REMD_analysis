#!/usr/bin/env python
"""
When running a replica exchange molecular dynamics (REMD) simulation, making each
replica start with the configuratin at the corresponding lambda values would shorten
the correlation time. For example, the initial configuration for the simulation of 
replica 40 should be a configuration having lambda values corresponding to replica 40.
Given the log file from a simulation which samples different intermediate states (EXE
or REMD simulation) and the corresponding tpr, xtc files, this Python code is able
to dump the configurations needed for a new REMD simulation.
"""

import argparse
import os
import sys
import copy
import numpy as np
import natsort


def initialize(args):
    parser = argparse.ArgumentParser(
        description='Extract .gro files at certain time frames from one trajectory.')
    parser.add_argument('-m',
                        '--method',
                        choices={'REMD', 'EXE'},
                        required=True,
                        help='Whether to extract configurations from an EXE simulation \
                            or an REMD simulation. Available options: REMD or EXE.')
    parser.add_argument('-nm',
                        '--mdp',
                        help='The .mdp file to use in the new REMD simulation.')
    parser.add_argument('-d',
                        '--dir',
                        required='REMD' in sys.argv[1:],  # True if 'REMD' is chosen (argv[0] is the script name)
                        help='The directory of where stateX folders are. This argument is only needed \
                            for the case of extracting configurations from an REMD simulation so\
                            that the program is able to find required files (.tpr, .trr and .log).')
    parser.add_argument('-l',
                        '--log',
                        nargs='+',
                        help='The file name of the log file of the formaer EXE or REMD \
                            simulation. Note that if REMD method is chosen, a log file \
                            from any replica should work and one might need multiple log\
                            files if the simulation was extended.')
    parser.add_argument('-f',
                        '--trj',
                        nargs='+',
                        help='The file name of the trajectory file(s) (.trr or .xtc). \
                            of either an EXE or REMD simulaiton. If the REMD method is \
                            chosen, multiple trjaectory files might be required (one for \
                            each replica), but the names of the trajectory files of different \
                            replicas are assumed to be the same.')
    parser.add_argument('-s',
                        '--tpr',
                        nargs='+',
                        help='The file name of the structure file (.tpr). of either an \
                            EXE or REMD simulation. If the REMD method is chosen, multiple \
                            trajectory files might be required (one fore each replica), but \
                            the names of the .tpr files of different replicas are assumed to \
                            be the same.')
    parser.add_argument('-p',
                        '--prefix',
                        help='The common prefix of the simulation files.')
    parser.add_argument('-M',
                        '--MPI',
                        default=False,
                        action='store_true',
                        help='Whether the MPI-version of GROMACS is used. "-m" \
                            specified means that MPI is used.')

    args_parse = parser.parse_args(args)

    # Below we deal with the case that any argument is not specified.

    # Try to find the mdp file in the current directory if not specified
    if args_parse.mdp is None:
        for file in os.listdir('.'):
            if file.endswith('.mdp'):
                args_parse.mdp = file
        if args_parse.mdp is None:
            print('Error: No mdp file for the new REMD simulation provided or found!')
            sys.exit()

    # The case that EXE method is chosen: search all the files in the current directory
    if args_parse.method == 'EXE':
        if args_parse.log is None:
            for file in os.listdir('.'):
                if file.endswith('.log'):
                    args_parse.log = file
            if args_parse.log is None:
                print('Error: No log file of an EXE or REMD simulation provided or found!')
                sys.exit()
        if args_parse.prefix is None:
            args_parse.prefix = args_parse.log.split('.')[0]  
        if args_parse.trj is None:
            if os.path.isfile('%s.trr' % args_parse.prefix) is True:
                args_parse.trj = '%s.trr' % args_parse.prefix
            elif os.path.isfile('%s.xtc' % args_parse.prefix) is True:
                args_parse.trj = '%s.xtc' % args_parse.prefix
            else:
                print('Error: No trajectory file (.xtc or .trr file) of the EXE simulation provided or found!')
                sys.exit()
        if args_parse.tpr is None:
            if os.path.isfile('%s.tpr' % args_parse.prefix) is True:
                args_parse.tpr = '%s.tpr' % args_parse.prefix
            else:
                print('Error: No structure file (.tpr file) provided or found!')
                sys.exit()

    # The case that REMD method is chosen: search the directory specified 
    if args_parse.method == 'REMD':
        # Note that here we assume that the file names of the trajectory files and the .tpr files
        # of different replicas of REMD are the same.
        if args_parse.log is None:
            # log file from any replica should be fine, so we'll just refer to the one(s) in state_0
            args_parse.log = find_REMD_files(args_parse.dir + '/state_0', '.log')
        if args_parse.trj is None:
            args_parse.trj = find_REMD_files(args_parse.dir + '/state_0', ['.trr', '.xtc'])
        if args_parse.tpr is None:
            args_parse.tpr = find_REMD_files(args_parse.dir + '/state_0', '.tpr')
        
        if args_parse.prefix is None:
            if isinstance(args_parse.log, list):
                args_parse.prefix = args_parse.log[0].split('.')[0]
            else:
                args_parse.prefix = args_parse.log.split('.')[0]

    return args_parse

def find_REMD_files(file_dir, ext):
    # file_dir: args.dir + 'state_0'
    # ext: file extension, '.log' for example
    # find all the files with the same extension in dir/state_X
    file_to_find = []

    for file in os.listdir(file_dir):
        if isinstance(ext, list):
            for i in ext:
                if file.endswith(i):
                    file_to_find.append(file)
        else:  # then should be a str
            if file.endswith(ext):
                file_to_find.append(file)
    
    if not file_to_find:
        print('Err:r No %s files provided or found!' % ext)
        sys.exit()
    else:
        file_to_find = natsort.natsorted(file_to_find, reverse=False)
        if '.xtc' in file_to_find[-1]:
            # for some reason, natsorted put filename.xtc behind filename.part0005.xtc
            # so we have to deal with this (no problems with other extension)
            if '00' not in file_to_find[-1]: 
                last_file = file_to_find[-1]  
                file_to_find.remove(last_file)     # remove the last element
                file_to_find.insert(0, last_file)  # insert back to the beginning

    return file_to_find


class ExtractState:
    def __init__(self, mdp, log):
        """
        This function get the information of the new lambda values by parsing the 
        mdp file for the new REMD simulation and the old lambda values by parsing
        the log file of the old REMD or EXE simulation. Note thate in the case that
        the configurations are extracted from an REMD simulation that had been
        extended, there might be multiple log file. If this is the case, we will
        just pick the last log file. Attributes like log_lambda and mdp_lambda will
        be added in the init function.
        """
        
        # info in mdp: for the new simulation
        self.mdp_lambda = []   # self.mdp_lambda = [coul[i], vdw[i], restr[i]]
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

        self.mdp_N = len(mdp_coul)  # number of states of the new REMD simulation
        if len(mdp_vdw) != self.mdp_N or len(mdp_restr) != self.mdp_N:
            print('Error: The number of lambda values for differnt interactiosn are not the same!')
            sys.exit()

        # combine all the lambda values in self.mdp_lambda

        for i in range(self.mdp_N):
            self.mdp_lambda.append([])
            self.mdp_lambda[i] = [mdp_coul[i], mdp_vdw[i], mdp_restr[i]]


        # info in log: out of the exisiting simulation
        self.log_lambda = []   # self.log_lambda = [coul[i], vdw[i], restr[i]]
        
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

    def get_extract_state(self):
        """
        This function compare the new and old lambda values to decide from which 
        states should extract to serve as the starting configurations of the new 
        REMD. Note that we always extract the structure at the last time frame, 
        which is the safest way if we want to ensure that the system has equilibrated.
        This function can be used for both EXE method and REMD method.
        """
        extract_state = np.zeros(self.mdp_N)
        # extract_state[i]: the configuration i was extracted from the state i in the exisiting simulation
        
        for i in range(len(self.mdp_lambda)):
            if self.mdp_lambda[i] in self.log_lambda:  # the states that remain the same
                extract_state[i] = self.log_lambda.index(self.mdp_lambda[i])
            else:
                if self.mdp_lambda[i][0] in self.log_coul:
                    # the case that coul lambda stays the same while other lambdas (basically lambda_restr) don't
                    extract_state[i] = self.log_coul.index(self.mdp_lambda[i][0])
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

        # extract_state = np.array(extract_state) + 1  # make the state index starts from 1 instead of 0

        return extract_state

    def get_extract_time(self, method, log, extract_state):
        # extract_time[i]: the time frame from which the confiuration i was extracted
        # Note that the list extract_state might look like:
        # [0, 0, 0, 1, 1,2, 2, 3, 4, 5, 5, 5, 6, ...]
        # So we will create a diction like {n: number of times that n repeats}
        # like: {1: 2, 2: 2, 3: 1, 4:1, 5: 3, ...}
        # The function will make sure that all the confiuguration are extracted from different time frame
        # even if the configurations are at the same states.

        if method == 'EXE':
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
        
        elif method == 'REMD':
            f = open(log, 'r')
            lines = f.readlines()
            f.close()

            lines.reverse()  # from this point, lines have been reversed
            for l in lines:
                # print(l)   # this will print from the bottom
                if 'Replica exchange at step' in l:
                    extract_time = l.split('time')[1]
                    break
            extract_time = int(extract_time.split('.')[0]) - 5  
            # make sure that there is a structure able to be extracted

        return extract_time


def main():
    args = initialize(sys.argv[1:])

    if args.method == 'EXE':
        print("Chosen method: extracting configurations from an EXE simulation.")
        print("EXE simulation output files from which the configurations will be extracted:")
        print("%s, %s, and %s." % (args.log, args.tpr, args.trj) )
        print("The .mdp file for the new REMD simulation is: %s\n" % args.mdp)
        print("Initializing the preparation of the .gro files ...")
        ES = ExtractState(args.mdp, args.log)
        extract_state = ES.get_extract_state()
        print("Identifying the time frames from which the configurations will be extracted ...")
        extract_time = ES.get_extract_time(args.method, args.log, extract_state)
        print('extract_time = ', extract_time)
        print("Start extracting configurations from the trajectory of the provided EXE simulation ...")
        
        for i in range(len(extract_time)):
            # note that we use os.popen intead of os.system to use pip to answer the prompts
            # raised by gromacs
            os.system("mkdir state_%s" % i)
            if args.MPI:
                p = os.popen("gmx_mpi trjconv -f %s -s %s -o state_%s/%s.gro -dump %s"
                            % (args.trj, args.tpr, i, args.prefix, extract_time[i]), 'w')
                p.write("0")

            if not args.MPI:
                p = os.popen("gmx trjconv -f %s -s %s -o state_%s/%s.gro -dump %s"
                            % (args.trj, args.tpr, i, args.prefix, extract_time[i]), 'w')
                p.write("0")  # 0: System

    elif args.method == 'REMD':
        print("Chosen method: extracting configurations from an REMD simulation.")
        print("The .mdp file for the new REMD simulation is: %s\n" % args.mdp)
        print("Initializing the preparation of the .gro files ...")
        log =  args.dir + 'state_0/' + args.log[-1]
        ES = ExtractState(args.mdp,log)
        extract_state = ES.get_extract_state() 
        print('The .gro files are going to be extracted from the following states of the provided REMD:')
        print('extract_state = ', extract_state)
        extract_time = ES.get_extract_time(args.method, log, extract_state)
        print('In the case of extracting configurations from an REMD simulation, we always extract from the last time frame.')
        print("Start extracting configurations from the trajectory of the provided REMD simulation ...")

        for i in range(len(extract_state)):
            s = int(extract_state[i])
            # note that we use os.popen intead of os.system to use pip to answer the prompts
            # raised by gromacs
            os.system("mkdir state_%s" % i)
            if args.MPI:
                print("Executing the command: gmx_mpi trjconv -f %sstate_%s/%s -s %sstate_%s/%s -o state_%s/%s.gro -dump %s" % (args.dir, s, args.trj[-1], args.dir, s, args.tpr[-1], i, args.prefix, extract_time))
                p = os.popen("gmx_mpi trjconv -f %sstate_%s/%s -s %sstate_%s/%s -o state_%s/%s.gro -dump %s"
                            % (args.dir, s, args.trj[-1], args.dir, s, args.tpr[-1], i, args.prefix, extract_time), 'w')
                p.write("0")   # 0: System
            if not args.MPI:
                print("Executing the command: gmx trjconv -f %sstate_%s/%s -s %sstate_%s/%s -o state_%s/%s.gro -dump %s" % (args.dir, s, args.trj[-1], args.dir, s, args.tpr[-1], i, args.prefix, extract_time))
                p = os.popen("gmx trjconv -f %sstate_%s/%s -s %sstate_%s/%s -o state_%s/%s.gro -dump %s"
                            % (args.dir, s, args.trj[-1], args.dir, s, args.tpr[-1], i, args.prefix, extract_time), 'w')
                p.write("0")   # 0: System

    # Might not reach here since GROMCAS does not send signal to Python after 
    # the trjconv commands are done. Press Ctrl + C to exit. (The work can be properly done though.)
    print("Complete!")  
