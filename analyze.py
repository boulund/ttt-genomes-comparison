#!/usr/bin/env python
# Fredrik Boulund 2014
# Analyze pneu mitis similarity

import argparse
from sys import argv, exit
from collections import namedtuple


def parse_mappings(mappings):
    infodict ={}
    fragment_tuple = namedtuple("Fragment_tuple", ["fragment", "target", "identity", "matches", "mismatches"])
    with open(mappings) as file:
        counter = 0
        for line in file:
            fragment_info = line.split()
            fragment = fragment_info[0]
            target = fragment_info[1]
            identity = float(fragment_info[2])
            matches = int(fragment_info[3])
            mismatches = int(fragment_info[4])
            ftuple = fragment_tuple(fragment, target, identity, matches, mismatches)
            try:
                infodict[fragment].append(ftuple)
            except KeyError:
                infodict[fragment] = [ftuple]
            counter += 1

        print "Found {} matched fragments with {} total matches.".format(len(infodict.keys()), counter)

    return infodict



def filter_hits(infodict, remove_noninformative=True, print_fragment="", print_fragments=False, matches=20, identity=100, mismatches=0):
    """ Filter out fragments that match the given criteria."""
    matches = int(matches)
    mismatches = int(mismatches)
    identity = float(identity)
    genomes = set()

    # Remove noninformative fragments
    informative_dict = {}
    if remove_noninformative:
        for fragment in infodict.iterkeys():
                hit_genomes = set()
                for infotuple in infodict[fragment]:
                    hit_genomes.add(infotuple.target)
                if len(hit_genomes) == 1:
                    informative_dict[fragment] = infodict[fragment]
    else:
        informative_dict = infodict
    print "Removed {} non-informative fragments. {} fragments remain.".format(len(infodict.keys())-len(informative_dict.keys()), len(informative_dict.keys()))
    
    # Filter hits based on user critera
    filtered_dict = {}
    for fragment in informative_dict.iterkeys():
        tuplelist = []
        for infotuple in informative_dict[fragment]:
            if infotuple.identity >= identity and\
                    infotuple.matches >= matches and\
                    infotuple.mismatches <= mismatches:
                tuplelist.append(infotuple)
        if tuplelist:
            [genomes.add(infotuple.target) for infotuple in tuplelist]
            filtered_dict[fragment] = tuplelist
    print "Filtered {} fragments based on critera. {} fragments remain.".format(len(informative_dict.keys())-len(filtered_dict.keys()), len(filtered_dict.keys()))

    if print_fragment:
        print filtered_dict[print_fragment]
    if print_fragments:
        print "{:<10} {:<10} {:<6} {:<4} {:<4}".format("Fragment", "target", "ID", "M", "m")
        for fragment, infotuples  in filtered_dict.iteritems():
            for infotuple in infotuples:
                print "{:<10} {:<10} {:<6} {:<4} {:<4}".format(infotuple.fragment, infotuple.target, infotuple.identity, infotuple.matches, infotuple.mismatches)

    # Summary printout
    bactdict = {}
    for genome in genomes:
        bactdict[genome] = 0

    for fragment in filtered_dict.keys():
        for infotuple in filtered_dict[fragment]:
            bactdict[infotuple.target] += 1

    print "{:<10} {:<5} {}".format("GENOME", "#", "PERCENTAGE")
    for bact in bactdict.keys():
        totalhits = float(sum(bactdict.values()))
        if totalhits == 0:
            print "{:<10} {:<5} {}".format(bact, bactdict[bact],  0)
        else:
            print "{:<10} {:<5} {}".format(bact, bactdict[bact],  100 * bactdict[bact] / totalhits)

    return filtered_dict




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze genome percentages in blat mapping results")
    parser.add_argument("mappings", metavar="BLAST8", 
            help="Mappings file in BLAST8 format")
    parser.add_argument("--id", metavar="I", dest="identity", default=100, type=int, 
            help="Minimum identity of mapped fragments [%(default)s]")
    parser.add_argument("--matches", metavar="M", dest="matches", default=20, type=int, 
            help="Minimum number of matches per fragment [%(default)s]")
    parser.add_argument("--mismatches", metavar="m", dest="mismatches", default=0, type=int, 
            help="Maximum number of mismatches per fragment [%(default)s]")
    parser.add_argument("-r", dest="remove", default=True, action="store_false", 
            help="Remove noninformative fragments [%(default)s]")
    parser.add_argument("-p", dest="print_fragments", default=False, action="store_true", 
            help="Print names of filtered fragments [%(default)s]")
    parser.add_argument("--print", dest="print_fragment", default="", 
            help="Print the hits of a single fragment identifier.")
    args = parser.parse_args()

    if len(argv) < 2:
        parser.print_help()
        exit()

    infodict = parse_mappings(args.mappings)
    filtered_dict = filter_hits(infodict, 
            remove_noninformative=args.remove,
            matches=args.matches,
            identity=args.identity,
            mismatches=args.mismatches,
            print_fragment=args.print_fragment,
            print_fragments=args.print_fragments)
