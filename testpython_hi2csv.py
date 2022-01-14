#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Python script to read data from a HIVAP output file of type  
# hivaperg.dat, extract cross-section data etc. and write 
# the data in a "compact format" to a file or to stdout.
#
# The script is compatible with output produced by the following two
# version of HIVAP:
# - version from 1990-06-20 with English text
# - version from 1994-06-07 with German text
#
# Note: the HIVAP run must be performed with DISC=3 in the HIVAP
# input file.
#
# Johan Nyberg, 2015-2018
#
# ToDo list:
# ----------
# - Use functions for
#   1) parsing the input file to determine projectile and target used
#   2) parsing shell effects table
#   3) parsing fusion cross section tables
#   4) parsing residue cross section tables
#   5) printing results in file
# - Write only specific residual nuclei, option -n "102sn 103sn 102In".
# - Extract and print spin-distributions of final nuclei, if available.
# - Is there info in the input file regarding particle multiplicities
#   for each exit channel? If so, one can get info from that regarding
#    the emission of alpha versus 2p2n.

# For use of the python 3.x print function in python >=2.6 :
from __future__ import print_function

# Script version:
__version__ = "2.3.7"

author_string = "Johan Nyberg, 2015-2018"

# Import modules

import argparse
#import io
#import os
import os.path
import re # Regular expressions
import sys
import string

# Script and python version strong (for option --version):
script_version = sys.argv[0] + " " + __version__
python_version = "python " + ".".join(map(str, sys.version_info[:3]))
myversion = script_version + ", " + python_version

def check_python_version():
  # Check that we have python version >= 2.7, else exit
  #print(sys.version_info[0] + sys.version_info[1])
  if sys.version_info[0] > 2: return
  if sys.version_info[0] >= 2 and sys.version_info[1] >= 7: return
  #print("Your python version is too old!")
  print("This script requires python version 2.7 or higher")
  print("Your version is: ", sys.version)
  sys.exit(1)

def get_args(argv=None):
  # This function decodes the command line arguments and returns the
  # arguments 
  #
  # Create an ArgumentParser object caller parser:
  parser = argparse.ArgumentParser(
    description = "Read HIVAP output file of type hivaperg.dat, \
    extract cross section data etc. from the file and write it \
    to an output file.", epilog = author_string)

  # Use the add_argument method to take the strings on the 
  # command line and turn them into objects:
  if (sys.version_info > (3, 4)):
    # Python version >= 3.4 code in this block.
    # The hivaperg.dat files produced by with the HIVAP version from
    # 1994-06-07 (with german text) contains a character on the 
    # third line (a german 'double s' character) that is not compatible
    # with utf-8 encoding -> python3 gives the error:
    #   UnicodeDecodeError: 'utf-8' codec can't decode byte 0xdf in 
    #   position 115: invalid continuation byte
    # This can be fixed by using
    # ... type=argparse.FileType('r', encoding='ISO-8859-1') ...
    # or
    # ... type=argparse.FileType('r', errors='replace') ...
    # in parser.add_argument(). However, the keyword arguments 
    # encoding and errors are not supported by python versions <3.4. 
    # Hence, the following ugly test of python version:
    parser.add_argument("infile", 
                        type=argparse.FileType(mode='r', 
                                               encoding='ISO-8859-1'),
                        #type=argparse.FileType(mode='r', 
                        # errors='replace'),
                        #type=argparse.FileType(mode='r'),
                        default=sys.stdin,
                        help="input file, use '-' for stdin")
  else:
    parser.add_argument("infile", 
                        type=argparse.FileType(mode='r'),
                        default=sys.stdin,
                        help="input file, use '-' for stdin")

  # The output file is opened with mode='w' (mode 'x' does not 
  # exist in python version < 3.3!). The best would be implement
  # a version that is compatible with python >=2.7 and that
  # prints an error and exits if the output file exists and
  # has an option -o/--overwrite that will allow for overwriting
  # an existing output file. I could not figure out how to do that
  # with argparse. 
  parser.add_argument("outfile", 
#                      type=argparse.FileType(mode='x'),
                      type=argparse.FileType(mode='w'),
                      default=sys.stdout,
                      help="output file, use '-' for stdout. \
                      Note: the file is silently overwritten  \
                      if it exists!")
  parser.add_argument("-r", "--reverse", action="store_true",
                      help="print table in reverse energy order")
  parser.add_argument("-v", "--verbose", action="store_true",
                      help="verbose output (for debugging)")
  parser.add_argument("--version", 
                      action='version', version=myversion,
                      help="print script and python version")

  # Require at least one argument:
  if len(sys.argv) < 2:
      parser.print_usage()
      sys.exit(1)

  # Return the decode command line arguments:
  return parser.parse_args(argv)

def main():

  # Check that we have good enough python version:
  check_python_version()

  # Get the command line arguments
  argvals = None             # init argv in case not testing
  #argvals = '6 2 -v'.split() #example of passing test params to parser
  args = get_args(argvals)

  #  sys.exit(0)

  # Verbose:
  if args.verbose:
    print("Verbosity turned on")
    print('Python version: ' + sys.version)

  # Read input data:
  if args.verbose: print("Start reading data.")
  #if args.verbose: print("args.infile = ", args.infile)
  data = args.infile.read()
  if args.verbose: print("End reading data.")

  # Test if output file exists. If so exit with an error:
  # This does not work, because outfile is opened in get_args
  # (new file created if it does not exist!)
  #if os.path.exists(args.outfile.name) == True:
  #  print("Output file", args.outfile.name, 
  #    "exists. Exiting script ...")
  #  sys.exit(1)

  # Create list my_list, which contains the lines in the input data, 
  # breaking at line boundaries:
  my_list = data.splitlines()
  set_list = my_list # set_list will be used to extract data from the
                     # shell effects table
  
  # Return length of the list:
  nrl = len(my_list)
  if args.verbose:
    print("{}{}".format("Nr of lines in input data = ", nrl))

  # Create and initialise variables and lists:

  # Element symbols:
  zmax = 118
  ElSym = ["" for i in range(0, zmax+1)]
  ElSym = ["n", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
           "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
           "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
           "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
           "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
           "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr",
           "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
           "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au",
           "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
           "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
           "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
           "Ds", "Rg", "Cn", "Nh", "Fl", "Ms", "Lv", "Ts", "Og"]

  # Element symbols as dictionary for easy lookup of which Z value
  # a given symbol has, e.g. dictElSym['Sn'] is 50:
  dictElSym = {}
  for i in range(0, zmax+1):
    dictElSym[ElSym[i]] = i

  #if args.verbose:
  #  for i in range(0, zmax):
  #    print(i, ElSym[i], dictElSym[ElSym[i]])

  # Define some max values of list indexes:
  maxnuc = 512 # max nr of nuclides
  maxe = 512   # max nr E* and Elab values
  maxz = 256   # max nr of Z values
  maxn = 256   # max nr of N values
  
  # Lists of nuclide names, Z, N, A and emitted particles:
  NucNam = [0 for i in range(maxnuc)] # example: '102Sn'
  Z = [0 for i in range(maxnuc)]      # example: 50
  N = [0 for i in range(maxnuc)]      # example: 52
  A = [0 for i in range(maxnuc)]      # example: 102
  EP = [0 for i in range(maxnuc)]     # example: '1a0p2n'
  # In NucNr[][] the nuclide nr (from 0 to maxnuc-1) is stored
  NucNr =  [[0.0 for z in range(maxz)] for n in range(maxn)]
  # In Itotope[] the nuclide nr (from 0 to maxnuc-1) is stored
  Isotope = [0 for i in range(maxnuc)]

  # Beam energy in lab system:
  Elab = [0.0 for i in range(maxe)]

  # Excitation energy of compound nucleus (CN):
  Estar = [0.0 for i in range(maxe)] 

  # Fusion cross section etc.
  sigfus = [0.0 for i in range(maxe)]
  sigsum = [0.0 for i in range(maxe)]
  sigrat = [0.0 for i in range(maxe)]
  sigvr = [0.0 for i in range(maxe)]
  sigfis = [0.0 for i in range(maxe)]
  lcrit = [0.0 for i in range(maxe)]
  bfis = [0.0 for i in range(maxe)]

  # Residual cross section as 2D array:
  sigres = [[0.0 for i in range(maxnuc)] for e in range(maxe)]
  
  # nrnuc = total nr of final nuclides
  # nri = nr of isotopes per element. Same for all elements.

  # Defince variables for the two different HIVAP versions:
  
  hivap1990 = \
  'THIS IS HIVAP VERSION HIVA 20 JUN 90'
  hivap1994 = \
  'Programm HIVAP,  Autor  W.Reisdorf   Version fuer "rzri6f" vom  7. 6. 1994'

  #if args.verbose:
  #  print(hivap1990
  #  print(hivap1994)
  #  print(my_list[0])

  if hivap1990 in my_list[0]: 
    ver=1990
    fusion_table = "PERCNT"
    #"ELAB   MEV/U   EXCIT      SIGFUS SIGEVA SIGFIS LCRIT PERCNT FIS"
    shell_effect_table = "EXP SHELLS"
    # EXP SHELLS  SHELL0=     .00
    barf = "BARF"
    #barfac = "BARFAC"
    # Output may contains a lot of such lines:
    # **** BARFIT CALLED WITH Z LESS THAN 19 OR GREATER THAN 111; 
    # or only this line
    # BARFAC=   1.000      BAR0=    .000
    start_cross_section_table = "CROSS SECTIONS"
    #cross_section_table = "CROSS SECTIONS"
    # CROSS SECTIONS (MB)
    end_cross_section_table_1 = "Z-DISTRIBUTION"
    end_cross_section_table_2 = "ENERGY"
    end_cross_section_table_3 = "THIS WAS HIVAP"
    #z_distribution_table = "Z-DISTRIBUTION"
    estar_cross_section_table = "EXCIT"
    # EXCIT     99TE     100TE     101TE     102TE     103TE     104TE
  elif hivap1994 in my_list[0]: 
    ver=1994
    fusion_table = "SI_fus"
    #"E_lab  A*MeV   E*/MeV  SI_fus  SI_VR  SI_Spalt L_krit b_Spalt"
    shell_effect_table = "Schaleneffekte"
    # Exp. Schaleneffekte,  SHELL0=   -1.20
    barf = "BARF"
    #barfac = "BARFAC"
    # Output may contains a lot of such lines:
    # **** BARFIT CALLED WITH Z LESS THAN 19 OR GREATER THAN 111; 
    # or only this line
    # BARFAC=   1.000      BAR0=    .000
    start_cross_section_table = "Querschnitte"
    #cross_section_table = "Querschnitte"
    # Querschnitte / mbarn
    end_cross_section_table_1 = "Z-Verteilung"
    end_cross_section_table_2 = "ENERGY"
    end_cross_section_table_3 = "Ende HIVAP"
    #z_distribution_table = "Z-Verteilung"
    # E*/MeV   103Te     104Te     105Te     106Te     107Te     108Te
    estar_cross_section_table = "E*/MeV"
  else:
    print("Unknown version of HIVAP input file. First line should \
    contain version info. See details in script. Exiting script ...")
    sys.exit(1)

  if args.verbose: print("Input file produced by HIVAP version from", ver)

  # Start parsing for projectile and target --------------------------

  # Extract ZT, AT, ZP and AP from a line of this type
  #   ZT= 13    AT  27    ZP 36    AP  78
  # and calculate the NucNam, Z, N and A of the compound nuclues.
  
  reaction_found = 0
  for lv in my_list:  # lv = list value
    if len(lv) == 0: # Skip all empty lines
      continue
    #if args.verbose: print(lv)
    if reaction_found == 1:
      break
    if "ZT=" in lv: # Line containing ZT, AT, ZP, AP
      #if args.verbose: print(lv)
      reaction_found = 1
      reaction = lv
      words = lv.split()
      nw = len (words)
      zt = int(words[1]); at = int(words[3]); nt = at - zt
      zp = int(words[5]); ap = int(words[7]); np = ap - zp
      tname = str(at) + ElSym[zt]
      pname = str(ap) + ElSym[zp]
  # End of lv for loop 

  # Find NucNam, Z, N and A of compound nucleus, obtained as
  # the nuclide with the largest Z and N:
  zcn = zp + zt; ncn = np + nt; acn = ap + at
  cnname = str(acn) + ElSym[zcn]
  #if args.verbose: print(cnname, zcn, ncn, acn)

  # Endparsing for projectile and target --------------------------
  
  # Start parsing shell effects table --------------------------

  # Extract the nuclide names of all nuclides produced. The names
  # are extracted from the table with the shell effects, which 
  # starts with the line containing the value in string 
  # shell_effect_table and ends with a line that contains the value in
  # string barf. The table itself contains two lines for each 
  # element. The first line contains the isotopes of the element, e.g.
  # Sn 105  Sn 106  Sn 107  Sn 108  Sn 109  Sn 110  Sn 111  Sn 112
  #
  # The second line contains the shell effect values in MeV for each 
  # isotope, e.g.
  # -8.00   -7.00   -6.20   -5.40   -4.81   -4.21   -3.79   -3.36

  # The variable barf_found is set = 1 when the data from the shell
  # effects table have been extracted, which means that a break of
  # the for loop below can be made.

  barf_found = 0
  nrnuc = 0 # counter of number of nuclides
  
  for l in range(0, nrl): # Loop through all lines.
    if barf_found == 1: break # We have parsed through the shell effects table
    #if args.verbose: print(l, my_list[l])
    #if "Exp. Schaleneffekte" in my_list[l]:
    if shell_effect_table in my_list[l]: # Found shell effect table
      # Copy the data in this table to the list set_list.
      # Remove all empty lines while doing this.
      # The table ends by a line containing BARFAC or BARFIT.
      nrlset = 0 # nr of lines in the extracted shell effect table
      for l2 in range(l+1, nrl):
        if len(my_list[l2]) == 0: # Skip all empty lines
          continue
        if(barf in my_list[l2]): 
          if args.verbose: print("Found BARF -> end of shell effect table")
          #if(my_list[l2] == "" or barf in my_list[l2]):
          #if args.verbose: print("Found empty line or BARF")
          barf_found = 1
          break
        set_list[nrlset] =  my_list[l2]
        nrlset = nrlset + 1
      # Now extract the data from the shell effect table list set_list:
      for l2 in range(0, nrlset, 2):
        if args.verbose: print(l2, set_list[l2])
        words = set_list[l2].split()
        nw = len (words)
        if args.verbose: print("nw = ", nw, "words ==", words)
        # nri = nr of isotopes per element. Same for all elements.
        # Used later when extracting cross section data/
        nri = int(nw/2)
        if args.verbose: print("nri = ", nri)
        for w in range(0, nw, 2):
          if args.verbose:
            print("w = ", w, ", ", end="")
            print("words[", w, "] = ", words[w], ", ", end="")
            print("words[", w+1, "] = ", words[w+1])
          s1 = words[w+1] + words[w] # Nuclide string e.g. '100Sn'
          # Extract the mass nr A, element symbol and atomic nr Z: 
          a1 = re.findall('\d+', s1) # A as list, e.g. ['100']
          a = int(a1[0])             # A values as int, e.g. 100
          es1 = re.findall('\D+', s1) # element symbol as list e.g. ['Sn']
          es = es1[0]                 # element symbol as string e.g. 'Sn'
          if es == "J": es = "I" # J is used instead of I for Iodine
          #print("s1 = ", s1, ", a = ", a, "es = ", es)
          z = dictElSym[es]           # Z, e.g. 50
          #if args.verbose:
          #  print("s1 = ", s1, ", a = ", a, "es = ", es, "z = ", z)
          #print("z = ", z)
          NucNam[nrnuc] = s1
          A[nrnuc] = int(a)
          Z[nrnuc] = int(z)
          N[nrnuc] = a-z
          NucNr[Z[nrnuc]][N[nrnuc]] = nrnuc
          nrnuc = nrnuc + 1
        # End of wfor loop through words on line
      # End of l2 for loop of shell effects table
    # End of shell effect table  

    # Create the list of emitted particles EP, e.g. 2p2n, 1a1p0n, etc.
    for i in range(0, nrnuc):
      nrn = ncn - N[i]
      nrp = zcn - Z[i] 
      if nrn == 2 and nrp == 2:
        nra = 0
      else:
        nra = int(min(nrn/2, nrp/2))
        nrn = int(nrn) - int(nra*2)
        nrp = int(nrp) - int(nra*2)
      EP[i] = str(nra) + "a" + str(nrp) + "p" + str(nrn) + "n"

  # End of l loop through all lines of my_list
  
  if (nrnuc == 0):
    print("Something wrong with the input file: could not find table ",
          "with shell effects (Exp. Schaleneffekte). ",
          "Exiting script...")
    sys.exit(1)

  if args.verbose: 
    print("Total nr of final nuclides =", nrnuc)
    print("Nr of isotopes per element =", nri)
    print("NucNam, Z, N, A, nra, nrp, nrn, EP:")
    for i in range(0, nrnuc):
      print(NucNam[i], Z[i], N[i], A[i], nra, nrp, nrn, EP[i])

  # End parsing shell effects table --------------------------

  # Start parsing fusion table --------------------------

  # Extract the fusion cross section etc.
  #
  # The beam energy, fusion cross section, l_crit etc are in 
  # tables containing the value of the string fusion_table 
  # in the header line. The tables end when a line containing the
  # value of string start_cross_section_table is found.
  #
  # Note that if there are more than 16 lab energy runs in the
  # HIVAP output file, there will be more than one fusion table!

  #dictEstar = {} # Dictionary of Estar for lookup of index to an Estar value
  # The above dict cannot be used because sometimes the Estar values are
  # the same in the HIVAP output file even though the ELAB values are different
  # (only one decimal is used in the file).
  
  ie = 0 # counter of energy value
  in_fusion_table = 0 # set = 1 if we are inside a fusion table
  
  for lv in my_list:  # lv = list value
    #if args.verbose: print(lv)
    if len(lv) == 0: # Skip all empty lines
      continue
    if fusion_table in lv: # Start of a fusion table
      in_fusion_table = 1
      continue
    if start_cross_section_table in lv: # End of fusion table
      in_fusion_table = 0
      continue
    if in_fusion_table == 1: # Start extracting data from a fusion table
      #if args.verbose:
      #  print(lv)
      words = lv.split()
      Elab[ie] = float(words[0])
      Estar[ie] = float(words[2])
      sigfus[ie] = float(words[3])
      sigvr[ie] = float(words[4])
      sigfis[ie] = float(words[5])
      lcrit[ie] = float(words[6])
      bfis[ie] = float(words[7]) #casout.f: bfis = sigfis/sigfus*100
      #dictEstar[Estar[ie]] = ie # This can be used to find the ie value for
      #                          # a given Estar value.
      #print(ie, Estar[ie])
      ie = ie + 1
      continue
    # End of if in_fusion_table
  # End of lv for loop
  nre = ie # Nr of energy values
    
  if args.verbose:
    print("Nr of energy values = ", nre)
    print("Elab = ", end = "")
    for i in range(nre):
      print(Elab[i], " ", end = "")
    print("")
    for i in range(nre):
      print(Estar[i], " ", end = "")
    print("")
    print("sigfus = ", end = "")
    for i in range(nre):
      print(sigfus[i], " ", end = "")
    print("")
    #print("dictEstar = ", dictEstar)

  # End parsing fusion tables --------------------------

  # Start parsing cross section tables --------------------------

  # The cross section data for the residual nuclei are located
  # between the lines containing the value of the strings 
  # start_cross_section_table and (end_cross_section_table_1 or
  # end_cross_section_table_2). The header of these tables look
  # like this:
  #   HIVAP version 1990:
  #   EXCIT    100In     101In     102In     103In     104In     105In
  #   HIVAP version 1994:
  #   E*/MeV   100In     101In     102In     103In     104In     105In
  # The variable in_xsects_table = 1 when these tables are processed.
    #
  # Note that if there are more than 16 lab energy runs in the
  # HIVAP output file, there will be more than one cross section table!
  
  iso = 0     # counter of isotope nr in each header line
  iso2 = 0    # index to first isotope in xsect table
  nt = 0      # counter of cross section table nr
  in_xsect_table = 0 # set = 1 if we are inside a cross section table

  ie1 = 0 # subcounter of energy value (reset to 0 for each new
          # res. nucl. cross section table
  ie = 0 # counter of energy values
  
  for lv in my_list:  # lv = list value
    #if args.verbose: print(lv)
    if len(lv) == 0: # Skip all empty lines
      continue
    if start_cross_section_table in lv: # Start of a cross section table
      in_xsect_table = 1
      ie = ie + ie1 # increment energy counter index
      continue
    if end_cross_section_table_1 in lv or \
       end_cross_section_table_2 in lv or \
       end_cross_section_table_3 in lv:       # End of cross section table
      iso = 0  # reset isotope counter
      iso2 = 0 # reset index to first isotope in xsect table
      nt = 0   # reset xsect table counter
      in_xsect_table = 0 # not parsing xsect tables anymore
      continue
    if in_xsect_table == 1: # Start extracting data from cross section tables
      #if args.verbose: print(lv)
      if estar_cross_section_table in lv: # header line of xsect table
        ie1 = 0
        #print("ie1 = ", ie1, ", ie = ", ie)
        words = lv.split()
        nrw = len(words)
        estarstr = words[0] # Check that first word = estar_cross_section_table
        if estarstr != estar_cross_section_table:
          print("Something is wrong with the first word of the header ",
                "line of the cross section table!")
          print("The header line is:")
          print(lv)
          print("The first word should be = ", estar_cross_section_table)
          print("Check consistency of input data file and script!")
          print("Exiting script...")
          sys.exit(1)
        # nri = nr of isotopes per line obtainefrom shell effects table
        if(nrw - 1 == nri):
          nrwstp = 1  # no space between mass number and element symbol
        else:
          nrwstp = 2  # no space between mass number and element symbol
        #if args.verbose: print("nrwstp = ", nrwstp)
        # Extract the index to the isotopes in the current xsect table.
        # These are saved in Isotope[iso] for use below when parsing the
        # data in the tables.
        for iw in range(1, nrw, nrwstp): # loop from word nr 1 to the end
          #if args.verbose: 
          #  print("iw = ", iw, "words[iw] = ", words[iw])
          if(nrwstp == 1):
            s1 = words[iw] # This word should be the nuclide name e.g. 100Sn
          else:
            s1 = words[iw] + words[iw+1] # join two words, e.g. 100 and Sn.
          # Extract the mass nr, element symbol and atomic nr: 
          a1 = re.findall('\d+', s1) # A as list, e.g. ['100']
          a = int(a1[0])             # A values as int, e.g. 100
          es1 = re.findall('\D+', s1) # element symbol as list e.g. ['Sn']
          es = es1[0]                 # element symbol as string e.g. 'Sn'
          if es == "J": es = "I" # J is used instead of I for Iodine
          #if args.verbose:
          #  print("s1 = ", s1, ", a = ", a, "es = ", es)
          z = dictElSym[es]           # Z, e.g. 50
          n = a - z                   # N
          Isotope[iso] = NucNr[z][n] # save the nuclide nr in Isotope[].
          #if args.verbose:
          #  print("s1 = ", s1, ", a = ", a, "z = ", z, "n = ", n,
          #        "iso = ", iso, "Isotope[iso] = ", Isotope[iso])
          iso = iso + 1
        # End iw for loop
        iso2 = nt * nri # index to first isotope in the next xsect table
        nt = nt + 1 # xsect table nr 
      else:
        #print(lv)
        # data lines (one per energy) of xsect tables:
        # E* (col 1) and xsects (cols 2-) 
        words = lv.split()
        nrw = len(words)
        for iw in range(0, nrw): # Extract E* and cross sections
          #if args.verbose: 
          #  print("iw = ", iw, "words[iw] = ", words[iw])
          if iw == 0: # This is the first column with the E* value
            estar = float(words[0])
            #ie = dictEstar[estar]
            #if args.verbose:
            #  print("estar = ", estar, "Estar[", ie1 + ie, "] = ", \
            #        Estar[ie1 + ie])
            if estar != Estar[ie1 + ie]:
              print("Error: E* values are out of synch for energy nr ie1+ie = ",
                    ie1 + ie)
              print("E* value from fusion cross section table Estar[ie1+ie] = ",
                    Estar[ie1 + ie])
              print("E* value from res. nuc. cross section table estar = ",
                    estar)
              print("Check consistency of input data file and script!")
              print("Line:")
              print(lv)
              print("Exiting script...")
              sys.exit(1)
          else: # cross section values in the remaining nri columns 
            # iso2 points to the first isotope in the current xsect table,
            # while iso2 + iw - 1 points to each isotope in this table.
            inuc = Isotope[iso2 + iw -1] 
            #if args.verbose:
            #  print("ie1 = ", ie1, ", ie = ", ie, "inuc = ", inuc)
            sigres[inuc][ie1 + ie] = float(words[iw])
            #if args.verbose:
            #  print("iw = ", iw, ", iso = ", iso, "iso2 = ", iso2,
            #         ", iso2 + iw - 1 = ", iso2 + iw - 1,
            #         ", inuc =", inuc, ", ie1 =", ie1, ie =", ie,
            #         ", Estar[ie1+ie] = ", Estar[ie1+ie],
            #         ", NucNam[inuc] =", NucNam[inuc],
            #         ", sigres[inuc][ie1+ie] =", sigres[inuc][ie1+ie])
          # End if/else w=0 or >0 
        # End iw for loop
        ie1 = ie1 + 1 # increment the energy counter
      # End if/else header line or data lines
      #if args.verbose: print("enr=", enr)
    # End if in_xsect_table = 1
  # End lv for loop
  nre2 = ie + ie1 # Nr of energy values

  # Check that nr of energy values from parsing the fusion table is
  # the same as for parsing the residue cross section tables. If not,
  # something is wrong!
  if (nre != nre2):
    print("Nr of energy values in the fusion table = ", nre, \
          "is != nr of energy values in the residual nuclei ", \
          "cross section table = ", nre2)
    print("Check consistency of input data file and script!")
    #print("Line:")
    #print(lv)
    print("Exiting script...")
    sys.exit(1)

  # End of parsing cross section tables --------------------------

  # Calculate SIGSUM and SIGSUM/SIGFUS ----------------------------

  for e in range(0, nre):
    sigsum[e] = 0.0
    for n in range(0, nrnuc):
      sigsum[e] = sigsum[e] + sigres[n][e]
      if (sigfus[e] > 0):
        sigrat[e] = sigsum[e] / sigfus[e]
      else:
        sigrat[e] = "undef"
    #if args.verbose: print(sigsum[e], sigrat[e])

  # Print results -------------------------------------------------
  
  # Print lines starting with hash:
  
  reaction = pname + " + " + tname + " -> " + cnname
  print("# HIVAP version: ", ver, file = args.outfile)
  print("# Reaction (beam + target-> CN): ", reaction, file = args.outfile)
  print("# Energies in MeV, cross sections in mb, ",
         "angular momenta in hbar", file = args.outfile)
  print("# SIGSUM = sum of residual nucleus cross sections. ",
         "SIGRAT = SIGSUM/SIGFUS", file = args.outfile)
  print("# ", file = args.outfile)

  # Print row with column nrs:
  print("{:>6}".format("#    0"), end='', file = args.outfile)
  print("{:>6}".format("     1"), end='', file = args.outfile)
  print("{:>10}".format("      2"), end='', file = args.outfile)
  print("{:>10}".format("      3"), end='', file = args.outfile)
  print("{:>10}".format("      4"), end='', file = args.outfile)
  print("{:>10}".format("      5"), end='', file = args.outfile)
  print("{:>10}".format("      6"), end='', file = args.outfile)
  print("{:>7}".format("     7"), end='', file = args.outfile)
  for n in range(0, nrnuc):
    cstr = str(n + 8)
    print("{:>13}".format(cstr), end='', file = args.outfile)
  print("", file = args.outfile)

  # Print table header lines:
  
  print("{:>6}".format("  ELAB"), end='', file = args.outfile)
  #print("{:>6}".format("ELAB/A"), end='', file = args.outfile)
  print("{:>6}".format("E*_CN"), end='', file = args.outfile)
  print("{:>10}".format(" SIGFUS"), end='', file = args.outfile)
  print("{:>10}".format(" SIGSUM"), end='', file = args.outfile)
  print("{:>10}".format(" SIGRAT"), end='', file = args.outfile)
  print("{:>10}".format("SIGEVAP"), end='', file = args.outfile)
  print("{:>10}".format(" SIGFIS"), end='', file = args.outfile)
  print("{:>7}".format(" LCRIT"), end='', file = args.outfile)
  for n in range(0, nrnuc):
    s1 =  NucNam[n]
    nucstr = " ".join(s1.split())
    nucstr =  EP[n] + "+" + " ".join(s1.split())
    print("{:>13}".format(nucstr), end='', file = args.outfile)
  print("", file = args.outfile)
  # This will print # followed by 69 spaces:
  #h = "#"
  #print("{:<69s}".format(h), end = "", file = args.outfile)
  #for n in range(0, nrnuc):
  #  s1 = EP[n]
  #  epstr = "".join(s1.split())
  #  print("{:>10}".format(epstr), end='', file = args.outfile)
  #  print("", file = args.outfile);

  print("{:>6}".format(" [MeV]"), end='', file = args.outfile)
  #print("{:>6}".format("MeV/A"), end='', file = args.outfile)
  print("{:>6}".format(" [MeV]"), end='', file = args.outfile)
  print("{:>10}".format("   [mb]"), end='', file = args.outfile)
  print("{:>10}".format("   [mb]"), end='', file = args.outfile)
  print("{:>10}".format("   [mb]"), end='', file = args.outfile)
  print("{:>10}".format("   [mb]"), end='', file = args.outfile)
  print("{:>10}".format("   [mb]"), end='', file = args.outfile)
  print("{:>7}".format("  hbar"), end='', file = args.outfile)
  for n in range(0, nrnuc):
    print("{:>13}".format("      [mb]"), end='', file = args.outfile)
  print("", file = args.outfile)
  
  #sys.exit(1)
  
  # Set the loop variables for printing with reverse order if
  # asked for:
  sta = 0;
  sto = nre
  stp = 1;
  if args.reverse:
    sta = nre - 1
    sto = -1
    stp = -1

  # Print data of the table:
  for e in range(sta, sto, stp):
    print("{:6.1f}".format(Elab[e]), end='', file = args.outfile)
    print("{:6.1f}".format(Estar[e]), end='', file = args.outfile)
    print("{:10.3e}".format(sigfus[e]), end='', file = args.outfile)
    print("{:10.3e}".format(sigsum[e]), end='', file = args.outfile)
    print("{:10.3e}".format(sigrat[e]), end='', file = args.outfile)
    print("{:10.3e}".format(sigvr[e]), end='', file = args.outfile)
    print("{:10.3e}".format(sigfis[e]), end='', file = args.outfile)
    print("{:7.1f}".format(lcrit[e]), end='', file = args.outfile)
    #print("{:7.2f}".format(bfis[e]), end='', file = args.outfile)
    for n in range(0, nrnuc):
      print("{0:13.2e}".format(sigres[n][e]), end='', file = args.outfile)
    print("", file = args.outfile)

# Execute function main:
main()  
