#!/usr/bin/python

import sys
import math
import os
import re
qcdbpkg_path = os.path.expanduser('~loriab/linux/qcdb')
sys.path.append(qcdbpkg_path)
import qcdb
#sys.path.append(os.path.abspath('./../qcdb'))
#print sys.path
import qcdb.basislist
sys.path.append(qcdbpkg_path + '/databases')


# load docstring info from database files (doesn't actually import database modules)
DBdocstrings = qcdb.dictify_database_docstrings()

# instructions
print """
 Welcome to imake-db.
    Just fill in the variables when prompted. 
    Hit ENTER to accept default.
    Strings should not be in quotes.
    Elements in arrays should be space-delimited.
    Nothing is case sensitive.
"""

# query database name
module_choices = dict(zip([x.upper() for x in DBdocstrings.keys()], DBdocstrings.keys()))

print 'n Choose your database.'
for item in module_choices.keys():
    print '    %-12s   %s' % ('[' + module_choices[item] + ']', DBdocstrings[module_choices[item]]['general'][0].lstrip(" |"))
print '\n'

user_obedient = False
while not user_obedient:
    temp = raw_input('    dbse = ').strip()
    if temp.upper() in module_choices.keys():
        db_name = module_choices[temp.upper()]
        user_obedient = True

# query database subset
subset_choices = dict(zip([x.upper() for x in DBdocstrings[db_name]['subset'].keys()], DBdocstrings[db_name]['subset'].keys()))

print '\n Choose your subset (multiple allowed).'
for key, val in DBdocstrings[db_name]['subset'].items():    
    print '    %-12s   %s' % ('[' + key + ']', val)
print '\n'

subset = []
user_obedient = False
while not user_obedient:
    temp = raw_input('    subset [all] = ').strip()
    ltemp = temp.split()
    if temp == "":
        user_obedient = True
    for item in ltemp:
        if item.upper() in subset_choices.keys():
            subset.append(subset_choices[item.upper()])
            user_obedient = True
        else:
            user_obedient = False
            subset = []
            break

# query qc program
print """
 Choose your quantum chemistry program.
    [qchem]
    [molpro]
    [psi4]         writes Q-Chem input files sufficient
    [nwchem]
    [xyz]          writes basic xyz files only
"""

user_obedient = False
while not user_obedient:
    temp = raw_input('    qcprog = ').strip()
    if temp.lower() in ['qchem', 'molpro', 'psi4', 'nwchem', 'xyz']:
        qcprog = temp.lower()
        user_obedient = True



# query basis set(s)
print """
 Choose your basis set (multiple allowed).
    e.g., aug-cc-pvdz or 6-31+G* or cc-pvtz may-cc-pvtz aug-cc-pvtz
"""

bases = []
user_obedient = False
while not user_obedient:
    temp = raw_input('    bases = ').strip()
    ltemp = temp.split()
    for item in ltemp:
        btemp = qcdb.basislist.corresponding_orbital(item)
        if btemp:
            bases.append(btemp)
            user_obedient = True
        else:
            print '    Basis set %s not recognized.' % (item)
            proceed = qcdb.query_yes_no('    Proceed anyway? =', False)
            if proceed:
                bases.append(item)
                user_obedient = True
            else:
                bases = []
                user_obedient = False
                break

# query directory prefix 
print """
 State your destination directory prefix.
"""

user_obedient = False
while not user_obedient:
    temp = raw_input('    dirprefix [try] = ').strip()
    if temp == "":
        dirprefix = 'try'
        user_obedient = True
    if temp.isalnum():
        dirprefix = temp
        user_obedient = True

# query memory
print """
 Choose your memory usage in MB.
"""

user_obedient = False
while not user_obedient:
    temp = raw_input('    memory [1600] = ').strip()
    if temp == "":
        memory = 1600
        user_obedient = True
    if temp.isdigit():
        memory = int(temp)
        user_obedient = True

# Load module for requested database
try: 
    database = __import__(db_name) 
except ImportError: 
    print('\nPython module for database %s failed to load\n\n' % (db_name)) 
    print('\nSearch path that was tried:\n') 
    print(", ".join(map(str, sys.path)))
    raise ValidationError("Python module loading problem for database " + str(db_name))
else:
    dbse = database.dbse
    HRXN = database.HRXN
    ACTV = database.ACTV
    RXNM = database.RXNM
    BIND = database.BIND
    TAGL = database.TAGL
    GEOS = database.GEOS
    try:
        DATA = database.DATA
    except AttributeError:
        DATA = {}


print """
        <<< SCANNED SETTINGS  SCANNED SETTINGS  SCANNED SETTINGS  SCANNED SETTINGS >>>

                          dbse = %s
                        subset = %s
                        qcprog = %s
                   HFUNCTIONAL = 
                         bases = %s
                     dirprefix = %s
                   memory [MB] = %d
                       usesymm = 
                       cast up =

        <<< SCANNED SETTINGS  DISREGARD RESULTS IF INAPPROPRIATE  SCANNED SETTINGS >>>

""" % (dbse, subset, qcprog, bases, dirprefix, memory)


# establish multiplicity hash table
mult = {
   1: "singlet",
   2: "doublet",
   3: "triplet",
   4: "quartet",
   5: "quintet",
   6: "sextet",
   7: "septet",
   8: "octet"}


# file extension
fext = 'xyz' if qcprog == 'xyz' else 'in'

# manipulate memory TODO: move to qcprog file
#if qcprog == 'molpro':
#    memory = int(math.ceil(memory / 8.0))

# merge and condense HRXN from subset
if len(subset) == 0:
    pass
else:
    temp = []
    for item in subset:
        if item == 'small':
            temp.append(database.HRXN_SM)
        elif item == 'large':
            temp.append(database.HRXN_LG)
        elif item == 'equilibrium':
            temp.append(database.HRXN_EQ)
        else:
            try:
                temp.append(getattr(database, item))
            except AttributeError:
                try:
                    temp.append(getattr(database, 'HRXN_' + item))
                except AttributeError:
                    raise ValidationError('Special subset \'%s\' not available for database %s.' % (item, db_name))
    HRXN = qcdb.drop_duplicates(temp)
#print 'HRXN', HRXN

# TODO: choose ACTV or merge ACTV
temp = []
for rxn in HRXN:
    temp.append(ACTV['%s-%s' % (dbse, rxn)])
HSYS = qcdb.drop_duplicates(temp)
#print 'HSYS', HSYS


methods = ['B3LYP-D', 'DF-MP2']
# commence the file-writing loop
tdir = '-'.join([dirprefix, dbse, qcprog])
try:
    os.mkdir(tdir)
except OSError:
    print 'Warning: directory %s already present.' % (tdir)

for basis in bases:
    basdir = qcdb.basislist.sanitize_basisname(basis)

    # TODO: handling basis sets skipped

    for method in methods:
        subdir = '-'.join([basdir, method])
        try:
            os.mkdir(tdir + '/' + subdir)
        except OSError:
            print 'Warning: directory %s/%s already present.' % (tdir, subdir)
    
        # TODO: forcing c1 symm skipped - still needed for xdm and molpro

#        ans = open('ssi.big', 'w')

        for system in HSYS:
            sfile = tdir + '/' + subdir + '/' + system + '.' + fext
            infile = open(sfile, 'w')

            GEOS[system].tagline = TAGL[system]
            GEOS[system].update_geometry()

            # write start of file and comment line
            if qcprog == 'qchem':
                infile.write('$comment\n%s %s\n$end\n\n' % (system, TAGL[system]))
            elif qcprog == 'molpro':
                infile.write('***, %s %s\n' % (system, TAGL[system]))
                infile.write('memory,%d,m\n' % (int(math.ceil(memory / 8.0))))
            elif qcprog == 'psi4':
                infile.write('# %s %s\n\n' % (system, TAGL[system]))
            elif qcprog == 'nwchem':
                infile.write('echo\ntitle "%s %s"\n\n' % (system, TAGL[system]))
            elif qcprog == 'xyz':
                pass

            # write molecule section
            if qcprog == 'xyz':
                infile.write(GEOS[system].save_string_xyz())
            elif qcprog == 'psi4':
                infile.write(GEOS[system].format_molecule_for_psi4())

            infile.close()

#            if re.search('dimer', system):
#                ans.write('%d\n' % (GEOS[system].natom()))
#        ans.close()

