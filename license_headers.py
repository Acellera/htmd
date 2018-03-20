# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import linecache
import sys
import difflib
from itertools import tee

# Do not use the triple quote block, as the script can back-fire and remove part of this string
license_header = ('# (c) 2015-2018 Acellera Ltd http://www.acellera.com\n'
                  '# All Rights Reserved\n'
                  '# Distributed under HTMD Software License Agreement\n'
                  '# No redistribution in whole or part\n'
                  '#\n'
                  )

# THIS IS A VERY DANGEROUS SCRIPT. USE WITH CAUTION.
# BY DEFAULT, IT WRITES THE OUTPUTS TO TEST DIRECTORY AND PRINTS DELTAS
# IT IS ADVISED TO CHECK THOSE DELTAS BEFORE WRITING DIRECTLY ON THE REPO (./)
outdir = './test_license_headers'
# outdir = './'
exclusions = ['./test_license_headers', './htmdlib', './htmd/data']

for root, dirs, files in os.walk('.'):
    for i, fname in enumerate([os.path.join(root, file) for file in files if not any(exclusion in root for exclusion in exclusions)]):
        try:
            if fname.endswith(('.py', '.sh')):
                # find all consecutive comment lines that include a flag string (exclude shebang lines that come before)
                headerlines = []
                flagstring = "License Agreement"
                with open(fname, 'r') as f:
                    for j, line in enumerate(f, 1):
                        if line.strip().startswith('#') and not line.strip().startswith('#!'):
                            if j == 1:
                                linestart = 1
                            if not linecache.getline(fname, j-1).strip().startswith('#') or \
                                    linecache.getline(fname, j-1).strip().startswith('#!'):
                                linestart = j
                            if not linecache.getline(fname, j+1).strip().startswith('#'):
                                lineend = j
                                tmp = {}
                                for k in range(linestart, lineend+1):
                                    tmp[k] = linecache.getline(fname, k)
                                for e in tmp:
                                    if flagstring in tmp[e]:
                                        headerlines.append(tmp)
                                        break
                print('File {} has {} contiguous headers'.format(fname, len(headerlines)))
                for n, h in enumerate(headerlines):
                    print('Header {}: {}'.format(n+1, ''.join(h[e] for e in h)))
                    print('Lines of the header {}'.format([e for e in h]))

                outfname = os.path.abspath(os.path.join(outdir, fname))
                if not os.path.exists(os.path.dirname(outfname)):
                    os.makedirs(os.path.dirname(outfname))
                # remove the lines
                with open(fname, 'r') as r:
                    file_lines = {}
                    for j, line in enumerate(r, 1):
                        file_lines[j] = line
                with open(outfname, 'w') as w:
                    for j in file_lines:
                        if j not in [e for h in headerlines for e in h]:
                            w.write(file_lines[j])
                # print license header after: triple quote block or shebang at the top. if starts with comment, insert
                # empty line. if file is empty, there will be no header.
                with open(outfname, 'r') as r:
                    file_lines = {}
                    triple = 0
                    shebang = False
                    empty = False
                    for j, line in enumerate(r, 1):
                        if j == 1 and line.strip().startswith('#!'):
                            shebang = True
                        if j == 1 and line.strip().startswith('#') and not line.strip().startswith('#!'):
                            empty = True
                        if triple == -1 and ('\"\"\"' in line or '\'\'\'' in line):
                            triple = j+1
                        if j == 1 and (line.strip().startswith('\"\"\"') or line.strip().startswith('\'\'\'')):
                            triple = -1
                            if len(line.strip()) != 3 and (line.strip().endswith('\"\"\"') or line.strip().endswith('\'\'\'')):
                                triple = 2
                        file_lines[j] = line
                with open(outfname, 'w') as w:
                    for j in file_lines:
                        if (shebang is False and triple == 0 and j == 1) or (shebang is True and j == 2) \
                                or (triple and j == triple):
                            w.write(license_header)
                        if j == 1 and empty is True:
                            w.write('\n')
                        w.write(file_lines[j])
                # make diff
                if outdir != './':
                    with open(fname, 'r') as orig:
                        with open(outfname, 'r') as modf:
                            diff = difflib.unified_diff(
                                orig.readlines(),
                                modf.readlines(),
                                fromfile=fname,
                                tofile=outfname,
                            )
                            diff1, diff2 = tee(diff)
                            if len([i for i in diff1]) == 0:
                                print('No changes')
                            else:
                                print('Delta:')
                                for line in diff2:
                                    sys.stdout.write(line)
                print(''.join('=' for _ in range(70)))
        except:
            print("Failed to process {}".format(fname))
