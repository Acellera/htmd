"""
    Extentions for PDB2PQR Suite

    This module provides various utilities for the PDB2PQR suite to be
    imported into other Python scripts.
    
    ----------------------------
   
    PDB2PQR -- An automated pipeline for the setup, execution, and analysis of
    Poisson-Boltzmann electrostatics calculations

    Copyright (c) 2002-2011, Jens Erik Nielsen, University College Dublin; 
    Nathan A. Baker, Battelle Memorial Institute, Developed at the Pacific 
    Northwest National Laboratory, operated by Battelle Memorial Institute, 
    Pacific Northwest Division for the U.S. Department Energy.; 
    Paul Czodrowski & Gerhard Klebe, University of Marburg.

    All rights reserved.

    Redistribution and use in source and binary forms, with or without modification, 
    are permitted provided that the following conditions are met:

        * Redistributions of source code must retain the above copyright notice, 
          this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright notice, 
          this list of conditions and the following disclaimer in the documentation 
          and/or other materials provided with the distribution.
        * Neither the names of University College Dublin, Battelle Memorial Institute,
          Pacific Northwest National Laboratory, US Department of Energy, or University
          of Marburg nor the names of its contributors may be used to endorse or promote
          products derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
    IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
    INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
    BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
    OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED 
    OF THE POSSIBILITY OF SUCH DAMAGE.

    ----------------------------
"""

import pkgutil

from optparse import OptionGroup, OptionConflictError, Option

_extList = [name for _, name, _ in pkgutil.iter_modules(__path__)]

extDict = {}

for extName in _extList:
    extDict[extName] = __import__(extName,globals(),locals(),[], 1)
    
def setupExtensionsOptions(parser):
    """
    Takes an instance of an OptionParser 
    and adds the options for all extensions
    
    If an extension adds it's own options, those
    options are put in their own group.
    """
    
    if len(extDict) == 0:
        return None
    
    firstGroup = OptionGroup(parser,"Extension options")
    groups = [firstGroup]
    
    for extName, extModule in list(extDict.items()):
        helpArg = {}
        if hasattr(extModule, 'usage'):
            helpArg['help'] = extModule.usage()
            
        extOption = Option('--' + extName, action='append_const', 
                             const = extName, dest = 'active_extensions', **helpArg )
        
        try:
            if hasattr(extModule, 'addExtensionOptions'):
                group = OptionGroup(parser, extName.capitalize() + " extension options")
                group.add_option(extOption)
                
                extModule.addExtensionOptions(group)
                
                if len(group.option_list) == 1:
                    opt = group.option_list[0]
                    group.remove_option(opt.get_opt_string())
                    firstGroup.add_option(opt)
                else:
                    groups.append(group)   
                
            else:
                firstGroup.add_option(extOption)
                
        except OptionConflictError as value:
            print('Error adding command line options for extension ' + extName + ' ' + '(' + str(value) + ')')
            
    
    for group in groups:
        parser.add_option_group(group)
    
    return groups

class extOutputHelper(object):
    """
    Simple class that makes writing to both file and output simple.
    """
    def __init__(self, routines, outfile):
        self.routines = routines
        self.outfile = outfile

    def write(self, s):
        self.routines.write(s)
        self.outfile.write(s)
