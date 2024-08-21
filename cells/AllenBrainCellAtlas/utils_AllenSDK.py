# Allen Institute Software License - This software license is the 2-clause BSD
# license plus a third clause that prohibits redistribution for commercial
# purposes without further permission.
#
# Copyright 2015-2017. Allen Institute. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# 3. Redistributions for commercial purposes are not permitted without the
# Allen Institute's written permission.
# For purposes of this license, commercial purposes is the incorporation of the
# Allen Institute's software into anything for which you will charge fees or
# other compensation. Contact terms@alleninstitute.org for commercial licensing
# opportunities.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
import logging
import os
from allensdk.model.biophys_sim.neuron.hoc_utils import HocUtils
from allensdk.core.nwb_data_set import NwbDataSet
import numpy as np
from pkg_resources import resource_filename  # @UnresolvedImport

PERISOMATIC_TYPE = "Biophysical - perisomatic"
ALL_ACTIVE_TYPE = "Biophysical - all active"

class Utils(HocUtils):
    '''A helper class for NEURON functionality needed for
    biophysical simulations.

    Attributes
    ----------
    h : object
        The NEURON hoc object.
    nrn : object
        The NEURON python object.
    neuron : module
        The NEURON module.
    '''

    _log = logging.getLogger(__name__)

    def __init__(self, description):
        self.update_default_cell_hoc(description)

        super(Utils, self).__init__(description)

    def update_default_cell_hoc(self, description, default_cell_hoc='cell.hoc'):
        ''' replace the default 'cell.hoc' path in the manifest with 'cell.hoc' packaged
        within AllenSDK if it does not exist '''

        hoc_files = description.data['neuron'][0]['hoc']
        try:
            hfi = hoc_files.index(default_cell_hoc)

            if not os.path.exists(default_cell_hoc):
                abspath_ch = resource_filename(__name__,
                                               default_cell_hoc)
                hoc_files[hfi] = abspath_ch

                if not os.path.exists(abspath_ch):
                    raise IOError("cell.hoc does not exist!")

                self._log.warning("Using cell.hoc from the following location: %s", abspath_ch)
        except ValueError as e:
            pass


    def generate_morphology(self, morph_filename):
        '''Load a swc-format cell morphology file.

        Parameters
        ----------
        morph_filename : string
            Path to swc.
        '''
        h = self.h

        swc = self.h.Import3d_SWC_read()
        swc.input(morph_filename)
        imprt = self.h.Import3d_GUI(swc, 0)

        h("objref this")
        imprt.instantiate(h.this)

        h("soma[0] area(0.5)")
        for sec in h.allsec():
            sec.nseg = 1 + 2 * int(sec.L / 40.0)
            if sec.name()[:4] == "axon":
                h.delete_section(sec=sec)
        h('create axon[2]')
        for sec in h.axon:
            sec.L = 30
            sec.diam = 1
            sec.nseg = 1 + 2 * int(sec.L / 40.0)
        h.axon[0].connect(h.soma[0], 0.5, 0.0)
        h.axon[1].connect(h.axon[0], 1.0, 0.0)

        h.define_shape()


    def load_cell_parameters(self):
        '''Configure a neuron after the cell morphology has been loaded.'''
        passive = self.description.data['passive'][0]
        genome = self.description.data['genome']
        conditions = self.description.data['conditions'][0]
        h = self.h

        h("access soma")

        # Set fixed passive properties
        for sec in h.allsec():
            sec.Ra = passive['ra']
            sec.insert('pas')
            for seg in sec:
                seg.pas.e = passive["e_pas"]

        for c in passive["cm"]:
            h('forsec "' + c["section"] + '" { cm = %g }' % c["cm"])

        # Insert channels and set parameters
        for p in genome:
            if p["section"] == "glob":  # global parameter
                h(p["name"] + " = %g " % p["value"])
            else:
                if p["mechanism"] != "":
                    h('forsec "' + p["section"] +
                      '" { insert ' + p["mechanism"] + ' }')
                h('forsec "' + p["section"] +
                  '" { ' + p["name"] + ' = %g }' % p["value"])

        # Set reversal potentials
        for erev in conditions['erev']:
            h('forsec "' + erev["section"] + '" { ek = %g }' % erev["ek"])
            h('forsec "' + erev["section"] + '" { ena = %g }' % erev["ena"])


class AllActiveUtils(Utils):

    def generate_morphology(self, morph_filename):
        '''Load a neurolucida or swc-format cell morphology file.

        Parameters
        ----------
        morph_filename : string
            Path to morphology.
        '''

        morph_basename = os.path.basename(morph_filename).decode()
        morph_extension = morph_basename.split('.')[-1]
        if morph_extension.lower() == 'swc':
            morph = self.h.Import3d_SWC_read()
        elif morph_extension.lower() == 'asc':
            morph = self.h.Import3d_Neurolucida3()
        else:
            raise Exception("Unknown filetype: %s" % morph_extension)

        morph.input(morph_filename)
        imprt = self.h.Import3d_GUI(morph, 0)

        self.h("objref this")
        imprt.instantiate(self.h.this)

        for sec in self.h.allsec():
            sec.nseg = 1 + 2 * int(sec.L / 40.0)

        self.h("soma[0] area(0.5)")
        axon_diams = [self.h.axon[0].diam, self.h.axon[0].diam]
        self.h.distance(sec=self.h.soma[0])
        for sec in self.h.allsec():
            if sec.name()[:4] == "axon":
                if self.h.distance(0.5, sec=sec) > 60:
                    axon_diams[1] = sec.diam
                    break
        for sec in self.h.allsec():
            if sec.name()[:4] == "axon":
                self.h.delete_section(sec=sec)
        self.h('create axon[2]')
        for index, sec in enumerate(self.h.axon):
            sec.L = 30
            sec.diam = axon_diams[index]

        for sec in self.h.allsec():
            sec.nseg = 1 + 2 * int(sec.L / 40.0)

        self.h.axon[0].connect(self.h.soma[0], 1.0, 0.0)
        self.h.axon[1].connect(self.h.axon[0], 1.0, 0.0)

        # make sure diam reflects 3d points
        self.h.area(.5, sec=self.h.soma[0])

    def load_cell_parameters(self):
        '''Configure a neuron after the cell morphology has been loaded.'''
        passive = self.description.data['passive'][0]
        genome = self.description.data['genome']
        conditions = self.description.data['conditions'][0]
        h = self.h

        h("access soma")

        # Set fixed passive properties
        for sec in h.allsec():
            sec.Ra = passive['ra']
            sec.insert('pas')
            # for seg in sec:
            #     seg.pas.e = passive["e_pas"]

        # for c in passive["cm"]:
        #     h('forsec "' + c["section"] + '" { cm = %g }' % c["cm"])

        # Insert channels and set parameters
        for p in genome:
            section_array = p["section"]
            mechanism = p["mechanism"]
            param_name = p["name"]
            param_value = float(p["value"])
            if section_array == "glob":  # global parameter
                h(p["name"] + " = %g " % p["value"])
            else:
                if hasattr(h, section_array):
                    if mechanism != "":
                        print('Adding mechanism %s to %s'
                              % (mechanism, section_array))
                        for section in getattr(h, section_array):
                            if self.h.ismembrane(str(mechanism),
                                                 sec=section) != 1:
                                section.insert(mechanism)

                    print('Setting %s to %.6g in %s'
                          % (param_name, param_value, section_array))
                    for section in getattr(h, section_array):
                        setattr(section, param_name, param_value)

        # Set reversal potentials
        for erev in conditions['erev']:
            erev_section_array = erev["section"]
            ek = float(erev["ek"])
            ena = float(erev["ena"])

            print('Setting ek to %.6g and ena to %.6g in %s'
                  % (ek, ena, erev_section_array))

            if hasattr(h, erev_section_array):
                for section in getattr(h, erev_section_array):
                    if self.h.ismembrane("k_ion", sec=section) == 1:
                        setattr(section, 'ek', ek)

                    if self.h.ismembrane("na_ion", sec=section) == 1:
                        setattr(section, 'ena', ena)
            else:
                print("Warning: can't set erev for %s, "
                      "section array doesn't exist" % erev_section_array)

        self.h.v_init = conditions['v_init']
        self.h.celsius = conditions['celsius']


def create_utils(description, model_type=None):
    ''' Factory method to create a Utils subclass.

    Parameters
    ----------
    description : Config instance
        used to initialize Utils subclass

    model_type : string
        Must be one of [PERISOMATIC_TYPE, ALL_ACTIVE_TYPE].  If none, defaults to PERISOMATIC_TYPE

    Returns
    -------
    Utils instance
    '''

    if model_type is None:
        try:
            model_type = description.data['biophys'][0]['model_type']
        except KeyError as e:
            logging.error("Could not infer model type from description")

    if model_type == PERISOMATIC_TYPE:
        return Utils(description)
    elif model_type == ALL_ACTIVE_TYPE:
        return AllActiveUtils(description)