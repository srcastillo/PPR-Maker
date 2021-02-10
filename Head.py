#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  28 2021

@author: Santiago Restrepo-Castillo & Gabriel Marti­nez-Galvez
"""

from opentrons import protocol_api

# Include experiment metadata
metadata = {
    'protocolName': 'FusR',
    'author': 'SRC & GMG',
    'description': 'This script results in the assembly of PPR proteins.',
    'apiLevel': '2.8'
}

# Set OT Labware

def run(protocol: protocol_api.ProtocolContext):

    # Labware
    R1 = protocol.load_labware('nest_96_wellplate_100ul_pcr_full_skirt', 4)
    R2 = protocol.load_labware('nest_96_wellplate_100ul_pcr_full_skirt', 5)
    R3 = protocol.load_labware('nest_96_wellplate_100ul_pcr_full_skirt', 6)
    MISC = protocol.load_labware('nest_96_wellplate_100ul_pcr_full_skirt', 2)
    tiprack_1 = protocol.load_labware('geb_96_tiprack_10ul', 1)
    tiprack_2 = protocol.load_labware('geb_96_tiprack_10ul', 3)
    tiprack_3 = protocol.load_labware('geb_96_tiprack_10ul', 7)
    tiprack_4 = protocol.load_labware('geb_96_tiprack_10ul', 8)

    # Pipettes
    p20 = protocol.load_instrument(
        'p20_single_gen2', 'right', tip_racks=
        [tiprack_1, tiprack_2, tiprack_3, tiprack_4])
