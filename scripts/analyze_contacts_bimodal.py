#!/usr/bin/env python3
"""
DEM Contact Analysis - Bimodal mode (AM_P + AM_S + SE, 3 types)
Same pipeline as standard, just different default type map.

Usage:
    python analyze_contacts_bimodal.py results/atoms.csv results/contacts.csv \
        -o ./results -t "1:AM_P,2:AM_S,3:SE" -s 1000
"""
import sys
import os

# Reuse standard analyze_contacts with different default
sys.path.insert(0, os.path.dirname(__file__))
from analyze_contacts import main as standard_main

if __name__ == '__main__':
    # Override default type map in sys.argv if not specified
    if '-t' not in sys.argv and '--type-map' not in sys.argv:
        sys.argv.extend(['-t', '1:AM_P,2:AM_S,3:SE'])
    standard_main()
